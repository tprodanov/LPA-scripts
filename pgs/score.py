#!/usr/bin/env python3

import itertools
import sys
from collections import namedtuple, deque
import pysam
from tqdm import tqdm
import numpy as np

from liftover import open


class PgsEntry:
    def __init__(self, rsid, chrom, pos, effect_allele, other_allele, effect, contigs):
        self.rsid = rsid
        self.chrom = chrom
        self.chrom_id = contigs.get(chrom, sys.maxsize)
        self.start = int(pos) - 1
        self.effect_allele = effect_allele
        self.other_allele = other_allele
        self.ref_allele = None
        self.effect = effect

        self.assumed_end = self.start + max(len(self.effect_allele), len(self.other_allele))

    @property
    def pos(self):
        return self.start + 1

    def obtain_ref_allele(self, reference):
        ref_seq = reference.fetch(self.chrom, self.start, self.assumed_end).upper()
        matches_effect = ref_seq.startswith(self.effect_allele)
        matches_other = ref_seq.startswith(self.other_allele)
        if matches_effect and matches_other:
            sys.stderr.write(f'WARN [{self.chrom}:{self.pos} "{self.rsid}"] Both alleles {self.effect_allele} '
                f'and {self.other_allele} match the reference {ref_seq}\n')
        elif matches_effect:
            self.ref_allele = self.effect_allele
            self.assumed_end = self.start + len(self.ref_allele)
        elif matches_other:
            self.ref_allele = self.other_allele
            self.assumed_end = self.start + len(self.ref_allele)
        else:
            sys.stderr.write(f'WARN [{self.chrom}:{self.pos} "{self.rsid}"] Both alleles {self.effect_allele} '
                f'and {self.other_allele} do not match the reference {ref_seq}\n')

    def __str__(self):
        return f'{self.rsid}\t{self.chrom}\t{self.pos}\t{self.effect_allele}'

    def __lt__(self, oth):
        return (self.chrom_id, self.chrom, self.start) < (oth.chrom_id, oth.chrom, oth.start)


class PGS:
    def __init__(self, filename, reference, contigs):
        self.reference = reference
        self.contigs = contigs
        self.f = open(filename)
        for line in self.f:
            if line.startswith('#'):
                continue
            colnames = line.rstrip('\n').split('\t')
            break
        self.rsid_col = colnames.index('rsID')
        self.chrom_col = colnames.index('chr_name')
        self.pos_col = colnames.index('chr_position')
        self.allele1_col = colnames.index('effect_allele')
        self.allele2_col = colnames.index('other_allele')
        self.effect_col = colnames.index('effect_weight')

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.f).rstrip('\n').split('\t')
        entry = PgsEntry(
            line[self.rsid_col], line[self.chrom_col],
            line[self.pos_col], line[self.allele1_col], line[self.allele2_col],
            line[self.effect_col], self.contigs)
        entry.obtain_ref_allele(self.reference)
        return entry


def process(pgs, vcf, out, log):
    cache = deque()
    for entry in tqdm(pgs):
        # sys.stderr.write(f'Processing {entry}\n')
        # while True:
        #     # sys.stderr.write(f'    Cache size {len(cache)}\n')
        #     if not cache:
        #         new_var = next(vcf, None)
        #         if new_var is None:
        #             break
        #         # sys.stderr.write(f'    Append {new_var.chrom}:{new_var.pos}\n')
        #         cache.append(new_var)

        #     var = cache[0]
        #     var_end = var.start + len(var.ref)
        #     # sys.stderr.write(f'    Check  {var.chrom}:{var.pos}-{var_end}\n')
        #     if entry.chrom_id < var.rid or entry.assumed_end <= var.start:
        #         # sys.stderr.write(f'    Break\n')
        #         break
        #     if var.rid < entry.chrom_id or var_end <= entry.start:
        #         cache.popleft()
        #         # sys.stderr.write(f'    Pop left\n')
        #         continue

        # vars = [var for var in cache if var.start <= entry.start and entry.assumed_end <= var_end]
        # sys.stderr.write(f'    Select {len(vars)} variants\n')

        if entry.chrom_id < sys.maxsize:
            vars = [var for var in vcf.fetch(entry.chrom, entry.start, entry.assumed_end)
                if var.start <= entry.start and entry.assumed_end <= var.start + len(var.ref)]
        else:
            vars = ()
        score_entry(entry, vcf.header.samples, vars, out, log)


def identify_allele(entry, var):
    padd_left = entry.start - var.start
    padd_right = var.start + len(var.ref) - entry.assumed_end
    for i, allele in enumerate(var.alleles):
        m = len(allele)
        if m > padd_left + padd_right and allele[padd_left : m - padd_right] == entry.effect_allele:
            return i
    return None


def separate_vars(vars):
    sep_vars = []
    for i, v in enumerate(vars):
        if i and sep_vars[-1][-1].start == v.start and sep_vars[-1][-1].ref == v.ref:
            sep_vars[-1].append(v)
        else:
            sep_vars.append([v])
    return sep_vars


def match_var(entry, vars, log):
    sep_vars = separate_vars(vars)
    best_i = None
    best_key = None
    best_allele_ixs = None

    for i, subvars in enumerate(sep_vars):
        ref_effect = False
        alt_effect = False
        allele_ixs = []
        for var in subvars:
            allele_ix = identify_allele(entry, var)
            log.write(f'\n{entry}\t{var.pos}\t{",".join(var.alleles)}\t')
            if allele_ix is not None:
                ref_effect |= allele_ix == 0
                alt_effect |= allele_ix != 0
                log.write(str(allele_ix))
            else:
                log.write('*')
            allele_ixs.append(allele_ix)

        if ref_effect and alt_effect:
            log.write('\tref&alt')
            continue
        elif not ref_effect and not alt_effect:
            continue

        key = (len(subvars), len(subvars[0].ref))
        if best_key is None or key < best_key:
            best_i = i
            best_key = key
            best_allele_ixs = allele_ixs

    if best_i is None:
        log.write('\tmismatch')
        return None
    log.write('\tuse')
    subvars = sep_vars[best_i]
    if len(sep_vars) > 1:
        u = subvars[0]
        log.write(f'[{u.pos},{u.ref}]')
    return subvars, best_allele_ixs


def score_entry(entry, samples, vars, out, log):
    if not vars:
        log.write(f'\n{entry}\t*\t*\tNA\tabsent')
    else:
        match = match_var(entry, vars, log)

    match = match_var(entry, vars, log)
    if match is None:
        log.write(f'\n{entry}\t*\t*\tNA\tabsent')
        out.write(f'{entry}\t{entry.effect}')
        out.write('\tNA' * len(samples))
        out.write('\n')
        return

    subvars, allele_ixs = match
    dose = [0] * len(samples)
    for var, ix in zip(subvars, allele_ixs):
        if ix is None:
            continue
        for i, (sample, gt) in enumerate(var.samples.items()):
            assert sample == samples[i]
            dose[i] += gt['GT'].count(ix)
    if allele_ixs[0] == 0:
        # alt_count = 2 * len(subvars) - ref_count
        # Dose = 2 - alt_count = 2 - 2 * len(subvars) + ref_count.
        c = 2 - 2 * len(subvars)
        dose = [max(0, d + c) for d in dose]
    dose = [min(d, 2) for d in dose]
    dose_str = '\t'.join(map(str, dose))
    out.write(f'{entry}\t{entry.effect}\t{dose_str}\n')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', metavar='FILE', required=True,
        help='Input PGS score for a correct reference genome.')
    parser.add_argument('-v', '--vcf', metavar='FILE', required=True,
        help='Multi-sample VCF file.')
    parser.add_argument('-r', '--reference', metavar='FILE', required=True,
        help='Genome reference.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-l', '--log', metavar='FILE', required=True,
        help='Output CSV file describing variant matching.')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    contigs = { contig: i for i, contig in enumerate(vcf.header.contigs) }
    reference = pysam.FastaFile(args.reference)
    sys.stderr.write('Loading PGS file\n')
    pgs = PGS(args.score, reference, contigs)
    entries = sorted(pgs)
    sys.stderr.write(f'Loaded {len(entries)} entries\n')

    out = open(args.output, 'w')
    samples_str = '\t'.join(vcf.header.samples)
    out.write(f'rsid\tchrom\tpos\teffect_allele\teffect\t{samples_str}\n')
    log = open(args.log, 'w')
    log.write('rsid\tchrom\tpos\teffect_allele\tvar_pos\tvar_alleles\teffect_ix\tstatus')

    sys.stderr.write('Processing VCF file\n')
    process(entries, vcf, out, log)
    log.write('\n')


if __name__ == '__main__':
    main()
