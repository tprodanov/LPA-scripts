#!/usr/bin/env python3

import liftover
import pysam
import gzip
import sys
from Bio.Seq import Seq


def open(filename, mode='r'):
    import builtins
    if filename == '-' or filename is None:
        return sys.stdin if mode == 'r' else sys.stdout

    if filename.endswith('.gz'):
        if 'b' not in mode:
            mode += 't'
        return gzip.open(filename, mode)
    return builtins.open(filename, mode)


class Columns:
    def __init__(self, line):
        # Convert to dictionary so that we can use `get` to fill missing columns with None.
        col_names = { col: i for i, col in enumerate(line.strip().split('\t')) }
        self.rsid = col_names.get('rsID')
        self.chrom = col_names.get('chrom_name') or col_names.get('hm_chr') or col_names.get('chr_name')
        self.pos = col_names.get('chr_position') or col_names.get('hm_pos')

        self.effect = col_names.get('effect_allele')
        self.other = col_names.get('other_allele')
        self.weight = col_names.get('effect_weight')


def convert_header(f_in, f_out, new_genome_name):
    for line in f_in:
        if line.startswith('#'):
            if line.startswith('#genome_build') or line.startswith('#HmPOS_build'):
                old_build = line.strip().split('=')[1]
                assert old_build in ('hg19', 'GRCh37'), f'Unexpected genome used in the PGS file ({old_build})'
                f_out.write(f'#genome_build={new_genome_name}\n')
            else:
                if line.startswith('#HmPOS_build'):
                    hm_build = line.strip().split('=')[1]
                    assert hm_build == 'GRCh37', f'Unexpected genome used in the PGS file ({hm_build})'
                f_out.write(line)

        else:
            f_out.write('##LIFT_OVER\n')
            f_out.write(f'#original_build={old_build}\n')
            return Columns(line)


def identify_ref_allele(reference, chrom, pos, effect_allele, other_allele):
    ref_seq = reference.fetch(chrom, pos - 1, pos + max(len(effect_allele), len(other_allele)))
    matches_effect = ref_seq.startswith(effect_allele)
    matches_other = ref_seq.startswith(other_allele)
    if matches_effect == matches_other:
        sys.stderr.write(f'[{chrom}:{pos}] Could not match reference sequence ({ref_seq}) against alleles '
            f'{effect_allele} and {other_allele}\n')
        return None
    return effect_allele if matches_effect else other_allele


def identify_variant(rsid, chrom, pos, allele1, allele2, variants):
    try:
        records = [rec for rec in variants.fetch(chrom, pos - 1, pos) if rec.pos == pos]
    except ValueError:
        return None, f'Could not fetch variants from {chrom}:{pos}'

    if len(records) == 1:
        return records[0], None
    for rec in records:
        if rec.id == rsid:
            return rec, None
    for rec in records:
        if allele1 in rec.alleles and (allele2 is None or allele2 in rec.alleles):
            return rec, None
    return None, f'Failed to identify appropriate variant ({len(records)} partial matches)'


# def process_wo_positions(f_in, f_out, header, converter, variants):
#     sys.stderr.write('Loading variants\n')
#     var_pos = {}
#     for rec in variants:
#         if rec.id:
#             var_pos[rec.id] = (rec.chrom, rec.pos, rec.alleles)

#     sys.stderr.write('Processing PGS list\n')
#     for line in f_in:
#         cols = line.rstrip('\n').split('\t')
#         rsid = cols[header.rsid]
#         weight = float(cols[header.weight])
#         effect_allele = cols[header.effect]
#         other_allele = cols[header.other] if header.other is not None else None

#         rec = var_pos.get(rsid)
#         if rec is None:
#             sys.stderr.write(f'ERROR [{rsid}] Variant absent from the VCF file\n')
#             continue

#         new_chrom = rec.chrom
#         new_pos = rec.pos

#         allele_i = rec.alleles.index(effect_allele) if effect_allele in rec.alleles else None
#         allele_j = rec.alleles.index(other_allele) if other_allele in rec.alleles else None
#         if allele_i is None or (allele_j is None and other_allele is not None):
#             sys.stderr.write(f'ERROR [{rsid} -> {new_chrom}:{new_pos}] '
#                 f'Failed to match alleles ({effect_allele} or {other_allele}) '
#                 f'to variant alleles ({", ".join(rec.alleles)}) [Possibly negative strand?]\n')
#             continue

# def check_ref():
    # ref_seq = reference.fetch(new_chrom, new_pos - 1,
    #     new_pos - 1 + max(len(effect_allele), len(other_allele))).upper()
    # effect_match = ref_seq.startswith(effect_allele)
    # other_match = ref_seq.startswith(other_allele)
    # if not effect_match and not other_match:
    #     sys.stderr.write(f'[{old_chrom}:{old_pos} "{rsid}" -> {new_chrom}:{new_pos}] '
    #         f'Both alleles ({effect_allele} and {other_allele}) do not match reference ({ref_seq})\n')
    #     continue
    # if effect_match and other_match:
    #     sys.stderr.write(f'[{old_chrom}:{old_pos} "{rsid}" -> {new_chrom}:{new_pos}] '
    #         f'Both alleles ({effect_allele} and {other_allele}) match reference ({ref_seq}) (still usable)\n')



def process(f_in, f_out, header, converter, variants):
    f_out.write('rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\t'
        'orig_rsID\torig_chr\torig_pos\torig_effect_allele\torig_other_allele\tliftover_strand\n')
    for line in f_in:
        cols = line.rstrip('\n').split('\t')
        rsid = cols[header.rsid] if header.rsid is not None else ''
        old_chrom = cols[header.chrom]
        old_pos = cols[header.pos]
        weight = cols[header.weight]
        try:
            old_pos = int(old_pos)
        except ValueError:
            sys.stderr.write('ERROR Missing or unexpected position ({})\n'.format(line.strip().replace('\t', '  ')))
            continue

        new_effect_allele = effect_allele = cols[header.effect]
        new_other_allele = other_allele = cols[header.other] if header.other is not None else None
        assert other_allele is None or ',' not in other_allele
        # ref_allele = identify_ref_allele(old_ref, old_chrom, old_pos, effect_allele, other_allele)

        converted = converter[old_chrom][old_pos]
        if len(converted) != 1:
            sys.stderr.write(f'ERROR [{old_chrom}:{old_pos} "{rsid}"] '
                f'Failed to lift-over to the new reference ({len(converted)} new positions)\n')
            continue

        new_chrom, new_pos, strand = converted[0]
        assert strand == '+' or strand == '-'
        if strand == '-':
            if other_allele is not None and len(effect_allele) != len(other_allele):
                sys.stderr.write(f'ERROR [{old_chrom}:{old_pos} "{rsid}" -> {new_chrom}:{new_pos}] '
                    f'Insert/deletion on negative strand\n')
                continue
            new_pos -= len(effect_allele) - 1
            new_effect_allele = str(Seq(effect_allele).reverse_complement())
            if other_allele is not None:
                new_other_allele = str(Seq(other_allele).reverse_complement())

        rec, msg = identify_variant(rsid, new_chrom, new_pos, new_effect_allele, new_other_allele, variants)
        if rec is None:
            sys.stderr.write(f'WARN [{old_chrom}:{old_pos} "{rsid}" -> {new_chrom}:{new_pos}] {msg}\n')
            new_rsid = ''
        else:
            new_rsid = rec.id
            allele_i = rec.alleles.index(new_effect_allele) if new_effect_allele in rec.alleles else None
            allele_j = rec.alleles.index(new_other_allele) \
                if new_other_allele is not None and new_other_allele in rec.alleles else None

            if allele_i is None or (other_allele is not None and allele_j is None):
                sys.stderr.write(f'WARN  [{old_chrom}:{old_pos} "{rsid}" -> {new_chrom}:{new_pos}] '
                    f'Failed to match alleles ({new_effect_allele} or {new_other_allele}) '
                    f'to variant alleles ({", ".join(rec.alleles)})\n')
            if allele_i != 0 and allele_j != 0:
                new_other_allele = rec.alleles[0]
            elif allele_i == 0 and other_allele is None:
                new_other_allele = rec.alleles[1]

        f_out.write(f'{new_rsid}\t{new_chrom}\t{new_pos}\t{new_effect_allele}\t{new_other_allele}\t{weight}\t'
            f'{rsid}\t{old_chrom}\t{old_pos}\t{effect_allele}\t{other_allele or "*"}\t{strand}\n')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input harmonized PGS file.')
    parser.add_argument('-o', '--output', metavar='FILE', default='-',
        help='Output PGS file.')
    parser.add_argument('-c', '--chain', metavar='FILE', required=True,
        help='Liftover chain.')
    # parser.add_argument('-r', '--old-ref', nargs=2, required=True,
    #     help='Original reference genome name and FASTA file.')
    parser.add_argument('-r', '--reference', default='CHM13v2.0',
        help='New reference genome name [%(default)s].')
    parser.add_argument('-v', '--variants', metavar='FILE', required=True,
        help='VCF file with variant definitions in the new reference genome.')
    args = parser.parse_args()

    converter = liftover.ChainFile(args.chain, one_based=True)
    variants = pysam.VariantFile(args.variants) if args.variants else None

    f_in = open(args.input)
    f_out = open(args.output, 'w')
    header = convert_header(f_in, f_out, args.reference)
    process(f_in, f_out, header, converter, variants)


if __name__ == '__main__':
    main()
