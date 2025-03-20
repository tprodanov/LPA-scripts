#!/usr/bin/env python3

import pysam
from collections import deque
from tqdm import tqdm


def var_matches(v, u):
    return v.chrom == u.chrom and v.pos == u.pos and v.alleles == u.alleles


def var_key(var):
    return (var.id, var.alts)


class Inputs:
    def __init__(self, filenames, region):
        self.vcfs = [pysam.VariantFile(fname) for fname in filenames]
        self.iters = [vcf.fetch(region=region) for vcf in self.vcfs]
        self.cache = [[] for it in self.iters]
        self.next = [next(it, None) for it in self.iters]

    def _load_next(self):
        it = self.iters[0]
        cache = self.cache[0]
        u = self.next[0]

        if u is None:
            assert all(u1 is None for u1 in self.next)
            return False

        cache.append(u)
        while True:
            v = next(it, None)
            if v is not None and u.chrom == v.chrom and u.pos == v.pos:
                cache.append(v)
            else:
                self.next[0] = v
                break
        cache.sort(key=var_key, reverse=True)
        n = len(cache) - 1

        for i, (it, cache) in enumerate(zip(self.iters[1:], self.cache[1:]), 1):
            cache.append(self.next[i])
            for _ in range(n):
                cache.append(next(it, None))
            self.next[i] = next(it, None)
            cache.sort(key=var_key, reverse=True)
        return True

    def __iter__(self):
        return self

    def __next__(self):
        if not self.cache[0]:
            if not self._load_next():
                raise StopIteration
        recs = [cache.pop() for cache in self.cache]
        assert all(rec is not None for rec in recs)
        return recs

    def debug(self):
        for vcf, u, cache in zip(self.vcfs, self.next, self.cache):
            print(f'    {vcf.header.samples[0]}')
            for v in cache[::-1]:
                print(f'        Cache: {v.chrom}:{v.pos}  {v.ref}  {",".join(v.alts)}')
            print(f'        Next:  {u.chrom}:{u.pos}  {u.ref}  {",".join(u.alts)}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output merged VCF file.')
    parser.add_argument('-r', '--region', metavar='STR',
        help='Merge specific region.')
    parser.add_argument('input', metavar='FILE', nargs='+',
        help='Input VCF files with the same variants by different samples.')
    args = parser.parse_args()

    inputs = Inputs(args.input, args.region)
    header = inputs.vcfs[0].header.copy()
    assert len(header.samples) == 1
    samples = { header.samples[0] }
    for f in inputs.vcfs[1:]:
        assert len(f.header.samples) == 1
        header.add_sample(f.header.samples[0])

    output = pysam.VariantFile(args.output, 'w', header=header)
    for recs in tqdm(inputs):
        rec0 = recs[0]
        assert all(var_matches(rec0, rec) for rec in recs[1:]), \
            f'Variant does not match at {rec0.chrom}:{rec0.pos}'
        new_rec = header.new_record()
        new_rec.chrom = rec0.chrom
        new_rec.pos = rec0.pos
        new_rec.id = rec0.id
        new_rec.alleles = rec0.alleles
        new_rec.qual = rec0.qual
        for filt in rec0.filter:
            new_rec.filter.add(filt)
        for k, v in rec0.info.items():
            new_rec.info[k] = v
        for i, rec in enumerate(recs):
            for k, v in rec.samples[0].items():
                new_rec.samples[i][k] = v
        output.write(new_rec)


if __name__ == '__main__':
    main()
