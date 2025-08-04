#!/usr/bin/env python3

from liftover import open
from tqdm import tqdm


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='FILE',
        help='Input scoring file.')
    parser.add_argument('-o', '--output', metavar='FILE', default='-',
        help='Output sum scores.')
    args = parser.parse_args()

    EFFECT_COL = 4
    SAMPLE_COL = 5
    with open(args.input) as inp:
        colnames = next(inp)
        assert not colnames.startswith('#')
        colnames = colnames.rstrip().split('\t')
        assert colnames[EFFECT_COL] == 'effect'
        samples = colnames[SAMPLE_COL:]
        scores = [0] * len(samples)

        for line in tqdm(inp):
            cols = line.rstrip().split('\t')
            effect = float(cols[EFFECT_COL])
            for i, dose in enumerate(cols[SAMPLE_COL:]):
                if dose != 'NA' and dose != 'nan':
                    scores[i] += effect * float(dose)

    with open(args.output, 'w') as out:
        out.write('sample\tscore\n')
        for sample, score in zip(samples, scores):
            out.write(f'{sample}\t{score:.10g}\n')


if __name__ == '__main__':
    main()
