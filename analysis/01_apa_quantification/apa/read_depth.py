# Write per-sample wig path and total mapped reads (from flagstat) for DaPars2.
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--sample_list', help='tab file; column 1 = sample ID')
parser.add_argument('--path_flagstat', default='./tmp', help='folder with <sample>.flagstat')
parser.add_argument('--path_wig', default='./wig', help='folder with <sample>.wig')
parser.add_argument('--output', help='output: wig path and read depth per sample')
args = parser.parse_args()


def extract_total_reads(flagstat_file):
    # first "mapped (" line (or "primary mapped") gives the mapped read count
    total_reads = '0'
    try:
        with open(flagstat_file, 'r') as f:
            for line in f:
                if ('mapped (' in line and 'primary' not in line) or 'primary mapped' in line:
                    total_reads = line.strip().split(' ')[0]
                    break
    except Exception as e:
        print(f'Error reading {flagstat_file}: {e}')
    return total_reads


samples = [line.strip().split('\t')[0] for line in open(os.path.abspath(args.sample_list), 'r')]
path_wig = os.path.join(os.path.abspath(args.path_wig), '')
path_flagstat = os.path.join(os.path.abspath(args.path_flagstat), '')

with open(args.output, 'w') as fho:
    for sample in samples:
        read_depth = extract_total_reads(path_flagstat + sample + '.flagstat')
        print('%s\t%s' % (path_wig + sample + '.wig', read_depth), file=fho)
