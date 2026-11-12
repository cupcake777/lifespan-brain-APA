# Write per-sample wig path and total mapped reads (from flagstat) for DaPars2.
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('--sample_list', help='the file contains all samples')
parser.add_argument('--path_flagstat', default='./tmp', help='the location of flagstat files')
parser.add_argument('--path_wig', default='./wig', help='the location of wig files')
parser.add_argument('--output', help='the final output file with read depth of each sample')
args = parser.parse_args()

# --- sample list and folders ---
selected_samples = [line.strip().split('\t')[0] for line in open(os.path.abspath(args.sample_list), 'r')]
path_wig = os.path.join(os.path.abspath(args.path_wig), '')
path_flagstat = os.path.join(os.path.abspath(args.path_flagstat), '')

# --- mapped reads per sample ---
with open(args.output, 'w') as fho:
    for sample in selected_samples:
        # first "mapped (" line (or "primary mapped" line) gives the mapped read count
        read_depth = '0'
        flagstat_file = path_flagstat + sample + '.flagstat'
        try:
            with open(flagstat_file, 'r') as f:
                for line in f:
                    if ('mapped (' in line and 'primary' not in line) or 'primary mapped' in line:
                        read_depth = line.strip().split(' ')[0]
                        break
        except Exception as e:
            print(f'Error reading {flagstat_file}: {e}')
        wig_location = path_wig + sample + '.wig'
        print('%s\t%s' % (wig_location, read_depth), file=fho)
