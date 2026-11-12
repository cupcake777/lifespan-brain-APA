# Write the DaPars2 configure file from the 3'UTR bed and the wig/read-depth table.
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--annotation_3utr', help="reference 3'UTR bed file")
parser.add_argument('--wigFile_depth', help='wig files and read depth table')
parser.add_argument('--coverage_threshold', type=str, default='10')
parser.add_argument('--threads', type=str, default='1')
parser.add_argument('--out_dir_prefix', type=str, default='Dapars2_out')
parser.add_argument('--out_file_prefix', type=str, default='Dapars2')
parser.add_argument('--out_config_name', type=str, default='Dapars2_running_configure.txt')
args = parser.parse_args()

wig_files = [line.strip().split('\t')[0] for line in open(os.path.abspath(args.wigFile_depth), 'r')]

with open(args.out_config_name, 'w') as fho:
    print("# Specify the reference of 3'UTR region", file=fho)
    print('\nAnnotated_3UTR=' + os.path.abspath(args.annotation_3utr), file=fho)
    print('\n# A comma separated list of wig files of all samples', file=fho)
    print('\nAligned_Wig_files=' + ','.join(wig_files), file=fho)
    print('\nOutput_directory=' + args.out_dir_prefix, file=fho)
    print('\nOutput_result_file=' + args.out_file_prefix, file=fho)
    print('\n# Specify Coverage threshold', file=fho)
    print('\nCoverage_threshold=' + args.coverage_threshold, file=fho)
    print('\n# Specify the number of threads to process the analysis', file=fho)
    print('\nNum_Threads=' + args.threads, file=fho)
    print('\n# Provide sequencing depth file for normalization', file=fho)
    print('\nsequencing_depth_file=' + os.path.abspath(args.wigFile_depth), file=fho)