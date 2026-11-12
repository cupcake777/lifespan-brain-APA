# Extract 3' UTR location (chr, start, end) per gene from the DaPars2 result table (column 4).
# Usage: python extract_3utr.py --dapars_res Dapars2_res.all_chromosomes.txt --output 3utr_location.txt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dapars_res', type=str, help='DaPars2 output file; merged all-chromosome file if run per chromosome')
parser.add_argument('--output', type=str, default='3utr_location.txt', help="specify 3' UTR location file")
args = parser.parse_args()

# --- load gene locations ---
gene2loc = {}
gene_order = []
fh = open(args.dapars_res, 'r')
fh.readline()
for line in fh.readlines():
    line = line.strip()
    w = line.split('\t')
    chrom, loc_code = w[3].split(':')
    start, end = loc_code.split('-')
    gene2loc[w[0]] = (chrom, start, end)
    gene_order.append(w[0])
fh.close()
print(len(gene_order), 'genes have been processed!')

# --- write ---
with open(args.output, 'w') as fho:
    print('Gene\tChr\tStart\tEnd', file=fho)
    for gene in gene_order:
        print('%s\t%s\t%s\t%s' % (gene, gene2loc[gene][0], gene2loc[gene][1], gene2loc[gene][2]), file=fho)
print('Done!')
