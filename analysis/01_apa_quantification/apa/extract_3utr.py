# Extract 3'UTR location (chr, start, end) per gene from a DaPars2 result table (column 4).
# Usage: python extract_3utr.py --dapars_res Dapars2_res.all_chromosomes.txt --output 3utr_location.txt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dapars_res', type=str, help='merged DaPars2 output (all chromosomes)')
parser.add_argument('--output', type=str, default='3utr_location.txt', help="3' UTR location file")
args = parser.parse_args()

gene2loc = {}
gene_order = []
with open(args.dapars_res, 'r') as fh:
    next(fh)
    for line in fh:
        w = line.strip().split('\t')
        chrom, loc_code = w[3].split(':')
        start, end = loc_code.split('-')
        gene2loc[w[0]] = (chrom, start, end)
        gene_order.append(w[0])
print(len(gene_order), 'genes have been processed!')

with open(args.output, 'w') as fho:
    print('Gene\tChr\tStart\tEnd', file=fho)
    for gene in gene_order:
        if gene in gene2loc:
            print('%s\t%s\t%s\t%s' % (gene, gene2loc[gene][0], gene2loc[gene][1], gene2loc[gene][2]), file=fho)
print('Done!')
