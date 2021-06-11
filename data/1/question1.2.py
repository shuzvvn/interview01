#!/usr/bin/python3

# question1.2.py

# Shu-Ting Cho <shutingcho@pitt.edu>
# 1.	請將 test.vcf.gz 以 Transcript 為單位轉成 TSV 格式(tab 分隔)
# 	a.	備註
# 		i.	VCF 格式中，會以 ; 作為註解的分隔號
# 		ii.	CSQ 註解中會以 , 作為 Transcript 的分隔號
# 		iii.	CSQ 註解中會以 | 作為 每個欄位的分隔號
# 		iv.	CSQ 欄位的標頭可從 VCF header 中，CSQ 描述中取得
# 		v.	略過沒有 CSQ 註解的部分
# 	b.	輸出欄位
# 		i.	CHROM
# 		ii.	POS
# 		iii.	REF
# 		iv.	ALT
# 		v.	FILTER
# 		vi.	Consequence
# 		vii.	SYMBOL
# 		viii.	Gene
# 		ix.	Feature
# 		x.	HGVSc
# 		xi.	HGVSp



# v2 2021/04/18 without library: vcf, csv 
# v1 2021/04/17



import sys
import argparse
import gzip


# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Convert .vcf.gz to .tsv by transcript.')
	parser.add_argument('-v', '--vcf', help='input file: .vcf.gz')
	parser.add_argument('-o', '--output', help='output file: .tsv')
	return parser.parse_args()

# to cover functions in vcf library
def vcf_line_reader(line):
	vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '	WT200617170161']
	words = line.split('\t')
	vcf_line_dict = {}
	for i in range(len(vcf_header)):
		vcf_line_dict[vcf_header[i]] = words[i]
	return vcf_line_dict

def vcf_get_CSQ_from_INFO(info):
	words = info.split(';')
	for word in words:
		if word[0:4] == 'CSQ=':
			record_INFO_CSQ = word[4:].split(',')
	return record_INFO_CSQ

# main
def main():
	args = parse_args()

	# read/open input and output
	vcf_file = gzip.open(args.vcf, 'r')

	output = open(args.output, 'w')

	# write header to output
	header = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'Consequence', 'SYMBOL', 'Gene', 'Feature', 'HGVSc', 'HGVSp']
	output.write('\t'.join(header) + '\n')

	# process line by line
	while True:
		line = vcf_file.readline()

		if not line:
			break
		else:
			line = line.decode('utf-8').rstrip()
		
			# ignore lines start with '#'
			if line[0] != '#':
				record = vcf_line_reader(line)

				# ignore position without CSQ info
				try:
					record_INFO_CSQ = vcf_get_CSQ_from_INFO(record['INFO'])
					for transcript in record_INFO_CSQ:
						words = transcript.split('|')

						# get info according to CSQ Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position ...
						row = [record['CHROM'], record['POS'], record['REF'], record['ALT'], record['FILTER'], words[1], words[3], words[4], words[6], words[10], words[11]]
						output.write('\t'.join(row) + '\n')
				except:
					print(record['CHROM'], record['POS'], 'has no CSQ info')
	output.close()

if __name__ == '__main__':
	main()