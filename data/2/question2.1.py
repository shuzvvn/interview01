#!/usr/bin/python3

# question2.1.py

# Shu-Ting Cho <shutingcho@pitt.edu>

# 2.	請根據 reference.csv.gz 的 probe 資訊，計算 input.tped 中，每個 probe 的統計結果，輸出成 TXT
# a.	備註
# i.	reference.csv.gz 中，Probe Set ID 欄位代表 Probe_id
# ii.	reference.csv.gz 中，Alt Allele 有多種可能時，會以 "/" 當分隔號表示
# iii.	input.tped 格式中，第2欄為 probe_id，第5欄開始為 sample call
# iv.	input.tped 中，sample call 如果為 "0 0" 代表 No-call
# b.	輸出欄位
# i.	probe_id
# 1.	轉大寫
# ii.	call_rate
# 1.	不是 No call 的樣本數 佔多少
# iii.	allele_count
# 1.	不是 No call 的樣本數
# iv.	allele_frequency
# 1.	alt allele 佔多少 %
# v.	heterogeneity
# 1.	alt allele == 1 的有多少
# vi.	homogeneity
# 1.	alt allele == 2 的有多少


import sys
import argparse
import gzip


# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Statistics of probes according to reference.')
	parser.add_argument('-r', '--reference', help='reference file: .csv.gz')
	parser.add_argument('-i', '--input', help='input file: .tped')
	parser.add_argument('-o', '--output', help='output file: .txt')
	return parser.parse_args()

# read ref to dict
def ref2dict(ref_file):
	with gzip.open(ref_file, 'r') as in_file_h:	
		# get index of Probe_id and Alt Allele from header
		line = in_file_h.readline()
		line = line.decode('utf-8').rstrip()
		headers = line.split(',')
		for header_i in range(len(headers)):
			if headers[header_i] == "Probe Set ID":
				probe_id_i = header_i # reference.csv.gz 中，Probe Set ID 欄位代表 Probe_id
			elif headers[header_i] == "Alt Allele":
				alt_allele_i = header_i
		ref_dict = {}
		while True:
			line = in_file_h.readline()
			if not line:
				break
			else:
				line = line.decode('utf-8').rstrip()
				words = line.split(',')
				probe_id = words[probe_id_i].upper() # 轉大寫
				alt_allele = words[alt_allele_i].split('/') # Alt Allele 有多種可能時，會以 "/" 當分隔號表示
				ref_dict[probe_id] = alt_allele
	return ref_dict

# main
def main():
	args = parse_args()

	# read/open input, ref and output
	ref_dict = ref2dict(args.reference)

	output = open(args.output, 'w')
	# write header to output
	header = ['probe_id', 'call_rate', 'allele_count', 'allele_frequency', 'heterogeneity', 'homogeneity']
	output.write('\t'.join(header) + '\n')

	count_in = 0
	with open(args.input, 'r') as in_file_h:
		# process line by line
		while True:
			line = in_file_h.readline()
			if not line:
				break
			else:
				count_in += 1
				words = line.rstrip().split('\t')
				probe_id_h = words[1].upper() # input.tped 格式中，第2欄為 probe_id, 轉大寫 
				sample_count = 0
				allele_count = 0
				heterogeneity = 0 # alt allele == 1 的有多少
				homogeneity = 0 # alt allele == 2 的有多少

				for i in range(4, len(words)): # 第5欄開始為 sample call
					sample_count += 1
					if words[i] != "0 0":
						allele_count += 1 # sample call 如果為 "0 0" 代表 No-call
						calls = words[i].split(' ')
						n_alt_allele = 0
						for call in calls: # check number of alt allele
							if call in ref_dict[probe_id_h]:
								n_alt_allele += 1
						if n_alt_allele == 1:
							heterogeneity += 1
						elif n_alt_allele == 2:
							homogeneity += 1

				# call_rate: 不是 No call 的樣本數 佔多少 %
				call_rate = '{:.0%}'.format(allele_count / sample_count)

				# allele_frequency: alt allele 佔多少 %
				allele_frequency = '{:.0%}'.format((heterogeneity + homogeneity) / allele_count)

				row = [probe_id_h, call_rate, str(allele_count), allele_frequency, str(heterogeneity), str(homogeneity)]
				output.write('\t'.join(row) + '\n')
	output.close()

if __name__ == '__main__':
	main()