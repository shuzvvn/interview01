#!/usr/bin/python3

# question4.2.py

# Shu-Ting Cho <shutingcho@pitt.edu>
# 4.	附檔中的bed檔為library prep kit的bed檔(hg38), 請找出clinical_ann_metadata.tsv 欄位Location為rsid的位點, 是否為在bed檔的含蓋範圍, 並以Level of Evidence 分別繪出各個bed檔含蓋的位點數(POS), 基因數(Gene), 藥物數(Related chamicals), 如下


# v2 2021/04/19	process by ranges in bed
# v1 2021/04/18


from datetime import datetime

# input data
tsv_data = './clinical_ann_metadata.tsv'

bed_dict = {}
bed_dict['IDT'] = './IDT_targets.bed'
bed_dict['roche'] = './ROCHE_KAPA_targets.bed'
bed_dict['V6'] = './Agilent_V6_targets.bed'
bed_dict['V7'] = './Agilent_V7_targets.bed'


# output data
out_file = './out_data.tsv'


# make empty summary dict for each bed
levels = ['1A', '1B', '2A', '2B', '3', '4']
infos = ['Chemicals', 'Gene', 'pos']

bed_summary_dict = {}

for info in infos:
	bed_summary_dict[info] = {}
	for bed in bed_dict:
		bed_summary_dict[info][bed] = {}
		for level in levels:
			bed_summary_dict[info][bed][level] = []


# open and read tsv data
tsv_data_h = open(tsv_data, 'r')
lines_data = tsv_data_h.readlines()

# get index for info needed:
header = lines_data[0].strip('\n').split('\t')
i_pos = header.index("Location")
i_level = header.index("Level of Evidence")
i_chemicals = header.index("Related Chemicals")
i_gene = header.index("Gene")


# read bed files line by line, store as list of postions
for bed in bed_dict:
	with open(bed_dict[bed], 'r') as in_file_h:
		n_line = 0
		print(datetime.now().time(), 'processing', bed)
		while True:
			line_bed = in_file_h.readline()
			n_line += 1
			if not line_bed:
				break
			else:
				words_bed = line_bed.strip('\n').split('\t')
				range_h = range(int(words_bed[1]), int(words_bed[2]))
				for line_data in lines_data[1:]:
					words_data = line_data.strip('\n').split('\t')
					if words_data[i_pos][:2] == 'rs':
						info_dict = {} # store infos in dict
						info_dict['pos'] = words_data[i_pos][2:]
						#print(datetime.now().time(), 'processing', info_dict['pos'])
						# check if position is in bed
						if int(info_dict['pos']) in range_h:
							print(datetime.now().time(), info_dict['pos'], 'in', bed, n_line)
							info_dict['level'] = words_data[i_level]
							info_dict['Chemicals'] = words_data[i_chemicals]
							info_dict['Gene'] = words_data[i_gene]
							for info in infos:
								bed_summary_dict[info][bed][level] += info_dict[info]

tsv_data_h.close()


# transform lists in bed_summary_dict to numbers of unique elements
print(datetime.now().time(), 'writing output')
with open(out_file, 'w') as out_file_h:
	for info in bed_summary_dict:
		for bed in bed_summary_dict[info]:
			for level in bed_summary_dict[info][bed]:
				value = len(set(bed_summary_dict[info][bed][level]))
				words = [info, level, bed, str(value)]
				out_file_h.write('\t'.join(words) + '\n')