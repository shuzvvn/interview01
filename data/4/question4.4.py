#!/usr/bin/python3

# question4_calculate.4.py

# Shu-Ting Cho <shutingcho@pitt.edu>
# 4.	附檔中的bed檔為library prep kit的bed檔(hg38), 請找出clinical_ann_metadata.tsv 欄位Location為rsid的位點, 是否為在bed檔的含蓋範圍, 並以Level of Evidence 分別繪出各個bed檔含蓋的位點數(POS), 基因數(Gene), 藥物數(Related chamicals), 如下

# v4 2021/04/19 store pos in dict for fast searching
# v3 2021/04/19 try with diff search approach
# v2 2021/04/19 in python
# v1 2021/04/18


#from datetime import datetime

# input data
tsv_data = './clinical_ann_metadata.tsv'

bed_dict = {}
bed_dict['IDT'] = './IDT_targets.bed'
bed_dict['roche'] = './ROCHE_KAPA_targets.bed'
bed_dict['V6'] = './Agilent_V6_targets.bed'
bed_dict['V7'] = './Agilent_V7_targets.bed'


# output data
out_file = './out_data.tsv'


# read bed files line by line, store as list of postions
pos_dict = {}

for kit in bed_dict:
	with open(bed_dict[kit], 'r') as in_file_h:
		#print(datetime.now().time(), 'processing', kit)
		while True:
			line = in_file_h.readline()
			if not line:
				break
			else:
				words = line.strip('\n').split('\t')
				for pos_h in range(int(words[1])+1, int(words[2])+1):
					try:
						pos_dict[pos_h].append(kit)
					except:
						pos_dict[pos_h] = [kit]


# make empty summary dict for each bed
levels = ['1A', '1B', '2A', '2B', '3', '4']
infos = ['Chemicals', 'Gene', 'pos']

bed_summary_dict = {}

for info in infos:
	bed_summary_dict[info] = {}
	for level in levels:
		bed_summary_dict[info][level] = {}
		for kit in bed_dict:
			bed_summary_dict[info][level][kit] = []


# process tsv data line by line
with open(tsv_data, 'r') as in_file_h:
	#print(datetime.now().time(), 'processing tsv data')
	# get index for info needed: 
	header = in_file_h.readline().strip('\n').split('\t')
	i_pos = header.index("Location")
	i_level = header.index("Level of Evidence")
	i_chemicals = header.index("Related Chemicals")
	i_gene = header.index("Gene")
	while True:
		line = in_file_h.readline()
		if not line:
			break
		else:
			words = line.strip('\n').split('\t')
			# check if location is rsID (rs[0-9]+)
			if words[i_pos][:2] == 'rs':
				pos = words[i_pos][2:]
				# check if position is in bed
				if int(pos) in pos_dict:
					for kit in set(pos_dict[int(pos)]):
						#print(datetime.now().time(), pos, 'in', kit)
						level = words[i_level]
						bed_summary_dict['Chemicals'][level][kit].append(words[i_chemicals])
						bed_summary_dict['Gene'][level][kit].append(words[i_gene])
						bed_summary_dict['pos'][level][kit].append(pos)


# transform lists in bed_summary_dict to numbers of unique elements
#print(datetime.now().time(), 'writing output')
with open(out_file, 'w') as out_file_h:
	# write header
	header = ['info', 'Evidence level', 'kit', 'Counts']
	out_file_h.write('\t'.join(header) + '\n')
	for info in bed_summary_dict:
		for level in bed_summary_dict[info]:
			for kit in bed_summary_dict[info][level]:
				value = len(set(bed_summary_dict[info][level][kit]))
				words = [info, level, kit, str(value)]
				out_file_h.write('\t'.join(words) + '\n')