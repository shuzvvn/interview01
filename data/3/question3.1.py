#!/usr/bin/python3

# question3.1.py

# Shu-Ting Cho <shutingcho@pitt.edu>
# 3.	請從A_gen.fasta檔案中分別輸出HLA:HLA25681 HLA:HLA05443 HLA:HLA26560 HLA:HLA26769 HLA:HLA26076中第251-270之間的序列 
# 並請以任意工具以柱狀圖比較5條序列中ACGT的個數


# v1 2021/04/18

# libraries
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


# input
in_file = './A_gen.fasta'

# seq need
id_list = ['HLA:HLA25681', 'HLA:HLA05443', 'HLA:HLA26560', 'HLA:HLA26769', 'HLA:HLA26076']


dict_h = {}
for record in SeqIO.parse(open(in_file, "r"), "fasta") :
	if record.id in id_list:
		seq_h = str(record.seq[250:270])
		print(record.id + '\t' + seq_h)
		countATCG = []
		for i in ['A', 'T', 'C', 'G']:
			countATCG.append(seq_h.count(i))
		dict_h[record.id] = countATCG



## show/compare ATCG count with grouped barplot

fig = plt.figure()

# set width of bars
barWidth = 0.4

# Set position of bar on X axis
rA = np.arange(len(dict_h)) * 2
rT = [x + barWidth for x in rA]
rC = [x + barWidth for x in rT]
rG = [x + barWidth for x in rC]

# set heights of bars
barsA = [dict_h[x][0] for x in dict_h]
barsT = [dict_h[x][1] for x in dict_h]
barsC = [dict_h[x][2] for x in dict_h]
barsG = [dict_h[x][3] for x in dict_h]

# Make the plot
plt.bar(rA, barsA, color='green', width=barWidth, edgecolor='white', label='A')
plt.bar(rT, barsT, color='red', width=barWidth, edgecolor='white', label='T')
plt.bar(rC, barsC, color='blue', width=barWidth, edgecolor='white', label='C')
plt.bar(rG, barsG, color='maroon', width=barWidth, edgecolor='white', label='G')

# Add xticks on the middle of the group bars
plt.xlabel('sequence ID', fontweight='bold')
plt.xticks([x + barWidth*1.5 for x in rA], list(dict_h), fontsize=8)
 
# Create legend & Show graphic
plt.legend()
plt.show()
fig.savefig('./out.pdf')