#!/usr/bin/env Rscript

# question4.1.R

# Shu-Ting Cho <shutingcho@pitt.edu>
# 4.	附檔中的bed檔為library prep kit的bed檔(hg38), 請找出clinical_ann_metadata.tsv 欄位Location為rsid的位點, 是否為在bed檔的含蓋範圍, 並以Level of Evidence 分別繪出各個bed檔含蓋的位點數(POS), 基因數(Gene), 藥物數(Related chamicals), 如下


# v1 2021/04/18


# read tables
data_tsv.df <- read.csv('/mnt/c/Users/vivia/project/test01/data/4/clinical_ann_metadata.tsv', sep="\t", header = T)

# list of bed.df
bed.df.list <- list()

bed.df.list[['IDT']] <- read.csv('/mnt/c/Users/vivia/project/test01/data/4/IDT_targets.bed', sep="\t", header = F)
bed.df.list[['roche']] <- read.csv('/mnt/c/Users/vivia/project/test01/data/4/ROCHE_KAPA_targets.bed', sep="\t", header = F)
bed.df.list[['V6']] <- read.csv('/mnt/c/Users/vivia/project/test01/data/4/Agilent_V6_targets.bed', sep="\t", header = F)
bed.df.list[['V7']] <- read.csv('/mnt/c/Users/vivia/project/test01/data/4/Agilent_V7_targets.bed', sep="\t", header = F)


# list of ranges for each bed
bed.range.list <- list()

for (bed in names(bed.df.list))
{
	bed.df <- bed.df.list[[bed]]
	regions <- c()
	for (row in rownames(bed.df))
	{
		regions <- c(regions, bed.df[row, 2]:bed.df[row, 3])
	}
	bed.range.list[[bed]] <- regions
	print(regions)
}


# summary for each bed
bed.summary.array <- array(character(), dim=c(6, 3, 4), dimnames = list(c('1A', '1B', '2A', '2B', '3', '4'), c('Chemicals', 'Gene', 'pos'), names(bed.df.list)))

bed.summary.array['1A', 'Chemicals', 'IDT'] <- c('chh')
bed.summary.array['1A', 'Chemicals', 'IDT'] <- c(bed.summary.array['1A', 'Chemicals', 'IDT'], 'ooo')


# read tsv data
for (row in rownames(data_tsv.df))
{
	# check if location is rsID
	if ( substr(data_tsv.df[row, 'Location'], start = 1, stop = 2) == 'rs' ) {
		position <- gsub('rs([0-9]+)', '\\1', data_tsv.df[row, 'Location'])
		# check if position is in bed
		for (bed in names(bed.range.list))
		{
			if (is.element(position, bed.range.list[[bed]])) {
				level <- data_tsv.df[row, 'Level of Evidence']
				bed.summary.list[[bed]][level, 'Chemicals'] <- paste(bed.summary.list[[bed]][level, 'Chemicals'], data_tsv.df[row, 'Related.Chemicals'], sep=',')
				bed.summary.list[[bed]][level, 'Gene'] <- paste(bed.summary.list[[bed]][level, 'Gene'], data_tsv.df[row, 'Gene'], sep=',')
				bed.summary.list[[bed]][level, 'pos'] <- paste(bed.summary.list[[bed]][level, 'pos'], position, sep=',')

			}
		}
	}
}

warnings()