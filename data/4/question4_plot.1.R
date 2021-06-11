#!/usr/bin/env Rscript

# question4_plot.1.R

# Shu-Ting Cho <shutingcho@pitt.edu>
# for plotting result from question4.4.py

# v1 2021/04/18


data <- read.csv('./out_data.tsv', sep="\t", header = T)

# Grouped and facet
library(ggplot2)
p <- ggplot(data, aes(fill=kit, y=Counts, x=Evidence.level)) + 
	geom_bar(position="dodge", stat="identity") + 
	scale_x_discrete(name="Evidence level") + 
	geom_text(aes(label=Counts), position=position_dodge(width=0.9), hjust=-0.1, size=3, inherit.aes = TRUE) + 
	expand_limits(y=c(0, 750)) + 
	coord_flip() + 
	facet_wrap(~info) + 
	theme_classic()

# write output pdf
pdf(file = './Rplot_bar.pdf', width=8, height=5)
p
dev.off() 