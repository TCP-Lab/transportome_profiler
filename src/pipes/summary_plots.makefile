#? Generate summary plots and other bits and bobs
#? 
#? This takes the output of running the heatmaps pipeline over and over
#? with all the different ranking methods and generates the final plots

.ONESHELL:

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@


