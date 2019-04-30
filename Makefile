figure-pepseq-labels/index.html: figure-pepseq-labels/index.Rmd
	cd figure-pepseq-labels && R -e 'rmarkdown::render("index.Rmd")'
pepseq9.labels.csv: pepseq9.labels.R
	R --vanilla < $<
figure-pepseq-example.png: figure-pepseq-example.R
	R --vanilla < $<
