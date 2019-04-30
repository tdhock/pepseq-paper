pepseq9.labels.csv: pepseq9.labels.R
	R --vanilla < $<
figure-pepseq-example.png: figure-pepseq-example.R
	R --vanilla < $<
