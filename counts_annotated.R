source("packages.R")

counts <- fread("counts_annotated.csv")
not.na <- counts[!is.na(`Accession (UniProt)`)]
setkey(not.na, `Accession (UniProt)`, `Peptide start`)
ends <- not.na[, list(
  first.end=min(`Peptide end`),
  last.end=max(`Peptide end`)
), by=list(`Accession (UniProt)`)]

fwrite(
  ends[order(last.end), list(
    `Accession (UniProt)`, last.end, `Accession (UniProt)`)],
  "pepseq9.txt",
  sep="\t",
  col.names=FALSE)
system("gzip pepseq9.txt")
not.na[, chromEnd := `Peptide end`]
not.na[, chromStart := chromEnd-1]

header.tmp <- paste(
  'track',
  'type=bedGraph',
  'db=pepseq9',
  'export=yes',
  'visibility=full',
  'maxSegments=20',
  'alwaysZero=on',
  'share=public',
  'graphType=points',
  'yLineMark=0',
  'yLineOnOff=on',
  'name=%s',
  'description="%s count"')
count.i.vec <- 3:41
for(col.i in count.i.vec){
  col.name <- names(not.na)[col.i]
  name.vec <- c("Accession (UniProt)", "chromStart", "chromEnd", col.name)
  header <- sprintf(header.tmp, gsub("[._]", "-", col.name), col.name)
  out <- not.na[, ..name.vec]
  dir.create("bedGraph", showWarnings=FALSE)
  out.bedGraph <- file.path("bedGraph", paste0(col.name, ".bedGraph.gz"))
  con <- gzfile(out.bedGraph, "w")
  writeLines(header, con)
  write.table(out, con, quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(con)
}

molt <- melt(not.na, measure.vars=3:41)
molt[, list(
  max.count=max(value, na.rm=TRUE)
), by=list(variable)][order(max.count)]


molt <- melt(not.na, measure.vars=3:11)
library(ggplot2)
molt[, acc := sub("_", "\n", `Accession (UniProt)`)]
molt[, var := sub("_", "\n", variable)]
ggplot()+
  theme_bw()+
  facet_grid(var ~ acc, scales="free", space="free")+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    chromEnd, value),
    shape=1,
    data=molt)+
  xlab("Peptide end")+
  ylab("Read counts")

ggplot()+
  theme_bw()+
  facet_grid(var ~ acc, scales="free", space="free")+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    chromEnd, value),
    shape=1,
    data=molt[acc=="Q9WLS3"])+
  xlab("Peptide end")+
  ylab("Read counts")
