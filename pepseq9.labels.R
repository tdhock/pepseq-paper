source("packages.R")

server <- "http://127.0.0.1:8080/"
download.fread <- function(u, f){
  if(!file.exists(f)){
    dir.create(dirname(f), showWarnings=FALSE)
    download.file(u, f)
  }
  fread(f)
}
dest.dir <- "pepseq9.labels"
profiles <- download.fread(
  paste0(server, "csv_profiles/"),
  "pepseq9.profile.list.csv")

match.dt <- namedCapture::df_match_variable(
  profiles, description=list(
    fullname=list(
      drb="[^.]+",
      "[.]",
      id="[^_ ]+"),
    cleaved=list(
      "[^ ]*"
    ), "?"))
pair.dt <- match.dt[description.cleaved %in% c("", "_cleaved")]
dcast(
  pair.dt,
  description.drb + description.id ~ description.cleaved)

uncleaved.dt <- pair.dt[description.cleaved == ""]

over.list <- list()
for(pair.i in 1:nrow(uncleaved.dt)){
  uncleaved.row <- uncleaved.dt[pair.i]
  name.list <- with(uncleaved.row, list(
    uncleaved=name,
    cleaved=paste0(name, "-cleaved")))
  dt.list <- list()
  for(sample.type in names(name.list)){
    sample.name <- name.list[[sample.type]]
    dt <- download.fread(
      paste0(server, "export/None/", sample.name, "/regions/csv/"),
      file.path(dest.dir, paste0(sample.name, ".csv")))
    new.list <- list(protein=dt$chromosome)
    for(col.name in c("min", "max", "annotation")){
      new.list[[paste0(sample.type, "_", col.name)]] <- dt[[col.name]]
    }
    new.dt <- do.call(data.table, new.list)
    setkeyv(
      new.dt,
      grep("annotation", names(new.list), invert=TRUE, value=TRUE))
    dt.list[[sample.type]] <- new.dt
  }
  join.dt <- foverlaps(dt.list$cleaved, dt.list$uncleaved, nomatch=0L)
  over.list[[pair.i]] <- data.table(
    uncleaved.name=uncleaved.row$description.fullname,
    join.dt)
}
over.dt <- do.call(rbind, over.list)

over.dt[, table(uncleaved_annotation)]
over.dt[uncleaved_annotation=="1breakpoint"]

(count.dt <- over.dt[, .(
  labels=.N
), by=list(
  uncleaved.name, protein, uncleaved_min, uncleaved_max, cleaved_annotation)])
print(count.dt[, table(labels, cleaved_annotation)])
bad <- count.dt[labels==1 & cleaved_annotation=="1breakpoint"]
if(nrow(bad)){
  print(bad)
  stop("only one 1breakpoint label on cleaved")
}

fwrite(over.dt, "pepseq9.labels.tall.csv")

(pepseq9.labels <- over.dt[, {
  label <- if(cleaved_annotation=="1breakpoint"){
    c("peakStart", "peakEnd")
  }else{
    "noPeaks"
  }
  data.table(
    label,
    labelStart=cleaved_min, labelEnd=cleaved_max)
}, by=list(
  uncleaved.name, protein, uncleaved_min, uncleaved_max, cleaved_annotation)])

fwrite(pepseq9.labels, "pepseq9.labels.csv")
