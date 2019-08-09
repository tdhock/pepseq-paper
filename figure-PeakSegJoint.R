library(namedCapture)
library(data.table)
library(animint2)

gz.vec <- Sys.glob("bedGraph/*.gz")
count.dt.list <- list()
for(file.i in seq_along(gz.vec)){
  gz <- gz.vec[[file.i]]
  one <- fread(cmd=paste("zcat", gz, "| grep O56266"), skip=1, col.names=c(
    "protein", "chromStart", "chromEnd", "count"))
  hla <- sub(".*/", "", sub(".bedGraph.gz", "", gz))
  count.dt.list[[file.i]] <- data.table(
    hla, one)
}
count.dt <- do.call(rbind, count.dt.list)

some <- namedCapture::df_match_variable(
  count.dt[grepl("DRB1.(14|08).0[13]", hla)],
  hla=list(
    nomatch.error=TRUE,
    fullname=list(
      drb="[^.]+",
      "[.]",
      id="[^_ ]+"),
    cleaved=list(
      "[^ ]*"
    ), "?"))[250 < chromEnd & chromEnd < 300]
some[, hla.Cleaved := ifelse(hla.cleaved=="_cleaved", "cleaved", "un-cleaved")]
some[, hla.panel := paste(hla.id, ifelse(hla.cleaved=="_cleaved", "cl", "un"))]
some[, sample.id := hla]
fit <- PeakSegJoint::PeakSegJointSeveral(some)
info <- PeakSegJoint::ConvertModelList(fit)

id.dt <- unique(some[, grep("hla", names(some), value=TRUE), with=FALSE])
seg.dt <- data.table(info$segments)[id.dt, on=list(sample.id=hla)]
big.models <- seg.dt[peaks==max(peaks)]
small.models <- seg.dt[peaks==0]
big.models[small.models, bkg.mean := i.mean, on=.(sample.id)]
status.dt <- seg.dt[, .(
  status=if(.N==3)"peak" else "bkg"
), by=.(sample.id, peaks)]
status.big <- status.dt[big.models, on=.(sample.id), allow.cartesian=TRUE]
status.big[, show.mean := ifelse(status=="peak", mean, bkg.mean)]

viz <- animint(
  out.dir="figure-PeakSegJoint",
  loss=ggplot()+
    ggtitle("Select joint peak model")+
    geom_point(aes(
      peaks, loss),
      data=info$loss)+
    geom_text(aes(
      max(peaks), max(loss[loss<Inf]),
      key=1,
      label=sprintf(
        "loss=%.2f peaks=%d",
        loss, peaks)),
      hjust=1,
      showSelected="peaks",
      data=info$loss)+
    theme_bw()+
    make_tallrect(info$loss, "peaks"),
  data=ggplot()+
    ggtitle("Selected model for DRB1, O56266 data")+
    geom_point(aes(
      chromEnd, count),
      data=some)+
    facet_grid(hla.panel ~ ., scales="free")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    xlab("Position in protein (amino acids)")+
    ylab("Count of aligned reads in cl=cleaved/un=uncleaved")+
    geom_segment(aes(
      chromStart+0.5, show.mean,
      key=chromStart,
      xend=chromEnd+0.5, yend=show.mean),
      color="green",
      showSelected="peaks",
      data=status.big),
  duration=list(peaks=1000),
  time=list(ms=2000, variable="peaks")
)
##animint2gist(viz)
