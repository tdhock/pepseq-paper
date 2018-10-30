source("packages.R")

pepseq <- fread("figure-pepseq-example.csv")

signal.name <- "HLA-DRB1*07:01:signal"
y.lab <- paste(signal.name)
pepseq$count <- as.integer(pepseq[[signal.name]])
pepseq[, chromStart := Position-1L]
pepseq[, chromEnd := Position]
proteins <- pepseq[, list(
  positions=max(Position),
  n.NA=sum(is.na(count))
  ), by=list(UniProt)][order(positions)]

some.pepseq <- pepseq[UniProt %in% proteins[50:60, UniProt] & !is.na(count)]
some.gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(UniProt ~ ., scales="free")+
  geom_point(aes(
    Position, count),
    shape=1,
    data=some.pepseq)

## Run peak detection on one sample up to 20 peaks.
one.pepseq <- pepseq[UniProt=="P01009" & !is.na(count)]
fit <- PeakSegPDPAchrom(one.pepseq, 20L)
max.feasible.peaks <- data.table(fit$loss)[feasible==TRUE, max(peaks)]
show.segments <- data.table(fit$segments)[peaks <= max.feasible.peaks+2]
show.changes <- show.segments[, data.table(
  position=chromStart[-1]+0.5,
  diff=diff(mean)
  ), by=list(peaks)]
show.changes[, constraint := ifelse(diff==0, "equality", "inequality")]

## Plot segment means of models from no peaks to max feasible peaks.
gg <- ggplot()+
  ggtitle("Equality constraints start at 11 peaks => Most complex model which still has alternating up and then down changes = 10 peaks")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(peaks ~ ., labeller=label_both)+
  ylab(y.lab)+
  geom_point(aes(
    Position, count),
    shape=1,
    data=one.pepseq)+
  geom_segment(aes(
    chromStart+0.5, mean,
    xend=chromEnd+0.5, yend=mean),
    color="green",
    data=show.segments)+
  scale_linetype_manual(values=c(equality="solid", inequality="dotted"))+
  geom_vline(aes(
    xintercept=position, linetype=constraint),
    color="green",
    data=show.changes)
png("figure-pepseq-example-mean.png", 12, 8, res=100, units="in")
print(gg)
dev.off()

## Plot peak regions of models from no peaks to max feasible peaks.
peak.y <- -25
gg <- ggplot()+
  ggtitle("Blue segments show peak regions (even-numbered segments); blue points show peak start position")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(peaks ~ ., labeller=label_both)+
  ylab(y.lab)+
  geom_point(aes(
    Position, count),
    shape=1,
    data=one.pepseq)+
  geom_segment(aes(
    chromStart+0.5, peak.y,
    xend=chromEnd+0.5, yend=peak.y),
    color="deepskyblue",
    size=1.5,
    data=show.segments[status=="peak"])+
  geom_point(aes(
    chromStart+0.5, peak.y),
    color="deepskyblue",
    shape=1,
    data=show.segments[status=="peak"])
png("figure-pepseq-example-peaks.png", 12, 8, res=100, units="in")
print(gg)
dev.off()

some.peaks <- some.pepseq[, {
  fit <- PeakSegPDPAchrom(.SD, 50L)
  max.feasible.peaks <- data.table(fit$loss)[feasible==TRUE, max(peaks)]
  data.table(fit$segments)[peaks==max.feasible.peaks & status=="peak"]
}, by=list(UniProt)]
max.dt <- some.pepseq[, list(
  max.count=max(count)
), by=list(UniProt)]
max.peaks <- some.peaks[max.dt, on=list(UniProt)]
(some.text <- some.peaks[, .SD[1], by=list(UniProt)])
max.peaks[, peak.y := -0.1*max.count]

gg <- some.gg+
  ggtitle("Most complex up-down model for several proteins")+
  geom_segment(aes(
    chromStart+0.5, peak.y,
    xend=chromEnd+0.5, yend=peak.y),
    color="deepskyblue",
    size=1.5,
    data=max.peaks)+
  geom_point(aes(
    chromStart+0.5, peak.y),
    color="deepskyblue",
    shape=1,
    data=max.peaks)

png("figure-pepseq-example.png", 12, 8, res=100, units="in")
print(gg)
dev.off()
 
