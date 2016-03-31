library(ggplot2)

dateString <- '2016-03-30-18-10-20'

setwd(sprintf('/home/develop/Dropbox/Experiments/Canon/FiducialTest/Test-%s', dateString))
data <- read.csv(sprintf('result-%s.csv', dateString))

## xstr, bystr, condstr correspond to y-axis, x-axis, and individual line in the final plot
aggregate.errors <- function(data, xstr, condstr, paramstr) {
    a <- aggregate(data[[condstr]], list(data[[condstr]]), mean)
    r <- data.frame()
    for (i in a$Group.1) {
        data2 <- data[data[[condstr]]==i,]
        m <- aggregate(data2[[paramstr]], list(data2[[xstr]]), mean)
        s <- aggregate(data2[[paramstr]], list(data2[[xstr]]), sd)
        t <- data.frame(x=m$Group.1, mean=m$x, sd=s$x)
        t$Condition <- sprintf("%s: %f", condstr, i)
        r <- rbind(r, t)
    }
    return (r)
}

plot.errors.bars <- function(r, xmin, xmax, ymin, ymax, title, xlabel, ylabel, xtick=10.0, ytick=1.0) {
    xmargin <- (xmax-xmin)*0.05
    ebwidth <- (xmax-xmin)*0.05
    p <- ggplot(r, aes(x=x, y=mean, group=Condition, color=Condition)) + geom_line() + geom_point() + geom_errorbar(aes(ymin=r$mean-r$sd, ymax=r$mean+r$sd), width=ebwidth, position=position_dodge(0.5)) + labs(title=title, x=xlabel, y=ylabel) + theme_classic() + scale_color_manual(values=c('#0072BE','#DA5319', '#EEB220', '#7E2F8E','#77AD30', '#4DBFEF', '#A3142F')) + scale_x_continuous(breaks=seq(xmin,xmax,xtick), limits=c(xmin-xmargin,xmax+xmargin)) +  scale_y_continuous(breaks=seq(ymin,ymax,ytick), limits=c(ymin,ymax))
    print(p)
}


ymin = 0.0
ymax = 10.0

## Remove NaN
data <- data[(!is.nan(data$FRE))&(!is.nan(data$TRE)),]

## Three comparisons:
##  1. Noise vs Number of fiducials (slice thickness = 2.0mm)
dataNoiseNFid <- data[data$Thickness==2,]
dataNoiseNfidFRE <- aggregate.errors(dataNoiseNFid, "Noise", "nFiducials", "FRE")
dataNoiseNfidTRE <- aggregate.errors(dataNoiseNFid, "Noise", "nFiducials", "TRE")

pdf("Noise-Fiducials-FRE.pdf")
plot.errors.bars(dataNoiseNfidFRE, 0, 0.5, ymin, ymax, "Noise vs FRE", "Signal-to-Noise Ratio", "FRE (mm)", 0.1, 1.0)
dev.off()
pdf("Noise-Fiducials-TRE.pdf")
plot.errors.bars(dataNoiseNfidTRE, 0, 0.5, ymin, ymax, "Noise vs TRE", "Signal-to-Noise Ratio", "TRE (mm)", 0.1, 1.0)
dev.off()

##  2. Slice thickness vs number of fiducials (noise = 0.2)
dataThickNfid <- data[data$Noise==0.2,]
dataThickNfidFRE <- aggregate.errors(dataThickNfid, "Thickness", "nFiducials", "FRE")
dataThickNfidTRE <- aggregate.errors(dataThickNfid, "Thickness", "nFiducials", "TRE")

pdf("Thickness-Fiducials-FRE.pdf")
plot.errors.bars(dataThickNfidFRE, 0, 5.0, ymin, ymax, "Slice Thickness vs FRE", "Thickness (mm)", "FRE (mm)", 1.0, 1.0)
dev.off()
pdf("Thickness-Fiducials-TRE.pdf")
plot.errors.bars(dataThickNfidTRE, 0, 5.0, ymin, ymax, "Slice Thickness vs TRE", "Thickness (mm)", "TRE (mm)", 1.0, 1.0)
dev.off()

##  3. Noise vs Slice thickness  (number of fiducils = 8)
dataNoiseThick <- data[data$nFiducials==8,]
dataNoiseThickFRE <- aggregate.errors(dataNoiseThick, "Noise", "Thickness", "FRE")
dataNoiseThickTRE <- aggregate.errors(dataNoiseThick, "Noise", "Thickness", "TRE")

pdf("Noise-Thickness-FRE.pdf")
plot.errors.bars(dataNoiseThickFRE, 0, 0.5, ymin, ymax, "Noise vs FRE", "Signal-to-Noise Ratio", "FRE (mm)", 0.1, 1.0)
dev.off()
pdf("Noise-Thickness-TRE.pdf")
plot.errors.bars(dataNoiseThickTRE, 0, 0.5, ymin, ymax, "Noise vs TRE", "Signal-to-Noise Ratio", "TRE (mm)", 0.1, 1.0)
dev.off()

