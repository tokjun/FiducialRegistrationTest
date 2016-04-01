library(ggplot2)

dateString <- '2016-03-31-17-11-42'

#setwd(sprintf('/home/develop/Dropbox/Experiments/Canon/FiducialTest/Test-%s', dateString))
setwd(sprintf('/Users/junichi/Dropbox/Experiments/Canon/FiducialTest/Test-%s', dateString))
data <- read.csv(sprintf('result-%s.csv', dateString))

## xstr, bystr, condstr correspond to y-axis, x-axis, and individual line in the final plot
aggregate.errors <- function(data, xstr, condstr, paramstr) {
    a <- aggregate(data[[condstr]], list(data[[condstr]]), mean)
    r <- data.frame()
    for (i in a$Group.1) {
        data2 <- data[data[[condstr]]==i,]
        m <- aggregate(data2[[paramstr]], list(data2[[xstr]]), mean)
        s <- aggregate(data2[[paramstr]], list(data2[[xstr]]), sd)
        min <- aggregate(data2[[paramstr]], list(data2[[xstr]]), min)
        max <- aggregate(data2[[paramstr]], list(data2[[xstr]]), max)
        t <- data.frame(x=m$Group.1, mean=m$x, sd=s$x, min=min$x, max=max$x)
        t$Condition <- sprintf("%s: %f", condstr, i)
        r <- rbind(r, t)
    }
    return (r)
}

plot.errors.bars.sd <- function(r, xmin, xmax, ymin, ymax, title, xlabel, ylabel, xtick=10.0, ytick=1.0, xlegendpos=0.2, ylegendpos=0.8) {
    xmargin <- (xmax-xmin)*0.05
    ebwidth <- (xmax-xmin)*0.3
    p <- ggplot(r, aes(x=x, y=mean, group=Condition))
    p <- p + geom_line(aes(linetype=Condition, colour=Condition)) + geom_point(aes(shape=Condition, colour=Condition), size=3)
    p <- p + geom_errorbar(aes(ymin=r$mean-r$sd, ymax=r$mean+r$sd, colour=Condition, linetype=Condition), width=ebwidth, position=position_dodge(0.1))
    p <- p + labs(title=title, x=xlabel, y=ylabel)
    p <- p + scale_color_manual(values=c('#0072BE','#DA5319', '#EEB220', '#7E2F8E','#77AD30', '#4DBFEF', '#A3142F'))
    p <- p + scale_x_continuous(breaks=seq(xmin,xmax,xtick), limits=c(xmin-xmargin,xmax+xmargin)) +  scale_y_continuous(breaks=seq(ymin,ymax,ytick), limits=c(ymin,ymax))
    p <- p + theme_bw() + theme(text = element_text(size=16), legend.key.size=unit(30, "pt"), legend.position=c(xlegendpos, ylegendpos))
    
    print(p)
}

plot.errors.bars <- function(r, xmin, xmax, ymin, ymax, title, xlabel, ylabel, xtick=10.0, ytick=1.0, xlegendpos=0.2, ylegendpos=0.8) {
    xmargin <- (xmax-xmin)*0.05
    ebwidth <- (xmax-xmin)*0.02
    p <- ggplot(r, aes(x=x, y=mean, group=Condition, color=Condition))
    p <- p + geom_line(aes(linetype=Condition, colour=Condition)) + geom_point(aes(shape=Condition, colour=Condition), size=3)
    p <- p + geom_errorbar(aes(ymin=r$min, ymax=r$max, colour=Condition, linetype=Condition), width=ebwidth, position=position_dodge(0.1))
    p <- p + labs(title=title, x=xlabel, y=ylabel)
    p <- p + scale_color_manual(values=c('#0072BE','#DA5319', '#EEB220', '#7E2F8E','#77AD30', '#4DBFEF', '#A3142F'))
    p <- p + scale_x_continuous(breaks=seq(xmin,xmax,xtick), limits=c(xmin-xmargin,xmax+xmargin)) +  scale_y_continuous(breaks=seq(ymin,ymax,ytick), limits=c(ymin,ymax))
    p <- p + theme_bw() + theme(text = element_text(size=16), legend.key.size=unit(30, "pt"), legend.position=c(xlegendpos, ylegendpos))
    print(p)
}

ymin = 0.0
ymax = 3.0

## Remove NaN
data <- data[(!is.nan(data$FRE))&(!is.nan(data$TRE)),]

## Convert noise level to %
data$Noise <- data$Noise*100

##  1. Noise vs Number of fiducials (slice thickness = 2.0mm)
dataNoiseNfid <- data[data$Thickness==2,]
dataNoiseNfidFRE <- aggregate.errors(dataNoiseNfid, "Noise", "nFiducials", "FRE")
dataNoiseNfidTRE <- aggregate.errors(dataNoiseNfid, "Noise", "nFiducials", "TRE")

pdf("Noise-Fiducials-FRE.pdf")
plot.errors.bars.sd(dataNoiseNfidFRE, 0, 50, ymin, ymax, "Noise vs FRE", "Noise Level (%)", "FRE (mm)", 10, 1.0)
dev.off()
pdf("Noise-Fiducials-TRE.pdf")
plot.errors.bars.sd(dataNoiseNfidTRE, 0, 50, ymin, ymax, "Noise vs TRE", "Noise Level (%)", "TRE (mm)", 10, 1.0)
dev.off()

##  2. Slice thickness vs number of fiducials (noise = 0.2)
dataThickNfid <- data[data$Noise==20,]
dataThickNfidFRE <- aggregate.errors(dataThickNfid, "Thickness", "nFiducials", "FRE")
dataThickNfidTRE <- aggregate.errors(dataThickNfid, "Thickness", "nFiducials", "TRE")

pdf("Thickness-Fiducials-FRE.pdf")
plot.errors.bars.sd(dataThickNfidFRE, 1.0, 4.0, ymin, ymax, "Slice Thickness vs FRE", "Thickness (mm)", "FRE (mm)", 1.0, 1.0)
dev.off()
pdf("Thickness-Fiducials-TRE.pdf")
plot.errors.bars.sd(dataThickNfidTRE, 1.0, 4.0, ymin, ymax, "Slice Thickness vs TRE", "Thickness (mm)", "TRE (mm)", 1.0, 1.0)
dev.off()

### 3. Noise vs Slice thickness  (number of fiducils = 8)
#dataNoiseThick <- data[data$nFiducials==8,]
#dataNoiseThickFRE <- aggregate.errors(dataNoiseThick, "Noise", "Thickness", "FRE")
#dataNoiseThickTRE <- aggregate.errors(dataNoiseThick, "Noise", "Thickness", "TRE")
#
#pdf("Noise-Thickness-FRE.pdf")
#plot.errors.bars.sd(dataNoiseThickFRE, 0, 0.5, ymin, ymax, "Noise vs FRE", "Noise Level (%)", "FRE (mm)", 0.1, 1.0)
#dev.off()
#pdf("Noise-Thickness-TRE.pdf")
#plot.errors.bars.sd(dataNoiseThickTRE, 0, 0.5, ymin, ymax, "Noise vs TRE", "Noise Level (%)", "TRE (mm)", 0.1, 1.0)
#dev.off()


## 4. Thickness-Time
ymin = 0.0
ymax = 7.0

dataThickNfidTime <- aggregate.errors(dataThickNfid, "Thickness", "nFiducials", "WallTime")
pdf("Thickness-Time-FRE.pdf")
plot.errors.bars.sd(dataThickNfidTime, 1.0, 4.0, ymin, ymax, "Slice Thickness vs Time", "Thickness (mm)", "Time (Wall Time) (s)", 1.0, 2.0, 0.8, 0.8)
dev.off()



## 5. Statistics:

# Synthesized images were processed in average XXX+/-XXX seconds/image. The overall success rate of automatic registration was XX%. The overall FRE and TRE were XXX+/-XXX mm and XXX+/-XXX for all fiducial frame configurations. FRE and TRE were greater when less fiducials were used (XX (n=5) vs XX (n=9); p=XXX) (Fig XXX). Under our experimental condition, the noise level did not impact the FRE/TRE (Figure XX A and B). However, higher noise level resulted in a lower success rate (XXX). The FRE and TRE were affected by the slice thickness especially when larger than 4mm (XX (1mm thickness) vs XX (4 mm thickness); p=XXX) (Figure XXX). The computation time was affected by the number of slices (XXX vs XXX; p=XXX), while it wasn’t by the number of fiducials (XXX). 
print(sprintf("Number of images: %d",  nrow(data)))
print(sprintf("Average process time: %f +/- %f", mean(data$WallTime), sd(data$WallTime)))
print(sprintf("Fiducial detection rate: %f %%", sum(data$Nfid)/sum(data$nFiducials) * 100))
print(sprintf("Success rate: %f %%", nrow(data[data$FRE<10.0,])/nrow(data) * 100))

print(sprintf("Overall FRE: %f +/- %f", mean(data$FRE), sd(data$FRE)))
print(sprintf("Overall TRE (@ 150 mm): %f +/- %f", mean(data$TRE), sd(data$TRE)))

dataThickNfid1 <- dataThickNfid[dataThickNfid$Thickness==1,]
dataThickNfid2 <- dataThickNfid[dataThickNfid$Thickness==2,]
dataThickNfid3 <- dataThickNfid[dataThickNfid$Thickness==3,]
dataThickNfid4 <- dataThickNfid[dataThickNfid$Thickness==4,]

testThickFRE2 <- t.test(dataThickNfid1$FRE, dataThickNfid2$FRE, alternative="two.sided")
testThickFRE3 <- t.test(dataThickNfid2$FRE, dataThickNfid3$FRE, alternative="two.sided")
testThickFRE4 <- t.test(dataThickNfid3$FRE, dataThickNfid4$FRE, alternative="two.sided")

print(sprintf("FRE for Thickness=1.0: %f +/- %f", mean(dataThickNfid1$FRE), sd(dataThickNfid1$FRE)))
print(sprintf("FRE for Thickness=2.0: %f +/- %f (p=%f)", mean(dataThickNfid2$FRE), sd(dataThickNfid2$FRE), testThickFRE2$p.value))
print(sprintf("FRE for Thickness=3.0: %f +/- %f (p=%f)", mean(dataThickNfid3$FRE), sd(dataThickNfid3$FRE), testThickFRE3$p.value))
print(sprintf("FRE for Thickness=4.0: %f +/- %f (p=%f)", mean(dataThickNfid4$FRE), sd(dataThickNfid4$FRE), testThickFRE4$p.value))

testThickTRE2 <- t.test(dataThickNfid1$TRE, dataThickNfid2$TRE, alternative="two.sided")
testThickTRE3 <- t.test(dataThickNfid2$TRE, dataThickNfid3$TRE, alternative="two.sided")
testThickTRE4 <- t.test(dataThickNfid3$TRE, dataThickNfid4$TRE, alternative="two.sided")

print(sprintf("TRE (@ 150 mm) for Thickness=1.0: %f +/- %f", mean(dataThickNfid1$TRE), sd(dataThickNfid1$TRE)))
print(sprintf("TRE (@ 150 mm) for Thickness=2.0: %f +/- %f (p=%f)", mean(dataThickNfid2$TRE), sd(dataThickNfid2$TRE), testThickTRE2$p.value))
print(sprintf("TRE (@ 150 mm) for Thickness=3.0: %f +/- %f (p=%f)", mean(dataThickNfid3$TRE), sd(dataThickNfid3$TRE), testThickTRE3$p.value))
print(sprintf("TRE (@ 150 mm) for Thickness=4.0: %f +/- %f (p=%f)", mean(dataThickNfid4$TRE), sd(dataThickNfid4$TRE), testThickTRE4$p.value))

dataNoiseNfid0 <- dataNoiseNfid[dataNoiseNfid$Noise==0,]
dataNoiseNfid1 <- dataNoiseNfid[dataNoiseNfid$Noise==10,]
dataNoiseNfid2 <- dataNoiseNfid[dataNoiseNfid$Noise==20,]
dataNoiseNfid3 <- dataNoiseNfid[dataNoiseNfid$Noise==30,]
dataNoiseNfid4 <- dataNoiseNfid[dataNoiseNfid$Noise==40,]
dataNoiseNfid5 <- dataNoiseNfid[dataNoiseNfid$Noise==50,]

testNoiseFRE1 <- t.test(dataNoiseNfid0$FRE, dataNoiseNfid1$FRE, alternative="two.sided")
testNoiseFRE2 <- t.test(dataNoiseNfid1$FRE, dataNoiseNfid2$FRE, alternative="two.sided")
testNoiseFRE3 <- t.test(dataNoiseNfid2$FRE, dataNoiseNfid3$FRE, alternative="two.sided")
testNoiseFRE4 <- t.test(dataNoiseNfid3$FRE, dataNoiseNfid4$FRE, alternative="two.sided")
testNoiseFRE5 <- t.test(dataNoiseNfid4$FRE, dataNoiseNfid5$FRE, alternative="two.sided")

print(sprintf("FRE for Noise= 0%%: %f +/- %f", mean(dataNoiseNfid0$FRE), sd(dataNoiseNfid0$FRE)))
print(sprintf("FRE for Noise=10%%: %f +/- %f (p=%f)", mean(dataNoiseNfid1$FRE), sd(dataNoiseNfid1$FRE), testNoiseFRE1$p.value))
print(sprintf("FRE for Noise=20%%: %f +/- %f (p=%f)", mean(dataNoiseNfid2$FRE), sd(dataNoiseNfid2$FRE), testNoiseFRE2$p.value))
print(sprintf("FRE for Noise=30%%: %f +/- %f (p=%f)", mean(dataNoiseNfid3$FRE), sd(dataNoiseNfid3$FRE), testNoiseFRE3$p.value))
print(sprintf("FRE for Noise=40%%: %f +/- %f (p=%f)", mean(dataNoiseNfid4$FRE), sd(dataNoiseNfid4$FRE), testNoiseFRE4$p.value))
print(sprintf("FRE for Noise=50%%: %f +/- %f (p=%f)", mean(dataNoiseNfid5$FRE), sd(dataNoiseNfid5$FRE), testNoiseFRE5$p.value))


testNoiseTRE1 <- t.test(dataNoiseNfid0$TRE, dataNoiseNfid1$TRE, alternative="two.sided")
testNoiseTRE2 <- t.test(dataNoiseNfid1$TRE, dataNoiseNfid2$TRE, alternative="two.sided")
testNoiseTRE3 <- t.test(dataNoiseNfid2$TRE, dataNoiseNfid3$TRE, alternative="two.sided")
testNoiseTRE4 <- t.test(dataNoiseNfid3$TRE, dataNoiseNfid4$TRE, alternative="two.sided")
testNoiseTRE5 <- t.test(dataNoiseNfid4$TRE, dataNoiseNfid5$TRE, alternative="two.sided")

print(sprintf("TRE for Noise= 0%%: %f +/- %f", mean(dataNoiseNfid0$TRE), sd(dataNoiseNfid0$TRE)))
print(sprintf("TRE for Noise=10%%: %f +/- %f (p=%f)", mean(dataNoiseNfid1$TRE), sd(dataNoiseNfid1$TRE), testNoiseTRE1$p.value))
print(sprintf("TRE for Noise=20%%: %f +/- %f (p=%f)", mean(dataNoiseNfid2$TRE), sd(dataNoiseNfid2$TRE), testNoiseTRE2$p.value))
print(sprintf("TRE for Noise=30%%: %f +/- %f (p=%f)", mean(dataNoiseNfid3$TRE), sd(dataNoiseNfid3$TRE), testNoiseTRE3$p.value))
print(sprintf("TRE for Noise=40%%: %f +/- %f (p=%f)", mean(dataNoiseNfid4$TRE), sd(dataNoiseNfid4$TRE), testNoiseTRE4$p.value))
print(sprintf("TRE for Noise=50%%: %f +/- %f (p=%f)", mean(dataNoiseNfid5$TRE), sd(dataNoiseNfid5$TRE), testNoiseTRE5$p.value))




dataNfid5 <- data[data$nFiducials==5,]
dataNfid6 <- data[data$nFiducials==6,]
dataNfid7 <- data[data$nFiducials==7,]
dataNfid8 <- data[data$nFiducials==8,]
dataNfid9 <- data[data$nFiducials==9,]

testNfidFRE6 <- t.test(dataNfid5$FRE, dataNfid6$FRE, alternative="two.sided")
testNfidFRE7 <- t.test(dataNfid5$FRE, dataNfid7$FRE, alternative="two.sided")
testNfidFRE8 <- t.test(dataNfid5$FRE, dataNfid8$FRE, alternative="two.sided")
testNfidFRE9 <- t.test(dataNfid5$FRE, dataNfid9$FRE, alternative="two.sided")

print(sprintf("FRE for nFiducials=5: %f +/- %f", mean(dataNfid5$FRE), sd(dataNfid5$FRE)))
print(sprintf("FRE for nFiducials=6: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid6$FRE), sd(dataNfid6$FRE), testNfidFRE6$p.value))
print(sprintf("FRE for nFiducials=7: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid7$FRE), sd(dataNfid7$FRE), testNfidFRE7$p.value))
print(sprintf("FRE for nFiducials=8: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid8$FRE), sd(dataNfid8$FRE), testNfidFRE8$p.value))
print(sprintf("FRE for nFiducials=9: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid9$FRE), sd(dataNfid9$FRE), testNfidFRE9$p.value))

testNfidTRE6 <- t.test(dataNfid5$TRE, dataNfid6$TRE, alternative="two.sided")
testNfidTRE7 <- t.test(dataNfid5$TRE, dataNfid7$TRE, alternative="two.sided")
testNfidTRE8 <- t.test(dataNfid5$TRE, dataNfid8$TRE, alternative="two.sided")
testNfidTRE9 <- t.test(dataNfid5$TRE, dataNfid9$TRE, alternative="two.sided")

print(sprintf("TRE for nFiducials=5: %f +/- %f", mean(dataNfid5$TRE), sd(dataNfid5$TRE)))
print(sprintf("TRE for nFiducials=6: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid6$TRE), sd(dataNfid6$TRE), testNfidTRE6$p.value))
print(sprintf("TRE for nFiducials=7: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid7$TRE), sd(dataNfid7$TRE), testNfidTRE7$p.value))
print(sprintf("TRE for nFiducials=8: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid8$TRE), sd(dataNfid8$TRE), testNfidTRE8$p.value))
print(sprintf("TRE for nFiducials=9: %f +/- %f (p=%f vs nFiducials=5)", mean(dataNfid9$TRE), sd(dataNfid9$TRE), testNfidTRE9$p.value))

