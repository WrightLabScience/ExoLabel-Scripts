SOURCE_DIR <- '~'

#########################
### SYNTHETIC RESULTS ###
#########################

## Note: Missing a lot of data here
THRESHOLD <- 4e4
in_file <- file.path(SOURCE_DIR,'SyntheticResults.RData')
load(in_file)
outfile <- basename(in_file)
outfile_supp <- gsub(".RData", "_Suppl.pdf", outfile)
outfile_main <- gsub(".RData", "_Main.pdf", outfile)
outfile_supp <- gsub("Data", "Results", outfile_supp)
outfile_main <- gsub("Data", "Results", outfile_main)
outfile_supp <- file.path(SOURCE_DIR, outfile_supp)
outfile_main <- file.path(SOURCE_DIR, outfile_main)

names(AllResults)
categories <- list(
  LP=c("LabelProp", "Speakeasy"),
  Modularity=c("FastGreedy", "Leiden", "Louvain", "LeadingEigen"),
  RandomWalk=c("Walktrap", "Infomap"),
  MCL=c("MCL1.4", "MCL2.0", "MCL4.0"),
  Other=c("Spinglass"),
  ExoLabel=c("ExoLabel00", "ExoLabel05", "ExoLabel10")
)
DATA_SKIP <- c("Spinglass")
CHANGE_SKIP <- c("Spinglass", "Speakeasy", "ExoLabel05", "ExoLabel10")
cols <- list(
  LP=c('#E0A608','#705304'), # label propagation
  Modularity=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
  RandomWalk=c('#000000', '#606060'),
  MCL=c('#E82B70','#C80B50','#A80030'),
  Other=c('#45A649'),
  ExoLabel=c("#824484", "#824484", "#824484")
)
pchs <- c(16,17,15,4)
lwd_val <- 0.5

AllResults <- AllResults[match(unlist(categories), names(AllResults))]
p <- !(names(AllResults) %in% CHANGE_SKIP)

#### SUPPLEMENTAL FIGURE ####
pdf(outfile_supp, height=4.3*1.5, width=4.3*1.5, onefile=TRUE)
par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
layout(matrix(c(1:6,8,8,7), nrow=3, byrow=TRUE))

cols <- unlist(lapply(seq_along(categories), \(i){
  cols[[i]][seq_along(categories[[i]])]
}))
pchs <- unlist(lapply(seq_along(categories), \(i){
  pchs[seq_along(categories[[i]])]
}))
#check_col <- "LabelProp"
check_col <- "ExoLabel00"
pctmasking <- colnames(AllResults[[1]])
nnodes <- GridParams[,"num_nodes"]

# data gets very wonky above this threshold
nnodes[nnodes > THRESHOLD] <- NA

nodecounts <- table(nnodes)
for(i in seq_len(6)){
  check <- GridParams[,"mu_topology"] >= GridParams[,"mu_weight"] & !is.na(AllResults[[check_col]][,i])
  plot(NULL, xlim=range(nnodes, na.rm = TRUE), ylim=c(0,1), main=paste0(pctmasking[i], " Masking"),
       ylab=ifelse(i%%3 == 1, "Mean AMI", ''),
       xlab=ifelse(i>3, "Number of Nodes", ''), log='x')

  for(j in seq_along(AllResults)){
    if(!p[j]) next
    v <- AllResults[[j]][,i]
    pts <- tapply(v[check], nnodes[check], mean, na.rm=TRUE)
    #pts <- tapply(v, nnodes, mean, na.rm=TRUE)
    prop_there <- tapply(v[check], nnodes[check], \(x) sum(!is.na(x)) / length(x))
    pts <- pts[prop_there >= 0.5]
    lines(y=pts, x=as.numeric(names(pts)), col=cols[j], lwd=lwd_val)
    points(y=pts, x=as.numeric(names(pts)), col=cols[j], pch=pchs[j], cex=0.5)
  }
}

## plot mean change vs. baseline
mean_change <- lapply(AllResults,
                      \(x){
                        v <- rep(0, 5)
                        for(i in seq_len(5)){
                          p <- !is.na(x[,1]) & x[,1] >= 0.01 & x[,i] >= 0.01
                          tmp <- (x[p,i+1] - x[p,1]) / x[p,1]
                          tmp[is.infinite(tmp)] <- 1
                          if(sum(is.na(tmp)) > 0.5 * length(tmp))
                            v[i] <- NA
                          else
                            v[i] <- weighted.mean(tmp, GridParams[p,3], na.rm=TRUE)*100
                            #v[i] <- mean(tmp, na.rm=TRUE)*100
                        }
                        v
                      })
#plot(0, axes=FALSE, ann=FALSE, frame.plot=FALSE, cex=0)
temp_n <- names(AllResults)
temp_n[which(temp_n == "ExoLabel00")] <- "ExoLabel"
names(AllResults) <- temp_n
p <- which(p)
plot(0, axes=FALSE, ann=FALSE, frame.plot=FALSE, cex=0)
legend("center", legend=names(AllResults)[p], lty=1,
       col=cols[p], pch=pchs[p], cex=1, bty='n', ncol=2, xpd=NA)
xlm <- c(0,50)
ylm <- c(-60,20)
plot(NULL, xlim=xlm, ylim=ylm, #ylim=range(unlist(mean_change)),
     main="", yaxt='n', xaxt='n',
     ylab="Weighted % Change in AMI", xlab="Percent Edges Masked")
axis(1, at=seq(xlm[1],xlm[2],10), labels = paste0(seq(xlm[1],xlm[2],10), '%'))
axis(2, at=seq(ylm[1],ylm[2],20), labels = paste0(seq(ylm[1],ylm[2],20), '%'))
for(i in seq_along(mean_change)){
  if(names(mean_change)[i] %in% CHANGE_SKIP) next
  lines(x=seq(0,50,10), y=c(0, mean_change[[i]]), col=cols[i], lwd=lwd_val*2, type='o', pch=pchs[i])
}
abline(h=0, lty=2)

dev.off(dev.list()['pdf'])

#### MAIN FIGURE ####

pdf(outfile_main, height=4.3, width=4.3*2, onefile=TRUE)
par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
#layout(matrix(c(1,3,2), nrow=1, byrow=TRUE), widths=c(1,1,0.33))
bench_col <- "ExoLabel00"
to_remove <- c("ExoLabel00", "ExoLabel05", "ExoLabel10", "Speakeasy")
pt_cex <- 0.75
HEAD_LEN <- 0.05
TXT_YOFF <- -0.8

p_use <- GridParams[,"mu_topology"] >= GridParams[,"mu_weight"] & !is.na(GridParams[,"num_nodes"])
all_amis <- vapply(AllResults, \(x) apply(AllResults[[bench_col]][p_use,] - x[p_use,], 2,
                                          mean,#weighted.mean,
                                          #median,
                                          #w=GridParams[p_use,"num_nodes"],
                                          na.rm=TRUE), numeric(6L))
all_amis_med <- vapply(AllResults, \(x) apply(AllResults[[bench_col]][p_use,] - x[p_use,], 2,
                                          median,
                                          na.rm=TRUE), numeric(6L))
all_amis_sds <- vapply(AllResults,\(x) apply(AllResults[[bench_col]][p_use,] - x[p_use,], 2,
                                             \(x)quantile(x,c(0.25,0.75),na.rm=TRUE)), numeric(12L))
to_remove <- match(to_remove, colnames(all_amis))
all_amis <- all_amis[,-to_remove]
all_amis_sds <- all_amis_sds[,-to_remove]
all_amis_med <- all_amis_med[,-to_remove]

o <- order(colMeans(all_amis))
nbars <- nrow(all_amis) * ncol(all_amis)
nc <- ncol(all_amis)
barcols <- cols[-to_remove][o]
barplot(c(all_amis[,o]),
        col=rep(paste0(cols[-to_remove][o], "C0"), each=6),
        space=c(1, rep(0,5)),
        ylim=c(-0.05,0.8), ylab="Median AMI Decrease from ExoLabel")
abline(h=0, lty=1)
# plot(NULL, xaxt='n', ylim=c(-1,1), xlim=c(0, nbars), xlab='',
#      ylab="AMI Decrease from ExoLabel", frame.plot=FALSE)
for(i in seq_len(ncol(all_amis))){
  cur_col <- colnames(all_amis)[o][i]
  for(j in seq_len(6)){
    ## Violin plots
    # r1 <- AllResults[[bench_col]][p_use,j]
    # r2 <- AllResults[[cur_col]][p_use,j]
    # #p <- !is.na(r1) & !is.na(r2)
    # #r1 <- r1[p]
    # #r2 <- r2[p]
    # #r1[is.na(r1)] <- 0
    # #r2[is.na(r2)] <- 0
    # r <- r1 - r2
    # r <- r[!is.na(r)]
    #
    # #den <- sm.density(r, 0.015, display='none', verbose=0)
    # # den <- sm.density(r, 0.015, model="Normal", display='none', verbose=0)
    # # yv <- den$eval.points
    # # xv_den <- den$estimate
    # den <- density(r)
    # xv_den <- den$y
    # yv <- den$x
    # xv_den <- xv_den / max(xv_den) * 0.5
    # crit_points <- quantile(r, c(0.025,0.5,0.975), na.rm=TRUE)
    #
    # xp <- 7*(i-1) + j - 1.5
    # polygon(y=c(yv,rev(yv)), x=c(xp+xv_den, xp-rev(xv_den)),
    #         col=barcols[i], border='black', lwd=0.5)
    # lines(y=rep(crit_points[2],2), x=xp+c(-0.25,0.25), col='black', lwd=2)
    #yv <- quantile(r, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
    #rect(xleft=xp[1], xright=xp[2], ybottom=yv[1], ytop=yv[2],
    #     col=barcols[i], border = 'black')
    #rect(xleft=xp[1], xright=xp[2], ybottom=yv[2], ytop=yv[3],
    #     col=barcols[i], border = 'black')

    ## confidence intervals
    lines(x=rep(7*(i-1)+0.5 + j, 2),
          #y=all_amis[,o][j,i] + (all_amis_sds[,o][2*j + c(-1,0),i]))
          y=all_amis_sds[,o][2*j + c(-1,0),i],
          lwd=2)
    points(x=7*(i-1)+0.5 + j, y=all_amis_med[,o][j,i])
  }
}
text(x=seq_len(ncol(all_amis))*7 - 3, y=0, colnames(all_amis)[o],
     xpd=NA, srt=0, adj=c(0.5,4), cex=0.75)
bar_width <- 5
start_pos <- 1 #0.25*nbars
ybar_start <- 0.45
ybar_inc <- 0.05
arr_y <- ybar_start - ybar_inc/2
rect(xleft=seq(0,5*bar_width,bar_width) + start_pos,
     xright=seq(0,5*bar_width,bar_width) + start_pos+bar_width,
     ybottom = ybar_start,
     ytop=seq(ybar_start+ybar_inc,ybar_start+ybar_inc*6,ybar_inc))
text(x=seq(0,5*bar_width,bar_width) + start_pos + bar_width/2,
     y=ybar_start, adj=c(0.5,-0.25),
     labels=paste0(seq(0,50,10), "%"), cex=0.5)
arrows(x0=start_pos, x1=6*bar_width+start_pos,
       y0=arr_y, y1=arr_y, length=0.05, lwd=2)
text(x=start_pos + 3*bar_width, y=arr_y,
     adj=c(0.5,1.5), labels="Edge Masking", cex=0.5)

X_ARR <- nbars+ncol(all_amis)+1
ARR_LEN <- c(0.4, 0.3)
arrows(x0=X_ARR, x1=X_ARR, y0=0, y1=ARR_LEN[2], length=HEAD_LEN, col='#45A649', lwd=2, xpd=NA)
# arrows(x0=X_ARR, x1=X_ARR, y0=0, y1=-1*ARR_LEN[1], length=HEAD_LEN, col='#D81B60', lwd=2, xpd=NA)
text(x=X_ARR, y=0, srt=270, col='#45A649', font=2, labels="ExoLabel better",
     adj=c(1.5,TXT_YOFF), cex=0.5, xpd=NA)
# text(x=X_ARR, y=0, srt=270, col='#D81B60', font=2, labels="ExoLabel worse",
#      adj=c(-0.1,TXT_YOFF), cex=0.5, xpd=NA)

legend(y=0.35, x=1, legend=c("Median", "Middle 50% Range"),
       lty=c(NA, 1), pch=c(1, NA), lwd=c(1,2), cex=0.75,
       pt.cex=1)

dev.off(dev.list()['pdf'])
