RealDataResults <- list()
GraphMeans <- list()
genome_levels <- c(50,
                   100,
                   200,
                   400)

PLOT_MAIN <- TRUE # False to plot supplemental
SOURCE_DIR <- '~'

if(PLOT_MAIN){
  PDF_SCALE <- 0.65
  PDF_DIM <- c(height=4.3*PDF_SCALE, width=4.3*3*PDF_SCALE)
  outfile <- file.path(dirname(SOURCE_DIR), "Figure2_main.pdf")
} else {
  PDF_SCALE <- 1.25
  PDF_DIM <- c(height=4.3*PDF_SCALE, width=4.3*1.5*PDF_SCALE)
  outfile <- file.path(dirname(SOURCE_DIR), "Figure2_suppl.pdf")
}

to_print <- character(0L)
for(i in seq_along(genome_levels)){
  fp <- file.path(SOURCE_DIR, paste0(genome_levels[i], "GenomeSetBlast_v5.RData"))
  if(!file.exists(fp)) stop("File for genome level ", genome_levels[i], " does not exist!")
  load(fp)
  RealDataResults[[paste(genome_levels[i], "Genomes")]] <- RealData
  GraphMeans[[i]] <- round(colMeans(GraphStats))

  to_print <- unique(c(to_print, do.call(rbind, RealData)$Alg))
  #GraphMeans[[i]] <- GraphStats
}


NUM_PLOTS <- length(RealDataResults)
ynam <- c("avg_PID")
PANELS_PER_PLOT <- length(ynam)
xnam <- c("mean_clusts")
PLOT_SCALE <- 1L

to_plot <- do.call(rbind, RealData)
to_print <- unique(to_plot$Alg)
cols <- c("LabelProp"='#E0A608',
          "Speakeasy"='#45A649', # label propagation
          "FastGreedy"='#1E88E5',
          "Louvain05"='#1E88E5',
          "Louvain10"='#1E88E5',
          "Louvain20"='#1E88E5',
          "Leiden05"='#0E68C5',
          "Leiden10"='#0E68C5',
          "Leiden20"='#0E68C5',
          "LeadingEigen"='#0048A5', # modularity
          "Walktrap2"='grey40',
          "Walktrap4"='grey40',
          "Walktrap6"='grey40',
          "Infomap"='gray60', # other random walks
          "DisjointSet"="black",
          "MCL I=1.4"='#E82B70',
          "MCL I=2"='#E82B70',
          "MCL I=3"='#E82B70',
          "MCL I=4"='#E82B70', # MCLs
          "Spinglass"='#E0A608') # spinglass, not used

EL_Cols <- c("#824484", NA, NA)
ExoLabel_SL <- c(0,25,50,100)
ExoLabel_Atten <- c(0, 0.5, 1.0)
el_algnames <- character(0)
for(i in seq_along(ExoLabel_Atten)){
  el_algnames <- c(el_algnames,
                   paste0("ExoLabel A=", sprintf("%.01f", ExoLabel_Atten[i]),
                          " SL=", ExoLabel_SL))
}
for(ii in seq_along(ExoLabel_Atten)){
  vel <- el_algnames[grepl(paste0("ExoLabel A=", sprintf("%.01f", ExoLabel_Atten[ii])), el_algnames)]
  for(elname in vel)
    cols[elname] <- EL_Cols[ii]
}
pchs <- c(16,17,15,18)
pchs <- unlist(lapply(rle(cols)$lengths, \(x) pchs[seq_len(x)]))

o <- match(to_print, names(cols))
cols <- cols[o]
pchs <- pchs[o]
CATEGORIES <- c("MCL",
                paste0("ExoLabel A=", sprintf("%.01f", ExoLabel_Atten)),
                "Louvain", "Leiden", "Walktrap")

if (PLOT_MAIN){
  pdf(outfile, width=PDF_DIM["width"], height=PDF_DIM["height"], onefile=TRUE)
  #ynam <- c("unweighted", "annot_nmi")
  LAYOUT_MAT <- matrix(c(1,2,3), nrow=1)
  layout(LAYOUT_MAT)
  par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
  xlims <- list(c(5,15),
                c(5,25),
                c(5,30),
                c(5,40))
} else {
  pdf(outfile, width=PDF_DIM["width"], height=PDF_DIM["height"], onefile=TRUE)
  #ynam <- c("unweighted", "annot_nmi")
  LAYOUT_MAT <- matrix(c(1,3,2,4,5,6), nrow=2)
  layout(LAYOUT_MAT)
  par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
  xlims <- list(c(5,15),
                c(5,25),
                c(5,30),
                c(5,40))
}

all_runtimes <- NULL
plotted_arrow_yet <- FALSE
for(j in seq_along(RealDataResults)){
  to_plot <- do.call(rbind, RealDataResults[[j]])
  for(i in seq_along(ynam)){
    xv <- to_plot[[xnam[i]]]
    yv <- to_plot[[ynam[i]]]

    xv <- tapply(xv, to_plot$Alg, mean)
    yv <- tapply(yv, to_plot$Alg, mean)
    runtimes <- tapply(to_plot$runtime, to_plot$Alg, mean)
    reorg <- match(to_print, names(xv))
    xv <- xv[reorg]
    yv <- yv[reorg]
    runtimes <- runtimes[reorg]
    if(is.null(all_runtimes)){
      all_runtimes <- runtimes
    } else {
      all_runtimes <- cbind(all_runtimes, runtimes)
    }
    if(PLOT_MAIN && j < 4){
      ## Main plot only includes the 400 genome case
      next
    } else if (PLOT_MAIN) {
      plot(y=yv, x=xv,
           xlab="Mean Cluster Size",
           ylab="Mean Within-Cluster PID",
           col=cols, pch=pchs,
           #xlim=xlims[[j]],
           xlim=c(0.8*min(xv),max(xv)),
           ylim=c(20,45),
           yaxt='n',
           main = "Accuracy on 400 Genomes")
    } else {
      plot(y=yv, x=xv,
           xlab=c("", "Mean Cluster Size")[(j%/%3)+1],
           ylab=c("", "Mean Within-Cluster PID")[j%%2+1],
           col=cols, pch=pchs,
           #xlim=xlims[[j]],
           xlim=c(0.8*min(xv),max(xv)),
           ylim=c(20,45),
           yaxt='n',
           main = paste(genome_levels[j], "Genomes"))
    }
    axis(2, at=seq(20,44,4), labels = paste0(seq(20,44,4), '%'))
    for(category in CATEGORIES){
      p <- which(grepl(category, to_print))
      cat_col <- cols[p][1]
      lines(y=yv[p[order(xv[p], yv[p])]], x=xv[p][order(xv[p], yv[p])], col=cat_col, lwd=1)
    }
    points(y=yv, x=xv, col=cols, pch=pchs)
    if(!plotted_arrow_yet){
      plotted_arrow_yet <- TRUE
      abs_arrow_mult <- 0.98
      yv_max <- 44L
      xv_max <- xv["DisjointSet"] * abs_arrow_mult
      arrow_scaling <- c(0.875, 0.9)
      arrows(x0=xv_max*arrow_scaling[1], x1=xv_max,
             y0=yv_max*arrow_scaling[2], y1=yv_max,
             col='#45A649', length=0.05, lwd=2)
      text(x=xv_max*arrow_scaling[1], y=yv_max*arrow_scaling[2],
           labels="Better Performance", col='#45A649', font=1, cex=1,
           srt=0, adj=c(0.65,1.1))
    }

    abline(h=yv["DisjointSet"], v=xv["DisjointSet"], lty=2)
    l_p <- which(pchs==16)
    # legend(c("topright", "topright", "bottomleft", "bottomright")[i],
    #        legend=to_print[l_p], pch=pchs[l_p], col=cols[l_p], cex=0.5, inset=0.025)
  }
}

## Runtimes
pos_disjoint_set <- which(rownames(all_runtimes) == "DisjointSet")
reorg <- c(seq_len(pos_disjoint_set-1), seq(pos_disjoint_set+1, nrow(all_runtimes)), pos_disjoint_set)
all_runtimes <- all_runtimes[reorg,]
# xv <- c(genome_levels, 4422)
xv <- genome_levels
pchs <- pchs[reorg]
cols <- cols[reorg]
GraphMeans <- do.call(rbind, GraphMeans)
#GraphMeans[4,] <- GraphMeans[4,]*11 ## Temporary until we get fixed data
xtick <- seq(1,14)
ytick <- seq(0,5)
to_include <- c(1,9,10,12,15,17,19,24)
runtimes <- all_runtimes[to_include,]
# runtimes <- cbind(runtimes, rep(NA_real_, nrow(runtimes)))
# runtimes[1,ncol(runtimes)] <- 54831
to_plot_x <- c(rep(xv, each=nrow(runtimes)), 4422)
to_plot_y <- c(runtimes, 54831)
# plot(y=c(runtimes), x=rep(xv, each=nrow(runtimes)),
plot(y=to_plot_y, x=to_plot_x,
     xlab="Number of Genomes",
     ylab="Runtime (sec.)",
     col=cols[to_include], pch=pchs[to_include],
     log='xy',
     #xlim=c(100000*xtick[1], 100000*xtick[length(xtick)]),
     ylim=c(10**ytick[1], 10**ytick[length(ytick)]),
     main = "Runtime", yaxt='n', xaxt='n')
#axis(1, at=1e5*xtick, labels=xtick)
axis(1, at=c(50, 100, 200, 400, 800, 1600, 4400))
axis(2, at=10**ytick, labels=parse(text=paste0(10, "^", ytick)))

y_line_marks <- c(1, 60, 60*60, 60*60*24, 60*60*24*7)
y_line_labs <- c("1 sec.", "1 min.", "1 hr.", "1 day", "1 week")
for(i in seq_along(y_line_marks)){
  abline(h=y_line_marks[i], lty=2)
  #text(x=100000, y=y_line_marks[i], labels = y_line_labs[i],
  text(x=5000, y=y_line_marks[i], labels = y_line_labs[i],
       adj=c(1,-0.20), font=2)
}

for(i in seq_len(nrow(runtimes))){
  yvr <- runtimes[i,]
  xvr <- xv
  if(i == 1){
    yvr <- c(yvr, 54831)
    xvr <- c(xvr, 4422)
  }
  lines(y=yvr, x=xvr, col=cols[to_include][i])
  points(y=yvr, x=xvr, col=cols[to_include][i], pch=pchs[to_include][i])
}

plot.new()
p <- which(!is.na(cols))
algnames <- rownames(all_runtimes[p,])
algnames <- c("ExoLabel (SL=0)", "ExoLabel (SL=25)", "ExoLabel (SL=50)", "ExoLabel (SL=100)",
              "FastGreedy", "Label Propagation",
              "Louvain (R=0.5)", "Louvain (R=1.0)", "Louvain (R=2.0)",
              "Leiden (R=0.5)", "Leiden (R=1.0)", "Leiden (R=2.0)",
              "InfoMap",
              "Walktrap (S=2)", "Walktrap (S=4)", "Walktrap (S=6)",
              "MCL (I=1.4)", "MCL (I=2.0)", "MCL (I=3.0)", "MCL (I=4.0)",
              "Naive Disjoint Set")
legend(x=-0.1, y=0.9, legend=algnames,
       pch=pchs[p], col=cols[p], xpd=NA, ncol=2, bty='n')

dev.off(dev.list()['pdf'])
