RealDataResults <- list()
GraphMeans <- list()
genome_levels <- c(50,
                   100,
                   200,
                   400)

PDF_SCALE <- 1.25
PDF_DIM <- c(height=4.3*PDF_SCALE, width=4.3*1.5*PDF_SCALE)
SOURCE_DIR <- '~'

for(i in seq_along(genome_levels)){
  fp <- file.path(SOURCE_DIR, paste0(genome_levels[i], "GenomeSetBlast_v5.RData"))
  if(!file.exists(fp)) stop("File for genome level ", genome_levels[i], " does not exist!")
  load(fp)
  RealDataResults[[paste(genome_levels[i], "Genomes")]] <- RealData
  GraphMeans[[i]] <- round(colMeans(GraphStats))

  to_print <- unique(c(to_print, do.call(rbind, RealData)$Alg))
  #GraphMeans[[i]] <- GraphStats
}
outfile <- file.path(SOURCE_DIR, "RealResults_full.pdf")

NUM_PLOTS <- length(RealDataResults)
ynam <- c("avg_PID")
PANELS_PER_PLOT <- length(ynam)
xnam <- c("mean_clusts")
PLOT_SCALE <- 1L

to_plot <- do.call(rbind, RealData)
to_print <- unique(to_plot$Alg)
cols <- c("LabelProp"='#45A649',
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
          "Spinglass"='#E0A608') # spinglass

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

pdf(outfile, width=PDF_DIM["width"], height=PDF_DIM["height"], onefile=TRUE)
#ynam <- c("unweighted", "annot_nmi")
LAYOUT_MAT <- matrix(c(1,3,2,4,5,6), nrow=2)
layout(LAYOUT_MAT)
par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
xlims <- list(c(5,15),
              c(5,25),
              c(5,30),
              c(5,35))

all_runtimes <- NULL
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
    plot(y=yv, x=xv,
         xlab=c("", "Mean Cluster Size")[(j%/%3)+1],
         ylab=c("", "Mean Within-Cluster PID")[j%%2+1],
         col=cols, pch=pchs,
         xlim=xlims[[j]],
         ylim=c(20,45),
         yaxt='n',
         main = paste(genome_levels[j], "Genomes"))
    axis(2, at=seq(20,44,4), labels = paste0(seq(20,44,4), '%'))
    for(category in CATEGORIES){
      p <- which(grepl(category, to_print))
      cat_col <- cols[p][1]
      lines(y=yv[p[order(xv[p], yv[p])]], x=xv[p][order(xv[p], yv[p])], col=cat_col, lwd=1)
    }
    points(y=yv, x=xv, col=cols, pch=pchs)

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
pchs <- pchs[reorg]
cols <- cols[reorg]
GraphMeans <- do.call(rbind, GraphMeans)
GraphMeans[4,] <- GraphMeans[4,]*11 ## Temporary until we get fixed data
xtick <- seq(1,14)
ytick <- seq(0,5)
plot(y=c(all_runtimes[-nrow(all_runtimes),]), x=rep(GraphMeans[,1], each=nrow(all_runtimes)-1),
     xlab="Number of Nodes (100,000s)",
     ylab="Runtime (sec.)",
     col=cols, pch=pchs,
     log='xy',
     xlim=c(100000*xtick[1], 100000*xtick[length(xtick)]),
     ylim=c(10**ytick[1], 10**ytick[length(ytick)]),
     main = "Runtime", xaxt='n', yaxt='n')
axis(1, at=1e5*xtick, labels=xtick)
axis(2, at=10**ytick, labels=parse(text=paste0(10, "^", ytick)))

for(i in seq_len(nrow(all_runtimes)-1)){
  lines(y=all_runtimes[i,], x=GraphMeans[,1], col=cols[i])
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
