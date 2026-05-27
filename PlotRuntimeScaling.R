sourcedir <- '~/Desktop/ExoBench'
if(!all(c("RuntimeData", "MemoryData") %in% ls()))
  load(file.path(sourcedir, "ExoBenchResults.RData"))

## height of the pdf, in inches
SCALE_FACTOR <- 4.3

## Global sizing
TEXT_CEX <- 0.85
PT_CEX <- 0.9

## Bounds of the plots.
## Note that the second value will be increased if it's not large enough
##  to show all datapoints.
RUNTIME_BOUNDS <- c(5e-1, 1e5)
MEM_BOUNDS <- c(5e1,1e6) * 1024 ## y-axis is in MB

#SKIP_ALGS <- c("HipMCL 1.4", "HipMCL 2.0")
SKIP_ALGS <- c("Spinglass", "Speakeasy")

cols <- c('#E0A608','#E0A608', # label propagation
          '#1E88E5','#0E68C5','#0048A5','#0028C5', # modularity
          'black', 'gray60', # other random walks
          '#E82B70','#C80B50','#A80030','#780000', # MCLs
          '#E0A608', # spinglass
          "#824484") # ExoLabel
pchs <- c(16,17,15,18)
plot_ids <- list(
  LabelPropagation=c("igraph-LP", "speakeasy-v2"),
  Modularity=c("igraph-Eigen", "igraph-FG", "igraph-Louvain", "igraph-Leiden"),
  RandomWalk=c("igraph-InfoMap", "igraph-Walktrap"),
  MCL=c("MCL1.4", "MCL2.0", "HipMCL1.4", "HipMCL2.0"),
  Other=c("igraph-Spinglass"),
  ExoLabel=c("ExoLabel")
)

pchs <- unlist(lapply(plot_ids, \(x) pchs[seq_along(x)]), use.names=FALSE)
p_reorg <- match(unlist(plot_ids, use.names=FALSE), colnames(RuntimeData))

PlotMemoryData <- MemoryData[,p_reorg]
PlotRuntimeData <- RuntimeData[,p_reorg]


new_labs <- c("Label Prop", "Speakeasy",
              "Leading Eigen", "FastGreedy", "Louvain", "Leiden",
              "InfoMap", "Walktrap",
              "MCL 1.4" , "MCL 2.0", "HipMCL 1.4", "HipMCL 2.0",
              "Spinglass",
              "ExoLabel")

pos_no_data <- which(apply(PlotRuntimeData, 2, \(x) all(is.na(x))))
names_no_data <- new_labs[pos_no_data]
SKIP_ALGS <- unique(c(SKIP_ALGS, names_no_data))

pdf(file.path(sourcedir, paste0("Figure_1bc", SKIP_CTR, ".pdf")),
    width=SCALE_FACTOR*1.5, height=SCALE_FACTOR*1.5, onefile=TRUE)
layout(matrix(c(1,4,2,3), ncol=2))
exponents <- rep(0L, length(new_labs))
par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0), cex.lab=TEXT_CEX, cex.axis=TEXT_CEX)
all_exponents <- list()
for(plotnum in seq_len(2)){
  d <- PlotRuntimeData
  r <- RUNTIME_BOUNDS
  if(plotnum == 2){
    r <- MEM_BOUNDS
    d <- PlotMemoryData
  }
  for(i in seq_len(ncol(d))){
    y <- d[,i]
    r[2] <- max(r[2], y, na.rm = TRUE)
  }
  xlim_v <- range(GraphStats[,2])
  plot(NULL,
       ylim=r,
       xlim=xlim_v,
       log='xy',
       xlab="Number of Nodes",
       ylab=c("Time (seconds)", "Memory (MB)")[plotnum],
       main=c("Runtime Scaling", "Memory Scaling")[plotnum],
       xaxt='n', yaxt='n')
  labx <- parse(text=paste0(10, "^", 2:8))
  axis(1, at=10**(2:8), labels=labx)
  ytick <- floor(log(r, 10))
  ytick <- 10**seq(ytick[1], ytick[2])
  if(plotnum == 1){
    laby <- log(ytick, base=10)
    laby <- parse(text=paste0(10, "^", laby))
    axis(2, at=ytick, labels=laby)
  } else {
    laby <- log(ytick / 1000, base=10)
    laby <- parse(text=paste0(10, "^", laby))
    axis(2, at=ytick, labels=laby)
  }
  if(plotnum==1){
    time_labs <- c("1 minute", "1 hour", "1 day", "1 week")
    time_vals <- c(60, 60*60, 24*60*60, 7*24*60*60)
    for(time_i in seq_along(time_labs)){
      line_y <- time_vals[time_i]
      text_val <- time_labs[time_i]
      abline(h=line_y, lty=2)
      text(x=xlim_v[2]*1.4, y=line_y, labels=text_val, adj=c(1,-0.3), cex=TEXT_CEX, font=2)
    }
  } else {
    for(mem_level in c(1, 16, 64, 256)){
      line_y <- mem_level * 1e6
      text_val <- paste(mem_level, "GB")
      abline(h=line_y, lty=2)
      text(x=xlim_v[2]*1.4, y=line_y, labels=paste(mem_level, "GB"), adj=c(1,-0.3), cex=TEXT_CEX, font=2)
    }
  }
  for(i in seq_len(ncol(d))){
    if(new_labs[i] %in% SKIP_ALGS) next
    y <- d[,i]
    x <- GraphStats[,2][!is.na(y)]
    y <- y[!is.na(y)]
    y[y == 0] <- r[1]
    o <- order(x)
    y <- y[o]
    x <- x[o]
    lines(x=x, y=y, col=cols[i], lwd=2, lty=1)
    points(x=x, y=y, pch=pchs[i], col=cols[i], cex=PT_CEX)

    if(length(x) >= 4){
      tp <- tail(y, 4)
      tpx <- tail(x, 4)
      exponents[i] <- lm(log(tp) ~ log(tpx))$coefficients[[2]]
    }
  }

  all_exponents[[plotnum]] <- exponents
}

##############
### LEGEND ###
##############

plot.new()
p_noprint <- match(SKIP_ALGS, new_labs)
leg_x <- -0.05
leg_y <- 0.9
x_off <- c(0.50, 0.3)
font_runtime <- rep(1, length(all_exponents[[1]][-p_noprint]))
font_memory <- rep(1, length(all_exponents[[2]][-p_noprint]))
font_runtime[which.min(all_exponents[[1]][-p_noprint])] <- 2
font_memory[which.min(all_exponents[[2]][-p_noprint])] <- 2
legend(x = leg_x, y=leg_y,
       col=cols[-p_noprint], lty=1, pch=pchs[-p_noprint], lwd=2,
       legend=new_labs[-p_noprint], cex=TEXT_CEX, bty='n', xpd=NA,
       title = "Algorithm\nName", title.font=2, title.adj=0.5, pt.cex=PT_CEX)
legend(x = leg_x + x_off[1], y=leg_y,
       legend=sprintf("%1.02f", all_exponents[[1]][-p_noprint]),
       cex=TEXT_CEX, bty='n', xpd=NA, text.font = font_memory,
       title = "Runtime\nScaling", title.font=2, title.adj=0.5, adj=0.3)
legend(x = leg_x + sum(x_off), y=leg_y,
       legend=sprintf("%1.02f", all_exponents[[2]][-p_noprint]),
       cex=TEXT_CEX, bty='n', xpd=NA, text.font = font_runtime,
       title = "Memory\nScaling", title.font=2, title.adj=0.5, adj=0.3)
