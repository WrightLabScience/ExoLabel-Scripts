# 305,615 nodes
# 4,407,504 edges
library(SynExtend)
library(igraph)
library(aricode)
set.seed(987L)

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 0){
  NUM_GENOMES <- as.integer(args[1])
} else {
  NUM_GENOMES <- 50L
}
VERSION_NUM <- 4
MCL_CUTOFF <- 100

source_dir <- file.path(paste0(NUM_GENOMES, 'GenomeSetBlast'))
pid_dir <- file.path(paste0(NUM_GENOMES, 'GenomeSet'))
out_dir <- "~"
out_file <- file.path(out_dir, paste0(basename(source_dir), "_v", VERSION_NUM, ".RData"))
plot_path <- file.path(out_dir, paste0(basename(source_dir), "_v", VERSION_NUM, "_results.pdf"))
all_files <- list.files(source_dir, full.names=TRUE)

DATA_COLNAMES <- c("unweighted", "weighted_byedge", "weighted_bynode", "weighted_bynm1",
                   "num_singletons", "mean_clusts_wsingle",
                   "unweighted_nozero", "weighted_nozero",
                   "runtime")
xnam <- rep("mean_clusts_wsingle", 4)
ynam <- c("unweighted", "weighted_bynm1", "unweighted_nozero", "weighted_nozero")
TO_PRINT <- seq_along(DATA_COLNAMES)
TO_PRINT <- c(6,7,8,1,2)
LAYOUT_MAT <- matrix(1:4, nrow=2)
PDF_DIM <- c(height=4.3*2, width=4.3*2)
PRINT_PREFIXES <- c("SIZE: ", "U: ", "W: ", "U0: ", "W0: ")
RealData <- list()

eval_clustering <- function(clusts, time=-1){
  ## renumber the clusters
  clusts[,2] <- match(clusts[,2], unique(clusts[,2]))
  all_counts <- table(clusts[,2])
  mean_clust_size <- mean(all_counts[all_counts > 1])
  num_single <- total_genes - sum(all_counts[all_counts != 1])
  mean_clust_wsingle <- mean(c(all_counts[all_counts>1], rep(1, num_single)))
  if(all(all_counts == 1)) return(c(0,0,num_single))

  clmap <- clusts[,2]
  names(clmap) <- clusts[,1]
  tmp_pid <- all_pids
  tmp_pid[,1] <- clmap[tmp_pid[,1]]
  tmp_pid[,2] <- clmap[tmp_pid[,2]]
  normalizer <- all_counts * (all_counts - 1) / 2

  p <- !is.na(tmp_pid[,1]) & !is.na(tmp_pid[,2]) & tmp_pid[,1] == tmp_pid[,2]
  tmp_pid <- tmp_pid[p,]
  tmp_nozero <- tapply(tmp_pid[,3], tmp_pid[,1], mean)
  avg_pid_nozero <- mean(tmp_nozero)
  avg_pid_nozero_weight <- weighted.mean(tmp_nozero, all_counts[names(tmp_nozero)])
  tmp_pid <- tapply(tmp_pid[,3], tmp_pid[,1], sum)
  weights <- all_counts[names(tmp_pid)]
  avg_pid <- numeric(length(all_counts))
  names(avg_pid) <- names(all_counts)
  avg_pid[names(tmp_pid)] <- tmp_pid
  avg_pid <- avg_pid[names(normalizer)]
  avg_pid <- avg_pid / normalizer
  clust_sizes <- all_counts[names(normalizer)]
  avg_pid[is.na(avg_pid)] <- 1

  c(mean(avg_pid[all_counts != 1]),
    weighted.mean(avg_pid, normalizer),
    weighted.mean(avg_pid, clust_sizes),
    weighted.mean(avg_pid, clust_sizes-1),
    num_single,
    mean_clust_wsingle,
    avg_pid_nozero,
    avg_pid_nozero_weight,
    time)
}

ExoLabel_SL <- c(0,25,50,100)
ExoLabel_Atten <- c(0, 0.5, 1.0)
el_algnames <- character(0)
for(i in seq_along(ExoLabel_Atten)){
  el_algnames <- c(el_algnames,
                   paste0("ExoLabel A=", sprintf("%.01f", ExoLabel_Atten[i]),
                          " SL=", ExoLabel_SL))
}

## igraph
igraph_algs <- list(
  cluster_fast_greedy,
  cluster_label_prop,
  \(x) cluster_louvain(x, resolution=0.5),
  cluster_louvain,
  \(x) cluster_louvain(x, resolution=2),
  \(x) cluster_leiden(x, objective_function="modularity", resolution=0.5),
  \(x) cluster_leiden(x, objective_function="modularity"),
  \(x) cluster_leiden(x, objective_function="modularity", resolution=2),
  cluster_infomap,
  cluster_leading_eigen,
  \(x) cluster_walktrap(x, steps=2),
  cluster_walktrap,
  \(x) cluster_walktrap(x, steps=6)
)
igraph_names <- c(
  "FastGreedy", "LabelProp",
  "Louvain05", "Louvain10", "Louvain20",
  "Leiden05", "Leiden10", "Leiden20",
  "Infomap",
  "LeadingEigen",
  "Walktrap2", "Walktrap4", "Walktrap6"
)

igraph_cutoffs <- c(100, 100, rep(Inf,3), rep(Inf,3), Inf, 10, rep(10,3))

GraphStats <- data.frame(NVert=numeric(length(all_files)), NEdges=numeric(length(all_files)))

igraph_keep <- which(igraph_cutoffs >= NUM_GENOMES)
igraph_names <- igraph_names[igraph_keep]
igraph_algs <- igraph_algs[igraph_keep]
tf <- tempfile()

for(fsi in seq_along(all_files)){
  cat('\n', basename(all_files[fsi]), ':\n')
  set.seed(fsi)
  cur_f <- all_files[fsi]
  cur_f_basename <- basename(cur_f)

  cat("Generating igraph graph...")
  igraph_gen_time <- Sys.time()
  all_pids <- read.delim(file.path(pid_dir, cur_f_basename), header=FALSE)
  g <- graph_from_data_frame(read.delim(cur_f, header=FALSE), directed=FALSE)
  total_genes <- length(V(g))
  igraph_gen_time <- difftime(Sys.time(), igraph_gen_time, units="secs")[[1]]
  cat("done.\n")

  ## ExoLabel
  ELRes <- matrix(NA_real_, nrow=length(ExoLabel_SL)*length(ExoLabel_Atten), ncol=length(DATA_COLNAMES))
  ctr <- 1L
  for(ii in seq_along(ExoLabel_Atten)){
    for(i in seq_along(ExoLabel_SL)){
      cat(el_algnames[ctr])
      t1 <- Sys.time()
      exoclust <- ExoLabel(cur_f, attenuation=ExoLabel_Atten[ii], add_self_loops=ExoLabel_SL[i], verbose=FALSE)
      elapsed <- difftime(Sys.time(), t1, units="secs")[[1]]
      cat(" (done)\n")
      ELRes[ctr,] <- eval_clustering(read.delim(exoclust$results, header=FALSE), elapsed)
      #cat(paste0("EL", ExoLabel_SL[i], " "))
      outp <- paste0(PRINT_PREFIXES, round(ELRes[ctr,TO_PRINT], 3))
      cat(outp, '\n')
      ctr <- ctr + 1L
    }
  }
  colnames(ELRes) <- DATA_COLNAMES
  ELRes <- as.data.frame(ELRes)
  ELRes$Alg <- el_algnames
  ELRes$File <- cur_f_basename
  RealData$ExoLabel <- rbind(RealData$ExoLabel, ELRes)
  #save(RealData, file=out_file)

  GraphStats[fsi,] <- c(total_genes, nrow(all_pids))
  igraph_res <- matrix(NA_real_, nrow=length(igraph_names)+1, ncol=length(DATA_COLNAMES))
  for(i in seq_along(igraph_names)){
    cat(igraph_names[i])
    t1 <- Sys.time()
    cl <- communities(igraph_algs[[i]](g))
    elapsed <- difftime(Sys.time(), t1, units="secs")[[1]] + igraph_gen_time
    cl <- data.frame(Vertex=unlist(cl, use.names=FALSE), Cluster=rep(seq_along(cl), times=lengths(cl)))
    cat(" (done)\n")
    igraph_res[i,] <- eval_clustering(cl, time=elapsed)
    outp <- paste0(PRINT_PREFIXES, round(igraph_res[i,TO_PRINT], 3))
    cat(outp, '\n')
  }

  ## disjoint sets analysis
  cat("Disjoint Set\n")
  t1 <- Sys.time()
  cl <- components(g)$membership
  cl <- data.frame(Vertex=names(cl), Cluster=cl)
  elapsed <- difftime(Sys.time(), t1)
  igraph_res[nrow(igraph_res),] <- eval_clustering(cl, time=elapsed)
  outp <- paste0(PRINT_PREFIXES, round(igraph_res[nrow(igraph_res),TO_PRINT], 3))
  cat(outp, '\n')
  rownames(igraph_res) <- c(igraph_names, "DisjointSet")
  colnames(igraph_res) <- DATA_COLNAMES
  igraph_res <- as.data.frame(igraph_res)
  igraph_res$Alg <- c(igraph_names, "DisjointSet")
  igraph_res$File <- cur_f_basename

  ## do speakeasy here too
  #cl <- speakeasyR::cluster(g)

  RealData$igraph <- rbind(RealData$igraph, igraph_res)
  #save(RealData, file=out_file)

  ## MCL
  if(NUM_GENOMES <= MCL_CUTOFF){
    infl_vals <- c(1.4, 2.0, 3.0, 4.0)
    mcl_res <- matrix(NA_real_, nrow=length(infl_vals), ncol=length(DATA_COLNAMES))
    for(i in seq_along(infl_vals)){
      cmd <- paste0("~/local/bin/mcl ", cur_f, " --abc -I ", infl_vals[i], " -o ", tf)
      cat(paste("MCL", infl_vals[i]))
      t1 <- Sys.time()
      system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
      res <- readLines(tf)
      elapsed <- difftime(Sys.time(), t1, units="secs")[[1]]
      res <- strsplit(res, '\t')
      res <- data.frame(Vertex=unlist(res, use.names=FALSE), Cluster=rep(seq_along(res), times=lengths(res)))
      cat(" (done)\n")
      mcl_res[i,] <- eval_clustering(res, time=elapsed)
      outp <- paste0(PRINT_PREFIXES, round(mcl_res[i,TO_PRINT],3))
      cat(outp, '\n')
      file.remove(tf)
    }
    colnames(mcl_res) <- DATA_COLNAMES
    mcl_res <- as.data.frame(mcl_res)
    mcl_res$Alg <- paste0("MCL I=", infl_vals)
    mcl_res$File <- cur_f_basename
    RealData$mcl <- rbind(RealData$mcl, mcl_res)
  }

  save(RealData, GraphStats, file=out_file)
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

  EL_Cols <- c("#824484", "#5f3260", "#341b35")
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
  CATEGORIES <- c("MCL", paste0("ExoLabel A=", sprintf("%.01f", ExoLabel_Atten)), "Louvain", "Leiden")

  pdf(plot_path, width=PDF_DIM["width"], height=PDF_DIM["height"], onefile=TRUE)
  #ynam <- c("unweighted", "annot_nmi")
  layout(LAYOUT_MAT)
  par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))
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
    plot(y=yv, x=xv,
         xlab=c("", "Mean Cluster Size")[i],
         ylab=ynam[i],
         col=cols, pch=pchs)
    for(category in CATEGORIES){
      p <- which(grepl(category, to_print))
      cat_col <- cols[p][1]
      lines(y=yv[p[order(xv[p], yv[p])]], x=xv[p][order(xv[p], yv[p])], col=cat_col, lwd=1)
    }
    points(y=yv, x=xv, col=cols, pch=pchs)
    if(i == length(ynam)){
      to_print <- paste0(to_print, " (", round(runtimes, 2), " sec.)")
    }
    abline(h=yv["DisjointSet"], v=xv["DisjointSet"], lty=2)
    l_p <- which(pchs==16)
    legend(c("topright", "topright", "bottomleft", "bottomright")[i],
           legend=to_print[l_p], pch=pchs[l_p], col=cols[l_p], cex=0.5, inset=0.025)
  }
  dev.off(dev.list()['pdf'])
}

warnings()
