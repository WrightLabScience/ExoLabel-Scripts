#nohup bash -c 'Rscript --vanilla ~/Desktop/ExoBench/BenchmarkAllAlgs.R' > '/home/eswright/Desktop/ExoBench/output_v9.txt' 2>&1 &
DEBUG <- FALSE
library(igraph)
MCL_PATH <- "mcl"
HIPMCL_PATH <- "/home/eswright/Desktop/ExoLabel/hipmcl/bin/hipmcl"
outdir <- '~/Desktop/ExoBench/'
TIME_CMD <- "/usr/bin/time" #gtime for mac
TIMEOUT_CMD <- "timeout"
STARTING_VERT_POW <- 3 # 10^pow
POW_INCREMENT <- 0.5
EDGES_PER_VERT <- 100
N_REPS <- 1
# time before skipping algorithm (in seconds, 60*60 = 1hr)
SKIP_VALS <- c(15, 60, 15*60, 60*60, 6*60*60, 12*60*60, 24*60*60, 48*60*60)
#SKIP_VALS <- c(5, 10, 15, 20, 100)
SKIP_CTR <- 1L
SKIP_MEM <- 400 * 1000000 # in KB, so 1M = 1GB
# number of times above SKIP_TIME before we disqualify an algorithm
NUM_FAILURES_TO_SKIP <- 1
if(DEBUG){
  outdir <- "~/Nextcloud/grad_school/PapersBeingSubmitted/ExoLabel/FigureGeneration/BenchmarkingRuntime/"
  HIPMCL_PATH <- "~/Documents/github/hipmcl/bin/hipmcl"
  MCL_PATH <- "~/local/bin/mcl"
  TIME_CMD <- "gtime"
  TIMEOUT_CMD <- "gtimeout"
  SKIP_MEM <- 16*1000000
  SKIP_VALS <- c(5,10,15,30,60,5*60)
}
PLOT_CALL <- paste("Rscript", file.path(outdir, "PlotRTBenchmarkResults.R"))

igraph_file <- tempfile(tmpdir=outdir)
igraph_file <- file.path(outdir, "igraph_file.tmp")
mcl_file <- tempfile(tmpdir=outdir)
exolabel_file <- tempfile(tmpdir=outdir)
exolabel_statfile <- tempfile(tmpdir=outdir)
generated_graphs <- character(0L)
PLOT_GENERATION_SCRIPT <- file.path(outdir, "PlotRTBenchmarkResults.R")
TEMP_GRAPH_DIR <- file.path(outdir, "GRAPHS")
if(!dir.exists(TEMP_GRAPH_DIR)) dir.create(TEMP_GRAPH_DIR)
outfile <- file.path(outdir, "ExoBenchResults.RData")


set.seed(555L)


algs <- c("igraph-LP", "igraph-Eigen", "igraph-FG",
          "igraph-Louvain", "igraph-InfoMap", "igraph-Walktrap",
          "igraph-Spinglass", "speakeasy-v2", "igraph-Leiden",
          "ExoLabel", "MCL1.4", "MCL2.0", "HipMCL1.4", "HipMCL2.0")
ids <- c("LP", "Eig", "FG", "Lou",
         "IM", "Walk", "Spin", "Speak", "Leid",
         "Exo", "M14", "M20", "HM14", "HM20")
if(!dir.exists(outdir)) dir.create(outdir)
if(file.exists(outfile)){
  load(outfile)
  nr <- nrow(RuntimeData)
  alg_control <- apply(RuntimeData, 2L, \(x){
    c(which(is.na(x)),1)[1] - 1
  })
  alg_control <- pmin(alg_control,
                     apply(MemoryData, 2L, \(x){
                       c(which(is.na(x) | x >= SKIP_MEM),1)[1]-1
                     }))
} else {
  GraphStats <- as.data.frame(matrix(nrow=0, ncol=5))
  colnames(GraphStats) <- c("Seed", "NVert", "NEdge", "Filepath", "ActualNodes")
  RuntimeData <- as.data.frame(matrix(nrow=0, ncol=length(algs)))
  colnames(RuntimeData) <- c(algs)
  MemoryData <- RuntimeData
  alg_control <- rep(0L, length(algs))
}
SKIP_VALS <- c(15, 60, 15*60, 60*60, 6*60*60, 12*60*60, 24*60*60, 48*60*60)
parse_output <- function(s){
  s <- strsplit(s, ' ')[[1]]
  rt <- s[1]
  mem <- as.numeric(s[2])
  rt <- strsplit(rt, ':')[[1]]
  rt <- rev(as.numeric(rt))
  rt <- sum(rt * (60^(seq(0, length(rt)-1L))))
  return(c(rt, mem))
}

build_igraph_command <- function(arg, path, extra_args=character(0L)){
  v <- paste0("suppressPackageStartupMessages(library(igraph));",
              "x <- communities(",
               arg,
               "(graph_from_data_frame(read.delim('",
               path,
               "', header = FALSE), directed=FALSE)",
              extra_args,
               "));",
              "writeLines(unlist(x, use.names=FALSE), con=\'", igraph_file, "\')")
  paste0('R -e "', v, '"')
}

build_speakeasy_command <- function(path, args=''){
  v <- paste0("suppressPackageStartupMessages(library(igraph));",
              "x <- ",
              "speakeasyR::cluster(graph_from_data_frame(read.delim('",
              path,
              "', header = FALSE), directed=FALSE)",
              args,
              "));",
              "writeLines(as.character(x), con=\'", igraph_file, "\')")
  paste0('R -e "', v, '"')
}

validate_igraph <- function(expected_num_nodes){
  if(!file.exists(igraph_file)) return(FALSE)
  x <- readLines(igraph_file)
  file.remove(igraph_file)
  if(length(x) != expected_num_nodes) return(FALSE)
  TRUE
}

validate_mcl <- function(expected_num_nodes){
  if(!file.exists(mcl_file)) return(FALSE)
  x <- readLines(mcl_file)
  file.remove(mcl_file)
  if(length(unlist(strsplit(x, ' '))) != expected_num_nodes) return(FALSE)
  TRUE
}

validate_exolabel <- function(dummyarg){
  if(!file.exists(exolabel_statfile)){
    return(FALSE)
  }
  res <- readRDS(exolabel_statfile)
  clustered <- read.delim(res$results, header=FALSE)
  file.remove(exolabel_statfile)
  file.remove(res$results)
  if(nrow(clustered) != res$graph_stats[1]) return(FALSE)
  TRUE
}

validation_func <- rep(
  c(validate_igraph, validate_exolabel, validate_mcl),
  times=c(9, 1, 4)
)

cmds <- c("cluster_label_prop", "cluster_leading_eigen", "cluster_fast_greedy",
          "cluster_louvain", "cluster_infomap", "cluster_walktrap", "cluster_spinglass")

tf <- tempfile(tmpdir=outdir)
cmdfile <- tempfile(tmpdir=outdir)
cmds <- vapply(cmds, \(x) build_igraph_command(x, tf), character(1L))
cmds <- c(cmds, build_speakeasy_command(tf))
cmds <- c(cmds, build_igraph_command("cluster_leiden", tf, ", objective_function=\'modularity\'"))
cmds <- c(cmds, paste0('R -e "suppressPackageStartupMessages(library(SynExtend));saveRDS(ExoLabel(\'',tf,
                       '\', attenuation=FALSE, verbose=FALSE, outfile=\'', exolabel_file, '\'), file=\'', exolabel_statfile ,'\')\"'))
cmds <- c(cmds, paste0(MCL_PATH, ' ', tf, " --abc -I 1.4 -o ", mcl_file))
cmds <- c(cmds, paste0(MCL_PATH, ' ', tf, " --abc -I 2.0 -o ", mcl_file))
cmds <- c(cmds, paste0(HIPMCL_PATH, " -M ", tf, " -I 1.4 -per-process-mem 15 -o ", mcl_file))
cmds <- c(cmds, paste0(HIPMCL_PATH, " -M ", tf, " -I 2.0 -per-process-mem 15 -o ", mcl_file))

run_cmd <- function(cmd, maxtime, nnodes, validation_func){
  options(scipen=999)
  cmd <- paste0("#!/bin/bash\n\n",
                TIME_CMD,
                " -f '%E %M' ",
                cmd)
  cat(cmd, file=cmdfile)
  system(paste("chmod +x", cmdfile))
  options(scipen=0)
  res <- suppressWarnings(system2(TIMEOUT_CMD, args=c(maxtime, "bash", cmdfile),
                                  stderr=TRUE, stdout=TRUE))
  if(!is.null(attr(res,'status')) && attr(res,'status') == 124 || !validation_func(nnodes)){
    return(c(NA_real_, NA_real_))
  }
  # runtime, then memory
  parse_output(res[length(res)])
}

while(SKIP_CTR <= length(SKIP_VALS)){
  SKIP_TIME <- SKIP_VALS[SKIP_CTR]
  ran_one <- FALSE
  for(i in seq_along(algs)){
    next_row <- alg_control[i] + 1
    if(alg_control[i] > 1 &&
       (is.na(RuntimeData[next_row-1, i]) ||
       (RuntimeData[next_row-1,i] >= SKIP_TIME) ||
       (MemoryData[next_row-1,i]) >= SKIP_MEM)){
      cat(gsub('.', '-', ids[i]), '')
      next
    }
    if(next_row > nrow(RuntimeData)){
      GraphStats[next_row,1] <- sample(100000, 1)
      GraphStats[next_row,2] <- floor(10^(STARTING_VERT_POW + POW_INCREMENT*(next_row-1)))
      GraphStats[next_row,3] <- floor(GraphStats[next_row,2] * EDGES_PER_VERT)
      GraphStats[next_row,4] <- file.path(TEMP_GRAPH_DIR, paste0("Graph_", next_row, "_", GraphStats[next_row,1], ".tsv"))
      GraphStats[next_row,5] <- NA_real_
    }
    seed <- GraphStats[next_row,1]
    nv <- GraphStats[next_row,2]
    ne <- GraphStats[next_row,3]
    fp <- GraphStats[next_row,4]
    expected_nodes <- GraphStats[next_row,5]
    cat(ids[i])
    if(fp == ""){
      ## catch graphs that are too large to generate
      next
    } else if(!file.exists(fp)){
      set.seed(seed)
      d <- try(as_data_frame(sample_gnm(nv, ne, directed=FALSE), what='edges'), silent=TRUE)
      if(!is(d, "data.frame")){
        GraphStats[next_row, 4] <- ""
        cat("! ")
        next
      }
      num_edges <- nrow(d)
      d <- cbind(d, round(runif(nrow(d)), 3))
      write.table(d, fp, append=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
      expected_nodes <- length(unique(c(d[,1], d[,2])))
      GraphStats[next_row,5] <- expected_nodes
      rm(d)
    } else if(is.na(expected_nodes)){
      d <- read.delim(fp, header=FALSE)
      expected_nodes <- length(unique(c(d[,1], d[,2])))
      GraphStats[next_row,5] <- expected_nodes
      rm(d)
    }
    cat(' ')
    file.copy(fp, tf, overwrite=TRUE)
    r <- run_cmd(cmds[i], SKIP_TIME, expected_nodes, validation_func[[i]])
    ran_one <- TRUE
    RuntimeData[next_row,i] <- r[1]
    MemoryData[next_row,i] <- r[2]
    alg_control[i] <- next_row
  }
  cat('\n')
  if(!ran_one){
    #source(PLOT_GENERATION_SCRIPT)
    system(PLOT_CALL, ignore.stdout=TRUE, ignore.stderr=TRUE)
    SKIP_CTR <- SKIP_CTR + 1L
    cat("***Time limit is now ", SKIP_VALS[SKIP_CTR], "***\n", sep='')
    nr <- nrow(RuntimeData)
    alg_control <- apply(RuntimeData, 2L, \(x){
      c(which(is.na(x)),1)[1] - 1
    })
    alg_control <- pmin(alg_control,
                   apply(MemoryData, 2L, \(x){
                     c(which(is.na(x) | x >= SKIP_MEM),1)[1]-1
                   }))
  }
  save(RuntimeData, MemoryData, GraphStats,
         STARTING_VERT_POW, POW_INCREMENT,
         SKIP_VALS, SKIP_CTR, SKIP_MEM,
         file=outfile)
}
