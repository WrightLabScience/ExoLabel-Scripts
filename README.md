# Scripts for ExoLabel

## Runtime Benchmark

The following scripts were used:

- BenchmarkAlgs_Runtime.R: Used to execute the runtime scaling analysis. Commands are dispatched from R, but are all run on the commandline using `time` or `gtime`.
- PlotRTBenchmarkResults.R: Used to create plots while BenchmarkAlgs_Runtime.R is running. This allowed us to see the progress of the benchmark without having to restart it for each iteration.
  * The specific graphs used for this benchmark are stored in another Zenodo archive (see source publication) due to space constraints.
  * Graph generation is seeded, so rerunning the script should produce the same graphs used in the publication without having to download graphs from the other archive.
- PlotRuntimeScaling.R: Generates Figure 1 in the source publication

Results from the above scripts are stored in ExoBenchResults.RData, which is
available in RealDataBenchmark.tar.gz from Zenodo at DOI https://doi.org/10.5281/zenodo.20652308.

## Real Data Benchmark

The Real Data Benchmark is shown in Fig. 2 and Fig. S2 of the main text.
Generating these plots requires the results of the benchmark, which is stored
in RealDataBenchmark.tar.gz from Zenodo at DOI https://doi.org/10.5281/zenodo.20419479. The following scripts are used:

- GenerateRealData.R: Creates the RData files in RealDataBenchmark.tar.gz by running each algorithm on disjoint sets of genomes and evaluating the results.
- PlotRealData.R: Generates Figure 2 and Supplemental Figure 2

The complete source set of scores from
reciprocal best BLAST hits is available in FullPairwiseBlastScores.tar.gz, and the
source assemblies are stored in Assemblies.tar.gz. All of these archies are available
from Zenodo at DOI https://doi.org/10.5281/zenodo.20419479.

## Synthetic Benchmark

PlotSynthetic.R generates Supplemental Figure 1 in the publication using the data
found in SyntheticBenchmark.tar.gz from Zenodo at DOI https://doi.org/10.5281/zenodo.20419479. The source data for this figure
are available in SyntheticDataSourceGraphs.tar.gz.
