# Run the whole module: the four pipeline steps, then the four figure scripts.
#
#   Rscript run_all.R            # everything
#   Rscript run_all.R figures    # figures only (needs data/, not raw_data/)
#
# Each script runs in its own R process, so nothing leaks from one step to the
# next. The GitHub version has no raw_data/, so only "figures" works there.

.root <- local({
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(hit)) {
    dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/"))
  } else {
    getwd()
  }
  while (!file.exists(file.path(d, "common.R")) && dirname(d) != d) d <- dirname(d)
  d
})

what <- commandArgs(trailingOnly = TRUE)
what <- if (!length(what)) "all" else what[1]

pipeline_scripts <- c(
  "pipeline/01_hypoxia_scores.R",
  "pipeline/02_mri_enhancement_stats.R",
  "pipeline/03_immune_hypoxia_stats.R",
  "pipeline/04_immune_mri_stats.R"
)
figure_scripts <- c(
  "figures/01_scores_by_group.R",
  "figures/02_scores_by_mri.R",
  "figures/03_immune_hypoxia_scatter.R",
  "figures/04_immune_by_mri.R"
)

scripts <- switch(
  what,
  all      = c(pipeline_scripts, figure_scripts),
  pipeline = pipeline_scripts,
  figures  = figure_scripts,
  stop("Unknown argument: ", what, ". Use all, pipeline or figures.")
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

for (s in scripts) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, c(shQuote(file.path(.root, s))))
  if (!identical(status, 0L)) {
    stop("Failed: ", s, " (exit status ", status, ")")
  }
}

message("\nAll done.")
