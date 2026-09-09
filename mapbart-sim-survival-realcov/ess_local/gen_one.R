#!/usr/bin/env Rscript
# gen_one.R --sc <n> --size_option <default|small> [--subsc <E>]
#          [--seed_rct <123>] [--seed_rwd <456>]
#
# Generate iteration 1 ONLY for one scenario, for the ess calibration to run on.
# data_gen_p10_v3.R has rm(list=ls()) and no CLI, so its hardcoded defaults can
# only be overridden by rewriting the assignment lines in its text (the same
# trick run_all.R uses).  --seed_rct / --seed_rwd default to run_all.R's values
# (123 / 456) and are applied the SAME way run_all.R does -- seed_rct rewrites
# the set.seed(<literal>) RCT seed-array generator, seed_rwd overrides the
# `seed_rwd <- ...` RWD-pool seed -- so _1 is identical to the full sweep's
# iteration 1 (which run_all.R then skips regenerating).
#
#   Rscript gen_one.R --sc 1 --size_option default    # n200, sc1E
#   Rscript gen_one.R --sc 3 --size_option small      # n30,  sc3E

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  k <- match(flag, args)
  if (is.na(k) || k == length(args)) return(default)
  args[k + 1]
}
sc_arg    <- get_arg("--sc")
size_arg  <- get_arg("--size_option")
subsc_arg <- get_arg("--subsc", "E")
seed_rct  <- get_arg("--seed_rct", "123")   # matches run_all.R seed_rct
seed_rwd  <- get_arg("--seed_rwd", "456")   # matches run_all.R seed_rwd
if (is.null(sc_arg) || is.null(size_arg))
  stop("usage: Rscript gen_one.R --sc <n> --size_option <default|small> [--subsc <E>] [--seed_rct <123>] [--seed_rwd <456>]")

data_gen <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-realcov/data_gen_p10_v3.R"

src <- readLines(data_gen)
src <- gsub("niter <- 500", "niter <- 1", src)
src <- sub("^sc <- .*",          sprintf("sc <- %s", sc_arg), src)
src <- sub("^subsc <- .*",       sprintf('subsc <- "%s"', subsc_arg), src)
src <- sub("^size_option <- .*", sprintf('size_option <- "%s"', size_arg), src)
# Seeds -- applied exactly as run_all.R's apply_overrides does:
src <- gsub("set\\.seed\\(\\s*\\d+\\s*\\)", sprintf("set.seed(%s)", seed_rct), src, perl = TRUE)
src <- sub("^seed_rwd <- .*", sprintf("seed_rwd <- %s", seed_rwd), src)
eval(parse(text = paste(src, collapse = "\n")))
