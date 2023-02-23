## if standard stan install did not work; run below
Sys.getenv("BINPREF")
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
install.packages("jsonlite", type = "source")
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

## check to see if stan is installed correctly
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\n CXX14FLAGS += -mtune=native -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2",
    file = M, sep = "\n", append = FALSE)
example(stan_model, package = "rstan", run.dontrun = TRUE)

## to install cmdstanr to speed up process, run below
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
install_cmdstan(cores = 8)
# CmdStan path set to: C:/Users/AmirNasrollahzadeh/OneDrive - Blend 360/Documents/.cmdstan/cmdstan-2.29.2