##### config files #####
# tar_config_set(script = "_targets_simulation.R",
#                store = "store_simulation",
#                project = "simulation",
#                reporter_make = "verbose",
#                shortcut = FALSE)
# tar_config_set(script = "_targets_scLANE_GLM.R",
#                store = "store_scLANE_GLM",
#                project = "scLANE_GLM",
#                reporter_make = "verbose",
#                shortcut = FALSE, 
#                inherits = "simulation")
# tar_config_set(script = "_targets_tradeSeq.R",
#                store = "store_tradeSeq",
#                project = "tradeSeq",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")

##### setup #####
setwd("/blue/rbacher/j.leary/repos/scLANE-Simulations/")
library(targets)
library(tarchetypes)

##### batch job submission #####
# sbatch -t 80:00:00 -c 1 --mem=100G -J scLANE_sim --account=rbacher --qos=rbacher-b --mail-type=END --mail-user=j.leary@ufl.edu --wrap="module load R; Rscript run.R"

##### simulation pipeline #####
Sys.setenv(TAR_PROJECT = "simulation")
tar_make_future(workers = 6)

##### scLANE GLM pipeline ####
Sys.setenv(TAR_PROJECT = "scLANE_GLM")
tar_make_future(workers = 6)

##### tradeSeq pipeline ####
Sys.setenv(TAR_PROJECT = "tradeSeq")
tar_make_future(workers = 6)
