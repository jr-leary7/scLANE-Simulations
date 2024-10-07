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
# tar_config_set(script = "_targets_scLANE_GEE.R",
#                store = "store_scLANE_GEE",
#                project = "scLANE_GEE",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_scLANE_GLMM.R",
#                store = "store_scLANE_GLMM",
#                project = "scLANE_GLMM",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_tradeSeq.R",
#                store = "store_tradeSeq",
#                project = "tradeSeq",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_Lamian_GEE.R",
#                store = "store_Lamian_GEE",
#                project = "Lamian_GEE",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_Lamian_GLMM.R",
#                store = "store_Lamian_GLMM",
#                project = "Lamian_GLMM",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_analysis.R",
#                store = "store_analysis",
#                project = "analysis",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")
# tar_config_set(script = "_targets_PseudotimeDE_GLM.R",
#                store = "store_PseudotimeDE_GLM",
#                project = "PseudotimeDE_GLM",
#                reporter_make = "verbose",
#                shortcut = FALSE,
#                inherits = "simulation")

##### setup #####
setwd("/blue/rbacher/j.leary/repos/scLANE-Simulations/")
library(targets)
library(tarchetypes)

##### batch job submission #####
# sbatch -t 80:00:00 --cpus-per-task=8 --ntasks=1 --mem=100G -J scLANE_sim --account=rbacher --qos=rbacher-b --mail-type=END --mail-user=j.leary@ufl.edu --wrap="module load R; Rscript run.R"
# sbatch -t 80:00:00 --cpus-per-task=25 --ntasks=1 --mem=312G -J scLANE_sim --account=rbacher --qos=rbacher-b --mail-type=END --mail-user=j.leary@ufl.edu --wrap="module load R; Rscript run.R"

##### simulation pipeline #####
# Sys.setenv(TAR_PROJECT = "simulation")
# tar_make_future(workers = 6)

##### scLANE GLM pipeline ####
# Sys.setenv(TAR_PROJECT = "scLANE_GLM")
# tar_make_future(workers = 25)

##### scLANE GLMM pipeline ####
Sys.setenv(TAR_PROJECT = "scLANE_GLMM")
tar_make_future(workers = 25)

##### PseudotimDE GLM pipeline ####
Sys.setenv(TAR_PROJECT = "PseudotimeDE_GLM")
tar_make_future(workers = 25)

##### tradeSeq pipeline ####
# Sys.setenv(TAR_PROJECT = "tradeSeq")
# tar_make_future(workers = 25)

##### scLANE GEE pipeline ####
# Sys.setenv(TAR_PROJECT = "scLANE_GEE")
# tar_make_future(workers = 25)

##### Lamian GEE pipeline #####
Sys.setenv(TAR_PROJECT = "Lamian_GEE")
tar_make_future(workers = 25)

##### Lamian GLMM pipeline #####
# Sys.setenv(TAR_PROJECT = "Lamian_GLMM")
# tar_make_future(workers = 6)

##### downstream analysis pipeline #####
Sys.setenv(TAR_PROJECT = "analysis")
tar_make_future(workers = 2)
