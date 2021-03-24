#################
# Project: Neanderthal niches
# Process: kuenm running (model calibration, final models, projections, post-modeling)
# Author: Marlon E. Cobos
# Date: 2020/12/30
#################

# packages
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(kuenm)) {
  devtools::install_github("marlonecobos/kuenm")
}
library(kuenm)
library(raster)

#####
# directory
setwd("YOUR/WORKING/DIRECTORY")


#####
# preparing data
sps <- c("MIS4_data", "MIS5a_data")

## reading data occs
m4 <- read.csv("MIS4_data/MIS4_coords_only_fixed_v2.csv")
m5 <- read.csv("MIS5a_data/MIS5a_coords_only_fixed_v2.csv")

## reading data vars
vm4 <- stack(list.files("MIS4_data/MIS4_climate/", pattern = ".asc$", full.names = T))
vm5 <- stack(list.files("MIS5a_data/MIS5a_climate/", pattern = ".asc$", full.names = T))

plot(vm4)
plot(vm5)

## check NAs
m4na <- extract(vm4, m4[, 2:3])
m5na <- extract(vm5, m5[, 2:3])

## check positions in geography
cols4 <- rep("blue", nrow(m4))
na4 <- which(is.na(m4na[, 1]))
cols4[na4] <- "red"

cols5 <- rep("blue", nrow(m5))
na5 <- which(is.na(m5na[, 1]))
cols5[na5] <- "red"

x11()
plot(vm4[[1]])
points(m4[, 2:3], col = cols4)

x11()
plot(vm5[[1]])
points(m5[, 2:3], col = cols5)


## splitting training and testing records
kuenm_occsplit(occ = "MIS4_data/MIS4_coords_only_fixed_v2.csv", train.proportion = 0.75,
               method = "random", save = TRUE, name = "MIS4_data/occ")

kuenm_occsplit(occ = "MIS5a_data/MIS5a_coords_only_fixed_v2.csv", train.proportion = 0.75,
               method = "random", save = TRUE, name = "MIS5a_data/occ")

## preparing M variables
kuenm_varcomb(var.dir = "MIS4_data/MIS4_climate", out.dir = "MIS4_data/M_variables",
              min.number = 2, in.format = "ascii", out.format = "ascii")

kuenm_varcomb(var.dir = "MIS5a_data/MIS5a_climate/", out.dir = "MIS5a_data/M_variables",
              min.number = 2, in.format = "ascii", out.format = "ascii")


#####
# model calibration
occ_joint <- paste0(sps, "/occ_joint.csv")
occ_tra <- paste0(sps, "/occ_train.csv")
M_var_dir <- paste0(sps, "/M_variables")
batch_cal <- paste0(sps, "/Candidate_models")
out_dir <- paste0(sps, "/Candidate_Models")
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- c("q", "lq", "lp", "qp", "lqp")
args <- NULL
maxent_path <- "C:/Maxent"
wait <- FALSE
run <- TRUE

occ_test <- paste0(sps, "/occ_test.csv")
out_eval <- paste0(sps, "/Calibration_results")
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE

for (i in 1:length(sps)) {
  ## candidate models
  kuenm_cal(occ.joint = occ_joint[i], occ.tra = occ_tra[i], M.var.dir = M_var_dir[i], batch = batch_cal[i],
            out.dir = out_dir[i], reg.mult = reg_mult, f.clas = f_clas, args = args,
            maxent.path = maxent_path, wait = wait, run = run)

  ## model evaluation and selection
  cal_eval <- kuenm_ceval(path = out_dir[i], occ.joint = occ_joint[i], occ.tra = occ_tra[i], occ.test = occ_test[i],
                          batch = batch_cal[i], out.eval = out_eval[i], threshold = threshold,
                          rand.percent = rand_percent, iterations = iterations, kept = kept,
                          selection = selection, parallel.proc = paral_proc)
}



#####
# final models
batch_fin <- paste0(sps, "/Final_models")
mod_dir <- paste0(sps, "/Final_Models")
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- TRUE
G_var_dir <- paste0(sps, "/G_variables")
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL

for (i in 1:length(sps)) {
  kuenm_mod(occ.joint = occ_joint[i], M.var.dir = M_var_dir[i], out.eval = out_eval[i], batch = batch_fin[i],
            rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir[i],
            out.format = out_format, project = project, G.var.dir = G_var_dir[i], ext.type = ext_type,
            write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,
            args = args, wait = wait1, run = run1)
}


#####
# extrapolation risks
sets_var <- list(c("Set_18", "Set_19"), c("Set_13", "Set_19"))
out_mop <- paste0(sps, "/MOP_results")
percent <- 5
paral <- FALSE

for (i in 1:length(sps)) {
  kuenm_mmop(G.var.dir = G_var_dir[i], M.var.dir = M_var_dir[i], sets.var = sets_var[[i]],
             out.mop = out_mop[i], percent = percent, parallel = paral)
}


#####
# MOP summary
smop_outdir <- paste0(sps, "/MOP_summary")
for (i in 1:length(sps)) {
  dir.create(smop_outdir[i])

  mops <- stack(list.files(out_mop[i], pattern = ".tif$", full.names = TRUE, recursive = TRUE))

  meam <- calc(mops, mean)
  minm <- calc(mops, min)

  writeRaster(meam, filename = paste0(smop_outdir[i], "/Mean_MOP.tif"), format = "GTiff")
  writeRaster(minm, filename = paste0(smop_outdir[i], "/Min_MOP.tif"), format = "GTiff")
}


#####
# model statistics
format <- "asc"
project <- TRUE
stats <- c("med", "range")
rep <- TRUE
ext_type <- c("E", "EC", "NE") # the type of extrapolation can be selected according to user requirements
out_dir <- paste0(sps, "/Final_Model_Stats")

for (i in 1:length(sps)) {
  sp_name <- as.character(read.csv(occ_joint[i])[1, 1])
  scenarios <- dir(dir(G_var_dir[i], full.names = TRUE)[1])

  kuenm_modstats(sp.name = sp_name, fmod.dir = mod_dir[i], format = format, project = project,
                 statistics = stats, replicated = rep, proj.scenarios = scenarios,
                 ext.type = ext_type, out.dir = out_dir[i])
}


#####
# projection changes
time_per <- c("MIS5a", "MIS4")
out_dir1 <-  paste0(sps, "/Projection_Changes")

for (i in 1:length(sps)) {
  kuenm_projchanges(occ = occ_joint[i], fmod.stats = out_dir[i], threshold = threshold,
                    time.periods = time_per[i], ext.type = ext_type,
                    out.dir = out_dir1[i])
}
