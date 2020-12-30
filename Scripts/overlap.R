#################
# Project: Neanderthal niches
# Process: Measuring niche overlap
# Author: Marlon E. Cobos
# Date: 2020/12/30
#################

# packages
if (!require(ellipsenm)) {
  devtools::install_github("marlonecobos/ellipsenm")
}

library(kuenm)
library(ellipsenm)
library(raster)
library(rgl)

#####
# directory
setwd("YOUR/WORKING/DIRECTORY")


#####
# pca
pcas <- kuenm_rpca(variables = "Variables/MIS5a", in.format = "ascii", var.scale = T,
                   write.result = T, out.format = "ascii", out.dir = "PCA", project = T,
                   proj.vars = "Variables/MIS4")


#####
# data for overlap analyses
MIS5a <- read.csv("MIS5a_coords_only_fixed_v2.csv")
MIS4 <- read.csv("MIS4_coords_only_fixed_v2.csv")

# raster layers of environmental data
vMIS5a <- pcas$PCRasters_initial[[1:3]]
vMIS4 <- pcas$PCRasters_MIS4[[1:3]]


#####
# preparing overlap objects to perform analyses
MIS5a_niche <- overlap_object(MIS5a, species =  "Site", longitude = "Long",
                              latitude = "Lat", method = "mve1", level = 95,
                              variables = vMIS5a)

MIS4_niche <- overlap_object(MIS4, species =  "Site", longitude = "Long",
                             latitude = "Lat", method = "mve1", level = 95,
                             variables = vMIS4)


#####
# niche overlap analysis with test of significance (MVE)
overlap_mve <- ellipsoid_overlap(MIS5a_niche, MIS4_niche, overlap_type = "back_union",
                                 significance_test = TRUE, replicates = 1000)

summary(overlap_mve)


#####
# plots
val_MIS5a <- na.omit(values(vMIS5a))
val_MIS4 <- na.omit(values(vMIS4))

lims <- apply(rbind(val_MIS5a, val_MIS4), 2, range)
xlim <- lims[, 1]; ylim <- lims[, 2]; zlim <- lims[, 3]

## overlap 3d
cols <- viridis::magma(10)[c(7, 3)]
plot_overlap(overlap_mve, legend = F)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1))

overlap_mve@variable_names
xat <- c(-6, -3, 0, 3)
yat <- c(-8, -4, 0, 4)
zat <- c(-3, -1.5, 0, 1.5)

# move the plot around and then save the view
view3d <- par3d()$userMatrix
saveRDS(view3d, "view3dov.rds")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "mve")
rgl.postscript("overlap_mve.pdf", fmt = "pdf")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
rgl.postscript("overlap_mve_MIS5a.pdf", fmt = "pdf")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS4, col = "gray75", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
rgl.postscript("overlap_mve_MIS4.pdf", fmt = "pdf")

## significance test
jpeg(filename = "significance_overlap_mve.jpg", width = 80, height = 80,
     units = "mm", res = 600)
par(mar = c(4.5, 4.5, 0.5, 0.5))
par(cex = 0.85)
hist(overlap_mve@significance_results$union_random$Niche_1_vs_2$overlap, breaks = seq(0, 0.8, 0.01),
     main = "", xlab = "Overlap", xlim = c(0, 0.8), ylim = c(0, 50))
abline(v = quantile(overlap_mve@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
       col = "darkgreen", lwd = 2, lty = 2)
abline(v = overlap_mve@union_overlap$overlap, col = "green", lwd = 2)
legend("topright", bty = "n", legend = c("Observed", "5% CL"),
       col = c("green", "darkgreen"), lty = c(1, 2), lwd = 2, cex = 0.8)
box(bty = "l")
dev.off()


#####
# preparing overlap objects to perform analyses (CVA)
MIS5a_niche <- overlap_object(MIS5a, species =  "Site", longitude = "Long",
                              latitude = "Lat", method = "covmat", level = 95,
                              variables = vMIS5a)

MIS4_niche <- overlap_object(MIS4, species =  "Site", longitude = "Long",
                             latitude = "Lat", method = "covmat", level = 95,
                             variables = vMIS4)


#####
# niche overlap analysis with test of significance (CVA)
overlap_cva <- ellipsoid_overlap(MIS5a_niche, MIS4_niche, overlap_type = "back_union",
                                significance_test = TRUE, replicates = 1000)

summary(overlap_cva)


#####
# plot
## overlap 3d
plot_overlap(overlap_cva, legend = F)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1))

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "cva")
rgl.postscript("overlap_cva.pdf", fmt = "pdf")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
rgl.postscript("overlap_cva_MIS5a.pdf", fmt = "pdf")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS4, col = "gray75", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "", ylab = "", zlab = "")
rgl.postscript("overlap_cva_MIS4.pdf", fmt = "pdf")

## significance test
jpeg(filename = "significance_overlap_cva.jpg", width = 80, height = 80,
     units = "mm", res = 600)
par(mar = c(4.5, 4.5, 0.5, 0.5))
par(cex = 0.85)
hist(overlap_cva@significance_results$union_random$Niche_1_vs_2$overlap, breaks = seq(0, 1, 0.01),
     main = "", xlab = "Overlap", xlim = c(0, 1), ylim = c(0, 30))
abline(v = quantile(overlap_cva@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
       col = "darkgreen", lwd = 2, lty = 2)
abline(v = overlap_cva@union_overlap$overlap, col = "green", lwd = 2)
legend("topright", bty = "n", legend = c("Observed", "5% CL"),
       col = c("green", "darkgreen"), lty = c(1, 2), lwd = 2, cex = 0.8)
box(bty = "l")
dev.off()

save(overlap_mve, overlap_cva, val_MIS5a, val_MIS4, file = "Will_Overlap.RData")

sink("Overlap_mve.txt")
summary(overlap_mve)
sink()

sink("Overlap_cva.txt")
summary(overlap_cva)
sink()


#####
# movies 3D
par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "mve")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "mve_MIS5a")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_mve, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS4, col = "gray75", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "mve_MIS4")


par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "cva")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS5a, col = "gray65", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "cva_MIS5a")

par3d(params = list(userMatrix = view3d))
plot_overlap(overlap_cva, legend = F, change_labels = T, niche_col = cols, data_col = cols)
plot3d(val_MIS4, col = "gray75", add = TRUE)
axes3d(xat = xat, yat = yat, zat = zat)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1, 1, 1),
           axes = FALSE, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3")
movie3d(spin3d(axis = c(0, 0, 1)), duration = 5, dir = getwd(), movie = "cva_MIS4")
