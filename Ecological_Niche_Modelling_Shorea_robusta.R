# Author: Kamana Pathak
# Date: 8/15/2025

## ============================================================
## Spatial Analysis & Species Distribution Modeling (SDM) in R
## ============================================================

#### Install and load required packages for SDM and ecological niche modeling

## 1. ENMeval - Evaluate ecological niche models, esp. MaxEnt
install.packages("ENMeval")

## 2. dplyr - Data manipulation (filter, select, mutate, summarize)
install.packages("dplyr")

## 3. geodata - Download climate, elevation, soil, occurrence datasets
install.packages("geodata")

## 4. terra - Modern raster and vector data handling (replacement for 'raster')
install.packages("terra")

## 5. corrplot - Visualize correlation matrices
install.packages("corrplot")

## 6. car & reshape - Regression support & data reshaping
if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car")
  install.packages("reshape")
}

## 7. ggplot2 - Data visualization with grammar of graphics
install.packages("ggplot2")

## 8. usdm - Uncertainty analysis for SDMs
install.packages("usdm")

## 9. ecospat - Spatial ecology and SDM evaluation tools
install.packages("ecospat")

## 10. tibble - Enhanced data frame handling
install.packages("tibble")

## 11. devtools - Install packages from GitHub
if (!require(devtools)) {
  install.packages("devtools")
}

## 12. kuenm - MaxEnt model calibration & evaluation (from GitHub)
if (!require(kuenm)) {
  devtools::install_github("marlonecobos/kuenm")
}

## 13. sp - Spatial vector data handling
install.packages("sp")

## 14. predicts - Environmental similarity metrics (MESS)
install.packages("predicts")

## Load all installed packages
library(ENMeval)
library(dplyr)
library(geodata)
library(terra)
library(corrplot)
library(car)
library(reshape)
library(ggplot2)
library(usdm)
library(ecospat)
library(tibble)
library(sp)
library(kuenm)
library(predicts)

## ============================================================
## Working directory setup
## ============================================================

getwd()  ## Check current working directory
setwd("D:/SDM_Shorea_robusta")  ## Set working directory

set.seed(100)  ## Set seed for reproducibility

## ============================================================
## Load and clean occurrence data
## ============================================================

bv <- read.csv(file ="Shorea_robusta.csv")
bv

### Check for duplicate records
any(duplicated(bv))

### Remove duplicates
occs <- bv[!duplicated(bv),]
occs

### Check again for duplicates
any(duplicated(occs))

### Check for NA values
any(is.na(occs))

### Remove NA occurrences
bv <- na.omit(occs)
bv

### Final check for NA
any(is.na(bv))

## ============================================================
## Download WorldClim bioclimatic variables
## ============================================================

bioclim <- geodata::worldclim_global(var = 'bio', res = 10, download = F,
                                     path = "D:/SDM_Shorea_robusta")
bioclim

### Crop raster to Nepal extent
bioclim_cropped <- crop(bioclim, ext(80, 88.2, 26.3, 30.5))
bioclim_cropped

### Convert raster to data frame for analysis
current_df <- as.data.frame(bioclim_cropped)
current_df <- na.omit(current_df)  ## Remove NA values
current_df

## ============================================================
## Multicollinearity check
## ============================================================

### Calculate variance inflation factor (VIF)
collin <- vifstep(bioclim_cropped, th=10)  ## Threshold = 10
collin

### Extract VIF results
vif <- collin@results
vif

### Save VIF results
write.csv(vif, "D:/SDM_Shorea_robusta/Vif_results.csv", row.names = FALSE)

### Get correlation matrix
corMatrix <- collin@corMatrix
corMatrix

### Exclude collinear variables
analogClimExcluded <- exclude(bioclim_cropped, collin)

### Keep only selected 8 bioclimatic variables
vars_to_keep <- c(
  "wc2.1_10m_bio_2", 
  "wc2.1_10m_bio_3", 
  "wc2.1_10m_bio_5", 
  "wc2.1_10m_bio_13", 
  "wc2.1_10m_bio_14", 
  "wc2.1_10m_bio_15", 
  "wc2.1_10m_bio_18", 
  "wc2.1_10m_bio_19"
)
analogClimExcluded <- analogClimExcluded[[vars_to_keep]]
print(names(analogClimExcluded))
nlyr(analogClimExcluded)  ## Should show 8 layers

### Stack layers
envs <- rast(analogClimExcluded)
plot(envs)

## ============================================================
## Prepare occurrence points
## ============================================================

### Fix coordinate columns
occs_fixed <- data.frame(
  Lon = as.numeric(occs$Lon),
  Lat = as.numeric(occs$Lat)
)

### Remove rows with NA coordinates
occs_fixed <- occs_fixed[!is.na(occs_fixed$Lon) & !is.na(occs_fixed$Lat), ]

### Extract raster cell numbers and remove duplicates
occs.cells <- raster::extract(envs[[1]], occs_fixed, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs_filtered <- occs_fixed[!occs.cellDups, ]
occs <- as.data.frame(occs_filtered)
occs

### Plot first environmental layer with occurrence points
plot(envs[[1]])
points(occs, col="red")

## ============================================================
## Extract environmental values at occurrence points
## ============================================================

occs.z <- raster::extract(envs, occs)

### Remove rows with NA
occs.na <- which(rowSums(is.na(occs.z)) > 0)
occs <- occs[-occs.na,]
occs.z <- na.omit(occs.z)

### Convert envs to terra SpatRaster
envs_terra <- rast(envs)
occs.z <- terra::extract(envs_terra, occs[, c("Lon", "Lat")])
na_rows <- which(rowSums(is.na(occs.z)) > 0)
if(length(na_rows) > 0){
  occs <- occs[-na_rows, ]
  occs.z <- occs.z[-na_rows, ]
}
if("ID" %in% colnames(occs.z)) {
  occs.z <- occs.z[, !(colnames(occs.z) == "ID")]
}

## ============================================================
## MESS analysis (Environmental similarity)
## ============================================================

occs.mess <- predicts::mess(envs_terra, occs.z)
plot(occs.mess, main = "Environmental similarity")
points(occs$Lon, occs$Lat, pch=20, col="red")

## ============================================================
## Create SpatialPoints and sf objects for occurrences
## ============================================================

occs.sp <- sp::SpatialPoints(occs)
proj4string(occs.sp) <- CRS("+proj=longlat +datum=WGS84")

occs.sf <- sf::st_as_sf(occs, coords = c("Lon","Lat"), crs = raster::crs(envs))
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

## ============================================================
## Buffer occurrence points and create background extent
## ============================================================

occs.buf <- sf::st_buffer(occs.sf, dist = 1000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>% 
  sf::st_transform((crs = raster::crs(envs)))
plot(envs[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

### Crop and mask environmental raster to study extent
envs.bg <- raster::crop(envs, occs.buf)
envs.bg <- raster::mask(envs.bg, occs.buf)
plot(envs.bg)
plot(envs.bg[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border ="blue", lwd =3, add = TRUE)

## ============================================================
## Generate background points
## ============================================================

bg <- dismo::randomPoints(envs.bg[[1]], n = 16000)
bg <- as.data.frame(bg)
colnames(bg) <- colnames(occs)
plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)

## ============================================================
## Partition occurrences for evaluation (block method)
## ============================================================

occs.z <- cbind(occs, raster::extract(envs, occs))
bg.z <- cbind(bg, raster::extract(envs, bg))

block <- get.block(occs, bg, orientation = "lat_lon")
table(block$occs.grp)

envs.bg <- terra::rast(envs.bg)
evalplot.grps(pts = occs, pts.grp = block$occs.grp, envs = envs.bg) + 
  ggtitle("Spatial block partitions: occurrences")

## ============================================================
## Run simple MaxNet model with linear feature class
## ============================================================

envs <- terra::rast(envs)
e.mx <- ENMevaluate(
  occs = occs,
  envs = envs,
  bg = bg,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(
    fc = c("L"),  # Linear features only
    rm = 1:5      # Regularization multipliers
  )
)
e.mx

## ============================================================
## Evaluate model performance and select best model
## ============================================================

e.mx.proc <- ENMevaluate(
  occs, envs, bg, 
  algorithm = "maxent",
  tune.args = list(fc = c("L"), rm = 1:5),
  partitions = "block",
  user.eval = proc
)

## Save evaluation results
write.csv(eval.results(e.mx), file = "Evaluation_resultss.csv")

## Select optimal model based on AICc and other metrics
res <- eval.results(e.mx)
write.csv(res, file = "D:/SDM_Shorea_robusta/Model_Evaluation_Result.csv", row.names = TRUE)

opt.aicc <- res %>% filter(delta.AICc == 0)
opt.seq <- res %>% filter(or.10p.avg == min(or.10p.avg)) %>% filter(auc.val.avg == max(auc.val.avg))
mod.seq <- eval.models(e.mx)[[opt.seq$tune.args]]

## Plot marginal response curves
plot(mod.seq, type = "cloglog")
jpeg("D:/SDM_Shorea_robusta/Cloglog_Response_Curves.jpeg", width = 8, height = 6, units = "in", res = 600)
plot(mod.seq, type = "cloglog", main = "Marginal Response Curves")
dev.off()

## ============================================================
## Save model predictions as raster maps
## ============================================================

tune_name <- as.character(opt.seq$tune.args)
pred.seq <- eval.predictions(e.mx)[[tune_name]]

## Save map with occurrence points
jpeg("D:/SDM_Shorea_robusta/Model_Prediction_With_Points.jpeg",
     width = 8, height = 6, units = "in", res = 600)
plot(pred.seq, zlim = c(0, 1), xlab = "Longitude", ylab = "Latitude", main = "Prediction with Occurrence Points")
points(eval.occs(e.mx), pch = 21, bg = eval.occs.grp(e.mx), cex = 0.7)
dev.off()

## Save map without occurrence points
jpeg("D:/SDM_Shorea_robusta/Best_Model_Prediction_Only.jpeg",
     width = 8, height = 6, units = "in", res = 600)
plot(pred.seq, zlim = c(0, 1), xlab = "Longitude", ylab = "Latitude", main = "Prediction Only")
dev.off()

## ============================================================
## Create map with countries, north arrow, scale, and color palette
## ============================================================

library(raster)
library(ggplot2)
library(sf)
library(ggspatial)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

pred.raster <- raster(pred.seq)
pred_df <- as.data.frame(pred.raster, xy = TRUE)
colnames(pred_df) <- c("lon", "lat", "suitability")
pred_df$lon <- as.numeric(pred_df$lon)
pred_df$lat <- as.numeric(pred_df$lat)
pred_df$suitability <- as.numeric(pred_df$suitability)

countries <- ne_countries(scale = "medium", returnclass = "sf")
centroids <- st_centroid(countries)
centroid_coords <- st_coordinates(centroids)
country_labels <- data.frame(lon = as.numeric(centroid_coords[, 1]),
                             lat = as.numeric(centroid_coords[, 2]),
                             name = countries$name)
country_labels <- subset(country_labels, lon >= 90 & lon <= 150 & lat >= 3 & lat <= 45)

cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
p <- ggplot() +
  geom_raster(data = pred_df, aes(x = lon, y = lat, fill = suitability)) +
  scale_fill_gradientn(colors = cols, name = "Suitability", limits = c(0, 1)) +
  geom_sf(data = countries, fill = NA, color = "black", size = 0.4) +
  geom_text_repel(data = country_labels, aes(x = lon, y = lat, label = name),
                  size = 3, color = "black", fontface = "bold") +
  labs(title = "Current", x = "Longitude", y = "Latitude") +
  scale_x_continuous(breaks = seq(90, 150, by = 15), labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(breaks = seq(5, 45, by = 10), labels = function(x) paste0(x, "°N")) +
  coord_sf(xlim = c(90, 150), ylim = c(3, 45), expand = FALSE, crs = st_crs(countries)) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(0.7, "cm"),
                         style = north_arrow_fancy_orienteering()) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = c(0.92, 0.6),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 11, face = "bold"))

p

## Save outputs
output_dir <- "D:/SDM_Shorea_robusta/Current_Prediction"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

jpeg_path <- file.path(output_dir, "Current_Best_Model_Prediction_Only_With_Countries.jpeg")
pdf_path  <- file.path(output_dir, "Current_Best_Model_Prediction_Only_With_Countries.pdf")

ggsave(jpeg_path, plot = p, width = 8, height = 6, units = "in", dpi = 600)
ggsave(pdf_path,  plot = p, width = 8, height = 6, units = "in", dpi = 600)
