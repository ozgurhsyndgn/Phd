## Update R & Packages ##

if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
}
updateR()

packs = as.data.frame(installed.packages(.libPaths()[1]), stringsAsFactors = F)

install.packages(packs$Package)


##biomod2
# setup environment ----

setwd("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future")

## install the latest release of biomod2

devtools::install_github('biomodhub/biomod2')

## load the required packages

library(devtools)
library(qpcR)
library(spThin)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(latticeExtra)
library(lattice)
library(corrplot)
library(RColorBrewer)
library(zoom)
library(ggpubr)
library(ggtext)
library(biomod2)
library(ggplot2)
library(viridis)
library(gridExtra)
library(rasterVis)
library(sf)
library(usdm)
library(tidyterra)
library(biomod2)
library(blockCV)
library(remotes)
library(terra)
library(ggcorrplot)

update.packages(checkbuilt=T)
remotes::install_github("rvalavi/blockCV")
pkgbuild::check_build_tools(debug = TRUE)


# Preprocessing -----------------------------------------------------------


castanea_occ <- read.csv('E:/Doktora/0_Thesis/2_Thesis_Data/0_Phd_Server/Working_Directory/DATA/Castanea_Locations_for_R/Centroids/castanea_centroids_spthin.csv')
castanea_shp <- st_read("E:/Doktora/0_Thesis/2_Thesis_Data/0_Phd_Server/Working_Directory/DATA/Castanea_Locations_for_R/Centroids/castanea_centroids_spthin.shp")
castanea_occ
#summary(castanea_occ)

# use an existing shapefile to restrict stack

extshp <- st_read("E:/Doktora/0_Thesis/2_Thesis_Data/TR_Shapefile/TUR_adm0.shp")
extshp_buffer <- st_read("E:/Doktora/0_Thesis/2_Thesis_Data/0_Phd_Server/Working_Directory/DATA/Castanea_Locations_for_R/Centroids/castanea_centroids_buffer.shp")


### Climatic Variables

worldclim <- rast(list.files(path = "E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/bio_tr/tr_split", full.names = TRUE))

worldclimcrop <- crop(worldclim,extshp,mask=T)

names(worldclimcrop)



## SUBSET ###

worldclimcrop_buffer <- crop(worldclimcrop ,extshp_buffer, mask=TRUE)

worldclimcrop_buffer_subset <- subset(worldclimcrop_buffer, c(14,2,4,7))

names(worldclimcrop_buffer_subset)

worldclimcrop_subset <- subset(worldclimcrop, c(14,2,4,7))

names(worldclimcrop_subset)

plot(worldclimcrop_subset)


# VIF / Plots -------------------------------------------------------------


vif_castanea<- vifstep(worldclimcrop_subset)
vif_castanea


vif <- vif_castanea@results
vif
corMatrix<-vif_castanea@corMatrix
corMatrix

worldclimcrop_excluded_buffer <- exclude(worldclimcrop_buffer,vif_castanea)

worldclimcrop_excluded <- exclude(worldclimcrop, vif_castanea)

plot(worldclimcrop_excluded_buffer)

####
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=FALSE)

levelplot(corMatrix,
          main = "Korelasyon Matrisi",
          xlab="",
          ylab="",
          col.regions = colorRampPalette(brewer.pal(4, "Blues"))(25),
          row.values = seq_len(nrow(corMatrix)),
          column.values = seq_len(ncol(corMatrix)),
          scales = list(x = list(rot = 26)),
          par.settings = list(axis.text = list(fontfamily = "Times New Roman", fontface="bold")),
          panel = function(...){
            panel.levelplot(...)
            coords <- expand.grid(x = 1:ncol(corMatrix), y = 1:nrow(corMatrix))
            panel.text(x = coords$x, y = coords$y,
                       labels = round(corMatrix, 2),
                       fontfamily = "Times New Roman", fontface="bold",
                       col = "black")
          })
ggcorrplot(corMatrix, hc.order = T, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_linedraw,
           colors = c("white",'#728FCE','#151B54'),
           lab = TRUE, lab_size = 6,lab_col='white', show.diag=F,
           legend.title = "", tl.cex=14)+
  theme(
    legend.key.width=unit(0.5, "cm"),
    legend.key.height=unit(1.5, "cm"),
    legend.title=element_text(size=40),
    legend.text=element_text(size=15)
  )



wc21_vif_data <- as.data.frame(vif_castanea@results)


ggplot(wc21_vif_data, aes(Variables, VIF, fill = VIF, label=round(VIF, digits = 2))) +
  scale_fill_distiller(palette="Blues",direction=1,type="div")+
  geom_col(aes(group = VIF),
           fill="#5f76b4",
           show.legend = FALSE,
           position = position_dodge2(), width = .4)+
  geom_text(aes(group = VIF),
            color = "white",
            position = position_stack(vjust = 0.5),
            size=5,
            fontface="bold",
            show.legend = FALSE)+
  theme_linedraw( base_size = 18)+
  labs(x = "DegiEkenler",
       y = "Varyans EiEme Degeri (VIF)")





#### PREDICTORS ##

worldclim_df <- as.data.frame(worldclimcrop_excluded, xy = TRUE, na.rm = TRUE)
worldclim_long_df <- 
  worldclim_df %>%
  pivot_longer(
    c(-x, -y),
    names_to = "variable",
    values_to = "value")
worldclim_long_df$variable<- factor(worldclim_long_df$variable,
                                    levels = c("bio4","bio10",
                                               "bio12","bio15"))
str(worldclim_long_df)

jpeg("2_predictors.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)



wc21_predictor_names <- c('bio4' = "bio4(sıcaklık mevsimselliği)",
                             'bio10' = "bio10(en sıcak çeyreğin ortalama sıcaklığı)",
                            'bio12' = "bio12(yıllık yağış)",
                              'bio15' = "bio15 (yağış mevsimselliği)")

plot1 <- ggarrange(nrow = 2, 
                   ncol = 2,
                   plotlist = lapply(split(worldclim_long_df, worldclim_long_df$variable),
                                     function(x){
                                       ggplot() + 
                                         geom_raster(data = x, 
                                                     aes(x =x, y = y, fill = value)) +
                                         scale_fill_distiller(palette= "Spectral") +
                                         facet_wrap(variable ~ .) +
                                         xlab("") +
                                         ylab("") +
                                         theme_bw() +
                                         theme(legend.title = element_blank()) +
                                         theme(strip.text.x = element_text(size = 12, face = "bold")) +
                                         theme(plot.title = element_text(hjust = 0.5))+
                                         theme(text = element_text(family = "Times New Roman"))
                                     }))

annotate_figure(plot1,
                bottom = text_grob("", face = "bold", 
                                   size = 16, family = "Times New Roman"),
                left = text_grob("", rot = 90, face = "bold", 
                                 size = 16, family = "Times New Roman"))

dev.off()





# BIOMOD_FormatingData ----------------------------------------------------

set.seed(42)
castanea_data <- 
  BIOMOD_FormatingData(
    resp.var = castanea_occ['species'],
    resp.xy = castanea_occ[, c('x', 'y')],
    expl.var = worldclimcrop_buffer_subset,
    resp.name = "castanea",
    PA.nb.rep = 10, PA.nb.absences = c(400,400,400,400,800,800,800,800,1200,1200),
    PA.strategy = "disk",
    PA.dist.min = 10000,
    na.rm = TRUE,
    filter.raster = TRUE
  )



castanea_data



# blockCV -----------------------------------------------------------------


castanea_dataPA1 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA1<-castanea_dataPA1|>filter(castanea_dataPA1$PA1==TRUE)|>dplyr::select(ID, PA1,x,y)

castanea_points_extractedPA1 <- sf::st_as_sf(castanea_dataPA1, coords = c("x", "y"), crs = 4326) ## for BLOCK-CV

castanea_points_extractedPA1

scv1 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA1,
  column = "PA1", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA2 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA2<-castanea_dataPA2|>filter(castanea_dataPA2$PA2==TRUE)|>dplyr::select(ID, PA2,x,y)

castanea_points_extractedPA2 <- sf::st_as_sf(castanea_dataPA2, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA2

scv2 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA2,
  column = "PA2", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA3 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA3<-castanea_dataPA3|>filter(castanea_dataPA3$PA3==TRUE)|>dplyr::select(ID,PA3,x,y)

castanea_points_extractedPA3 <- sf::st_as_sf(castanea_dataPA3, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA3

scv3 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA3,
  column = "PA3", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA4 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA4<-castanea_dataPA4|>filter(castanea_dataPA4$PA4==TRUE)|>dplyr::select(ID,PA4,x,y)

castanea_points_extractedPA4 <- sf::st_as_sf(castanea_dataPA4, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA4

scv4 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA4,
  column = "PA4", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA5 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA5<-castanea_dataPA5|>filter(castanea_dataPA5$PA5==TRUE)|>dplyr::select(ID,PA5,x,y)

castanea_points_extractedPA5 <- sf::st_as_sf(castanea_dataPA5, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA5

scv5 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA5,
  column = "PA5", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA6 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA6<-castanea_dataPA6|>filter(castanea_dataPA6$PA6==TRUE)|>dplyr::select(ID,PA6,x,y)

castanea_points_extractedPA6 <- sf::st_as_sf(castanea_dataPA6, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA6

scv6 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA6,
  column = "PA6", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA7 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA7<-castanea_dataPA7|>filter(castanea_dataPA7$PA7==TRUE)|>dplyr::select(ID,PA7,x,y)

castanea_points_extractedPA7 <- sf::st_as_sf(castanea_dataPA7, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA7

scv7 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA7,
  column = "PA7", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA8 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA8<-castanea_dataPA8|>filter(castanea_dataPA8$PA8==TRUE)|>dplyr::select(ID,PA8,x,y)

castanea_points_extractedPA8 <- sf::st_as_sf(castanea_dataPA8, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA8

scv8 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA8,
  column = "PA8", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA9 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA9<-castanea_dataPA9|>filter(castanea_dataPA9$PA9==TRUE)|>dplyr::select(ID,PA9,x,y)

castanea_points_extractedPA9 <- sf::st_as_sf(castanea_dataPA9, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA9

scv9 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA9,
  column = "PA9", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

castanea_dataPA10 <- castanea_data@PA.table |> tibble::rownames_to_column( "ID") |> cbind(castanea_data@coord)

castanea_dataPA10<-castanea_dataPA10|>filter(castanea_dataPA10$PA10==TRUE)|>dplyr::select(ID,PA10,x,y)

castanea_points_extractedPA10 <- sf::st_as_sf(castanea_dataPA10, coords = c("x", "y"), crs = 4326)

castanea_points_extractedPA10

scv10 <- blockCV::cv_spatial(
  x = castanea_points_extractedPA10,
  column = "PA10", # the response column (binary or multi-class)
  r = worldclimcrop_buffer_subset,
  k = 4, # number of folds
  hexagon = FALSE,
  rows_cols = c(4,9),
  #size = 150000.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  #iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
)

##

## SPATIAL AUTOCORELATION ##

castanea_range <- cv_spatial_autocor(
  x = castanea_shp, # species data
  column = "species", # column storing presence-absence records (0s and 1s)
  plot = FALSE
)

castanea_range$range

castanea_DataSplitTable <- scv1$biomod_table


nK <- 4
nPseudo <- 1

CVTablePA1 <- matrix(NA, nrow = nrow(scv1$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA1) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA1[,c(i*nK-(nK-1)):c(i*nK)] <- scv1$biomod_table
  colnames(CVTablePA1)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA1 <-as.data.frame(cbind(CVTablePA1,castanea_dataPA1$ID))
head(CVTablePA1)


###

CVTablePA2 <- matrix(NA, nrow = nrow(scv2$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA2) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA2[,c(i*nK-(nK-1)):c(i*nK)] <- scv2$biomod_table
  colnames(CVTablePA2)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA2 <-as.data.frame(cbind(CVTablePA2,castanea_dataPA2$ID))
head(CVTablePA2)


###

CVTablePA3 <- matrix(NA, nrow = nrow(scv3$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA3) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA3[,c(i*nK-(nK-1)):c(i*nK)] <- scv3$biomod_table
  colnames(CVTablePA3)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA3 <-as.data.frame(cbind(CVTablePA3,castanea_dataPA3$ID))
head(CVTablePA3)


###

CVTablePA4 <- matrix(NA, nrow = nrow(scv4$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA4) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA4[,c(i*nK-(nK-1)):c(i*nK)] <- scv4$biomod_table
  colnames(CVTablePA4)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA4 <-as.data.frame(cbind(CVTablePA4,castanea_dataPA4$ID))
head(CVTablePA4)

###

CVTablePA5 <- matrix(NA, nrow = nrow(scv5$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA5) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA5[,c(i*nK-(nK-1)):c(i*nK)] <- scv5$biomod_table
  colnames(CVTablePA5)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA5 <-as.data.frame(cbind(CVTablePA5,castanea_dataPA5$ID))
head(CVTablePA5)

###

CVTablePA6 <- matrix(NA, nrow = nrow(scv6$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA6) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA6[,c(i*nK-(nK-1)):c(i*nK)] <- scv6$biomod_table
  colnames(CVTablePA6)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA6 <-as.data.frame(cbind(CVTablePA6,castanea_dataPA6$ID))
head(CVTablePA6)


###

CVTablePA7 <- matrix(NA, nrow = nrow(scv7$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA7) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA7[,c(i*nK-(nK-1)):c(i*nK)] <- scv7$biomod_table
  colnames(CVTablePA7)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA7 <-as.data.frame(cbind(CVTablePA7,castanea_dataPA7$ID))
head(CVTablePA7)

###

CVTablePA8 <- matrix(NA, nrow = nrow(scv8$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA8) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA8[,c(i*nK-(nK-1)):c(i*nK)] <- scv8$biomod_table
  colnames(CVTablePA8)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA8 <-as.data.frame(cbind(CVTablePA8,castanea_dataPA8$ID))
head(CVTablePA8)

###

CVTablePA9 <- matrix(NA, nrow = nrow(scv9$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA9) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA9[,c(i*nK-(nK-1)):c(i*nK)] <- scv9$biomod_table
  colnames(CVTablePA9)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA9 <-as.data.frame(cbind(CVTablePA9,castanea_dataPA9$ID))
head(CVTablePA9)

###

CVTablePA10 <- matrix(NA, nrow = nrow(scv10$biomod_table), ncol = nK*nPseudo)
colnames(CVTablePA10) <- paste0("x", seq_len(nK*nPseudo))
for(i in 1:nPseudo){
  CVTablePA10[,c(i*nK-(nK-1)):c(i*nK)] <- scv10$biomod_table
  colnames(CVTablePA10)[c(i*nK-(nK-1)):c(i*nK)] <- paste0("_PA",i, "_RUN", seq_len(nK))
}
CVTablePA10 <-as.data.frame(cbind(CVTablePA10,castanea_dataPA10$ID))
head(CVTablePA10)

###


colnames(CVTablePA1) <- c("_PA1_RUN1","_PA1_RUN2","_PA1_RUN3","_PA1_RUN4","ID")
colnames(CVTablePA2) <- c("_PA2_RUN1","_PA2_RUN2","_PA2_RUN3","_PA2_RUN4","ID")
colnames(CVTablePA3) <- c("_PA3_RUN1","_PA3_RUN2","_PA3_RUN3","_PA3_RUN4","ID")
colnames(CVTablePA4) <- c("_PA4_RUN1","_PA4_RUN2","_PA4_RUN3","_PA4_RUN4","ID")
colnames(CVTablePA5) <- c("_PA5_RUN1","_PA5_RUN2","_PA5_RUN3","_PA5_RUN4","ID")
colnames(CVTablePA6) <- c("_PA6_RUN1","_PA6_RUN2","_PA6_RUN3","_PA6_RUN4","ID")
colnames(CVTablePA7) <- c("_PA7_RUN1","_PA7_RUN2","_PA7_RUN3","_PA7_RUN4","ID")
colnames(CVTablePA8) <- c("_PA8_RUN1","_PA8_RUN2","_PA8_RUN3","_PA8_RUN4","ID")
colnames(CVTablePA9) <- c("_PA9_RUN1","_PA9_RUN2","_PA9_RUN3","_PA9_RUN4","ID")
colnames(CVTablePA10) <- c("_PA10_RUN1","_PA10_RUN2","_PA10_RUN3","_PA10_RUN4","ID")


CVPA1_2<-CVTablePA1|>full_join(CVTablePA2,by=c("ID"))
CVPA3_4<-CVTablePA3|>full_join(CVTablePA4,by=c("ID"))
CVPA5_6<-CVTablePA5|>full_join(CVTablePA6,by=c("ID"))
CVPA7_8<-CVTablePA7|>full_join(CVTablePA8,by=c("ID"))
CVPA9_10<-CVTablePA9|>full_join(CVTablePA10,by=c("ID"))


CVPA1_2_3_4<-CVPA1_2|>full_join(CVPA3_4,by=c("ID"))
CVPA5_6_7_8<-CVPA5_6|>full_join(CVPA7_8,by=c("ID"))
CVPA5_6_7_8_9_10<- CVPA5_6_7_8|>full_join(CVPA9_10,by=c("ID"))

CVPA1_2_3_4_5_6_7_8_9_10<-CVPA1_2_3_4|>full_join(CVPA5_6_7_8_9_10,by=c("ID"))
CVPA1_2_3_4_5_6_7_8_9_10_s <- CVPA1_2_3_4_5_6_7_8_9_10[with(CVPA1_2_3_4_5_6_7_8_9_10, order(as.numeric(ID))),]
CVPA1_2_3_4_5_6_7_8_9_10_s_r <- CVPA1_2_3_4_5_6_7_8_9_10_s %>% remove_rownames %>% column_to_rownames(var="ID")
CVPA1_2_3_4_5_6_7_8_9_10_s_r <- dplyr::mutate_all(CVPA1_2_3_4_5_6_7_8_9_10_s_r, as.logical)

CVPA1_2_3_4_5_6_7_8_9_10_s_r

CVPA1_2_3_4_5_6_7_8_9_10_s_r_m <- data.matrix(CVPA1_2_3_4_5_6_7_8_9_10_s_r)



# BIOMOD_Tuning -----------------------------------------------------------


## formatted object summary
castanea_data

## plot of selected pseudo-absences

plot(castanea_data)

library(doParallel)
cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
time.seq <- system.time(opt.t <- bm_ModelingOptions(data.type = 'binary',
                            bm.format = castanea_data,
                            models = c('RF','GBM','GLM'),
                            strategy = 'default',
                            calib.lines=CVPA1_2_3_4_5_6_7_8_9_10_s_r_m))
stopCluster(cl)

opt.t

ModelsTable

####

unregister_dopar <- function() {
  tuned.rf <- foreach:::.foreachGlobals
  rm(list=ls(name=tuned.rf), pos=tuned.rf)
}
unregister_dopar()

library(doParallel)
gl <- makeCluster(6)
doParallel::registerDoParallel(gl)
time.seq <- system.time(tuned.rf <- bm_Tuning(model = 'RF',
                        tuning.fun = 'rf', ## see in ModelsTable
                        do.formula = TRUE,
                        bm.options = opt.t@options$RF.binary.randomForest.randomForest,
                        bm.format = castanea_data,
                        calib.lines = CVPA1_2_3_4_5_6_7_8_9_10_s_r_m))
stopCluster(gl)

tuned.rf

###

unregister_dopar <- function() {
  tuned.gbm <- foreach:::.foreachGlobals
  rm(list=ls(name=tuned.gbm), pos=tuned.gbm)
}
unregister_dopar()


library(doParallel)
bl <- makeCluster(6)
doParallel::registerDoParallel(bl)
time.seq <- system.time(tuned.gbm <- bm_Tuning(model = 'GBM',
                                               tuning.fun = 'gbm', ## see in ModelsTable
                                               do.formula = TRUE,
                                               bm.options = opt.t@options$GBM.binary.gbm.gbm,
                                               bm.format = castanea_data,
                                               calib.lines = CVPA1_2_3_4_5_6_7_8_9_10_s_r_m))
stopCluster(bl)

###

unregister_dopar <- function() {
  tuned.mars <- foreach:::.foreachGlobals
  rm(list=ls(name=tuned.mars), pos=tuned.mars)
}
unregister_dopar()

library(doParallel)
el <- makeCluster(6)
doParallel::registerDoParallel(el)
time.seq <- system.time(tuned.glm <- bm_Tuning(model = 'GLM',
                                               tuning.fun = 'glm', ## see in ModelsTable
                                               do.formula = TRUE,
                                               bm.options = opt.t@options$GLM.binary.stats.glm,
                                               bm.format = castanea_data,
                                               calib.lines = CVPA1_2_3_4_5_6_7_8_9_10_s_r_m))

stopCluster(el)



user.val <- list(RF.binary.randomForest.randomForest = tuned.rf,
                 GBM.binary.gbm.gbm = tuned.gbm,
                 GLM.binary.stats.glm = tuned.glm)
user.val


#opt.d <- bm_ModelingOptions(data.type = 'binary',
 #                           models = c("GBM","RF"),
  #                          bm.format = castanea_data,
   #                         strategy = 'bigboss',
    #                        user.base = "bigboss")
#opt.d

###

unregister_dopar <- function() {
  myOptions <- foreach:::.foreachGlobals
  rm(list=ls(name=myOptions), pos=myOptions)
}
unregister_dopar()


library(doParallel)
dl <- makeCluster(6)
doParallel::registerDoParallel(dl)
time.seq <- system.time(myOptions <- bm_ModelingOptions(data.type = 'binary',
                                                        bm.format = castanea_data,
                                                        models = c("RF","GBM","GLM"),
                                                        strategy = 'tuned',
                                                        user.val = user.val,
                                                        calib.lines = CVPA1_2_3_4_5_6_7_8_9_10_s_r_m))
stopCluster(dl)


myOptions



# BIOMOD_Modeling ---------------------------------------------------------

#### Single models

### Merge the parallelizations ###

unregister_dopar <- function() {
  myBiomodModelOut <- foreach:::.foreachGlobals
  rm(list=ls(name=myBiomodModelOut), pos=myBiomodModelOut)
}
unregister_dopar()


# Model single models

myBiomodModelOut <- BIOMOD_Modeling(bm.format = castanea_data,
                                    modeling.id = 'AllModels',
                                    models = c("RF","GBM","GLM"),
                                    bm.options = myOptions,
                                    CV.strategy = "user.defined",
                                    CV.user.table = CVPA1_2_3_4_5_6_7_8_9_10_s_r,
                                    var.import = 2,
                                    metric.eval = c('TSS','ROC',"BOYCE"),
                                    #prevalence = 0.5,
                                    scale.models=FALSE,
                                    CV.do.full.models = FALSE,
                                    do.progress= TRUE,
                                    seed.val = 42,
                                    models.pa = list("GBM" = c("PA1","PA2","PA3","PA4","PA5","PA6","PA7","PA8","PA9","PA10"),
                                                     "RF" = c("PA1","PA2","PA3","PA4","PA5","PA6","PA7","PA8","PA9","PA10"),
                                                     "GLM" = c("PA1","PA2","PA3","PA4","PA5","PA6","PA7","PA8","PA9","PA10")))
#models.pa = list("GBM" = c("PA1","PA2","PA3","PA4"),
 #                "RF" = c("PA1","PA2","PA3","PA4"),
  #               "MARS" = c("PA1","PA2","PA3","PA4")))

myBiomodModelOut


# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)
write.csv(myBiomodModelOut@models.evaluation@val,"evaluation_scores.csv", row.names=TRUE)


x11()
# Represent evaluation scores & variables importance
jpeg("9_validation_scores.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
dev.off()
jpeg("10_validation_box.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'PA'), dataset = 'validation')
dev.off()
jpeg("11_exp_variables_run.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'), dataset = 'validation')
dev.off()
jpeg("12_exp_variables_box.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)


##### WC21 Validation Scores ####

bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'), dataset = 'validation')


wc21_validation_plots <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, 
                                            group.by = c('algo', 'algo'), dataset = 'validation')

wc21_validation_plots_df <- as.data.frame(wc21_validation_plots$tab)

write.csv(wc21_validation_plots_df, "wc21_validation_plots_df.csv")


wc21_validation_plots_df_plot <- ggplot(wc21_validation_plots_df, aes(x = metric.eval, y = validation, fill = metric.eval)) + 
                                  geom_boxplot(width=0.35) +
                                  scale_fill_brewer(palette="Dark2")+
                                  labs(x="", y="validasyon skoru")+
                                  theme_bw(base_size= 16)+
                                  theme(strip.text.x = element_text(size = 16, colour = "black",face="bold"),
                                        legend.position = "bottom",legend.text = element_text(size=16))+
                                  guides(fill=guide_legend("Model Degerlendirme Yontemleri"))+
                                  #scale_fill_manual(values=c("#2CB63B", "orange", "#E4E71C"))+
                                  facet_grid(. ~ algo)






#### WC21 Expl. Variables ####

wc21_exp_var_plots<- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
                                      

wc21_exp_var_plots_df <- as.data.frame(wc21_exp_var_plots$tab)

write.csv(wc21_exp_var_plots_df, "wc21_exp_var_plots_df.csv")


wc21_exp_var_plots_df_plot <- ggplot(wc21_exp_var_plots_df, aes(x = expl.var, y = var.imp, fill = expl.var)) + 
                                geom_boxplot(width=0.35) +
                                scale_fill_brewer(palette="Dark2")+
                                labs(x="", y="deDiEken onemi")+
                                theme_bw(base_size= 16)+
                                theme(strip.text.x = element_text(size = 16, colour = "black",face="bold"),
                                      legend.position = "bottom",legend.text = element_text(size=16))+
                                #scale_fill_manual(values=c("#2CB63B", "orange", "#E4E71C"))+
                                guides(fill=guide_legend("DeDiEkenler"))+
                                facet_grid(. ~ algo)


ggarrange(wc21_validation_plots_df_plot,wc21_exp_var_plots_df_plot, labels = c("(a)","(b)"),font.label = list(size = 22),
          ncol= 2, nrow= 1, common.legend = F, xscale)



# Represent response curves

jpeg("12_response_curves_all_runs.jpeg",
     units = "px", width = 2500, height = 1500, res = 300)
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = "all",
                      fixed.var = 'median')
dev.off()



# BIOMOD_Projections ------------------------------------------------------

### Project models

#### Single models


# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  #nb.cpu = 4,
                                  proj.name = 'Current',
                                  new.env = worldclimcrop_subset,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = FALSE,
                                  seed.val = 42)


#output.format = '.tif',


jpeg("13_projection.jpeg",
     units = "px", width = 3200, height = 3200, res= 300)
plot(myBiomodProj)
dev.off()


# BIOMOD_EnsembleModeling / EnsembleForecasting -------------------------------------------------


myBiomodEM_algo <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                           models.chosen = 'all',
                                           em.by = 'algo',
                                           metric.select = c('TSS'),
                                           metric.select.thresh = c(0.75),
                                           metric.select.dataset = 'validation',
                                           var.import = 2,
                                           metric.eval = c('TSS', 'ROC','BOYCE'),
                                           em.algo = c('EMmean'),
                                           EMwmean.decay = 'proportional')

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.75),
                                      metric.select.dataset = 'validation',
                                      var.import = 2,
                                      metric.eval = c('TSS', 'ROC','BOYCE'),
                                      em.algo = c('EMmean'),
                                      EMwmean.decay = 'proportional')


myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

write.csv(myBiomodEM@models.evaluation@val,"worldclim_ensemble_evaluation_scores_4lyrs.csv", row.names=TRUE)
write.csv(myBiomodEM@variables.importance@val,"worldclim_variable_importance_4lyrs.csv", row.names=TRUE)

# Represent evaluation scores & variables importance

jpeg("17_variable_importance_box.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)
bm_PlotVarImpBoxplot(bm.out = myBiomodEM_algo, group.by = c('algo', 'expl.var', 'full.name'))
dev.off()

# Represent response curves
jpeg("18_response_curves_varimp4.jpeg",
     units = "px", width = 2000, height = 2000, res = 300)
bm_PlotResponseCurves(bm.out = myBiomodEM_algo, 
                      models.chosen = get_built_models(myBiomodEM_algo),
                      fixed.var = 'median')
dev.off()


#Try to manipulate s4 and enhance response curves

myBiomodEM_algo_curves <- bm_PlotResponseCurves(bm.out = myBiomodEM_algo, 
                                                models.chosen = get_built_models(myBiomodEM_algo),
                                                fixed.var = 'median')

myBiomodEM_algo_curves_df <- myBiomodEM_algo_curves$tab


write.csv(myBiomodEM_algo_curves_df, "myBiomodEM_algo_curves_df.csv")

myBiomodEM_algo_curves_df <- read.csv("myBiomodEM_algo_curves_df.csv")


myBiomodEM_algo_curves_df$pred.name <- as.factor(myBiomodEM_algo_curves_df$pred.name)
myBiomodEM_algo_curves_df$expl.name <- factor(myBiomodEM_algo_curves_df$expl.name,
                                              levels = c("bio4",
                                                         "bio10",
                                                         "bio12",
                                                         "bio15"))


myBiomodEM_algo_curves_df_plot<- ggplot(data = myBiomodEM_algo_curves_df, aes(x = expl.val, 
                                             y = pred.val,
                                             group = pred.name,
                                             color = pred.name)) +
                                      geom_line(linewidth=0.65, linetype = 1) +
                                      facet_wrap(.~expl.name, 
                                      scales="free_x", 
                                      ncol = 2) +
                                      scale_colour_discrete(
                                      limits = c("castanea_EMmeanByTSS_mergedData_mergedRun_GBM",
                                                  "castanea_EMmeanByTSS_mergedData_mergedRun_GLM",
                                                  "castanea_EMmeanByTSS_mergedData_mergedRun_RF"),
                                                  labels = c("GBM", "GLM", "RF")) +
                                      theme_bw() +
                                      xlab("") +
                                      ylab("") +
                                      theme(legend.position = "bottom") +
                                      theme(legend.text = element_text())+
                                      theme(legend.title=element_blank())

#Try to manipulate s4 and enhance response curves -- ENSEMBLE

myBiomodEM_curves <- bm_PlotResponseCurves(bm.out = myBiomodEM, 
                                                models.chosen = get_built_models(myBiomodEM),
                                                fixed.var = 'median')

myBiomodEM_curves_df <- myBiomodEM_curves$tab


write.csv(myBiomodEM_curves_df, "myBiomodEM1_curves_df.csv")

myBiomodEM_curves_df <- read.csv("myBiomodEM1_curves_df.csv")


myBiomodEM_curves_df$pred.name <- as.factor(myBiomodEM_curves_df$pred.name)
myBiomodEM_curves_df$expl.name <- factor(myBiomodEM_curves_df$expl.name,
                                              levels = c("bio4",
                                                         "bio10",
                                                         "bio12",
                                                         "bio15"))


myBiomodEM_curves_df_plot <- ggplot(data = myBiomodEM_curves_df, aes(x = expl.val, 
                                             y = pred.val,
                                             group = pred.name,
                                             color = pred.name)) +
                                    geom_line(linewidth=0.7, linetype = 1) +
                                    facet_wrap(.~expl.name, 
                                               scales="free_x", 
                                               ncol = 2) +
                                    scale_colour_discrete(
                                      limits = c("castanea_EMmeanByTSS_mergedData_mergedRun_mergedAlgo"),
                                      labels = c("Birle??ik Model (Ensemble)")) +
                                    theme_bw() +
                                    xlab("") +
                                    ylab("") +
                                    theme(legend.position = "bottom") +
                                    theme(legend.text = element_text())+
                                    theme(legend.title=element_blank())



################################# GGARRANGE - Curves #########

ggarrange(myBiomodEM_algo_curves_df_plot, myBiomodEM_curves_df_plot,
          ncol= 2, nrow= 1, common.legend = F)


#Project ensemble models (from single projections)
myBiomodEMProj_algo_buffer <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo,
                                                         #nb.cpu = 4,
                                                         proj.name = 'CurrentEM_buffer_algo',
                                                         new.env = worldclimcrop_buffer_subset,
                                                         models.chosen = 'all',
                                                         metric.binary = 'all',
                                                         metric.filter = 'all')

myBiomodEMProj_EM_buffer <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                       #nb.cpu = 4,
                                                       proj.name = 'CurrentEM_buffer',
                                                       new.env = worldclimcrop_buffer_subset,
                                                       models.chosen = 'all',
                                                       metric.binary = 'all',
                                                       metric.filter = 'all')


# Project ensemble models (building single projections)
myBiomodEMProj_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo,
                                                  #nb.cpu = 4,
                                                  proj.name = 'CurrentEM_algo',
                                                  new.env = worldclimcrop_subset,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')

# Project ensemble models (building single projections)
myBiomodEMProj_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                #nb.cpu = 4,
                                                proj.name = 'CurrentEM',
                                                new.env = worldclimcrop_subset,
                                                models.chosen = 'all',
                                                metric.binary = 'all',
                                                metric.filter = 'all')

##########################################################################################
                                            

myBiomodEMProj

plot(myBiomodEMProj_algo)
plot(myBiomodEMProj_EM)
plot(myBiomodEMProj_algo_buffer)
plot(myBiomodEMProj_EM_buffer)


# BIOMOD_SpeciesRangeChange (SRC) -----------------------------------------

display.brewer.all()

# wc21_IPSL_50_26 ---------------------------------------------------------

# Load environmental variables extracted from worldclim

wc21_IPSL_50_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/IPSL/26/wc21_26_ipsl_41_60.tif")
names(wc21_IPSL_50_261)
names(wc21_IPSL_50_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_IPSL_50_261 <-
  subset(wc21_IPSL_50_261, c(4,10,12,15))
wc21_IPSL_50_26 <- crop(wc21_IPSL_50_261, extshp,mask=T)

names(wc21_IPSL_50_26)

# Project onto future conditions
castanea_models_proj_2050_IPSL26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                      proj.name = '2050_IPSL26',
                                                      new.env = wc21_IPSL_50_26,
                                                      models.chosen = 'all',
                                                      metric.binary = 'TSS',
                                                      build.clamping.mask = FALSE)


castanea_models_proj_2050_IPSL26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                  bm.proj = castanea_models_proj_2050_IPSL26_tr,
                                                                  proj.name = '2050_IPSL26_EM',
                                                                  models.chosen = 'all',
                                                                  metric.binary = 'all',
                                                                  metric.filter = 'all')

castanea_models_proj_2050_IPSL26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                     bm.proj = castanea_models_proj_2050_IPSL26_tr,
                                                                     proj.name = '2050_IPSL26_algo',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')


plot(castanea_models_proj_2050_IPSL26_tr_EM)
plot(castanea_models_proj_2050_IPSL26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_IPSL26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_IPSL26_EM/proj_2050_IPSL26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_IPSL26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_IPSL26)

myBiomodRangeSize_2050_IPSL26$Compt.By.Models

plot(myBiomodRangeSize_2050_IPSL26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "white")


writeRaster(myBiomodRangeSize_2050_IPSL26$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_IPSL26.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_IPSL26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_IPSL_50_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_IPSL_50_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/IPSL/85/wc21_85_ipsl_41_60.grd")
names(wc21_IPSL_50_851)
names(wc21_IPSL_50_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_IPSL_50_851 <-
  subset(wc21_IPSL_50_851, c(4,10,12,15))
wc21_IPSL_50_85 <- crop(wc21_IPSL_50_851, extshp,mask=T)


names(wc21_IPSL_50_85)

# Project onto future conditions
castanea_models_proj_2050_IPSL85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2050_IPSL85',
                                                         new.env = wc21_IPSL_50_85,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2050_IPSL85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2050_IPSL85_tr,
                                                                     proj.name = '2050_IPSL85_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2050_IPSL85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2050_IPSL85_tr,
                                                                       proj.name = '2050_IPSL85_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2050_IPSL85_tr_EM)
plot(castanea_models_proj_2050_IPSL85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_IPSL85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_IPSL85_EM/proj_2050_IPSL85_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_IPSL85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_IPSL85)

myBiomodRangeSize_2050_IPSL85$Compt.By.Models

plot(myBiomodRangeSize_2050_IPSL85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2050_IPSL85$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_IPSL85.tif",overwrite=T)


# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_IPSL85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)



# wc21_IPSL_70_26 ---------------------------------------------------------

wc21_IPSL_70_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/IPSL/26/IPSL_26_61_80.tif")
names(wc21_IPSL_70_261)
names(wc21_IPSL_70_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_IPSL_70_261 <-
  subset(wc21_IPSL_70_261, c(4,10,12,15))
wc21_IPSL_70_26 <- crop(wc21_IPSL_70_261, extshp,mask=T)

names(wc21_IPSL_70_26)

# Project onto future conditions
castanea_models_proj_2070_IPSL26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2070_IPSL26',
                                                         new.env = wc21_IPSL_70_26,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2070_IPSL26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2070_IPSL26_tr,
                                                                     proj.name = '2070_IPSL26_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2070_IPSL26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2070_IPSL26_tr,
                                                                       proj.name = '2070_IPSL26_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2070_IPSL26_tr_EM)
plot(castanea_models_proj_2070_IPSL26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_IPSL26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_IPSL26_EM/proj_2070_IPSL26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2070_IPSL26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_IPSL26)

myBiomodRangeSize_2070_IPSL26$Compt.By.Models

plot(myBiomodRangeSize_2070_IPSL26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2070_IPSL26$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_IPSL26.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_IPSL26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_IPSL_70_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_IPSL_70_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/IPSL/85/wc21_85_ipsl_61_80.grd")
names(wc21_IPSL_70_851)
names(wc21_IPSL_70_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_IPSL_70_851 <-
  subset(wc21_IPSL_70_851, c(4,10,12,15))
wc21_IPSL_70_85 <- crop(wc21_IPSL_70_851, extshp,mask=T)


names(wc21_IPSL_70_85)

# Project onto future conditions
castanea_models_proj_2070_IPSL85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2070_IPSL85',
                                                         new.env = wc21_IPSL_70_85,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2070_IPSL85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2070_IPSL85_tr,
                                                                     proj.name = '2070_IPSL85_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2070_IPSL85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2070_IPSL85_tr,
                                                                       proj.name = '2070_IPSL85_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2070_IPSL85_tr_EM)
plot(castanea_models_proj_2070_IPSL85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_IPSL85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_IPSL85_EM/proj_2070_IPSL85_EM_castanea_ensemble_TSSbin.tif")


# Compute differences
myBiomodRangeSize_2070_IPSL85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_IPSL85)

myBiomodRangeSize_2070_IPSL85$Compt.By.Models

plot(myBiomodRangeSize_2070_IPSL85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2070_IPSL85$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_IPSL85.tif",overwrite=T)


# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_IPSL85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_MPI_50_26 ----------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_MPI_50_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/MPI/26/MPI_26_41_60.tif")
names(wc21_MPI_50_261)
names(wc21_MPI_50_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_MPI_50_261 <-
  subset(wc21_MPI_50_261, c(4,10,12,15))
wc21_MPI_50_26 <- crop(wc21_MPI_50_261, extshp,mask=T)

names(wc21_MPI_50_26)

# Project onto future conditions
castanea_models_proj_2050_MPI26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2050_MPI26',
                                                         new.env = wc21_MPI_50_26,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2050_MPI26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2050_MPI26_tr,
                                                                     proj.name = '2050_MPI26_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2050_MPI26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2050_MPI26_tr,
                                                                       proj.name = '2050_MPI26_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2050_MPI26_tr_EM)
plot(castanea_models_proj_2050_MPI26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_MPI26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_MPI26_EM/proj_2050_MPI26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_MPI26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_MPI26)

myBiomodRangeSize_2050_MPI26$Compt.By.Models

plot(myBiomodRangeSize_2050_MPI26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2050_MPI26$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_MPI26.tif",overwrite=T)


# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_MPI26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_MPI_50_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_MPI_50_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/MPI/85/MPI_85_41_60.tif")
names(wc21_MPI_50_851)
names(wc21_MPI_50_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_MPI_50_851 <-
  subset(wc21_MPI_50_851, c(4,10,12,15))
wc21_MPI_50_85 <- crop(wc21_MPI_50_851, extshp,mask=T)


names(wc21_MPI_50_85)

# Project onto future conditions
castanea_models_proj_2050_MPI85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2050_MPI85',
                                                         new.env = wc21_MPI_50_85,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2050_MPI85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2050_MPI85_tr,
                                                                     proj.name = '2050_MPI85_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2050_MPI85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2050_MPI85_tr,
                                                                       proj.name = '2050_MPI85_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2050_MPI85_tr_EM)
plot(castanea_models_proj_2050_MPI85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_MPI85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_MPI85_EM/proj_2050_MPI85_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_MPI85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_MPI85)

myBiomodRangeSize_2050_MPI85$Compt.By.Models

plot(myBiomodRangeSize_2050_MPI85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2050_MPI85$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_MPI85.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_MPI85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)



# wc21_MPI_70_26 ---------------------------------------------------------

wc21_MPI_70_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/MPI/26/MPI_26_61_80.tif")
names(wc21_MPI_70_261)
names(wc21_MPI_70_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_MPI_70_261 <-
  subset(wc21_MPI_70_261, c(4,10,12,15))
wc21_MPI_70_26 <- crop(wc21_MPI_70_261, extshp,mask=T)

names(wc21_MPI_70_26)

# Project onto future conditions
castanea_models_proj_2070_MPI26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2070_MPI26',
                                                         new.env = wc21_MPI_70_26,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2070_MPI26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2070_MPI26_tr,
                                                                     proj.name = '2070_MPI26_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2070_MPI26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2070_MPI26_tr,
                                                                       proj.name = '2070_MPI26_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2070_MPI26_tr_EM)
plot(castanea_models_proj_2070_MPI26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_MPI26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_MPI26_EM/proj_2070_MPI26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2070_MPI26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_MPI26)

myBiomodRangeSize_2070_MPI26$Compt.By.Models

plot(myBiomodRangeSize_2070_MPI26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2070_MPI26$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_MPI26.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_MPI26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_MPI_70_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_MPI_70_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/MPI/85/wc21_85_mpi_61_80.grd")
names(wc21_MPI_70_851)
names(wc21_MPI_70_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                             'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_MPI_70_851 <-
  subset(wc21_MPI_70_851, c(4,10,12,15))
wc21_MPI_70_85 <- crop(wc21_MPI_70_851, extshp,mask=T)


names(wc21_MPI_70_85)

# Project onto future conditions
castanea_models_proj_2070_MPI85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                         proj.name = '2070_MPI85',
                                                         new.env = wc21_MPI_70_85,
                                                         models.chosen = 'all',
                                                         metric.binary = 'TSS',
                                                         build.clamping.mask = FALSE)


castanea_models_proj_2070_MPI85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                     bm.proj = castanea_models_proj_2070_MPI85_tr,
                                                                     proj.name = '2070_MPI85_EM',
                                                                     models.chosen = 'all',
                                                                     metric.binary = 'all',
                                                                     metric.filter = 'all')

castanea_models_proj_2070_MPI85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                       bm.proj = castanea_models_proj_2070_MPI85_tr,
                                                                       proj.name = '2070_MPI85_algo',
                                                                       models.chosen = 'all',
                                                                       metric.binary = 'all',
                                                                       metric.filter = 'all')


plot(castanea_models_proj_2070_MPI85_tr_EM)
plot(castanea_models_proj_2070_MPI85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_MPI85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_MPI85_EM/proj_2070_MPI85_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2070_MPI85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_MPI85)

myBiomodRangeSize_2070_MPI85$Compt.By.Models

plot(myBiomodRangeSize_2070_MPI85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2070_MPI85$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_MPI85.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_MPI85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)

# wc21_GFDL_50_26 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_GFDL_50_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/GFDL/26/GFDL_26_41_60.tif")
names(wc21_GFDL_50_261)
names(wc21_GFDL_50_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                            'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_GFDL_50_261 <-
  subset(wc21_GFDL_50_261, c(4,10,12,15))
wc21_GFDL_50_26 <- crop(wc21_GFDL_50_261, extshp,mask=T)

names(wc21_GFDL_50_26)

# Project onto future conditions
castanea_models_proj_2050_GFDL26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                        proj.name = '2050_GFDL26',
                                                        new.env = wc21_GFDL_50_26,
                                                        models.chosen = 'all',
                                                        metric.binary = 'TSS',
                                                        build.clamping.mask = FALSE)


castanea_models_proj_2050_GFDL26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                    bm.proj = castanea_models_proj_2050_GFDL26_tr,
                                                                    proj.name = '2050_GFDL26_EM',
                                                                    models.chosen = 'all',
                                                                    metric.binary = 'all',
                                                                    metric.filter = 'all')

castanea_models_proj_2050_GFDL26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                      bm.proj = castanea_models_proj_2050_GFDL26_tr,
                                                                      proj.name = '2050_GFDL26_algo',
                                                                      models.chosen = 'all',
                                                                      metric.binary = 'all',
                                                                      metric.filter = 'all')


plot(castanea_models_proj_2050_GFDL26_tr_EM)
plot(castanea_models_proj_2050_GFDL26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_GFDL26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_GFDL26_EM/proj_2050_GFDL26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_GFDL26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_GFDL26)

myBiomodRangeSize_2050_GFDL26$Compt.By.Models

plot(myBiomodRangeSize_2050_GFDL26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2050_GFDL26$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_GFDL26.tif",overwrite=T)


# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_GFDL26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_GFDL_50_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_GFDL_50_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2041_2060/GFDL/85/wc21_85_gfdl_41_60.grd")
names(wc21_GFDL_50_851)
names(wc21_GFDL_50_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                            'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_GFDL_50_851 <-
  subset(wc21_GFDL_50_851, c(4,10,12,15))
wc21_GFDL_50_85 <- crop(wc21_GFDL_50_851, extshp,mask=T)


names(wc21_GFDL_50_85)

# Project onto future conditions
castanea_models_proj_2050_GFDL85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                        proj.name = '2050_GFDL85',
                                                        new.env = wc21_GFDL_50_85,
                                                        models.chosen = 'all',
                                                        metric.binary = 'TSS',
                                                        build.clamping.mask = FALSE)


castanea_models_proj_2050_GFDL85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                    bm.proj = castanea_models_proj_2050_GFDL85_tr,
                                                                    proj.name = '2050_GFDL85_EM',
                                                                    models.chosen = 'all',
                                                                    metric.binary = 'all',
                                                                    metric.filter = 'all')

castanea_models_proj_2050_GFDL85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                      bm.proj = castanea_models_proj_2050_GFDL85_tr,
                                                                      proj.name = '2050_GFDL85_algo',
                                                                      models.chosen = 'all',
                                                                      metric.binary = 'all',
                                                                      metric.filter = 'all')


plot(castanea_models_proj_2050_GFDL85_tr_EM)
plot(castanea_models_proj_2050_GFDL85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2050_GFDL85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2050_GFDL85_EM/proj_2050_GFDL85_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2050_GFDL85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2050_GFDL85)

myBiomodRangeSize_2050_GFDL85$Compt.By.Models

plot(myBiomodRangeSize_2050_GFDL85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")


writeRaster(myBiomodRangeSize_2050_GFDL85$Diff.By.Pixel,file = "./myBiomodRangeSize_2050_GFDL85.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2050_GFDL85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)



# wc21_GFDL_70_26 ---------------------------------------------------------

wc21_GFDL_70_261 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/GFDL/26/GFDL_26_61_80.tif")
names(wc21_GFDL_70_261)
names(wc21_GFDL_70_261) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                            'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_GFDL_70_261 <-
  subset(wc21_GFDL_70_261, c(4,10,12,15))
wc21_GFDL_70_26 <- crop(wc21_GFDL_70_261, extshp,mask=T)

names(wc21_GFDL_70_26)

# Project onto future conditions
castanea_models_proj_2070_GFDL26_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                        proj.name = '2070_GFDL26',
                                                        new.env = wc21_GFDL_70_26,
                                                        models.chosen = 'all',
                                                        metric.binary = 'TSS',
                                                        build.clamping.mask = FALSE)


castanea_models_proj_2070_GFDL26_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                    bm.proj = castanea_models_proj_2070_GFDL26_tr,
                                                                    proj.name = '2070_GFDL26_EM',
                                                                    models.chosen = 'all',
                                                                    metric.binary = 'all',
                                                                    metric.filter = 'all')

castanea_models_proj_2070_GFDL26_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                      bm.proj = castanea_models_proj_2070_GFDL26_tr,
                                                                      proj.name = '2070_GFDL26_algo',
                                                                      models.chosen = 'all',
                                                                      metric.binary = 'all',
                                                                      metric.filter = 'all')


plot(castanea_models_proj_2070_GFDL26_tr_EM)
plot(castanea_models_proj_2070_GFDL26_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_GFDL26 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_GFDL26_EM/proj_2070_GFDL26_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2070_GFDL26 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_GFDL26)

myBiomodRangeSize_2070_GFDL26$Compt.By.Models

plot(myBiomodRangeSize_2070_GFDL26$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "#99FFFF")



writeRaster(myBiomodRangeSize_2070_GFDL26$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_GFDL26.tif",overwrite=T)



# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_GFDL26$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)


# wc21_GFDL_70_85 ---------------------------------------------------------


# Load environmental variables extracted from worldclim

wc21_GFDL_70_851 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/Worldclim/2061_2080/GFDL/85/GFDL_85_61_80.tif")
names(wc21_GFDL_70_851)
names(wc21_GFDL_70_851) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10',
                            'bio11', 'bio12' , 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')
wc21_GFDL_70_851 <-
  subset(wc21_GFDL_70_851, c(4,10,12,15))
wc21_GFDL_70_85 <- crop(wc21_GFDL_70_851, extshp,mask=T)


names(wc21_GFDL_70_85)

# Project onto future conditions
castanea_models_proj_2070_GFDL85_tr <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                                        proj.name = '2070_GFDL85',
                                                        new.env = wc21_GFDL_70_85,
                                                        models.chosen = 'all',
                                                        metric.binary = 'TSS',
                                                        build.clamping.mask = FALSE)


castanea_models_proj_2070_GFDL85_tr_EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                                    bm.proj = castanea_models_proj_2070_GFDL85_tr,
                                                                    proj.name = '2070_GFDL85_EM',
                                                                    models.chosen = 'all',
                                                                    metric.binary = 'all',
                                                                    metric.filter = 'all')

castanea_models_proj_2070_GFDL85_tr_algo <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_algo, 
                                                                      bm.proj = castanea_models_proj_2070_GFDL85_tr,
                                                                      proj.name = '2070_GFDL85_algo',
                                                                      models.chosen = 'all',
                                                                      metric.binary = 'all',
                                                                      metric.filter = 'all')


plot(castanea_models_proj_2070_GFDL85_tr_EM)
plot(castanea_models_proj_2070_GFDL85_tr_algo)



# Load current and future binary projections
CurrentProj <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_CurrentEM/proj_CurrentEM_castanea_ensemble_TSSbin.tif")
FutureProj_2070_GFDL85 <- rast("E:/Doktora/0_Thesis/2_Thesis_Data/Raster/0_STACK/worldclim_current_future/castanea/proj_2070_GFDL85_EM/proj_2070_GFDL85_EM_castanea_ensemble_TSSbin.tif")



# Compute differences
myBiomodRangeSize_2070_GFDL85 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2070_GFDL85)

myBiomodRangeSize_2070_GFDL85$Compt.By.Models

plot(myBiomodRangeSize_2070_GFDL85$Diff.By.Pixel,
     col = c("#FF6666", "goldenrod1", "#E0E0E0","#006600"),
     colNA = "white")


writeRaster(myBiomodRangeSize_2070_GFDL85$Diff.By.Pixel,file = "./myBiomodRangeSize_2070_GFDL85.tif",overwrite=T)


# Calculate the SRC in hectar

myDiff <- myBiomodRangeSize_2070_GFDL85$Diff.By.Pixel
map.cellsize <- cellSize(myDiff, unit = "ha")
dfres <- data.frame()
for (thisvalue in c(-1,1,0,-2)) {
  tmp <- classify(myDiff, matrix(c(thisvalue, 1), ncol = 2), others = 0)*map.cellsize
  for(thismodel in 1:nlyr(tmp)){
    dfres <- rbind(dfres,
                   data.frame("model" = names(myDiff)[thismodel],
                              value = thisvalue,
                              area = sum(values(tmp[[thismodel]]), na.rm = TRUE))
    ) 
  }
}
show(dfres)
