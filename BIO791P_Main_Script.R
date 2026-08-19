#### BIO791P Research Project Code ####
## Louis James Boettger - 2025/26


# Install Packages (if not installed already) -----------------------------

# remove comment to install package
# install.packages("readxl")

# Load Packages -----------------------------------------------------------

invisible(lapply(c(
  "readxl", "elevatr", "terra", "sf", "dplyr", "wdpar", "geodata",
  "ggplot2", "maptiles", "tidyterra", "patchwork", 
  "rnaturalearth", "rnaturalearthdata", "vegan", "iNEXT", 
  "tidyverse", "scales", "tidyr", "FSA", "viridis",
  "corrplot", "rstatix", "car"), library, character.only = TRUE))


# Making Functions --------------------------------------------------------

#### crtl+alt+T to run entire section ####

# function for shannon diversity index
shannon <- function(x){
  total_abundance <- colSums(x)
  output <- diversity(total_abundance, index = 'shannon')
  print(output)
}
# function for simpson diversity index
simpson <- function(x){
  total_abundance <- colSums(x)
  output <- diversity(total_abundance, index = 'simpson')
  print(output)
}

# function to output both shannon and simpson indices at once
diversityindices <- function(x){
  shannon(x)
  simpson(x)
}
# function formatting incidence data for species accumulation 
freq_format <- function(x) {
  x <- ifelse(x != 0, 1, 0)
  num_days <- nrow(x)
  incidence_freq <- colSums(x)
  output <- c(num_days, incidence_freq)
}

# species accumulation plot function
speciesaccum_plot <- function(data,
                              line_colour = "",
                              fill_colour = "",
                              endpoint = NULL) {
  result <- iNEXT(
    data,
    q = 0,
    datatype = "incidence_freq",
    endpoint = endpoint
  )
  ggiNEXT(
    result,
    type = 1,
    color.var = "Order.q"
  ) +
    labs(
      x = "Days",
      y = "Species Richness"
    ) +
    coord_cartesian(xlim = c(0, 25)) +
    scale_color_manual(values = line_colour) +
    scale_fill_manual(values = fill_colour)
}



# Load Data ---------------------------------------------------------------

all_data <- read_xlsx("eDNA_Sampling.xlsx")
all_data$Location <- as.factor(all_data$Location)

indices_all <- read_xlsx("Indices_Scotland.xlsx")
indices_all$Location <- as.factor(indices_all$Location)

hvy_metal<- read_xlsx("Heavy_Metal.xlsx")
hvy_metal$Location <- as.factor(hvy_metal$Location)



# Site Maps ---------------------------------------------------------------

# download UK reserves 
uk_wdpa <- wdpa_fetch("United Kingdom")

# retrieve NNR
craigellachie <- uk_wdpa |>
  filter(grepl("Craigellachie", NAME)) |> 
  filter(DESIG=='National Nature Reserve')

craigellachie

# convert to shapefile
craigellachie_sf <- st_as_sf(craigellachie)

# reproject to national grid
craigellachie_sf <- st_transform(craigellachie_sf, 27700)

# plot outline of reserve
craigellachie_map <- ggplot() +
  geom_sf(data = craigellachie_sf, colour = "black")

craigellachie_map

#### Satellite Map ####
craigellachie_wgs <- st_transform(craigellachie_sf, 4326)

# add a buffer
craigellachie_buffer <- st_buffer(craigellachie_wgs, dist = 300)

# get satellite tiles
craigellachie_tiles <- get_tiles(craigellachie_buffer,
  provider = "Esri.WorldImagery", crop = TRUE, zoom = 16)

craigellachie_satellite_map <- ggplot() +
  geom_spatraster_rgb(data = craigellachie_tiles) +
  geom_sf(data = craigellachie_wgs, fill = NA,
          colour = "paleturquoise3", linewidth = 1) +
  labs(x = "Longitude", y = "Latitude", title = "C") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

craigellachie_satellite_map


sampling_points <- all_data |>
  filter(Sample_Type == "Soil") |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

craigellachie_sat_points <- ggplot()+
  geom_spatraster_rgb(data = craigellachie_tiles) +
  geom_sf(data = craigellachie_wgs, fill = NA, colour = "paleturquoise3", 
          linewidth = 1)+
  geom_sf(data = sampling_points, 
          colour = "maroon2", size = 2, shape = 18)+
  labs(title = "C", x= "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

craigellachie_sat_points


camera_audio_labels <- sampling_points |> 
  filter(Location %in% c(3, 18))

camera_labels <- sampling_points |> 
  filter(Location %in% c(13, 20))

audiomoth_labels <- sampling_points |> 
  filter(Location == "19")

metal_labels <- sampling_points |> 
  filter(Location %in% c(1, 5, 24, 4))

metal_camera <- sampling_points |> 
  filter(Location == "8")

metal_camera_audio <- sampling_points |> 
  filter(Location == "9")


sat_equip_map <- ggplot()+
  geom_spatraster_rgb(data = craigellachie_tiles) +
  
  geom_sf(data = craigellachie_wgs, fill = NA, colour = "white", 
          linewidth = 1)+
  
  geom_sf(data = sampling_points, colour = "steelblue3", size = 1)+
  
  geom_sf(data = camera_labels, colour = "cyan3", size = 1.5) +
  geom_sf(data = camera_audio_labels, colour = "red3", size = 1.5) +
  geom_sf(data = audiomoth_labels, colour = "palegreen3", size = 1.5) +
  geom_sf(data = metal_labels, colour = "magenta3", size = 1.5) +
  geom_sf(data = metal_camera, colour = "orange", size = 1.5) +
  geom_sf(data = metal_camera_audio, colour = "pink2", size = 1.5) +
  
  
  
  geom_sf_text(data = audiomoth_labels, aes(label = Location), colour = "white", 
               size = 1.5, fontface = "bold") +
  geom_sf_text(data = camera_labels, aes(label = Location), colour = "white",
               size = 1.5, fontface = "bold") +
  geom_sf_text(data = camera_audio_labels, aes(label = Location), colour = "white",
               size = 1.5, fontface = "bold") +
  geom_sf_text(data = metal_camera, aes(label = Location), colour = "white",
               size = 1.5, fontface = "bold") +
  geom_sf_text(data = metal_camera_audio, aes(label = Location), colour = "white",
               size = 1.5, fontface = "bold") +
  geom_sf_text(data = metal_labels, aes(label = Location), colour = "white",
               size = 1.5, fontface = "bold") +
  
  labs(title = "A", x= "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sat_equip_map



#### Elevation Map ####
craigellachie_elev <- get_elev_raster(
  locations = craigellachie_sf,
  z = 13,              
  clip = "locations")

craigellachie_elev <- rast(craigellachie_elev)

craigellachie_elev_map <- ggplot() +
  geom_spatraster(data = craigellachie_elev, na.rm = F) +
  scale_fill_viridis_c(name = "Elevation (m)", na.value= "NA") +
  geom_sf(data = craigellachie_sf, fill = NA, colour = "black") +
  labs(x = "Longitude", y = "Latitude", 
       title = "B")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

craigellachie_elev_map


# both maps side by side
craigellachie_satellite_map | craigellachie_elev_map


satellite_elevation_maps <- sat_equip_map | craigellachie_elev_map
ggsave(
  filename = "satellite_elevation_map.png",
  plot = satellite_elevation_maps,
  width = 8,
  height = 6,
  dpi = 300)


# reserve location within Scotland/UK
# UK administrative boundaries
uk <- gadm(country = "GBR", level = 1, path = tempdir())

# Convert to sf
uk_sf <- st_as_sf(uk)

# Extract Scotland
scotland <- uk_sf |>
  filter(NAME_1 == "Scotland")

# GPS point
craigellachie_gps_sf <- st_as_sf(
  data.frame(lon = -3.840780, lat = 57.191720),
  coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = uk_sf, fill = "grey92", color = "black", linewidth = 0.2) +
  geom_sf(data = craigellachie_gps_sf, color = "red", size = 2) +
  theme_minimal()

ggplot() +
  geom_sf(data = scotland, fill = "grey95", color = "black") +
  geom_sf(data = craigellachie_gps_sf, color = "red", size = 3) +
  geom_sf_text(data = craigellachie_gps_sf, aes(label = "Craigellachie"), nudge_y = 0.03) +
  coord_sf(xlim = c(-4.8, -2.5), ylim = c(56.7, 57.7), expand = FALSE) +
  theme_minimal()

cairngorms_np <- uk_wdpa |> filter(NAME == "Cairngorms", DESIG == "National Park")

# convert to shapefile
cairngorms_sf <- st_as_sf(cairngorms_np)

# reproject to national grid
cairngorms_sf <- st_transform(cairngorms_sf, 27700)

# plot outline of reserve
cairngorms_map <- ggplot() +
  geom_sf(data = cairngorms_sf, colour = "black")
cairngorms_map


cairngorms_map_wide <- ggplot()+
  geom_sf(data = uk_sf, fill = "grey95", color = "black") +
  geom_sf(data = cairngorms_sf, fill = "grey35", colour = "black")+
  geom_sf(data = craigellachie_sf, fill = "cyan2", colour = "cyan2")+
  coord_sf(xlim = c(-8, 0), ylim = c(55, 59), expand = F)+
  labs(x= "Longitude", y= "Latitude")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cairngorms_map_wide

cairngorms_map_close <- ggplot()+
  geom_sf(data = uk_sf, fill = "grey95", color = "black") +
  geom_sf(data = cairngorms_sf, fill = "mediumorchid2", colour = "black")+
  geom_sf(data = craigellachie_sf, fill = "paleturquoise3", colour = "black")+
  coord_sf(xlim = c(-4.5, -2.78), ylim = c(56.7, 57.45), expand = F)+
  labs(x= "Longitude", y= "Latitude", title = "B")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cairngorms_map_close


combined_map<-(cairngorms_map_wide|cairngorms_map_close)/(craigellachie_sat_points|craigellachie_elev_map)
combined_map

cairngorms_map_wide |(craigellachie_elev_map/sat_equip_map)

# Bird Species Accumulations ----------------------------------------------

birds_all <- read_excel('All_Combined80.xlsx')

common_species <- birds_all |>
  count(`Common name`, name = "Detections") |>
  filter(Detections >= 5) |>
  pull(`Common name`)

birds_filtered <- birds_all |>
  filter(`Common name` %in% common_species)


# separate raw data by location
locations <- split(birds_filtered, birds_filtered$Location)

# reformat data to have 1 row per sampling day
wide_locations <- lapply(locations, function(df) {
  df |>
    mutate(SamplingDay = paste(Year, Month, Day, sep = "-"),value = 1) |>
    select(SamplingDay, `Common name`, value) |>
    distinct() |>
    pivot_wider(names_from = `Common name`, values_from = value, values_fill = 0) |>
    select(-SamplingDay)
})

# use format function to format data for iNEXT, e.g. presence/absence only
birds_3_format <- freq_format(wide_locations[["3"]])
birds_9_format <- freq_format(wide_locations[["9"]])
birds_18_format <- freq_format(wide_locations[["18"]])
birds_19_format <- freq_format(wide_locations[["19"]])

# run iNEXT species accumulation graphs per location
birds_3_speciesaccum <- speciesaccum_plot(
  birds_3_format, line_colour = "cyan2", fill_colour = "cyan3")+
  theme_bw()+ 
  theme(legend.position = "none")+ 
  labs(x="",y="Species Richness", subtitle="A")+ 
  coord_cartesian(ylim=c(0,30))
birds_3_speciesaccum

birds_9_speciesaccum <- speciesaccum_plot(
  birds_9_format, line_colour = "magenta2", fill_colour = "magenta3")+
  theme_bw()+ 
  theme(legend.position = "none")+ 
  labs(x="",y="", subtitle = "B")+ 
  coord_cartesian(ylim=c(0,30))
birds_9_speciesaccum

birds_18_speciesaccum <- speciesaccum_plot(
  birds_18_format, 
  line_colour = "aquamarine3", fill_colour = "aquamarine4")+
  theme_bw()+ 
  theme(legend.position = "none", axis.title.x = element_text(hjust = -0.25))+ 
  labs(x = "Days", y="", subtitle = "C")+ 
  coord_cartesian(ylim=c(0,30))

birds_18_speciesaccum

birds_19_speciesaccum <- speciesaccum_plot(birds_19_format, 
                                           line_colour = "gold", fill_colour = "gold3")+ 
  theme_bw()+ theme(legend.position = "none") + labs(x="",y="", subtitle = "D")+ 
  coord_cartesian(ylim=c(0,30)) 
birds_19_speciesaccum

bird_accumulation_locations <- birds_3_speciesaccum | birds_9_speciesaccum |birds_18_speciesaccum | birds_19_speciesaccum
bird_accumulation_locations 

ggsave(
  filename = "bird_accumulation_locations.png",
  plot = bird_accumulation_locations,
  width = 12,
  height = 6,
  dpi = 300)



# combine location data into 1 list for an iNEXT with all accumulations
bird_sites <- list(
  "Location 3" = birds_3_format,
  "Location 9" = birds_9_format,
  "Location 18" = birds_18_format,
  "Location 19" = birds_19_format)

# result for each curve
birds_result <- iNEXT(bird_sites, q = 0, datatype = "incidence_freq", endpoint = 30)

# accumulation with 4 curves, 1 for each location
ggiNEXT(birds_result, type = 1, color.var = "Assemblage") +
  labs( x = "Sampling days", y = "Species richness")


# combine all accumulation data for a single accumulation curve
birds_combined <- birds_all |>
  mutate(SamplingDay = paste(Year, Month, Day, sep = "-"), value = 1) |>
  select(SamplingDay, `Common name`, value) |>
  distinct() |>
  pivot_wider(names_from = `Common name`, values_from = value, values_fill = 0) |>
  select(-SamplingDay)

# format combined data to 1 row per day & presence/absence only
birds_combined_format <- freq_format(birds_combined)

# plot combined data
speciesaccum_plot(birds_combined_format, line_colour = "black",
  fill_colour = "grey70")


# Bird Community PCoAs & Stats --------------------------------------------

#### Presence/absence community matrix ####

# cleaning data, and filtering per location per day
birds_pa <- birds_filtered |>
  mutate(SamplingDay = paste(Year, Month, Day, sep = "-")) |>
  mutate(value = 1) |>
  distinct(Location, SamplingDay, `Common name`, .keep_all = TRUE) |>
  select(Location, SamplingDay, `Common name`, value) |>
  pivot_wider(names_from = `Common name`, values_from = value, values_fill = 0)


# detection frequency community matrix
birds_abundance <- birds_filtered |>
  mutate(SamplingDay = paste(Year, Month, Day, sep = "-")) |>
  count(Location, SamplingDay, `Common name`) |>
  pivot_wider(names_from = `Common name`, values_from = n, values_fill = 0)

birds_pa_metadata <- birds_pa |>
  select(Location, SamplingDay)

birds_pa_matrix <- birds_pa |>
  select(-Location, -SamplingDay)


# Bray-Curtis dissimilarity
birds_pa_dist <- vegdist(birds_pa_matrix, method = "bray")

birds_pa_pcoa <- cmdscale(
  birds_pa_dist,
  eig = TRUE,
  k = 2)

birds_pa_variance <-
  round(
    birds_pa_pcoa$eig / sum(birds_pa_pcoa$eig[birds_pa_pcoa$eig > 0]) * 100, 1)

birds_pa_scores <- data.frame(
  birds_pa_metadata,
  PCoA1 = birds_pa_pcoa$points[,1],
  PCoA2 = birds_pa_pcoa$points[,2])

# ensure locations are factors not numeric
birds_pa_scores$Location <- factor(birds_pa_scores$Location)

# assign location colours
location_cols <- c(
  "18" = "aquamarine4",
  "19" = "gold",
  "3"  = "cyan2",
  "9"  = "magenta2")

#### P/A PCoA: Convex Hulls ####

# create hulls
birds_pa_hulls <- birds_pa_scores |>
  group_by(Location) |>
  slice(chull(PCoA1, PCoA2))

# PCoA plot
birds_pa_plot <-
  ggplot(birds_pa_scores,
         aes(PCoA1, PCoA2, colour = Location, fill = Location)) +
  geom_polygon(data = birds_pa_hulls, alpha = 0.2, linewidth = 0) +
  geom_point(size = 3) +
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = paste0("PCoA1 (", birds_pa_variance[1], "%)"),
    y = paste0("PCoA2 (", birds_pa_variance[2], "%)"))

birds_pa_plot

#### P/A PCoA: Ordination Ellipse ####
birds_pa_ellipse <- ggplot(birds_pa_scores,
       aes(x = PCoA1, y = PCoA2, colour = Location, fill = Location)) +
  
  stat_ellipse(geom = "polygon",
               alpha = 0.15,
               level = 0.95,      
               linewidth = 0.75) +
  
  geom_point(size = 3) +
  
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = paste0("PCoA1 (", birds_pa_variance[1], "%)"),
       y = paste0("PCoA2 (", birds_pa_variance[2], "%)"),
       subtitle = "A")+
  theme(legend.position = "none")

birds_pa_ellipse

#### Abundance PCoA: Convex Hulls####

birds_abundance_metadata <- birds_abundance |>
  select(Location, SamplingDay)

birds_abundance_matrix <- birds_abundance |>
  select(-Location, -SamplingDay)

birds_abundance_dist <-
  vegdist(birds_abundance_matrix, method = "bray")

birds_abundance_pcoa <-cmdscale(birds_abundance_dist, eig = TRUE, k = 2)

birds_abundance_variance <-
  round(birds_abundance_pcoa$eig /
      sum(birds_abundance_pcoa$eig[birds_abundance_pcoa$eig > 0]) * 100, 1)

birds_abundance_scores <-
  data.frame(birds_abundance_metadata,
    PCoA1 = birds_abundance_pcoa$points[,1],
    PCoA2 = birds_abundance_pcoa$points[,2])

birds_abundance_scores$Location <- factor(birds_abundance_scores$Location)


# creating convex hulls for abundance
birds_abundance_hulls <- birds_abundance_scores |>
  group_by(Location) |>
  slice(chull(PCoA1, PCoA2))

# plotting PCoA with convex hulls
birds_abundance_plot <-
  ggplot(birds_abundance_scores,
         aes(PCoA1, PCoA2, colour = Location, fill = Location)) +
  
  geom_polygon(data = birds_abundance_hulls, alpha = 0.2, linewidth = 0) +
  
  geom_point(size = 3) +
  
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  
  theme_bw() +
  
  labs(x = paste0("PCoA1 (", birds_abundance_variance[1], "%)"),
       y = paste0("PCoA2 (", birds_abundance_variance[2], "%)"))

birds_abundance_plot

#### Abundance PCoA: Ordination Ellipse ####
birds_abundance_ellipse <- ggplot(birds_abundance_scores,
       aes(x = PCoA1, y = PCoA2, colour = Location, fill = Location)) +
  
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,     
    linewidth = 0.75) +
  
  geom_point(size = 3) +
  
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = paste0("PCoA1 (", birds_abundance_variance[1], "%)"),
    y = paste0("PCoA2 (", birds_abundance_variance[2], "%)"),
    subtitle = "B")

birds_abundance_ellipse


#### PCoA Plots Combined ####

# convex hull plots side by side
birds_pa_plot | birds_abundance_plot

#ellipse plots side by side
bird_community_pcoas_ellipse <- birds_pa_ellipse | birds_abundance_ellipse

ggsave(filename = "bird_community_pcoas_ellipse.png",
  plot = bird_community_pcoas_ellipse,
  width = 12,
  height = 6,
  dpi = 300)

#### Community Stats ####

## presence/absence ##

# p/a Bray-Curtis permanova
birds_pa_permanova <-adonis2(birds_pa_matrix ~ Location, 
                             data = birds_pa_metadata, 
                             method = "bray", permutations = 999)

birds_pa_permanova

# p/a dispersion test
birds_pa_dispersion <- betadisper(
    birds_pa_dist,
    birds_pa_metadata$Location)

anova(birds_pa_dispersion)
TukeyHSD(birds_pa_dispersion)

## abundance ##

# abundance Bray-Curtis permanova
birds_abundance_permanova <- adonis2(birds_abundance_matrix ~ Location,
                                     data = birds_abundance_metadata,
                                     method = "bray", permutations = 999)

birds_abundance_permanova

# abundance dispersion test
birds_abundance_dispersion <- betadisper(birds_abundance_dist,
    birds_abundance_metadata$Location)

anova(birds_abundance_dispersion)


# Bird Diversity Indices --------------------------------------------------

bird_diversity_daily <- birds_filtered |>
  mutate(SamplingDay = paste(Year, Month, Day, sep = "-")) |>
  count(Location, SamplingDay, `Common name`) |>
  group_by(Location, SamplingDay) |>
  summarise(
    Shannon = diversity(n, index = "shannon"),
    Simpson = diversity(n, index = "simpson"),
    .groups = "drop")

bird_diversity_daily$Location <- as.factor(bird_diversity_daily$Location)

bird_shannon_boxp <- ggplot(bird_diversity_daily,
       aes(x = Location, y = Shannon, fill = Location)) +
  geom_boxplot(alpha = 0.65) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = "", y = "Shannon Diversity Index", subtitle = "A")+
  theme(legend.position = "none")

bird_shannon_boxp

bird_simpson_boxp <- ggplot(bird_diversity_daily,
       aes(x = Location, y =  Simpson, fill = Location)) +
  geom_boxplot(alpha=0.65) +
  geom_jitter(width = 0.15,alpha = 0.5, size = 2) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = "Location", y = "Simpson Diversity Index", subtitle = "B")+
  theme(axis.title.x = element_text(hjust = -0.135))

bird_simpson_boxp

bird_diversity_indices_boxp <- bird_shannon_boxp | bird_simpson_boxp
bird_diversity_indices_boxp

ggsave(filename = "bird_diversity_indices_boxp.png",
  plot = bird_diversity_indices_boxp,
  width = 12,
  height = 6,
  dpi = 300)


kruskal.test(Shannon ~ Location, data = bird_diversity_daily)

dunnTest(Shannon ~ factor(Location), 
         data = bird_diversity_daily, method = "holm")

kruskal.test(Simpson ~ Location, data = bird_diversity_daily)

dunnTest(Simpson ~ factor(Location),
         data = bird_diversity_daily, method = "holm")

# Acoustic Indices: PCoAs & Stats -----------------------------------------

# splitting indices to locations per day, using mean of each index
acoustic_daily <- indices_all |> 
  mutate(SamplingDay = substr(Recording, 1, 8)) |> 
  group_by(Location, SamplingDay) |> 
  summarise(across(c(ACI.tot, ACI.min, ADI, AEve, Bio, Bio.sd,
        H, Hs, Ht, NDSI, NDSI.anthro, NDSI.bio), mean, na.rm = TRUE), 
        n_recordings = n(), .groups = "drop")

#### Making PCoA ####
acoustic_metadata <- acoustic_daily |>
  select(Location, SamplingDay)

# matrix of acoustic indices per location per day
acoustic_matrix <- acoustic_daily |>
  select(ACI.tot, ADI, AEve, Bio, NDSI.anthro, NDSI.bio)

# re scaling, different indices have different scales e.g. 1-2, >2000, 0-1
acoustic_matrix_scaled <- scale(acoustic_matrix)
acoustic_dist <- dist(acoustic_matrix_scaled, method = "euclidean")

acoustic_pcoa <- cmdscale(acoustic_dist, eig = TRUE, k = 2)

acoustic_variance <-
  round(acoustic_pcoa$eig /
      sum(acoustic_pcoa$eig[acoustic_pcoa$eig > 0]) * 100, 1)

acoustic_scores <- data.frame(acoustic_metadata,
  PCoA1 = acoustic_pcoa$points[,1],
  PCoA2 = acoustic_pcoa$points[,2])

acoustic_scores$Location <- factor(acoustic_scores$Location)

#### Acoustic Indices PCoA: Convex Hulls ####

# creating convex hulls
acoustic_hulls <- acoustic_scores |>
  group_by(Location) |>
  slice(chull(PCoA1, PCoA2))

# PCoA plot
ggplot(acoustic_scores,
       aes(PCoA1, PCoA2, colour = Location, fill = Location)) +
  
  geom_polygon(data = acoustic_hulls, alpha = 0.2, linewidth = 0) +
  
  geom_point(size = 3) +
  
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = paste0("PCoA1 (", acoustic_variance[1], "%)"),
    y = paste0("PCoA2 (", acoustic_variance[2], "%)"))


#### Acoustic PCoA: Ordination Ellipse ####


# PCoA Plot
PCoA_Indices_Ellipses <- ggplot(acoustic_scores,
       aes(x = PCoA1, y = PCoA2, colour = Location, fill = Location)) +
  
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,      
    linewidth = 0.75) +
  
  geom_point(size = 3) +
  
  scale_colour_manual(values = location_cols) +
  scale_fill_manual(values = location_cols)+
  theme_bw() +
  labs(x = paste0("PCoA1 (", acoustic_variance[1], "%)"),
       y = paste0("PCoA2 (", acoustic_variance[2], "%)"))

PCoA_Indices_Ellipses

ggsave(filename = "PCoA_Indices_Ellipses.png",
  plot = PCoA_Indices_Ellipses,
  width = 8,
  height = 6,
  dpi = 300)


#### Acoustic Indices Stats ####

# multivariate euclidean permanova
acoustic_permanova <- adonis2(acoustic_matrix_scaled ~ Location, data = acoustic_metadata,
                              method = "euclidean", permutations = 999)

acoustic_permanova

# indices dispersion test
acoustic_dispersion <- betadisper(acoustic_dist, acoustic_metadata$Location)

anova(acoustic_dispersion)

# post-hoc test
TukeyHSD(acoustic_dispersion)

plot(acoustic_dispersion)


# correlations between indices (checking redundancy)
acoustic_correlations <-
  cor(acoustic_matrix_scaled, acoustic_pcoa$points[,1:2])

round(acoustic_correlations, 2)

corrplot(cor(acoustic_matrix_scaled), method="circle")



# reformat to long form data for desired indices
acoustic_daily_long <- acoustic_daily |> 
  select(Location, SamplingDay, ACI.tot, ADI, AEve, 
         Bio, NDSI.anthro, NDSI.bio) |> 
  pivot_longer(cols = -c(Location, SamplingDay),
    names_to = "Index",
    values_to = "Value")

# make location a factor instead of numeric
acoustic_daily_long$Location <- as.factor(acoustic_daily_long$Location)


# kruskal-wallis tests for all indices
acoustic_kruskal <- acoustic_daily_long |>
  group_by(Index) |>
  summarise(statistic = kruskal.test(Value ~ Location)$statistic,
    p = kruskal.test(Value ~ Location)$p.value, 
    df= kruskal.test(Value ~Location)$parameter,
    .groups = "drop")  |> 
  mutate(p_adj = p.adjust(p, method = "holm"))

acoustic_kruskal


# mean values for each index
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_NDSI_anthro = mean(NDSI.anthro, na.rm = TRUE))
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_Bio= mean(Bio, na.rm = TRUE))
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_ADI= mean(ADI, na.rm = TRUE))
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_ACI.tot= mean(ACI.tot, na.rm = TRUE))
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_NDSI.bio= mean(NDSI.bio, na.rm = TRUE))
acoustic_daily |>
  group_by(Location) |>
  summarise(mean_AEve= mean(AEve, na.rm = TRUE))

# Acoustic Indices: Other Plots -------------------------------------------

# all desired indices, boxplots per location per day
acoustic_indices_boxp <- acoustic_daily_long |> 
  ggplot(aes(x = factor(Location), y = Value, fill = Location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  facet_wrap(~Index, scales = "free_y", nrow =2) +
  scale_fill_manual(values = location_cols)+
  labs(x = "Location")+
  theme(legend.position = "bottom")+
  theme_bw()

acoustic_indices_boxp

ggsave(
  filename = "acoustic_indices_boxp.png",
  plot = acoustic_indices_boxp,
  width = 8,
  height = 6,
  dpi = 300)



# Heavy Metal Plots -------------------------------------------------------

hvy_metal<- read_xlsx("Heavy_Metal.xlsx")
hvy_metal$Location <- as.factor(hvy_metal$Location)

heavy_metal_avgconc_bar <- hvy_metal |> 
  filter(Units == "avg_conc_mg/g" & Measurement >0) |> 
  ggplot(aes(Location, Measurement)) +
  geom_col(fill = "aquamarine3") +
  facet_wrap(~Heavy_Metal, scales = "free_y") +
  scale_y_continuous(labels = label_number())+
  ylab("Mean Concentration (mg/g)")+
  theme_bw()

heavy_metal_avgconc_bar

ggsave(filename = "heavy_metal_avgconc_bar.png",
  plot = heavy_metal_avgconc_bar,
  width = 10,
  height = 6,
  dpi = 300)


# Heavy Metal Stats -------------------------------------------------------

hvy_metal_distA9 <- hvy_metal |> 
  left_join(
    all_data |> 
      select(Location, Dist_A9), by = "Location")



hvy_metal_distA9_kruskal <- hvy_metal_distA9 |>
  filter(Units == "mg/g") |>
  group_by(Heavy_Metal) |>
  kruskal_test(Measurement ~ Dist_A9)

hvy_metal_distA9_kruskal



assumption_results <- hvy_metal |>
  filter(Units == "mg/g") |>
  group_by(Heavy_Metal) |>
  group_modify(~{
    model <- aov(Measurement ~ factor(Location), data = .x)
    
    data.frame(shapiro_p = shapiro.test(residuals(model))$p.value,
      levene_p = leveneTest(
        Measurement ~ factor(Location),
        data = .x)$`Pr(>F)`[1])
  })

assumption_results


hvy_metal_kruskal <- hvy_metal |>
  filter(Units == "mg/g") |>
  group_by(Heavy_Metal) |>
  kruskal_test(Measurement ~ Location)

hvy_metal_kruskal

hvy_metal_dunn <- hvy_metal |>
  filter(Units == "mg/g") |>
  group_by(Heavy_Metal) |>
  dunn_test(
    Measurement ~ Location,
    p.adjust.method = "holm")

hvy_metal_dunn

hvy_metal_dunn_sig <- hvy_metal_dunn |> 
  filter(p.adj < 0.05)
hvy_metal_dunn_sig

print(hvy_metal_dunn_sig, n = Inf)


# Heavy Metal PCoA --------------------------------------------------------


#### creating the pcoa ####
metal_matrix <- hvy_metal |>
  filter(Units == "avg_conc_mg/g") |>
  select(Location, Heavy_Metal, Measurement) |>
  pivot_wider(
    names_from = Heavy_Metal,
    values_from = Measurement)

metal_metadata <- metal_matrix |>
  select(Location)

metal_values <- metal_matrix |>
  select(-Location)

# standardising the  metals
metal_values_scaled <- scale(metal_values)


metal_dist <- dist(
  metal_values_scaled,
  method = "euclidean")


metal_pcoa <- cmdscale(
  metal_dist,
  eig = TRUE,
  k = 2)

metal_variance <- round(
  metal_pcoa$eig /
    sum(metal_pcoa$eig[metal_pcoa$eig > 0]) * 100, 1)


metal_scores <- data.frame(metal_metadata,
  PCoA1 = metal_pcoa$points[, 1],
  PCoA2 = metal_pcoa$points[, 2])

metal_scores$Location <- factor(metal_scores$Location)

#### plotting the PCoA ####
PCoA_Metals <- ggplot(
  metal_scores,
  aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(colour = Location), size = 4) +
  geom_text(aes(label = Location),
    vjust = -1, size = 4) +
  theme_bw() +
  labs(x = paste0("PCoA1 (", metal_variance[1], "%)"),
    y = paste0("PCoA2 (", metal_variance[2], "%)"))+
  theme(legend.position = "none")

PCoA_Metals


ggsave(
  filename = "PCoA_Metals.png",
  plot = PCoA_Metals,
  width = 8,
  height = 6,
  dpi = 300)



# Combined Summaries ------------------------------------------------------

acoustic_summary <- acoustic_daily |>
  group_by(Location) |>
  summarise(
    Mean_ACI = mean(ACI.tot, na.rm = TRUE),
    Mean_ADI = mean(ADI, na.rm = TRUE),
    Mean_AEve = mean(AEve, na.rm = TRUE),
    Mean_Bio = mean(Bio, na.rm = TRUE),
    Mean_NDSI_anthro = mean(NDSI.anthro, na.rm = TRUE),
    Mean_NDSI_bio = mean(NDSI.bio, na.rm = TRUE),
    .groups = "drop")

bird_diversity_summary <- bird_diversity_daily |>
  group_by(Location) |>
  summarise(
    Bird_Shannon = mean(Shannon, na.rm = TRUE),
    Bird_Simpson = mean(Simpson, na.rm = TRUE),
    .groups = "drop")

heavy_metal_summary <- hvy_metal |>
  filter(Units == "mg/g") |>
  group_by(Location, Heavy_Metal) |>
  summarise(Mean_Concentration = mean(Measurement, na.rm = TRUE), 
            .groups = "drop") |>
  pivot_wider(
    names_from = Heavy_Metal,
    values_from = Mean_Concentration,
    names_prefix = "Mean_")

combined_summary <- acoustic_summary |>
  full_join(bird_diversity_summary, by = "Location") |>
  full_join(heavy_metal_summary, by = "Location") |>
  full_join(eDNA_location_summary, by = "Location") |>
  arrange(as.numeric(as.character(Location)))

combined_summary

# eDNA & Acoustics & Metals Integrated ------------------------------------


#### Acoustic Indices PCoA: Profile & eDNA ####
acoustic_location <- acoustic_daily |>
  group_by(Location) |>
  summarise(
    ACI.tot = mean(ACI.tot, na.rm = TRUE),
    ADI = mean(ADI, na.rm = TRUE),
    AEve = mean(AEve, na.rm = TRUE),
    Bio = mean(Bio, na.rm = TRUE),
    NDSI.anthro = mean(NDSI.anthro, na.rm = TRUE),
    NDSI.bio = mean(NDSI.bio, na.rm = TRUE),
    .groups = "drop")

acoustic_location_matrix <- acoustic_location |>
  select(ACI.tot, ADI, AEve, Bio, NDSI.anthro, NDSI.bio)

acoustic_location_scaled <- scale(acoustic_location_matrix)

acoustic_location_dist <- dist(acoustic_location_scaled, method = "euclidean")

acoustic_location_pcoa <- cmdscale(acoustic_location_dist, eig = TRUE, k = 2)

acoustic_location_variance <- round(
  acoustic_location_pcoa$eig /
    sum(acoustic_location_pcoa$eig[
      acoustic_location_pcoa$eig > 0]) * 100, 1)

acoustic_location_scores <- data.frame(Location = acoustic_location$Location,
  PCoA1 = acoustic_location_pcoa$points[, 1],
  PCoA2 = acoustic_location_pcoa$points[, 2])



acoustic_edna_pcoa <- acoustic_location_scores |> 
  left_join(
    combined_summary |> 
      select(Location, eDNA_Shannon, eDNA_Simpson, eDNA_Species_Richness),
    by = "Location")

acoustic_edna_pcoa_plot <- ggplot(
  acoustic_edna_pcoa,
  aes(x = PCoA1, y = PCoA2)) +
  
  geom_point(
    aes(colour = eDNA_Shannon, size = eDNA_Species_Richness),
    alpha = 0.9) +
  
  geom_text(
    aes(label = Location), vjust = -1.2, size = 4) +
  
  scale_colour_viridis_c(
    name = "eDNA Shannon\nDiversity") +
  scale_size_continuous(
    name = "eDNA Species\nRichness") +
  theme_bw() +
  labs(x = paste0("PCoA1 (",acoustic_location_variance[1],"%)"),
    y = paste0("PCoA2 (",acoustic_location_variance[2],"%)"),
    subtitle = "A")+
  theme(legend.position = "none")

acoustic_edna_pcoa_plot

ggsave(
  filename = "acoustic_edna_pcoa_plot.png",
  plot = acoustic_edna_pcoa_plot,
  width = 10,
  height = 8,
  dpi = 300)

#### Metal & eDNA ####
metal_edna <- metal_scores |>
  left_join(
    eDNA_location_summary |>
      select(Location, eDNA_Shannon, eDNA_Simpson, eDNA_Species_Richness),
    by = "Location")

PCoA_Metals_eDNA <- ggplot(
  metal_edna,
  aes(x = PCoA1, y = PCoA2)) +
  
  geom_point(
    aes(colour = eDNA_Shannon, size = eDNA_Species_Richness),
    alpha = 0.9) +
  
  geom_text(aes(label = Location),
    vjust = -1.2, size = 4) +
  scale_size_continuous(name = "eDNA Species\nRichness") +
  scale_colour_viridis_c(name = "eDNA Shannon\nDiversity",
    breaks = c(0.5, 1.0, 1.5, 2.0),
    limits = c(0.5, 2.0)) +
  theme_bw() +
  
  labs(x = paste0("PCoA1 (",metal_variance[1],"%)"),
    y = paste0("PCoA2 (",metal_variance[2],"%)"),
    subtitle = "B")+
  theme(legend.position = "right")

PCoA_Metals_eDNA


ggsave(
  filename = "PCoA_Metals_eDNA.png",
  plot = PCoA_Metals_eDNA,
  width = 10,
  height = 8,
  dpi = 300)


acoustic_metal_edna_pcoa <- acoustic_edna_pcoa_plot | PCoA_Metals_eDNA
acoustic_metal_edna_pcoa

ggsave(
  filename = "acoustic_metal_edna_pcoa.png",
  plot = acoustic_metal_edna_pcoa,
  width = 12,
  height = 8,
  dpi = 300)




acoustic_means_long <- combined_summary |>
  rename(
    `Mean ACI` = Mean_ACI,
    `Mean ADI` = Mean_ADI,
    `Mean AEve` = Mean_AEve,
    `Mean Bio` = Mean_Bio,
    `Mean NDSI.anthro` = Mean_NDSI_anthro,
    `Mean NDSI.bio` = Mean_NDSI_bio) |>
  pivot_longer(
    cols = c(
      `Mean ACI`, `Mean ADI`, `Mean AEve`,
      `Mean Bio`, `Mean NDSI.anthro`, `Mean NDSI.bio`),
    names_to = "Acoustic_Index",
    values_to = "Acoustic_Value")


eDNA_diversity_long <- combined_summary |>
  pivot_longer(
    cols = c(
      Bird_Shannon, Bird_Simpson,
      eDNA_Shannon, eDNA_Simpson,
      eDNA_Species_Richness),
    names_to = "Diversity_Index",
    values_to = "Diversity_Value")

combined_summary_long <- acoustic_means_long |>
  left_join(eDNA_diversity_long, by = "Location")


acoustic_eDNA_diversity_facet <- ggplot(
  acoustic_means_long,
  aes(x = Acoustic_Value,y = eDNA_Shannon)) +
  geom_point(aes(colour = Location),size = 3, na.rm = TRUE) +
  geom_smooth(method = "lm", se = T, na.rm = TRUE) +
  facet_wrap(~ Acoustic_Index, scales = "free") +
  theme_bw() +
  labs(x = "Acoustic Index Value", y = "eDNA Shannon Diversity")+
  scale_colour_manual(values = location_cols)

acoustic_eDNA_diversity_facet

ggsave(
  filename = "acoustic_eDNA_diversity_facet.png",
  plot = acoustic_eDNA_diversity_facet,
  width = 12,
  height = 8,
  dpi = 300)


metal_means_long <- combined_summary |>
  rename(
    `Mean Cd2288` = Mean_Cd2288,
    `Mean Cu3247` = Mean_Cu3247,
    `Mean Cr2835` = Mean_Cr2835,
    `Mean Pb2203` = Mean_Pb2203,
    `Mean Fe2599` = Mean_Fe2599,
    `Mean Zn2025` = Mean_Zn2025) |>
  pivot_longer(
    cols = c(
      `Mean Cd2288`, `Mean Cu3247`, `Mean Cr2835`,
      `Mean Pb2203`, `Mean Fe2599`, `Mean Zn2025`),
    names_to = "Heavy_Metal",
    values_to = "Metal_Concentration")


metal_eDNA_diversity_facet <- metal_means_long |> 
  filter(Location %in% c(1, 4, 5, 8, 9, 24)) |> 
           ggplot(aes(x = Metal_Concentration, y = eDNA_Shannon)) +
  geom_point(aes(colour = Location),size = 3, na.rm = TRUE) +
  geom_smooth(method = "lm", se = T, na.rm = TRUE) +
  facet_wrap(~ Heavy_Metal, scales = "free") +
  theme_bw() +
  labs(x = "Metal Concentration (mg/g)", y = "eDNA Shannon Diversity")+
  scale_colour_brewer(palette="Set2")

metal_eDNA_diversity_facet

ggsave(
  filename = "metal_eDNA_diversity_facet.png",
  plot = metal_eDNA_diversity_facet,
  width = 12,
  height = 8,
  dpi = 300)


metal_eDNA_spearman <- metal_means_long |>
  filter(Location %in% c(1, 4, 5, 8, 9, 24)) |> 
  group_by(Heavy_Metal) |>
  summarise(
    spearman = cor(Metal_Concentration, eDNA_Shannon, method = "spearman"),
    p_value = cor.test(Metal_Concentration, eDNA_Shannon, method = "spearman")$p.value)
metal_eDNA_spearman

acoustic_eDNA_spearman <- acoustic_means_long |>
  filter(Location %in% c(3, 9, 18, 19)) |> 
  group_by(Acoustic_Index) |>
  summarise(
    spearman = cor(Acoustic_Value, eDNA_Shannon, method = "spearman"),
    p_value = cor.test(Acoustic_Value, eDNA_Shannon, method = "spearman")$p.value)
acoustic_eDNA_spearman
