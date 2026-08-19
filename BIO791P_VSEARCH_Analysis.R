

# Load Packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(vegan)
library(tibble)



# Load & Prep Data --------------------------------------------------------

# excel sheet with location & barcode key
eDNA_Sampling <- read_xlsx(
  "eDNA_Sampling.xlsx") |> 
  mutate(
    Barcode = as.character(Barcode))


# vsearch output files
vsearch_files <- list.files(
  "Blast_Output/All_Taxa",
  pattern = "\\.txt$",
  full.names = TRUE)


# function for reading in and formatting vsearch output
read_vsearch <- function(file) {
  
  dat <- read_tsv(
    file,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  dat |> 
    mutate(Barcode = str_extract(
        basename(file),
        "(?<=barcode)[0-9]+"),
      Barcode = str_remove(
        Barcode,
        "^0+"),
      
      OTU = X1,
      OTU_ID = str_remove(
        X1,
      ";size=[0-9]+$"),
      OTU_size = as.numeric(
        str_extract(
          X1,
          "(?<=size=)[0-9]+")),
      
      Reference = X2,
      
      Percent_identity = as.numeric(
        X3),
      
      Taxon_raw = X13) |> 
    
    select(
      Barcode,
      OTU_ID,
      OTU_size,
      Reference,
      Percent_identity,
      Taxon_raw)
}


vsearch_all <- map_dfr(
  vsearch_files,
  read_vsearch)



# make sure barcodes have matched properly
barcode_check <- tibble(
  Barcode = sort(
    unique(vsearch_all$Barcode))) |> 
  left_join(eDNA_Sampling |> 
      select(Barcode, Location), by = "Barcode")

print(barcode_check, n = Inf)



# Taxonomy Functions ------------------------------------------------------


extract_taxa <- function(x) {
  
  if (
    length(x) == 0 ||
    is.na(x) ||
    x == ""
  ) {
    return(character(0))
  }
  
  x |> 
    str_split("\\|") |> 
    unlist() |> 
    str_trim() |> 
    unique()
}


# extract genus name from species column

extract_genus <- function(x) {
  
  taxa <- extract_taxa(x)
  
  # Remove unidentified annotations
  
  taxa <- taxa[!str_detect(taxa, regex("unidentified", ignore_case = TRUE))]

  if (length(taxa) == 0) {
    return(character(0))
  }
  
  genera <- str_extract(
    taxa,
    "^[A-Za-z]+")
  
  genera <- unique(
    genera[!is.na(genera) & genera != ""])
  
  genera
}



# extract species name & clean; remove uncertain species & other annotations

extract_species <- function(x) {
  
  taxa <- extract_taxa(x)
  
  if (length(taxa) == 0) {
    return(character(0))
  }
  
  taxa <- taxa[!str_detect(taxa,
      regex(
        "unidentified|_nr\\.?$|_sp\\.?$|_cf\\.?$",ignore_case = TRUE))]
  
  if (length(taxa) == 0) {
    return(character(0))
  }
  
  species <- taxa[
    str_detect(
      taxa,
      "^[A-Za-z]+_[A-Za-z]+$")]
  
  unique(species)
}



# taxonomy for each OTU; identify ambiguous, and unidentified etc

resolve_otu <- function(df) {
  
  best_id <- max(
    df$Percent_identity,
    na.rm = TRUE
  )
  
  best <- df |> 
    filter(
      Percent_identity == best_id)
  
  genera <- unique(
    unlist(
      lapply(
        best$Taxon_raw,
        extract_genus)))
  
  genera <- genera[!is.na(genera) & genera != ""]

  
  if (length(genera) == 0) {
    
    return(
      tibble(
        Genus = NA_character_,
        Species = NA_character_,
        Assignment_level = "Unidentified",
        Best_identity = best_id))
  }
  
  if (length(genera) > 1) {
    
    return(
      tibble(
        Genus = NA_character_,
        Species = NA_character_,
        Assignment_level = "Ambiguous",
        Best_identity = best_id))
  }
  
  genus <- genera[1]
  
  species <- unique(
    unlist(
      lapply(
        best$Taxon_raw,
        extract_species)))
  
  species <- species[!is.na(species) & species != ""]
  
  if (length(species) == 1) {
    
    return(
      tibble(
        Genus = genus,
        Species = species[1],
        Assignment_level = "Species",
        Best_identity = best_id))
  }
  
  if (length(species) > 1) {
    
    return(
      tibble(
        Genus = genus,
        Species = NA_character_,
        Assignment_level = "Genus",
        Best_identity = best_id))
  }
  
  tibble(
    Genus = genus,
    Species = NA_character_,
    Assignment_level = "Genus",
    Best_identity = best_id)
}



# Resolve Taxonomy --------------------------------------------------------

otu_taxonomy <- vsearch_all |>
  group_by(
    OTU_ID
  ) |>
  group_modify(
    ~ resolve_otu(.x)
  ) |>
  ungroup()



# add taxonomy to vsearch data
vsearch_all <- vsearch_all |> 
  left_join(
    otu_taxonomy,
    by = "OTU_ID")


# table for OTUs
otu_table <- vsearch_all |> 
  distinct(
    Barcode,
    OTU_ID,
    .keep_all = TRUE)



# filter out scat samples, only keeping soil samples
otu_soil <- otu_table  |> 
  left_join(
    eDNA_Sampling |>
      select(Barcode, Sample_Type), by = "Barcode") |>
  filter(Sample_Type == "Soil") |>
  select(-Sample_Type)


# summary for taxonomy
taxonomic_summary <- otu_soil |>
  distinct(OTU_ID, .keep_all = TRUE) |>
  count(Assignment_level, sort = TRUE)

taxonomic_summary


# taxpnomic abundance

taxonomic_abundance <- otu_soil |>
  distinct(OTU_ID, .keep_all = TRUE) |>
  group_by(Assignment_level) |>
  summarise(n_OTUs = n(),
    total_reads = sum(OTU_size, na.rm = TRUE),
    .groups = "drop") |>
  arrange(desc(total_reads))

taxonomic_abundance


sort(table(vsearch_all$Taxon_raw),decreasing = T)

# Abundance & P/A Matrices ------------------------------------------------


#### species abundance matrix ####

eDNA_species <- otu_soil |>
  filter(
    Assignment_level == "Species") |>
  group_by(Barcode, Species) |>
  summarise(Abundance = sum(OTU_size, na.rm = TRUE),
            .groups = "drop") |>
  pivot_wider(
    names_from = Species,
    values_from = Abundance,
    values_fill = 0)


#### species presence/absence matrix ####

eDNA_species_pa <- eDNA_species |>
  mutate(across(-Barcode, ~ ifelse(. > 0, 1, 0)))


#### genus abundance matrix ####

eDNA_genus <- otu_soil |>
  filter(
    Assignment_level %in% c("Species", "Genus")) |>
  group_by(Barcode, Genus) |>
  summarise(Abundance = sum(OTU_size, na.rm = TRUE),
    .groups = "drop") |>
  pivot_wider(
    names_from = Genus,
    values_from = Abundance,
    values_fill = 0)


#### genus p/a matrix ####

eDNA_genus_pa <- eDNA_genus |>
  mutate(across(-Barcode, ~ ifelse(
        . > 0, 1, 0)))



# Species & Genus Summaries -----------------------------------------------

#### species summary ####

species_summary <- otu_soil |>
  filter(Assignment_level == "Species") |>
  group_by(Species) |>
  summarise(n_OTUs = n_distinct(OTU_ID),
    total_reads = sum(OTU_size, na.rm = TRUE),
    .groups = "drop") |>
  arrange(desc(total_reads))

species_summary


#### genus summary ####

genus_summary <- otu_soil |>
  filter(Assignment_level %in% c(
      "Species",
      "Genus")) |>
  group_by(Genus) |>
  summarise(n_OTUs = n_distinct(OTU_ID),
    total_reads = sum(OTU_size, na.rm = TRUE),
    .groups = "drop") |>
  arrange(desc(total_reads))

genus_summary


#### amiguous OTUs (OTUs with multiple top IDs)

ambiguous_otus <- otu_soil |>
  filter(Assignment_level == "Ambiguous") |>
  distinct(OTU_ID, .keep_all = TRUE) |>
  arrange(desc(OTU_size))

print(ambiguous_otus, n = Inf)


# checking specific ambiguous OTU

vsearch_all |>
  filter(OTU_ID == "OTU_85555") |>
  select(OTU_ID, OTU_size, Reference,
    Percent_identity, Taxon_raw, Genus,
    Species, Assignment_level)



# Species Abundance PCoA --------------------------------------------------

#### making the pcoa ####
species_abundance_matrix <- eDNA_species |>
  column_to_rownames(
    "Barcode")

eDNA_metadata <- eDNA_Sampling |>
  filter(Sample_Type == "Soil") |>
  filter(Barcode %in% rownames(species_abundance_matrix)) |>
  mutate(
    Location = factor(Location),
    Habitat = factor(Habitat),
    Trees = factor(Trees),
    Sub_Habitat = factor(Sub_Habitat),
    Elevation_Group = factor(Elevation_Group)) |>
  distinct(Barcode, .keep_all = TRUE) |>
  column_to_rownames("Barcode")

species_abundance_dist <- vegdist(
  species_abundance_matrix, method = "bray")

species_abundance_pcoa <- cmdscale(
  species_abundance_dist, eig = TRUE, k = 2)

species_abundance_variance <- round(
  species_abundance_pcoa$eig /
    sum(
      species_abundance_pcoa$eig[
        species_abundance_pcoa$eig > 0]) * 100,1)

species_abundance_scores <- data.frame(
  Barcode = rownames(species_abundance_pcoa$points),
  PCoA1 = species_abundance_pcoa$points[, 1],
  PCoA2 = species_abundance_pcoa$points[, 2]) |>
  left_join(eDNA_Sampling,by = "Barcode") |>
  mutate(
    Location = factor(Location),
    Sample_Type = factor(Sample_Type),
    Habitat = factor(Habitat),
    Trees = factor(Trees),
    Sub_Habitat = factor(Sub_Habitat))

#### species abundance PCoA plot ####
ggplot(
  species_abundance_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Dist_A9_group,
    fill = Dist_A9_group)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linewidth = 0.75) +
  geom_point(
    size = 3) +
  theme_bw() +
  labs(
    colour = "Dist_A9_group",
    fill = "Dist_A9_group",
    x = paste0("PCoA1 (",species_abundance_variance[1],"%)"),
    y = paste0("PCoA2 (",species_abundance_variance[2],"%)"))



# Species P/A PCoA --------------------------------------------------------

#### making the PCoA ####
species_pa_matrix <- eDNA_species_pa |>
  column_to_rownames(
    "Barcode")

species_pa_dist <- vegdist(
  species_pa_matrix,
  method = "bray")

species_pa_pcoa <- cmdscale(
  species_pa_dist,
  eig = TRUE,
  k = 2)

species_pa_variance <- round(
  species_pa_pcoa$eig /
    sum(
      species_pa_pcoa$eig[species_pa_pcoa$eig > 0]) *100, 1)

species_pa_scores <- data.frame(
  Barcode = rownames(
    species_pa_pcoa$points),
  PCoA1 = species_pa_pcoa$points[, 1],
  PCoA2 = species_pa_pcoa$points[, 2]) |>
  left_join(
    eDNA_Sampling,
    by = "Barcode") |>
  mutate(
    Location = factor(Location),
    Sample_Type = factor(Sample_Type),
    Habitat = factor(Habitat),
    Trees = factor(Trees),
    Sub_Habitat = factor(Sub_Habitat))

#### species p/a PCoA plot ####
ggplot(
  species_pa_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Dist_A9_group,
    fill = Dist_A9_group)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linewidth = 0.75) +
  geom_point(
    size = 3) +
  theme_bw() +
  labs(
    colour = "Dist_A9_group",
    fill = "Dist_A9_group",
    x = paste0("PCoA1 (",species_pa_variance[1],"%)"),
    y = paste0("PCoA2 (",species_pa_variance[2],"%)"))



# Genus Abundance PCoA ----------------------------------------------------

#### making the PCoA ####
genus_abundance_matrix <- eDNA_genus |>
  column_to_rownames(
    "Barcode")

genus_abundance_dist <- vegdist(
  genus_abundance_matrix,
  method = "bray")

genus_abundance_pcoa <- cmdscale(
  genus_abundance_dist, eig = TRUE, k = 2)

genus_abundance_variance <- round(
  genus_abundance_pcoa$eig /
    sum(genus_abundance_pcoa$eig[genus_abundance_pcoa$eig > 0]) *100, 1)

genus_abundance_scores <- data.frame(
  Barcode = rownames(
    genus_abundance_pcoa$points),
  PCoA1 = genus_abundance_pcoa$points[, 1],
  PCoA2 = genus_abundance_pcoa$points[, 2]) |>
  left_join(
    eDNA_Sampling,
    by = "Barcode") |>
  mutate(
    Location = factor(Location),
    Sample_Type = factor(Sample_Type),
    Habitat = factor(Habitat),
    Trees = factor(Trees),
    Sub_Habitat = factor(Sub_Habitat))

#### genus abundance PCoA ####
ggplot(
  genus_abundance_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Dist_A9_group,
    fill = Dist_A9_group)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linewidth = 0.75) +
  geom_point(
    size = 3) +
  theme_bw() +
  labs(
    colour = "Dist_A9_group",
    fill = "Dist_A9_group",
    x = paste0("PCoA1 (",genus_abundance_variance[1],"%)"),
    y = paste0("PCoA2 (",genus_abundance_variance[2],"%)"))



# Genus P/A PCoA ----------------------------------------------------------


#### making the PCoA ####
genus_pa_matrix <- eDNA_genus_pa |>
  column_to_rownames(
    "Barcode")

genus_pa_dist <- vegdist(
  genus_pa_matrix,
  method = "bray")

genus_pa_pcoa <- cmdscale(
  genus_pa_dist, eig = TRUE, k = 2)

genus_pa_variance <- round(
  genus_pa_pcoa$eig / sum(genus_pa_pcoa$eig[genus_pa_pcoa$eig > 0]) *100, 1)

genus_pa_scores <- data.frame(
  Barcode = rownames(
    genus_pa_pcoa$points),
  PCoA1 = genus_pa_pcoa$points[, 1],
  PCoA2 = genus_pa_pcoa$points[, 2]) |>
  left_join(
    eDNA_Sampling,
    by = "Barcode") |>
  mutate(
    Location = factor(Location),
    Sample_Type = factor(Sample_Type),
    Habitat = factor(Habitat),
    Trees = factor(Trees),
    Sub_Habitat = factor(Sub_Habitat))

#### genus p/a PCoA plot ####
genus_pa_distA9_PCoA <- ggplot(genus_pa_scores,
  aes(x = PCoA1,y = PCoA2,
    colour = Dist_A9_group,
    fill = Dist_A9_group)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linewidth = 0.75) +
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("magenta3", "gold2","cyan3"))+
  scale_colour_manual(values = c("magenta3","gold2","cyan3"))+
  labs(colour = "Distance to A9", fill = "Distance to A9",
    x = paste0("PCoA1 (",genus_pa_variance[1],"%)"),
    y = paste0("PCoA2 (",genus_pa_variance[2],"%)"))

genus_pa_distA9_PCoA

ggsave("genus_pa_distA9_PCoA.png",
       genus_pa_distA9_PCoA,
       width = 8,
       height = 6,
       dpi = 300)


# Taxonomic Assignment Summary --------------------------------------------


assignment_summary <- otu_soil |>
  distinct(OTU_ID, .keep_all = TRUE) |>
  group_by(Assignment_level) |>
  summarise(n_OTUs = n(), total_reads = sum(OTU_size, na.rm = TRUE),
    .groups = "drop") |>
  arrange(desc(total_reads))

assignment_summary



# Summary Stats -----------------------------------------------------------

total_hits <- sum(otu_soil$OTU_size, na.rm = TRUE)
total_hits


n_samples <- vsearch_all |>
  filter(!Barcode %in% c("12", "4")) |>
  summarise(n = n_distinct(Barcode)) |>
  pull(n)

n_samples

# mean hits per soil sample
mean_hits_per_sample <- otu_soil |>
  group_by(Barcode) |>
  summarise(hits = sum(OTU_size, na.rm = TRUE)) |>
  summarise(mean_hits = mean(hits)) |>
  pull(mean_hits)

mean_hits_per_sample






soil_species_pa <- eDNA_species_pa |>
  filter(!Barcode %in% c("12","4"))

soil_richness <- soil_species_pa |>
  mutate(richness = rowSums(across(-Barcode))) |>
  select(Barcode, richness)

mean_richness <- mean(soil_richness$richness)

mean_richness


soil_richness |> slice_max(richness, n = 1)
soil_richness |> slice_min(richness, n = 1)

# mean OTUs per soil sample
mean_otus_per_sample <- otu_soil |>
  filter(!Barcode %in% c("12","4")) |>
  group_by(Barcode) |>
  summarise(n_otus = n_distinct(OTU_ID)) |>
  summarise(mean_otus = mean(n_otus)) |>
  pull(mean_otus)
mean_otus_per_sample


# OTU singleton stats (OTU size = 1)
singleton_stats <- otu_soil |>
  distinct(OTU_ID, .keep_all = TRUE) |>
  summarise(
    total_otus = n(),
    singletons = sum(OTU_size == 1),
    percent_singletons = 100 * singletons / total_otus
  )

singleton_stats



# eDNA Diversity Indices --------------------------------------------------

eDNA_diversity <- eDNA_species |>
  column_to_rownames("Barcode") |>
  as.data.frame()

eDNA_shannon_simpson <- tibble(
  Barcode = rownames(eDNA_diversity),
  eDNA_Shannon = vegan::diversity(eDNA_diversity,index = "shannon"),
  eDNA_Simpson = vegan::diversity(eDNA_diversity,index = "simpson"),
  eDNA_Species_Richness = specnumber(eDNA_diversity))



eDNA_diversity_location <- eDNA_shannon_simpson |>
  left_join(eDNA_Sampling |>
      select(Barcode, Location), by = "Barcode") |>
  filter(!Barcode %in% c("12", "4"))


eDNA_location_summary <- eDNA_diversity_location |>
  select(
    Location,
    eDNA_Shannon,
    eDNA_Simpson,
    eDNA_Species_Richness)

eDNA_location_summary


# eDNA Diversity & Richness Plots -----------------------------------------

ggplot(eDNA_location_summary, aes(x = eDNA_Shannon)) +
  geom_histogram(
    binwidth = 0.25,
    boundary = 0,
    colour = "black",
    fill = "mediumorchid3") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 11, 1))+
  scale_x_continuous(breaks = seq(0, 2.25, 0.25))+
  labs(
    x = "eDNA Shannon Diversity Index",
    y = "Number of locations")



# Stats -------------------------------------------------------------------

# adding environmental data - A9 distance
eDNA_diversity_location_A9 <- eDNA_diversity_location |>
  left_join(
    genus_pa_scores |> select(Barcode, Dist_A9),
    by = "Barcode"
  )

# linear models for both diversity indices
summary(lm(eDNA_Shannon ~ Dist_A9, eDNA_diversity_location_A9))
summary(lm(eDNA_Simpson ~ Dist_A9, eDNA_diversity_location_A9))

# plotting diversity indices against proximity to A9
simpson_a9_edna <- eDNA_diversity_location_A9 |> 
  ggplot(aes(x = Dist_A9, y= eDNA_Simpson))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x = "Distance to A9", y = "Simpson Index", subtitle = "B")

simpson_a9_edna

shannon_a9_edna <- eDNA_diversity_location_A9 |> 
  ggplot(aes(x = Dist_A9, y= eDNA_Shannon))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x = "Distance to A9", y = "Shannon Index", subtitle = "A")

shannon_a9_edna

# visualising plots side by side
shannon_a9_edna | simpson_a9_edna


# read file with species & associated Orders
order_taxonomy <- read_xlsx("Genus_Orders.xlsx")


genus_species_counts <- species_summary |>
separate_wider_delim(Species,
                     delim = "_",
                     names = c("Genus", "species")) |>
  count(Genus, name = "n_species")

genus_species_counts



genus_species_counts <- genus_species_counts |>
  left_join(order_taxonomy, by = "Genus")

# bar chart showing distribution of species within orders and genera
genus_species_order_plot <- ggplot(genus_species_counts,
       aes(x = reorder(Genus, n_species),
           y = n_species,
           fill = Order)) +
  geom_col() +
  coord_flip() +
  labs(
  x = "Genus",
  y = "Number of Species") +
  scale_fill_brewer(palette="Set3")+
  theme_bw() +
  theme(legend.position = "right")

genus_species_order_plot


ggsave("genus_species_order_plot.png",
       genus_species_order_plot,
       width = 8,
       height = 3,
       dpi = 300)
