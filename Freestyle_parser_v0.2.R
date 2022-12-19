#!/usr/bin/env Rscript



# #___ dummy data for testing
# setwd("/media/lucaz/DATA/Seafile/Laboratorio/myscripts/GitHub/Freestyle_parser")
# raw_data <- "toy_data"
# # raw_data <- "toy_data/Output software @IgG.xlsx"
# annotation_db <- "annotation_db.csv"
# config_file <- "config_file.tsv"
# #___ for testing



#! parse command line argument - python style
if (! "optparse" %in% rownames(installed.packages())) {
  cat("Installing missing package optparse\n\n")
  install.packages("optparse")
}
suppressMessages(suppressWarnings(library("optparse")))

option_list <- list(
  make_option(c("-r", "--raw_data"), type="character", default=NULL, 
              help="path to a folder containing multiple Excel files or to a single Excel file", metavar="character"),
  make_option(c("-a", "--annotation_db"), type="character", default=NULL, 
              help="path to a 2-columns *.csv* file containing the list of glycoform names and related glycan masses", metavar="character"),
  make_option(c("-c", "--config_file"), type="character", default=NULL, 
              help="path to a 2-columns *.csv* file containing tolerance_value, overlap_label_range, and any analyzed peptide name with related mass", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



#! parsing function
Freestyle_parser <- function(raw_data, annotation_db, config_file) {
  
  ##! libraries
  required_libs <- c("readxl","tidyverse","ggplot2","ggrepel")
  missing_libs <- required_libs[!required_libs %in% rownames(installed.packages())]
  if (length(missing_libs) > 1) {
    cat("Installing missing packages: ", paste(missing_libs, collapse = ", "), "\n\n")
    sapply(missing_libs, function(i) install.packages(i))
  }
  invisible(suppressMessages(suppressWarnings(lapply(required_libs, library, character.only = TRUE))))
  
  
  
  ##! read-in raw data
  if (file_test("-d", raw_data)) {
    raw_tables <- list.files(raw_data, pattern = ".xlsx$|.xls$", 
                             full.names = T, recursive = T)
    raw_table <- lapply(raw_tables, function(i) {
      if (grepl(".xlsx$", i)) {
        i_tab <- read_xlsx(i)
        i_name <- gsub(".xls$|.xlsx$", "", basename(i))
        i_name <- gsub(" ", "_", i_name)
        i_target_peptide <- gsub(".*@", "", i_name)
        # ??? '@' can be changed, or even be specified e.g. in config file. But target peptide name must be contained between fixed and unique delimiter (e.g. @ and file extension)
        i_tab <- data.frame(Sample=i_name, Target.Peptide=i_target_peptide, i_tab) 
        
      } else {
        i_tab <- read_xls(i)
        i_name <- gsub(".xls$|.xlsx$", "", basename(i))
        i_name <- gsub(" ", "_", i_name)
        i_target_peptide <- gsub(".*@", "", i_name)
        i_tab <- data.frame(Sample=i_name, Target.Peptide=i_target_peptide, i_tab) 
        
      }
      
    })
    raw_table <- do.call(rbind, raw_table)
    
  } else if (file_test("-f", raw_data)) {
    if (grepl(".xlsx$", raw_data)) {
      raw_table <- read_xlsx(raw_data)
      tab_name <- gsub(".xls$|.xlsx$", "", basename(raw_data))
      tab_name <- gsub(" ", "_", tab_name)
      target_peptide <- gsub(".*@", "", tab_name)
      raw_table <- data.frame(Sample=tab_name, Target.Peptide=target_peptide, raw_table) 
      
    } else if (grepl(".xls$", raw_data)) {
      raw_table <- read_xls(raw_data)
      tab_name <- gsub(".xls$|.xlsx$", "", basename(raw_data))
      tab_name <- gsub(" ", "_", tab_name)
      target_peptide <- gsub(".*@", "", tab_name)
      raw_table <- data.frame(Sample=tab_name, Target.Peptide=target_peptide, raw_table) 
      
    } else {
      stop(cat("The value provided for 'raw_data' (", raw_data,
               ") doesn't point to a valid directory nor input file"))
      
    }
    
  } else {
    stop(cat("The value provided for 'raw_data' (", raw_data,
             ") doesn't point to a valid directory nor input file"))
    
  }
  
  
  
  ##! read-in config file
  conf_file <- read.csv(config_file, h=F, row.names = 1) %>% 
    t() %>% as.data.frame()
  
  if (colnames(conf_file)[1] != "tolerance_value") {
    stop("Config file is missing the parameter 'tolerance_value'. Add a row like:\n\n",
    "tolerance_value\t5\n\n\n")
  }
  
  if (colnames(conf_file)[2] != "overlap_label_range") {
    stop("Config file is missing the parameter 'overlap_label_range'. Add a row like:\n\n",
         "overlap_label_range\t50\n\n\n")
  }

  missing_peptide <- unique(raw_table$Target.Peptide[!raw_table$Target.Peptide %in% names(conf_file)])
  if (length(missing_peptide) > 0) {
    stop("In the config file, the following peptide masses are missing:\n\n", 
         paste(missing_peptide, collapse = ", "),
         "\n\nPlease add the missing entries required for processing the following samples:\n\n", 
         paste(unique(raw_table$Sample[raw_table$Target.Peptide %in% missing_peptide]), collapse = ", "),
         "\n\n\n")
  }
  
  
  
  ##! read-in annotation database
  ann_db <- read.csv(annotation_db, h=T)
  colnames(ann_db) <- c("Glycoform", "Glycan_Mass") # to ensure compliance with downstream code
  
  # strick check (too strickt?) as, in the annotaion step, the ppm is calculated on 'Monoisotopic.Mass' and not on 'Glycan_Mass'
  Glycoform_issues <- ann_db %>%
    mutate(glycan_ppm = Glycan_Mass/1E6*conf_file$tolerance_value,
           upper_threshold = Glycan_Mass + glycan_ppm,
           lower_threshold = Glycan_Mass - glycan_ppm
    ) %>% 
    rowwise() %>% mutate(similar_Glycoform = with(ann_db[ann_db$Glycoform != Glycoform, ], 
                                                  paste0(Glycoform[Glycan_Mass >= lower_threshold & 
                                                                     Glycan_Mass <= upper_threshold], collapse = ", "))) %>%
    filter(similar_Glycoform != "")
  
  if (nrow(Glycoform_issues) > 0) {
    stop(paste0("In the list of annotated Glycoforms, there are entries with a mass difference",
                " lower that the tollerance value provided (", conf_file$tolerance_value,"ppm).\n",
                "I.e., it is not possible to correctly discriminate between them.\n",
                "Affected Glycoforms are:\n\n"),
         apply(Glycoform_issues, 1, function(i) paste(i[1], "\twith\t", i[6], "\n")),
         "\n\n\n")
  }
  
  
  
  ##! filter table
  raw_table_filt <- raw_table %>%
    filter(!is.na(No.)) %>% # remove peak details
    mutate(across(3:10, as.numeric)) %>%
    select(Sample, Target.Peptide, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity) # keep only desired cols
  
  
  
  ##! annotate peaks
  peaks_ann <- raw_table_filt %>%
    mutate(target_pept_mass = as.numeric(conf_file[Target.Peptide]),
           measurement_tolerance = Monoisotopic.Mass/1E6*conf_file$tolerance_value # ??? ppm should be calculated on the Monoisotopic.Mass, right? as it accounts for the instrument measurement error
    ) %>%
    rowwise() %>% mutate(Glycoform = {
      ###! look if peak's mass matches the mass of the target peptide + a glycan
      i_ann = ann_db$Glycoform[target_pept_mass + ann_db$Glycan_Mass >= Monoisotopic.Mass - measurement_tolerance &
                                 target_pept_mass + ann_db$Glycan_Mass <= Monoisotopic.Mass + measurement_tolerance]
      
      if (identical(i_ann, character(0))) {
        ###! look if peak's mass matches the mass of the not-glycosylated target peptide 
        if (target_pept_mass >= Monoisotopic.Mass - measurement_tolerance &
            target_pept_mass <= Monoisotopic.Mass + measurement_tolerance) { "not-glycosylated" 
        } else { NA } ###! peaks without annotation
        
      } else { i_ann }
      
    }) %>%
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>% 
    mutate(Rel.Abundance_ann = if_else(!is.na(Glycoform), Sum.Intensity/sum(Sum.Intensity[!is.na(Glycoform)])*100, NA_real_),
           Rel.Abundance_highest = Sum.Intensity/max(Sum.Intensity)
           )
  
  
  ##! output tables
  peaks_ann_all <- peaks_ann %>%
    select(Sample, Target.Peptide, Glycoform, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity)
  write.table(peaks_ann_all, "Peaks_annotated_all.tsv",
              col.names = T, row.names = F, quote = F, sep = "\t")
  
  peaks_ann_clean <- peaks_ann %>%
    filter(!is.na(Glycoform)) %>%
    select(Sample, Target.Peptide, Glycoform, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity, Rel.Abundance_ann)
  write.table(peaks_ann_clean, "Peaks_annotated_clean.tsv",
              col.names = T, row.names = F, quote = F, sep = "\t")
  
  
  
  ##! plot from clean table
  min_mass_diff <- conf_file$overlap_label_range
  z <- 1
  peaks_ann_extra <- peaks_ann %>%
    ##! find peak clusters - consider only annotated peaks!
    filter(!is.na(Glycoform)) %>%
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>%
    mutate(too_close_sx = if_else(abs(Monoisotopic.Mass - lag(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
           too_close_dx = if_else(abs(Monoisotopic.Mass - lead(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
    ) %>%
    rowwise() %>% mutate(
      dummy_cluster = {
        if ((too_close_sx == F & too_close_dx == T) | 
            all(too_close_sx, too_close_dx)) {
          z
        } else {
          z <- z + 1
          z-1
        } 
      }) %>% 
    group_by(Sample, dummy_cluster) %>% 
    mutate(drop_label = ifelse(Sum.Intensity == max(Sum.Intensity), F, T),
           label_filt = ifelse(drop_label == T, NA, Glycoform),
           mass_filt = ifelse(drop_label == T, NA, round(Monoisotopic.Mass, 3)),
           label_filt_wmass = ifelse(drop_label == T, NA, paste0(label_filt, "\n", sprintf("%0.3f", mass_filt)))) %>%
    ungroup() %>%
    select(-c(too_close_sx, too_close_dx, dummy_cluster)) %>%
    ###! merge back with full table
    left_join(peaks_ann, .) %>%
    replace_na(list(drop_label = T)) %>%
    ###! find peak clusters - any peaks!
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>%
    mutate(too_close_sx2 = if_else(abs(Monoisotopic.Mass - lag(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
           too_close_dx2 = if_else(abs(Monoisotopic.Mass - lead(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
    ) %>%
    rowwise() %>% mutate(
      dummy_cluster2 = {
        if ((too_close_sx2 == F & too_close_dx2 == T) | 
            all(too_close_sx2, too_close_dx2)) {
          z
        } else {
          z <- z + 1
          z-1
        } 
      }) %>% 
    group_by(Sample, dummy_cluster2) %>% 
    mutate(drop_label2 = ifelse(Sum.Intensity == max(Sum.Intensity) &
                                  drop_label != F, F, T),
           mass_filt2 = ifelse(drop_label2 == T, NA, round(Monoisotopic.Mass, 3)))
  
  p1 <- peaks_ann_extra %>% filter(!is.na(Glycoform)) %>%
    ggplot(aes(Monoisotopic.Mass, Rel.Abundance_highest)) +
    geom_segment(aes(xend = Monoisotopic.Mass, yend = 0), linewidth = 1, lineend = "butt") +
    geom_label_repel(aes(label = label_filt_wmass), size = 2.3, force = 5,
                     direction = "y", nudge_y = .1,
                     force_pull = 0, # do not pull toward data points
                     min.segment.length = 0, segment.color = "red", 
                     segment.linetype = 2, segment.size = 0.3) +
    coord_cartesian(clip = "off") + 
    # scale_y_continuous(expand = expansion(mult = c(0, 1))) +
    facet_grid(Sample~., scales = "free") +
    labs(y = "Intensity (realtive to highest peak)") +
    theme_minimal() + theme(text = element_text(size = 15))
  ggsave("Peaks_annotated_clean.pdf", p1, dpi = 300, 
         width = 16, height = 3 * length(unique(peaks_ann_extra$Sample)))
  
  p2 <- ggplot(peaks_ann_extra, aes(Monoisotopic.Mass, Rel.Abundance_highest)) +
    geom_segment(aes(xend = Monoisotopic.Mass, yend = 0), linewidth = 1, lineend = "butt") +
    geom_label_repel(aes(label = mass_filt2), size = 1.3, force = 5,
                     direction = "y", nudge_y = .1,
                     force_pull = 0, # do not pull toward data points
                     min.segment.length = 0, segment.color = "red", 
                     segment.linetype = 2, segment.size = 0.3) +
    geom_label_repel(aes(label = label_filt_wmass), size = 2.3, force = 5,
                     direction = "y", nudge_y = .1,
                     force_pull = 0, # do not pull toward data points
                     min.segment.length = 0, segment.color = "red", 
                     segment.linetype = 2, segment.size = 0.3) +
    coord_cartesian(clip = "off") + 
    # scale_y_continuous(expand = expansion(mult = c(0, 1))) +
    facet_grid(Sample~., scales = "free") +
    xlim(500, 6500) + #??? restrict window???
    labs(y = "Intensity (realtive to highest peak)") +
    theme_minimal() + theme(text = element_text(size = 15))
  ggsave("Peaks_annotated_all.pdf", p2, dpi = 300, 
         width = 16, height = 3 * length(unique(peaks_ann_extra$Sample)))
  
}



#! execute
Freestyle_parser(raw_data = opt$raw_data, 
                 annotation_db = opt$annotation_db, 
                 config_file = opt$config_file)