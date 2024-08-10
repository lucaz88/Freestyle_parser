#!/usr/bin/env Rscript



# #___ dummy data for testing
# setwd("/media/lucaz/DATA/Seafile/Laboratorio/myscripts/GitHub/Freestyle_parser")
# # raw_data <- "toy_data"
# # raw_data <- "toy_data/Output software @IgG.xlsx"
# raw_data <- "toy_adducts"
# annotation_db <- "annotation_db.csv"
# # config_file <- "config_file.csv"
# config_file <- "config_file_adducts.csv"
# out_dir <- "res_adducts"
# #___ for testing



#! load libraries
required_libs <- c("optparse","readxl","tidyverse","ggplot2","ggrepel")
invisible(lapply(required_libs, library, character.only = TRUE))




#! parse command line argument - python style
option_list <- list(
  make_option(c("-r", "--raw_data"), type="character", default=NULL,
              help="path to a folder containing multiple Excel files or to a single Excel file", metavar="character"),
  make_option(c("-c", "--config_file"), type="character", default=NULL, 
              help="path to a 2-columns *.csv* file containing tolerance_value, overlap_label_range, and any analyzed peptide name with related mass", metavar="character"),
  make_option(c("-a", "--annotation_db"), type="character", default=NULL, 
              help="path to a 2-columns *.csv* file containing the list of glycoform names and related glycan masses", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="path to a folder for the output files", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



#! parsing function
Freestyle_parser <- function(raw_data, config_file, annotation_db, out_dir) {
  
  dir.create(out_dir, recursive = T, showWarnings = F)
  
  
  
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
               ") doesn't point to a valid directory nor input file\n\n\n"))
      
    }
    
  } else {
    stop(cat("The value provided for 'raw_data' (", raw_data,
             ") doesn't point to a valid directory nor input file\n\n\n"))
    
  }
  
  
  
  ##! read-in config file
  conf_file <- read.csv(config_file, h=F, row.names = 1) %>% 
    t() %>% as.data.frame()
  
  if (colnames(conf_file)[1] != "tolerance_value") {
    stop("Config file is missing the parameter 'tolerance_value'. Add a row such as:\n\n",
         "tolerance_value,5\n\n\n")
  }
  
  if (colnames(conf_file)[2] != "overlap_label_range") {
    stop("Config file is missing the parameter 'overlap_label_range'. Add a row such as:\n\n",
         "overlap_label_range,50\n\n\n")
  }
  
  if (colnames(conf_file)[3] != "mass_range") {
    stop("Config file is missing the parameter 'mass_range'. Add a row such as:\n\n",
         "mass_range,1000:NA\n\n\n")
  }
  
  if (colnames(conf_file)[4] != "modification_name") {
    stop("Config file is missing the parameter 'modification_name'. Add a row such as:\n\n",
         "modification_name,NH4_adduct\n\n\n")
  }
  
  if (colnames(conf_file)[5] != "modification_mass") {
    stop("Config file is missing the parameter 'modification_mass'. Add a row such as:\n\n",
         "modification_mass,18.03382\n\n\n")
  }
  
  conf_file <- conf_file %>% 
    mutate_at(vars(-c("mass_range", "modification_name")), as.numeric)
  
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
      
    })
  
  
  
  ##! annotate adducts
  if (!is.na(conf_file$modification_name)) {
    cat("\n\nAdduct annotation requested.", 
        "\nNB: it is possible to check for one type of adduct (occurring once) at the time.\n",
        "Repeat the analysis editing the config file to check for different adducts.\n\n")
    
    ###! create adducts table
    ann_db <- ann_db %>%
      mutate(Adduct = paste0(Glycoform, "_", conf_file$modification_name),
             Adduct_Mass = Glycan_Mass + conf_file$modification_mass)
    
    ###! run annotation
    peaks_ann <- peaks_ann %>%
      rowwise() %>% mutate(Adduct = {
        ###! look if peak's mass matches the mass of the target peptide + a glycan
        i_adduct = ann_db$Adduct[target_pept_mass + ann_db$Adduct_Mass >= Monoisotopic.Mass - measurement_tolerance &
                                   target_pept_mass + ann_db$Adduct_Mass <= Monoisotopic.Mass + measurement_tolerance]
        
        if (identical(i_adduct, character(0))) { NA ###! peaks that are not adducts
        } else { i_adduct }
        
      }) %>%
      group_by(Sample) %>% 
      rowwise() %>% mutate(same_peak = {
        i_peak = as.character(gsub(paste0("_", conf_file$modification_name), "",
                                   na.omit(c(Glycoform, Adduct))))
        if (length(i_peak) > 0) { i_peak } else { NA }
      })
    
  }
  
  
  
  ##! calculate peak abundance
  peaks_ann <- peaks_ann %>%
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>% 
    mutate(
      Rel.Abundance_ann = if_else(!is.na(Glycoform), Sum.Intensity/sum(Sum.Intensity[!is.na(Glycoform)])*100, NA_real_),
      Rel.Abundance_highest = Sum.Intensity/max(Sum.Intensity)*100
    ) %>% 
    arrange(Monoisotopic.Mass, .by_group = T)
  
  if (!is.na(conf_file$modification_name)) {
    peaks_ann <- peaks_ann %>% 
      group_by(Sample) %>%
      mutate(
        tot_smpl_rel_abund = sum(Rel.Abundance_highest[!is.na(same_peak)])
      ) %>%
      group_by(Sample, same_peak) %>%
      mutate(
        Rel.Abundance_w.adduct = if_else(!is.na(same_peak), sum(Rel.Abundance_highest)/tot_smpl_rel_abund*100, NA_real_),
      ) %>% 
      group_by(Sample) %>% arrange(Monoisotopic.Mass, .by_group = T)
  }
  
  
  
  ##! output tables
  if (is.na(conf_file$modification_name)) {
    peaks_ann_all <- peaks_ann %>%
      select(Sample, Target.Peptide, Glycoform, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity, Rel.Abundance_ann)
    write.table(peaks_ann_all, file.path(out_dir, "Peaks_annotated_all.tsv"),
                col.names = T, row.names = F, quote = F, sep = "\t")
    
    peaks_ann_clean <- peaks_ann %>%
      filter(!is.na(Glycoform)) %>%
      select(Sample, Target.Peptide, Glycoform, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity, Rel.Abundance_ann)
    write.table(peaks_ann_clean, file.path(out_dir, "Peaks_annotated_clean.tsv"),
                col.names = T, row.names = F, quote = F, sep = "\t")
    
  } else {
    peaks_ann_all <- peaks_ann %>%
      select(Sample, Target.Peptide, Glycoform, Adduct, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity, Rel.Abundance_ann, Rel.Abundance_w.adduct)
    write.table(peaks_ann_all, file.path(out_dir, "Peaks_annotated_all.tsv"),
                col.names = T, row.names = F, quote = F, sep = "\t")
    
    peaks_ann_clean <- peaks_ann %>%
      filter(!is.na(same_peak)) %>%
      select(Sample, Target.Peptide, Glycoform, Adduct, Monoisotopic.Mass, Number.of.Charge.States, RT.Range, Sum.Intensity, Rel.Abundance_ann, Rel.Abundance_w.adduct)
    write.table(peaks_ann_clean, file.path(out_dir, "Peaks_annotated_clean.tsv"),
                col.names = T, row.names = F, quote = F, sep = "\t")
  }
  
  
  
  ##! find peak clusters - consider only annotated peaks!
  min_mass_diff <- conf_file$overlap_label_range
  peaks_ann_extra <- peaks_ann %>%
    filter(!is.na(Glycoform)) %>%
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>%
    mutate(too_close_sx = if_else(c(min_mass_diff, diff(Monoisotopic.Mass)) < min_mass_diff, T, F),
           # too_close_dx = if_else(abs(Monoisotopic.Mass - lead(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
           dummy_cluster = if_else(!too_close_sx, 1, 0),
           dummy_cluster = cumsum(dummy_cluster)
           )%>% 
    group_by(Sample, dummy_cluster) %>% 
    mutate(drop_label = ifelse(Sum.Intensity == max(Sum.Intensity), F, T),
           label_filt = ifelse(drop_label == T, NA, Glycoform),
           mass_filt = ifelse(drop_label == T, NA, round(Monoisotopic.Mass, 3)),
           label_filt_wmass = ifelse(drop_label == T, NA, paste0(label_filt, "\n", sprintf("%0.3f", mass_filt)))
           ) %>%
    ungroup() %>%
    select(-c(too_close_sx, dummy_cluster)) %>%
    ###! merge back with full table
    left_join(peaks_ann, .) %>%
    replace_na(list(drop_label = T)) %>%
    ###! find peak clusters - any peaks!
    group_by(Sample) %>% 
    arrange(Monoisotopic.Mass, .by_group = T) %>%
    mutate(too_close_sx2 = if_else(c(min_mass_diff, diff(Monoisotopic.Mass)) < min_mass_diff, T, F),
           # too_close_dx2 = if_else(abs(Monoisotopic.Mass - lead(Monoisotopic.Mass, default = 0)) <= min_mass_diff, T, F),
           dummy_cluster2 = if_else(!too_close_sx2, 1, 0),
           dummy_cluster2 = cumsum(dummy_cluster2)
    ) %>%
    group_by(Sample, dummy_cluster2) %>% 
    mutate(drop_label2 = ifelse(Sum.Intensity == max(Sum.Intensity) &
                                  drop_label != F, F, T),
           mass_filt2 = ifelse(drop_label2 == T, NA, round(Monoisotopic.Mass, 3)))
  
  
  
  ##! plot from full table
  if (is.na(conf_file$modification_name)) {
    p1 <- ggplot(peaks_ann_extra, aes(Monoisotopic.Mass, Rel.Abundance_highest)) + 
      geom_segment(aes(xend = Monoisotopic.Mass, yend = 0), 
                   linewidth = 1, lineend = "butt")
  } else {
    p1 <- ggplot(peaks_ann_extra, aes(Monoisotopic.Mass, Rel.Abundance_highest)) + 
      geom_segment(aes(xend = Monoisotopic.Mass, yend = 0, color = !is.na(Adduct)), 
                   linewidth = 1, lineend = "butt") +
      scale_color_manual(values = c("black", "orange")) +
      labs(color = "Adduct")
  }
  p1 <- p1 + geom_label_repel(aes(label = mass_filt2), size = 1.3, force = 5,
                              direction = "y", nudge_y = .1, na.rm=T,
                              force_pull = 0, # do not pull toward data points
                              min.segment.length = 0, segment.color = "red", 
                              segment.linetype = 2, segment.size = 0.3) +
    geom_label_repel(aes(label = label_filt_wmass), size = 2.3, force = 5,
                     direction = "y", nudge_y = .1, na.rm=T,
                     force_pull = 0, # do not pull toward data points
                     min.segment.length = 0, segment.color = "red", 
                     segment.linetype = 2, segment.size = 0.3) +
    coord_cartesian(clip = "off") + 
    # scale_y_continuous(expand = expansion(mult = c(0, 1))) +
    facet_wrap(~Sample, ncol = 1, scales = "free") +
    # ylim(0, max(peaks_ann_extra$Rel.Abundance_highest[!is.na(peaks_ann_extra$Glycoform)])) + #??? not working with faceted plots, cut out labels
    xlim(as.numeric(gsub(":.*", "", conf_file$mass_range)),
         as.numeric(gsub(".*:", "", conf_file$mass_range))) +
    labs(y = "Intensity (realtive to highest peak)") +
    theme_minimal() + theme(text = element_text(size = 15), legend.position = "top")
  ggsave(file.path(out_dir, "Peaks_annotated_all.pdf"), p1, dpi = 300, units = "cm", limitsize = F, scale = 2,
         width = 20, height = 5.5 * length(unique(peaks_ann_extra$Sample)) + .5)
  
  
  
  ##! plot from clean table
  if (is.na(conf_file$modification_name)) {
    p2 <- peaks_ann_extra %>% filter(!is.na(Glycoform)) %>%
      ggplot(aes(Monoisotopic.Mass, Rel.Abundance_highest)) + 
      geom_segment(aes(xend = Monoisotopic.Mass, yend = 0), 
                   linewidth = 1, lineend = "butt")
  } else {
    p2 <- peaks_ann_extra %>% filter(!is.na(same_peak)) %>%
      ggplot(aes(Monoisotopic.Mass, Rel.Abundance_highest)) + 
      geom_segment(aes(xend = Monoisotopic.Mass, yend = 0, color = !is.na(Adduct)), 
                   linewidth = 1, lineend = "butt") +
      scale_color_manual(values = c("black", "orange")) +
      labs(color = "Adduct")
  }
  p2 <- p2 + geom_label_repel(aes(label = label_filt_wmass), size = 2.3, force = 5,
                              direction = "y", nudge_y = .1, na.rm=T,
                              force_pull = 0, # do not pull toward data points
                              min.segment.length = 0, segment.color = "red", 
                              segment.linetype = 2, segment.size = 0.3) +
    coord_cartesian(clip = "off") +
    # scale_y_continuous(expand = expansion(mult = c(0, 1))) +
    facet_wrap(~Sample, ncol = 1, scales = "free") +
    labs(y = "Intensity (realtive to highest peak)") +
    theme_minimal() + theme(text = element_text(size = 15), legend.position = "top")
  ggsave(file.path(out_dir, "Peaks_annotated_clean.pdf"), p2, dpi = 300, units = "cm", limitsize = F, scale = 2,
         width = 16, height = 3.5 * length(unique(peaks_ann_extra$Sample)) + .5)
  
}



#! execute
Freestyle_parser(raw_data = opt$raw_data, 
                 config_file = opt$config_file,
                 annotation_db = opt$annotation_db, 
                 out_dir = opt$out_dir)