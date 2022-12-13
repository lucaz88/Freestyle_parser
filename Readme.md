# R function to parse Freestyle Excel outputs

## Input specifications:
* raw_data: path to a folder containing multiple Excel files, or to a specific file;
           filenames of each Excel HAS to include, unambiguously, the peptide name: *@pep_name.xls*
*  annotation_db: use '_' instead of spaces, use '@' to separate name types e.g. "H5N2@Man5"
*  config_file: 2-columns tab separated list in which:
              - 1st row - 'tolerance_value' followed by value (in ppm) setting the stringency of peaks' annotation
              - 2nd row - 'overlap_label_range' followed by value (in Dalton) defining range within which only the label of highest intensity peak with be plotted
              - following rows contain the names of target peptides followed by their masses (in Dalton).
                All peptides relevant for the analyzed samples must be included

## How to run: