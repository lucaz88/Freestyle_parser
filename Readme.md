# R function to parse Freestyle Excel outputs

## Required inputs:
* ***raw_data***: path to a folder containing multiple Excel files (e.g. Project_01/Freestyle_outputs), or path pointing to a single Excel file (e.g. Project_01/Freestyle_outputs/Outputfile[]()@IgM_site4.xlsx). Note that the target peptide name has to be appended **unambiguously** to the filename of each Excel file (e.g. *@peptide_X.xls*)
* ***annotation_db***: path to a 2-columns *.csv* file containing the list of glycoform names and related glycan masses. Note that the first row in the file contain the column names. For the glycoform names, use '_' instead of spaces, and '@' to separate name types, e.g. "H5N2@Man5"

    File example:

       glycoform,Glycan_Mass
       N1@single-GlcNAc_(Gn),203.0794
       H2N2@MU_,730.2644
       H2N2X@MUX_,862.3067
       H3N2@MM_,892.3172

* ***config_file***: path to a 2-columns *.csv* file containing:
    * 1st row - 'tolerance_value' followed by the value (in ppm) setting the stringency of peaks' annotation
    * 2nd row - 'overlap_label_range' followed by the value (in Dalton) defining a range within which only the highest intensity peak will be labbeled in the plot
    * 3nd row - 'mass_range' followed by the values' range (in Dalton:Dalton) of min and max Monoisotopic mass that will be plotted. Use NA to allows for automatic selection of the value/s (e.g. 750:NA)
    * *n*rows containing names of any analyzed peptides followed by their masses (in Dalton)
    
    File example:
    
       tolerance_value,5
       overlap_label_range,50
       mass_range,NA:NA
       IgG,1189.512
       IgM_site1,1284.618

## Set-up the environment:
* R needs to be installed on the computer. Please follow the instruction provided under *Download and Install R* at https://cran.r-project.org/
* Install the necessary packages:
    * Start an R session
    * Type in the following code and follow the prompted instructions:

          required_libs <- c("optparse","readxl","tidyverse","ggplot2","ggrepel")
          missing_libs <- required_libs[!required_libs %in% rownames(installed.packages())]
          if (length(missing_libs) > 1) {
            cat("Installing missing packages: ", paste(missing_libs, collapse = ", "), "\n\n")
            sapply(missing_libs, function(i) install.packages(i))
          }
          invisible(lapply(required_libs, library, character.only = TRUE))

* Fetch the tool repository
    * in **Windows**:
        - Almost at the top of this GitHub page, click *"<> Code"*
        - Click on *"Download ZIP"*
        - Save and extract the file where you prefer (e.g. C:\Users\Luca\Tools)

    * in **Linux**:

        open a terminal and type:

          sudo apt update
          sudo apt install git
          cd /home/user/tools # or to whatever folder in which you want to save the tool
          git clone https://github.com/lucaz88/Freestyle_parser.git

## How to run:

For **Windows** system:

* Press the Windows Start button on the screen or keyboard
* Type in "Command Prompt"
* Left click on Command Prompt
* Move into the repository folder `cd C:\Users\Luca\Tools\Freestyle_parser-main` in which toy files are provided to test the tool
* Type `"C:\Program Files\R\R-4.2.2\bin\Rscript.exe" Freestyle_parser_v0.3.R -r toy_data -a annotation_db.csv -c config_file.csv`
\# update the paths to where you actually installed R and downloaded the repo

For **Linux** system:

open a terminal and type:

    cd /home/user/tools/Freestyle_parser
    Rscript Freestyle_parser_v0.3.R -r toy_data -a annotation_db.csv -c config_file.csv

**To run you own sample simply `cd` into the project folder with your data and provide the related values the strings for `-r`, `-a` and/or `-c` arguments.**