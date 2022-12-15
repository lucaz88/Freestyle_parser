# R function to parse Freestyle Excel outputs

## Input specifications:
* ***raw_data***: path to a folder containing multiple Excel files (e.g. Project_01/Freestyle_outputs), or path pointing to a single Excel file (e.g. Project_01/Freestyle_outputs/Outputfile@IgM_site4.xlsx). Note that the target peptide name has to be appended **unambiguously** to the filename of each Excel file (e.g. *@peptide_X.xls*)
*  ***annotation_db***: path to a 2-columns *.csv* file containing the list of glycoform names and related glycan masses. Note that the first row in the file contain the column names. For the glycoform names, use '_' instead of spaces, and '@' to separate name types e.g. "H5N2@Man5"
Example:

    glycoform,Glycan_Mass
    N1@single-GlcNAc_(Gn),203.0794
    H2N2@MU_,730.2644
    H2N2X@MUX_,862.3067
    H3N2@MM_,892.3172

*  ***config_file***: path to a 2-columns *.csv* file containing:
    *1st row - 'tolerance_value' followed by the value (in ppm) setting the stringency of peaks' annotation
    *2nd row - 'overlap_label_range' followed by the value (in Dalton) defining a range within which only the highest intensity peak will be labbeled in the plot
    **n* rows containing names of any analyzed peptides followed by their masses (in Dalton).
    Example:
    
        tolerance_value,5
        overlap_label_range,50
        IgG,1189.512
        IgM_site1,1284.618

## Setting up the environment:
* R needs to be installed on the computer. Please follow the instruction provided under *Download and Install R* at https://cran.r-project.org/
* fetch the tool repository.
For Windows system:

    ???

For Linux system

open a terminal and type:

    sudo apt update
    sudo apt install git
    cd ~ # or to whatever folder in which you want to save the tool
    git clone https://github.com/lucaz88/Freestyle_parser.git

## How to run:

For Windows system:

    ???

For Linux system:

open a terminal and type:

    Rscript Freestyle_parser_v0.2.R -r toy_data -a annotation_db.csv -c config_file.csv 
