#' @author Lei Li
#' @title mms2plot
#' @description Visualization of multiple MS/MS spectra for groups of modified 
#' and non-modified peptides
#' @export mms2plot
#' @param id_table_path File path of the table that contains the information of 
#'        peptide-spectrum matches (PSMs) and number labels. The PSM information 
#'        is referred to as the output file msms.txt from Maxquant.
#' @param mod_xml_path File path of the modification file. The Maxquant 
#'        mdification xml file modifications.xml can be directly used.
#' @param par_filepath File path of of the parameter batch table that includes
#'        the parameter xml file and the fragment mass tolerance (ppm). The 
#'        parameter file format is referred to as par.xml in Maxquant.
#' @param output_path Folder Path that stores the output image files.
#' @param min_intensity_ratio minimal percentage threshold of MS2 intensity,
#'        compared with the highest intensity. (default=0.01).
#' @param pdf_width The width of a single PSM figure area in inches. The area
#'        includes both plot region and outer margin area. The default is 3.35
#'        for a single column. The width is 7 for double-column.
#' @param pdf_height The height of a single PSM figure area in inches. The
#'        default is pdf_width/2.4.
#' @param xmai Margin of the figure in number of inches for x axis.
#'        (default=pdf_width*0.15/3.35).
#' @param ymai Margin of the figure in number of inches for y axis.
#'        (default=pdf_width*0.15/3.35).
#' @param ppm The threshold of mass error in parts per million(ppm):
#'        (exactMass-accurateMass)/exactMass*1E6. (default=20).
#' @param y_ion_col y ion color. The default is red.
#' @param b_ion_col b ion color. The default is blue.
#' @param peaks_col color of all M/Z peaks. The default is grey.
#' @param ymax The height of the plot region relative to the default.
#'        (default=1.6).
#' @param info_height The height of MS2 annotation. (default=1.5)
#' @param peptide_height The height of peptide sequence annotation in the plot
#'        region, relative to the default. (default=1.3)
#' @param mod_height The height of modification annotation relative to the
#'        location where peptide sequence is annotation. (default=0.07).
#' @param len_annoSpace The length of b/y ion annotation segments.(default=0.1).
#' @param lwd line width relative to the default. (default=pdf_width/3.35).
#' @param cex A numerical value giving the amount by which plotting text and
#'        symbols is magnified relative to the default.(default=pdf_width/3.35).
#' @param show_letterBY Logical: should "b" or "y" characters are shown on the
#'        peak annotation? The default is FALSE.
#'
#' @return No value is returned.
#' @import xml2
#' @import gsubfn
#' @importFrom MSnbase readMSData
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines par plot segments text
#' @importFrom data.table fread
#' @note See vignettes for more details
#' @examples
#' ###
#' general_path = system.file( package = "mms2plot",dir = "extdata" )
#' setwd( general_path )
#' ###################################
#' # Generate mms2plot for label-free data
#' lf_path = system.file( package = "mms2plot",dir = "extdata/label_free" )
#' # id_table_path expands the Maxqaunt output msms.txt by adding "label" column
#' id_table_path = dir( lf_path, "msms_labelfree.txt", full.names = TRUE )
#' mod_xml_path = dir( general_path, "modifications.xml", full.names = TRUE )
#' # par_filepath contains par.xml with full file path and PPM cutoff
#' par_filepath = dir( general_path, "par_batch.txt", full.names = TRUE )
#' output_path = general_path
#' #mms2plot( id_table_path, mod_xml_path, par_filepath, output_path) #not run
#' 
#' ###################################
#' # Generate mms2plot for TMT labelling
#' TMT_path = system.file( package = "mms2plot",dir = "extdata/TMT" )
#' # id_table_path expands the Maxqaunt output msms.txt by adding "label" column
#' id_table_path = dir( TMT_path, "msms_TMT.txt", full.names = TRUE )
#' mod_xml_path = dir( general_path, "modifications.xml", full.names = TRUE )
#' # par_filepath contains par.xml with full file path and PPM cutoff
#' par_filepath = dir( general_path, "par_batch.txt", full.names = TRUE )
#' output_path = general_path
#' #mms2plot( id_table_path, mod_xml_path, par_filepath, output_path) #not run
#' 
#' #####################################
#' # Generate mms2plot for SILAC labelling
#' SILAC_path = system.file( package = "mms2plot",dir = "extdata/silac" )
#' id_table_path = dir( SILAC_path, "msms_SILAC.txt", full.names = TRUE )
#' mod_xml_path = dir( general_path, "modifications.xml", full.names = TRUE )
#' par_filepath = dir( general_path, "par_batch.txt", full.names = TRUE )
#' output_path = general_path
#' #mms2plot( id_table_path, mod_xml_path,par_filepath,output_path ) #not run
#' 
#' #####################################
#' # Generate mms2plot for dimethyl labelling
#' dim_path = system.file(package="mms2plot",dir="extdata/Dimethyl_Labelling")
#' id_table_path = dir( dim_path, "msms_dim.txt", full.names = TRUE )
#' mod_xml_path = dir( general_path, "modifications.xml", full.names = TRUE )
#' par_filepath = dir( general_path, "par_batch.txt", full.names = TRUE )
#' output_path = general_path
#' #mms2plot(id_table_path, mod_xml_path, par_filepath, output_path) #not run
#' 
#roxygen2::roxygenise()
# rm(list=ls())
# .libPaths( c( .libPaths(), "E:/mms2plot/MMS2plot/inst/extdata") )
# library(xml2)
# library(MSnbase)
# library(data.table)
# library(DescTools)  # MixColor
# library(gsubfn)
# 
# PPM_denominator=1E6
# 
# source("R/plot_mirror_or_group.R")
# source("R/add_mod_aa.R")
# source("R/psm_calculation.R")
# source("R/plot_components.R")
# 
# setwd("e:/lei_package4/MMS2plot/")

mms2plot <- function(id_table_path,
                    mod_xml_path,
                    par_filepath,
                    output_path,
                    min_intensity_ratio=0.01,
                    pdf_width=3.35,
                    pdf_height=pdf_width/2.4,
                    xmai = 0.15*pdf_width/3.35,
                    ymai = 0.3*pdf_width/3.35,
                    ppm=20,
                    y_ion_col="red",
                    b_ion_col="blue",
                    peaks_col = "grey",
                    ymax = 1.6,
                    info_height = 1.5,
                    peptide_height = 1.3,
                    mod_height = 0.07,
                    len_annoSpace = 0.1,
                    lwd=1*pdf_width/3.35,
                    cex=1*pdf_width/3.35,
                    show_letterBY=FALSE){
    srt <- 0
    #browser()
    if(! file.exists(output_path)) {
        output_path = paste(output_path, "/", sep="")
        if(! file.exists(output_path)){
            stop(paste("The output dictionary [", output_path, "] does NOT \
                exist!"))
        }
    }
    
    if(! file.exists(par_filepath)){
        stop(paste0("The file '", par_filepath, "' does NOT exist."))
    }
    
    par_files<-data.table::fread(par_filepath, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
    #browser()
    par_ppm <- data.table::rbindlist(apply(par_files, 1, readpar_ppm))
    #browser()
    
    input_table <- data.table::fread(id_table_path, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
    
    # check the columns in id_table_path
    msms_file = base::basename(id_table_path)
    if(! any(colnames(input_table) %in% "Raw file")){
        stop(paste0("The column 'Raw file' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "label")){
        stop(paste0("The column 'label' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Scan number")){
        stop(paste0("The column 'Scan number' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Sequence")){
        stop(paste0("The column 'Sequence' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Modifications")){
        stop(paste0("The column 'Modifications' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Modified sequence")){
        stop(paste0("The column 'Modified sequence' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Gene Names")){
        stop(paste0("The column 'Gene Names' does not exist in ", msms_file, "!"))
    }else if(! any(colnames(input_table) %in% "Charge")){
        stop(paste0("The column 'Charge' does not exist in ", msms_file, "!"))
    }
    input_table$base_rawFile <- basename(input_table$`Raw file`)
    #browser()
    input_table <- check_input_table(input_table, id_table_path, par_ppm, par_filepath)
    #browser()
    lapply(unique(input_table$`Raw file`), drawms2plot_samerawfile, input_table,
        mod_xml_path, output_path, par_ppm, min_intensity_ratio, pdf_width,
        pdf_height, xmai, ymai, y_ion_col, b_ion_col, peaks_col, ymax,
        peptide_height, info_height, mod_height, len_annoSpace, lwd, cex,
        show_letterBY, srt) # call for individual raw_files
    invisible(gc())
}


# load("data/aa_mw_table.rda")
# load("data/atom_mw_table.rda")
# ##save(aa_mw_table, atom_mw_table, PPM_denominator, file = "data/data.rda")
# ##load aa_mw and atom_mw files
# ##aa_mw_table <-   data.table::fread("inst/extdata/AA_MW.txt", sep = "\t",
# ##                                   fill = TRUE, header = TRUE)
# ##atom_mw_table <- data.table::fread("inst/extdata/atom_MW.txt", sep = "\t",
# ##                                   fill = TRUE, header = TRUE)
# 
# 
# mod_xml_path = "inst/extdata/modifications.xml"
# par_filepath = "inst/extdata/par_batch_test.txt"
# id_table_path = "inst/extdata/TMT/msms_TMT_test.txt"
# id_table_path = "inst/extdata/Dimethyl_Labelling/msms_dim_test.txt"
# id_table_path = "inst/extdata/silac/msms_SILAC_test.txt"
# #id_table_path = "inst/extdata/label_free/msms_labelfree_test.txt"
# #
# output_path = "d:"
# mms2plot(id_table_path=id_table_path, mod_xml_path=mod_xml_path,
#        par_filepath=par_filepath, output_path="d:", pdf_width=7, show_letterBY=T)
# 
#library(BiocCheck)
#BiocCheck::BiocCheck("e:/r_packages/mms2plot")
#BiocCheck::BiocCheck("E:/lei_packages3/MMS2plot")

# id_table_path = "D:/test/input/identification.txt"
# mod_xml_path = "D:/test/input/modifications.xml"
# mqpar_filepath = "D:/test/input/parameter_batch.txt"
# output_path = "D:/test"
# mms2plot(id_table_path,mod_xml_path,mqpar_filepath,output_path)
