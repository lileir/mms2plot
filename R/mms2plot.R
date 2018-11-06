#' @author Lei Li
#' @title mms2plot
#' @description Visualization of multiple MS/MSs for (un)modified peptides
#' @export mms2plot
#' @param id_table_path File path name of a table that contains MS2 information
#'        of identified (un)modified peptides plus a group labelling. The format
#'        of MS2 information is referred to as the output file ms2.txt from
#'        the Maxquant search software.
#' @param par_xml_path Xml file path of parameters for modifications and
#'        labelling. The file is the Maxquant parameter file modifications.xml.
#' @param mqpar_filepath File path name that includes a list of parameter files
#'        for search engines. The parameter file format is referred to as
#'        mqpar.xml in Maxquant.
#' @param output_path path of output files.
#' @param min_intensity_ratio minimal percentage threshold of MS2 intensity,
#'        compared with the highest intensity. (default=0.01).
#' @param pdf_width The width of a single PSM figure area in inches. The area
#'        includes both plot region and outer margin area. The default is 3.35
#'        for a single column. The width is 7 for double-column.
#' @param pdf_height The height of a single PSM figure area in inches. The
#'        default is pdf_width/2.4.
#' @param xmai margin of the figure in number of inches for x axis.
#'        (default=pdf_width*0.15/3.35).
#' @param ymai margin of the figure in number of inches for y axis.
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
#' \dontrun{
#' # Generate mms2plot for TMT labelling
#' # id_table_path expands the Maxqaunt output msms.txt by adding "label" column
#' id_table_path = "inst/extdata/TMT/msms_TMT.txt"
#' par_xml_path = "inst/extdata/modifications.xml" # Maxquant modifications.xml
#' # mqpar_filepath contains mqpar.xml with full file path and PPM cutoff
#' mqpar_filepath = "inst/extdata/mqpar_batch.txt"
#' output_path = "inst/extdata"
#' mms2plot(id_table_path, par_xml_path, mqpar_filepath, output_path)
#'
#' # Generate mms2plot for SILAC labelling
#' id_table_path = "inst/extdata/silac/msms_SILAC.txt"
#' par_xml_path = "inst/extdata/modifications.xml"
#' mqpar_filepath = "inst/extdata/mqpar_batch.txt"
#' output_path = "inst/extdata"
#' mms2plot(id_table_path,par_xml_path,mqpar_filepath,output_path,pdf_width=7)
#'
#' # Generate mms2plot for dimethyl labelling
#' id_table_path = "inst/extdata/Dimethyl_Labelling/msms_dim.txt"
#' par_xml_path = "inst/extdata/modifications.xml"
#' mqpar_filepath = "inst/extdata/mqpar_batch.txt"
#' output_path = "inst/extdata"
#' mms2plot( id_table_path, par_xml_path, mqpar_filepath, output_path )
#' }
#rm(list=ls())
#.libPaths( c( .libPaths(), "D:/Rpackages_tmp") )
#library(xml2)
#library(MSnbase)
#library(data.table)
#library(DescTools)  # MixColor
#library(gsubfn)

mms2plot <-function(id_table_path,
                    par_xml_path,
                    mqpar_filepath,
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
    mqpar_files<-data.table::fread(mqpar_filepath, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)

    mqpar_ppm <- data.table::rbindlist(apply(mqpar_files, 1, readMQPar_ppm))

    input_table <- data.table::fread(id_table_path, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)

    input_table$base_rawFile <- basename(input_table$`Raw file`)
    #browser()
    input_table <- check_input_table(input_table, id_table_path, mqpar_ppm)
    #browser()
    lapply(unique(input_table$`Raw file`), drawms2plot_samerawfile, input_table,
        par_xml_path, output_path, mqpar_ppm, min_intensity_ratio, pdf_width,
        pdf_height, xmai, ymai, y_ion_col, b_ion_col, peaks_col, ymax,
        peptide_height, info_height, mod_height, len_annoSpace, lwd, cex,
        show_letterBY, srt) # call for individual raw_files
    invisible(gc())
}

#rm(list=ls())
#load("data/data.rda")
#save(aa_mw_table, atom_mw_table, PPM_denominator, file = "data/data.rda")
# load aa_mw and atom_mw files
# aa_mw_table <-   data.table::fread("inst/extdata/AA_MW.txt", sep = "\t",
#                                    fill = TRUE, header = TRUE)
# atom_mw_table <- data.table::fread("inst/extdata/atom_MW.txt", sep = "\t",
#                                    fill = TRUE, header = TRUE)
# PPM_denominator=1E6
#
#source("R/plot_mirror_or_group.R")
#source("R/add_mod_aa.R")
#source("R/psm_calculation.R")
#source("R/plot_components.R")

#mqpar_filepath = "inst/extdata/mqpar_batch.txt"
#par_xml_path = "inst/extdata/modifications.xml"
#id_table_path = "inst/extdata/TMT/msms_TMT.txt"
#id_table_path = "inst/extdata/Dimethyl_Labelling/msms_dim.txt"
#id_table_path = "inst/extdata/silac/msms_SILAC.txt"

#mms2plot(id_table_path=id_table_path, par_xml_path=par_xml_path,
#          mqpar_filepath=mqpar_filepath, output_path="d:", pdf_width=7)

#library(BiocCheck)
#BiocCheck::BiocCheck("e:/r_packages/mms2plot")
#BiocCheck::BiocCheck("e:/test/mms2plot")

