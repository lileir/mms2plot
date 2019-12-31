# readpar_ppm function for individual par.xml files and ppms
# output:
#      rawfile  variableModifications             isobaricLabels       labelMods
#   test.mzXML                Oxid_11                         NA              NA
#  test2.mzxml          Oxidation (M)           TMT6plex-Nter127              NA
#
# comma(,) separae filenames, different types of labelling or reporter ion
# semicolon(;) used in maxquant to separate labelled aa in the same label types
readpar_ppm <- function(par_filename) {
    #browser()
#    ion_type <-par_filename["ion_type"]  
    ion_type = switch(par_filename["ion_type"],  # y or z or yz for b/y or c/z or b/y/c/z
       "y" = par_filename["ion_type"],
       "z" = par_filename["ion_type"],
       "yz" = par_filename["ion_type"],
       stop(paste0("'", par_filename["ion_type"], "' is not included as the aceptable labels (y, z or yz) in the 'ion_type' column for '", par_filename["par_path"], "'!"))
    )
    
    ppm <- par_filename["ppm"]
    xmlread <- xml2::read_xml(par_filename["par_path"])
    filePaths <-  xml2::xml_find_all(xmlread, "//filePaths") # for files
    variableModifications <-
        xml2::xml_find_all(xmlread, "//variableModifications")
    isobaricLabels <-
        xml2::xml_find_all(xmlread, "//isobaricLabels") #reporter ion: TMT/iTRAQ
    labelMods <- xml2::xml_find_all(xmlread, "//labelMods") #labelling: Arg10
    fixedModifications <- xml2::xml_find_all(xmlread, "//fixedModifications")

    filePaths_to_generalpaths = gsub("\\", "/", xml2::xml_text(xml2::xml_children(filePaths)), fixed = T)
    filePaths <- basename(filePaths_to_generalpaths)
    variableModifications <-
        xml2::xml_text(xml2::xml_children(variableModifications))
    isobaricLabels <- xml2::xml_text(xml2::xml_children(isobaricLabels))
    labelMods <- xml2::xml_text(xml2::xml_children(labelMods))
    fixedModifications <- xml2::xml_text(xml2::xml_children(fixedModifications))

    if (all(nchar(filePaths) > 0)) {
        # connect multiple input_files by comma
        rawfile <- paste(filePaths, collapse = ",")

        if (all(nchar(variableModifications) == 0) ||
            length(variableModifications) == 0) {
            variableModifications <- NA
        } else{
            variableModifications <- paste(variableModifications,collapse = ",")
        }

        if (all(nchar(isobaricLabels) == 0)  || length(isobaricLabels) == 0) {
            isobaricLabels <- NA
        } else{
            isobaricLabels <- paste(isobaricLabels, collapse = ",")
        }

        if (all(nchar(labelMods) == 0) || length(labelMods) == 0) {
            labelMods <- NA
        } else{
            labelMods <- paste(labelMods, collapse = ",")
        }

        if(all(nchar(fixedModifications) == 0)||length(fixedModifications)==0) {
            #Arg10;Lys8
            fixedModifications <- NA
        } else{
            fixedModifications <- paste(fixedModifications, collapse = ",")
        }

        output <- data.table::data.table(rawfile, variableModifications,
            isobaricLabels, labelMods, fixedModifications, ppm, ion_type)
#        browser()
        return(output)
    } else{
        stop("No raw file is included in this par.xml!")
    }
    #browser()
    return(NULL)
}



# check if the input_table has the expected format
check_input_table<-function(input_table, id_table_path, par_ppm, par_filepath ){
    #browser()
    
    if(nrow(input_table) == 0){ stop(paste("The file", id_table_path,
        "is empty! Please see the example file. \
        [note:stopped in check_input_table()].")) }
    #print(par_ppm)
    #browser()
    unique_rawfiles <-
        base::unique(tools::file_path_sans_ext(basename(input_table$`Raw file`)))
    parfile <-
        base::unique(tools::file_path_sans_ext(unlist(strsplit(par_ppm$rawfile, ","))))
    
    if(! all(unique_rawfiles %in% parfile)){
        rawfiles=base::setdiff(unique_rawfiles, parfile)
        rawfiles_collapse = base::paste(rawfiles, collapse = "/")
        stop(paste0("The raw file(s) '",rawfiles_collapse, "' in ", id_table_path, " are not included in par file(s) of ",par_filepath))
    }

    col_check <- c("Raw file","Scan number","Sequence","Modifications",
        "Gene Names","label")

    # check if specific columns exist in the input_table
    col_not_exist <- col_check[!col_check %in% colnames(input_table)]
    if(length(col_not_exist) > 0 ){
        stop(paste("The column(s) '", paste(col_not_exist,collapse = " and ") ,
            "' is(are) not existent! [note:stopped in check_input_table()].",
            sep = "" ))
    }

    # check if each label group only contains a single raw file names
    #browser()
    unique_rawFile <-
        input_table[, unique(input_table$`Raw file`), by = input_table$label]
    #if(nrow(unique_rawFile) != nrow(unique(unique_rawFile[,1]))) {
        #stop("MS2 IDs in each group should be derived from the same MS file! \
            #[note:stopped in check_input_table()].")}

    input_table$`Modified sequence` <-
        gsub("_","", input_table$`Modified sequence`)
    input_table_unmod <-
        subset(input_table, nchar(input_table$`Modified sequence`) == 0)
    input_table_unmod$`Modified sequence` = input_table_unmod$Sequence
    input_table_mod <-
        subset(input_table, nchar(input_table$`Modified sequence`)>0)
    # remove modificaiton to check sequence identity
    mod_remove_brackets <-
        gsub("\\(\\w\\w\\)", "", input_table_mod$`Modified sequence`)
    if( !all(input_table_mod$Sequence == mod_remove_brackets) ){
        stop("Sequences in the 'Sequence' column are different from those in \
            the 'Modified sequence' column ! \
            [note:stopped in check_input_table()].")
    }

    # modify modseq column by adding related unmodified peptides in empty rows
    input_table$`Modified sequence` <-
        mapply(function(x,y){ifelse((is.na(y) || nchar(y) == 0), x, y)},
            input_table$Sequence, input_table$`Modified sequence`)
    return(input_table)
}


# drawms2plot_samerawfile function
# get unique MS2 number & read info for each raw file using get_mzIntensity
# call plot_mms2 for each label type
# aa_mw_mod_table <- list_aaMwModTable_ppm[[1]]
# ppm <- list_aaMwModTable_ppm[[2]]
drawms2plot_samerawfile <- function(MS2FileName, input_table,  mod_xml_path,
    output_path, par_ppm, min_intensity_ratio, pdf_width, pdf_height,
    xmai, ymai, y_ion_col, b_ion_col, peaks_col, ymax, peptide_height,
    info_height, mod_height, len_annoSpace, lwd, cex, show_iontype=TRUE, srt){

    # extract from input_table MS2 info from the same file
    input_table_sameRawFile <-
        input_table[input_table$base_rawFile == basename(MS2FileName)]
    # Processing modification.xml files
    # And read site, title, composition and merge into aa_mw_table.
    # add mod_aa to the table, labelling data are annotated by group flag
    #browser()
    list_aaMwModTable_ppm<-add_mod_aa(mod_xml_path, basename(MS2FileName),
        mms2plot::aa_mw_table, par_ppm)

    aa_mw_mod_table = list_aaMwModTable_ppm[[1]] ## aa plus modaa
    ppm = list_aaMwModTable_ppm[[2]]  ##  ppm
    ion_type = list_aaMwModTable_ppm[[3]] # y, z or yz ion
    #browser()

    # unique MS2 scan_number from the extract MS2 info
    scan_number <- unique(input_table_sameRawFile$`Scan number`)

    print(paste("Reading the raw MS file: ", MS2FileName, "... ..."))
    MS2s_frFile <-
        MSnbase::readMSData(MS2FileName, msLevel=2, mode="onDisk",verbose=FALSE)
    print("Reading ... ... completed!")

    # call get_ms2info function to extract MS2 M/Z, intensities and etc.
    mzIntensity_list <-  lapply(scan_number, get_ms2info, MS2s_frFile)
    # Garbage Collection for MS2s_frFile
    rm(MS2s_frFile);  invisible(gc())
    #browser()
    # change list as data.frame, each row contain one MS2 info
    mzIntensity <- do.call(rbind, mzIntensity_list)
    scannumber_charge = base::unique(base::subset(input_table_sameRawFile, select=c("Scan number", "Charge")))
    ### charge retrieved from MS raw file is replaced by the charge determined by users
    mzIntensity_tmp  = base::merge(mzIntensity, scannumber_charge, by="Scan number")
    if(nrow(mzIntensity_tmp) == nrow(mzIntensity) ){
        mzIntensity = mzIntensity_tmp
    }else{
        stop("There is an error in adding charge from user input file! Please contact the developer!")
    }
    # merge input_table and mzIntensity
    input_table_sameRawFile <-
        data.table::setDT(mzIntensity)[input_table_sameRawFile,on="Scan number"]
    #browser()
    tmp<-by(input_table_sameRawFile, input_table_sameRawFile$label, plot_mms2,
        output_path, aa_mw_mod_table, min_intensity_ratio, pdf_width,
        pdf_height, xmai, ymai, ppm, ion_type, y_ion_col,
        b_ion_col, peaks_col, ymax, peptide_height, info_height, mod_height,
        len_annoSpace, lwd, cex, show_iontype, srt)
    invisible(gc())
}

# Get MS2 m/z and intensities
get_ms2info <- function(scan_number, ms2_samefile){
    #browser()
    MS1table <- Biobase::fData(ms2_samefile)
    MS1_specific <-
        subset(MS1table, MS1table$acquisitionNum == scan_number)$filterString
    if( length(MS1_specific) ==0){
        stop( paste("The scan_number [", scan_number, "] is not found in the \
            raw MS file! [note:stopped in the function get_ms2info].", sep="") )
    }
    # mzML: get ms1 mz; NA for mzXML
    MS1_mz <- as.numeric(utils::tail(unlist(strsplit(unlist(
        strsplit(MS1_specific, "@"))[1] , " ")),n=1))

    ms2_info <- ms2_samefile[[which(Biobase::fData(ms2_samefile)$acquisitionNum
        == scan_number)]]  # fData(ms2_samefile) shows all the MS2 info
    ms2_info_table <- data.table::data.table("Scan number"=scan_number,
        max_intensity = max( ms2_info@intensity ), `Retention time`=ms2_info@rt,
        `m/z` = ifelse(is.na(MS1_mz), ms2_info@precursorMz, MS1_mz),
        #Charge = ms2_info@precursorCharge, # replace the charge from the raw file by the charge determined by the user 
        Monoisotopicmz = ms2_info@precursorMz,
        mz=paste(round( ms2_info@mz, digits = 3), collapse=";" ),
        intensity = paste( round( ms2_info@intensity, digits=3 ), collapse=";"))
    #browser()
    return(ms2_info_table)
}


# plot_mms2 function
plot_mms2 <- function(input_table, output_path, aa_mw_mod_table,
    min_intensity_ratio, pdf_width, pdf_height, xmai, ymai, ppm, ion_type, y_ion_col,
    b_ion_col, peaks_col, ymax, peptide_height, info_height, mod_height,
    len_annoSpace, lwd, cex, show_iontype, srt){
    #browser()
    if(nrow(input_table) == 2){ # for mirror plot
        plot_mirror(input_table, output_path, aa_mw_mod_table,
            min_intensity_ratio, pdf_width, pdf_height, xmai, ymai, ppm, ion_type,
            y_ion_col, b_ion_col, peaks_col, ymax, peptide_height, info_height,
            mod_height, len_annoSpace, lwd, cex, show_iontype, srt)
    }else{ # for aligned plot
        plot_group(input_table, output_path, aa_mw_mod_table,
            min_intensity_ratio, pdf_width, pdf_height, xmai, ymai, ppm, ion_type, 
            y_ion_col, b_ion_col, peaks_col, ymax, peptide_height, info_height,
            mod_height, len_annoSpace, lwd, cex, show_iontype, srt)
    }
    invisible(gc())
}


# the main function to draw mirror_plot
plot_mirror <- function(input_table, output_path, aa_mw_mod_table,
    min_intensity_ratio, pdf_width, pdf_height, xmai, ymai, ppm, ion_type, y_ion_col,
    b_ion_col, peaks_col, ymax, peptide_height, info_height, mod_height,
    len_annoSpace, lwd, cex, show_iontype, srt){
    
    #browser()
    outputFilename <- paste(input_table$base_rawFile[1], input_table$label[1],
        input_table$Sequence[1], "mirror", sep = "_")
    outputFilename <- paste(output_path, outputFilename, sep="/")
    outputFilename <- paste(outputFilename,"pdf",sep=".")

    # set intensity value range between 0 to 1
    mz_intensity_percent <- get_intensity_perc(input_table, min_intensity_ratio)
    two_mz_intensity_percent <- do.call(rbind, mz_intensity_percent)
    # calculate therotical b/y ions for the given peptides and ppm is considered
    #browser()
    AA_mzs <- mapply(calculate_aa_mzs,input_table$`Modified sequence`,
        input_table$Charge, input_table$Monoisotopicmz,
        MoreArgs=list(ppm, ion_type, aa_mw_mod_table), SIMPLIFY = FALSE )
    #browser()
    # calcualte PSM for each MS2 plot iteratively
    PSMs <- mapply(find_matchedIons, AA_mzs, mz_intensity_percent,
        MoreArgs=list(ion_type, b_ion_col, y_ion_col), SIMPLIFY = FALSE )
    two_PSMs <- do.call(rbind, PSMs)
    #browser()
    grDevices::graphics.off()
    pdf(file=outputFilename, width=pdf_width, height=pdf_height*2)

    old_options <- options(scipen=22) # max of fixed notation: 1e22
    on.exit(options(old_options), add = TRUE)
    old_par <- par()
    on.exit(options(old_options), add = TRUE)

    options(scipen=22)

    # Set the outer margin area. Note that this parameter applies for the entire
    # figure and altering it while drawing figures may cause problems.
    par(oma=c(1,1,1,0), mai=c(xmai,ymai,0,ymai), cex=cex, lwd=lwd)
    max_mz<-max(unlist(lapply(mz_intensity_percent, `[[`, 1)))+100 # max mz
    #########################################################################

    graphics::plot(two_mz_intensity_percent$mz,
        two_mz_intensity_percent$intensity_perc, type = "h",las = 1,
        xlab = "", ylab = "", xlim = c(0, max_mz), xaxt="n", yaxt="n", xaxs="i",
        yaxs="i", ylim = c(-1*ymax, ymax), col = peaks_col, cex = cex, lwd=0.5,
        frame.plot=FALSE)
    graphics::box(lwd=lwd/2)

    # label x/y axis
    label_axis(input_table$max_intensity, xlab="m/z", ylab="Intensity", max_mz,
        cex, lwd)
    abline(h = 0, lwd=0.5)
    # draw peaks, extension lines and ionLabel
    draw_peak_ionL(two_PSMs, lwd, len_annoSpace, srt, show_iontype)
    ## draw PSM annotation (peptide, mod, b ions and y ions)
    #browser()
    mapply(draw_psmanno, AA_mzs, PSMs, max_mz, peptide_height, mod_height,
        len_annoSpace, y_ion_col, b_ion_col, lwd=lwd/2)
    #browser()
    ## MS2 information ( file name, RT, Scan number, m/z, charge and gene name)
    mapply(draw_ms2generalinfo, input_table$`Retention time`,
        input_table$`Scan number`, input_table$`m/z`, input_table$Charge,
        input_table$`Gene Names`, PSMs, info_height)
    # title outside oma
    graphics::mtext("m/z", side=1, line=0, cex=0.5*cex, outer=TRUE)
    graphics::mtext("Intensity", side=2, line=0, cex=0.5*cex, outer=TRUE)
    rawFileName <- paste("File:", input_table$base_rawFile[1])
    graphics::mtext(rawFileName, side=3, line=0, cex=0.33*cex, outer=TRUE)

    print(paste("The pdf file '", outputFilename, "' was generated.", sep=""))
    dev.off()
    #browser()
}

plot_group <-function(input_table, output_path, aa_mw_mod_table,
    min_intensity_ratio, pdf_width, pdf_height, xmai, ymai, ppm, ion_type, y_ion_col,
    b_ion_col, peaks_col, ymax, peptide_height, info_height, mod_height,
    len_annoSpace, lwd, cex, show_iontype, srt){

    #browser()
    outputFilename <- paste(input_table$base_rawFile[1], input_table$label[1],
        input_table$Sequence[1], "align", sep = "_")
    outputFilename <- paste(output_path, outputFilename, sep="/")
    outputFilename <- paste(outputFilename,"pdf",sep=".")
    # set intensity value range between 0 to 1
    mz_intensity_percent <- get_intensity_perc(input_table, min_intensity_ratio)
    # calculate therotical b/y ions for the given peptides and ppm is considered
    AA_mzs <- mapply(calculate_aa_mzs,input_table$`Modified sequence`,
        input_table$Charge, input_table$Monoisotopicmz,
        MoreArgs=list(ppm, ion_type, aa_mw_mod_table), SIMPLIFY = FALSE )
    #browser()
    # calcualte PSM for each MS2 plot iteratively
    PSMs <- mapply(find_matchedIons, AA_mzs, mz_intensity_percent,
        MoreArgs=list(ion_type, b_ion_col, y_ion_col), SIMPLIFY = FALSE )
    #browser()
    #max_mz<-max(input_table$`m/z` * input_table$Charge)
    max_mz<-max(unlist(lapply(mz_intensity_percent, `[[`, 1)))+100 # max mz
    ##  draw pictures ########################

    grDevices::graphics.off()
    grDevices::pdf(file=outputFilename, width=pdf_width,
        height=pdf_height*nrow(input_table))

    options(scipen=22) # max of fixed (not exponential) notation: 1e22

    # Set the outer margin area. Note that this parameter applies for the entire
    #   figure and altering it while drawing figures may cause problems.
    nplot <- nrow(input_table)
    par(mfrow=c(nplot, 1), oma=c(2,1,1,0),mai=c(0, ymai,0,ymai),cex=cex,lwd=lwd)

    # Curtain layering
    mapply(plot_group_individual, mz_intensity_percent, AA_mzs, PSMs,
        input_table$max_intensity, input_table$base_rawFile,
        input_table$`Retention time`, input_table$`Scan number`,
        input_table$`m/z`, input_table$Charge, input_table$`Gene Names`,
        MoreArgs=list(y_ion_col, b_ion_col, peaks_col, ymax, lwd,
        peptide_height, mod_height, len_annoSpace, srt, cex, max_mz,
        info_height, show_iontype))

    # title outside oma
    graphics::mtext("m/z", side=1, line=1, cex=0.5*cex, outer=TRUE)
    graphics::mtext("Intensity", side=2, line=0, cex=0.5*cex, outer=TRUE)
    rawFileName <- paste("File:", input_table$base_rawFile[1])
    graphics::mtext(rawFileName, side=3, line=0, cex=0.33*cex, outer=TRUE)
    # x label
    axis(side = 1, at=seq(0, max_mz,by=200),lwd=0.5*lwd,tck=-0.02,labels=FALSE)
    graphics::mtext(side = 1, text=seq(0, max_mz, by=200),
        at=seq(0, max_mz, by=200), cex=0.33*cex)

    print(paste("The pdf file '", outputFilename, "' was generated.", sep=""))
    dev.off()
    #browser()
}
