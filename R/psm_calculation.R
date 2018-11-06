#attach mod to AA (including before or after AA) and calculate theoretical m/z
#of b ions and y ions with charge 1 or 2
calculate_aa_mzs <-function(seq, charge, Monoisotopicmz, ppm, aa_mw_mod_table){

    #browser()
    AA_mzs <- list()
    if(length(unique(aa_mw_mod_table$labelmod_group))>1){ # aa with label mod
        AA_mzs<-by(aa_mw_mod_table, aa_mw_mod_table$labelmod_group,
            calculate_AAmz_individual_label, seq, charge, Monoisotopicmz, ppm,
            flag="labelmod", simplify=FALSE)
    }else if(length(unique(aa_mw_mod_table$reporterion_group))>1 ||
        grepl("plex|TMT", unique(aa_mw_mod_table$reporterion_group))){
        #browser()
        AA_mzs<-by(aa_mw_mod_table, aa_mw_mod_table$reporterion_group,
            calculate_AAmz_individual_label, seq, charge, Monoisotopicmz, ppm,
            flag="reporterion", simplify = FALSE)
    }else{ # aa with label free
        #browser()
        AA_mzs<-calculate_AAmz_individual_label(aa_mw_mod_table, seq, charge,
            Monoisotopicmz, ppm )
        AA_mzs<-list(AA_mzs)
    }
    #browser()
    # remove empty elements in the list AA_mz
    AA_mzs_final <- AA_mzs[lengths(AA_mzs) != 0]

    if(length(AA_mzs_final) == 1 ){
        return(AA_mzs_final[[1]])
    }else{
        # browser()
        stop("more than 1 labelling matches. Please contact the developer! \
            [note:stopped in calculate_aa_mzs()].")
    }
}

# for label-free, labelled or reportor ions (e.g. TMT)
calculate_AAmz_individual_label <-function(aa_mw_mod_table, seq, charge,
    Monoisotopicmz, ppm, flag){

    H_weight <- subset(mms2plot::atom_mw_table,
        mms2plot::atom_mw_table$Element == "H")$Monoisotopic
    H2O_weight <- 2* H_weight + subset(mms2plot::atom_mw_table,
        mms2plot::atom_mw_table$Element == "O")$Monoisotopic

    #seq<-"AAAVLPVLDLAQR"
    #charge <- 3
    ms1_mzThreshold <- 0.5 # therotical ms1 m/z minus measured ms1 m/z
    PPM_denominator <- 1E6
    mz_range <- ppm/PPM_denominator
    #browser()
    # using calculateFragments for testing ##########################
    #if(charge==2){
    #  a=calculateFragments(seq)
    #}else{
    #  a=calculateFragments(seq)
    #  b=calculateFragments(seq, z=2)
    #}
    ##########################################

    seq_origin <- seq
    mod_N <- "" # mod at prot N terminal, set as empty string, nchar(mod_N) == 0
    if(grepl("^\\(", seq)){
        if( substr(seq, 6,6) == "(" ) { # the first AA is modified at both sides
            mod_N<-substr(seq[1], 1, 9)
            seq<-substring(seq, 10)
        }else{  # the first AA is modified at left side
            mod_N<-substr(seq[1], 1, 5)
            seq<-substring(seq, 6)
        }
    }
    seq_sep <- strsplit(unlist(strsplit(seq, "\\(\\w\\w\\)")), "")
    mod <- unlist(strsplit(seq, "[A-Z]+"))
    tmp <- paste(mod,collapse ="")
    AA <- vector()
    if(nchar(tmp) != 0 ){ # contain mod
        mod <- c(mod[2:length(mod)],"")

        # add mod after AA
        add_mod_C <-function(subseq,mod){
            subseq[length(subseq)] <- paste(subseq[length(subseq)], mod, sep="")
            return(subseq)
        }
        AA<-unlist(mapply(add_mod_C, seq_sep, mod))

        if(nchar(mod_N) > 0 ){AA <- c(mod_N, AA)}
    }else{
        AA <- unlist(seq_sep)
        if(nchar(mod_N)>0){ AA <- c(mod_N, AA) }
    }
    #calculate theoretical m/z of b ions and y ions
    AA_seq <- data.table::data.table(aa_varmod=AA, index=seq_len(length(AA)))

    AA_mz <- data.table::setDT(aa_mw_mod_table)[AA_seq, on="aa_varmod"]
    AA_mz_NA<-subset(AA_mz, is.na(AA_mz$weight)) #extract rows with NA weight
    if(nrow(AA_mz_NA)>0){
        mod_NA_weight <- paste(as.character(AA_mz_NA$AA), collapse=" & ")
        stop(paste("The modification(s) [ ", mod_NA_weight, " ] in '",
            seq_origin, "' is not included in the modification.xml file! \
            [note:stopped in calculate_AAmz_individual_label()].", sep=""))
    }

    # for reportor ion (e.g. TMT)
    if(flag == "reporterion"){
        AA_mz_final_start <-
            subset(AA_mz, AA_mz$reporterion %in% "anyNterm" & AA_mz$index == 1)
        if(nrow(AA_mz_final_start)>0){
            AA_mz_final_others <-
                subset(AA_mz, is.na(AA_mz$reporterion)      & AA_mz$index != 1)
            AA_mz <- rbind(AA_mz_final_start, AA_mz_final_others)
        }
    }else if(flag == "labelmod"){
        AA_mz_final_start <-
            subset(AA_mz, AA_mz$labelmod %in% "anyNterm" & AA_mz$index == 1)
        if(nrow(AA_mz_final_start)>0){
            AA_mz_final_others <-
                subset(AA_mz, is.na(AA_mz$labelmod)      & AA_mz$index != 1)
            AA_mz <- rbind(AA_mz_final_start, AA_mz_final_others)
        }
    }
    AA_mz_b <- AA_mz[with(AA_mz, order(AA_mz$index)),]
    AA_mz_b$mz_b <- cumsum(AA_mz_b$weight)
    AA_mz_b$mz_b <- AA_mz_b$mz_b + H_weight

    # sort and calculate MW for y ion for
    AA_mz_y <- AA_mz[with(AA_mz, order(-AA_mz$index)),]
    AA_mz_y$mz_y <- cumsum(AA_mz_y$weight)
    AA_mz_y$mz_y <- AA_mz_y$mz_y + H2O_weight + H_weight

    AA_mz <- merge(AA_mz_b, AA_mz_y)
    #browser()
    # add N-terminal and C-terminal weight minus one original H_weight atom
    mz_ideal <- sum(sum( AA_mz$weight) + H2O_weight + H_weight*charge )/charge
    # If ideal MS1 Monoisotopic mz is not similar to measured Monoisotopic mz,
    # this type of labelling does not fit.
    if(abs(mz_ideal-Monoisotopicmz) > ms1_mzThreshold){
        return(NULL)
    } else{
        # sort and calculate MW for b ion
        AA_mz_b <- AA_mz[with(AA_mz, order(AA_mz$index)),]
        AA_mz_b$mz_b <- cumsum(AA_mz_b$weight)
        AA_mz_b$mz_b <- AA_mz_b$mz_b + H_weight

        # sort and calculate MW for y ion for
        AA_mz_y <- AA_mz[with(AA_mz, order(-AA_mz$index)),]
        AA_mz_y$mz_y <- cumsum(AA_mz_y$weight)
        AA_mz_y$mz_y <- AA_mz_y$mz_y + H2O_weight + H_weight

        AA_mz <- merge(AA_mz_b, AA_mz_y)
        AA_mz$charge <- 1

        if(charge >2){ # only calculate the b/y ions with charge of 2
            charge <- 2
            AA_mz2 <- AA_mz
            # adjusted, as AA_mz2_b$mz_b already add H_weight for charge 1
            AA_mz2$mz_b <- (AA_mz$mz_b + H_weight)/charge
            # adjusted, as AA_mz2_y$mz_y already add H2O_weight+H_weight for 1+
            AA_mz2$mz_y <- (AA_mz$mz_y + H_weight)/charge
            AA_mz2$charge <- 2
            AA_mz <- rbind(AA_mz, AA_mz2)
        }

        AA_mz$mz_b_min <- (1-mz_range)*AA_mz$mz_b
        AA_mz$mz_b_max <- (1+mz_range)*AA_mz$mz_b
        AA_mz$mz_y_min <- (1-mz_range)*AA_mz$mz_y
        AA_mz$mz_y_max <- (1+mz_range)*AA_mz$mz_y

        AA_mz <- AA_mz[order(AA_mz$charge,AA_mz$index)]
        #browser()
        return(AA_mz)
    }
}


# Find matched b/y ions from MS2 peaks
test_individualIon<-function(mz_intensity_percent, AA_mz, b_ion_col, y_ion_col){
    #browser()
    # intensity
    MIP <- data.frame(as.list(mz_intensity_percent))
    ion_info<-data.table::data.table()
    # search for b ion match
    ion <- subset(AA_mz, MIP$mz > AA_mz$mz_b_min & MIP$mz < AA_mz$mz_b_max)
    # mz matched and mz is smaller than the mw of the AA sequence (using index)
    if(nrow(ion) == 1 && ion$index < max(AA_mz$index)){
        ion_info <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
            "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
            "abs_intensity_prc_ext" = ifelse(abs(MIP$intensity_perc)>0.30,
            abs(MIP$intensity_perc),
            stats::runif(sum(abs(MIP$intensity_perc)<0.30),min=0.30, max=0.65)),
            "ionLabel" = paste("b", ion$index, paste(rep("+", ion$charge),
                collapse=""), sep=""),
            "ion" = paste("b", ion$index, sep=""), "col" = b_ion_col,
            "direction" = ifelse(MIP$intensity_perc>0, 1, -1))
    }
    ion <- NULL

    # search for y ion match
    ion <- subset(AA_mz, MIP$mz > AA_mz$mz_y_min & MIP$mz < AA_mz$mz_y_max)
    # mz matched and mz is smaller than the mw of the AA sequence (using index)
    # index is for b index. calculate y ion by minux b index
    if(nrow(ion) == 1 && max(AA_mz$index)-ion$index+1 < max(AA_mz$index) ){
        yion <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
            "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
            "abs_intensity_prc_ext" = ifelse(abs(MIP$intensity_perc)>0.7,
            abs(MIP$intensity_perc),
            stats::runif(sum(abs(MIP$intensity_perc)<0.7), min=0.7, max=1)),
            "ionLabel" = paste("y", max(AA_mz$index)-ion$index+1,
                paste(rep("+", ion$charge),collapse=""), sep=""),
            "ion" = paste("y", max(AA_mz$index)-ion$index+1, sep=""),
            "col" = y_ion_col,"direction" = ifelse(MIP$intensity_perc>0, 1, -1))
        if(nrow(ion_info)>0){
            ion_info <- rbind(ion_info, yion)
        }else{
            ion_info <- yion
        }
    }
    return(ion_info)
}

# find matched ions from MS2
find_matchedIons<-function(AA_mz, mz_intensity_percent, b_ion_col, y_ion_col){
    #browser()
    # find b/y ions with info
    psm <- apply( mz_intensity_percent, 1, test_individualIon, AA_mz, b_ion_col,
        y_ion_col)
    psm <- data.table::rbindlist(psm) # can remove NULL elements
    if(nrow(psm)<4){stop("The matched ion peaks are limited (<4). Please check \
        the ppm threshold. [note:stopped in the function find_matchedIons].")}

    #browser()
    # Keep the mz with the largest intensity if multiple mzs match same b/y ion
    # find the max per row and retain the columns with the max
    abs_intensity_prc_ext_max <-
        stats::aggregate(intensity~ionLabel, data=psm, max)
    psm<-merge(abs_intensity_prc_ext_max, psm)

    ####### test #################
    #  y_psm[1,2] <- b_psm[1,2]
    #  y_psm[1,3] <- b_psm[1,3]
    #  y_psm[1,4] <- b_psm[1,4]
    ################################

    #  yb_psm <- rbind(b_psm, y_psm)
    #  yb_psm$ion <- unlist(strsplit(yb_psm$ionLabel, split="\\+|\\("))

    # If one m/z corresponds to multi b/y ion, these ions are labelled together
    # original ionLabel: label for b/y ion annotation;
    # ionlabel_peak: labelling of peaks for b/y ions
    ionLabel_comb <-
        stats::aggregate(ionLabel~abs_intensity_prc_ext+mz+intensity_perc,
            psm, paste, collapse=",")
    psm<-merge(ionLabel_comb, psm,
        by=c("abs_intensity_prc_ext", "mz","intensity_perc"),
        suffixes = c("_peak", ""))
    if(any(grepl(",", psm$ionLabel_peak))){
        psm[grepl(",", psm$ionLabel_peak),]$col =
            DescTools::MixColor(y_ion_col, b_ion_col)
    }
    #browser()
    return(psm)
}

# add intensity_perc column as values/max &
# (intensity_perc * -1 for downMS2 in mirrorplot)
get_intensity_perc<-function(input_table, min_intensity_ratio){

    # funtion to retrieve mz, intensity and intensity_perc
    #for each row of input_table
    get_peakswoNoise<-function(input_table, min_intensity_ratio){
        max_intensity <- as.numeric(input_table["max_intensity"])
        peaks<- data.table::data.table(
            mz = as.numeric(unlist(strsplit(input_table["mz"],";"))),
            intensity=as.numeric(unlist(strsplit(input_table["intensity"],";")))
            )
        # remove intensity < 1% of max_intensity
        peaks_wo_Noise<-
            subset(peaks, peaks$intensity > max_intensity * min_intensity_ratio)
        peaks_wo_Noise$intensity_perc <- peaks_wo_Noise$intensity/max_intensity
        return(peaks_wo_Noise)
    }
    peaks <- apply(input_table, 1, get_peakswoNoise,  min_intensity_ratio)
    if(length(peaks) == 2) {  # for plot_mirror
        peaks[[2]]$intensity_perc <- peaks[[2]]$intensity_perc * -1
    }
    return(peaks)
}
