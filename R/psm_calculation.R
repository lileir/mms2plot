#attach mod to AA (including before or after AA) and calculate theoretical m/z
#of b ions and y ions with charge 1 or 2
calculate_aa_mzs <- function(seq, charge, Monoisotopicmz, ppm, ion_type, aa_mw_mod_table){

    #browser()
    aa_mw_mod_table$bcyz = "" # add ion_type
    AA_mzs <- list()
    if(length(unique(aa_mw_mod_table$labelmod_group))>1){ # mod or SILAC
        AA_mzs<-by(aa_mw_mod_table, aa_mw_mod_table$labelmod_group,
            calculate_AAmz_individual_label, seq, charge, Monoisotopicmz, ppm, 
            ion_type, flag="labelmod")#, simplify=FALSE)
    }else if(length(unique(aa_mw_mod_table$reporterion_group))>1 ||
            grepl("plex|TMT", unique(aa_mw_mod_table$reporterion_group))){
        #browser()
        AA_mzs<-by(aa_mw_mod_table, aa_mw_mod_table$reporterion_group,
            calculate_AAmz_individual_label, seq, charge, Monoisotopicmz, ppm, 
            ion_type, flag="reporterion")#, simplify = FALSE)
    }else{ # aa with label free
        #browser()
        AA_mzs<-calculate_AAmz_individual_label(aa_mw_mod_table, seq, charge,
            Monoisotopicmz, ppm, ion_type, flag="")#, simplify = FALSE )
        #browser()
        AA_mzs<-list(AA_mzs)
    }
    #browser()
    # remove empty elements in the list AA_mz
    AA_mzs_final <- AA_mzs[lengths(AA_mzs) != 0]
    return(AA_mzs_final)
}

# AAmz for label-free, labelled or reportor ions (e.g. TMT) for b/y ions
calculate_AAmz_individual_label <-function(aa_mw_mod_table, seq, charge,
    Monoisotopicmz, ppm, ion_type, flag){

    H_weight <- subset(mms2plot::atom_mw_table,
        mms2plot::atom_mw_table$Element == "H")$Monoisotopic
    H2O_weight <- 2 * H_weight + subset(mms2plot::atom_mw_table,
        mms2plot::atom_mw_table$Element == "O")$Monoisotopic
    NH3_weight <- 3 * H_weight + subset(mms2plot::atom_mw_table,
        mms2plot::atom_mw_table$Element == "N")$Monoisotopic
    #seq<-"AAAVLPVLDLAQR"
    #charge <- 3
    ms1_mzThreshold <- 1#0.5 # therotical ms1 m/z minus measured ms1 m/z
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
    #browser()
    AA <- vector()
    if(nchar(tmp) != 0 ){ # contain mod
        mod <- c(mod[2:length(mod)],"")

        # add mod after AA
        add_mod_C <-function(subseq,mod){
            subseq[length(subseq)] <- paste(subseq[length(subseq)], mod, sep="")
            return(subseq)
        }
        #browser()
        if(length(seq_sep) == length(mod)){ # mod in the middle of the seq
            AA<-as.vector(unlist(mapply(add_mod_C, seq_sep, mod))) # AA as a AA vector
        }else if( length(mod) - length(seq_sep) == 1  ){ # mod at end of the seq
            if(any(mod == "")){
                mod = mod[1:length(mod)-1]
                AA<-as.vector(unlist(mapply(add_mod_C, seq_sep, mod))) # AA as a AA vector
            } else{
                stop("Please contact developer. The function calculate_AAmz_individual_label is wrong.")
            }   
        }else{
            stop("Please contact developer. The function calculate_AAmz_individual_label is wrong.")
        }
        #browser()
        if(nchar(mod_N) > 0 ){AA <- c(mod_N, AA)}
    }else{
        AA <- unlist(seq_sep)
        if(nchar(mod_N)>0){ AA <- c(mod_N, AA) }
    }
    #calculate theoretical m/z of b ions and y ions
    AA_seq <- data.table::data.table(aa_varmod=AA, index=seq_len(length(AA)))
    #browser()
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
    #browser()
    
    # for charge == 1
    AA_mz$charge <- 1
    AA_mz_inc <- AA_mz[with(AA_mz, order(AA_mz$index)),]  # from N-terminal to C-terminal
    
    AA_mz_inc$mz_b <- cumsum(AA_mz_inc$weight) + H_weight
    AA_mz_inc$mz_c <- cumsum(AA_mz_inc$weight) + H_weight + NH3_weight

    # sort and calculate MW for y ion for
    AA_mz_dec <- AA_mz[with(AA_mz, order(-AA_mz$index)),]
    AA_mz_dec$mz_y <- cumsum(AA_mz_dec$weight) + H2O_weight + H_weight
    AA_mz_dec$mz_z0 <- cumsum(AA_mz_dec$weight) + H2O_weight + H_weight - NH3_weight + H_weight
    AA_mz_dec$mz_z1 <- cumsum(AA_mz_dec$weight) + H2O_weight + H_weight - NH3_weight + H_weight * 2
    AA_mz = merge(AA_mz_inc, AA_mz_dec) # for b/c/y/z ions
    #browser()
    if(any(grepl("neutralloss", colnames(AA_mz)))){
        if( sum(as.numeric(AA_mz$neutralloss) != 0) != 1){AA_mz$neutralloss <- NULL} # without neutral loss
        else{ # with neutral loss
            AA_mz$neutralloss = as.numeric(AA_mz$neutralloss)
            AA_mz$neutralloss = AA_mz$weight - AA_mz$neutralloss
            AA_mz_loss_inc <- AA_mz[with(AA_mz, order(AA_mz$index)),]  # from N-terminal to C-terminal
            
            AA_mz_loss_inc$mz_b_loss <- cumsum(AA_mz_loss_inc$neutralloss) + H_weight
            AA_mz_loss_inc$mz_c_loss <- cumsum(AA_mz_loss_inc$neutralloss) + H_weight + NH3_weight
            
            AA_mz_loss_inc$mz_b_loss <- (AA_mz_loss_inc$mz_b_loss != AA_mz_loss_inc$mz_b) * AA_mz_loss_inc$mz_b_loss
            AA_mz_loss_inc$mz_c_loss <- (AA_mz_loss_inc$mz_c_loss != AA_mz_loss_inc$mz_c) * AA_mz_loss_inc$mz_c_loss
            
            # sort and calculate MW for y ion for
            AA_mz_loss_dec <- AA_mz[with(AA_mz, order(-AA_mz$index)),]
            AA_mz_loss_dec$mz_y_loss <- cumsum(AA_mz_loss_dec$neutralloss) + H2O_weight + H_weight
            AA_mz_loss_dec$mz_z0_loss <- cumsum(AA_mz_loss_dec$neutralloss) + H2O_weight + H_weight - NH3_weight + H_weight
            AA_mz_loss_dec$mz_z1_loss <- cumsum(AA_mz_loss_dec$neutralloss) + H2O_weight + H_weight - NH3_weight + H_weight * 2
            
            AA_mz_loss_dec$mz_y_loss <- (AA_mz_loss_dec$mz_y_loss != AA_mz_loss_dec$mz_y) * AA_mz_loss_dec$mz_y_loss
            AA_mz_loss_dec$mz_z0_loss <- (AA_mz_loss_dec$mz_z0_loss != AA_mz_loss_dec$mz_z0) * AA_mz_loss_dec$mz_z0_loss
            AA_mz_loss_dec$mz_z1_loss <- (AA_mz_loss_dec$mz_z1_loss != AA_mz_loss_dec$mz_z1) * AA_mz_loss_dec$mz_z1_loss
            
            AA_mz = merge(AA_mz_loss_inc, AA_mz_loss_dec) # for b/c/y/z ions
        }
    }
    if(charge >2){ # only calculate the b/y ions with charge of 2
        AA_mz2 = AA_mz
        charge <- 2
        # adjusted, as AA_mz2_b$mz_b already add H_weight for charge 1
        AA_mz2$mz_b <- (AA_mz$mz_b + H_weight)/charge
        # adjusted, as AA_mz2_y$mz_y already add H2O_weight+H_weight for 1+
        AA_mz2$mz_y <- (AA_mz$mz_y + H_weight)/charge
        AA_mz2$mz_c <- (AA_mz$mz_c + H_weight)/charge
        AA_mz2$mz_z0 <- (AA_mz$mz_z0 + H_weight)/charge
        AA_mz2$mz_z1 <- (AA_mz$mz_z1 + H_weight)/charge
        AA_mz2$charge <- charge
        if( any(grepl("neutralloss", colnames(AA_mz)))){ # colnames contains neutral loss
            #browser()
            # adjusted, as AA_mz_loss2_b$mz_b_loss already add H_weight for charge 1
            AA_mz2$mz_b_loss <- (AA_mz$mz_b_loss + H_weight)/charge
            # adjusted, as AA_mz_loss2_y$mz_y already add H2O_weight+H_weight for 1+
            AA_mz2$mz_y_loss <- (AA_mz$mz_y_loss + H_weight)/charge
            AA_mz2$mz_c_loss <- (AA_mz$mz_c_loss + H_weight)/charge
            AA_mz2$mz_z0_loss <- (AA_mz$mz_z0_loss + H_weight)/charge
            AA_mz2$mz_z1_loss <- (AA_mz$mz_z1_loss + H_weight)/charge
        }
        AA_mz <- rbind(AA_mz, AA_mz2)
    }
    if( any(grepl("neutralloss", colnames(AA_mz)))){ # colnames contains neutral loss
        AA_mz$mz_b_loss_min <- (1-mz_range)*AA_mz$mz_b_loss
        AA_mz$mz_b_loss_max <- (1+mz_range)*AA_mz$mz_b_loss
        AA_mz$mz_c_loss_min <- (1-mz_range)*AA_mz$mz_c_loss
        AA_mz$mz_c_loss_max <- (1+mz_range)*AA_mz$mz_c_loss
        AA_mz$mz_y_loss_min <- (1-mz_range)*AA_mz$mz_y_loss
        AA_mz$mz_y_loss_max <- (1+mz_range)*AA_mz$mz_y_loss
        AA_mz$mz_z0_loss_min <- (1-mz_range)*AA_mz$mz_z0_loss
        AA_mz$mz_z0_loss_max <- (1+mz_range)*AA_mz$mz_z0_loss
        AA_mz$mz_z1_loss_min <- (1-mz_range)*AA_mz$mz_z1_loss
        AA_mz$mz_z1_loss_max <- (1+mz_range)*AA_mz$mz_z1_loss
    }
    AA_mz$mz_b_min <- (1-mz_range)*AA_mz$mz_b
    AA_mz$mz_b_max <- (1+mz_range)*AA_mz$mz_b
    AA_mz$mz_c_min <- (1-mz_range)*AA_mz$mz_c
    AA_mz$mz_c_max <- (1+mz_range)*AA_mz$mz_c
    AA_mz$mz_y_min <- (1-mz_range)*AA_mz$mz_y
    AA_mz$mz_y_max <- (1+mz_range)*AA_mz$mz_y
    AA_mz$mz_z0_min <- (1-mz_range)*AA_mz$mz_z0
    AA_mz$mz_z0_max <- (1+mz_range)*AA_mz$mz_z0
    AA_mz$mz_z1_min <- (1-mz_range)*AA_mz$mz_z1
    AA_mz$mz_z1_max <- (1+mz_range)*AA_mz$mz_z1
    #browser()

    if(ion_type == "y"){ # for b/y ions
        col_wo_cz = ! grepl("^mz_z|^mz_c", colnames(AA_mz)) # remove c/z ion annotations
        
        AA_mz <- subset(AA_mz, select = col_wo_cz)
    }else if(ion_type == "z"){ # for c/z ions
         col_wo_by = ! grepl("^mz_b|^mz_y", colnames(AA_mz)) # remove b/y ion annotations
        AA_mz <- subset(AA_mz, select = col_wo_by)
    }
    AA_mz <- AA_mz[order(AA_mz$charge,AA_mz$index)]
    #browser()
    return(AA_mz)
}

test_Ions <- function(AA_mz, mz_intensity_percent, ion_type, b_ion_col, y_ion_col){
    #browser()
    psm <- apply( mz_intensity_percent,1, test_individualIon, AA_mz, ion_type,
            b_ion_col, y_ion_col ) 
    psm <- data.table::rbindlist(psm) # can remove NULL elements
    if(nrow(psm)==0){stop(paste("The number of matched ion peaks is ", nrow(psm), " for the peptide '", 
        paste(AA_mz[[1]]$aa, collapse="") ,"'. Please check identification results from the search engine.", sep=""))}
    else if(nrow(psm)<4){warning(paste("The number of matched ion peaks is ", nrow(psm), " for the peptide '", 
        paste(AA_mz[[1]]$aa, collapse="") ,"', which is really limited. Please check 
        the ppm threshold or identification results from the search engine.", sep=""))}
    #browser()
    return(psm)
}

# Find matched b/y ions from MS2 peaks
test_individualIon<-function( mz_intensity_percent, AA_mz, ion_type, b_ion_col, y_ion_col){
    #browser()
    # intensity
    MIP <- data.frame(as.list(mz_intensity_percent))
    ion_info<-data.table::data.table()
    
    if(grepl("y", ion_type ) ){
        # search for b ion match
        ion <- subset(AA_mz, MIP$mz > AA_mz$mz_b_min & MIP$mz < AA_mz$mz_b_max)
        # mz matched and mz is smaller than the mw of the AA sequence (using index)
        if(nrow(ion) == 1 && ion$index < max(AA_mz$index)){
            #browser()
            b_max = 0.55
            b_min = abs(MIP$intensity_perc)
            if(abs(MIP$intensity_perc)>b_max){
                b_max = 1.01 # in case b_min = 100%
            }
            if(b_min < 0.1){
                b_min = b_min + 0.1    
            }
            abs_intensity_prc_ext = stats::runif(1,min=b_min, max=b_max)
            #browser()
            ion_info <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                "ionLabel" = paste("b", ion$index, paste(rep("+", ion$charge),
                    collapse=""), sep=""),
                "ion" = paste("b", ion$index, sep=""), "col" = b_ion_col,
                "direction" = ifelse(MIP$intensity_perc>0, 1, -1))
        }
        ion <- NULL
        #browser()
        # search for y ion match
        ion <- subset(AA_mz, MIP$mz > AA_mz$mz_y_min & MIP$mz < AA_mz$mz_y_max)
        # mz matched and mz is smaller than the mw of the AA sequence (using index)
        # index is for b index. calculate y ion by minux b index
        if(nrow(ion) == 1 && max(AA_mz$index)-ion$index+1 < max(AA_mz$index) ){
            #browser()
            y_max = 1.01
            y_min = abs(MIP$intensity_perc)
            if(y_min < 0.3){
                y_min = y_min + 0.2    
            }
            abs_intensity_prc_ext = stats::runif(1,min=y_min, max=y_max)
            #browser()
            yion <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                "abs_intensity_prc_ext" = abs_intensity_prc_ext,
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
        
        # for neutral loss
        if(any(grepl("loss", colnames(AA_mz)))){
            #browser()
            ion <- subset(AA_mz, MIP$mz > mz_b_loss_min & MIP$mz < AA_mz$mz_b_loss_max)
            # mz matched and mz is smaller than the mw of the AA sequence (using index)
            if(nrow(ion) == 1 && ion$index < max(AA_mz$index)){
                #browser()
                b_max = 0.55
                b_min = abs(MIP$intensity_perc)
                if(abs(MIP$intensity_perc)>b_max){
                    b_max = 1.01 # in case b_min = 100%
                }
                if(b_min < 0.1){
                    b_min = b_min + 0.1    
                }
                abs_intensity_prc_ext = stats::runif(1,min=b_min, max=b_max)
                #browser()
                bion_loss_info <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                    "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                    "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                    "ionLabel" = paste("b", ion$index, paste(rep("+", ion$charge),
                        collapse=""),"*", sep=""),
                    "ion" = paste("b", ion$index, sep=""), "col" = b_ion_col,
                    "direction" = ifelse(MIP$intensity_perc>0, 1, -1))
                if(nrow(ion_info)>0){
                    ion_info <- rbind(ion_info, bion_loss_info)
                }else{
                    ion_info <- bion_loss_info
                }
            }
            ion <- NULL
            
            # search for y ion match
            ion <- subset(AA_mz, MIP$mz > AA_mz$mz_y_loss_min & MIP$mz < AA_mz$mz_y_loss_max)
            #browser()
            # mz matched and mz is smaller than the mw of the AA sequence (using index)
            # index is for b index. calculate y ion by minux b index
            if(nrow(ion) == 1 && max(AA_mz$index)-ion$index+1 < max(AA_mz$index) ){
                #browser()
                y_max = 1.01
                y_min = abs(MIP$intensity_perc)
                if(y_min < 0.3){
                    y_min = y_min + 0.2    
                }
                abs_intensity_prc_ext = stats::runif(1,min=y_min, max=y_max)
                #browser()
                yion_loss_info <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                    "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                    "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                    "ionLabel" = paste("y", max(AA_mz$index)-ion$index+1,
                        paste(rep("+", ion$charge),collapse=""), "*",  sep=""),
                    "ion" = paste("y", max(AA_mz$index)-ion$index+1, sep=""),
                    "col" = y_ion_col,"direction" = ifelse(MIP$intensity_perc>0, 1, -1))
                if(nrow(ion_info)>0){
                    ion_info <- rbind(ion_info, yion_loss_info)
                }else{
                    ion_info <- yion_loss_info
                }
            }
            
        }
    }
    #browser()
    if(grepl("z", ion_type ) ){
        # search for b ion match
        ion <- subset(AA_mz, MIP$mz > AA_mz$mz_c_min & MIP$mz < AA_mz$mz_c_max)
        # mz matched and mz is smaller than the mw of the AA sequence (using index)
        if(nrow(ion) == 1 && ion$index < max(AA_mz$index)){
            #browser()
            c_max = 0.55
            c_min = abs(MIP$intensity_perc)
            if(abs(MIP$intensity_perc)>c_max){
                c_max = 1.01
            }
            if(c_min < 0.1){
                c_min = c_min + 0.1    
            }
            abs_intensity_prc_ext = stats::runif(1,min=c_min, max=c_max)
            #browser()
            cion <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                "ionLabel" = paste("c", ion$index, paste(rep("+", ion$charge),
                    collapse=""), sep=""),
                "ion" = paste("c", ion$index, sep=""), "col" = b_ion_col,
                "direction" = ifelse(MIP$intensity_perc>0, 1, -1))
                        
            if(nrow(ion_info)>0){
                ion_info <- rbind(ion_info, cion)
            }else{
                ion_info <- cion
            }        
        }
        ion <- NULL
    
        # search for z0 ion match
        ion <- subset(AA_mz, MIP$mz > AA_mz$mz_z0_min & MIP$mz < AA_mz$mz_z0_max)
        # mz matched and mz is smaller than the mw of the AA sequence (using index)
        # index is for b index. calculate y ion by minux b index
        if(nrow(ion) == 1 && max(AA_mz$index)-ion$index+1 < max(AA_mz$index) ){
            z_max = 1.01
            z_min = abs(MIP$intensity_perc)
            if(z_min < 0.3){
                z_min = z_min + 0.2    
            }
            abs_intensity_prc_ext = stats::runif(1,min=z_min, max=z_max)
            #browser()
            z0ion <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                "ionLabel" = paste("z0", max(AA_mz$index)-ion$index+1,
                    paste(rep("+", ion$charge),collapse=""), sep=""),
                "ion" = paste("z0", max(AA_mz$index)-ion$index+1, sep=""),
                "col" = y_ion_col,"direction" = ifelse(MIP$intensity_perc>0, 1, -1))
            
                        
            if(nrow(ion_info)>0){
                ion_info <- rbind(ion_info, z0ion)
            }else{
                ion_info <- z0ion
            }
        }

        # search for z1 ion match
        ion <- subset(AA_mz, MIP$mz > AA_mz$mz_z1_min & MIP$mz < AA_mz$mz_z1_max)
        # mz matched and mz is smaller than the mw of the AA sequence (using index)
        # index is for b index. calculate y ion by minux b index
        if(nrow(ion) == 1 && max(AA_mz$index)-ion$index+1 < max(AA_mz$index) ){
            #browser()
            z_max = 1.01
            z_min = abs(MIP$intensity_perc)
            if(z_min < 0.3){
                z_min = z_min + 0.2    
            }
            abs_intensity_prc_ext = stats::runif(1,min=z_min, max=z_max)
            #browser()
            z1ion <- data.table::data.table("mz"=MIP$mz, "index" = ion$index,
                "intensity"=MIP$intensity, "intensity_perc" =  MIP$intensity_perc,
                "abs_intensity_prc_ext" = abs_intensity_prc_ext,
                "ionLabel" = paste("z'", max(AA_mz$index)-ion$index+1,
                    paste(rep("+", ion$charge),collapse=""), sep=""),
                "ion" = paste("z'", max(AA_mz$index)-ion$index+1, sep=""),
                "col" = y_ion_col,"direction" = ifelse(MIP$intensity_perc>0, 1, -1))
            
            if(nrow(ion_info)>0){
                ion_info <- rbind(ion_info, z1ion)
            }else{
                ion_info <- z1ion
            }
        }
    }
    return(ion_info)
}

# find matched ions from MS2
find_matchedIons<-function(AA_mz, mz_intensity_percent, ion_type, b_ion_col, y_ion_col){
    #browser()
    if( length(AA_mz) == 1 ){ ## list(AA_mz) == 1
        psm <- apply( mz_intensity_percent,1,test_individualIon, AA_mz[[1]], ion_type,
            b_ion_col, y_ion_col ) 
        #browser()
        psm <- data.table::rbindlist(psm) # can remove NULL elements
        if(nrow(psm)==0){stop(paste("The number of matched ion peaks is ", nrow(psm), " for the peptide '", 
            paste(AA_mz[[1]]$aa, collapse="") ,"'. Please check identification results from the search engine.", sep=""))}
        else if(nrow(psm)<4){warning(paste("The number of matched ion peaks is ", nrow(psm), " for the peptide '", 
            paste(AA_mz[[1]]$aa, collapse="") ,"', which is really limited. Please check 
            the ppm threshold or identification results from the search engine.", sep=""))}
    }else{ ## list(AA_mz) >= 2; SILAC go this way
        #browser()
        psm_list <- lapply(AA_mz, test_Ions, mz_intensity_percent, ion_type, b_ion_col, y_ion_col ) 
        nrow_psm_df=unlist(lapply(psm_list, nrow)) # count rows of data.frame in list
        psm_max_row_df_pos=which(nrow_psm_df == max(nrow_psm_df)) # find the index of max 
        psm = psm_list[[psm_max_row_df_pos]]
        #psms <- mapply(test_Ions, AA_mz, MoreArgs = list( mz_intensity_percent, b_ion_col, y_ion_col ) )
    }

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
    #browser()
    ionLabel_comb <-
        stats::aggregate(ionLabel~mz+intensity_perc,
            psm, paste, collapse=",")
    psm<-merge(ionLabel_comb, psm,
        by=c("mz","intensity_perc"),
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
get_intensity_perc <- function(input_table, min_intensity_ratio){
    #browser()
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
