# label x/y axis
label_axis<-function(max_intensity, xlab, ylab, max_mz, cex,lwd){
    if(length(max_intensity) == 2){ # for mirror_plot
        if( max_intensity[1] > 0){ # quarter of the max_intensity
            axis_y_positive <- seq(0, max_intensity[1], length.out = 5)
            axis_y_negative <- seq(max_intensity[2], 0, length.out = 5)
        }else{# quarter of the max_intensity
            axis_y_positive <- seq(max_intensity[1], 0, length.out = 5)
            axis_y_negative <- seq(0, max_intensity[2], length.out = 5)
        }
        # remove duplicate "0" and keep 4 digits
        axis_y <- signif(c(axis_y_negative, axis_y_positive[-1]), digits=4)
        # abs( separate -1 to 1 with the length of 9 )
        axis_y_percent <- scales::percent(abs(seq(-1, 1, length.out = 9)))
        axis(side = 2,las = 1, at = seq(-1,1, length.out = 9),
            labels = as.integer(axis_y), tck=-0.01, mgp=c(0, 0.2, 0),
            cex.axis=0.33, lwd=0.5*lwd)
        axis(side = 4,las = 1, at = seq(-1,1, length.out = 9),
            labels = axis_y_percent,     tck=-0.01, mgp=c(0, 0.2, 0),
            cex.axis=0.33, lwd=0.5*lwd)
        axis(side = 1, at=seq(0, max_mz, by=200), lwd=0.5*lwd, tck=-0.01,
            labels=FALSE)
        # mgp=c(0,0.1,0) NOT work for dist. betw. label and x-axis. mtext used.
        graphics::mtext(side = 1, text=seq(0, max_mz, by=200),
            at=seq(0, max_mz, by=200), cex=0.33*cex) #
    }else{ # for group_plot
        axis_y <- seq(0, max_intensity, length.out = 5)
        axis_y_percent <- scales::percent(seq(0, 1, length.out = 5))
        axis(side = 2,las = 1, at = seq(0, 1, length.out = 5),
            labels = as.integer(axis_y), tck=-0.02, mgp=c(0, 0.2, 0),
            cex.axis=0.33, lwd=0.5*lwd)
        axis(side = 4,las = 1, at = seq(0, 1, length.out = 5),
            labels = axis_y_percent, tck=-0.02, mgp=c(0, 0.2, 0),
            cex.axis=0.33, lwd=0.5*lwd)
    }
}

generate_ionLabel <- function(ionLabel_peak, show_iontype){
    #browser()
    
#    ionLabel_peaks_split = unlist(strsplit(ionLabel_peak), ",")
#    for(i in 1:length(ionLabel_peaks_split)){
        # ionLabel_peak = ionLabel_peaks_split[i]
        # ion_name <- substr(ionLabel_peak, 1, regexpr(pattern ='\\+', ionLabel_peak)-1)
        # number_plus <- substring(ionLabel_peak, regexpr(pattern ='\\+',ionLabel_peak))
        # if(number_plus == "+"){number_plus = ""}
        # else if(number_plus == "+*"){number_plus="*"}
        # else if(number_plus == "++"){number_plus="2+"}
        # else if(number_plus == "++*"){number_plus="2+*"}
        # 
        # if(grepl("^z0", ion_name )){
        #     ion_number = substring(ion_name, 3)
        #     if(show_iontype){
        #         label = as.expression(bquote('z'^o*.(ion_number)^.(number_plus)))
        #     }else{
        #         label = as.expression(bquote(.(ion_number)^.(number_plus)))
        #     }
        # }else{
        #     if(show_iontype){
        #         label = as.expression(bquote(.(ion_name)^.(number_plus)))
        #     }else{
        #         ion_name = substring(ion_name, 2)
        #         label = as.expression(bquote(.(ion_name)^.(number_plus)))
        #     }
        # }
#    }
    show_iontype=T
 #   ionLabel_peak = "z08++,y7+"
 #   ionLabel_peak = "z08++"
    ionLabel_peak = gsub("z0", "w", ionLabel_peak)
    if(grepl(",", ionLabel_peak)){    
        ionLabel_peak = unlist(strsplit(ionLabel_peak, ","))
    }
    ionLabel <- substr(ionLabel_peak, 1, 1)
    ionLabel_style  <- ifelse(ionLabel == "w", "o", "")
    ionLabel = gsub("w", "z", ionLabel)
    ionLabel_number <- substr(ionLabel_peak, 2, regexpr(pattern ='\\+', ionLabel_peak)-1)
    number_plus <- substring(ionLabel_peak, regexpr(pattern ='\\+',ionLabel_peak))
    number_plus = gsub("^\\+\\+", "2+", number_plus)
    number_plus = gsub("^\\+", "", number_plus)
    
    if(length(ionLabel) == 1){
        label = as.expression(bquote(.(ionLabel)^.(ionLabel_style)*.(ionLabel_number)^.(number_plus)))
    }else if(length(ionLabel) == 2){
        label = as.expression(bquote(paste(.(ionLabel[1])^.(ionLabel_style[1])*.(ionLabel_number[1])^.(number_plus[2]),",",.(ionLabel[2])^.(ionLabel_style[2])*.(ionLabel_number[2])^.(number_plus[2]),sep="")))
        
#        label = as.expression(bquote(.(ionLabel[1])^.(ionLabel_style[1])*.(ionLabel_number[1])^.(number_plus[2]),.(ionLabel[2])^.(ionLabel_style[2])*.(ionLabel_number[2])^.(number_plus[2])))
    }

    # if(length(ionLabel_peak) == 1){
    #     if(grepl("^z0", ion_name )){
    #      ion_number = substring(ion_name, 3)
    #         if(show_iontype){
    #             label = as.expression(bquote('z'^o*.(ion_number)^.(number_plus)))
    #         }else{
    #             label = as.expression(bquote(.(ion_number)^.(number_plus)))
    #         }
    #     }else{
    #         if(show_iontype){
    #             label = as.expression(bquote(.(ion_name)^.(number_plus)))
    #         }else{
    #             ion_name = substring(ion_name, 2)
    #             label = as.expression(bquote(.(ion_name)^.(number_plus)))
    #         }
    #     }
    # }
    
    return(label)
}
# peaks annotation  plot(rnorm(30), xlab = as.expression(bquote(.(a)^th*.(b))))
####################
draw_peak_ionL<-function(psm, lwd, len_annoSpace, srt,show_iontype){
    #browser()
    len_annoSpace = 0.05 #  len_annoSpace / 2   distance between peak and label
    ion_unique<-unique(psm, select=c("mz", "abs_intensity_prc_ext", "direction",
        "intensity_perc","ionLabel_peak","col"))
    lines(ion_unique$mz, ion_unique$intensity_perc, type = "h",las=1,
        col=ion_unique$col, lwd = lwd/2)
    segments(ion_unique$mz,ion_unique$intensity_perc, ion_unique$mz,
        ion_unique$abs_intensity_prc_ext*ion_unique$direction,
        col=ion_unique$col, lwd = lwd/2, lty=3) # draw peak extention segments
    labels = mapply(generate_ionLabel, ion_unique$ionLabel_peak, show_iontype)

    text( ion_unique$mz, (ion_unique$abs_intensity_prc_ext+len_annoSpace)
        *ion_unique$direction,
        labels <- labels, cex = 0.3, col=ion_unique$col,adj=0.5,srt=0)
    #browser()
}

# generate and draw PSM labels
draw_psmanno<-function(AA_mz, PSM, max_mz, peptide_height, mod_height,
    len_annoSpace, y_ion_col, b_ion_col, lwd){
    #browser()
    AA_mz = AA_mz[[1]] # change a list of data.frame to data.frame
    # for peptide annotation, only consider aa with charge 1 for printing
    AA_mz <- subset(AA_mz, AA_mz$charge == 1)
    direction <- PSM$direction[1] # draw MS2 on upplot (1) or downplot (-1)
    peptide <- c("-", as.character(AA_mz$aa_varmod), "-")
    peptide_list <- vector("list", length(peptide)*2-1) # set peptide + two dash
    peptide_list_b <- peptide_list # preset NULL for b ions
    peptide_list_y <- peptide_list # preset NULL for y ions
    peptide_list_c <- peptide_list # preset NULL for c ions
    peptide_list_z0 <- peptide_list # preset NULL for z0 ions
    peptide_list_z1 <- peptide_list # preset NULL for z1 ions
    # set peptide
    peptide_list[c(TRUE,FALSE)] <- as.list(peptide)
    peptide_list[c(FALSE,TRUE)] <- as.list(".")
    # set b ions
    peptide_list_b[c(FALSE,TRUE)] <-
        paste("b", seq(0,length(peptide)-2,by=1),sep="")
    # set y ions
    peptide_list_y[c(FALSE,TRUE)] <-
        paste("y", rev(seq(0,length(peptide)-2,by=1)), sep="")

    # set c ions
    peptide_list_c[c(FALSE,TRUE)] <-
        paste("c", seq(0,length(peptide)-2,by=1),sep="")
    # set z ions
    peptide_list_z0[c(FALSE,TRUE)] <-
        paste("z0", rev(seq(0,length(peptide)-2,by=1)), sep="")
    peptide_list_z1[c(FALSE,TRUE)] <-
        paste("z'", rev(seq(0,length(peptide)-2,by=1)), sep="")

    # set the width of AA sequence in the plot
    x_quarter<-seq(0,max_mz, length.out = 20)[c(3,18)]  # seq from 5/20 to 17/20
    AA_pos <- seq(x_quarter[1], x_quarter[2], length.out=length(peptide_list))
    pos_start_dash <- x_quarter[1]  # x-axis value "-" at the front of peptides
    pos_end_dash   <- x_quarter[2] # x-axis value "-" at the end of peptides

    # integrate as data.table
    PSMlabel<-data.table::data.table(AA_pos=AA_pos, peptide=peptide_list,
        bion=peptide_list_b, yion=peptide_list_y, cion = peptide_list_c, 
        z0ion = peptide_list_z0, z1ion = peptide_list_z1)
    #browser()

    apply(PSMlabel, 1, print_modpeptide,  peptide_height, mod_height,
        pos_start_dash, pos_end_dash, direction)# print AA letters for labelling

    # draw b ions between AA letters
    PSManno_bion <- subset(PSMlabel, PSMlabel$bion %in% PSM$ion)$AA_pos
    index <- match(PSManno_bion,PSMlabel$AA_pos)-1
    bion_xsmall <- (PSMlabel$AA_pos[index] + PSManno_bion)/2
    PSManno_bsmall <-
        substring(subset(PSMlabel, PSMlabel$bion %in% PSM$ion)$bion,2)

        # draw y ions between AA letters
    PSManno_yion <- subset(PSMlabel, PSMlabel$yion %in% PSM$ion)$AA_pos
    index <- match(PSManno_yion,PSMlabel$AA_pos)+1
    yion_xlarge <- (PSMlabel$AA_pos[index] + PSManno_yion)/2
    PSManno_ysmall <-
        substring(subset(PSMlabel, PSMlabel$yion %in% PSM$ion)$yion,2)

    # draw c ions between AA letters
    PSManno_cion <- subset(PSMlabel, PSMlabel$cion %in% PSM$ion)$AA_pos
    index <- match(PSManno_cion,PSMlabel$AA_pos)-1
    cion_xsmall <- (PSMlabel$AA_pos[index] + PSManno_cion)/2
    PSManno_csmall <-
        substring(subset(PSMlabel, PSMlabel$cion %in% PSM$ion)$cion,2)
    
    # draw z0 ions between AA letters
    PSManno_z0ion <- subset(PSMlabel, PSMlabel$z0ion %in% PSM$ion)$AA_pos
    index <- match(PSManno_z0ion,PSMlabel$AA_pos)+1
    z0ion_xlarge <- (PSMlabel$AA_pos[index] + PSManno_z0ion)/2
    PSManno_z0small <-
        substring(subset(PSMlabel, PSMlabel$z0ion %in% PSM$ion)$z0ion,3) # removing z'

    # draw z1 ions between AA letters
    PSManno_z1ion <- subset(PSMlabel, PSMlabel$z1ion %in% PSM$ion)$AA_pos
    index <- match(PSManno_z1ion,PSMlabel$AA_pos)+1
    z1ion_xlarge <- (PSMlabel$AA_pos[index] + PSManno_z1ion)/2
    PSManno_z1small <-
        substring(subset(PSMlabel, PSMlabel$z1ion %in% PSM$ion)$z1ion,3) # removing z'

    #browser()
    if(direction ==1){
        if(length(PSManno_bion)>0){ # bion
            segments(PSManno_bion, rep(peptide_height, length(PSManno_bion)),
                PSManno_bion,rep(peptide_height-len_annoSpace,length(PSManno_bion)),
                col=b_ion_col, lwd=lwd)
            segments(bion_xsmall, rep(peptide_height-len_annoSpace,
                length(PSManno_bion)), PSManno_bion,
                rep(peptide_height-len_annoSpace, length(PSManno_bion)),
                col=b_ion_col, lwd=lwd)
    
            text((bion_xsmall+PSManno_bion)/2, rep(peptide_height-len_annoSpace,
                length(PSManno_bion)), PSManno_bsmall, cex = 0.25, adj=c(0.5, 1.3),
                col=b_ion_col)
        }
        if(length(PSManno_yion) > 0 ) {  # y ion
            segments(PSManno_yion, rep(peptide_height, length(PSManno_yion)),
                PSManno_yion,rep(peptide_height+len_annoSpace,length(PSManno_yion)),
                col=y_ion_col, lwd=lwd)
            segments(yion_xlarge,
                rep(peptide_height+len_annoSpace, length(PSManno_yion)),
                PSManno_yion, rep(peptide_height+len_annoSpace,
                length(PSManno_yion)), col=y_ion_col, lwd=lwd)
    
            text((yion_xlarge+PSManno_yion)/2,
                rep(peptide_height+len_annoSpace, length(PSManno_yion)),
                PSManno_ysmall, cex = 0.25, adj=c(0.5, -0.3),col=y_ion_col)
        }
        if(length(PSManno_cion)>0){ # c ion
            segments(PSManno_cion, rep(peptide_height, length(PSManno_cion)),
                PSManno_cion,rep(peptide_height-len_annoSpace,length(PSManno_cion)),
                col=b_ion_col, lwd=lwd)
            segments(cion_xsmall, rep(peptide_height-len_annoSpace,
                length(PSManno_cion)), PSManno_cion,
                rep(peptide_height-len_annoSpace, length(PSManno_cion)),
                col=b_ion_col, lwd=lwd)
    
            text((cion_xsmall+PSManno_cion)/2, rep(peptide_height-len_annoSpace,
                length(PSManno_cion)), PSManno_csmall, cex = 0.25, adj=c(0.5, 1.3),
                col=b_ion_col)
        }
        if(length(PSManno_z0ion) > 0 ) { #z0 ion
            segments(PSManno_z0ion, rep(peptide_height, length(PSManno_z0ion)),
                PSManno_z0ion,rep(peptide_height+len_annoSpace,length(PSManno_z0ion)),
                col=y_ion_col, lwd=lwd)
            segments(z0ion_xlarge,
                rep(peptide_height+len_annoSpace, length(PSManno_z0ion)),
                PSManno_z0ion, rep(peptide_height+len_annoSpace,
                length(PSManno_z0ion)), col=y_ion_col, lwd=lwd)
    
            text((z0ion_xlarge+PSManno_z0ion)/2,
                rep(peptide_height+len_annoSpace, length(PSManno_z0ion)),
                PSManno_z0small, cex = 0.25, adj=c(0.5, -0.3),col=y_ion_col)
        }
        if(length(PSManno_z1ion) > 0 ) { # z1 ion
            segments(PSManno_z1ion, rep(peptide_height, length(PSManno_z1ion)),
                PSManno_z1ion,rep(peptide_height+len_annoSpace,length(PSManno_z1ion)),
                col=y_ion_col, lwd=lwd)
            segments(z1ion_xlarge,
                rep(peptide_height+len_annoSpace, length(PSManno_z1ion)),
                PSManno_z1ion, rep(peptide_height+len_annoSpace,
                length(PSManno_z1ion)), col=y_ion_col, lwd=lwd)
    
            text((z1ion_xlarge+PSManno_z1ion)/2,
                rep(peptide_height+len_annoSpace, length(PSManno_z1ion)),
                PSManno_z1small, cex = 0.25, adj=c(0.5, -0.3),col=y_ion_col)
        }
    }else{
        if(length(PSManno_bion)>0){
            segments(PSManno_bion, rep(peptide_height*-1, length(PSManno_bion)),
                PSManno_bion,
                rep((peptide_height+len_annoSpace)*-1, length(PSManno_bion)),
                col=b_ion_col, lwd=lwd)
            segments(bion_xsmall,
                rep((peptide_height+len_annoSpace)*-1, length(PSManno_bion)),
                PSManno_bion, rep((peptide_height+len_annoSpace)*-1,
                length(PSManno_bion)), col=b_ion_col, lwd=lwd)
    
            text((bion_xsmall+PSManno_bion)/2,rep((peptide_height+len_annoSpace)*-1,
                length(PSManno_bion)), PSManno_bsmall, cex = 0.25, adj=c(0.5, 1.3),
                col=b_ion_col)
        }
        if(length(PSManno_yion) > 0 ) {
            segments(PSManno_yion, rep(peptide_height*-1, length(PSManno_yion)),
                PSManno_yion, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_yion)), col=y_ion_col, lwd=lwd)
            segments(yion_xlarge, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_yion)), PSManno_yion,
                rep((peptide_height-len_annoSpace)*-1, length(PSManno_yion)),
                col=y_ion_col, lwd=lwd)
    
            text((yion_xlarge+PSManno_yion)/2,rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_yion)), PSManno_ysmall, cex = 0.25, adj=c(0.5, -0.3),
                col=y_ion_col)
        }
        if(length(PSManno_cion)>0){
            segments(PSManno_cion, rep(peptide_height*-1, length(PSManno_cion)),
                PSManno_cion,
                rep((peptide_height+len_annoSpace)*-1, length(PSManno_cion)),
                col=b_ion_col, lwd=lwd)
            segments(cion_xsmall,
                rep((peptide_height+len_annoSpace)*-1, length(PSManno_cion)),
                PSManno_cion, rep((peptide_height+len_annoSpace)*-1,
                length(PSManno_cion)), col=b_ion_col, lwd=lwd)
    
            text((cion_xsmall+PSManno_cion)/2,rep((peptide_height+len_annoSpace)*-1,
                length(PSManno_cion)), PSManno_csmall, cex = 0.25, adj=c(0.5, 1.3),
                col=b_ion_col)
        }
        if(length(PSManno_z0ion) > 0 ) {
            segments(PSManno_z0ion, rep(peptide_height*-1, length(PSManno_z0ion)),
                PSManno_z0ion, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z0ion)), col=y_ion_col, lwd=lwd)
            segments(z0ion_xlarge, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z0ion)), PSManno_z0ion,
                rep((peptide_height-len_annoSpace)*-1, length(PSManno_z0ion)),
                col=y_ion_col, lwd=lwd)
    
            text((z0ion_xlarge+PSManno_z0ion)/2,rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z0ion)), PSManno_z0small, cex = 0.25, adj=c(0.5, -0.3),
                col=y_ion_col)
        }
        if(length(PSManno_z1ion) > 0 ) {
            segments(PSManno_z1ion, rep(peptide_height*-1, length(PSManno_z1ion)),
                PSManno_z1ion, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z1ion)), col=y_ion_col, lwd=lwd)
            segments(z1ion_xlarge, rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z1ion)), PSManno_z1ion,
                rep((peptide_height-len_annoSpace)*-1, length(PSManno_z1ion)),
                col=y_ion_col, lwd=lwd)
    
            text((z1ion_xlarge+PSManno_z1ion)/2,rep((peptide_height-len_annoSpace)*-1,
                length(PSManno_z1ion)), PSManno_z1small, cex = 0.25, adj=c(0.5, -0.3),
                col=y_ion_col)
        }
    }
}

# print identified peptide
print_modpeptide <- function(PSMlabel, peptide_height, mod_height,
    pos_start_dash, pos_end_dash, direction){
    #browser()
    peptide_height <- peptide_height * direction
    mod_height  <- mod_height *  direction
    if (PSMlabel$peptide == ".") {
        return(NULL)# ignore
    } else if (grepl("^[A-Z]", PSMlabel$peptide) ||
        PSMlabel$peptide == "-") {
        # AA at the front
        if (nchar(PSMlabel$peptide) == 1) {
            # without modification
            text( PSMlabel$AA_pos, peptide_height, PSMlabel$peptide,
                cex = 0.5, adj = 0.5 )
        } else{
            # with modification at the end
            AA <- substr(PSMlabel$peptide, 1, 1)
            AA_endMod <-
                substr(PSMlabel$peptide, 3, nchar(PSMlabel$peptide) - 1)
            text( PSMlabel$AA_pos, peptide_height, bquote(.(AA)),
                cex = 0.5, adj = 0.5 )
            if (direction > 0) {
                text( PSMlabel$AA_pos, peptide_height + mod_height,
                    bquote(""[.(AA_endMod)]), cex = 0.5, adj = 0.5 )
            } else{
                text( PSMlabel$AA_pos, peptide_height - mod_height,
                    bquote(""[.(AA_endMod)]), cex = 0.5, adj = 0.5 )
            }
        }
    } else if (grepl("^\\(", PSMlabel$peptide)) {
        # mod at the start
        #browser()
        bracket_start <- gregexpr("(", PSMlabel$peptide, fixed = TRUE)[[1]]
        if (length(bracket_start) == 1) {
            # only one mod which is at the start
            AA_frontMod <-
                substr(PSMlabel$peptide, 2, nchar(PSMlabel$peptide) - 2)
            AA <- substring(PSMlabel$peptide, nchar(PSMlabel$peptide))
            if (direction > 0) {
                text( pos_start_dash, peptide_height + mod_height,
                    bquote(""[.(AA_frontMod)]), cex = 0.5, adj = 0.5 )
                text( PSMlabel$AA_pos, peptide_height, bquote(.(AA)),
                    cex = 0.5, adj = 0.5 )
            }else{
                text( pos_start_dash, peptide_height - mod_height,
                    bquote(""[.(AA_frontMod)]), cex = 0.5, adj = 0.5 )
                text(PSMlabel$AA_pos, peptide_height, bquote(.(AA)), cex = 0.5,
                    adj = 0.5 )
            }
        } else if (length(bracket_start) == 2) {
            # two mod: one at the front and another at the end
            AA_frontMod <- substr(PSMlabel$peptide, 2, bracket_start[2] - 3)
            AA <-
                substr(PSMlabel$peptide, bracket_start[2]-1, bracket_start[2]-1)
            AA_endMod <-
                substr(PSMlabel$peptide, bracket_start[2] + 1,
                    nchar(PSMlabel$peptide) - 1)
            text(PSMlabel$AA_pos, peptide_height, bquote(.(AA)), cex = 0.5,
                adj = 0.5)
            if (direction > 0) {
                text( pos_start_dash, peptide_height + mod_height,
                    bquote(""[.(AA_frontMod)]), cex = 0.5, adj = 0.5 )
                text( PSMlabel$AA_pos, peptide_height + mod_height,
                    bquote(""[.(AA_endMod)]), cex = 0.5, adj = 0.5 )
            } else{
                text( pos_start_dash, peptide_height - mod_height,
                    bquote(""[.(AA_frontMod)]), cex = 0.5, adj = 0.5 )
                text(PSMlabel$AA_pos,peptide_height - mod_height,
                    bquote(""[.(AA_endMod)]), cex = 0.5, adj = 0.5)
            }
        } else{
            stop("CAN NOT BE THREE MODS! [note:stopped in print_modpeptide()].")
        }
    } else{
        stop("CAN NOT BE THREE MODS! [note:stopped in print_modpeptide()].")
    }
}


# write information about file, protein, MS2 retention time, charge and etc.
draw_ms2generalinfo<-function(rt, scan, mz, charge, gene, PSM, info_height){
    direction <- PSM$direction[1]
    rt <- paste("RT (min):", ifelse(rt<480, rt, round(rt/60, 1)))
    scan <- paste("Scan:", scan)
    mz <- paste("m/z:", round(mz, 2))
    charge <- paste("Charge:", charge)
    gene <- paste("Gene:", gene)
    text(25, info_height*direction, paste(scan, mz, charge, rt, gene, sep="; "),
        pos=4, offset=0, cex = 0.33)
}

#  plot_group_individual
plot_group_individual<-function(mz_intensity_percent, AA_mz,  PSM,
    max_intensity, base_rawFile, rt, scan, mz, charge, gene, y_ion_col,
    b_ion_col, peaks_col, ymax, lwd, peptide_height, mod_height, len_annoSpace,
    srt, cex, max_mz, info_height, show_iontype){
    #browser()
    graphics::plot(mz_intensity_percent$mz, mz_intensity_percent$intensity_perc,
        type = "h",las = 1, xlab = "", ylab = "", xlim = c(0, max_mz), xaxt="n",
        yaxt="n", xaxs="i", yaxs="i",  # w/o using axes=FALSE as want x/ylim
        ylim = c(0, ymax), col = peaks_col, cex = cex, lwd=0.5,
        frame.plot=FALSE) # max_intensity：1. give space of 0.7 for annotation
    graphics::box(lwd=lwd/2)

    label_axis(max_intensity, xlab="", ylab="", max_mz, cex, lwd)
    # draw and color peak extension line
    draw_peak_ionL(PSM, lwd, len_annoSpace, srt, show_iontype)
    ## draw PSM annotation (peptide, mod, b ions and y ions)
    draw_psmanno( AA_mz, PSM, max_mz, peptide_height, mod_height, len_annoSpace,
        y_ion_col, b_ion_col, lwd=lwd/2 )
    ## MS2 information ( file name, RT, Scan number, m/z, charge and gene name)
    draw_ms2generalinfo( rt, scan, mz, charge, gene, PSM, info_height)
}
