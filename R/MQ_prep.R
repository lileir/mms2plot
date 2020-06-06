#' Conversion of Maxquant search result file (i.e. msms.txt) to the mms2plot-required identification file. 
#' @export MQ_prep
#' 
#' @param config_table_path File path of the table that contains the raw MS filepath, 
#'        the type of modifiation, mass shift and the tolerance of the mass shift. 
#'        Please see teh user guide for the details.
#' @param msms_path File path of the MaxQuant search result file (i.e. msms.txt).
#' @param output_file File path of of the parameter batch table that includes
#'        the parameter xml file and the fragment mass tolerance (ppm). The 
#'        parameter file format is referred to as par.xml in Maxquant.
#' @param output_path a character string naming a output file.
#' @return No value is returned.
#' 
#' @examples
#' \dontrun{
#' general_path = system.file( package = "mms2plot",dir = "extdata" )
#' setwd( general_path )
#' config_table_path = 'extdata/prep/user_table_forMSMS.txt'
#' msms_path = 'extdata/prep/MaxQuant/msms.txt'
#' output_file = 'extdata/prep/MaxQuant/conversion/identification.txt'
#' MQ_prep(config_table_path,msms_path,output_file)
#' }
#' 
#rm(list = ls())
#gc()
#library(data.table)
#library(stringr)
#
MQ_prep = function(config_table_path,msms_path,output_file){
    
    user_table = data.table::fread(config_table_path,data.table = FALSE)
    msms_read_0 = data.table::fread(msms_path,data.table = FALSE)
    ions_match = msms_read_0[,Find_cols(msms_read_0,rex = '[Mm]atches')]
    count_ions_butNOTneu = stringr::str_count(ions_match,';') + 1 - stringr::str_count(ions_match,'-')
    msms_read = msms_read_0[which(count_ions_butNOTneu > 6),]
    
    raws_user = basename(stringr::str_split(user_table$rawfilepath,'\\.',simplify = TRUE)[,1])###Raw(s) in usr's table
    raws = unique(msms_read$`Raw file`)###Raw(s) in msms.txt
    
    if(any(!raws %in% raws_user)){
        stop(paste0('Can not find the full path of "',
                    raws[!raws %in% raws_user],'.mzML" ! Please check the user table.'))
    }
    
    user_table_1 = user_table[raws_user %in% raws,]
    raws_user_1 = basename(stringr::str_split(user_table_1$rawfilepath,'\\.',simplify = TRUE)[,1])
    ls = lapply(raws,Score_fil,df = data.frame(msms_read))
    bind_score = data.frame(do.call(rbind,ls),check.names = FALSE)
    cols_fix = c('[Rr]aw.?[Ff]ile','[Ss]can.?[Nn]umber','[Mm]odifications?',
                 '[Ss]equence','[Mm]odified.?[Ss]equence','[Cc]harge')
    
    if(any(stringr::str_detect(colnames(bind_score),'^[Gg]ene.?[Nn]ames?$'))){
        cols_select = c(cols_fix,'[Gg]ene.?[Nn]ames?')
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = bind_score))
        sel = bind_score[,num_cols]
    }else if(any(stringr::str_detect(colnames(bind_score),'^[Pp]roteins?$'))){
        cols_select = c(cols_fix,'[Pp]roteins?')
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = bind_score))
        sel = bind_score[,num_cols]
    }else{
        cols_select = cols_fix
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = bind_score))
        sel = bind_score[,num_cols]
        sel$`Gene Names` = NA
    }
    colnames(sel) = c('Raw file','Scan number','Modifications','Sequence',
                      'Modified sequence','Charge','Gene Names')
    ###---Only pairs/groups---###
    msms_unmod = sel[with(sel,`Modifications` == 'Unmodified'),]
    msms_mod = sel[with(sel,`Modifications` != 'Unmodified'),]
    
    if(all(nrow(msms_unmod)!=0,nrow(msms_mod)!=0)){
        key_unmod = unique(paste(msms_unmod$`Raw file`,msms_unmod$Sequence,msms_unmod$Charge))
        key_mod = unique(paste(msms_mod$`Raw file`,msms_mod$Sequence,msms_mod$Charge))
        key_intersect = intersect(key_unmod,key_mod)
        
        if(length(key_intersect) != 0){
            keep_unmod = msms_unmod[with(msms_unmod,paste(`Raw file`,`Sequence`,`Charge`) %in% key_intersect),]
            keep_mod = msms_mod[with(msms_mod,paste(`Raw file`,`Sequence`,`Charge`) %in% key_intersect),]
            bind_groups = rbind(keep_mod,keep_unmod)
            keys = paste(bind_groups$`Raw file`,bind_groups$Sequence,bind_groups$Charge)
        }else{
            stop('There are no pairs/groups in msms.txt')
        }
    }else{
        stop('There are no pairs/groups in msms.txt')
    }
    #Add lables to msms.txt by raw_sequence_charge...
    label_ls = lapply(seq(length(key_intersect)),Add_label,df = bind_groups,keys = keys)
    label = do.call(rbind,label_ls)
    raw_ls = lapply(X = raws,Sub_raw,df1 = user_table_1,raws,df2 = label)
    bind = do.call(rbind,raw_ls)
    bind$`Modified sequence` = stringr::str_remove_all(bind$`Modified sequence`,'_')
    if(!dir.exists(dirname(output_file))){dir.create(dirname(output_file),recursive = TRUE)}
    utils::write.table(bind,output_file,quote = F,row.names = F,sep = '\t')
    print("The conversion is successfully complete.")
}

# config_table_path = 'extdata/user_table_forMSMS.txt'
# msms_path = 'extdata/MaxQuant/msms.txt'# output_file = 'extdata/MaxQuant/conversion/identification.txt'
# 
# 
# ####Run---MQ_prep
# MQ_prep(config_table_path,msms_path,output_file)
# 


