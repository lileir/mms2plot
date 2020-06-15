#' Conversion of Maxquant search result file (i.e. msms.txt) to the mms2plot-required identification file. 
#' @export MQ_prep
#' 
#' @param config_table_path File path of the table that contains four columns:the raw MS filepath, 
#'        the type of modifiation, mass shift and the tolerance of the mass shift. 
#'        Please see the user guide for the details.
#' @param msms_path File path of the MaxQuant search result file (i.e. msms.txt).
#' @param output_path a character string naming a output file as the identification txt file for mms2plot.
#' @return No value is returned.
#' 
#' @examples
#' general_path = system.file( package = "mms2plot",dir = "extdata" )
#' setwd( general_path )
#' config_table_path = 'prep/user_table_forMSMS.txt'
#' msms_path = 'prep/MaxQuant/msms.txt'
#' output_path = 'prep/MaxQuant/conversion/identification.txt'
#' MQ_prep(config_table_path,msms_path,output_path)
#'
#' 
#rm(list = ls())
#gc()
#library(data.table)
#library(stringr)

MQ_prep = function(config_table_path,msms_path,output_path){
    
    user_table = data.table::fread(config_table_path,data.table = FALSE)
    msms_read = data.table::fread(msms_path,data.table = FALSE)
    raws = unique(msms_read$`Raw file`)###Raw(s) in msms.txt
    
    if(stringr::str_detect(basename(user_table$rawfilepath),'\\.mzML')){
        raws_user = stringr::str_remove(basename(user_table$rawfilepath),'\\.mzML')###Raw(s) in usr's table
    }else{
        stop(paste0('Can not find the full path of raw file ...mzML ! Please check the user table.'))
    }
    
    if(any(!raws %in% raws_user)){
        stop(paste0('Can not find the full path of "',
                    raws[!raws %in% raws_user],'.mzML" ! Please check the user table.'))
    }
    ###>>>>>
    user_table_1 = user_table[raws_user %in% raws,]
    ###>>>>>    
    cols_fix = c('[Rr]aw.?[Ff]ile','[Ss]can.?[Nn]umber','[Mm]odifications?',
                 '[Ss]equence','[Mm]odified.?[Ss]equence','[Cc]harge')
    
    if(any(stringr::str_detect(colnames(msms_read),'^[Gg]ene.?[Nn]ames?$'))){
        cols_select = c(cols_fix,'[Gg]ene.?[Nn]ames?')
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = msms_read))
        sel = msms_read[,num_cols]
    }else if(any(stringr::str_detect(colnames(msms_read),'^[Pp]roteins?$'))){
        cols_select = c(cols_fix,'[Pp]roteins?')
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = msms_read))
        sel = msms_read[,num_cols]
    }else{
        cols_select = cols_fix
        num_cols = unlist(lapply(X = cols_select,FUN = Find_cols,df = msms_read))
        sel = msms_read[,num_cols]
        sel$`Gene Names` = NA
    }
    colnames(sel) = c('Raw file','Scan number','Modifications','Sequence',
                      'Modified sequence','Charge','Gene Names')
    
    keys = paste(sel$`Raw file`,sel$Sequence,sel$Charge,sep = '_')
    
    label_ls = lapply(seq(length(unique(keys))),Add_label,df = sel,keys = keys)
    label = do.call(rbind,label_ls)
    raw_ls = lapply(X = raws,Sub_raw,df1 = user_table_1,raws,df2 = label)
    bind = do.call(rbind,raw_ls)
    bind$`Modified sequence` = stringr::str_remove_all(bind$`Modified sequence`,'_')
    if(!dir.exists(dirname(output_path))){dir.create(dirname(output_path),recursive = TRUE)}
    utils::write.table(bind,output_path,quote = F,row.names = F,sep = '\t')
    print("The conversion is successfully complete.")
}

# ###For test >>>
# general = 'E:/mingliya/mms2plot-master/inst/extdata'
# setwd(general)
# config_table_path = 'prep/user_table_forMSMS.txt'
# msms_path = 'prep/MaxQuant/msms.txt'
# output_path = 'prep/MaxQuant/test.txt'
# # ####Run---MQ_prep
# MQ_prep(config_table_path,msms_path,output_path)
