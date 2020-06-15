#' Conversion of the research result with the mzid format to the mms2plot-required identification file. 
#' @export MZID_prep
#' 
#' @param config_table_path file path of the table that contains five columns: the raw MS filepath, 
#'        the type of modifiation, mass shift, the tolerance of the mass shift and the search result MZID filepath. 
#'        Please see teh user guide for the details.
#' @param output_path a character string naming a output file as the identification txt file for mms2plot.
#' @return No value is returned.
#' 
#' @examples
#' general_path = system.file( package = "mms2plot",dir = "extdata" )
#' setwd( general_path )
#' config_table_path = 'prep/user_table_forMSMS.txt'
#' path_user_table = 'prep/user_table_forMZID.txt'
#' output_path = 'prep/comet/conversion/identification.txt'
#' MZID_prep(path_user_table,output_path)#,dir_mzid)#label.by.Raw_seq_charge)
#' 
#' 
#' 
# rm(list=ls())
# gc()
# if(!requireNamespace('stringr')) install.packages('stringr')
# library(stringr)
# if(!requireNamespace('mzID')) install.packages('mzID')
# library(mzID)
# if(!requireNamespace('dplyr')) install.packages('dplyr')
# library(dplyr)
# 
# options(stringsAsFactors = FALSE)
# options(digits = 15)

MZID_prep = function(path_user_table,output_path){
  user_table = data.table::fread(path_user_table)#read user table
  result_df = lapply(X = seq(nrow(user_table)),FUN = Convert_by_usr,table = user_table)#,paths_mzid = paths_mzid)
  
  df = do.call(rbind,result_df)
  #browser()
  
  df_mod = df[with(df,Modifications == 'Unmodified'),]
  df_unmod = df[with(df,Modifications != 'Unmodified'),]
  
  key_unmod = unique(paste(df_unmod$`Raw file`,df_unmod$Sequence,df_unmod$Charge))
  key_mod = unique(paste(df_mod$`Raw file`,df_mod$Sequence,df_mod$Charge))
  key_intersect = intersect(key_unmod,key_mod)
  
  if(length(key_intersect) > 0){
      keep_unmod = df_unmod[with(df_unmod,paste(`Raw file`,`Sequence`,`Charge`) %in% key_intersect),]
      keep_mod = df_mod[with(df_mod,paste(`Raw file`,`Sequence`,`Charge`) %in% key_intersect),]
      bind_groups = rbind(keep_mod,keep_unmod)
  }else{
      bind_groups = df
  }
  
  keys = paste(bind_groups$`Raw file`,bind_groups$Sequence,bind_groups$Charge)
  df_sub_ls = lapply(seq(length(key_intersect)),Add_label,df = bind_groups,keys = keys)
  df_labeled = do.call(rbind,df_sub_ls)
  
  if(!dir.exists(dirname(output_path))){dir.create(dirname(output_path),recursive = TRUE)}
  utils::write.table(df_labeled,output_path,quote = F,row.names = F,sep = '\t')
  print("The conversion is successfully complete.")
}

# path_user_table = 'inst/extdata/prep/user_table_forMZID.txt'
# output_path = 'inst/extdata/prep/comet/conversion/identification.txt'
# ####Run---MZID_prep
# MZID_prep(path_user_table,output_path)
# path_user_table = 'E:\\mingliya\\project\\mms2plot\\search_result\\comet\\user_table_forMZID.txt'
# output_path = paste(dirname(path_user_table),'test_out.txt',sep = '/')
# MZID_prep(path_user_table,output_path)
