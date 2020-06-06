####***AddModLabel:add labels by raw_sequence***#########
Add_label = function(i,df,keys){
    df_sub = df[keys == unique(keys)[i],]
    df_sub$label = i
    return(df_sub)
}


Score_fil = function(raw,df){
    col_num = Find_cols(df,rex = '[Rr]aw.?file')
    df_1 = df[which(df[,col_num] == raw),]
    col_num_1 = Find_cols(df_1,rex = '[Ss]core')
    score = as.numeric(df_1[,col_num_1])
    df_2 = df_1[score >= stats::fivenum(score)[2],]
}


Sub_raw = function(raw,df1,df2,raws_user){
    user_sub = df1[which(raws_user == raw),]
    ###---***for one raw in raws***---###
    label_sub = df2[with(df2,`Raw file` == raw),]
    label_sub$`Raw file` = user_sub$rawfilepath
    #browser()
    label_sub$`Modified sequence` = gsub('\\(\\w+\\s?\\(\\w\\)\\)|\\(\\w+\\)',
                                         paste0('(',tolower(user_sub$modification),')'),
                                         label_sub$`Modified sequence`)
    return(label_sub)
}

Find_cols = function(df,rex){grep(paste0('^',rex,'$'),colnames(df))}

###***AddModLabel:add modification symbols into sequence***####
AddModLabel = function(i,seqs,poss,modlabel){
    seq = seqs[i]
    pos = as.numeric(poss[[i]])
    split = unlist(strsplit(seq,''))
    mod.res = split[pos]
    split[pos] = paste0(mod.res,modlabel)
    modseq = paste(split,collapse = '')
    return(modseq)
}

#Get indicies beyond the mw range of user-specified modification and fixed modification
Ind_modtype_err =  function(MW,flo,cei,c_flo,c_cei){
    mw = as.numeric(MW)
    ind_within = which(flo < mw & mw < cei | c_flo < mw & mw < c_cei)
}

#get the mw and position within range of fixed modification
Ind_mod_fix =  function(MW,c_flo,c_cei){
    mw = as.numeric(MW)
    ind_within = which(c_flo < mw & mw < c_cei)
}

#paste mws and positions
Ind_mod_fix1 = function(MW_POS,c_flo,c_cei){
    mw = as.numeric(unlist(stringr::str_extract_all(MW_POS,"[^\\(]\\d+\\.?\\d+[^\\)]")))
    pos = stringr::str_remove_all(string = unlist(stringr::str_extract_all(MW_POS,'\\(\\d+\\)')),pattern = '[\\(\\)]')
    ind_within = which(c_flo < mw & mw < c_cei)
    mw_1 = mw[-ind_within]
    pos_1 = pos[-ind_within]
    pas = paste(paste(mw_1,paste0('(',pos_1,')')),collapse=';')
    return(pas)
}


Convert_by_usr <- function(i,table){
    # i is the row number of table
    # table is user_table(dataframe)
    user_table = table[i,]
    raw_path = user_table$rawfilepath
    path_mzid_user = user_table$mzidpath
    ######******MZID --> flatten file
    mzid = mzID::mzID(path_mzid_user)
    fla = mzID::flatten(mzid)
    fla = fla[with(fla,rank == 1),]#comet,for MSGFPlus,almost every rank was '1'
    # browser()
    ######******filter modification mass shift in flatten file
    inv.col = c("[Ss]tart","[Ee]nd","[Pp]re","[Pp]ost","[Aa]ccession","[Dd]escription")
    num_cols = unlist(lapply(X = inv.col,FUN = Find_cols,df = fla))
    
    if(length(num_cols) > 0){
        fla = fla[!duplicated(fla[,-num_cols]),]
    }
    
    mod = subset(fla,fla$modified == TRUE)
    
    unmod = subset(fla,fla$modified == FALSE)
    
    ######******Extract all of mw and midified positions
    if( nrow(mod) > 0 ){
        ###***range of mw,only for user_table with single row
        accu_usr = user_table$mass_shift_tolerance
        mw_usr = user_table$mass_shift
        mw_floor = as.numeric(mw_usr) - as.numeric(accu_usr)
        mw_ceiling = as.numeric(mw_usr) + as.numeric(accu_usr)
        fix_c_mw = 57.021464
        c_fixed_floor = as.numeric(fix_c_mw) - as.numeric(accu_usr)
        c_fixed_ceiling = as.numeric(fix_c_mw) + as.numeric(accu_usr)
        
        mw_pos=mod$modification#mass shift and position
        extra_mw = stringr::str_extract_all(mw_pos,"[^\\(]\\d+\\.?\\d+[^\\)]")
        extra_pos = mapply(FUN = stringr::str_remove_all,
                           string = stringr::str_extract_all(mw_pos,'\\(\\d+\\)'),pattern = '[\\(\\)]')
        #####***Save mw within range [mw-variance,mw+variance]
        #remove mw not within modification specified by user and fixed modification
        inds_correct_ls = lapply(FUN = Ind_modtype_err,X = extra_mw,
                                 flo = mw_floor,cei = mw_ceiling,
                                 c_flo = c_fixed_floor,c_cei = c_fixed_ceiling)
        count_mw = unlist(lapply(extra_mw,length))
        count_correct = unlist(lapply(inds_correct_ls,length))
        inds_correct = which(count_mw == count_correct)
        mod_user_1 = mod[inds_correct,]#subset within mw range and fixed modification mw range
        if(nrow(mod_user_1)>0){
            mw_pos_1=mod_user_1$modification
            extra_mw_1 = stringr::str_extract_all(mw_pos_1,"[^\\(]\\d+\\.?\\d+[^\\)]")
            inds_fix_within = lapply(FUN=Ind_mod_fix,X = extra_mw_1,c_flo = c_fixed_floor,c_cei = c_fixed_ceiling)
            count_fix = unlist(lapply(inds_fix_within,length))
            inds_fix = which(count_fix != 0)
            ###^^^
            if(length(inds_fix) > 0){
                mod_fix_within = mod_user_1[inds_fix]#fix within
                mod_nonfix = mod_user_1[-inds_fix,]#no fix
                
                mw_pos_2=mod_fix_within$modification
                extra_mw_2 = stringr::str_extract_all(mw_pos_2,"[^\\(]\\d+\\.?\\d+[^\\)]")
                count_mw_2 = unlist(lapply(extra_mw_2,length))
                #When a peptide has only fixed modification,
                #the modification column will be changed to NA, modified to FALSE, and merged into the unmod subset
                ind_only_fix = which(count_fix[count_fix != 0] == count_mw_2)
                
                if(length(ind_only_fix) > 0){
                    mod_mul_types = mod_fix_within[-ind_only_fix,]#subset containing fixed and variable modifications simultaneously
                    mod_only_fix = mod_fix_within[ind_only_fix,]#subset has only fixed modification
                    
                    
                    mod_only_fix$modified = FALSE
                    mod_only_fix$modification = NA
                    
                    mw_pos_clean = unlist(lapply(X=mod_mul_types$modification,c_flo=c_fixed_floor,c_cei=c_fixed_ceiling,FUN=Ind_mod_fix1))
                    mod_mul_types$modification = mw_pos_clean
                    
                    mod_1 = rbind(mod_mul_types,mod_nonfix)#merge into mod subset
                    unmod_1 = rbind(mod_only_fix,unmod)#merge into unmod subset
                }else{
                    mod_mul_types = mod_fix_within
                    mw_pos_clean = unlist(lapply(X=mod_mul_types$modification,c_flo=c_fixed_floor,c_cei=c_fixed_ceiling,FUN=Ind_mod_fix1))
                    mod_mul_types$modification = mw_pos_clean
                    mod_1 = mod_mul_types
                    unmod_1 = unmod
                }
            }else{
                mod_1 = mod_user_1
                unmod_1 = unmod
                seq_intsect = intersect(unmod_1$pepseq,mod_1$pepseq)
                ind_sub = which(mod_1$pepseq %in% seq_intsect)
                mod_1_sub = mod_1[ind_sub,]
                unmod_sub = unmod_1[unmod_1$pepseq %in% seq_intsect,]
                ######******Extract mw and modified position from flatten file******######
                mw_pos_mod_1 = mod_1_sub$modification
                MWs = stringr::str_extract_all(mw_pos_mod_1,"[^\\(]\\d+\\.?\\d+[^\\)]")
                POSs = mapply(FUN = stringr::str_remove_all,
                              string = stringr::str_extract_all(mw_pos_mod_1,'\\(\\d+\\)'),pattern = '[\\(\\)]')
                ######******Get modified sequence with modification labels such as "(ox)"******######
                SEQ = mod_1_sub$pepseq
                modlabel = paste0('(',tolower(user_table$modification),')')
                
                Modified_seq = mapply(AddModLabel,seq(nrow(mod_1_sub)),MoreArgs = list(seqs=SEQ,poss=POSs,modlabel=modlabel))#modified sequences
                ######******Get the dataframe required by mms2plot******######
                unmod_sub$`Modified sequence` = unmod_sub$pepseq
                mod_1_sub$`Modified sequence` = Modified_seq
                unmod_sub$Modifications = 'Unmodified'
                mod_1_sub$Modifications = user_table$modification
                bind_mod_unmod = rbind(mod_1_sub,unmod_sub)
                bind_mod_unmod$`Raw file` = raw_path
                
                cols_pick = c('[Pp]roteins?','[Gg]ene.?names?','[Dd]escription','[Aa]ccession')
                
                num_cols = unlist(lapply(X = cols_pick,FUN = Find_cols,df = bind_mod_unmod))
                
                if(length(num_cols) > 0){
                    sel = bind_mod_unmod[,num_cols]
                    if(length(num_cols) > 1){
                        strpaste = apply(X = sel,MARGIN = 1,FUN = paste,collapse = '_')
                    }else{
                        strpaste = sel
                    }
                    GNs = stringr::str_remove(stringr::str_extract(strpaste,'GN=\\w+'),'GN=')
                    GNs[which(is.na(GNs))] = stringr::str_extract(strpaste[which(is.na(GNs))],'\\w{2}\\|\\w+\\|\\w+')
                    GNs[which(is.na(GNs))] = 'Uncharacterized'
                }else{GNs = 'Uncharacterized'}
                
                cols_select = c('Raw file','[Aa]cquisitionnum','Modifications','[Pp]epseq',
                                'Modified sequence','[Cc]hargestate')
                num_col = unlist(lapply(X = cols_select,FUN = Find_cols,df = bind_mod_unmod))
                df_mms2plot = cbind(dplyr::select(bind_mod_unmod,num_col),GNs)
                colnames(df_mms2plot) = c('Raw file','Scan number','Modifications','Sequence',
                                          'Modified sequence','Charge','Gene Names')
                return(df_mms2plot)
            }
            ###$$$  
        }else{
            stop('The mass shift exceeds the threshold , please check the input!!')
        }
    }else{
        stop('In the search results,there are no PSMs with variable modification that you are interested in...Please check the search process!')
    }  
}