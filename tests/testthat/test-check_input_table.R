context("check input table")

test_that("check input table and path", {
  mqpar_files=readRDS(file = "mqpar_files.rds")
  mqpar_ppm=readRDS(file = "mqpar_ppm.rds")
  input_table=readRDS(file = "input_table.rds")
  par_xml_path=readRDS(file = "par_xml_path.rds")

  #saveRDS("modifications.xml", file = "par_xml_path.rds")

  list_aaMwModTable_ppm=readRDS(file = "list_aaMwModTable_ppm.rds")

  input_table$`Modified sequence` <-
    gsub("_","", input_table$`Modified sequence`)

  # run check_input_table()
  input_table_check <- check_input_table(input_table,id_table_path, mqpar_ppm)
  # test for check_input_table
  expect_equal(input_table, input_table_check)

  # # test for add_mod_aa
  # 
  #   MS2FileName<-unique(input_table$`Raw file`)[1]
  #   list_aaMwModTable_ppm_ama<-add_mod_aa(par_xml_path, basename(MS2FileName),
  #       mms2plot::aa_mw_table, mqpar_ppm)
  #   list_aaMwModTable_ppm=readRDS(file = "list_aaMwModTable_ppm.rds")
  #   expect_equal(list_aaMwModTable_ppm, list_aaMwModTable_ppm_ama)
})
