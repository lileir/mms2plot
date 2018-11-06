## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ---- echo=TRUE, results='asis'------------------------------------------
getwd()
id_table_path = "../inst/extdata/silac/msms_SILAC.txt"
input_table <- data.table::fread(id_table_path, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
knitr::kable(input_table)

## ---- echo=TRUE, results='asis'------------------------------------------
getwd()
mqpar_filepath = "../inst/extdata/mqpar_batch.txt"
mqpar_files<-data.table::fread(mqpar_filepath, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
knitr::kable(mqpar_files)

## ----pressure, echo=FALSE, fig.cap="", out.width = '100%'----------------
knitr::include_graphics("annotation.png")

