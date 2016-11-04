# Viva 450K
# started by Allan Just
# importing Illumina 450K data
##############################

# mostly methylumi pipeline - start with iDat files

library(Hmisc)
library(methylumi)
library(wateRmelon)
library(ggplot2)
library(data.table)
library(tools) # for md5sum

importnotload <- T

#Path to data - will change later!
datadir <- "/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data"

## after running linux command to pull out all of the idat files not in base folders and move just those
# find IDATs_final/ -maxdepth 2 -name "*.idat" -exec cp {} idat_base_official/ \;
# /net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/IDATs_final
# [goa306@rclogin10 IDATs_final]$ ls | wc -l
# 117
# /net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/idat_base_official
# [goa306@rclogin10 idat_base_official]$ ls | wc -l
# 2804


##########################
if(importnotload){
  # barcodes from idat filenames
  barcodes <- list.files(path = file.path(datadir, "idat_base_official"), pattern = "idat$")
  barcodes <- unique(substring(barcodes, 1, 17))
  
  
  Sys.time()  
  system.time(raw.idat.all <- methylumIDAT(barcodes = barcodes, 
                                           parallel = F, oob = T, idatPath=file.path(datadir, "idat_base_official")))
  rm(barcodes)
  save(raw.idat.all, file = file.path(datadir, "imported", paste0("raw_idat_all_official.RData")))
} else if(!exists("raw.idat.all")) {
  load(file.path(datadir, "imported", paste0("raw_idat_all_official.RData")))
}

Sys.time()
# system.time(raw.idat.pf <- pfilter(raw.idat.all))


# cleanup a bit
rm(importnotload)
gc()


# End of file