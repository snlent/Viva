 # Viva 450K

# Preprocessing and normalization of methylation data objects

##############################

# library(Hmisc)
library(methylumi)
library(wateRmelon)
# library(ggplot2)
# library(RColorBrewer)
# library(plyr)
# library(matrixStats)
library(data.table)

# showplots <- F
# saveplots <- F

load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/datbgdb_official_CB_2015-12-16.Rdata")

# within-array normalization with BMIQ
##############################

# I was able to run BMIQ with all samples, with no problem
dat.bg.db.cord.bmiq <- BMIQ(dat.bg.db.cord)

save(dat.bg.db.cord.bmiq, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", "datbgdbbmiq_official_CB_", Sys.Date(), ".Rdata"))

