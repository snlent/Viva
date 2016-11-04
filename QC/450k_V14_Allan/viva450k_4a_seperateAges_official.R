
# Seperate out the stages of samples: Cord blood, Y3, Y7
# Otherwise, BMIQ-ing everything all together take a million Gs and hours - it times out/gets killed.

# this was doen interactively, not as a batch job.

library(methylumi)
library(data.table)


load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/datbgdb_official_2015-12-14.Rdata")


# save out CB
dat.bg.db.cord <- dat.bg.db[, pData(dat.bg.db)$COLLECTION =="Proband Delivery"]
dim(dat.bg.db.cord)
# Features  Samples
# 485577      549

save(dat.bg.db.cord, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", 
                                   "datbgdb_official_CB_", Sys.Date(), ".Rdata"))

rm(dat.bg.db.cord)


# DO THE REST LATER - IT TAKES INCREDIBLY LONG TO SAVE OUT EACH FILE


# save out Y3
dat.bg.db.Y3 <- dat.bg.db[, pData(dat.bg.db)$COLLECTION =="Proband Age 3"]
dim(dat.bg.db.Y3)

save(dat.bg.db.Y3, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", 
                                   "datbgdb_official_Y3_", Sys.Date(), ".Rdata"))

rm(dat.bg.db.Y3)


# save out Y7
dat.bg.db.Y7 <- dat.bg.db[, pData(dat.bg.db)$COLLECTION =="Proband Age 7"]
dim(dat.bg.db.Y7)

save(dat.bg.db.Y7, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", 
                                 "datbgdb_official_Y7_", Sys.Date(), ".Rdata"))

rm(dat.bg.db.Y7)
