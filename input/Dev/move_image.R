

args <- commandArgs()
print(args)

ss <- as.numeric(args[3])
# ss <- 24

dir1 <- list.files('/projectnb/planet/PLSP/Img_raw/rawImage',full.names=T)
dir2 <- list.files('/projectnb/planet/PLSP/Img_raw/rawImage')
desd <- '/projectnb/modislc/users/mkmoon/Planet/rawImage/rawImage/'


desD <- paste0(desd,dir2[ss])
if (!dir.exists(desD)) {dir.create(desD)}

system(paste0('mv -v ',dir1[ss],'/* ',desD))