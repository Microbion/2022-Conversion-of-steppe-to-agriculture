pcnm_vpa <- function(otu, meta){
    require(SoDA)
    require(adespatial)
    require(vegan)
    xy <- geoXY(meta$lat, meta$lon, unit=1000)
    geo.trans <- data.frame(xy)
    otu.h.det <- resid(lm(as.matrix(otu)~., data=geo.trans))
    otu.d1 <- dist(otu.h.det)
    otu.correlog <- mantel.correlog(otu.d1, XY=geo.trans, nperm=99)
    geo.d1 <- dist(geo.trans)
    PCNM.auto <- dbmem(geo.d1,MEM.autocor="positive")
    PCNM.auto <- as.data.frame(PCNM.auto)
    undet.PCNM.rda <- rda(otu, PCNM.auto)
    undet.PCNM.R2a <- RsquareAdj(undet.PCNM.rda)$adj.r.squared
    undet.PCNM.fwd <- forward.sel(otu, as.matrix(PCNM.auto),adjR2thresh=undet.PCNM.R2a,nperm=99)
    #the number of significant PCNM
    nb.sig.PCNM <- nrow(undet.PCNM.fwd)
    # Identity of significant PCNMs in increasing order
    PCNM.sign <- sort(undet.PCNM.fwd$order)
    #Write significant PCNMs to a new object
    PCNM.red <- PCNM.auto[,c(PCNM.sign)]
    #Now that we have selected the best PCNM variables we run a new PCNM analysis with the 5 significant PCNM variables
    PCNM.rda2 <- rda(otu.h.det~., data=PCNM.red)
    fwd.R2a <- RsquareAdj(PCNM.rda2)$adj.r.squared
    XY.fwd <- forward.sel(otu, as.matrix(meta[,c("lat", "lon")]), adjR2thresh=fwd.R2a,nperm=99)
    XY.sign <- sort(XY.fwd$order)
    XY.red <- meta[,c("lat", "lon")][,c(XY.sign)]
    env.rda <- rda(otu, subset(meta, select=c(pH,TC.TN,TC,salt,AP,NO3_N,AK)))
    env.R2a <- RsquareAdj(env.rda)$adj.r.squared
    env.fwd <- forward.sel(otu,subset(meta, select=c(pH,TC.TN,TC,salt,AP,NO3_N,AK)), adjR2thresh=env.R2a, nperm=99 )
    env.sign <- sort(env.fwd$order)
    env.red <- subset(meta, select=c(pH,TC.TN,TC,salt,AP,NO3_N,AK))[,c(env.sign)]
    list(env.red, XY.red, PCNM.red)
}