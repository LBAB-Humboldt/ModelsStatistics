
# 1. Revisar que las proyecciones coinciden
CheckProjection <- function(r, proj) {
      if(is.na(testN2@crs@projargs)){
    testN2@crs@projargs <- proj
    cat('sin proyeccion')
  }
    if (testN2@crs@projargs != proj){
    if (testN2@crs@projargs == "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
    { testN2@crs@projargs <- proj
    cat('sin datum')
    } 
    else {
      testN2 <- projectRaster(testN2, crs = proj)
      compareRaster(testN2, clcBrick, extent = TRUE, rowcol = TRUE, crs = TRUE)
      cat('Resampling ... \n' )
    }
  }
}


# 2. Calcular el tamanno de ocurrencia de la espcie km2
EstimateRangeSize <- function(r, area.raster) {
  outArea_fx <- cellStats( r * area.raster, stat = "sum")
  return(outArea_fx)
}



# 3. Minimo poligono convexo para Raster
PredictMcpRast <- function(r, area.raster, rast) {
  cells <- Which(r == 1, cells = T)
  pts <-  xyFromCell(r, cell = cells) #Calcula el centro del pixel
  ch <- convHull(pts) #Predice que la especie este dentro de determinado convex hull con base en los puntos de entrenamiento
  ch.pred <- predict(ch, r)
  ch.area <- cellStats(ch.pred * area.raster, "sum")
  
  if (rast == T) {
    writeRaster(ch.pred, paste0(sppFolder,'/mcp_Rast.tif'), overwrite=TRUE)
  }
 
  return (list (ch.area = ch.area, ch = ch.pred))
}



# 4. Minimo poligono convexo para Registros de especies
PredictMcpRec <- function (r, Proj_col, rast, s, sppPath, proj) {
  recc <- read.csv(as.character(sppPath$Rec_file[s]), as.is = T)
  lon_lat <- matrix(nrow = nrow(recc), ncol = 2)
  lon_lat[,1] <- recc$lon
  lon_lat[,2] <- recc$lat
  in.pts <- SpatialPoints(lon_lat, proj4string = CRS(proj))
  cells <- unique(cellFromXY(r, in.pts))
  if(length(cells) > 2){ 
    ch_pol <- convHull(lon_lat)
    ch_pol@polygons@proj4string@projargs <- proj
    ch_proj <- spTransform(ch_pol@polygons, CRS("+init=epsg:3116")) ##EPSG 3116 codigo Magna origen Bogota
    ch_proj <- spTransform(ch_pol@polygons, Proj_col)
    MCPrec_km2 <- gArea(ch_proj) / 1000000
    return (list (MCP_recp = ch_proj, MCPrec_km2 = MCPrec_km2)) 
  } else {
    return (list (MCP_recp = NA, MCPrec_km2 = NA))
  }
}


# 5. Habitat CLC
HabitatAreaClc<- function (r, area.raster, clcBrick, clcKey, rast){
  r <- r * area.raster 
  clc_area <- rep( 0, nlayers(clcBrick))
  clc_area <- data.frame(matrix(0,nrow=1,ncol=nlayers(clcBrick)))
  colnames(clc_area) <- gsub('.tif', '', clcTif)
  for(l in 1:nlayers(clcBrick)){
    clcSp <- ((r * clcBrick[[l]] /100)) # Se divide por 100 ya que los valores de clcBrick estan en porcentaje, r contien los valores del raster con respecto al tama?o del pixel correspondiente
    clc_area[1,l] <- cellStats(clcSp, stat ="sum", na.rm = TRUE)
    }
  return(clc_area)
}


# 6. Historico perdida de bosque
ForestLoss <- function(r, area.raster, f90, f00, f05, f10, f12, rast=TRUE) { 
  r <- r *area.raster
  com90 <- r * f90 # Los valores f90, f00, f05, f10, f12, corresponden al area en km2 de bosque en el pixel. El rango es de 0 a 1, siendo 1km2 el area total del pixel.
  com00 <- r * f00
  com05 <- r * f05
  com10 <- r * f10
  com12 <- r * f12
  
  forestArea <- data.frame(y1990 = cellStats(com90, sum),
                             y2000 = cellStats(com00, sum),
                             y2005 = cellStats(com05, sum),
                             y2010 = cellStats(com10, sum),
                             y2012 = cellStats(com12, sum))
  return(forestArea)
}


# 7. Escenarios 2030
Escenarios2030 <- function (r, area.raster, con2030_res,  des2030_res, his2030_res, sppFolder, rast) {
  inRaster <- r *area.raster
  escenarios2030 <- data.frame (con_2030 = NA, des_20302 = NA, his_2030 = NA) 
  
  #con2030
  inRaster_2030c <- con2030_res * inRaster
  escenarios2030[1,1] <- cellStats (inRaster_2030c, stat = "sum")
  
  #des2030
  inRaster_2030d <- des2030_res * inRaster
  escenarios2030[1,2] <- cellStats (inRaster_2030d, stat = "sum")
  
  #his2030
  inRaster_2030h <- his2030_res * inRaster
  escenarios2030[1,3] <-cellStats (inRaster_2030h, stat = "sum")
  
  return (escenarios2030)
}


# 8. Areas protegidas
AllProtectedFigures <- function(r, ap, pn, sc, ot, area.raster, rast){
  r <- r * area.raster
  spOccEx <- EstimateRangeSize( r = testN2, area.raster = area.raster) ##.............. Deber?a dividirse por el area dentro del poligono
  # Los calulos del la interseccion entre el modelo y el area protegida se dividen por cien porque los valores de pn est?n en porcentaje rango de 1 a 100
  # Se divide por cien porque los valores de pn est?n en porcentaje rango de 1 a 100
  intAp <- (ap * r /100)
  areaAp <- sum(intAp[], na.rm = TRUE)
  intPn <- (pn * r /100)
  areaPn <- sum(intPn[], na.rm = TRUE) #Suma de las probabilidades de ocurrencia dentro de 'pn' para el calculo de las estadisticas
  intSc <- (sc * r / 100)
  areaSc <- sum(intSc[], na.rm = TRUE)
  intOt <- (ot * r /100)
  areaOt <- sum(intOt[], na.rm = TRUE)
  repr <- data.frame(protArea = areaAp*100/spOccEx, #Proporcion del area de distribucion entre el area protegida y el total
                     natPark =  areaPn*100/spOccEx,
                     civilRes = areaSc*100/spOccEx,
                     othFigur = areaOt*100/spOccEx)
  
  if (rast==T) {
    rastBrc <- stack(intAp, intPn, intSc, intOt)
    print('TIFF protected areas')
    writeRaster(rastBrc, paste0(sppFolder,'/sinap.tif'), overwrite=TRUE)
  }
  return(repr)
}



# AMENAZAS 
# 9.human foot print hfp

HumanFootPrint <- function (r, hfp, rast){
  nhfp <- hfp * r
  sumVal  <- sum(nhfp[], na.rm = TRUE)
  minVal  <- min(nhfp[], na.rm = TRUE)
  maxVal  <- max(nhfp[], na.rm = TRUE)
  meanVal <- mean(nhfp[], na.rm = TRUE)
  varVal  <- var(nhfp[], na.rm = TRUE)
  cvVal   <- cv(nhfp[], na.rm = TRUE)
  size <- length(which(!is.na(nhfp[])))
  
  dfResult <- data.frame(size, sumVal, meanVal, minVal, maxVal, varVal, cvVal)  
  
  if (rast == T) {
    writeRaster(nhfp, paste0(sppFolder,'/hfp.tif'), overwrite=TRUE)
  }
    return(df = dfResult)
  
}

# 10. Vias/Roads
RoadsLength <- function(r, rds, rast){
  nrds <- rds * r
  size <- length(which(!is.na(nrds[])))
  if (rast == T) {
    writeRaster(nrds, paste0(sppFolder,'/Roads.tif'), overwrite=TRUE)
  }
      return(size)
  
}

# 11. Titulos mineros
TitulosMineros <- function(r, area.raster, ntif, rast) {
  nTit <- (titMin * r * area.raster /100)
  #size <- length(which(!is.na(nTit[]))) ---Evaluar si se requiere
  Tit_area  <- cellStats(nTit, stat = "sum", na.rm = TRUE)
  
  if (rast == T) {
    print ('Raster titulos mineros')
    writeRaster(nTit, paste0(sppFolder,'/Tit_min.tif'), overwrite=TRUE)
  }
  return(Tit_area)
}


