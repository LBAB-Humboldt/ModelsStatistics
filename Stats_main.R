##  ---------------- ESTAD?STICAS MODELOS DE DISTRIBUCION--------------------------------
## Este script genera 78 arhivos (outputs) por modelo que corresponden a:
## 54 TIFF de coberturas de la tierra nivel 3, rango de 0 a 1 (km2)
## 5 TIFF de Bosque para 1990, 2000, 2005, 2010, 2012, rango de 0 a 1 (km2)
## 1 TIFF Human Foot Print, ragno 0 a ??????
## 1 Shapefile (8 archivos) de minConvexPolygon_pcs
## 2 TIFF; (1) Minning, rango 0 a 1 (km2); (2) Roads, rango 0 a 1
## 4 TIFF ProtArea (1) _allProtArea; (2) _othFigures; (3)_NatPark; (4)_civilRe; con rangos de 0 a 100
## 4 CSV; StatA; StatB; StatC; StatD



# Cargar librerias
library(rgdal)
library(raster)
library(dismo)
library(rgeos)

# Cargar capas y scripts
source('//192.168.11.113/Lab_biogeografia2/Estadisticas/stats_SDM.R') # Ajustar ruta segun las carpertas del PC

# Existe un script 'Generate_Rdata.R' para generar los siguientes archivos Rdata
load('//192.168.11.113/Lab_biogeografia2/Estadisticas/clc_n3.RData') # Archivos localizados en la NAS, ajustar si estan en otro directorio
load('//192.168.11.113/Lab_biogeografia2/Estadisticas/forestTifs.RData')
load('//192.168.11.113/Lab_biogeografia2/Estadisticas/Escenarios_2030.RData')
load('//192.168.11.113/Lab_biogeografia2/Estadisticas/protAreaTifs.RData') 
load('//192.168.11.113/Lab_biogeografia2/Estadisticas/rds&Min.RData') 


# Lista de modelos que van a ser evaluados, en este caso Zamias
sppPath <- read.csv('//192.168.11.113/Lab_biogeografia2/Modelos/Zamias/finalRasterList.csv', sep=",", header = TRUE, as.is = T)
sppRecords <- ('//192.168.11.113/Lab_biogeografia2/Modelos/Zamias/21112016_Zamias2/') # Diretorio de los archivos csv con los valores de registro o puntos geograficos 'records'
sppPath$Rec_file <- paste0(sppRecords, sppPath$Species,'.csv')
sppName <- paste0(gsub('.tif', '', basename(as.character(sppPath[, 1]))))


# Set proyecciones para los calculos
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Proj_col <- "+proj=tmerc +lat_0=4.596200417 +lon_0=-74.07750791700001 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +units=m +no_defs"


# Tablas de resultados
R_size <- data.frame()
metadata.all <- data.frame()
header.met <- read.csv('//192.168.11.113/Lab_biogeografia2/Estadisticas/Indice_metadatos.csv')
threats.all <- data.frame()


for(s in 21:nrow(sppPath)){

# Direccion de salida outputs
sppFolder <- paste0('C:/Biomodelos/Grilla/Zamias_2504/', gsub('.tif', '', basename(as.character(sppPath[s, 2])))) # Establecer el directorio de salida
dir.create(sppFolder, recursive = TRUE)
  
  
#Generar capas solo para el primer registro
if (s==1) {rast = TRUE} else {rast = FALSE}
testN2 <- raster(as.character(sppPath[s, 2]))
cat(s)

## Funciones calculos estadisticas
# 1. Revisar que las proyecciones coinciden
CheckProjection(r = testN2, proj=proj)

area.raster <- area(testN2)

# 2. Calcular el tamanno de ocurrencia de la espcie km2
OccExtension <- EstimateRangeSize( r = testN2, area.raster = area.raster)

# 3. Minimo poligono convexo para Raster
MCPRast <- PredictMcpRast( r = testN2, area.raster = area.raster, rast = rast)

# 4. Minimo poligono convexo para Registros de especies
MCPRec <- PredictMcpRec(r = testN2, Proj_col = Proj_col, proj = proj, rast = rast, s = s,  sppPath = sppPath)

# 5. Habitat CLC
Habitat <- HabitatAreaClc(r = testN2, area.raster = area.raster, clcBrick = clcBrick, clcKey = clcKey, rast = rast)
HabitatPercent <- Habitat*100/sum(Habitat) 

# 6. Historico perdida de bosque
ForestArea <- ForestLoss( r = testN2, area.raster = area.raster, f90 = f90, f00 = f00, f05 = f05, f10 = f10, f12 = f12, rast = TRUE) 
ForestPercentage <- (ForestArea / OccExtension * 100)

# 7. Escenarios 2030
EscenariosBosque <- Escenarios2030 (r = testN2, area.raster = area.raster, con2030_res = con2030_res, des2030_res = des2030_res, his2030_res = his2030_res, rast = TRUE, sppFolder = sppFolder)
EscenariosBosquePorcentaje <- (EscenariosBosque / OccExtension * 100)

# 8. Areas protegidas
ProtAreaRep <- AllProtectedFigures( r = testN2, area.raster = area.raster, ap = ap, pn = pn, sc = sc, ot = ot, rast = rast)

# Amenazas 
# 9.Human foot print
HFPrint <- HumanFootPrint (r = testN2, hfp = hfp, rast = rast)

# 10. Vias/Roads
RdsLength <- RoadsLength (r = testN2, rds = rds, rast = rast)

# 11. Titulos mineros
MinTitArea <- TitulosMineros (r = testN2, area.raster = area.raster, ntif = ntif, rast = rast)

sppModel <- gsub('.tif', '', basename(as.character(sppPath[s, 2]))) # Conserva el nombre de la especie con información del modelo, sí es concenso o nivel 2

#Out data frames
metadata <- cbind.data.frame(sppModel, HabitatPercent, ForestPercentage, EscenariosBosquePorcentaje, MCPRast$ch.area, MCPRec$MCPrec_km2, OccExtension, ProtAreaRep)
threats <- cbind.data.frame(sppModel, HFPrint$size, HFPrint$sumVal, HFPrint$meanVal, RdsLength, MinTitArea)
metadata.all <- rbind(metadata.all, metadata)
threats.all <- rbind(threats.all, threats)

}

colnames(metadata.all) <- colnames(header.met)
write.csv(metadata.all, 'C:/Biomodelos/Grilla/Registros/Zamias_metadata.csv', row.names = T)
write.csv(threats.all, 'C:/Biomodelos/Grilla/Registros/Zamias_threats.csv', row.names = T)




