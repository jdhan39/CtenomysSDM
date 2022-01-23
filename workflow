#Bash Work Flow
#Species of interest: Ctenomys opimus

#Pull raster files from:
#https://opentopography.s3.sdsc.edu/minio/dataspace/OTDS.012020.4326.1/raster/slope/
#Slope raster:
#https://opentopography.s3.sdsc.edu/dataspace/OTDS.012020.4326.1/raster/slope/slope_90M_s60w180.tar.gz
#Geographic region of interest: according to the MERIT tile system
#s30w090
#Citation for species polygons: 

cd /media/sf_LVM_shared/project/GBIF

#Pulling species data from GBIF
wget https://api.gbif.org/v1/occurrence/download/request/0249938-200613084148143.zip

#Selecting the polygons and creating a new .shp file
cd /media/sf_LVM_shared/project/species_data
ogr2ogr -sql "SELECT * FROM Ctenomys_distributions WHERE BINOMIAL = 'Ctenomys opimus'" \
opimus_poly.shp Ctenomys_distributions.shp

#Downloading, organizing, and cropping topographic raster files
cd /media/sf_LVM_shared/project

ulx=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\("," "); print int($2 -1)}')
uly=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\)"," "); print int($6 +1)}')
lrx=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\("," "); print int($5 +1)}')
lry=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\)"," "); print int($3 -1)}')

baseurl=https://opentopography.s3.sdsc.edu/dataspace/OTDS.012020.4326.1/raster
var_list=(pcurv elev-stdev roughness)

gridcode_list=(s30w090)

for var in ${var_list[@]}; do
	mkdir $var #mkdir for topographic element
	cd $var

	for gridcode in ${gridcode_list[@]}; do
		mkdir $gridcode
		cd $gridcode
		filename=${var}_90M_${gridcode}.tar.gz
		rm -f ${filename} # make sure this only removes the tar ball
		wget ${baseurl}/${var}/${filename}
		tar -xzf $filename
		cd ..
	done
	
	gdalbuildvrt -overwrite ${var}_90M.vrt ./*/*.tif
	gdal_translate -projwin $ulx $uly $lrx $lry -co COMPRESS=DEFLATE -co ZLEVEL=9 \
	${var}_90M.vrt ${var}_90M_${gridcode}crop.tif
	cd ..

done

#Downloading climate data

ulx=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\("," "); print int($2 -1)}')
uly=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\)"," "); print int($6 +1)}')
lrx=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\("," "); print int($5 +1)}')
lry=$(ogrinfo -al -so species_data/ctenomys_group.shp ctenomys_group | grep Extent | awk \
'{ gsub("\\)"," "); print int($3 -1)}')

mkdir climate_data
cd /media/sf_LVM_shared/project/climate_data

mkdir tmax tmin prec

# download data from https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V1
for MM in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
	#download tmin data
	cd tmin
	gdal_translate -projwin $ulx $uly $lrx $lry -co COMPRESS=DEFLATE -co ZLEVEL=9 \
	/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/\
	tmin/CHELSA_tmin10_${MM}_1979-2013_V1.2_land.tif CHELSA_tmin10_${MM}_1979-2013_V1.2_land_crop.tif
	cd ..
	#download precpitation data
	cd prec
	gdal_translate -projwin $ulx $uly $lrx $lry -co COMPRESS=DEFLATE -co ZLEVEL=9 \
	/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/\
	prec/CHELSA_prec_${MM}_V1.2_land.tif CHELSA_prec_${MM}_V1.2_land_crop.tif
	cd ..
	#download tmax data
	cd tmax
	gdal_translate -projwin $ulx $uly $lrx $lry -co COMPRESS=DEFLATE -co ZLEVEL=9 \
	/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/\
	tmax/CHELSA_tmax10_${MM}_1979-2013_V1.2_land.tif CHELSA_tmax10_${MM}_1979-2013_V1.2_land_crop.tif
	cd ..
done

#Calculating Mean and Standard Deviation for Climate Data
cd /media/sf_LVM_shared/project
# calculate mean and stdev for tmin10,

gdalbuildvrt -overwrite -separate CHELSA_tmin10_1979-2013_V1.2_land_crop.vrt \
climate_data/tmin/CHELSA_tmin10_*_1979-2013_V1.2_land_crop.tif \

pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32768 -f --mean -i \
CHELSA_tmin10_1979-2013_V1.2_land_crop.vrt -o climate_data/tmin/\
CHELSA_tmin10_1979-2013_V1.2_land_crop_mean.tif
pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32768 -f --stdev -i \
CHELSA_tmin10_1979-2013_V1.2_land_crop.vrt -o climate_data/tmin/\
CHELSA_tmin10_1979-2013_V1.2_land_crop_stdev.tif

gdal_edit.py -a_nodata -32768 climate_data/tmin/CHELSA_tmin10_1979-2013_V1.2_land_crop_mean.tif
gdal_edit.py -a_nodata -32768 climate_data/tmin/CHELSA_tmin10_1979-2013_V1.2_land_crop_stdev.tif

# calculate mean and stdev for tmax,

gdalbuildvrt -overwrite -separate CHELSA_tmax10_1979-2013_V1.2_land_crop.vrt \
climate_data/tmax/CHELSA_tmax10_*_1979-2013_V1.2_land_crop.tif

pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32768 -f --mean -i \
CHELSA_tmax10_1979-2013_V1.2_land_crop.vrt -o climate_data/tmax/\
CHELSA_tmax10_1979-2013_V1.2_land_crop_mean.tif
pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32768 -f --stdev -i \
CHELSA_tmax10_1979-2013_V1.2_land_crop.vrt -o climate_data/tmax/\
CHELSA_tmax10_1979-2013_V1.2_land_crop_stdev.tif

gdal_edit.py -a_nodata -32768 climate_data/tmax/CHELSA_tmax10_1979-2013_V1.2_land_crop_mean.tif
gdal_edit.py -a_nodata -32768 climate_data/tmax/CHELSA_tmax10_1979-2013_V1.2_land_crop_stdev.tif

# calculate mean and stdev for prec,

gdalbuildvrt -overwrite -separate CHELSA_prec_V1.2_land_crop.vrt \
climate_data/prec/CHELSA_prec_*_V1.2_land_crop.tif

pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32767 -f --mean -i \
CHELSA_prec_V1.2_land_crop.vrt -o climate_data/prec/CHELSA_prec_V1.2_land_crop_mean.tif
pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -nodata -32767 -f --stdev -i \
CHELSA_prec_V1.2_land_crop.vrt -o climate_data/prec/CHELSA_prec_V1.2_land_crop_stdev.tif

gdal_edit.py -a_nodata -32767 climate_data/prec/CHELSA_prec_V1.2_land_crop_mean.tif
gdal_edit.py -a_nodata -32767 climate_data/prec/CHELSA_prec_V1.2_land_crop_stdev.tif

#-------------------------------------------------------------------------------------
#Modeling in R

library(ggplot2)
library(rgdal)
library(raster)
library(rgeos)
library(reshape)
library(rasterVis)
library(dismo)
library(InformationValue)
library(mgcv)
library(randomForest)
library(e1071)
library(caret)
library(maptools)
library(tidyr)

set.seed(4532)

setwd('/media/sf_LVM_shared/project/GBIF')
gbif.occ <- read.delim(file = "gbif_genus_occ.csv")
head(gbif.occ, n=2)

# Working df of only opimus with all N/A coordinates removed.
sp.raw <- gbif.occ[gbif.occ$species == 'Ctenomys opimus',]
sp.coor <- sp.raw[,c('decimalLatitude', 'decimalLongitude')]
sp.df.org <- sp.coor[!is.na(sp.coor$decimalLatitude),]
head(sp.df.org)

#Reordering to put longtitude first
sp.df <- sp.df.org[,c(2,1)]
colnames(sp.df) <- c("lon","lat")
head(sp.df)
nrow(sp.df)

#Enviornmental data
setwd('/media/sf_LVM_shared/project')

#Topographic Data
elev <- raster('./elev-stdev/elev-stdev_90M_s30w090crop.tif')
pcurv <- raster('./pcurv/pcurv_90M_s30w090crop.tif')
ruff<- raster('./roughness/roughness_90M_s30w090crop.tif')

#Climate Data: mean only
tmin <- raster('./climate_data/tmin/CHELSA_tmin10_1979-2013_V1.2_land_crop_mean.tif')
tmax<- raster('./climate_data/tmax/CHELSA_tmax10_1979-2013_V1.2_land_crop_mean.tif')
prec <- raster('./climate_data/prec/CHELSA_prec_V1.2_land_crop_mean.tif')

#Adding polygons
sp.poly <-readOGR("/media/sf_LVM_shared/project/species_data", "opimus_poly")

plot(tmin)
points(sp.df$lon, sp.df$lat, col = "red", cex = .3)
plot(sp.poly, add=TRUE)

# Stacking and resampling
topo <- stack(c(ruff, pcurv, elev))
climate <- stack(c(tmin,tmax,prec))
r.topo <- resample(topo, climate, method='bilinear')

# Making complete enviornmental stack
layers <- stack(r.topo, climate)

# Renaming layers
names(layers) <- c("Roughness", "Pcurv", "Elev", "Mean.MIN.Temp", "Mean.MAX.Temp", "Precipitation")
names(layers)

# Vizualizing Stack layers
options(repr.plot.width=18, repr.plot.height=11)
plot(layers,cex.axis=2,cex.lab=2,cex.main=2.5,legend.args=list(text=NULL, side=1, cex.lab = 3, line=2.3))

#Identifying Correlations in data
options(repr.plot.width=18, repr.plot.height=16)
pairs(layers,maxpixels=1000,cex=0.5)

#Dropping elev
layers <- dropLayer(layers, 3)
names(layers)

#Generating Absence Data

# scaling layers
predictors <-c("Roughness",'Pcurv', "Mean.MIN.Temp", "Mean.MAX.Temp", "Precipitation")
sc.layers <- scale(layers[[predictors]])

#Generate background points
back.xy <- randomPoints(sc.layers, p=sp.df, n=230)
colnames(back.xy) <- c("lon","lat")

options(repr.plot.width=9, repr.plot.height=9)
plot(ruff)
points(back.xy)

#extract GIS data
pres.idpv <- raster::extract(sc.layers, sp.df)
back.idpv <- raster::extract(sc.layers, back.xy)

#link data
df.pres <- data.frame(sp.df, pres.idpv, pres=1)
df.back <- data.frame(back.xy, back.idpv, pres=0)

#remove any potential NAs
df.pres <- df.pres[complete.cases(df.pres),]
df.back <- df.back[complete.cases(df.back),]

# merge together
df.all <- rbind(df.pres, df.back)
head(df.all)

#Scaling Data
df.all[,predictors] <- scale(df.all[,predictors])

head(df.all, n=2)

#Pivot
df.all.piv <- tidyr::pivot_longer(df.all, -c(lon, lat, pres), names_to="predictors")
# Ordering by perdictor and coordinates
df.all.po <- df.all.piv[with(df.all.piv, order(predictors)),]
head(df.all.po)

ggplot(df.all.po,aes(x=value,y=pres))+facet_wrap(~predictors)+
	geom_point()+
	stat_smooth(method = "lm", formula = y ~ x + I(x^2), col="red")+
	geom_smooth(method="gam",formula=y ~ s(x, bs = "cs")) +
	theme(text = element_text(size = 20))

#Model Fitting: GLM

summary(df.all)

# Training-Test split
df.all$grp <- kfold(df.all,5)
#gives count of unique values
table(df.all$grp)

mdl.glm <- glm(pres~Roughness+Pcurv*I(Pcurv^2)+Mean.MIN.Temp+Mean.MAX.Temp+Precipitation, \
family=binomial(link=logit), data=subset(df.all,grp!=1))

summary(mdl.glm)

#ROC curve 
pred.glm1 <- predict(mdl.glm,df.all[which(df.all$grp!=1),predictors],type="response")
pred.glm2 <- predict(mdl.glm,df.all[which(df.all$grp==1),predictors],type="response")

plotROC(df.all[which(df.all$grp!=1),"pres"],pred.glm1)

plotROC(df.all[which(df.all$grp==1),"pres"],pred.glm2)

#GLM OutMapping
p1 <- raster::predict(sc.layers, mdl.glm,type="response")

options(repr.plot.width=8, repr.plot.height=9)
gplot(p1)+geom_tile(aes(fill=value))+
	scale_fill_gradientn(
		colours= c("dodgerblue4","green","yellow","orange","firebrick"),
		na.value = "transparent")+
	geom_polygon(aes(x=long,y=lat,group=group),
							data=fortify(sp.poly),fill="transparent",col="black", size = 1.3)+
geom_point(aes(x = lon, y = lat), data = subset(df.all,pres==1),shape=21, fill='snow1', col="black",size=3)+
coord_equal()

#Cross Validation
ctrl <- trainControl(method = "cv", number = 10)
	mdl.glm.cv <- train(as.factor(pres)~Roughness+Pcurv+Mean.MIN.Temp+Mean.MAX.Temp+Precipitation, \
	family=binomial(link=logit), data = df.all, method = "glm", trControl = ctrl, \
	metric='Accuracy')
	
	mdl.glm.cv
	
#Predictive mapping
pred.glm.resp <- raster::predict(sc.layers, mdl.glm.cv, type='raw')
plot(pred.glm.resp)

pred.mdl.glm.prob1 <- raster::predict(sc.layers,mdl.glm.cv,type="prob" , index=1)
plot(pred.mdl.glm.prob1)

pred.mdl.glm.prob2 <- raster::predict(sc.layers,mdl.glm.cv ,type="prob" , index=2)
plot(pred.mdl.glm.prob2)

