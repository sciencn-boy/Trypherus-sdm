###---------------SDM------------------###
#--------------download packages------------------#
#install.packages(c("dismo", "terra", "ggplot2", "data.table", "raster", "rJava", "psych", "ENMeval",'sdm','geodata','usdm','data.table'))
library(dismo)
library(dplyr)
library(tidyr)
library(mapview)
library(raster)
library(geodata)
library(sdm)
library(usdm)
library(data.table)
setwd("E:/Trypherus")  #Setting up the working directory
spg <- read.csv("Trypherus.csv")
sp <- read.csv("Trypherus.csv")
spg <- spg[,-1]
head(spg)
colnames(spg)[1] <- "lon"
colnames(spg)[2] <- "lat"
spg$species <- 1
head(spg)
class(spg)
spg <- spg%>%drop_na()
nrow(spg)
coordinates(spg) <- c('lon','lat')
class(spg)

#-----------------read data---------------------#
vars <- c(1:21)      #BIO20：elevation
bio <-stack(sprintf("worldclim/Current/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d.tif", vars))
names(bio) <- paste0("bio",vars)

#------------exclude-----------#
ex <- raster::extract(bio,spg)
v <- vifstep(ex)
bioc <- exclude(bio,v)
bioc

# vars_bioc <- c(2,3,4,8,9,10,15,17,18,20) #maxent
vars_bioc <- c(3,7,8,9,10,14,15,18,20,21)    #vif
bioc<-stack(sprintf("worldclim/Current/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d.tif", vars_bioc))
names(bioc) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE","HII")
#future
bioc_126_50 <- stack(sprintf("worldclim/Future/ssp126_2050/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_126_50) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE","HII")
bioc_126_70 <- stack(sprintf("worldclim/Future/ssp126_2070/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_126_70) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE","HII")
bioc_585_50 <- stack(sprintf("worldclim/Future/ssp585_2050/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_585_50) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE","HII")
bioc_585_70 <- stack(sprintf("worldclim/Future/ssp585_2070/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_585_70) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE","HII")

#-------------model----------------#
#sdm，use'glm','gam','rf','maxent'；sub，boot，0.25，3
d <- sdmData(species~.,spg,predictors = bioc,bg=list(method='gRandom',n=1000))
m <- sdm(species~., d, methods = c('glm','gam','rf','maxent'),replication=c('sub','boot'),
         test.percent=25,n=10,parallelSetting=list(ncore=4,method='parallel'))
glm <- sdm(species~., d, methods = 'glm',replication=c('sub','boot'),
           test.percent=25,n=10,parallelSetting=list(ncore=4,method='parallel'))
gam <- sdm(species~., d, methods = 'gam',replication=c('sub','boot'),
           test.percent=25,n=10,parallelSetting=list(ncore=4,method='parallel'))
rf <- sdm(species~., d, methods = 'rf',replication=c('sub','boot'),
          test.percent=25,n=10,parallelSetting=list(ncore=4,method='parallel'))
maxent <- sdm(species~., d, methods ='maxent',replication=c('sub','boot'),
              test.percent=25,n=10,parallelSetting=list(ncore=4,method='parallel'))

# predict
p_glm <- predict(glm,bioc)
terra::writeRaster((p_glm), filename="glm.tif",overwrite=TRUE)
p_gam <- predict(gam,bioc)
terra::writeRaster((p_gam), filename="gam.tif",overwrite=TRUE)
p_rf <- predict(rf,bioc)
terra::writeRaster((p_rf), filename="rf.tif",overwrite=TRUE)
p_maxent <- predict(maxent,bioc)
terra::writeRaster((p_maxent), filename="maxent.tif",overwrite=TRUE)

#---------ensemble-------------------#
#en1 <- ensemble(m,vir,filename = "en.img",setting=list(method='weighted',stat='tss',opt=2))
en1 <- ensemble(m,bioc,setting=list(method='weighted',stat='tss',opt=2))  #ensemble
plot(en1)
terra::writeRaster((en1), filename="ensemble.tif",overwrite=TRUE)

#126_50
bioc_126_50 <- stack(sprintf("worldclim/Future/ssp126_2050/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_126_50) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE")
#
p_glm <- predict(glm,bioc_126_50)
terra::writeRaster((p_glm), filename="fu126_50_glm.tif",overwrite=TRUE)
p_gam <- predict(gam,bioc_126_50)
terra::writeRaster((p_gam), filename="fu126_50_gam.tif",overwrite=TRUE)
p_rf <- predict(rf,bioc_126_50)
terra::writeRaster((p_rf), filename="fu126_50_rf.tif",overwrite=TRUE)
p_maxent <- predict(maxent,bioc_126_50)
terra::writeRaster((p_maxent), filename="fu126_50_maxent.tif",overwrite=TRUE)

en_fu126_50 <- ensemble(m,bioc_126_50,setting=list(method='weighted',stat='tss',opt=2))
plot(en_fu126_50)
terra::writeRaster((en_fu126_50), filename="ensemble_fu126_50.tif",overwrite=TRUE)

#126_70

bioc_126_70 <- stack(sprintf("worldclim/Future/ssp126_2070/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_126_70) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE")
# 
p_glm <- predict(glm,bioc_126_70)
terra::writeRaster((p_glm), filename="fu126_70_glm.tif",overwrite=TRUE)
p_gam <- predict(gam,bioc_126_70)
terra::writeRaster((p_gam), filename="fu126_70_gam.tif",overwrite=TRUE)
p_rf <- predict(rf,bioc_126_70)
terra::writeRaster((p_rf), filename="fu126_70_rf.tif",overwrite=TRUE)
p_maxent <- predict(maxent,bioc_126_70)
terra::writeRaster((p_maxent), filename="fu126_70_maxent.tif",overwrite=TRUE)

en_fu126_70 <- ensemble(m,bioc_126_70,setting=list(method='weighted',stat='tss',opt=2))
plot(en_fu126_70)
terra::writeRaster((en_fu126_70), filename="ensemble_fu126_70.tif",overwrite=TRUE)

#585_50

bioc_585_50 <- stack(sprintf("worldclim/Future/ssp585_2050/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_585_50) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE")
# 
p_glm <- predict(glm,bioc_585_50)
terra::writeRaster((p_glm), filename="fu585_50_glm.tif",overwrite=TRUE)
p_gam <- predict(gam,bioc_585_50)
terra::writeRaster((p_gam), filename="fu585_50_gam.tif",overwrite=TRUE)
p_rf <- predict(rf,bioc_585_50)
terra::writeRaster((p_rf), filename="fu585_50_rf.tif",overwrite=TRUE)
p_maxent <- predict(maxent,bioc_585_50)
terra::writeRaster((p_maxent), filename="fu585_50_maxent.tif",overwrite=TRUE)

en_fu585_50 <- ensemble(m,bioc_585_50,setting=list(method='weighted',stat='tss',opt=2))
plot(en_fu585_50)
terra::writeRaster((en_fu585_50), filename="ensemble_fu585_50.tif",overwrite=TRUE)

#585_70

bioc_585_70 <- stack(sprintf("worldclim/Future/ssp585_2070/wc2.1_2.5m_60/wc2.1_2.5m_bio_%d_fu.tif",vars_bioc))
names(bioc_585_70) <- c("BIO3","BIO7","BIO8","BIO9","BIO10","BIO14","BIO15","BIO18","ELE")
# 
p_glm <- predict(glm,bioc_585_70)
terra::writeRaster((p_glm), filename="fu585_70_glm.tif",overwrite=TRUE)
p_gam <- predict(gam,bioc_585_70)
terra::writeRaster((p_gam), filename="fu585_70_gam.tif",overwrite=TRUE)
p_rf <- predict(rf,bioc_585_70)
terra::writeRaster((p_rf), filename="fu585_70_rf.tif",overwrite=TRUE)
p_maxent <- predict(maxent,bioc_585_70)
terra::writeRaster((p_maxent), filename="fu585_70_maxent.tif",overwrite=TRUE)

en_fu585_70 <- ensemble(m,bioc_585_70,setting=list(method='weighted',stat='tss',opt=2))
plot(en_fu585_70)
terra::writeRaster((en_fu585_70), filename="ensemble_fu585_70.tif",overwrite=TRUE)


#-------------#
ch <- en2 - en1
cl2 <- colorRampPalette(c('red','orange','white','gray','green','blue'))
plot(ch,col=cl2(200))

#-----------------------
df <- as.data.frame(d)    
head(df)
#d1 <- as.data.frame(d)
df <- data.frame(species=df$species,coordinates(d@info@coords))  
head(df)
xy <- as.matrix(df[,c('lon','lat')])
head(xy)
p <- raster::extract(en1,xy)    
head(p)
nrow(df)
length(p)
ev <- evaluates(df$species,p)
ev@statistics
ev@threshold_based   

th <- ev@threshold_based$threshold[2]  

pa1 <- raster(en1)

pa1[] <- ifelse(en1[] >= th, 1 , 0)
plot(pa1)

pa2 <- raster(en2)

pa2[] <- ifelse(en2[] >= th, 1 , 0)
plot(pa2)

chp <- pa2-pa1
plot(chp,col=c('red','gray','blue'))  

#-----------------------
rcurve(m,id=7:12)  
getVarImp(m,id=1)  
plot(getVarImp(m))