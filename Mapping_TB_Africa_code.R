#here you need to keep only the countries where you have data and remove the other ones
#because if not the model will predict strange values in Sahara etc.!

#simple R-INLA model as template
#1.step install the packages you need

###############user specifications#########################################
#define number of samples 
nn = 10000#(10000 for final results, 150 for fast results)
aggfactor <- 5#factor of aggregation for final maps
popt <- 5 #threshold for population to mask areas
# unit: persons per square kilometer.
mu <-0.025#we do not expect more than 25 per 1,000 prev per cell, to be discussed with experts
# define the unit for prevalence map
prevunit <- 1000#if 1000: values are in per 1,000; if 100, values are in %.
#note that only prevalence maps will be in per 1,000, counts will remain untouched
#Run various model alternatives
#with population filter
popfilter <- FALSE
#removing Mozambique (outlier)
mozout <- FALSE
#if allrunis TRUE: run all models (with/out popfilter and with/out Mozambique)
allrun <- FALSE
#if allrun is FALSE it will run one model based on mozout and popilter choices
###############end user specifications#########################################

#INLA used to fit Bayesian models
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer", "ggmap", "rgdal", "rgeos","maptools", "tmap","gtools","fmsb")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
library(INLA)

#2. step: define the paths
path_input <- paste0(getwd(),"/INPUT")


#loop
if(allrun==TRUE){
  combinations <-expand.grid(popfilter=c(TRUE,FALSE),mozout=c(TRUE,FALSE))
}else{
  combinations <-data.frame(popfilter=popfilter,mozout=mozout)
}
#intial data loading


#study area
Africa <- rnaturalearth::ne_countries()
Africa <- Africa[Africa$continent=="Africa",]
#plot(Africa, main="Adm. Boundaries Africa Level 0")

#TB data 
TB0 <- read_excel(paste0(path_input,"/TB_Africa_v2.xlsx"),trim_ws = TRUE)
TB0$Number_examined <- as.numeric(TB0$Number_examined)
TB0$`Total TB cases` <- as.numeric(TB0$`Total TB cases`)
TB0$Latitude <- as.numeric(TB0$Latitude)
TB0$Longitude <- as.numeric(TB0$Longitude)

#############OPTION TO BE DISCUSSED##################
#keep all data except if country-level
TB0 <- TB0[TB0$`Lowest admin level (0-4)`>0,]
#############OPTION TO BE DISCUSSED##################
#remove potential mistakes
TB0$Number_examined <- ifelse(TB0$Number_examined<0,0,TB0$Number_examined)# non neg response
TB0$`Total TB cases` <- round(TB0$`Total TB cases`,0)# count response
TB0$Number_examined <- round(TB0$Number_examined,0)# count denominator


#hist(TB$`Lowest admin level (0-4)`)
#remove potential NA in country 
TB0 <-TB0[complete.cases(TB0[,c("ISO3")]),]
#*************Temp from worldclim*******************************************************************************************
T0 <- raster::getData("worldclim",var="bio",res=5) # note that in WorldClim bio1 is mean annual Temp and bio 12 is precpt
T0 <- T0[[c(1)]]
#*************Prec from worldclim*****************************************************************************************
P0 <- raster::getData("worldclim",var="bio",res=5) # note that in WorldClim bio12 is the annual mean precpt
P0 <- P0[[c(12)]]
#**************altitude*******************************************************************************************************
alt0<-geodata::elevation_global(res=5,path=path_input)
alt0<-raster::raster(alt0)
#*************travel time to health facility from MAP*********************************************************************************************
ahf0 <- raster(paste0(path_input,"/2020_walking_only_travel_time_to_healthcare.tif"))
#**********************population density**********************************************************************************************
popden0 <-raster(paste0(path_input,"/gpw-v4-population-density_2000.tif"))

##########################################################################
#loop across all model specifications#####################################
for (k in 1:nrow(combinations))
{
mycomb <- combinations[k,]
subfiles <- ifelse(mycomb$popfilter & mycomb$mozout, "NOMOZ/FILTER",
                        ifelse(mycomb$popfilter & !mycomb$mozout, "ALL/FILTER",
                               ifelse(!mycomb$popfilter & mycomb$mozout, "NOMOZ/NOFILTER",
                                      "ALL/NOFILTER")))
path_output <- file.path(getwd(), "OUTPUT", subfiles)
#3. remove Mozambique or not
if (mycomb$mozout==TRUE){
  TB <- subset(TB0,ISO3!= "MOZ")
} else {TB <-TB0}


#write.csv(TB,"TB.csv")
#select countries in the study area based on worldbank ISO 3
#identify countries to be removed 
outc <- setdiff(unique(Africa$wb_a3), unique(TB$ISO3))#vector 1,2 as arguments
outc <- outc[complete.cases(outc)]
#remove iteratively countries not in TB data
myarea <- Africa
for (i in 1:length(outc)){
myarea = subset(myarea,wb_a3 != outc[i])
}
#plot(Africa);plot(myarea,add=TRUE,col="red")
#check if the right countries have been selected
if(isFALSE(unique(sort(myarea$wb_a3))==unique(sort(TB$ISO3))))
  stop("Error: countries of TB do not match study area!")
## crop and mask
T1 <- mask(crop(T0, extent(myarea)),myarea)
## Check that it worked
#plot(T1, main="Temperature");plot(myarea, add=TRUE, lwd=1)
## crop and mask
P1 <- mask(crop(P0, extent(myarea)),myarea)
## Check that it worked
#plot(P1, main="Precipitation")
alt <- mask(crop(alt0, extent(myarea)),myarea)
#plot(alt, main="Altitude")
# problem to access it
# acc <- malariaAtlas::getRaster(surface="A global map of travel time to cities to assess inequalities in accessibility in 2015",shp=Africa)
# raster::raster(acc,file=paste0(path_input,"/accAfrica.tif"))
# acc <- raster::writeRaster(acc,paste0(path_input,"/accAfrica.tif"))
acc <- raster::raster(paste0(path_input,"/accAfrica.tif"))

# plot(acc)
ahf <- mask(crop(ahf0, extent(myarea)),myarea)
#plot(ahf)
popden <- mask(crop(popden0, extent(myarea)),myarea)
#plot(popden)
#*******stack rasters using stack(raster1,raster2,...)************************************************
# as the variables (T3,P3,alt,Access,popden) are from different sources and dimensiton*so it needs to b resampleed before stack************
#T3 <- resample(T3,alt)
acc <- resample(acc,alt)
popden <- resample(popden,alt)
#P3 <- resample(P3,alt, method='ngb')
ahf <- resample(ahf,alt)
# rs <- stack(T1,P1,alt,ahf,popden)
# names(rs) <- c("TMP","PCP","ALT","ACH","POP")
rs <- stack(T1,P1,alt,ahf,acc,popden)
names(rs) <- c("TMP","PCP","ALT","ACH","ACC","POP")

#*************x-mean/sd to standardize the covariates rs ****************************************
#rs2 <-scale(rs)
#quantile normalization (rank values and make correspond to normally distributed data)
norm<-list()
new_var<-list()
new_var_full<-list()
st2<-list()
x<-list()
n<-list()
for (i in 1:nlayers(rs)){
  # linear is my raster
  st2[[i]] <- rs[[i]]
  x[[i]] <- getValues(rs[[i]])
  n[[i]] <- length(na.omit(x[[i]]))
  norm[[i]] <- qnorm(seq(0.0, 1, length.out = n[[i]] + 2)[2:(n[[i]] + 1)])
  new_var[[i]] <- norm[[i]][base::rank(na.omit(x[[i]]))]
  new_var_full[[i]] <- rep(NA, length(x[[i]]))
  new_var_full[[i]][!is.na(x[[i]])] <- new_var[[i]]
  values(st2[[i]]) <- new_var_full[[i]]
}
rs2<-stack(st2)#

#mapping covariates
rsmap <- rs
#names(rsmap) <- c("TMP","PCP","ALT","ACH","logPOP")
#names(rsmap) <- c("TMP","PCP","ALT","ACH","ACC","POP")

#log pop for mapping purposes
rsmap[["POP"]] <- log(subset(rsmap, 'POP')+1)
names(rsmap[["POP"]]) <- "logPOP"

af.layer <- function() {
  plot(Africa, add=TRUE)}
plot.new()
pdf(file=paste0(path_output,"/pdf/Covariates.pdf"))
plot(rsmap,col=viridis(10),nc=2,addfun=af.layer)
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#save data preparation output here
save(rs2,rs,myarea,TB,popden, path_input,path_output,Africa,file=paste0(path_output,"/dataprep.Rdata"))
#load(file=paste0(path_output,"/dataprep.Rdata"))

#********************************************************
pt <-TB
pt <- data.frame(TB=pt$`Total TB cases`,Nscreen = pt$Number_examined, x=pt$Longitude, y=pt$Latitude)
pt <- pt[complete.cases(pt),]

#hist(pt$TB)

#**********combine the coordinate**********************
xy <- cbind(pt$x, pt$y)
#extract covariate at point coordinates
covariate_all <-data.frame(raster::extract(rs2, xy))
covariate_z <- data.frame(covariate_all)
#apply VIF to subset data
source("INPUT/vif.R")
#The function uses three arguments. The first is a matrix or data frame of the explanatory variables,
#the second is the threshold value to use for retaining variables, and the third is a logical argument 
#indicating if text output is returned as the stepwise selection progresses. The output indicates the VIF
#values for each variable after each stepwise comparison. The function calculates the VIF values for all 
#explanatory variables, removes the variable with the highest value, and repeats until all VIF values are
#below the threshold. The final output is a list of variable names with VIF values that fall below the threshold. 
keep.dat <-vif_func(in_frame=covariate_z,thresh=5,trace=F)
#temperature and altitude are also highly correlated (so remove altitude)
rsmap <- rsmap[[keep.dat]]
rs <- rs[[keep.dat]]
rs2 <- rs2[[keep.dat]]
covariate_z <- covariate_z[paste0(keep.dat)]

#correlation
# jnk=raster::layerStats(covariate_z, 'pearson', na.rm=TRUE)
# corr_matrix=jnk$'pearson correlation coefficient'
corr_matrix <- cor(covariate_z,use="complete.obs")
write.csv(corr_matrix,paste0(path_output,"/csv/covariate_correlation.csv"))


#pop and acc highly correlated so we remove pop
#ach and acc highly correlated so we remove ach
covariate_z$POP <- NULL
covariate_z$ACH <- NULL

#*********combine all data together**
ptcov<-cbind(pt,covariate_z)
#summary(ptcov)
#######MODELLING PART STARTS HERE#####################################################################
######################################################################################################
#define mesh object
bdry <- INLA::inla.sp2segment(myarea)
bdry$loc <- INLA::inla.mesh.map(bdry$loc)
mesh<-INLA::inla.mesh.2d(loc=xy, boundary=bdry, max.edge=c(0.6,4),offset=c(0.6,4),cutoff = 0.6)#to be fine tuned
mesh$n
#save mesh plot
plot.new()
pdf(file=paste0(path_output,"/pdf/mesh.pdf"))
plot(Africa);plot(mesh, add=TRUE); lines(bdry, col=1);plot(Africa,add=TRUE,border="darkgreen",lwd=2.5);
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#define RINLA objects
A = inla.spde.make.A(mesh=mesh, loc=as.matrix(xy));dim(A)   #A matrix
spde = inla.spde2.pcmatern(mesh, prior.range=c(10,0.1),prior.sigma = c(1,0.1))#basic spde object with default priors
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
Y=ptcov$TB
maxprev <- max(Y/ptcov$Nscreen)
#stack
stk = inla.stack(data=list(Y=ptcov$TB, n=ptcov$Nscreen),A=list(A,1),effects=list(
  list(z.field=1:spde$n.spde),list(#z.intercept=rep(1,length(Y)),
                                   covariate=covariate_z)),tag="est.z")
covn <- ncol(covariate_z)

###model selection (other criteria can be used)

mypb <- txtProgressBar(min = 0, max = covn, initial = 0, width = 150, style = 3) 
for(i in 1:covn){
 # f1 <- as.formula(paste0("Y ~ -1 + z.intercept + f(z.field, model=spde) + ", paste0(colnames(covariate_z)[1:i], collapse = " + ")))
  f1 <- as.formula(paste0("Y ~  f(z.field, model=spde) + ", paste0(colnames(covariate_z)[1:i], collapse = " + ")))
  model1<-inla(f1, #the formula
               data=inla.stack.data(stk,spde=spde),  #the data stack
               family= 'binomial',   #which family the data comes from
               Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
               control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
               #control.compute = list(mlik = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               control.compute = list(waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               control.fixed = list(prec=1000,prec.intercept=.001),
               verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
  # model_selection <- if(i==1){rbind(c(model = paste(colnames(covariate_z)[1:i]),mlik = model1$mlik[1]))
  #   }else{rbind(model_selection,c(model = paste(colnames(covariate_z)[1:i],collapse = " + "),mlik=model1$mlik[1]))
  model_selection <- if(i==1){rbind(c(model = paste(colnames(covariate_z)[1:i]),waic = model1$waic$waic))
  }else{rbind(model_selection,c(model = paste(colnames(covariate_z)[1:i],collapse = " + "),waic=model1$waic$waic))
        
       }
  setTxtProgressBar(mypb, i, title = "Model fit completed", label = i)
}
model_selection <- data.frame(model_selection)#provides the predictive performance for each investigated model
model_selection[,2] <- as.numeric(model_selection[,2])
write.csv(model_selection,paste0(path_output,"/csv/predperf_models.csv"))
#select the best model based on marginal likelihood
modsel <- model_selection[which.min(model_selection[,2]),]

#6. Run the selected model based on the results of the model selection
formula.z = as.formula(paste("Y ~  f(z.field, model=spde) +",modsel$model))

mymodel<-inla(formula.z, #the formula
            data=inla.stack.data(stk,spde=spde),  #the data stack
            family= 'binomial',   #which family the data comes from
            Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
            control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
           control.fixed = list(prec=1000,prec.intercept=.001),
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
#***********************************
#summary(mymodel)
#improve cpo computation (optional, takes about 10-15mn for 143 cases)
if(mymodel$ok==FALSE){
mymodel = inla.cpo(mymodel, force=TRUE)
}
#save coefficient correlation 
modeloutput <- data.frame(round(mymodel$summary.fixed,3))
write.csv(modeloutput,paste0(path_output,"/csv/finalmodeloutput.csv"))
#plot
plot.new()
pdf(file=paste0(path_output,"/pdf/validitycheck1.pdf"),width=8,height=6,onefile = TRUE)
plot(mymodel)
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#one useful diagnostic plot for PIT
#The PIT is the probability of a new response less than the observed response using a model based on the rest of the data. 
#We'd expect the PIT values to be uniformly distributed if the model assumptions are correct.
#check http://julianfaraway.github.io/brinla/examples/chicago.html for info
pit <- mymodel$cpo$pit
n <- length(Y)
uniquant <- (1:n)/(n+1)
plot(uniquant, sort(pit), xlab="uniform quantiles", ylab="Sorted PIT values")
# plot.new()
# pdf(file=paste0(path_output,"/pdf/validitycheck2.pdf"),width=8,height=6,onefile = TRUE)
# plot(uniquant, sort(pit), xlab="uniform quantiles", ylab="Sorted PIT values")
# abline(0,1)
# if(!is.null(dev.list())) dev.off()

#Report the plot with the logit transform
library(gtools)
plot.new()
pdf(file=paste0(path_output,"/pdf/validitycheck2.pdf"),width=8,height=6,onefile = TRUE)
plot(logit(uniquant), logit(sort(pit)), xlab="uniform quantiles (logit scale)", ylab="Sorted PIT values (logit scale)", main="")
abline(0,1)
if(!is.null(dev.list())) dev.off()
if(!is.null(dev.list())) dev.off()

#6. step: make predictions
#get predictive locations based on covariate (take anyone, we took the first one)
mask<-rs2[[1]]; NAvalue(mask)=-9999
mask<-aggregate(mask, fact=aggfactor)#for computational reasons we make predictions at coarser level
#create population mask as well for filtering areas where pop < threshold
popmask<-aggregate(popden, fact=aggfactor)#for computational reasons we make predictions at coarser level
popmask[popmask < popt] <- NA
popmask[!is.na(popmask),] <- 1
#plot(popmask)

#code to make maps
pred_val<-getValues(mask)
w<-is.na(pred_val)
index<-1:length(w)
index<-index[!w]
pred_locs<-xyFromCell(mask,1:ncell(mask))
pred_locs<-pred_locs[!w,]
colnames(pred_locs)<-c('longitude','latitude')
locs_pred <- pred_locs
#Mapping between meshes and continuous space
A.pred = inla.spde.make.A(mesh=mesh, loc=locs_pred)
mycovnames <- mymodel$names.fixed
mycovnames <- mycovnames[-1]
covariates_pred <- rs2[[mycovnames]]
covariates_pred = data.frame(extract(covariates_pred, pred_locs))
names(covariates_pred) <- mycovnames

#use link that account for max prevalence from the sample
#to avoid extreme values in the prediction of TB prev
linkfun <- function(x){
  mu/(1+exp(-x))
}
#sampling the posterior
set.seed(999)
samp = inla.posterior.sample(nn, mymodel)
pred = matrix(NA,nrow=dim(A.pred)[1], ncol=nn)
k = dim(covariates_pred)[2] ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('z.field',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('(Intercept)',rownames(samp[[i]]$latent)),]
  beta = NULL
  for (j in 1:k){
    beta[j] = samp[[i]]$latent[grep(names(covariates_pred)[j],rownames(samp[[i]]$latent)),]
  }
  #compute beta*covariate for each covariate
  linpred<-list()
  for (j in 1:k){
    linpred[[j]]<-beta[j]*covariates_pred[,j]
  }
  linpred<-Reduce("+",linpred)
  lp = intercept + linpred + drop(A.pred%*%field)
  ## Predicted values
   pred[,i] = linkfun(lp) #for transformation from 0 to mu
 # pred[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_2.5 <- apply(pred, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_med = apply(pred, 1, function(x) quantile(x, probs=c(0.5), na.rm=TRUE))
pred_sd = apply(pred, 1, sd)
pred_25pct = apply(pred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_75pct = apply(pred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
IQR = pred_75pct - pred_25pct
pred_mean <- apply(pred, 1, function(x) mean(x, na.rm = TRUE))
pred_975 <- apply(pred, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
#Saving predictive maps as raster files
predinput=list(pred_2.5,pred_med,pred_sd,pred_25pct,pred_75pct,IQR,pred_mean,pred_975)
prednames=as.list(c("LUI","prob_med","Prob_sd","Prob_q25","Prob_q75","Prob_IQR","Prob_mean","UUI"))
mainnames=as.list(c("PR 2.5th","PR med","PR sd","PR 25th pct","PR 75th pct","IQR","mean","PR 97.5th"))
#destfile <- paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- predinput[[j]]
  out[[j]] = setValues(mask, pred_val)
  #mask with population density
  if (mycomb$popfilter==TRUE){
  out[[j]] = out[[j]]* popmask * prevunit} else {
    out[[j]] = out[[j]]* prevunit} 
  writeRaster(out[[j]], paste0(path_output,'/',"prevalence/",prednames[[j]],'.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=paste0(path_output,'/',"prevalence"),pattern =".tif$", full.names=TRUE)
#extract raster data
rasterls<-list()
for (i in 1:length(raster.list)){
  rasterls[[i]]<-raster::raster(raster.list[[i]])
}
b <- raster::brick(rasterls)
#add informative text before the graphs
t1 = "Predicted median prevalence (per 1,000)"
t2 = "Predicted IQR (per 1,000)"
t3 = "Predicted Q25 (per 1,000)"
t4 = "Predicted Q75 (per 1,000)"
t5 = "Predicted sd (per 1,000)"
allt <- list(t1,t2,t3,t4,t5)

#scale bar
# Function to add a scalebar to a base-graphics plot
myScalebar = function(units_label, yadj=2) {
  # Get plot coordinates
  pc = par("usr") 
  # Position scale line between last two major x-axis tick marks
  # and 2/10th, (0.15) of the total y-range above the lower y-axis coordinate
  lines(c(floor(pc[2]-1),floor(pc[2])),     
        rep(pc[3] + 0.15*(pc[4] - pc[3]), 2))
  # Place the units label at the midpoint of and just below the scale line
  text(x=mean(c(floor(pc[2]-1), floor(pc[2]))), 
       y=pc[3] + 0.15*(pc[4] - pc[3]),
       label=units_label, adj=c(0.5, yadj))
}
# Function to create a color palette based on standard deviation cuts
  create_breaks <- function(raster_data,n) {
    breaks <- c(min(values(raster_data),na.rm=TRUE))
    for (i in seq(1, n)) {
      breaks <- c(breaks, -i * sd(values(raster_data),na.rm=TRUE), i * sd(values(raster_data),na.rm=TRUE))
    }
    breaks <- c(breaks, max(values(raster_data),na.rm=TRUE))
    return(breaks)
  }

# Create a color palette with n standard deviation cuts
mybreaks <- create_breaks(b[[2]],10)

plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Prediction_med_probability.pdf"))
  plot(Africa,main=allt[1],lwd=1.5,col="white",border="grey25")
  plot(myarea,lwd=2,col="grey75",border="grey75",add=TRUE)
  plot(b[[2]],add=TRUE,col=viridis(length(mybreaks)+1),cut=mypal,legend=FALSE)
  plot(b[[2]],add=TRUE,col=viridis(length(mybreaks)+1),cut=mypal,legend.only=TRUE, horizontal = TRUE)
  plot(Africa,lwd=1.5,border="grey25",col='transparent',add=TRUE)
  myScalebar("1 km")
  if(!is.null(dev.list())) dev.off() 
  if(!is.null(dev.list())) dev.off()
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Prediction_sd_probability.pdf"))
plot(Africa,main=allt[5],lwd=1.5,col="white",border="grey25")
plot(myarea,lwd=2,col="grey75",border="grey75",add=TRUE)
plot(b[[5]],add=TRUE,col=viridis(length(mybreaks)+1),cut=mypal,legend=FALSE)
plot(b[[5]],add=TRUE,col=viridis(length(mybreaks)+1),cut=mypal,legend.only=TRUE, horizontal = TRUE)
plot(Africa,lwd=1.5,border="grey25",col='transparent',add=TRUE)
myScalebar("1 km")
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()

#computing the lower bound, med, and upper bound estimation of TB count
#simple way (which would need some improvement for finer assessment,
#e.g. including catchment area etc.)
#computing population size from population density for each pixel
#popdensity<- subset(rs2, 'POP')
#computing population size from population density for each pixel
popdensity <- aggregate(popden, fact=aggfactor)#for computational reasons we make predictions at coarser level
if (mycomb$popfilter==TRUE){
popdensity <- popdensity*popmask} 
#aggregate(popdensity, fact=4)#for computational reasons we make predictions at coarser level
#pop density is in pop per sq.km, so compute area of cell in sq.km

cellarea <- area(popdensity, na.rm=FALSE, weights=FALSE)#area per cell is about 13km
#compute population size as pop.size [pop] = pop.density [pop/sq.km] * area [sq.km]
popsize <- popdensity * cellarea
#putting lower bound, med, and upper bound of TB prevalence estimation together
TBprev <-raster::stack(b[[1]], b[[2]],b[[4]],b[[5]],b[[6]],b[[7]],b[[8]])
TBcount <- TBprev * popsize / prevunit #divide by prevunit to go back to prevalence unit (%)
TBcount <- round(TBcount)
names(TBcount) <- c("LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")
#add informative text before the graphs
c1 = "Predicted lower TB cases"
c2 = "Predicted median TB cases"
c3 = "Predicted upper TB cases"
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Prediction_logmed(count).pdf"))
plot(Africa,main=c2,lwd=2,col="white")
plot(myarea,lwd=2,col="grey55",add=TRUE)
plot(log(TBcount[[2]]),add=TRUE,col=viridis(10),legend=FALSE)
plot(log(TBcount[[2]]),add=TRUE,col=viridis(10),legend.only=TRUE, horizontal = TRUE)
plot(Africa,add=TRUE,lwd=2)
myScalebar("1 km")
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#save maps as raster files
writeRaster(TBcount[[2]], paste0(path_output,'/','tif/TBcount_median.tif'), overwrite=TRUE)
writeRaster(TBcount[[1]], paste0(path_output,'/','tif/TBcount_lower.tif'), overwrite=TRUE)
writeRaster(TBcount[[3]], paste0(path_output,'/','tif/TBcount_upper.tif'), overwrite=TRUE)

#check with official statistics if the estimation is not too far from real estimates for the country
casesest <- data.frame(cellStats(TBcount, sum))
write.csv(casesest,paste0(path_output,"/csv/estimatedcases.csv"))

# #make spplot which has a better output than basic plot
# logTBcount<- log(TBcount)
# names(logTBcount) <- c("logTBcases(low)",
#                        "logTBcases(median)","logTBcases(upper)")
# 
# # Plot with layer-list passed to `sp.layout`
# plot.new()
# pdf.options()
# pdf(file=paste0(path_output,'/',"pdf/Prediction_log(count)_spplot.pdf"))
# print(spplot(logTBcount,sp.layout = Africa.layer,
#         first = TRUE))
# if(!is.null(dev.list())) dev.off() 
# if(!is.null(dev.list())) dev.off()
#extract TB counts per admin regions in Ethiopia

#adm1
Africa.layer <- list("sp.polygons", Africa, col = "black")
adm1poly <- list()
mycountries <- unique(myarea$iso_a3)
for (i in 1:length(mycountries)){
  #adm1poly[[i]] <- raster::getData("GADM", country=mycountries[[i]], path=path_input,level=1)
  # after downlowded you can run the line below:
  adm1poly[[i]] <- readRDS(paste0(path_input, "/gadm36_",mycountries[[i]],"_1_sp.rds"))
}
#put polygons together
adm1poly <- do.call(rbind,adm1poly)

#keep admin1 polygons in the study area
Af1TB <-extract(TBcount, adm1poly,fun="sum",na.rm=TRUE,sp=TRUE)
Af1TB@data <- Af1TB@data[c("GID_0","NAME_0","GID_1","NAME_1","LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")]
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm1_cases.pdf"))
#print(spplot(Af1TB,"median.TB.cases",sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#save the shapefile
library(dplyr)
Af1TB@data <- Af1TB@data %>%
  mutate_if(is.character, ~iconv(., from = "", to = "UTF-8"))

#library(rgdal)
raster::shapefile(Af1TB,file=paste0(path_output,'/',"shp/TBcases_adm1"),overwrite=TRUE)
#rgdal::writeOGR(Af1TB, dsn = paste0(path_output,'/',"shp"), layer = "TBcases_adm1", driver = "ESRI Shapefile", encoding = "UTF-8",overwrite_layer = TRUE)
# Af1TB_sf <- sf::st_as_sf(Af1TB)
# sf::st_write(Af1TB_sf, paste0(path_output,'/',"shp/TBcases_adm1_v2.shp"), encoding = "UTF-8")

#plot the median prevalence at adm-1 level
#keep admin1 polygons in the study area
names(TBprev) <- c("LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")
Af1TBprev <- raster::extract(TBprev, adm1poly,fun="median",na.rm=TRUE,sp=TRUE)
Af1TBprev@data <- Af1TBprev@data[c("GID_0","GID_1", names(TBprev))]
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm1_prevalence.pdf"))
#print(spplot(Af1TBprev, sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm1_median_prevalence.pdf"))
#print(spplot(Af1TBprev,"medianTBprev", sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()

#save the shapefile
raster::shapefile(Af1TBprev,file=paste0(path_output,'/',"shp/TBprevalence_adm1"), overwrite=TRUE)
#rgdal::writeOGR(Af1TBprev, paste0(path_output,'/',"shp"),"TBprevalence_adm1", 
#               driver = "ESRI Shapefile",overwrite_layer=TRUE) 


#adm0
adm0poly <- list()
mycountries <- unique(myarea$iso_a3)
for (i in 1:length(mycountries)){
  #  adm0poly[[i]] <- raster::getData("GADM", country=mycountries[[i]], path=path_input,level=0)
  # after downlowded you can run the line below:
  adm0poly[[i]] <- readRDS(paste0(path_input, "/gadm36_",mycountries[[i]],"_0_sp.rds"))
}
#put polygons together
adm0poly <- do.call(rbind,adm0poly)

#plot the median prevalence at adm-0 level
Af0TBprev <- raster::extract(TBprev, adm0poly,fun="median",na.rm=TRUE,sp=TRUE)
Af0TBprev@data <- Af0TBprev@data[c("GID_0", names(TBprev))]
#plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm0_prevalence.pdf"))
#print(spplot(Af0TBprev, sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#save the shapefile
raster::shapefile(Af0TBprev,file=paste0(path_output,'/',"shp/TBprevalence_adm0"),overwrite=TRUE)
# rgdal::writeOGR(Af0TBprev, paste0(path_output,'/',"shp"),"TBprevalence_adm0",
#          driver = "ESRI Shapefile",overwrite_layer=TRUE,layer_options = "ENCODING=UTF-8")
#plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm0_median_prevalence.pdf"))
#print(spplot(Af0TBprev,"medianTBprev", sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()

#same at ADM0
Af0TB <-extract(TBcount, adm0poly,fun="sum",na.rm=TRUE,sp=TRUE)
Af0TB@data <- Af0TB@data[c("GID_0","NAME_0","LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")]
plot.new()
pdf.options()
pdf(file=paste0(path_output,"/pdf/Africa_adm0_cases.pdf"))
#print(spplot(Af0TB,"median.TB.cases",sp.layout = Africa.layer))
if(!is.null(dev.list())) dev.off() 
if(!is.null(dev.list())) dev.off()
#save the shapefile
raster::shapefile(Af0TB,file=paste0(path_output,'/',"shp/TBcases_adm0"),overwrite=TRUE)

#save csv output
#ADM0
#prevalence
prevadm0 <- data.frame(Af0TBprev@data)
prevadm0 <- merge(prevadm0,adm0poly@data,by='GID_0') 
prevadm0$NAME_0 <- NULL
#count
countadm0 <- data.frame(Af0TB@data)
#merge prevalence and count in the table
prevadm0 <- merge(prevadm0,countadm0,by="GID_0")
write.csv(prevadm0,paste0(path_output,"/csv/prevadm0.csv"))

#ADM1
#prevalence
prevadm1 <- data.frame(Af1TBprev@data)
prevadm1 <- merge(prevadm1,adm1poly@data,by=c('GID_1')) 
prevadm1 <-prevadm1[c("GID_0.x","NAME_0","GID_1","NAME_1","LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")]
colnames(prevadm1)[colnames(prevadm1) == 'GID_0.x'] <- 'GID_0'

#count
countadm1 <- data.frame(Af1TB@data)
countadm1 <- countadm1[c("GID_1","LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")]

#merge prevalence and count in the table
prevadm1 <- merge(prevadm1,countadm1,by="GID_1")
write.csv(prevadm1,paste0(path_output,"/csv/prevadm1.csv"))
}

###exploratory data analysis########
eda <- read.csv('C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/INPUT/TB_Africa_v2.csv') 
# Load libraries
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
library(ggplot2)
# Extract African country boundaries
africa <- ne_countries(continent = "Africa", returnclass = "sf")

eda_sf <- st_as_sf(eda, coords = c("Longitude", "Latitude"), crs = 4326)
eda_sf$Prev <- as.numeric(gsub("%", "", eda_sf$Prev)) / 100
eda_sf$Prev <-eda_sf$Prev*100


# Specify the path where you want to save the shapefile
output_path <- "C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/OUTPUT/shp/eda_sf_shapefile.shp"

# Write the sf object to a shapefile
st_write(eda_sf, output_path)
####scatter plot outlier adm0
all0<-read.csv("C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/OUTPUT/ALL/NOFILTER/csv/prevadm0.csv")
nomoz0<-read.csv("C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/OUTPUT/NOMOZ/NOFILTER/csv/prevadm0.csv")
combined_data <- merge(all0,nomoz0, by = "GID_0", suffixes = c("_1", "_2"))

ggplot(combined_data, aes(x = medianTBprev.x_1, y = medianTBprev.x_2)) +
  geom_point(size = 3) +  # Increase the size of the dots
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # Fit a linear model line
  labs(x = "Model including outlier", y = "Model excluding outlier", 
       title = "Scatter plot of Country level prevalence of TB") +
  theme(
    text = element_text(size = 12)  # Set the font size of all text elements to 12
  )
####scatter plot outlier adm1
all1<-read.csv("C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/OUTPUT/ALL/NOFILTER/csv/prevadm1.csv")
nomoz1<-read.csv("C:/Users/20734635/Desktop/Phd works/TB Africa_Alemneh/Alemneh files/TB Africa - uncertainity - 5km/OUTPUT/NOMOZ/NOFILTER/csv/prevadm1.csv")

combined_data1 <- merge(all1,nomoz1, by = "GID_0", suffixes = c("_1", "_2"))

ggplot(combined_data1, aes(x = medianTBprev.x_1, y = medianTBprev.x_2)) +
  geom_point(size = 3) +  # Increase the size of the dots
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") +  # Add a diagonal line
  labs(x = "Model including outlier", y = "Model excluding outlier", 
       title = "Scatter plot of Region level prevalence of TB") +
  theme(
    text = element_text(size = 12)  # Set the font size of all text elements to 12
  )
p1 <- ggplot(combined_data1, aes(x = medianTBprev.x_1, y = medianTBprev.x_2)) +
  geom_point(size = 3) +  # Increase the size of the dots
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") +  # Add a diagonal line
  labs(x = "Model including outlier", y = "Model excluding outlier", 
       title = "Scatter plot of Region-level Prevalence of TB") +
  theme(
    text = element_text(size = 12)  # Set the font size of all text elements to 12
  ) +
  facet_grid( ~ GID_0)  # Facet by country (GID_0) on rows, and by admin1 region (NAME_1_1) on columns

# Display the plot
plot(p1)


mse <- mean((combined_data1$medianTBprev.x_1 - combined_data1$medianTBprev.x_2)^2)
       # Print the MSE value
     print(paste("Mean Squared Error (MSE):", mse))
     #######COMPARING COVARIATE DISTRIBUTION, GMRF AND PREVALENCE RESULTS
     spatial_field_samples <- inla.posterior.sample(nn, mymodel)
     spatial_field_values <- sapply(1:nn, function(i) spatial_field_samples[[i]]$latent[grep('z.field', rownames(spatial_field_samples[[i]]$latent)),])
    
     A2.pred <- inla.spde.make.A(mesh = mesh, loc = pred_locs)
     predicted_spatial_field <- A2.pred %*% spatial_field_values
     spatial_field_mean <- apply(predicted_spatial_field, 1, mean)
     transformed <- linkfun(spatial_field_mean)
      # Aggregate the posterior samples to obtain mean spatial random effects
    
     spatial_field_data <- data.frame(
       longitude = pred_locs[, 1],  # Assuming longitude is in the first column of pred_locs
       latitude = pred_locs[, 2],   # Assuming latitude is in the second column of pred_locs
       transfoemed_mean = transformed  # The mean of the estimated GMRF
     )
     # Plot the estimated GMRF
     ggplot() +
       geom_tile(data = spatial_field_data, aes(x = longitude, y = latitude, fill = transfoemed_mean)) +
       scale_fill_viridis_c() +
       theme_minimal() +
       labs(title = "Estimated Spatial Random Effect (GMRF)")
     
     pred2 = matrix(NA,nrow=dim(A.pred)[1], ncol=nn)
     for (i in 1:nn){
       field = samp[[i]]$latent[grep('z.field',rownames(samp[[i]]$latent)),]
       lp =  drop(A.pred%*%field)
       ## Predicted values
       pred[,i] = linkfun(lp) #for transformation from 0 to mu
       # pred[,i] = plogis(lp) #for a bernoulli likelihood
       # pred[,i] = exp(lp)#for a poisson likelihood
     }
     pred_mean_field <- apply(pred, 1, function(x) mean(x, na.rm = TRUE))
     predinput=list(pred_mean_field)
     prednames=as.list(c("Prob_mean_field"))
     mainnames=as.list(c("mean"))
     #destfile <- paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket
     
     out<-list()
     for (j in 1:length(predinput))
     {
       pred_val[!w] <- predinput[[j]]
       out[[j]] = setValues(mask, pred_val)
       #mask with population density
       if (mycomb$popfilter==TRUE){
         out[[j]] = out[[j]]* popmask * prevunit} else {
           out[[j]] = out[[j]]* prevunit} 
       writeRaster(out[[j]], paste0(path_output,'/',"prevalence/",prednames[[j]],'.tif'), overwrite=TRUE)
     }
     library(raster)
     
