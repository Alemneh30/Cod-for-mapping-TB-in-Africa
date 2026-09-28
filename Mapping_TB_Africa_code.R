# Code to generate the main results 
# Default parameters when this script is run independently --------------------

if (!exists("nn"))            nn <- 10000
if (!exists("mainaggfactor")) mainaggfactor <- 2
if (!exists("mu"))            mu <- 0.025
if (!exists("prevunit"))      prevunit <- 1000
if (!exists("popt"))          popt <- 5
if (!exists("allrun"))        allrun <- FALSE
# run all model specs or only the main spec.
if (allrun) {
  combinations <- expand.grid(popfilter = c(FALSE, TRUE),mozout = c(FALSE, TRUE))
} else {
  # Main specification only
  combinations <- data.frame(popfilter = FALSE,mozout = FALSE)
}
###############end user specifications#########################################

#INLA download#########################################################
#Note that using recent INLA versions may lead to different results 
#To replicate the work, you can use INLA 24.05.10 from tgz (archives). For mac (M processor):
#install.packages("~/Downloads/INLA_24.05.10.tgz",repos = NULL, type = "mac.binary")#from https://inla.r-inla-download.org/R/stable
# get fmesher compatible with INLA 24.05.10 from: https://cran.ms.unimelb.edu.au/src/contrib/Archive/fmesher/
# install.packages(
#   "~/Downloads/fmesher_0.1.7.tar.gz",
#   repos = NULL,
#   type = "source"
# )
###############end INLA download#########################################

#INLA used to fit Bayesian models
list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","viridis", "geodata", "rnaturalearth", "malariaAtlas", "readxl","ggplot2",
                      "RColorBrewer", "ggmap", "tmap","fmsb","geodata","sf","rnaturalearthdata",
                      "ggrepel","dplyr","gtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
library(INLA)

#2. step: define the input and output paths
# To facilitate reproduciblity, all folders are created in the current working directory.
path_input <- paste0(getwd(),"/INPUT")
output_path <- paste0(getwd(),"/OUTPUT")
dir.create(path_input,  recursive = TRUE, showWarnings = FALSE)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


#loop
if(allrun==TRUE){
  combinations <-expand.grid(popfilter=c(TRUE,FALSE),mozout=c(TRUE,FALSE))
  #folder for Mozambique robustness
   dir.create(file.path(output_path,"robustness_Mozambique"), recursive = TRUE, showWarnings = FALSE)
}else{
  combinations <-data.frame(popfilter=popfilter,mozout=mozout)
}
#intial data loading
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
TB0 <- TB0[complete.cases(TB0[,c("ISO3")]),]

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
# create the output directories if they don't exist
output_dirs <- c(path_output, file.path(path_output, "shp"),
  file.path(path_output, "pdf"), file.path(path_output, "csv"),
  file.path(path_output, "prevalence"),file.path(path_output, "tif") )
invisible(lapply(output_dirs, dir.create,recursive = TRUE, showWarnings = FALSE))

#Robustness test: whether to remove Mozambique or not
if (mycomb$mozout==TRUE){
  TB <- subset(TB0,ISO3!= "MOZ")
} else {TB <-TB0}

#select countries in the study area based on worldbank ISO 3
#identify countries to be removed 
#study area
Africa <- rnaturalearth::ne_countries()
Africa <- Africa[Africa$continent=="Africa",]
outc <- setdiff(unique(Africa$wb_a3), unique(TB$ISO3))#vector 1,2 as arguments
outc <- outc[complete.cases(outc)]
#remove iteratively countries not in TB data
myarea <- Africa
for (i in 1:length(outc)){
myarea = subset(myarea,wb_a3 != outc[i])
}
#check if the right countries have been selected
if(isFALSE(unique(sort(myarea$wb_a3))==unique(sort(TB$ISO3))))
  stop("Error: countries of TB do not match study area!")
#for raster operations put myarea into a SPDF object
myareaplot <- myarea

# load orginal covariates at 10m resolution (from the original study)
# covariates are: T1=temperature, P1=precipitation, alt=altitude, ahf=accessibility to health care, acc=accessibility to cities, popden=population density
load(file.path(path_input, "originalcovariates.Rdata"))
rs <- raster::stack(T1, P1, alt, ahf, acc, popden)
#rs <- raster::crop(rs, raster::extent(myarea))
#rs <- raster::mask(rs, myarea)
names(rs) <- c("TMP", "PCP", "ALT", "ACH", "ACC", "POP")

#*************x-mean/sd to standardize the covariates rs ****************************************
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
rs2 <- stack(st2)#

#TB data
pt <-TB
pt <- data.frame(TB=pt$`Total TB cases`,Nscreen = pt$Number_examined, x=pt$Longitude, y=pt$Latitude)
pt <- pt[complete.cases(pt),]

#extract covariate at point coordinates
xy <- cbind(pt$x, pt$y)
covariate_all <- data.frame(raster::extract(rs2, xy))
covariate_z <- data.frame(covariate_all)
#apply VIF to subset data
source(paste0(path_input,"/vif.R"))
#The function uses three arguments. The first is a matrix or data frame of the explanatory variables,
#the second is the threshold value to use for retaining variables, and the third is a logical argument 
#indicating if text output is returned as the stepwise selection progresses. The output indicates the VIF
#values for each variable after each stepwise comparison. The function calculates the VIF values for all 
#explanatory variables, removes the variable with the highest value, and repeats until all VIF values are
#below the threshold. The final output is a list of variable names with VIF values that fall below the threshold. 
keep.dat <- vif_func(in_frame=covariate_z,thresh=5,trace=F)
#temperature and altitude are also highly correlated (so remove altitude)
rs <- rs[[keep.dat]]
rs2 <- rs2[[keep.dat]]
covariate_z <- covariate_z[paste0(keep.dat)]

#correlation
corr_matrix <- cor(covariate_z,use="complete.obs")
write.csv(corr_matrix,paste0(path_output,"/csv/covariate_correlation.csv"))
#pop and acc highly correlated so we remove pop
#ach and acc highly correlated so we remove ach
covariate_z$POP <- NULL
covariate_z$ACH <- NULL

#combine all data together
ptcov <- cbind(pt,covariate_z)
#summary(ptcov)

#######MODELLING #####################################################################
######################################################################################
bdry <- INLA::inla.sp2segment(myarea)
bdry$loc <- INLA::inla.mesh.map(bdry$loc)
mesh<-INLA::inla.mesh.2d(loc=xy, boundary=bdry, max.edge=c(0.6,4),offset=c(0.6,4),cutoff = 0.6)#to be fine tuned
#mesh$n
#save mesh plot
tv <- mesh$graph$tv
mesh_segments <- do.call(rbind, lapply(seq_len(nrow(tv)), function(i) {v <- tv[i,]; rbind(data.frame(x=mesh$loc[v[1],1],y=mesh$loc[v[1],2],xend=mesh$loc[v[2],1],yend=mesh$loc[v[2],2]), data.frame(x=mesh$loc[v[2],1],y=mesh$loc[v[2],2],xend=mesh$loc[v[3],1],yend=mesh$loc[v[3],2]), data.frame(x=mesh$loc[v[3],1],y=mesh$loc[v[3],2],xend=mesh$loc[v[1],1],yend=mesh$loc[v[1],2]))}))
myarea_sf <- st_as_sf(myarea)
country_labels <- st_point_on_surface(myarea_sf) |> st_coordinates() |> as.data.frame() |> mutate(ISO3=myarea_sf$iso_a3)
africa_sf <- sf::st_as_sf(Africa)
#plot
p_mesh <- ggplot() +
  geom_sf(data=africa_sf,fill="grey85",colour="grey45",linewidth=.7) +
  geom_segment(data=mesh_segments,aes(x=x,y=y,xend=xend,yend=yend),colour="navy",linewidth=.15,alpha=.45) +
  #geom_sf(data=myarea_sf,fill=NA,colour="navy",linewidth=.7) +
  geom_label_repel(data=country_labels,aes(x=X,y=Y,label=ISO3),size=3,fontface="bold",fill="white",label.size=.15,box.padding=.25,point.padding=.05,min.segment.length=0,seed=123) +
  coord_sf(expand=FALSE) +
   scale_y_continuous(breaks=seq(-40,40,20)) +
  scale_x_continuous(breaks=seq(-20,40,20)) +
  theme_classic()+ xlab("")+ ylab("")+
  theme(plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=8),panel.grid=element_line(colour="grey90",linewidth=.3))
#save plot
ggsave(file.path(path_output,"pdf","mesh.pdf"),p_mesh,width=9,height=8)

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
  f1 <- as.formula(paste0("Y ~  f(z.field, model=spde) + ", paste0(colnames(covariate_z)[1:i], collapse = " + ")))
  model1<-inla(f1, #the formula
               data=inla.stack.data(stk,spde=spde),  #the data stack
               family= 'binomial',   #which family the data comes from
               Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
               control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
              # control.compute = list(mlik = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               control.compute = list(waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              # inla.mode = "classic",
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
#select the best model based on marginal likelihood
modsel <- model_selection[which.min(model_selection[,2]),]#min for waic
#modsel <- model_selection[which.max(model_selection[,2]),]#max for mlik

#spatial field only
f0 <- as.formula(paste0("Y ~ f(z.field, model=spde)"))
stk0 = inla.stack(data=list(Y=ptcov$TB, n=ptcov$Nscreen),A=list(A),effects=list(
  list(z.field=1:spde$n.spde)))
model0<-inla(f0, #the formula
             data=inla.stack.data(stk0,spde=spde),  #the data stack
             family= 'binomial',   #which family the data comes from
             Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
             control.predictor=list(A=inla.stack.A(stk0),compute=TRUE),  #compute gives you the marginals of the linear predictor
             #control.compute = list(mlik = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             control.compute = list(waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             control.fixed = list(prec.intercept=.001),
             #   inla.mode = "classic",
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
df0 <- data.frame(model = "Spatial field (no covariates)",waic = model0$waic$waic)
#df0 <- data.frame(model = "Spatial field (no covariates)",mlik = model0$mlik[1])
df0[,2] <- as.numeric(df0[,2])
model_selection <- rbind(model_selection,df0)

#model with intercept only
f0 <- as.formula(paste0("Y ~  -1 + z.intercept"))
stk0 = inla.stack(data=list(Y=ptcov$TB, n=ptcov$Nscreen),A=list(A),effects=list(
  list(z.intercept=1)))
model0<-inla(f0, #the formula
             data=inla.stack.data(stk0),  #the data stack
             family= 'binomial',   #which family the data comes from
             Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
             control.predictor=list(A=inla.stack.A(stk0),compute=TRUE),  #compute gives you the marginals of the linear predictor
             #control.compute = list(mlik = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             control.compute = list(waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             control.fixed = list(prec.intercept=.001),
             #   inla.mode = "classic",
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
df0 <- data.frame(model = "Intercept only",waic = model0$waic$waic)
#df0 <- data.frame(model = "Intercept only",mlik = model0$mlik[1])
df0[,2] <- as.numeric(df0[,2])
model_selection <- rbind(model_selection,df0)
model_selection <- model_selection[order(model_selection$waic, decreasing = TRUE), ]
#model_selection <- model_selection[order(model_selection$mlik, decreasing = TRUE), ]
#save model selection results
write.csv(model_selection,paste0(path_output,"/csv/Supplementary_Table2.csv"))

#Run the selected model based on the results of the model selection
formula.z = as.formula(paste("Y ~  f(z.field, model=spde) +",modsel$model))

mymodel<-inla(formula.z, #the formula
            data=inla.stack.data(stk,spde=spde),  #the data stack
            family= 'binomial',   #which family the data comes from
            Ntrials = n,      #this is specific to binomial as we need to tell it the number of examined
            control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            control.fixed = list(prec=1000,prec.intercept=.001),
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
#summary(mymodel)
#improve cpo computation (optional, takes about 10-15mn for 143 cases)
if(mymodel$ok==FALSE){
mymodel = inla.cpo(mymodel, force=TRUE)
}
#save coefficient correlation 
modeloutput_OR <- mymodel$summary.fixed
# Remove mode and kld and sd
modeloutput_OR <- modeloutput_OR[, 
  !names(modeloutput_OR) %in% c("mode", "kld", "sd")]
# Transform log-odds to odds ratios
modeloutput_OR[] <- lapply(modeloutput_OR, exp)
# Round it to 3 decimals
modeloutput_OR <- round(modeloutput_OR, 3)

library(dplyr)
Table2 <- modeloutput_OR %>%
  mutate(
    Covariate = rownames(modeloutput_OR),
    `Odds Ratio (95% CrI)` = sprintf(
      "%.2f (%.2f, %.2f)", mean, `0.025quant`, `0.975quant`)
  ) %>%
  select(Covariate, `Odds Ratio (95% CrI)`) %>%
  arrange(Covariate == "(Intercept)")
rownames(Table2) <- NULL
#save Table 2
write.csv(Table2,paste0(path_output,"/csv/Table2.csv"))

#validity check plot
#The PIT is the probability of a new response less than the observed response using a model based on the rest of the data. 
#We'd expect the PIT values to be uniformly distributed if the model assumptions are correct.
#check http://julianfaraway.github.io/brinla/examples/chicago.html for info
pit <- mymodel$cpo$pit
n <- length(Y)
uniquant <- (1:n)/(n+1)
#Report the plot with the logit transform
pit_df <- data.frame(x = logit(uniquant),y = logit(sort(pit)))
#plot
p_pit <- ggplot(pit_df, aes(x=x, y=y)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=.5) +
  geom_point(size=1.25, shape=21, fill="white", colour="black", stroke=0.5)+
  labs(x="uniform quantiles (logit scale)",y="Sorted PIT values (logit scale)") +
  scale_x_continuous(breaks=seq(-4,4,2)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    plot.title=element_text(face="bold"),
    plot.caption=element_text(hjust=0, size=8),
    panel.grid.minor=element_blank()
  )
  #saveplot
ggsave(file.path(path_output,"pdf","PITvalidity.pdf"),p_pit, width=6, height=4)

#Predictions (mapping) at prediction locations
#Use coarser grid cells to account for spatial uncertainty and to reduce computational burden
#aggregation factor is defined by the user at the beginning of the code
mask <- rs2[[1]]
raster::NAvalue(mask) <- -9999
mask <- raster::aggregate(mask, fact = mainaggfactor)
#create population mask as well for filtering areas where pop < threshold
popmask <- aggregate(popden, fact=mainaggfactor)#for computational reasons we make predictions at coarser level

#filter by population density if requested by the user
if(mycomb$popfilter==TRUE){
popmask[popmask < popt] <- NA
}
popmask[!is.na(popmask),] <- 1
#plot(popmask)

#code to make maps
pred_val <- getValues(mask)
w <- is.na(pred_val)
index <- 1:length(w)
index <- index[!w]
pred_locs <- xyFromCell(mask,1:ncell(mask))
pred_locs <- pred_locs[!w,]
colnames(pred_locs)<-c('longitude','latitude')
locs_pred <- pred_locs
#Mapping between meshes and continuous space
A.pred = inla.spde.make.A(mesh=mesh, loc=locs_pred)
mycovnames <- mymodel$names.fixed
mycovnames <- mycovnames[-1]
covariates_pred <- rs2[[mycovnames]]
covariates_pred = data.frame(raster::extract(covariates_pred, pred_locs))
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
  #predicted values
   pred[,i] = linkfun(lp) #for transformation from 0 to mu
 # pred[,i] = plogis(lp) #for a bernoulli likelihood
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
#save raster files of prevalence
out<-list()
for (j in 1:length(predinput)){
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
  rasterls[[i]]<-raster::raster(raster.list[[i]])}
b <- raster::brick(rasterls)

#computing the lower bound, med, and upper bound estimation of TB count
#computing population size from population density for each pixel of the predictive map
popdensity <- raster::aggregate(popden, fact=mainaggfactor)#for computational reasons we make predictions at coarser level

#population filter if requested by the user
if (mycomb$popfilter==TRUE){
popdensity <- popdensity*popmask} 

#pop density is in pop per sq.km, so compute area of cell in sq.km
cellarea <- raster::area(popdensity, na.rm=FALSE, weights=FALSE)#area per cell 
#compute population size as pop.size [pop] = pop.density [pop/sq.km] * area [sq.km]
popsize <- popdensity * cellarea
#putting lower bound, med, and upper bound of TB prevalence estimation together
TBprev <-raster::stack(b[['LUI']], b[['prob_med']],b[['Prob_q25']],b[['Prob_q75']],b[['Prob_IQR']],b[['Prob_mean']],b[['UUI']])
TBcount <- TBprev * popsize / prevunit #divide by prevunit to go back to prevalence unit 
TBcount <- round(TBcount)
names(TBcount) <- c("LITBcount","medianTBcount", "lowqTBcount","uppqTBcount","IQRTBcount","meanTBcount","UITBcount")
#save maps as raster files
writeRaster(TBcount[['medianTBcount']], paste0(path_output,'/','tif/TBcount_median.tif'), overwrite=TRUE)
writeRaster(TBcount[['meanTBcount']], paste0(path_output,'/','tif/TBcount_mean.tif'), overwrite=TRUE)
writeRaster(TBcount[['lowqTBcount']], paste0(path_output,'/','tif/TBcount_lower.tif'), overwrite=TRUE)
writeRaster(TBcount[['uppqTBcount']], paste0(path_output,'/','tif/TBcount_upper.tif'), overwrite=TRUE)

#generate table of global estimates
casesest <- data.frame(cellStats(TBcount, sum))
table_global <- data.frame(
  Area = "Global", `Mean estimated TB cases (95% UI)` = paste0(
    formatC(casesest["meanTBcount", 1], format = "d", big.mark = ","),
    " (", formatC(casesest["LITBcount", 1], format = "d", big.mark = ","),
    "–",formatC(casesest["UITBcount", 1], format = "d", big.mark = ","),")"),
  Median = formatC( casesest["medianTBcount", 1],format = "d",big.mark = ","),
  IQR = formatC(casesest["IQRTBcount", 1],format = "d", big.mark = ","),
  check.names = FALSE
)
#save global estimates
write.csv(table_global,paste0(path_output,"/csv/globalestimates.csv"))

#extract predicted values at adm1 level
Africa.layer <- list("sp.polygons", Africa, col = "black")
adm1poly <- list()
mycountries <- unique(myarea$iso_a3)
for (i in 1:length(mycountries)){
  adm1poly[[i]] <- readRDS(paste0(path_input, "/gadm36_",mycountries[[i]],"_1_sp.rds"))
}
#put polygons together
adm1poly <- do.call(rbind,adm1poly)

#extract predicted count TB at adm1 level
Af1TB <- raster::extract(TBcount, adm1poly,fun="sum",na.rm=TRUE,sp=TRUE)
Af1TB@data <- Af1TB@data[c("GID_0","NAME_0","GID_1","NAME_1","LITBcount","medianTBcount", "lowqTBcount","uppqTBcount","IQRTBcount","meanTBcount","UITBcount")]
#save the shapefile
library(dplyr)
Af1TB@data <- Af1TB@data %>%
  mutate_if(is.character, ~iconv(., from = "", to = "UTF-8"))
# save the shapefile
raster::shapefile(Af1TB,file=paste0(path_output,'/',"shp/TBcases_adm1"),overwrite=TRUE)

#extract predicted prevalence TB at adm1 level
names(TBprev) <- c("LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")
Af1TBprev <- raster::extract(TBprev, adm1poly,fun="mean",na.rm=TRUE,sp=TRUE)
Af1TBprev@data <- Af1TBprev@data[c("GID_0","GID_1", names(TBprev))]
#save the shapefile
raster::shapefile(Af1TBprev,file=paste0(path_output,'/',"shp/TBprevalence_adm1"), overwrite=TRUE)

#extract predicted values at adm0 level
adm0poly <- list()
mycountries <- unique(myarea$iso_a3)
for (i in 1:length(mycountries)){
  adm0poly[[i]] <- readRDS(paste0(path_input, "/gadm36_",mycountries[[i]],"_0_sp.rds"))
}
#put polygons together
adm0poly <- do.call(rbind,adm0poly)

#extract predicted prevalence at adm0 level
Af0TBprev <- raster::extract(TBprev, adm0poly,fun="mean",na.rm=TRUE,sp=TRUE)
Af0TBprev@data <- Af0TBprev@data[c("GID_0", names(TBprev))]
#save the shapefile
raster::shapefile(Af0TBprev,file=paste0(path_output,'/',"shp/TBprevalence_adm0"),overwrite=TRUE)

#extract predicted count TB at adm0 level
Af0TB <- raster::extract(TBcount, adm0poly,fun="sum",na.rm=TRUE,sp=TRUE)
Af0TB@data <- Af0TB@data[c("GID_0","NAME_0","LITBcount","medianTBcount", "lowqTBcount","uppqTBcount","IQRTBcount","meanTBcount","UITBcount")]
#save the shapefile
raster::shapefile(Af0TB,file=paste0(path_output,'/',"shp/TBcases_adm0"),overwrite=TRUE)

#save predicted prevalence and count at adm0 and adm1 level in csv files
#prevalence at adm0 level
prevadm0 <- data.frame(Af0TBprev@data)
prevadm0 <- merge(prevadm0,adm0poly@data,by='GID_0') 
prevadm0$NAME_0 <- NULL
#count at adm0 level
countadm0 <- data.frame(Af0TB@data)
#merge prevalence and count in the table
prevadm0 <- merge(prevadm0,countadm0,by="GID_0")
write.csv(prevadm0,paste0(path_output,"/csv/prev_and_count_adm0.csv"))

#Table 1 (Estimated TB cases by country)
library(dplyr)
table1 <- prevadm0 %>%
  select(Country = NAME_0,LITBcount,medianTBcount,
  meanTBcount,IQRTBcount, UITBcount) %>%
  arrange(desc(meanTBcount)) %>%
  mutate(
    `Estimated mean TB cases (95% UI)` = paste0(
      formatC(meanTBcount, format = "d", big.mark = ","),
      " (",formatC(LITBcount, format = "d", big.mark = ","), "–",
      formatC(UITBcount, format = "d", big.mark = ","), ")"),
    Median = formatC( medianTBcount, format = "d", big.mark = ","),
    IQR = formatC(IQRTBcount, format = "d", big.mark = ",")) %>%
  select( Country, `Estimated mean TB cases (95% UI)`, Median, IQR )
# Save Table 1
write.csv(table1,paste0(path_output, "/csv/Table1.csv"),row.names = FALSE)

#Table 1b (Estimated TB prevalence by country)
table1b <- prevadm0 %>%
  select(Country = NAME_0,LITBprev,medianTBprev, meanTBprev, IQRTBprev,UITBprev ) %>%
  arrange(desc(meanTBprev)) %>%
  mutate(
    `Estimated mean TB prevalence per 1,000 (95% UI)` = paste0(
      formatC(meanTBprev, format = "f", digits = 2), " (",
      formatC(LITBprev, format = "f", digits = 2), "–",
      formatC(UITBprev, format = "f", digits = 2), ")"),
    Median = formatC( medianTBprev, format = "f", digits = 2),
    IQR = formatC(IQRTBprev, format = "f", digits = 2 )) %>%
  select(Country, `Estimated mean TB prevalence per 1,000 (95% UI)`, Median,IQR )

# Save Table 1b
write.csv(table1b,paste0(path_output, "/csv/Table1b.csv"),row.names = FALSE)

#extract predicted values at adm1 level and save in csv files
#prevalence
prevadm1 <- data.frame(Af1TBprev@data)
prevadm1 <- merge(prevadm1,adm1poly@data,by=c('GID_1')) 
prevadm1 <-prevadm1[c("GID_0.x","NAME_0","GID_1","NAME_1","LITBprev","medianTBprev", "lowqTBprev","uppqTBprev","IQRTBprev","meanTBprev","UITBprev")]
colnames(prevadm1)[colnames(prevadm1) == 'GID_0.x'] <- 'GID_0'

#count
countadm1 <- data.frame(Af1TB@data)
countadm1 <- countadm1[c("GID_1","LITBcount","medianTBcount", "lowqTBcount","uppqTBcount","IQRTBcount","meanTBcount","UITBcount")]

#merge prevalence and count in the table
prevadm1 <- merge(prevadm1,countadm1,by="GID_1")
write.csv(prevadm1,paste0(path_output,"/csv/prev_and_count_adm1.csv"))

#Analysing the spatial field (GMRF) and plotting the estimated spatial random effect 
spatial_field_samples <- inla.posterior.sample(nn, mymodel)
spatial_field_values <- sapply(1:nn, function(i) spatial_field_samples[[i]]$latent[grep('z.field', rownames(spatial_field_samples[[i]]$latent)),])
A2.pred <- inla.spde.make.A(mesh = mesh, loc = pred_locs)
predicted_spatial_field <- A2.pred %*% spatial_field_values
spatial_field_mean <- apply(predicted_spatial_field, 1, mean)
transformed <- linkfun(spatial_field_mean)
# Aggregate the posterior samples to obtain mean spatial random effects
spatial_field_data <- data.frame(
  longitude = pred_locs[, 1],  # longitude is in the first column of pred_locs
  latitude = pred_locs[, 2],   # latitude is in the second column of pred_locs
  transformed_mean = transformed  # The mean of the estimated GMRF
)
# Create raster from prediction grid
z <- raster::getValues(mask)
z[!w] <- transformed
spatial_field_raster <- raster::setValues(mask, z)

# save raster
raster::writeRaster(
  spatial_field_raster,
  paste0(path_output, "/prevalence/spatial_field_mean.tif"),
  overwrite=TRUE)
}
##########################################################################
#end loop across all model specifications#################################
###Analysis of effects by removing Mozambique from the analysis
if(allrun == TRUE){
# load the CSV files for both scenarios (with and without Mozambique)
all0 <- read.csv(paste0(output_path,"/ALL/NOFILTER/csv/prev_and_count_adm0.csv"))
nomoz0 <- read.csv(paste0(output_path,"/NOMOZ/NOFILTER/csv/prev_and_count_adm0.csv"))
combined_data <- merge(all0,nomoz0, by = "GID_0", suffixes = c("_1", "_2"))

#adm-0 level of analysis
#compare mean TB prevalence when Mozambique is included vs excluded from the analysis
rmse <- sqrt(mean(
  (combined_data$meanTBprev_1 - combined_data$meanTBprev_2)^2,
  na.rm=TRUE
))
# prepare data
plot_prev_adm0 <- combined_data %>%
  dplyr::select(GID_0,
   meanTBprev_1, LITBprev_1, UITBprev_1,
   meanTBprev_2, LITBprev_2, UITBprev_2) %>%
  tidyr::pivot_longer(
    cols=-GID_0,
    names_to=c(".value","model"),
    names_pattern="(meanTBprev|LITBprev|UITBprev)_([12])"
  ) %>%
  mutate(model=factor(model, levels=c("1","2"),
                      labels=c("Mean value using model\nincluding Mozambique",
                      "Mean value using model\nexcluding Mozambique")))

p0_prev <- ggplot(plot_prev_adm0, aes(x=meanTBprev, y=reorder(GID_0,meanTBprev), shape=model)) +
  geom_errorbarh(aes(xmin=LITBprev,xmax=UITBprev), height=.15,
                 position=position_dodge(width=.5), linewidth=.5) +
  geom_point(position=position_dodge(width=.5), size=2.5) +
  annotate("label", x=Inf, y=-Inf,
           label=paste0("RMSE = ", round(rmse,3), " per 1,000"),
           hjust=1.1, vjust=-.5, size=3.5) +
  labs(
    x="ADM-0 estimated mean TB prevalence per 1,000 (95% UI)",
    y="",
    shape=NULL,
    title=""
  ) +
  theme_bw() +
  theme(text=element_text(size=11), legend.position="bottom",
        panel.grid.minor=element_blank())
#save
ggsave(p0_prev, filename=paste0(output_path,"/robustness_Mozambique/rob_Moz_mean_prev_adm0.pdf"), width=5.5, height=5)
#adm-1 level of analysis
#compare mean TB prevalence when Mozambique is included vs excluded from the analysis
all1 <- read.csv(paste0(output_path,"/ALL/NOFILTER/csv/prev_and_count_adm1.csv"))
nomoz1 <- read.csv(paste0(output_path,"/NOMOZ/NOFILTER/csv/prev_and_count_adm1.csv"))  
combined_data1 <- merge(all1, nomoz1, by = c("GID_0", "NAME_1"), suffixes = c("_1", "_2"))
# ADM1: simple comparison of mean TB prevalence -----------------------------

rmse_adm1 <- sqrt(mean(
  (combined_data1$meanTBprev_1 - combined_data1$meanTBprev_2)^2,
  na.rm=TRUE
))

pearson_adm1 <- cor(combined_data1$meanTBprev_1,combined_data1$meanTBprev_2,
  use="complete.obs",method="pearson")

p1 <- ggplot(combined_data1, aes(x=meanTBprev_2, y=meanTBprev_1)) +
  geom_point(size=.8, alpha=.4) +
  geom_abline(intercept=0, slope=1, colour="black", linetype="dashed", linewidth=.5) +
  annotate("label", x=Inf, y=-Inf,
           label=paste0("Pearson r = ", round(pearson_adm1,3),
                        "\nRMSE = ", round(rmse_adm1,3), " per 1,000"),
           hjust=1.05, vjust=-.5, size=3.5) +
  labs(
    x="Model excluding Mozambique\nADM-1 estimated mean TB prevalence per 1,000",
    y="Model including Mozambique\nADM-1 estimated mean TB prevalence per 1,000",
    title=""
  ) +
  coord_equal() +
  theme_bw() +
  theme( panel.grid.minor=element_blank())

ggsave(p1,
       filename=paste0(output_path,"/robustness_Mozambique/rob_Moz_mean_prev_adm1.pdf"),
       width=5.5, height=5)

# ADM0 TB counts ---------------------------------------------------------------

rmse_count <- sqrt(mean(
  (combined_data$meanTBcount_1 - combined_data$meanTBcount_2)^2,
  na.rm=TRUE
))

# prepare data
plot_count_adm0 <- combined_data %>%
  dplyr::select(GID_0,
    meanTBcount_1, LITBcount_1, UITBcount_1,
    meanTBcount_2, LITBcount_2, UITBcount_2) %>%
  tidyr::pivot_longer(
    cols=-GID_0,
    names_to=c(".value","model"),
    names_pattern="(meanTBcount|LITBcount|UITBcount)_([12])"
  ) %>%
        mutate(model=factor(model, levels=c("1","2"),
        labels=c("Mean value using model\nincluding Mozambique",
       "Mean value using model\nexcluding Mozambique")))

p0_count <- ggplot(plot_count_adm0,
                   aes(x=meanTBcount, y=reorder(GID_0,meanTBcount), shape=model)) +
  geom_errorbarh(aes(xmin=LITBcount,xmax=UITBcount), height=.15,
                 position=position_dodge(width=.5), linewidth=.5) +
  geom_point(position=position_dodge(width=.5), size=2.5) +
  annotate("label", x=Inf, y=-Inf,
           label=paste0("RMSE = ", format(round(rmse_count), big.mark=","), " cases"),
           hjust=1.1, vjust=-.5, size=3.5) +
  scale_x_log10(labels=scales::label_comma()) +
  labs(
    x="ADM-0 Estimated mean TB cases (95% UI; log10 scale)",
    y="",
    shape=NULL,
    title=""
  ) +
  theme_bw() +
  theme(text=element_text(size=11), legend.position="bottom",
        panel.grid.minor=element_blank())

ggsave(
  p0_count,
  filename=paste0(output_path,"/robustness_Mozambique/rob_Moz_mean_count_adm0.pdf"),
  width=5.5
)

# -----------------------------------------------------------------------------
# TB COUNT: ADM1
# -----------------------------------------------------------------------------

rmse_count_adm1 <- sqrt(mean(
  (combined_data1$meanTBcount_1 - combined_data1$meanTBcount_2)^2,
  na.rm=TRUE
))

pearson_count_adm1 <- cor(combined_data1$meanTBcount_1,combined_data1$meanTBcount_2,
  use="complete.obs",method="pearson")

p1_count <- ggplot(combined_data1, aes(x=meanTBcount_2, y=meanTBcount_1)) +
  geom_point(size=.8, alpha=.4) +
  geom_abline(intercept=0, slope=1, colour="black", linetype="dashed", linewidth=.5) +
  annotate("label", x=Inf, y=-Inf,
           label=paste0("Pearson r = ", round(pearson_count_adm1,3),
                        "\nRMSE = ", format(round(rmse_count_adm1), big.mark=","), " cases"),
           hjust=1.05, vjust=-.5, size=3.5) +
  labs(
    x="Model excluding Mozambique\nADM-1 Estimated mean TB cases",
    y="Model including Mozambique\nADM-1 Estimated mean TB cases",
    title=""
  ) +
  coord_equal() +
  theme_bw() +
  theme( panel.grid.minor=element_blank())

ggsave(
  p1_count,
  filename=paste0(output_path,"/robustness_Mozambique/rob_Moz_mean_count_adm1.pdf"),
  width=4, height=4.5)

}
# End of the code########################################################################
