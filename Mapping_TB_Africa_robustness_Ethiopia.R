# =============================================================================
# ROBUSTNESS: INFLUENCE OF ETHIOPIAN OBSERVATIONS ON TB CASE COUNTS
# Fits the same model with all data and excluding Ethiopia, then compares
# estimated TB case counts outside Ethiopia at ADM0 and ADM1 levels.
# =============================================================================

# 1. SETTINGS -----------------------------------------------------------------

if (!exists("nn")) nn <- 10000
if (!exists("mainaggfactor")) mainaggfactor <- 2
if (!exists("mu")) mu <- 0.025
if (!exists("prevunit")) prevunit <- 1000
set.seed(999)

pkgs <- c("INLA","raster","rnaturalearth","readxl","ggplot2","dplyr","ggrepel")
invisible(lapply(pkgs, require, character.only=TRUE))

path_input <- file.path(getwd(),"INPUT")
path_output <- file.path(getwd(),"OUTPUT","robustness_Ethiopia")
dir.create(file.path(path_output,"csv"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(path_output,"pdf"),recursive=TRUE,showWarnings=FALSE)


# 2. DATA ---------------------------------------------------------------------

TB0 <- readxl::read_excel(file.path(path_input,"TB_Africa_v2.xlsx"),trim_ws=TRUE)
TB0$Number_examined <- round(as.numeric(TB0$Number_examined))
TB0$`Total TB cases` <- round(as.numeric(TB0$`Total TB cases`))
TB0$Latitude <- as.numeric(TB0$Latitude); TB0$Longitude <- as.numeric(TB0$Longitude)
TB0 <- TB0[TB0$`Lowest admin level (0-4)` > 0 & complete.cases(TB0$ISO3),]
TB0$Number_examined <- pmax(TB0$Number_examined,0)

Africa <- rnaturalearth::ne_countries()
Africa <- Africa[Africa$continent=="Africa",]

load(file.path(path_input,"originalcovariates.Rdata"))
rs <- raster::stack(T1,P1,alt,ahf,acc,popden)
names(rs) <- c("TMP","PCP","ALT","ACH","ACC","POP")


# 3. QUANTILE-NORMALISE COVARIATES --------------------------------------------

st2 <- lapply(seq_len(raster::nlayers(rs)),function(i) {
  x <- raster::getValues(rs[[i]]); n <- sum(!is.na(x))
  z <- qnorm(seq(0,1,length.out=n+2)[2:(n+1)])
  new <- rep(NA_real_,length(x)); new[!is.na(x)] <- z[base::rank(x[!is.na(x)])]
  r <- rs[[i]]; raster::values(r) <- new; r
})
rs2 <- raster::stack(st2); names(rs2) <- names(rs)


# 4. SELECT FINAL COVARIATES ON FULL DATA -------------------------------------

pt_all <- data.frame(TB=TB0$`Total TB cases`,Nscreen=TB0$Number_examined,x=TB0$Longitude,y=TB0$Latitude)
pt_all <- pt_all[complete.cases(pt_all),]; xy_all <- cbind(pt_all$x,pt_all$y)
cov_all <- data.frame(raster::extract(rs2,xy_all))

source(file.path(path_input,"vif.R"))
keep.dat <- vif_func(in_frame=cov_all,thresh=5,trace=FALSE)
keep.dat <- setdiff(keep.dat,c("POP","ACH"))
rs2 <- rs2[[keep.dat]]

# Same WAIC forward selection as original analysis
cov_z_all <- data.frame(raster::extract(rs2,xy_all)); names(cov_z_all) <- keep.dat
myarea_all <- Africa[Africa$wb_a3 %in% unique(TB0$ISO3),]
bdry_all <- INLA::inla.sp2segment(myarea_all); bdry_all$loc <- INLA::inla.mesh.map(bdry_all$loc)
mesh_all <- INLA::inla.mesh.2d(loc=xy_all,boundary=bdry_all,max.edge=c(.6,4),offset=c(.6,4),cutoff=.6)
A_all <- INLA::inla.spde.make.A(mesh=mesh_all,loc=xy_all)
spde_all <- INLA::inla.spde2.pcmatern(mesh_all,prior.range=c(10,.1),prior.sigma=c(1,.1))

stk_all <- INLA::inla.stack(data=list(Y=pt_all$TB,n=pt_all$Nscreen),A=list(A_all,1),
  effects=list(list(z.field=seq_len(spde_all$n.spde)),list(covariate=cov_z_all)),tag="est.z")

waic <- sapply(seq_along(keep.dat),function(i) {
  f <- as.formula(paste("Y ~ f(z.field, model=spde_all) +",paste(keep.dat[1:i],collapse=" + ")))
  m <- INLA::inla(f,data=INLA::inla.stack.data(stk_all,spde=spde_all),family="binomial",Ntrials=pt_all$Nscreen,
    control.predictor=list(A=INLA::inla.stack.A(stk_all),compute=TRUE),control.compute=list(waic=TRUE),
    control.fixed=list(prec=1000,prec.intercept=.001),verbose=FALSE)
  m$waic$waic
})

final_covariates <- keep.dat[seq_len(which.min(waic))]
message("Final covariates held fixed in both models: ",paste(final_covariates,collapse=", "))


# 5. COMMON PREDICTION AND POPULATION GRIDS -----------------------------------

mask <- raster::aggregate(rs2[[1]],fact=mainaggfactor); raster::NAvalue(mask) <- -9999
pred_val_template <- raster::getValues(mask); w <- is.na(pred_val_template)
pred_locs <- raster::xyFromCell(mask,seq_len(raster::ncell(mask)))[!w,,drop=FALSE]
colnames(pred_locs) <- c("longitude","latitude")

covariates_pred <- data.frame(raster::extract(rs2[[final_covariates]],pred_locs))
names(covariates_pred) <- final_covariates

# Population density -> population per prediction cell
popdensity <- raster::aggregate(popden,fact=mainaggfactor)
cellarea <- raster::area(popdensity,na.rm=FALSE,weights=FALSE)
popsize <- popdensity * cellarea

linkfun <- function(x) mu/(1+exp(-x))


# 6. ADMINISTRATIVE BOUNDARIES -------------------------------------------------

countries <- unique(TB0$ISO3)

adm0poly <- do.call(rbind,lapply(countries,function(cc)
  readRDS(file.path(path_input,paste0("gadm36_",cc,"_0_sp.rds")))))

adm1list <- lapply(countries,function(cc) {
  f <- file.path(path_input,paste0("gadm36_",cc,"_1_sp.rds"))
  if (file.exists(f)) readRDS(f) else NULL
})
adm1poly <- do.call(rbind,adm1list[!vapply(adm1list,is.null,logical(1))])


# 7. FIT + PREDICT TB COUNTS ---------------------------------------------------

fit_predict <- function(TB,scenario) {

  message("Running: ",scenario," (n = ",nrow(TB),")")

  pt <- data.frame(TB=TB$`Total TB cases`,Nscreen=TB$Number_examined,x=TB$Longitude,y=TB$Latitude)
  pt <- pt[complete.cases(pt),]; xy <- cbind(pt$x,pt$y)
  cov_z <- data.frame(raster::extract(rs2[[final_covariates]],xy)); names(cov_z) <- final_covariates

  myarea <- Africa[Africa$wb_a3 %in% unique(TB$ISO3),]
  bdry <- INLA::inla.sp2segment(myarea); bdry$loc <- INLA::inla.mesh.map(bdry$loc)
  mesh <- INLA::inla.mesh.2d(loc=xy,boundary=bdry,max.edge=c(.6,4),offset=c(.6,4),cutoff=.6)
  A <- INLA::inla.spde.make.A(mesh=mesh,loc=xy)
  spde <- INLA::inla.spde2.pcmatern(mesh,prior.range=c(10,.1),prior.sigma=c(1,.1))

  formula.z <- as.formula(paste("Y ~ f(z.field, model=spde) +",paste(final_covariates,collapse=" + ")))
  stk <- INLA::inla.stack(data=list(Y=pt$TB,n=pt$Nscreen),A=list(A,1),
    effects=list(list(z.field=seq_len(spde$n.spde)),list(covariate=cov_z)),tag="est.z")

  model <- INLA::inla(formula.z,data=INLA::inla.stack.data(stk,spde=spde),family="binomial",Ntrials=pt$Nscreen,
    control.predictor=list(A=INLA::inla.stack.A(stk),compute=TRUE),control.compute=list(config=TRUE),
    control.fixed=list(prec=1000,prec.intercept=.001),verbose=FALSE)

  A.pred <- INLA::inla.spde.make.A(mesh=mesh,loc=pred_locs)
  set.seed(999); samp <- INLA::inla.posterior.sample(nn,model)
  pred <- matrix(NA_real_,nrow=nrow(A.pred),ncol=nn)

  for (i in seq_len(nn)) {
    latent <- samp[[i]]$latent; rn <- rownames(latent)
    field <- latent[grep("z.field",rn),]
    intercept <- latent[grep("\\(Intercept\\)",rn),][1]
    beta <- vapply(final_covariates,function(v) latent[grep(v,rn,fixed=TRUE),][1],numeric(1))
    lp <- intercept + as.vector(as.matrix(covariates_pred)%*%beta) + drop(A.pred%*%field)
    pred[,i] <- linkfun(lp)
  }

  pred_lower <- apply(pred,1,quantile,probs=.025,na.rm=TRUE)
  pred_mean <- rowMeans(pred,na.rm=TRUE)
  pred_upper <- apply(pred,1,quantile,probs=.975,na.rm=TRUE)

  make_raster <- function(v) {
    z <- pred_val_template; z[!w] <- v
    raster::setValues(mask,z)*prevunit
  }

  # Prevalence raster needed only internally to calculate counts
  TBprev <- raster::stack(make_raster(pred_lower),make_raster(pred_mean),make_raster(pred_upper))
  names(TBprev) <- c("lower","mean","upper")

  # Convert prevalence per 1,000 to estimated TB cases
  TBcount <- round(TBprev * popsize / prevunit)
  names(TBcount) <- c("lower","mean","upper")

  # Aggregate counts by SUM, not mean
  a0 <- raster::extract(TBcount,adm0poly,fun=sum,na.rm=TRUE,sp=TRUE)
  a1 <- raster::extract(TBcount,adm1poly,fun=sum,na.rm=TRUE,sp=TRUE)

  adm0 <- data.frame(GID_0=a0@data$GID_0,NAME_0=a0@data$NAME_0,
                     lower=a0@data$lower,count=a0@data$mean,upper=a0@data$upper)

  adm1 <- data.frame(GID_0=a1@data$GID_0,GID_1=a1@data$GID_1,NAME_1=a1@data$NAME_1,
                     lower=a1@data$lower,count=a1@data$mean,upper=a1@data$upper)

  rm(samp,pred,TBprev,TBcount); gc()
  list(adm0=adm0,adm1=adm1)
}


# 8. FIT BOTH SCENARIOS --------------------------------------------------------

result_all <- fit_predict(TB0,"All observations")
TB_noeth <- subset(TB0,ISO3!="ETH")
message("Ethiopian observations excluded: ",nrow(TB0)-nrow(TB_noeth))
result_noeth <- fit_predict(TB_noeth,"Ethiopia excluded")


# 9. ADM0 TB-COUNT COMPARISON --------------------------------------------------

adm0 <- merge(result_all$adm0,result_noeth$adm0,by="GID_0",suffixes=c("_all","_noeth")) %>%
  filter(GID_0!="ETH") %>%
  mutate(difference=count_noeth-count_all,percent_change=100*difference/count_all)

rmse0 <- sqrt(mean((adm0$count_all-adm0$count_noeth)^2,na.rm=TRUE))
pearson0 <- cor(adm0$count_all,adm0$count_noeth,use="complete.obs")
spearman0 <- cor(adm0$count_all,adm0$count_noeth,use="complete.obs",method="spearman")

write.csv(adm0,file.path(path_output,"csv","Ethiopia_exclusion_TBcount_ADM0.csv"),row.names=FALSE)
# Position correlation box on log-scale plot
x_box0 <- quantile(adm0$count_all, .999, na.rm=TRUE)
y_box0 <- quantile(adm0$count_noeth, .009, na.rm=TRUE)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

p0 <- ggplot(adm0,aes(count_all,count_noeth)) +
  geom_abline(intercept=0,slope=1,linetype="dashed") + geom_point(size=2) +
    geom_text_repel(aes(label=NAME_0_all), size=2.8, seed=999, box.padding=.8,
                  segment.alpha=.3, point.padding=.3, min.segment.length=0, max.overlaps=Inf) +
  scale_x_log10(labels=scales::label_comma()) + scale_y_log10(labels=scales::label_comma()) +
   annotate(
    "label",
       x=x_box0, y=y_box0,
    label=paste0(
      "Pearson r = ", round(pearson0,3),
      "\nSpearman rho = ", round(spearman0,3)
    ),
    hjust=1, vjust=0, size=3.5
  )  +
  labs(title="",
         x="Global ADM-0 estimated mean TB cases\nModel including Ethiopia data [log scale]",
       y="Global ADM-0 estimated mean TB cases\nModel excluding Ethiopia data [log scale]",
       caption="") +
  theme_bw()

ggsave(file.path(path_output,"pdf","Ethiopia_exclusion_TBcount_ADM0.pdf"),p0,width=6.5,height=5)


# 10. ADM1 TB-COUNT COMPARISON -------------------------------------------------

adm1 <- merge(result_all$adm1,result_noeth$adm1,by="GID_1",suffixes=c("_all","_noeth")) %>%
  filter(GID_0_all!="ETH") %>%
  mutate(difference=count_noeth-count_all,percent_change=100*difference/count_all)

rmse1 <- sqrt(mean((adm1$count_all-adm1$count_noeth)^2,na.rm=TRUE))
pearson1 <- cor(adm1$count_all,adm1$count_noeth,use="complete.obs")
spearman1 <- cor(adm1$count_all,adm1$count_noeth,use="complete.obs",method="spearman")

write.csv(adm1,file.path(path_output,"csv","Ethiopia_exclusion_TBcount_ADM1.csv"),row.names=FALSE)

# Position correlation box on log-scale plot
x_box1 <- quantile(adm1$count_all, .999, na.rm=TRUE)
y_box1 <- quantile(adm1$count_noeth, .009, na.rm=TRUE)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

p1 <- ggplot(adm1,aes(count_all,count_noeth)) +
  geom_abline(intercept=0,slope=1,linetype="dashed") + geom_point(size=.8,alpha=.4) +
  scale_x_log10(labels=scales::label_comma()) + scale_y_log10(labels=scales::label_comma()) +
   annotate(
    "label",
    x=x_box1, y=y_box1,
    label=paste0(
      "Pearson r = ", round(pearson1,3),
      "\nSpearman rho = ", round(spearman1,3)
    ),
    hjust=1, vjust=0, size=3.5
  )  +
  labs(title="",
       x="Global ADM-1 estimated mean TB cases\nModel including Ethiopia data [log scale]",
       y="Global ADM-1 estimated mean TB cases\nModel excluding Ethiopia data [log scale]",
       caption="") +
  theme_bw()

ggsave(file.path(path_output,"pdf","Ethiopia_exclusion_TBcount_ADM1.pdf"),p1,width=6.5,height=5)


# 11. SUMMARY -----------------------------------------------------------------

summary_results <- data.frame(
  level=c("ADM0","ADM1"),RMSE_cases=c(rmse0,rmse1),
  Pearson=c(pearson0,pearson1),Spearman=c(spearman0,spearman1),
  median_abs_percent_change=c(median(abs(adm0$percent_change),na.rm=TRUE),median(abs(adm1$percent_change),na.rm=TRUE)),
  max_abs_percent_change=c(max(abs(adm0$percent_change),na.rm=TRUE),max(abs(adm1$percent_change),na.rm=TRUE))
)

write.csv(summary_results,file.path(path_output,"csv","Eth_exc_TBcount_summary.csv"),row.names=FALSE)
print(summary_results)

adm0_ranked <- adm0 %>%
  select(Country=NAME_0_all,count_with_Ethiopia=count_all,count_without_Ethiopia=count_noeth,percent_change) %>%
  arrange(desc(abs(percent_change)))

write.csv(adm0_ranked,file.path(path_output,"csv","Eth_exc_TBcount_ADM0_ranked_changes.csv"),row.names=FALSE)
print(adm0_ranked)

message("Done. TB-count influence analysis of excluding data from Ethiopia saved in: ",path_output)
