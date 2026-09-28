# Prior and prediction-grid robustness analysis for 'Mapping TB in Africa'
# Runs FILTER or NOFILTER configuration
# Prior robustness: 8 prior combinations.

# Resolution robustness uses:
#   central estimate = posterior mean
#   uncertainty = 95% posterior UI (2.5%, 97.5%)
############################ USER SETTINGS ####################################
if (!exists("nn"))            nn <- 10000
if (!exists("mainaggfactor")) mainaggfactor <- 2
if (!exists("mu"))            mu <- 0.025
if (!exists("prevunit"))      prevunit <- 1000
if (!exists("popt"))          popt <- 5
# this analysis is carried out for the main spec (no population filter)
# but it can also be done with pop filter
popfilter <- FALSE 
set.seed(999)
# Prediction-grid resolution:
#   aggfactor 1  = 10 arcmin
#   aggfactor 2  = 20 arcmin
#   aggfactor 3  = 30 arcmin 
#   aggfactor 5  = 50 arcmin
#   aggfactor 10 = 100 arcmin
# range of aggregation to be assessed
aggfactor_grid <- c(1, 2, 3, 5, 10)# units in arcminutes

###############################################################################
# Packages
pkgs <- c("INLA", "raster", "rnaturalearth", "readxl", "ggplot2", "dplyr", "tidyr")
invisible(lapply(pkgs, require, character.only = TRUE))
# Paths
path_input <- file.path(getwd(), "INPUT")
path_output <- file.path(getwd(), "OUTPUT", "ALL", ifelse(popfilter, "FILTER", "NOFILTER"), "robustness")
dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path_output, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path_output, "pdf"), recursive = TRUE, showWarnings = FALSE)
# -----------------------------------------------------------------------------
# 1. DATA PREPARATION: ALL observations, NO population filter
# -----------------------------------------------------------------------------
TB <- readxl::read_excel(file.path(path_input, "TB_Africa_v2.xlsx"), trim_ws = TRUE)
TB$Number_examined <- round(as.numeric(TB$Number_examined))
TB$`Total TB cases` <- round(as.numeric(TB$`Total TB cases`))
TB$Latitude <- as.numeric(TB$Latitude)
TB$Longitude <- as.numeric(TB$Longitude)
TB <- TB[TB$`Lowest admin level (0-4)` > 0, ]
TB$Number_examined <- ifelse(TB$Number_examined < 0, 0, TB$Number_examined)
TB <- TB[complete.cases(TB[, "ISO3"]), ]
Africa <- rnaturalearth::ne_countries()
Africa <- Africa[Africa$continent == "Africa", ]
outc <- setdiff(unique(Africa$wb_a3), unique(TB$ISO3))
outc <- outc[complete.cases(outc)]
myarea <- Africa
for (i in seq_along(outc)) myarea <- subset(myarea, wb_a3 != outc[i])

# Covariates used in the original analysis
load(file.path(path_input, "originalcovariates.Rdata"))
rs <- raster::stack(T1, P1, alt, ahf, acc, popden)
# rs <- raster::crop(rs, raster::extent(myarea))
# rs <- raster::mask(rs, myarea)
names(rs) <- c("TMP", "PCP", "ALT", "ACH", "ACC", "POP")

# Quantile normalisation
st2 <- vector("list", raster::nlayers(rs))
for (i in seq_len(raster::nlayers(rs))) {
  x <- raster::getValues(rs[[i]])
  n <- sum(!is.na(x))
  z <- qnorm(seq(0, 1, length.out = n + 2)[2:(n + 1)])
  new <- rep(NA_real_, length(x))
  new[!is.na(x)] <- z[base::rank(x[!is.na(x)])]
  st2[[i]] <- rs[[i]]
  raster::values(st2[[i]]) <- new
}
rs2 <- raster::stack(st2)
names(rs2) <- names(rs)
pt <- data.frame(TB = TB$`Total TB cases`, Nscreen = TB$Number_examined, x = TB$Longitude, y = TB$Latitude)
pt <- pt[complete.cases(pt), ]
xy <- cbind(pt$x, pt$y)
covariate_z <- data.frame(raster::extract(rs2, xy))
source(file.path(path_input, "vif.R"))
keep.dat <- vif_func(in_frame = covariate_z, thresh = 5, trace = FALSE)
rs2 <- rs2[[keep.dat]]
covariate_z <- covariate_z[keep.dat]

# Same exclusions as the original analysis
covariate_z$POP <- NULL
covariate_z$ACH <- NULL
rs2 <- rs2[[names(covariate_z)]]
ptcov <- cbind(pt, covariate_z)
Y <- ptcov$TB
ntrials <- ptcov$Nscreen
# -----------------------------------------------------------------------------
# 2. MESH
# -----------------------------------------------------------------------------
bdry <- INLA::inla.sp2segment(myarea)
bdry$loc <- INLA::inla.mesh.map(bdry$loc)
mesh <- INLA::inla.mesh.2d(loc = xy, boundary = bdry, max.edge = c(0.6, 4), offset = c(0.6, 4), cutoff = 0.6)
A <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(xy))
# -----------------------------------------------------------------------------
# 3. ORIGINAL PREDICTION GRID AND POPULATION GRID
# -----------------------------------------------------------------------------
mask <- rs2[[1]]
raster::NAvalue(mask) <- -9999
mask <- raster::aggregate(mask, fact = mainaggfactor)
pred_val_template <- raster::getValues(mask)
w <- is.na(pred_val_template)
pred_locs <- raster::xyFromCell(mask, seq_len(raster::ncell(mask)))[!w, , drop = FALSE]
colnames(pred_locs) <- c("longitude", "latitude")
A.pred <- INLA::inla.spde.make.A(mesh = mesh, loc = pred_locs)
# Population for case-count calculations; NOFILTER only
popdensity <- raster::aggregate(popden, fact = mainaggfactor)

# if popfilter is TRUE, then apply population filter to the population grid
if(popfilter==TRUE) {
  popmask <- popdensity
  popmask[popmask < popt] <- NA
  popmask[!is.na(popmask)] <- 1
  popdensity <- popdensity * popmask
}

cellarea <- raster::area(popdensity, na.rm = FALSE, weights = FALSE)
popsize <- popdensity * cellarea
# Administrative boundaries needed for Tables 1 and 1b
mycountries <- unique(myarea$iso_a3)
adm0poly <- lapply(mycountries, function(cc) {
  readRDS(file.path(path_input, paste0("gadm36_", cc, "_0_sp.rds")))
})
adm0poly <- do.call(rbind, adm0poly)
# -----------------------------------------------------------------------------
# 4. INLA
# -----------------------------------------------------------------------------
formula.z <- as.formula(paste("Y ~ f(z.field, model = spde) +", paste(names(covariate_z), collapse = " + ")))
# Eight prior combinations
prior_grid <- expand.grid(cov_prec = c(1000, 10000), range0 = c(10, 5), sigma0 = c(1, 2), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  arrange(desc(cov_prec == 1000), desc(range0 == 10), desc(sigma0 == 1)) %>%
  mutate(prior_id = paste0("P", row_number()), scenario = sprintf("C%s_R%s_S%s", cov_prec, range0, sigma0), original = cov_prec == 1000 & range0 == 10 & sigma0 == 1)
prior_key <- prior_grid %>%
  transmute(Prior = prior_id, `Covariate precision (tau_beta)` = cov_prec, `rho: P(range < rho) = 0.1` = range0, `sigma: P(SD > sigma) = 0.1` = sigma0, Original = ifelse(original,
    "Yes", ""))
prior_note <- paste0("Prior specifications: ", paste0(prior_grid$prior_id, ": tau_beta=", prior_grid$cov_prec, ", rho=", prior_grid$range0, ", sigma=", prior_grid$sigma0, ifelse(prior_grid$original,
  " (original)", ""), collapse = "; "), "\n", "tau_beta = covariate prior precision; ", "rho denotes P(range < rho) = 0.1; ", "sigma denotes P(SD > sigma) = 0.1.")
write.csv(prior_grid, file.path(path_output, "csv", "prior_scenarios.csv"), row.names = FALSE)
write.csv(prior_key, file.path(path_output, "csv", "prior_specifications.csv"), row.names = FALSE)
all_table1 <- list()
all_table1b <- list()
all_table2 <- list()
all_global <- list()
linkfun <- function(x) mu/(1 + exp(-x))
# -----------------------------------------------------------------------------
# 5. FIT THE 8 PRIOR SCENARIOS
# -----------------------------------------------------------------------------
for (s in seq_len(nrow(prior_grid))) {
  pr <- prior_grid[s, ]
  message("Running prior scenario ", s, "/", nrow(prior_grid), ": ", pr$scenario)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(pr$range0, 0.1), prior.sigma = c(pr$sigma0, 0.1))
  stk <- INLA::inla.stack(data = list(Y = Y, n = ntrials), A = list(A, 1), effects = list(list(z.field = seq_len(spde$n.spde)), list(covariate = covariate_z)), tag = "est.z")
  mymodel <- INLA::inla(formula.z, data = INLA::inla.stack.data(stk, spde = spde), family = "binomial", Ntrials = ntrials, control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
    control.compute = list(config = TRUE), control.fixed = list(prec = pr$cov_prec, prec.intercept = 0.001), verbose = FALSE)
  #improve cpo computation (optional, takes about 10-15mn for 143 cases)
 if(mymodel$ok==FALSE){
 mymodel = inla.cpo(mymodel, force=TRUE)
 }
  # ---- Table 2: odds ratios and 95% CrI -------------------------------
  fixed <- mymodel$summary.fixed
  table2_num <- data.frame(Covariate = rownames(fixed), mean = exp(fixed$mean), lower = exp(fixed$`0.025quant`), upper = exp(fixed$`0.975quant`), scenario = pr$scenario, prior_id = pr$prior_id,
    original = pr$original, row.names = NULL)
  all_table2[[s]] <- table2_num
  # ---- Posterior prediction -------------------------------------------
  mycovnames <- mymodel$names.fixed[-1]
  covariates_pred <- data.frame(raster::extract(rs2[[mycovnames]], pred_locs))
  names(covariates_pred) <- mycovnames
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
  }
  # samp <- INLA::inla.posterior.sample(nn, mymodel)
  # pred <- matrix(NA_real_, nrow = nrow(A.pred), ncol = nn)
  # for (i in seq_len(nn)) {
  #   latent <- samp[[i]]$latent
  #   rn <- rownames(latent)
  #   field <- latent[grep("z.field", rn), ]
  #   intercept <- latent[grep("\\(Intercept\\)", rn), ][1]
  #   beta <- vapply(mycovnames, function(v) latent[grep(v, rn, fixed = TRUE), ][1], numeric(1))
  #   lp <- intercept + as.vector(as.matrix(covariates_pred) %*% beta) + drop(A.pred %*% field)
  #   pred[, i] <- linkfun(lp)
  # }
  pred_lower <- apply(pred, 1, quantile, probs = 0.025, na.rm = TRUE)
  pred_median <- apply(pred, 1, median, na.rm = TRUE)
  pred_q25 <- apply(pred, 1, quantile, probs = 0.25, na.rm = TRUE)
  pred_q75 <- apply(pred, 1, quantile, probs = 0.75, na.rm = TRUE)
  pred_iqr <- pred_q75 - pred_q25
  pred_mean <- rowMeans(pred, na.rm = TRUE)
  pred_upper <- apply(pred, 1, quantile, probs = 0.975, na.rm = TRUE)
  make_raster <- function(x) {
    z <- pred_val_template
    z[!w] <- x
    raster::setValues(mask, z) * prevunit
  }
  if(popfilter == TRUE) TBprev <- raster::stack(
  make_raster(pred_lower) * popmask,
  make_raster(pred_median) * popmask,
  make_raster(pred_q25) * popmask,
  make_raster(pred_q75) * popmask,
  make_raster(pred_iqr) * popmask,
  make_raster(pred_mean) * popmask,
  make_raster(pred_upper) * popmask
) else TBprev <- raster::stack(
  make_raster(pred_lower),
  make_raster(pred_median),
  make_raster(pred_q25),
  make_raster(pred_q75),
  make_raster(pred_iqr),
  make_raster(pred_mean),
  make_raster(pred_upper)
)
names(TBprev) <- c("LITBprev", "medianTBprev", "lowqTBprev", "uppqTBprev", "IQRTBprev", "meanTBprev", "UITBprev")
TBcount <- round(TBprev * popsize / prevunit)
names(TBcount) <- c("LITBcount", "medianTBcount", "lowqTBcount", "uppqTBcount", "IQRTBcount", "meanTBcount", "UITBcount")
  
  # ---- Global estimates -----------------------------------------------
  casesest <- data.frame(raster::cellStats(TBcount, sum))
  global_num <- data.frame(Area = "Global", mean = casesest["meanTBcount", 1], lower = casesest["LITBcount", 1], upper = casesest["UITBcount", 1], median = casesest["medianTBcount",
    1], IQR = casesest["IQRTBcount", 1], scenario = pr$scenario, prior_id = pr$prior_id, original = pr$original)
  all_global[[s]] <- global_num
 
  # ---- Country Tables 1 and 1b ---------------------------------------
  Af0TBprev <- raster::extract(TBprev, adm0poly, fun = median, na.rm = TRUE, sp = TRUE)
  Af0TB <- raster::extract(TBcount, adm0poly, fun = sum, na.rm = TRUE, sp = TRUE)
  table1_num <- data.frame(Country = Af0TB@data$NAME_0, mean = Af0TB@data$meanTBcount, lower = Af0TB@data$LITBcount, upper = Af0TB@data$UITBcount, median = Af0TB@data$medianTBcount,
    IQR = Af0TB@data$IQRTBcount, scenario = pr$scenario, prior_id = pr$prior_id, original = pr$original)
  all_table1[[s]] <- table1_num
  table1b_num <- data.frame(Country = Af0TBprev@data$NAME_0, mean = Af0TBprev@data$meanTBprev, lower = Af0TBprev@data$LITBprev, upper = Af0TBprev@data$UITBprev, median = Af0TBprev@data$medianTBprev,
    IQR = Af0TBprev@data$IQRTBprev, scenario = pr$scenario, prior_id = pr$prior_id, original = pr$original)
  all_table1b[[s]] <- table1b_num
  rm(mymodel, samp, pred, TBprev, TBcount)
  gc()
}
# -----------------------------------------------------------------------------
# 6. SAVE PRIOR ROBUSTNESS RESULTS
# -----------------------------------------------------------------------------
table1_all <- dplyr::bind_rows(all_table1)
table1b_all <- dplyr::bind_rows(all_table1b)
table2_all <- dplyr::bind_rows(all_table2)
global_all <- dplyr::bind_rows(all_global)
write.csv(table1_all, file.path(path_output, "csv", "rob_Table1.csv"), row.names = FALSE)
write.csv(table1b_all, file.path(path_output, "csv", "rob_Table1b.csv"), row.names = FALSE)
write.csv(table2_all, file.path(path_output, "csv", "rob_Table2.csv"), row.names = FALSE)
write.csv(global_all, file.path(path_output, "csv", "rob_globalestimates.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------
# 7. PRIOR ROBUSTNESS PLOTS
# -----------------------------------------------------------------------------
plot_country_robustness <- function(dat, title, xlab, filename, log_scale = FALSE) {
  ref <- dat %>%
    filter(original) %>%
    select(Country, ref = mean)
  dd <- dat %>%
    left_join(ref, by = "Country")
  p <- ggplot(dd, aes(x = mean, y = prior_id)) + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) + 
  geom_point(aes(shape = original), size = 2) + 
     scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Mean value using model with\noriginal prior", "FALSE" = "Mean value using model with\nalternative priors"),
    name = ""
  )+
  geom_vline(aes(xintercept = ref),
    linetype = "dashed", linewidth = 0.3) + facet_wrap(~Country, scales = "free_x") + 
    labs(title = title, x = xlab, y = "Prior specification", 
    caption = "") + theme_bw() +  guides(shape=guide_legend(nrow=1, byrow=TRUE)) +
 theme(
  strip.background=element_blank(),
  strip.text.x=element_text( size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=6),
   axis.text.y=element_text(size=6),
  plot.caption=element_text(hjust=0, size=7.5)
)
    #caption = if (log_scale)
    #paste0(prior_note, " X-axis is log10 transformed.") else prior_note) + theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 7.5))
  if (log_scale)
    p <- p + scale_x_log10(labels=scales::label_comma())
  ggsave(file.path(path_output, "pdf", filename), p, width = 11, height = 9)
}
plot_country_robustness(table1_all, "", "ADM-0 estimated mean TB cases (95% UI), log10 scale",
 "rob_Table1.pdf", log_scale = TRUE)
plot_country_robustness(table1b_all, "", "ADM-0 estimated mean TB prevalence per 1,000 (95% UI)", 
"rob_Table1b.pdf")

ref2 <- table2_all %>%
  filter(original) %>%
  select(Covariate, ref = mean)
t2plot <- table2_all %>%
  left_join(ref2, by = "Covariate")
p2 <- ggplot(t2plot, aes(x = mean, y = prior_id)) + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) + 
geom_point(aes(shape = original), size = 2) + 
     scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Mean value using model with\noriginal prior", "FALSE" = "Mean value using model with\nalternative priors"),
    name = ""
  )+
geom_vline(aes(xintercept = ref),
  linetype = "dashed", linewidth = 0.3) + geom_vline(xintercept = 1, linetype = "dotted") + facet_wrap(~Covariate, scales = "free_x") + 
  labs(title = "",
  x = "Estimated odds ratio (95% UI)", y = "Prior specification", shape = "Original priors") + 
  theme_bw() + theme(
  plot.caption = element_text(hjust = 0, size = 7.5))+
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text( size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=8),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)

ggsave(file.path(path_output, "pdf", "rob_Table2.pdf"), p2, width = 8, height = 6)

refg <- global_all %>%
  filter(original) %>%
  pull(mean)
#Plot
pg <- ggplot(global_all, aes(x = mean, y = prior_id)) + 
geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15) + 
geom_point(aes(shape = original), size = 2) + 
     scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Mean value using model with\noriginal prior", "FALSE" = "Mean value using model with\nalternative priors"),
    name = ""
  )+ geom_vline(xintercept = refg,
  linetype = "dashed", linewidth = 0.4) + scale_x_log10(labels = scales::label_comma()) + 
  labs(title = "", x = "Global estimated mean TB cases (95% UI; log10 x-axis)",
  y = "Prior specification", shape = "Original priors") + 
  theme_bw() +
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text(face="bold", size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=8),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)
ggsave(file.path(path_output, "pdf", "rob_globalestimates.pdf"), pg, width = 5, height = 6)
# -----------------------------------------------------------------------------
# 8. ROBUSTNESS TO PREDICTION-GRID RESOLUTION
# Main specification only:
# aggfactor = 1,2, 3, 5, 10
#
# Central estimate = posterior mean
# Uncertainty = 95% UI
# -----------------------------------------------------------------------------
resolution_results <- list()
resolution_adm0_results <- list()
# ---- Fit main specification once -------------------------------------------
spde_main <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(10, 0.1), prior.sigma = c(1, 0.1))
stk_main <- INLA::inla.stack(data = list(Y = Y, n = ntrials), A = list(A, 1), effects = list(list(z.field = seq_len(spde_main$n.spde)), list(covariate = covariate_z)), tag = "est.z")
model_main <- INLA::inla(formula.z, data = INLA::inla.stack.data(stk_main, spde = spde_main), family = "binomial", Ntrials = ntrials, control.predictor = list(A = INLA::inla.stack.A(stk_main),
  compute = TRUE), control.compute = list(config = TRUE), control.fixed = list(prec = 1000, prec.intercept = 0.001), verbose = FALSE)
# Same posterior samples for all resolutions
set.seed(999)
samp_main <- INLA::inla.posterior.sample(nn, model_main)
mycovnames <- model_main$names.fixed[-1]
# ---- Loop over resolutions --------------------------------------------------
for (aggfactor_res in aggfactor_grid) {
  message("Running resolution robustness: aggfactor = ", aggfactor_res)
  # ---- Prediction grid ------------------------------------------------------
  mask_res <- rs2[[1]]
  raster::NAvalue(mask_res) <- -9999
  mask_res <- raster::aggregate(mask_res, fact = aggfactor_res)
  pred_val_template_res <- raster::getValues(mask_res)
  w_res <- is.na(pred_val_template_res)
  pred_locs_res <- raster::xyFromCell(mask_res, seq_len(raster::ncell(mask_res)))[!w_res, , drop = FALSE]
  colnames(pred_locs_res) <- c("longitude", "latitude")
  A.pred_res <- INLA::inla.spde.make.A(mesh = mesh, loc = pred_locs_res)
  # ---- Population grid ------------------------------------------------------
  # popdensity_res <- raster::aggregate(popden, fact = aggfactor_res)
  # cellarea_res <- raster::area(popdensity_res, na.rm = FALSE, weights = FALSE)
  # popsize_res <- popdensity_res * cellarea_res

  # ---- Population grid ------------------------------------------------------
  popdensity_res <- raster::aggregate(popden, fact = aggfactor_res)
  popmask_res <- popdensity_res 
  if(popfilter == TRUE) {
    popmask_res[popmask_res < popt] <- NA
    popmask_res[!is.na(popmask_res)] <- 1
    popdensity_res <- popdensity_res * popmask_res
  }
  cellarea_res <- raster::area(popdensity_res, na.rm = FALSE, weights = FALSE)
  popsize_res <- popdensity_res * cellarea_res

  # ---- Covariates at prediction locations -----------------------------------
  covariates_pred_res <- data.frame(raster::extract(rs2[[mycovnames]], pred_locs_res))
  names(covariates_pred_res) <- mycovnames
  # ---- Posterior prediction -------------------------------------------------
  pred_res <- matrix(NA_real_, nrow = nrow(A.pred_res), ncol = nn)
  for (i in seq_len(nn)) {
    latent <- samp_main[[i]]$latent
    rn <- rownames(latent)
    field <- latent[grep("z.field", rn), ]
    intercept <- latent[grep("\\(Intercept\\)", rn), ][1]
    beta <- vapply(mycovnames, function(v) latent[grep(v, rn, fixed = TRUE), ][1], numeric(1))
    lp <- intercept + as.vector(as.matrix(covariates_pred_res) %*% beta) + drop(A.pred_res %*% field)
    pred_res[, i] <- linkfun(lp)
  }
  # ---- Posterior summaries --------------------------------------------------
  # Mean is the central estimate; UI remains the 2.5% and 97.5% quantiles.
  pred_lower_res <- apply(pred_res, 1, quantile, probs = 0.025, na.rm = TRUE)
  pred_mean_res <- rowMeans(pred_res, na.rm = TRUE)
  pred_upper_res <- apply(pred_res, 1, quantile, probs = 0.975, na.rm = TRUE)
  # ---- Convert prevalence to raster ----------------------------------------
  make_raster_res <- function(x) {
    z <- pred_val_template_res
    z[!w_res] <- x
    raster::setValues(mask_res, z) * prevunit
  }
  TBprev_res <- raster::stack(make_raster_res(pred_lower_res),  make_raster_res(pred_mean_res), make_raster_res(pred_upper_res) )
  if(popfilter == TRUE) {
    TBprev_res <- TBprev_res * popmask_res
  } 
  names(TBprev_res) <- c("lower", "mean", "upper")
  # ---- Estimated TB counts --------------------------------------------------
  TBcount_res <- round(TBprev_res * popsize_res/prevunit)
  names(TBcount_res) <- c("lower", "mean", "upper")
  # ---- Global estimates -----------------------------------------------------
  global_count_res <- raster::cellStats(TBcount_res, sum, na.rm = TRUE)
  global_population_res <- raster::cellStats(popsize_res, sum, na.rm = TRUE)
  global_prev_res <- global_count_res/global_population_res * prevunit
  # ---- ADM0 estimates -------------------------------------------------------
  adm0_prev_res <- raster::extract(TBprev_res, adm0poly, fun = median, na.rm = TRUE, sp = TRUE)
  adm0_count_res <- raster::extract(TBcount_res, adm0poly, fun = sum, na.rm = TRUE, sp = TRUE)
  # ---- Store global results -------------------------------------------------
  resolution_results[[as.character(aggfactor_res)]] <- data.frame(aggfactor = aggfactor_res, arcmin = aggfactor_res * 10, global_prev_lower = as.numeric(global_prev_res["lower"]),
    global_prev_mean = as.numeric(global_prev_res["mean"]), global_prev_upper = as.numeric(global_prev_res["upper"]), global_count_lower = as.numeric(global_count_res["lower"]),
    global_count_mean = as.numeric(global_count_res["mean"]), global_count_upper = as.numeric(global_count_res["upper"]))
  # ---- Store ADM0 results ---------------------------------------------------
  resolution_adm0_results[[as.character(aggfactor_res)]] <- data.frame(Country = adm0_prev_res@data$NAME_0, aggfactor = aggfactor_res, arcmin = aggfactor_res * 10, prevalence_lower = adm0_prev_res@data$lower,
    prevalence = adm0_prev_res@data$mean, prevalence_upper = adm0_prev_res@data$upper, count_lower = adm0_count_res@data$lower, count = adm0_count_res@data$mean, count_upper = adm0_count_res@data$upper)
  rm(mask_res, pred_locs_res, A.pred_res, pred_res, TBprev_res, TBcount_res)
  gc()
}
# -----------------------------------------------------------------------------
# 9. COMBINE AND COMPARE RESOLUTION RESULTS
# -----------------------------------------------------------------------------
resolution_global <- dplyr::bind_rows(resolution_results)
resolution_adm0 <- dplyr::bind_rows(resolution_adm0_results)
# Main specification: aggfactor = mainaggfactor
ref_global <- resolution_global %>%
  filter(aggfactor == mainaggfactor) %>%
  select(ref_prevalence = global_prev_mean, ref_count = global_count_mean)
resolution_global <- resolution_global %>%
  mutate(ref_prevalence = ref_global$ref_prevalence, ref_count = ref_global$ref_count, prevalence_percent_change = (global_prev_mean - ref_prevalence)/ref_prevalence * 100, count_percent_change = (global_count_mean -
    ref_count)/ref_count * 100)
ref_adm0 <- resolution_adm0 %>%
  filter(aggfactor == mainaggfactor) %>%
  select(Country, ref_prevalence = prevalence, ref_count = count)
resolution_adm0 <- resolution_adm0 %>%
  left_join(ref_adm0, by = "Country") %>%
  mutate(prevalence_percent_change = (prevalence - ref_prevalence)/ref_prevalence * 100, count_percent_change = (count - ref_count)/ref_count * 100)
# -----------------------------------------------------------------------------
# 10. SAVE RESOLUTION ROBUSTNESS TABLES
# -----------------------------------------------------------------------------
write.csv(resolution_global, file.path(path_output, "csv", "rob_resolution_global.csv"), row.names = FALSE)
write.csv(resolution_adm0, file.path(path_output, "csv", "rob_resolution_ADM0.csv"), row.names = FALSE)
resolution_summary <- resolution_global %>%
  select(aggfactor, arcmin, global_prev_lower, global_prev_mean, global_prev_upper, global_count_lower, global_count_mean, global_count_upper, prevalence_percent_change, count_percent_change)
write.csv(resolution_summary, file.path(path_output, "csv", "rob_resolution_summary.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------
# 11. RESOLUTION ROBUSTNESS PLOTS
# -----------------------------------------------------------------------------
resolution_adm0$arcmin_label <- factor(resolution_adm0$aggfactor, levels = as.factor(aggfactor_grid), labels = paste0(aggfactor_grid * 10, " arcmin"))
resolution_global$arcmin_label <- factor(resolution_global$aggfactor, levels = as.factor(aggfactor_grid)  , labels = paste0(aggfactor_grid * 10, " arcmin"))
# ADM0 prevalence 
p_res_prev <- ggplot(resolution_adm0, aes(x = prevalence, y = arcmin_label)) + geom_errorbarh(aes(xmin = prevalence_lower, xmax = prevalence_upper), height = 0.15) + geom_point(size = 2) +
  geom_vline(data = resolution_adm0 %>%
    filter(aggfactor == mainaggfactor), aes(xintercept = prevalence), linetype = "dashed", linewidth = 0.4) + facet_wrap(~Country, scales = "free_x") + labs(title = "",
  x = "Estimated TB prevalence per 1,000 (95% UI)", y = "Prediction-grid resolution", caption = "Points and horizontal lines show posterior means and 95% uncertainty intervals. Dashed line = main specification.") +
  theme_bw() +
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text(size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=8),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)
ggsave(file.path(path_output, "pdf", "rob_res_ADM0_prev.pdf"), p_res_prev, width = 9, height = 8)
# ADM0 TB counts 
p_res_count <- ggplot(resolution_adm0, aes(x = count, y = arcmin_label)) + geom_errorbarh(aes(xmin = count_lower, xmax = count_upper), height = 0.15) + geom_point(size = 2) + 
geom_vline(data = resolution_adm0 %>%
  filter(aggfactor == mainaggfactor), aes(xintercept = count), linetype = "dashed", linewidth = 0.4) + facet_wrap(~Country, scales = "free_x") + scale_x_log10(labels = scales::label_comma()) +
  labs(title = "", x = "Estimated TB cases (95% UI; log10 x-axis)", y = "Prediction-grid resolution", caption = "Points and horizontal lines show posterior means and 95% uncertainty intervals. Dashed line = main specification. X-axis is log10 transformed.") +
   theme_bw() +
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text(size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=7),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)
ggsave(file.path(path_output, "pdf", "rob_res_ADM0_count.pdf"), p_res_count, width = 10, height = 8)
# Global prevalence 
global_ref_prev <- resolution_global %>%
  filter(aggfactor == mainaggfactor) %>%
  pull(global_prev_mean)
p_res_global_prev <- ggplot(resolution_global, aes(x = global_prev_mean, y = arcmin_label)) + geom_errorbarh(aes(xmin = global_prev_lower, xmax = global_prev_upper), height = 0.15) +
  geom_point(size = 2.5) + geom_vline(xintercept = global_ref_prev, linetype = "dashed", linewidth = 0.4) + 
  labs(title = "",
  x = "Global estimated TB prevalence per 1,000 (95% UI)", y = "Prediction-grid resolution", 
  caption = "Points and horizontal lines show posterior means and 95% uncertainty intervals.\nDashed line = main specification.") +
   theme_bw() +
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text(size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=8),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)
ggsave(file.path(path_output, "pdf", "rob_res_global_prev.pdf"), p_res_global_prev, width = 6, height = 4)
# ---- Global TB counts --------------------------------------------------------
global_ref_count <- resolution_global %>%
  filter(aggfactor == mainaggfactor) %>%
  pull(global_count_mean)
p_res_global_count <- ggplot(resolution_global, aes(x = global_count_mean, y = arcmin_label)) + geom_errorbarh(aes(xmin = global_count_lower, xmax = global_count_upper), height = 0.15) +
  geom_point(size = 2.5) + geom_vline(xintercept = global_ref_count, linetype = "dashed", linewidth = 0.4) + scale_x_log10(labels = scales::label_comma()) + 
  labs(title = "",
  x = "Global estimated TB cases (95% UI; log10 x-axis)", y = "Prediction-grid resolution", 
  caption = "Points and horizontal lines show posterior means and 95% uncertainty intervals\nDashed line = main specification. X-axis is log10 transformed.") +
  theme_bw() +
   theme(
  strip.background=element_blank(),
  strip.text.x=element_text(size=10, hjust=0),
  strip.text = element_text(margin=margin(l=0)),
  legend.position="bottom",
  legend.direction="horizontal",
  axis.text.x=element_text(size=8),
   axis.text.y=element_text(size=8),
  plot.caption=element_text(hjust=0, size=7.5)
)
ggsave(file.path(path_output, "pdf", "rob_res_global_count.pdf"), p_res_global_count, width = 5, height = 4)
# -----------------------------------------------------------------------------
# 12. FINAL
# -----------------------------------------------------------------------------
message("Done. Robustness outputs are in: ", path_output)
# End of the code

