=============================================================================
# Analysis:
# Comparing reported results with re-analysis at ADM-0
# =============================================================================
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
###############################################################################
# Packages
###############################################################################
pkgs <- c( "INLA", "raster", "rnaturalearth", "readxl", "ggplot2", "dplyr","ggrepel" )
invisible(lapply(pkgs, require, character.only = TRUE))
###############################################################################
# Paths
###############################################################################
path_input <- file.path(getwd(), "INPUT")
path_output <- file.path(getwd(), "OUTPUT", "ALL", ifelse(popfilter, "FILTER", "NOFILTER"), "reported_vs_reanalysis")
dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path_output, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path_output, "pdf"), recursive = TRUE, showWarnings = FALSE)
###############################################################################
# 1. DATA PREPARATION
###############################################################################
TB <- readxl::read_excel( file.path(path_input, "TB_Africa_v2.xlsx"), trim_ws = TRUE )
TB$Number_examined <- round(as.numeric(TB$Number_examined))
TB$`Total TB cases` <- round(as.numeric(TB$`Total TB cases`))
TB$Latitude <- as.numeric(TB$Latitude)
TB$Longitude <- as.numeric(TB$Longitude)
TB <- TB[TB$`Lowest admin level (0-4)` > 0, ]
TB$Number_examined <- ifelse( TB$Number_examined < 0, 0, TB$Number_examined )
TB <- TB[complete.cases(TB[, "ISO3"]), ]
###############################################################################
# Africa map
###############################################################################
Africa <- rnaturalearth::ne_countries()
Africa <- Africa[Africa$continent == "Africa", ]
outc <- setdiff( unique(Africa$wb_a3), unique(TB$ISO3) )
outc <- outc[complete.cases(outc)]
myarea <- Africa
for (i in seq_along(outc)) {
  myarea <- subset(myarea, wb_a3 != outc[i])
}
###############################################################################
# Covariates
###############################################################################
load(file.path(path_input, "originalcovariates.Rdata"))
rs <- raster::stack( T1, P1, alt, ahf, acc, popden )
names(rs) <- c( "TMP", "PCP", "ALT", "ACH", "ACC", "POP" )
###############################################################################
# Quantile normalisation
###############################################################################
st2 <- vector( "list", raster::nlayers(rs) )
for (i in seq_len(raster::nlayers(rs))) {
  x <- raster::getValues(rs[[i]])
  n <- sum(!is.na(x))
  z <- qnorm( seq(0, 1, length.out = n + 2)[2:(n + 1)] )
  new <- rep(NA_real_, length(x))
  new[!is.na(x)] <-
    z[base::rank(x[!is.na(x)])]
  st2[[i]] <- rs[[i]]
  raster::values(st2[[i]]) <- new
}
rs2 <- raster::stack(st2)
names(rs2) <- names(rs)
###############################################################################
# Observation locations
###############################################################################
pt <- data.frame( TB = TB$`Total TB cases`, Nscreen = TB$Number_examined, x = TB$Longitude, y = TB$Latitude )
pt <- pt[complete.cases(pt), ]
xy <- cbind(pt$x, pt$y)
covariate_z <- data.frame( raster::extract(rs2, xy) )
source(file.path(path_input, "vif.R"))
keep.dat <- vif_func( in_frame = covariate_z, thresh = 5, trace = FALSE )
rs2 <- rs2[[keep.dat]]
covariate_z <- covariate_z[keep.dat]
# Same exclusions as the original analysis
covariate_z$POP <- NULL
covariate_z$ACH <- NULL
rs2 <- rs2[[names(covariate_z)]]
ptcov <- cbind( pt, covariate_z )
Y <- ptcov$TB
ntrials <- ptcov$Nscreen
###############################################################################
# 2. INLA MESH
###############################################################################
bdry <- INLA::inla.sp2segment(myarea)
bdry$loc <- INLA::inla.mesh.map( bdry$loc )
mesh <- INLA::inla.mesh.2d( loc = xy, boundary = bdry, max.edge = c(0.6, 4), offset = c(0.6, 4), cutoff = 0.6 )
A <- INLA::inla.spde.make.A( mesh = mesh, loc = as.matrix(xy) )
###############################################################################
# 3. PREDICTION GRID
###############################################################################
mask <- rs2[[1]]
raster::NAvalue(mask) <- -9999
mask <- raster::aggregate( mask, fact = mainaggfactor )
pred_val_template <- raster::getValues(mask)
w <- is.na(pred_val_template)
pred_locs <- raster::xyFromCell( mask, seq_len(raster::ncell(mask)) )[!w, , drop = FALSE]
colnames(pred_locs) <- c( "longitude", "latitude" )
A.pred <- INLA::inla.spde.make.A( mesh = mesh, loc = pred_locs )
###############################################################################
# 4. POPULATION GRID
###############################################################################
popdensity <- raster::aggregate( popden, fact = mainaggfactor )
# if popfilter is TRUE, then apply population filter to the population grid
if(popfilter==TRUE) {
  popmask <- popdensity
  popmask[popmask < popt] <- NA
  popmask[!is.na(popmask)] <- 1
  popdensity <- popdensity * popmask
}

cellarea <- raster::area(popdensity, na.rm = FALSE, weights = FALSE)
popsize <- popdensity * cellarea
###############################################################################
# 5. ADMINISTRATIVE POLYGONS
###############################################################################
mycountries <- unique(myarea$iso_a3)
# ADM0
adm0poly <- lapply(
  mycountries,
  function(cc) {
    readRDS( file.path( path_input, paste0("gadm36_", cc, "_0_sp.rds") ) )
  }
)
adm0poly <- do.call( rbind, adm0poly )
###############################################################################
# 6. ORIGINAL MODEL ONLY
###############################################################################
spde <- INLA::inla.spde2.pcmatern( mesh, prior.range = c(10, 0.1), prior.sigma = c(1, 0.1) )
formula.z <- as.formula( paste( "Y ~ f(z.field, model = spde) +", paste( names(covariate_z), collapse = " + " ) ) )
stk <- INLA::inla.stack( data = list( Y = Y, n = ntrials ), A = list( A, 1 ), effects = list( list( z.field = seq_len(spde$n.spde) ), list( covariate = covariate_z ) ), tag = "est.z" )
mymodel <- INLA::inla( formula.z, data = INLA::inla.stack.data( stk, spde = spde ), family = "binomial", Ntrials = n, control.predictor = list( A = INLA::inla.stack.A(stk), compute = TRUE ), 
control.compute = list( config = TRUE ), control.fixed = list( prec = 1000, prec.intercept = 0.001 ), verbose = FALSE )
###############################################################################
# 7. POSTERIOR PREDICTION
###############################################################################
mycovnames <- mymodel$names.fixed[-1]
covariates_pred <- data.frame( raster::extract( rs2[[mycovnames]], pred_locs ) )
names(covariates_pred) <- mycovnames
set.seed(999)
samp <- INLA::inla.posterior.sample( nn, mymodel )
pred <- matrix( NA_real_, nrow = nrow(A.pred), ncol = nn )
linkfun <- function(x) {
  mu / (1 + exp(-x))
}
for (i in seq_len(nn)) {
  latent <- samp[[i]]$latent
  rn <- rownames(latent)
  field <- latent[
    grep("z.field", rn),
    ]
  intercept <- latent[
    grep("\\(Intercept\\)", rn),
    ][1]
  beta <- vapply(
    mycovnames,
    function(v) {
      latent[
        grep( v, rn, fixed = TRUE ),
      ][1]
    },
    numeric(1)
  )
  lp <-
    intercept +
    as.vector( as.matrix(covariates_pred) %*% beta ) +
    drop( A.pred %*% field )
  pred[, i] <- linkfun(lp)
}
###############################################################################
# predictions 
###############################################################################
pred_2.5 <- apply(pred, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_med = apply(pred, 1, function(x) quantile(x, probs=c(0.5), na.rm=TRUE))
pred_sd = apply(pred, 1, sd)
pred_25pct = apply(pred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_75pct = apply(pred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
IQR = pred_75pct - pred_25pct
pred_mean <- apply(pred, 1, function(x) mean(x, na.rm = TRUE))
pred_975 <- apply(pred, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
###############################################################################
# Convert vector to raster
###############################################################################
make_raster <- function(x) {
  z <- pred_val_template
  z[!w] <- x
  raster::setValues( mask, z ) * prevunit
}
TBprev_mean <- make_raster( pred_mean )
 if(popfilter == TRUE) {
  TBprev_mean <- TBprev_mean * popmask
}
names(TBprev_mean) <- "TB_prev_mean"

TBprev_median <- make_raster( pred_median )
 if(popfilter == TRUE) {
  TBprev_median <- TBprev_median * popmask
}
names(TBprev_median) <- "TB_prev_median"

# Lower 95% UI
TBprev_lower <- make_raster(pred_lower)
if(popfilter == TRUE) TBprev_lower <- TBprev_lower * popmask
names(TBprev_lower) <- "TB_prev_lower"

# Upper 95% UI
TBprev_upper <- make_raster(pred_upper)
if(popfilter == TRUE) TBprev_upper <- TBprev_upper * popmask
names(TBprev_upper) <- "TB_prev_upper"

# Posterior IQR
TBprev_IQR <- make_raster(IQR)
if(popfilter == TRUE) TBprev_IQR <- TBprev_IQR * popmask
names(TBprev_IQR) <- "TB_prev_IQR"

TBprev_updated <- raster::stack(TBprev_mean, TBprev_lower, TBprev_upper)
TBcount_updated <- round( TBprev_updated * popsize / prevunit )
names(TBcount_updated) <- c('TB_count_mean','TB_count_lower','TB_count_upper')

###############################################################################
# adm0
adm0poly <- adm0poly[, c("GID_0", "NAME_0")]
adm0_count_mean <- raster::extract( TBcount_updated[['TB_count_mean']], adm0poly, fun = sum, na.rm = TRUE, sp = TRUE )

# Create ADM0 comparison data frame
adm0_compare <- data.frame(
  GID_0 = adm0_count_mean@data$GID_0,
  Country = adm0_count_mean@data$NAME_0,
  count_mean = adm0_count_mean@data$TB_count_mean
)

#comparison on what was reported in the paper with updated values
#values from Table 1 (paper, median)
reported_count <- c(
  "Nigeria" = 460247,
  "Mozambique" = 120622,
  "Ghana" = 111828,
  "Kenya" = 97763,
  "South Africa" = 88147,
  "Uganda" = 74571,
  "Ethiopia" = 68838,
  "Cameroon" = 63127,
  "Zimbabwe" = 45031,
  "Chad" = 41320,
  "Zambia" = 24726,
  "Mauritania" = 11603,
  "Rwanda" = 2207,
  "Guinea-Bissau" = 1952
)

adm0_compare$reported_count <- reported_count[adm0_compare$Country]
# ---------------------------------------------------------------------------
# compared reported values with updated values
# ---------------------------------------------------------------------------
# Correlations
cor_spearman <- cor(
  adm0_compare$reported_count,
  adm0_compare$count_mean,
  use="complete.obs",
  method="spearman"
)

cor_pearson <- cor(
  adm0_compare$reported_count,
  adm0_compare$count_mean,
  use="complete.obs",
  method="pearson"
)

x_box <- quantile(adm0_compare$count_mean, .99, na.rm=TRUE)
y_box <- quantile(adm0_compare$reported_count, .01, na.rm=TRUE)

#plot ADM-0
p1 <- ggplot(adm0_compare, aes(x=count_mean, y=reported_count)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  geom_text_repel(aes(label=Country), size=2.8, seed=999, box.padding=.8,
                  segment.alpha=.3, point.padding=.3, min.segment.length=0, max.overlaps=Inf) +
  annotate("label", x=x_box, y=y_box,
           label=paste0("Pearson r = ", round(cor_pearson,3),
                        "\nSpearman rho = ", round(cor_spearman,3)),
           hjust=1, vjust=0, size=3) +
  scale_x_log10(labels=scales::label_comma()) +
  scale_y_log10(labels=scales::label_comma()) +
  labs(
    title="",
      x="Updated global ADM-0 mean TB counts [log scale]",
      y="Global ADM-0 mean TB counts [log scale]\nas reported in Table 1"
  ) +
  theme_bw()
#save updated Table 1
ggsave(file.path(path_output,"pdf","Updated_Table1.pdf"),p1, width=7, height=5.5)


# ADM1: comparison

adm1list <- lapply(countries,function(cc) {
  f <- file.path(path_input,paste0("gadm36_",cc,"_1_sp.rds"))
  if (file.exists(f)) readRDS(f) else NULL
})
adm1poly <- do.call(rbind,adm1list[!vapply(adm1list,is.null,logical(1))])

# Updated ADM1 estimates
updated_count <- raster::extract(
  TBcount_updated, adm1poly,
  fun=sum, na.rm=TRUE, sp=TRUE
)

adm1_compare <- data.frame(
  Admin_0=updated_count@data$NAME_0,
  GID_1=updated_count@data$GID_1,
  Admin_1=updated_count@data$NAME_1,
  updated_count_mean=updated_count@data$TB_count_mean,
  updated_count_lower=updated_count@data$TB_count_lower,
  updated_count_upper=updated_count@data$TB_count_upper
)

# Load reported ADM1 values
# Change filename if necessary
reported_adm1 <- read.csv(
  file.path(path_input, "supdata3.csv")
)

# Join reported and updated estimates
adm1_compare <- adm1_compare %>%
  dplyr::left_join(
    reported_adm1 %>%
      dplyr::select(
        Admin_0, Admin_1,
        reported_count=median.TB.cases,
        reported_count_lower=TB.cases.lower.estimate,
        reported_count_upper=TB.cases.upper.estimate
      ),
    by=c("Admin_0","Admin_1"),
  relationship="many-to-one"
  )

# Keep complete observations 
adm1_compare_plot <- adm1_compare %>%
  dplyr::filter(
    !is.na(reported_count),
    !is.na(updated_count_mean),
    !is.na(reported_count_lower),
    !is.na(updated_count_lower),
    !is.na(reported_count_upper),
    !is.na(updated_count_upper)
  )
# Pearson correlation f
cor_spearman1 <- cor(
  adm1_compare$reported_count,
  adm1_compare$updated_count_mean,
  use="complete.obs",
  method="spearman"
)

cor_pearson1 <- cor(
  adm1_compare$reported_count,
  adm1_compare$updated_count_mean,
  use="complete.obs",
  method="pearson"
)
# correlation by country
cor_by_country <- adm1_compare_plot %>%
  dplyr::group_by(Admin_0) %>%
  dplyr::summarise(
    pearson_r=cor(
      reported_count,
      updated_count_mean,
      use="complete.obs",
      method="pearson"
    ),
    .groups="drop"
  )
# Plot
p <- ggplot(adm1_compare_plot,
            aes(x=updated_count_mean, y=reported_count)) +
  geom_point(size=.8, alpha=.4) +
  geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=.4) +

  # Pearson r in top-left of each facet
  geom_label(
    data=cor_by_country,
    aes(x=-Inf, y=Inf, label=paste0("Pearson r = ",round(pearson_r,3))),
    inherit.aes=FALSE,
    hjust=-.05, vjust=1.15, size=3
  ) +

  facet_wrap(~Admin_0, scales="fixed",nrow=3)+
  coord_equal()+
  # Log axes, but display full values with commas
  scale_x_log10(
  labels=scales::label_comma(),
  expand=expansion(mult=c(.08,.12))
) +
scale_y_log10(
  labels=scales::label_comma(),
  expand=expansion(mult=c(.08,.12))
)+
  labs(
    x="Updated global ADM-1 mean TB counts [log scale]",
    y="Global ADM-1 TB counts [log scale]\nas reported in Sup. Data 3",
     caption=paste0(
      "Pearson's r = ", round(cor_pearson1,3),
      "\nSpearman's rho = ", round(cor_spearman1,3),
      "\n(based on raw counts)"
    )
  ) +
  theme_bw() + theme(
    strip.background=element_blank(),
    strip.text=element_text(size=10),
    axis.text.x=element_text(size=6),
   axis.text.y=element_text(size=6),
    plot.caption=element_text(size=7, hjust=1),
    panel.grid.minor=element_blank()
  )
#save ADM-1 plot
ggsave(
  file.path(path_output,"pdf","Updated_ADM1.pdf"),
  p, width=8.5, height=6.25
)

# Save correlations for ADM0 and ADM1
cordf <- data.frame(
  adm=c("ADM0","ADM1"),
  pearson=c(cor_pearson,cor_pearson1),
  spearman=c(cor_spearman,cor_spearman1),
  quantities="reported_vs_updated"
)

write.csv(
  cordf,
  file.path(path_output,"csv","correlation.csv"),
  row.names=FALSE
)
# DONE
###############################################################################
message( "Done. Outputs saved in: ", path_output )
