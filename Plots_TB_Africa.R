# Replication code for the paper "Mapping TB in Africa" by T. K. et al. (2024)
if (!exists("allrun"))        allrun <- FALSE
# To make all plots (main and SI): all four model configurations =============
# run all model specs or only the main spec.
if (allrun) {
  combinations <- expand.grid(popfilter = c(FALSE, TRUE),mozout = c(FALSE, TRUE))
} else {
  # Main specification only
  combinations <- data.frame(popfilter = FALSE,mozout = FALSE)
}

# Packages ===============================================================

list.of.packages <- c(
  "raster", "rnaturalearth", "sf", "rnaturalearthdata",
  "dplyr", "scales", "ggplot2", "ggspatial", "grid"
)

new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))


# Paths ==================================================================

path_input <- file.path(getwd(), "INPUT")
output_path <- file.path(getwd(), "OUTPUT")

dir.create(path_input, recursive = TRUE, showWarnings = FALSE)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Africa =================================================================

Africa_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
Africa_sf <- Africa_sf[Africa_sf$continent == "Africa", ]

# Common scale bar =======================================================

sy <- -37
sx <- -15.5
kmdeg <- 111.32 * cos(sy * pi / 180)

major_km <- c(0, 1000, 2000, 4000)
minor_km <- c(500, 1500, 2500, 3000, 3500)

major_x <- sx + major_km / kmdeg
minor_x <- sx + minor_km / kmdeg

major_ticks <- data.frame(
  x = major_x, xend = major_x,
  y = sy, yend = sy + 1.5
)

minor_ticks <- data.frame(
  x = minor_x, xend = minor_x,
  y = sy, yend = sy + .9
)


# Continuous colour palette =============================================

tb_cols <- c(
  "#59A0CB",  # blue
  "#87B2BE",
  "#A4C4B3",
  "#C1D6A7",
  "#DBE99A",
  "#F3FA8D",  # yellow
  "#F2D777",
  "#EEB663",
  "#EA934F",
  "#E46D3E",
  "#D84B40"   # red
)

# Discrete colour palette ================================================

agg_cols <- c(
  "#2166B1", "#D2DCB0", "#F2F86D",
  "#EFB951", "#E5743B", "#A32F35"
)


# =======================================================================
# COMMON MAP ELEMENTS
# =======================================================================

add_map_elements <- function(p) {

  p +
    annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(.30, "in"), pad_y = unit(.18, "in"),
      height = unit(.3, "in"), width = unit(.3, "in"),
      style = north_arrow_fancy_orienteering(
      line_col = "black",
     fill = c("black", "white")
)
    ) +
    annotate(
      "segment",
      x = major_x[1], xend = major_x[4],
      y = sy, yend = sy,
      linewidth = .55
    ) +
    geom_segment(
      data = major_ticks,
      aes(x, y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = .55
    ) +
    geom_segment(
      data = minor_ticks,
      aes(x, y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = .45
    ) +
    annotate(
      "text",
      x = major_x,
      y = sy + 2.2,
      label = c("0", "1,000", "2,000", "4,000"),
      size = 3.1
    ) +
    annotate(
      "text",
      x = major_x[4] + 2,
      y = sy + 2.2,
      label = "kilometers",
      hjust = 0,
      size = 3
    )
}


# =======================================================================
# CONTINUOUS RASTER MAP FUNCTION
# =======================================================================

plot_TB_map <- function(r, legend_title, filename, path_output) {

  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("longitude", "latitude", "value")

zmin  <- min(df$value, na.rm = TRUE)
zmax  <- max(df$value, na.rm = TRUE)
zmean <- mean(df$value, na.rm = TRUE)
zsd   <- sd(df$value, na.rm = TRUE)

# ArcGIS Standard Deviation stretch (2 SD)
stretch_min <- max(zmin, zmean - 2 * zsd)
stretch_max <- min(zmax, zmean + 2 * zsd)

# Transform to 0-1 display scale and clamp tails
df$stretch <- scales::rescale(
  df$value,
  from = c(stretch_min, stretch_max),
  to = c(0, 1)
)

df$stretch <- pmax(0, pmin(1, df$stretch))

  p <- ggplot() +
    geom_sf(
      data = Africa_sf,
      fill = "#E6E6E6",
      colour = "#777777",
      linewidth = .3
    ) +
    geom_raster(
      data = df,
      aes(longitude, latitude, fill = stretch)
    ) +
    geom_sf(
      data = Africa_sf,
      fill = NA,
      colour = "#777777",
      linewidth = .3
    ) +
    scale_fill_gradientn(
  colours = tb_cols,
  values = seq(0, 1, length.out = length(tb_cols)),
  limits = c(0, 1),
  guide = "none"
)+
    coord_sf(
      xlim = c(-20, 55),
      ylim = c(-38.5, 39),
      expand = FALSE,
      datum = NA
    ) +
    theme_void() +
    theme(
      panel.background = element_rect(
        fill = "#C6E7FC",
        colour = "black",
        linewidth = .65
      ),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5)
    )

  p <- add_map_elements(p)

  # Custom continuous legend
  lx1 <- -17
  lx2 <- -13.2
  ly1 <- -22
  ly2 <- -8
  nleg <- 200

  leg <- data.frame(
  x = (lx1 + lx2) / 2,
  y = seq(ly1, ly2, length.out = nleg),
  stretch = seq(0, 1, length.out = nleg)
)

  p <- p +
    geom_tile(
  data = leg,
  aes(x, y, fill = stretch),
  width = lx2 - lx1,
  height = (ly2 - ly1) / (nleg - 1),
  inherit.aes = FALSE
) +
    annotate(
      "rect",
      xmin = lx1, xmax = lx2,
      ymin = ly1, ymax = ly2,
      fill = NA,
      colour = "grey65",
      linewidth = .25
    ) +
    annotate(
      "text",
      x = lx1, y = ly2 + 4.5,
      label = legend_title,
      hjust = 0,
      size = 4.5
    ) +
    annotate(
      "text",
      x = lx2 + .6, y = ly2,
      label = sprintf("%.2f", zmax),
      hjust = 0, vjust = 1,
      size = 3.2
    ) +
    annotate(
      "text",
      x = lx2 + .6, y = ly1,
      label = sprintf("%.2f", zmin),
      hjust = 0, vjust = 0,
      size = 3.2
    ) +
    annotate(
  "rect",
  xmin = lx1, xmax = lx2,
  ymin = ly1 - 5, ymax = ly1 - 1.8,
  fill = "#E6E6E6",
  colour = "grey65",
  linewidth = .25
) +
annotate(
  "text",
  x = lx2 + .6, y = ly1 - 3.4,
  label = "no data",
  hjust = 0,
  size = 3.2
)

  ggsave(
    file.path(path_output, "pdf", paste0(filename, ".pdf")),
    p, width = 7, height = 7.2, units = "in"
  )

  ggsave(
    file.path(path_output, "pdf", paste0(filename, ".png")),
    p, width = 7, height = 7.2, units = "in", dpi = 400
  )

  p
}


# =======================================================================
# AGGREGATED ADM0 / ADM1 MAP FUNCTION
# =======================================================================

plot_agg_map <- function(
  x, var, breaks, legend_title, filename, path_output,
  border_width = .3, digits = 2) {

  zmin <- min(x[[var]], na.rm = TRUE)
  zmax <- max(x[[var]], na.rm = TRUE)

  breaks <- sort(unique(as.numeric(breaks)))
  breaks <- breaks[breaks > zmin & breaks < zmax]

  if (!length(breaks))
    stop("No valid internal class breaks for ", var)

  brks <- c(zmin, breaks, zmax)

  fmt <- function(z)
    formatC(z, format = "f", digits = digits, big.mark = ",")

  labels <- paste0(
    fmt(head(brks, -1)),
    " - ",
    fmt(tail(brks, -1))
  )

  x$map_class <- cut(
    x[[var]],
    breaks = brks,
    labels = labels,
    include.lowest = TRUE,
    right = TRUE
  )

  class_cols <- colorRampPalette(agg_cols)(length(labels))
  cols <- setNames(class_cols, labels)

  p <- ggplot() +

    # Africa background = no data
    geom_sf(
      data = Africa_sf,
      aes(colour = "no data"),
      fill = "#E6E6E6",
      linewidth = .3
    ) +

    # Estimated ADM0/ADM1 values
    geom_sf(
      data = x,
      aes(fill = map_class),
      colour = "#777777",
      linewidth = border_width
    ) +

    # Country borders
    geom_sf(
      data = Africa_sf,
      fill = NA,
      colour = "#777777",
      linewidth = .35,
      show.legend = FALSE
    ) +

    # Statistical classes
    scale_fill_manual(
  values = cols,
  drop = FALSE,
  na.translate = FALSE,
  name = legend_title
)+

    # Separate no-data legend
    scale_colour_manual(
      values = c("no data" = "#777777"),
      name = NULL,
      guide = guide_legend(
        order = 2,
        override.aes = list(
          fill = "#E6E6E6",
          colour = "#777777",
          linewidth = .3
        )
      )
    ) +

    guides(
      fill = guide_legend(order = 1)
    ) +

    coord_sf(
      xlim = c(-20, 55),
      ylim = c(-38.5, 39),
      expand = FALSE,
      datum = NA
    ) +

    theme_void() +

    theme(
      panel.background = element_rect(
        fill = "#C6E7FC",
        colour = "black",
        linewidth = .65
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.position = c(.05, .1),
      legend.justification = c(0, 0),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(.38, "cm"),
      legend.key.spacing.y = unit(.15, "cm"),
      legend.background = element_blank(),
      plot.margin = margin(5)
    )

  p <- add_map_elements(p)

  ggsave(
    file.path(path_output, "pdf", paste0(filename, ".pdf")),
    p, width = 7, height = 7.2, units = "in"
  )

  ggsave(
    file.path(path_output, "pdf", paste0(filename, ".png")),
    p, width = 7, height = 7.2, units = "in", dpi = 400
  )

  p
}

# =======================================================================
# RUN ALL FOUR CONFIGURATIONS
# =======================================================================

for (k in seq_len(nrow(combinations))) {

  popfilter <- combinations$popfilter[k]
  mozout <- combinations$mozout[k]

  subfiles <- ifelse(
    popfilter & mozout, "NOMOZ/FILTER",
    ifelse(
      popfilter & !mozout, "ALL/FILTER",
      ifelse(
        !popfilter & mozout, "NOMOZ/NOFILTER",
        "ALL/NOFILTER"
      )
    )
  )

  path_output <- file.path(output_path, subfiles)
  dir.create(file.path(path_output, "pdf"),
             recursive = TRUE, showWarnings = FALSE)

  message("\nProcessing: ", subfiles)


  # =====================================================================
  # RASTER PREVALENCE MAPS
  # =====================================================================

  raster.list <- list.files(
    file.path(path_output, "prevalence"),
    pattern = "\\.tif$",
    full.names = TRUE
  )

  b <- raster::brick(lapply(raster.list, raster::raster))

  p_mean <- plot_TB_map(
    b[["Prob_mean"]],
    "TB prevalence mean\n(per 1,000)",
    "TB_prevalence_mean",
    path_output
  )

  p_IQR <- plot_TB_map(
    b[["Prob_IQR"]],
    "Predicted TB prevalence\n(IQR per 1,000)",
    "TB_prevalence_IQR",
    path_output
  )
  
  p_sd <- plot_TB_map(
    b[["Prob_sd"]],
    "SD of predicted TB\nprevalance (per 1,000)",
    "TB_prevalence_sd",
    path_output
  )

  p_gmrf <- plot_TB_map(
    b[["spatial_field_mean"]],
    "mean spatial random field\n(GMRF)",
    "TB_prevalence_GMRF",
    path_output
  )

  # =====================================================================
  # LOAD ADM0 / ADM1 PREVALENCE
  # =====================================================================

  Af0TB_sf <- sf::st_read(
    file.path(path_output, "shp", "TBprevalence_adm0.shp"),
    quiet = TRUE
  )

  Af1TBprev_sf <- sf::st_read(
    file.path(path_output, "shp", "TBprevalence_adm1.shp"),
    quiet = TRUE
  )

  Af0TB_sf <- sf::st_transform(Af0TB_sf, 4326)
  Af1TBprev_sf <- sf::st_transform(Af1TBprev_sf, 4326)


  # ADM0 mean prevalence -------------------------------------------------

  p_adm0 <- plot_agg_map(
    Af0TB_sf,
    "meanTBprev",
    breaks = c(1.95, 2.68, 3.68, 4.01, 4.70),
    legend_title = "TB mean prevalence (per 1,000)",
    filename = "TB_mean_prevalence_ADM0",
    path_output = path_output
  )
 
  # ADM0 IQR prevalence -------------------------------------------------
 
    p_adm0 <- plot_agg_map(
    Af0TB_sf,
    "IQRTBprev",
   IQR_breaks <- c( 2.00, 3.20, 4.16, 4.50,  6.50),
    legend_title = "TB IQR prevalence\n(per 1,000)",
    filename = "TB_IQR_prevalence_ADM0",
    path_output = path_output
  )
  # ADM1 mean prevalence -------------------------------------------------

  p_adm1 <- plot_agg_map(
    Af1TBprev_sf,
    "meanTBprev",
    breaks = c(1.85, 2.68, 3.68, 4.77, 5.84),
    legend_title = "TB mean prevalence (per 1,000)",
    filename = "TB_mean_prevalence_ADM1",
    path_output = path_output,
    border_width = .15
  )


  # ADM1 prevalence IQR --------------------------------------------------

  IQR_breaks <- stats::quantile(
    Af1TBprev_sf$IQRTBprev,
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_adm1_IQR <- plot_agg_map(
    Af1TBprev_sf,
    "IQRTBprev",
    breaks =  c(1.6, 2.99, 4.19, 5.46, 7.74),
    legend_title = "TB prevalence IQR",
    filename = "TB_prevalence_IQR_ADM1",
    path_output = path_output,
    border_width = .15
  )


  # =====================================================================
  # LOAD ADM0 / ADM1 TB CASES
  # =====================================================================

  Af0TB <- sf::st_read(
    file.path(path_output, "shp", "TBcases_adm0.shp"),
    quiet = TRUE
  )

  Af1TB <- sf::st_read(
    file.path(path_output, "shp", "TBcases_adm1.shp"),
    quiet = TRUE
  )

  Af0TB <- sf::st_transform(Af0TB, 4326)
  Af1TB <- sf::st_transform(Af1TB, 4326)


  # Shapefile field names are truncated to 10 characters:
  # medianTBcount -> medianTBco
  # meanTBcount   -> meanTBcoun


  # ADM0 median TB cases -------------------------------------------------

  adm0_median_breaks <- stats::quantile(
    Af0TB[["medianTBco"]],
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_count_adm0_median <- plot_agg_map(
    Af0TB,
    "medianTBco",
    breaks = adm0_median_breaks[2:6],
    legend_title = "Median TB cases",
    filename = "TB_cases_median_ADM0",
    path_output = path_output,
    digits = 0
  )


  # ADM0 mean TB cases ---------------------------------------------------

  adm0_mean_breaks <- stats::quantile(
    Af0TB[["meanTBcoun"]],
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_count_adm0_mean <- plot_agg_map(
    Af0TB,
    "meanTBcoun",
    breaks = adm0_mean_breaks[2:6],
    legend_title = "Mean TB cases",
    filename = "TB_cases_mean_ADM0",
    path_output = path_output,
    digits = 0
  )


  # ADM0 IQR TB cases ----------------------------------------------------

  adm0_IQR_breaks <- stats::quantile(
    Af0TB[["IQRTBcount"]],
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_count_adm0_IQR <- plot_agg_map(
    Af0TB,
    "IQRTBcount",
    breaks = adm0_IQR_breaks[2:6],
    legend_title = "TB cases IQR",
    filename = "TB_cases_IQR_ADM0",
    path_output = path_output,
    digits = 0
  )


  # ADM1 mean TB cases ---------------------------------------------------

  adm1_mean_breaks <- stats::quantile(
    Af1TB[["meanTBcoun"]],
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_count_adm1_mean <- plot_agg_map(
    Af1TB,
    "meanTBcoun",
    breaks = adm1_mean_breaks[2:6],
    legend_title = "Mean TB cases",
    filename = "TB_cases_mean_ADM1",
    path_output = path_output,
    border_width = .15,
    digits = 0
  )


  # ADM1 IQR TB cases ----------------------------------------------------

  adm1_IQR_breaks <- stats::quantile(
    Af1TB[["IQRTBcount"]],
    probs = seq(0, 1, length.out = 7),
    na.rm = TRUE,
    names = FALSE
  )

  p_count_adm1_IQR <- plot_agg_map(
    Af1TB,
    "IQRTBcount",
    breaks = adm1_IQR_breaks[2:6],
    legend_title = "TB cases IQR",
    filename = "TB_cases_IQR_ADM1",
    path_output = path_output,
    border_width = .15,
    digits = 0
  )
}
