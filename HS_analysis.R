# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "MBA","reshape2","magrittr","RColorBrewer")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw(base_size=18))
library(gganimate)
library(RColorBrewer)
library(magrittr)
library(MBA)
library(reshape2)

# Colour palette for the temp-oxy index
idx_pal <- c("#252525",brewer.pal(n = 11, name = "RdYlGn")[-c(6,11)])

# Load and plot data ---------------------------------------------------

# Define site order
site.order <- c("5km Center", "Polly's Point 2", "Hohm Island 2", "River/Outfall") # Order of site names
site.labels <- c("5km Center", "Polly's Pt.", "Hohm Isle", "Estuary") # Shorter site names

#Load the raw data
hs <- read.csv("Catalyst_HarbourSurvey_2022_all.csv", na.strings = c("")) %>% 
  select(-pH) %>% 
  filter(!depth %in% c("River","River:")) %>% 
  mutate(date = date %>% paste("-2022",sep = "") %>% as.Date(format = "%d-%b-%Y"), #specify year = 2022
         julian = date %>% format("%j"),
         across(c(depth,julian), as.integer),
         # Add lat/long data for each site
         lat = case_when(
           site == "5km Center" ~ 49.201260,
           site == "Polly's Point 2" ~ 49.215729,
           site == "Hohm Island 2" ~ 49.227984,
           site == 'River/Outfall' ~ 49.247812
         ),
         lon = case_when(
           site == "5km Center" ~ -124.820737,
           site == "Polly's Point 2" ~ -124.821066,
           site == "Hohm Island 2" ~ -124.823543,
           site == 'River/Outfall' ~ -124.816844
         ),
         site = factor(site, levels = site.order,labels = site.labels)) %>% 
  filter(!between(julian, 91,121)) # Remove April surveys conducted using handheld probe on HS advice

# set date range for plots.
max_jul <- max(hs$julian) 
min_jul <- max_jul - 40



# Create function that interpolates values between surveys
intp_fn <- function(station, var) {
  df <- hs %>% 
    filter(site == station, !is.na({{var}})) %>% 
    select(julian,depth,{{var}}) %>% 
    mba.surf(no.X = max_jul, no.Y = max_jul, extend = T,
             # Parameters to control the stretch along the x axis
             n = 1, m = 2, # Higher n:m ratio means more x stretch, less y stretch
             # Set bounding box to ensure data are interpolated over the same date range for all sites
             b.box = c(min(hs$julian),max(hs$julian),min(.$depth),max(.$depth))) 
  
  dimnames(df$xyz.est$z) <- list(df$xyz.est$x, df$xyz.est$y)
  
  melt(df$xyz.est$z, varnames = c('julian', 'depth')) %>% 
    mutate(julian = round(julian, 0),
           date = as.Date(julian, origin = "2021-12-31"),
           depth = round(depth, 0)) %>% 
    group_by(date, julian, depth) %>% 
    summarize(value = mean(value)) %>% 
    ungroup()
}

# Apply function to all sites in the data and combine into single dfs
sites <- purrr::set_names(levels(as.factor(hs$site))) # List of site names


#Temperature
temps <- map_df(sites, ~intp_fn(.x, temp),.id = "site") %>% 
  mutate(site = factor(site, levels = site.labels))

#DO
do <- map_df(sites, ~intp_fn(.x, do_mgl),.id = "site") %>% 
  mutate(site = factor(site, levels = site.labels))

#Salinity
sal <- map_df(sites, ~intp_fn(.x, salinity),.id = "site") %>% 
  mutate(site = factor(site, levels = site.labels))


# Raster plot for time series data (based on temp data, tweaks needed for DO and salinity)
ts_p_fn <- function(df, hs_var, min_julian, max_julian) {
  df %>% 
    filter(between(julian, min_julian, max_julian)) %>% # Truncate interpolated data to desired time period
    ggplot(aes(x = date, y = depth)) +
    facet_grid(site~., scales = "free_y", space = "free_y") +
    geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
    geom_point(data = hs %>% 
                 filter(between(julian, min_julian, max_julian)) %>% # Truncate raw data
                 drop_na(contains(hs_var)), 
               aes(x = date, y = depth),
               size = 0.2, alpha = 0.15, shape = 8) +
    scale_fill_distiller(palette="RdYlGn", name = "Temperature (°C)",
                         labels = as.integer, breaks = c(5,10,15,20), 
                         limits = c(min(temps$value), max(temps$value)),
                         na.value = "transparent") +
    geom_contour(aes(z = value), binwidth = 2, colour = "black", alpha = 0.2) +
    #geom_contour(aes(z = value), breaks = 5, colour = "black") +
    labs(y = "Depth (m)", x = NULL) +
    scale_y_reverse(expand = c(0,0), labels = as.integer) +
    scale_x_date(expand = c(0,0), breaks = "2 weeks", date_labels = "%d %b") +
    theme(strip.background = element_rect(fill = "white"),
          legend.position = c(0.68,0.62),
          legend.direction = "horizontal",
          legend.justification = c("left", "bottom"),
          legend.background = element_rect(colour = "black",fill = alpha("white",0.75)),
          panel.spacing.y = unit(1, "lines"),
          axis.title.y = element_text(vjust = 2))
}


# Render plots for each variable
ts_p_fn(temps, "temp", min_jul, max_jul) # Temperature

ts_p_fn(do, "do_mgl", min_jul, max_jul) + #DO (plus appropriate legend tweaks)
  scale_fill_distiller(palette="RdYlGn", trans = "reverse", name = "DO (mg/L)",
                       labels = as.integer, breaks = c(3,6,9,12),
                       limits = c(max(do$value),min(do$value)))

sal_p <- ts_p_fn(sal, "salinity", min_jul, max_jul) + #Salinity (plus appropriate legend tweaks)
  scale_fill_distiller(direction = 1, name = "Salinity (ppt)",
                       limits = c(min(sal$value), max(sal$value))) 

sal_p[["layers"]][[3]][["stat_params"]][["binwidth"]] <- 5 #5ppm contours

sal_p

  

# Holding conditions index ------------------------------------------------

# What values correspond to the mean and +/- 1 SD of the two variables?
hs %>% 
  summarize(across(temp:do_mgl, 
                   .fns = list(mean=mean,sd=sd), 
                   .names = "{.fn}_{.col}", 
                   na.rm = T))

# We want the index to reflect moderate temps as 15C not 11.5C and 
# moderate do as 5ppm not 6.65 ppm, so we opt to manually enter the
# center and scale values in the index calculation:

# Create z-scored variables for temp and do and combine into holding index
hs %<>% 
  mutate(temp_cs = -scale(temp, center = 15, scale = 3),
         do_cs = scale(do_mgl, center = 5, scale = 2)) %>% 
  #mutate(temp_cs = -temp_cs) %>%   # Flip the sign on the temp index
  rowwise() %>% 
  mutate(hold_idx = case_when(  # calculate combined index
    do_mgl < 5 ~ mean(c(3*do_cs, temp_cs)),
    between(do_mgl,5, 6) ~ mean(c(2*do_cs, temp_cs)),
    do_mgl > 6 ~ mean(c(do_cs, temp_cs))
  )) %>%
  ungroup() %>% #Strip rowwise nature
  # Add categorical index based on Howard Stiff's recommendations
  mutate(idx_cat = case_when(
    temp < 12 & do_mgl > 4 ~ 5,
    between(temp, 12,16) & do_mgl > 4 ~ 4,
    between(temp, 16,18) & do_mgl > 4 ~ 3,
    temp < 18 & between(do_mgl,3,4) ~ 2,
    between(temp, 18,24) | between(do_mgl,2,3) ~ 1,
    temp > 24 | do_mgl < 2 ~ 0 # Catastrophic levels
    ))




# Interpolate time series of the continuous index
idx <- map_df(sites, ~intp_fn(.x, hold_idx), .id = "site") %>% 
  mutate(site = factor(site, levels = site.labels))


# Plot from the continuous index 
idx_p <- ts_p_fn(idx, "hold_idx", min_jul,max_jul) + 
  scale_fill_distiller(palette="RdYlGn", trans = "reverse", name = "Index",
                       labels = as.integer,
                       limits = c(max(idx$value),min(idx$value)))

idx_p[["layers"]][[3]][["aes_params"]][["binwidth"]] <- 1 # Change contour spacing

idx_p





# Plot from the categorical index
idx_c <- map_df(sites, ~intp_fn(.x, idx_cat), .id = "site") %>% 
  # Constrain values to [5,0]
  # Equation below from: https://stats.stackexchange.com/a/281164
  mutate(value = ((5-0)*((value - min(value))/(max(value)-min(value))))+0,
         site = factor(site, levels = site.labels)) 


idxcat_p <- ts_p_fn(idx_c, "idx_cat", min_jul, max_jul) + 
  # Custom colour scale based on the index and Howard Stiff's suggestions
  scale_fill_gradientn(colours = idx_pal,
                       name = "Temp-oxy\nindex",
                       limits = c(0,5))

idxcat_p[["layers"]][[3]][["stat_params"]][["binwidth"]] <- 1 # Change contour spacing

idxcat_p






# Add annotation for heat dome (from 2021, but left in case other annotations useful)
ann_text <- data.frame(date = as.Date("2021-07-02"), depth = 15,
                       site = factor("Polly's Pt.", levels = site.labels),
                       text = "Heat dome strikes\nsouthern BC")

idxcat_p +
  annotate("rect", fill = "white", alpha = 0.4,
           xmin = as.Date("2021-06-25"),xmax = as.Date("2021-07-01"),
           ymin = -Inf, ymax = Inf) +
  geom_label(data = ann_text, aes(label = text),
             alpha = 0.75,nudge_x = 10,label.r = unit(0, "lines"))


# Cross-section of the Inlet from most recent survey ----------------------

# Translate site names into their physical distances from the river mouth
hs %<>% 
  mutate(dist_km = case_when(
    site == "5km Center" ~ 5,
    site == "Hohm Isle" ~ 2.24,
    site == "Polly's Pt." ~ 3.45,
    site == "Estuary" ~ 0),
    # Shift the measurements at river and 5km mark for plotting
    plotting_shift = case_when(
      dist_km == 5 ~ 4.95, 
      dist_km == 0 ~ 0.05, 
      TRUE ~ dist_km
      ))

# Create function that interpolates values between sites
intp_dist_fn <- function(var, day_of_year) {
  df <- hs %>% 
    filter(julian == day_of_year) %>% 
    #Add some rows to extend interpolated values down to 60 m. 
    #Assumes linear water column below bottom deepest measurement at 5km site
    add_row(depth = c(55,60), site = rep("5km Center",2)) %>% 
    filter(!(site %in% site.order[-1] & is.na(.data[[var]]))) %>% #.data[[]] for when variable names are given in quotes
    arrange(site, depth) %>% 
    group_by(site) %>% 
    fill(c(salinity:dist_km), .direction = "down") %>% 
    ungroup() %>% 
    # Now keep only the columns needed for interpolation
    select(dist_km,depth,.data[[var]]) %>% 
    mba.surf(no.X = 100, no.Y = 100, extend = T,
             n = 6, m = 1) # Parameters to control the stretch along the x axis
  
  dimnames(df$xyz.est$z) <- list(df$xyz.est$x, df$xyz.est$y)
  
  melt(df$xyz.est$z, varnames = c('dist_km', 'depth')) %>% 
    mutate(depth = round(depth, 0))
}

# Apply function to all sites in the data and combine into single dfs
vars <- purrr::set_names(names(select(hs, c(salinity:do_mgl,idx_cat)))) # List of variable names

# All variables interpolated for every survey of 2021
#d1 <- crossing(vars,unique(hs$julian)) %>% rename(julian = 2) #df that aligns vars and days

# pass df1 to map2 to get all combos of variable and day.
#xs_vars <- map2_df(d1$vars,d1$julian,~intp_dist_fn(.x,.y)) # Causes R studio to crash
#xs_vars <- map_df(vars, ~intp_dist_fn(.x, max(hs$julian)), .id = "var")

#d2 <- filter(d1, vars == "idx_cat") 
#idx_ts <- map2(d2$vars, d2$julian, ~intp_dist_fn(.x,.y)) 

# Manually input bathymetric data from web application: https://data.chs-shc.ca/map
bathy <- data.frame(dist_km = c(-0.1,-0.1,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,4.75,4.8,4.85, 5.01),
                    depth = c(61,10,12, 15, 18, 20, 24, 30, 33, 33, 35, 35,46,55,61))

# Plotting function
xs_p_fn <- function(df) {
  df %>% 
    ggplot(aes(x = dist_km, y = depth)) +
    #facet_wrap(~site, scales = "free_y", strip.position = "right", ncol = 1) +
    geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
    geom_point(data = hs %>%
                 ## Update below when complete surveys become available
                 filter(julian == max(julian)) %>%
                 drop_na(salinity:do_sat), 
               aes(x = plotting_shift, y = depth),
               size = 0.2, alpha = 0.15, shape = 8) +
    scale_fill_distiller(palette="RdYlGn", name = "Temperature (°C)",
                         labels = as.integer, breaks = c(5,10,15,20), 
                         limits = c(min(temps$value), max(temps$value))) +
    geom_contour(aes(z = value), binwidth = 2, colour = "black", alpha = 0.2) +
    #geom_contour(aes(z = value), breaks = 5, colour = "black") +
    geom_polygon(data = bathy, fill = "grey80", colour = "black") +
    labs(y = "Depth (m)", x = "Distance from river mouth (km)") +
    coord_cartesian(xlim = c(0,5), expand = F) +
    scale_y_reverse(labels = as.integer) +
    theme(strip.background = element_rect(fill = "white"),
          legend.position = c(0.02,0.2),
          legend.direction = "horizontal",
          legend.justification = c("left", "bottom"),
          legend.background = element_rect(colour = "black"),
          panel.spacing.y = unit(1, "lines"),
          axis.title.y = element_text(vjust = 2))
}

# Plot for temperature
#xs_p_fn(xs_vars, "temp")

# Plot for temp-oxy index
idx_xs_p <- xs_p_fn(intp_dist_fn("idx_cat",max(hs$julian)) %>% 
          # Squish the interpolated index values into [0,5]
          #filter(var == "idx_cat") %>% 
          mutate(value = (5-0)*((value - min(value))/(max(value)-min(value)))+0)) +
  coord_cartesian(ylim = c(60, 0), xlim = c(0,5), expand = F) + # Chop it off at 25m because 5km site not surveyed any deeper.
  # Custom colour scale based on the index and Howard Stiff's suggestions
  scale_fill_gradientn(colours = idx_pal,
                       name = "Temp-oxy\nindex", limits = c(0,5))

idx_xs_p[["layers"]][[3]][["stat_params"]][["binwidth"]] <- 1 # Change contour spacing

idx_xs_p



## Need to figure out similar fix to the scale limits as I applied above for the index. 
## Maybe split xs_vars into separate df for each measure

# DO plot
ts_p_fn(xs_vars, "do_mgl") + 
  scale_fill_distiller(palette="RdYlGn", trans = "reverse", name = "DO (mg/L)",
                       labels = as.integer, breaks = c(3,6,9,12),
                       #limits = c(max(do$value),min(do$value)))

#Salinity plot
sal_p <- ts_p_fn(xs_vars, "salinity") + 
  scale_fill_distiller(direction = 1, name = "Salinity (ppt)",
                       #limits = c(min(sal$value), max(sal$value))) 

sal_p[["layers"]][[3]][["stat_params"]][["binwidth"]] <- 5 #5ppm contours

sal_p

