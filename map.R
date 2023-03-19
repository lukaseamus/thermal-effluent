##########################################################################################
## Project: Impacts of thermal effluent on Posidonia oceanica and associated macrofauna ##
## Script purpose: Visualisation of site map                                            ##
## Author: Luka Seamus Wright                                                           ##
##########################################################################################

#### 1.  Spatial data ####
#### 1.1 Sites ####
sites <- data.frame(site = c("KS", "KA", "KB", "KC", "KD",
                             "ZS", "ZA", "ZB", "ZC", "ZD"),
                    lat = c(35.84058, 35.84042, 35.84063, 35.84058,
                            35.84072, 35.83592, 35.83648, 35.83623,
                            35.83603, 35.83605),
                    lon = c(14.56192, 14.56268, 14.5627, 14.56223,
                            14.56208, 14.55981, 14.56013, 14.56045,
                            14.5603, 14.56018))

sites$site <- factor(sites$site, levels = c("ZS", "ZA", "ZB", "ZC", "ZD",
                                            "KS", "KA", "KB", "KC", "KD"))

# KS = Kbira source
# KA = Kbira site A
# KB = Kbira site B
# KC = Kbira site C
# KD = Kbira site D
# ZS = Zghira source
# ZA = Zghira site A
# ZB = Zghira site B
# ZC = Zghira site C
# ZD = Zghira site D

transects <- data.frame(line = c(rep("K1", 2), rep("K2", 2), rep("K3", 2), 
                                 rep("K4", 2), rep("K5", 2), rep("Z1", 2),
                                 rep("Z2", 2), rep("Z3", 2), rep("Z4", 2),
                                 rep("Z5", 2)),
                        lat = c(35.84058, 35.8414269, 35.84058, 35.8412704,
                                35.84058, 35.84103063, 35.84058, 35.8407365,
                                35.84058, 35.84042349, 35.83592, 35.83681783,
                                35.83592, 35.83680757, 35.83592, 35.83665827,
                                35.83592, 35.83630088, 35.83592, 35.83591999),
                        lon = c(14.56192, 14.56229858, 14.56192, 14.56263149,
                                14.56192, 14.56287858, 14.56192, 14.56301006,
                                14.56192, 14.56301005, 14.55981, 14.55971353,
                                14.55981, 14.5600022, 14.55981, 14.56044484,
                                14.55981, 14.56081311, 14.55981, 14.56091681))

# K1 = Kbira transect 1
# K2 = Kbira transect 2
# K3 = Kbira transect 3
# K4 = Kbira transect 4
# K5 = Kbira transect 5
# Z1 = Zghira transect 1
# Z2 = Zghira transect 2
# Z3 = Zghira transect 3
# Z4 = Zghira transect 4
# Z5 = Zghira transect 5

  
#### 1.2 Base map ####
require(sf)
crop <- c(xmin = 14.1766, xmax = 14.575, ymin = 35.786, ymax = 36.09)
base <- read_sf("~/PATH/land-polygons-complete-4326/land_polygons.shp") %>% st_crop(crop)
seagrass <- read_sf("~/PATH/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC013014-Seagrasses-Py-v7_1.shp")
seagrass <- st_transform(seagrass, crs = 4326)
seagrass <- st_make_valid(seagrass)
seagrass <- seagrass %>% st_crop(crop)
sf_use_s2(FALSE)

#### 2.  Visualisation ####
#### 2.1 Set theme ####
require(ggplot2)
mytheme <- theme(plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill = NA, size = 1),
                 axis.line = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_text(size = 12, face = "bold"),
                 legend.background = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 2.2 Plot country map ####
malta <- ggplot() +
  geom_sf(seagrass, mapping = aes(), colour = NA, fill = "#81a512") +
  geom_sf(base, mapping = aes(), colour = NA, fill = "#b5b8ba") +
  geom_rect(aes(xmin = 14.5588, xmax = 14.5725, ymin = 35.834, ymax = 35.844), 
            fill = NA, colour = "#000000", size = 0.5, linejoin = "mitre") +
  mytheme +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank())
malta


#### 2.3 Plot site map ####
site <- ggplot() +
  geom_sf(base, mapping = aes(), colour = NA, fill = "#b5b8ba") +
  geom_line(data = transects, mapping = aes(x = lon, y = lat, group = line),
            size = 0.5, colour = "#cdd0d1") +
  geom_point(data = sites, mapping = aes(x = lon, y = lat, fill = site, shape = site),
             size = 1.7, colour = "#b5b8ba", stroke = 0.5) +
  geom_rect(mapping = aes(xmin = 14, xmax = 14.5588, ymin = 35, ymax = 36),
            fill = "#ffffff", colour = "#000000") +
  scale_fill_manual(values = c(rep("#f5a54a", 5), rep("#6ea4be", 5)),
                    labels = c("Source", "Plot A", "Plot B", "Plot C", "Plot D",
                               "Source", "Plot A", "Plot B", "Plot C", "Plot D"),
                    guide = guide_legend(ncol = 2, title = "  Il-Ħofra       Il-Ħofra \n  ż-Żgħira      l-Kbira")) +
  scale_shape_manual(values = rep(c(21, 22, 23, 24, 25), 2),
                     labels = c("Source", "Plot A", "Plot B", "Plot C", "Plot D",
                                "Source", "Plot A", "Plot B", "Plot C", "Plot D"),
                     guide = guide_legend(ncol = 2, title = "  Il-Ħofra       Il-Ħofra \n  ż-Żgħira      l-Kbira")) +
  theme(legend.position = c(0.18, 0.25)) +
  scale_y_continuous(breaks = seq(35.834, 35.844, by = 0.0025), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(14.5525, 14.5725, by = 0.005), expand = c(0, 0)) +
  coord_sf(xlim = c(14.5525, 14.5725), ylim = c(35.834, 35.844)) +
  mytheme
site

#### 2.4 Plot combined map ####
require(cowplot)
combined <- ggdraw(site) +
            draw_plot(malta, width = 0.48, height = 0.48,
                      x = 0.035, y = 0.49)
combined
# dimensions: 4.5 x 8.3 in

#### 3.  Clean up ####
detach(package:sf)
detach(package:ggplot2)
detach(package:cowplot)
rm(list = ls())
graphics.off()
cat("\014")
