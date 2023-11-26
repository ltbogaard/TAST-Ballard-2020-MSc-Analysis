## TAST MSC CLEANING AND EDA SCRIPT
## Laura Bogaard

# Load data and libraries
library(tidyverse)
library(reshape2)
library(readr)
library(useful)
library(nlme)
library(mgcv)
library(lubridate)

# set plotting themes
library(extrafont)
loadfonts(device = "win")
windowsFonts(Times = windowsFont("TT Times New Roman"))
library(ggplot2); theme_set(theme_bw(base_size=12, base_family='Times New Roman')+
                              theme(
                                legend.title = element_text(size = 10),
                                legend.text = element_text(size = 10), 
                                plot.title = element_text(hjust = 0.5, size = 14),
                                legend.key.size = unit(5, "mm"), 
                                axis.text.x = element_text(size = 12),
                                axis.text.y = element_text(size = 12),
                                axis.title.x = element_text(size = 12),
                                axis.title.y = element_text(size = 12)))


#set wd 
setwd("~/Desktop/TASTBALLARD/TAST WD/Dissertation")
#######################################################
# load data and adjust data types
#### distance refers to the dataset where each row is a surfacing
distance <- read.csv("Ballard_distance+fishcts.csv", header = T)
distance$time <- parse_time(distance$time, format = "%T")
distance$treatment <- factor(distance$treatment)
distance <- unite(distance, col = "blockid", block:block_ns, sep = "", remove = FALSE, na.rm = FALSE)
head(distance)
# Seperate time into hms so we can use hour of day as a factor variable 
distance <- distance %>%
  separate(time, c("hour", "minute","second"),
         sep = ":", remove = TRUE)
distance$foraging <- ifelse(distance$foraging == "Yes", 1, 0) #make foraging binary
distance <- distance%>% filter(platform_distance <= 250)

#### pred refers to the dataset where each row is a surveying session
pred <- read.csv("Ballard_Predation.csv")
pred$session_id<- factor(pred$Session_ID)


pred$Date <- as.Date(pred$Date, format = "%Y-%m-%d")
pred$Start_Time <- parse_time(pred$Start_Time, format = "%T")
pred$End_Time <- parse_time(pred$End_Time, format = "%T")

# Seperate time into hms 
pred <- pred %>%
    separate(Start_Time, c("hour", "minute","second"), 
             sep = ":", remove = F)

  
#### fish refers to the dataset where each row is a date which corresponds # of fish
fish <- read_csv("fishpassage_surveyonly.csv", 
                 col_types = cols(date = col_date(format = "%m/%d/%y"),
                             treatment = col_factor(levels = c("on", "off"))))


################################################################################
#### RL data EDA + wrangling 
# RL refers to the dataset with distances and calculated RLs
RL <- data.frame(read_csv("ballard_rl.csv"))
 
#check                
head(RL)
tail(RL)
RL[19,]

#weird, remove last 6 rows of NAs
RL <- slice(RL, 1:(n()-6))

#rename distance columns
rldist <- as.numeric(RL$platform_distance) #distance from recording location to platform
location <- factor(RL$location) #whether recording was taken at edge of canal or centre
mean_rl <- as.numeric(RL$mean_rl) #mean received level in dB re 1uPa, assuming spherical spreading
plot(rldist, mean_rl)

# build simple linear model
rlmodemp <- lm(mean_rl ~ I(20*log10(rldist)), data = RL) 
rlmod_sp <- function(distance){182 - 20*log10(distance)}
rlmod_cy <- function(distance){182 - 10*log10(distance)}


summary(rlmod)
rlmod$coefficients
summary(rlmod)$adj.r.squared

mutate(RL_theoretical = as.numeric(182 - (20*log10(distance)))) %>% 
  mutate(RL_emperical = c(predict(rlmod, newdata = data.frame(rldist = distance))))

# For each value of distance, I can get the value of RL estimated by the model, 
library(ggplot2)
install.packages("ggnewscale")
library(ggnewscale)
# now do this in gg plot
RL_plot <- ggplot(data = RL) +
  geom_point(mapping = aes(x = rldist, y = mean_rl)) +
  ggtitle("Transmission Loss of TAST signal") + 
  labs(x = "Distance (m)", y = "Received level dB re 1µPa", linetype = "Model") +
  stat_smooth(mapping = aes(x = rldist, y = mean_rl, linetype = "Emperical"),
              colour = "black", method = "lm", formula = y ~ I(20*log(x)),  se=T, level=0.95) +
  stat_function(fun = function(rldist){182 - 20*log10(rldist)}, aes(linetype = "Spherical")) +
  stat_function(fun = function(rldist){182 - 10*log10(rldist)}, aes(linetype = "Cylindrical")) +
  theme_light() +
  theme(legend.position = c(0.85, 0.8), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.key.size = unit(12, "mm"), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(0, 5000, by = 50)) +
  scale_y_continuous(breaks = seq(100, 180, by = 10)) +
  scale_linetype_manual(values=c( "Cylindrical" ='dashed', "Spherical" = 'dotted', "Emperical" = 'solid')) 



RL_plot

# Boxplot for all distances
distance %>% 
  ggplot(mapping = aes(y = platform_distance, x = treatment, fill = treatment, alpha = 0.5)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("blue", "darkgreen"), name = "Treatment", labels = c("On", "Off")) +
  ylab("Distance") +
  xlab("Treatment") +
  theme_classic()+
  theme(legend.position="none")

# # plot using polar coords, interesting but not informative
# library(plotrix)
# plt.lns <- filter(valid_bearings, platform_distance <= 250)
# angles <- valid_bearings$platform_bearing
# polar.plot(plt.lns$platform_distance, polar.pos=angles, labels="", rp.type = "s")

# Overlapped Hist for Distance by Treatment (Final Figure for Report)
valid_bearings %>% 
  ggplot(aes(x = tide_dist, fill = treatment)) + 
  geom_histogram(binwidth = 10, color = "white", position = "dodge", alpha = 0.5) + 
  #scale_alpha_discrete(name = "Fish Species", range = c(0.4, 1), labels = c("Coho", "Chinook")) + 
  scale_fill_manual(name = "Treatment", values = c("blue", "red"), labels = c("Off", "On")) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  xlim(0, 250) +
  #scale_x_continuous(breaks = seq(0, 250, 50)) +
  theme_bw() +
  ylab("Frequency") + 
  xlab("Distance from TAST (m)") 

nrow(clean_dist2[treatment == "ON",]) + nrow(clean_dist2[treatment == "OFF",])

valid_bearings[treatment == "ON",]

########### ########### ########### ########### ########### ########### ########### 

################################################################################

### PREDATION DATA EDA ###

species <- pred$Phase != "Baseline"

#Boxplot comparing pred rate and salmon species/phase with baseline
# TODO: get rid of baseline somehow 
pred %>%
  ggplot(aes (y = pred_rate_per_hour, x = Phase !="Baseline", fill = Phase)) +
  geom_boxplot() +
  scale_fill_manual(name = "Salmon Species", values = c("white", "lightgrey", "darkgrey"), 
                    labels = c("Baseline", "Chinook", "Coho")) +
  theme_classic() +
  ylab("Predation Rate") +
  xlab("Salmon Species")

################## 
# FISH DATA EDA
#chinook by treatment
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = chinook_day, fill = treatment)) + geom_bar(stat = "identity")
#coho by treatment-
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = coho_day, fill = treatment)) + geom_bar(stat = "identity")

# Chinook and Coho by Treatment (FINAL for Report)
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  melt(id.vars = c("time", "treatment"), measure.vars = c("coho_day", "chinook_day")) %>% 
  ggplot(aes(x = time, y = value, fill = variable, alpha = treatment)) + 
  geom_col() + 
  scale_alpha_discrete(name = "Treatment", range = c(1, 0.4), labels = c("On", "Off")) + 
  scale_fill_discrete(name = "Fish Species", labels = c("Coho", "Chinook")) +
  theme_classic() +
  ylab("Fish count") + 
  xlab("Date") 

###################################### 
# WRANGLE DATA INTO PRESENCE ABSENCE
###################################### 
# get counts of surfacings within 10m by 10m bins for each session
library(dplyr)
grid_counts <- valid_bearings %>%  
  mutate(
    cut_x = cut(x, breaks = seq(0, 250, by = 10), include.lowest = T),
    cut_y = cut(y, breaks = seq(0, 250, by = 10), include.lowest = T),
  ) %>%
  count(cut_x, cut_y) %>% mutate(n_bin = n())

# make this a function 
gridcount <- function(data){
  data %>%  
  mutate(
    cut_x = cut(x, breaks = seq(0, 250, by = 10), include.lowest = T),
    cut_y = cut(y, breaks = seq(0, 250, by = 10), include.lowest = T),
    distance = 
  ) %>%
    count(cut_x, cut_y) %>% 
    mutate(n_bin = n())
}

# plot all data
grid_count_plot <- ggplot(valid_bearings, aes(x, y)) + 
  stat_bin_2d(binwidth = 10) + 
  geom_point(color = "white", cex = 0.25)
grid_count_plot


# plot on 
gridcount_ON <- filter(valid_bearings, treatment == "ON")
gridon <- gridcount(gridcount_ON)
gridcount_ON <- ggplot(gridcount_ON, aes(x, y)) + 
  stat_bin_2d(binwidth = 10) + 
  geom_point(color = "white", cex = 0.25) +
  xlim(0, 250) +
  labs(title = "TAST ON")

gridcount_ON

#plot off
gridcount_OFF <- filter(valid_bearings, treatment == "OFF")
gridoff <- gridcount(gridcount_OFF)
gridcount_OFF <- ggplot(gridcount_OFF, aes(x, y)) + 
  stat_bin_2d(binwidth = 10, show.legend = F) + 
  geom_point(color = "white", cex = 0.25)+
  labs(title = "TAST OFF")
gridcount_OFF


# plot on and off side by side
require(gridExtra)
grid.arrange(gridcount_OFF, gridcount_ON, ncol = 2)

### TODO: how can I add counts for multi-individual surfacings?????

# rearrange pred data set so binary presence/absence for each 10x10m square 
library(reshape2)

#list all unique Session IDs to loop through for valid bearings
session_IDs <- unique(valid_bearings$session_id)

# create empty vector to fill with loop output
datalist = list()

# This loop reates a unique ID for each session + x + y combo with an observation
for (ID in session_IDs) {
  # subset the particular survey session
  df_sub <- subset(valid_bearings, session_id == ID)
  # for that session, list x y square for each observation
  a <- gridcount(df_sub)
  # create unique ID for each session + x + y  combo
  # b <- expand.grid(session_IDs, unique(a$cut_x), unique(a$cut_y))
  # save in new column of a
  datalist[[ID]] <- paste(a$cut_x, a$cut_y, ID, sep = "_")
}
datalist
# turn list into a data frame
datalisst <- plyr::ldply(datalist, rbind)
#group by session make each row an observation
uniqueIDs <- datalisst %>%
  pivot_longer(!.id, names_to = "count", values_to = "unique_ID")

# expand grid to create list of every possible combo of session ID and square
head(grid_counts)

all_combo <- expand.grid(unique(grid_counts$cut_x), unique(grid_counts$cut_y), session_IDs)
head(all_combo)

#create 4th column with unique ID
all_combo$Var4 <- paste(all_combo$Var1, all_combo$Var2, all_combo$Var3, sep = "_")

# Create PA column
#loop through all combo, if value exists in unique ids, put a 1, if not, put a 0
Df1 <- all_combo
Df2 <- data.frame(uniqueIDs)

Df1$PA <- as.integer(Df1$Var4 %in% Df2$unique_ID)
head(Df1)

# add column for treatment
library(stringr)
Df1$treatment <- ifelse(str_sub(Df1$Var3, 1, 1) == "T", paste("ON"), paste("OFF"))

# add column for square ID
Df1$square_ID <- paste(Df1$Var1, Df1$Var2, sep = "_")

# fix "Var" columns
colnames(Df1) <- c("cut_x", "cut_y", "session_id", "unique_id", "p_a", "treatment", "square_ID")

# add a group size and pred rate column by matching session ID
pred$session_id <- pred$Session_ID
pred3 <- pred %>%
  select(session_id, Seal_Freq, pred_rate_per_hour, Observer1, Date, Julian_day_from_start, hour, dur_hr)
Df1 <- left_join(Df1, pred3, by = "session_id")

# add fish count 
dist2 <- distance%>%
  select(session_id, fish_day)
Df1 <- left_join(Df1, dist2, by = "session_id")

## DF1 is a data frame where P/A is recorded for each square during each survey
Df1$treatment <- as.factor(Df1$treatment)
head(valid_bearings)

# make x column from cut_x so we are dealing with a single value for each square
xsq <- regmatches(Df1$cut_x ,gregexpr("(?<=[\\[\\(,])[0-9.]+(?=[\\]\\),])", Df1$cut_x, perl = TRUE))
xsq <- sapply(xsq, function(r) r[2])
Df1$xsq <- as.numeric(xsq)

# make y column # using regular expressiongs (YUCK) 
ysq <- regmatches(Df1$cut_y ,gregexpr("(?<=[\\[\\(,])[0-9.]+(?=[\\]\\),])", Df1$cut_y, perl = TRUE))
ysq <- sapply(ysq, function(x) x[2])
Df1$ysq <- as.numeric(ysq)
head(Df1)

# okay now for every PA = 1 bring in Distance and Bearing rl x y etc
PA1 <- Df1%>%
  filter(p_a == 1)
head(PA1)
PA1 <- distinct(PA1, unique_id, .keep_all = TRUE)

# add this new rl information to DF1
################################################################################
########### okay bring in clean distance data to validate + help from abinand
################################################################################

#filter clean distance data for only days after 02 sep or jul day 39
clean_dist2 <- read.csv("cleandist.csv") 
valid_bearings <-  filter(clean_dist2, jul_day > 39)

numon <- length(which(valid_bearings$treatment == "ON"))

# make sure df is in ascending order by datetime so 1st index matches refcoord
valid_bearings$datetime <- as_datetime(paste(valid_bearings$Date, valid_bearings$Time, sep = " "))
valid_bearings <- arrange(valid_bearings, datetime)

## great package for spatial stuff
library(sf)
library(sp)
library(rgdal)

## Read in shape file
poly <- st_read("BSA3-polygon.shp")

## convert from lat-long to UTM
poly1 <- st_transform(poly,crs = "+proj=utm +zone=10")

## get coordinates
poly_coord <- st_coordinates(poly1)

## get our reference coordinates 
origin <- c(-122.39772,47.664821)
ref <- c(-122.398166, 47.665465) 
ref_pts <- data.frame(x = c(origin[1],ref[1]), y = c(origin[2],ref[2]))

##convert to spatial object
coordinates(ref_pts) <- ~x+y

##convert to UTM as done with polygon
ref_pts <- ref_pts %>% st_as_sf() %>% st_set_crs(st_crs(poly)) %>% st_transform("+proj=utm +zone=10") %>% st_coordinates()

##Redefine the origin as in our dataset
poly_coord <- poly_coord %>% data.frame() %>% mutate(X = X - ref_pts[1,1], Y = Y - ref_pts[1,2]) %>% as.matrix()

####Check distances match
ref_coord <- c(ref_pts[2,]-ref_pts[1,])
sqrt(sum(ref_coord^2)) 
valid_bearings$tide_dist[1]

#### There seems to be a small discrepancy of about 2cm. Check if this will be 
# an issue, otherwise looks good.

#### We've got the shapefile to the same scale as our data but they aren't 
# aligned yet. Need to rotate one of the two to make them match

#### Rotating the polygon to the dataset, might want to swap this at some 
# point if I want to make plots for the report for geographic consistency

###Target vector to align to
target_coord <- c(valid_bearings$x[1],valid_bearings$y[1])

###Find angle to rotate by
theta <- acos(sum(ref_coord*target_coord)/(sqrt(sum(ref_coord^2))*sqrt(sum(target_coord^2))))

###Create a rotation matrix
rotmat <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow = 2,byrow = T)

###Rotate the polygon by matrix multiplication
as.matrix(poly_coord)
poly_coord <- poly_coord[,1:2]%*%rotmat %>% data.frame()

### Convert the coordinates to a polygon again
polygon <- poly_coord %>%
  st_as_sf(coords = c("X1", "X2")) %>%
   mutate(geometry = st_combine(geometry)) %>% 
  st_cast("POLYGON")

###Original polygon with origin and reference coordinate
ggplot()+
  geom_sf(data = poly)+
  geom_point(aes(x = origin[1],y = origin[2]))+
  geom_point(aes(x = ref[1],y =ref[2]))

### Rotated polygon in our new coordinate system. looks ok to me
ggplot()+
  geom_sf(data = polygon)+
  geom_point(aes(x = 0,y = 0))+
  geom_point(aes(x = target_coord[1], y= target_coord[2]))

## Check other points ----  
valid_bearings$x <- as.numeric(valid_bearings$x)
valid_bearings$y <- as.numeric(valid_bearings$y)

### TODO: most points line up. clean up last of weird points when Katie get the angles to me 
# plot of final points layout 
new <- ggplot()+
  geom_sf(data = polygon)+
  geom_point(data = valid_bearings, aes(x = x, y = y, color = blockid))+
  geom_text(data = valid_bearings, position=position_jitter(width=5,height=5), aes(x = x, y = y, label = ID))
new
old <- ggplot()+
  geom_sf(data = polygon)+
  geom_point(data = valid_bearings, aes(x = x, y = y, color = blockid))+
  geom_text(data = valid_bearings, position=position_jitter(width=5,height=5), aes(x = x, y = y, label = ID))
old

################################################################################
## exploring angle errors by observer

asila <- valid_bearings%>% filter(tide_dist <= 250)%>%filter(obs == "AB")

asila_plot <- ggplot()+
  geom_sf(data = polygon)+
  geom_point(data = asila, aes(x = x,y = y))+
  labs(title = "Asila")+ 
  theme(plot.title = element_text(size=10))

amina <- valid_bearings%>% filter(tide_dist <= 250)%>%filter(obs == "AC")

amina_plot <- ggplot()+
  geom_sf(data = polygon)+
  geom_point(data = amina, aes(x = x,y = y))+
  labs(title = "Amina")+ 
  theme(plot.title = element_text(size=10))

drea <- valid_bearings%>% filter(tide_dist <= 250)%>%filter(obs == "AMB")

drea_plot <- ggplot()+
  geom_sf(data = polygon)+
  geom_point(data = drea, aes(x = x,y = y))+
  labs(title = "Andrea")+ 
  theme(plot.title = element_text(size=10))

laura <- valid_bearings%>% filter(tide_dist <= 250)%>%filter(obs == "LB")

laura_plot <- ggplot() +
  geom_sf(data = polygon) +
  geom_point(data = laura, aes(x = x, y = y)) +
  labs(title = "Laura") + 
  theme(plot.title = element_text(size=10))

cowplot::plot_grid(laura_plot, amina_plot, asila_plot, drea_plot)
################################################################################
### Create grid. I've chosen the centres to be at multiples of 10. Just seems easier 
### you can change the offset to c(-10,-40) if you want the centres at multiples 
### of 5 and grid edges to be at multiples of 10

gr <- st_make_grid(polygon,cellsize = 10,offset=c(-10,-40))%>% st_sf()

##assign id based on centre
gr$id <- apply(st_coordinates(st_centroid(gr)),1,paste,collapse = ',')

#### Select grids that intersects with the polygon
gr1 <- gr[st_intersects(polygon,gr,sparse = F),]

ggplot()+
  geom_sf(data = polygon)+
  geom_sf(data = gr1, fill = NA) #ask abinand about this 

##### We don't want grid areas that are outside so instead clip grids and select above certain area? 
##### I've retained any grid above 1m square, if this is too small still, you can change it to something larger

gr2 <- st_intersection(gr,polygon)

gr2 <- gr2[st_area(gr2)>1,] %>% st_as_sf() #### Ask steve about this 


###Some of the grids are not exactly sqaures at the point. Use centroids of the shapes as the centres for grid values
centres <- st_centroid(gr2)

### calculate spatial covariates 
gr2 <- gr2 %>% 
  mutate(distance = as.numeric(st_distance(centres,st_point(x = c(0,0))))) %>% 
  mutate(RL_theoretical = as.numeric(182 - (20*log10(distance)))) %>% 
  mutate(Emperical_RL = c(predict(rlmod, newdata = data.frame(rldist = distance))))
## RL PLOT
##### 
# plot emperical RL over study area 
RLPLOT <- ggplot()+
  geom_sf(data = gr2,aes(fill = Emperical_RL))+
  theme_classic()+
 # points(0,0, pch= "*", col="green", cex=3)+ dont work 
  theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10), 
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.key.size = unit(5, "mm"), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(x = "x position (metres)", title = "Empirically Modelled Transmission Loss in Study Area", y = "y position (metres)", fill = "Emperical RL (dB re 1µPa)")

RLPLOT

RL_plot <- ggplot(data = RL) +
  geom_point(mapping = aes(x = rldist, y = mean_rl)) +
  ggtitle("Transmission Loss of TAST signal") + 
  labs(x = "Distance (m)", y = "Received level dB re 1µPa", linetype = "Model") +
  stat_smooth(mapping = aes(x = rldist, y = mean_rl, linetype = "Empirical"),
              colour = "black", method = "lm", formula = y ~ I(20*log(x)),  se=T, level=0.95) +
  stat_function(fun = function(rldist){182 - 20*log10(rldist)}, aes(linetype = "Spherical")) +
  stat_function(fun = function(rldist){182 - 10*log10(rldist)}, aes(linetype = "Cylindrical")) +
  theme_light() +
  theme(legend.position = c(0.85, 0.8), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.key.size = unit(12, "mm"), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(0, 5000, by = 50)) +
  scale_y_continuous(breaks = seq(100, 180, by = 10)) #+
  #scale_linetype_manual(values=c( "Cylindrical" ='dashed', "Spherical" = 'dotted', "Empirical" = 'solid')) 



RL_plot



# # add an RL column to massive distance data set
# distance <-  distance %>%
#   mutate(RL_emperical = c(predict(rlmod, newdata = data.frame(rldist = tide_dist))))
#####
####Once you have your data cleaned, there are two ways to go ahead.

####1. Since our grids are intervals of 10, we can arithmetically assign the respective grid ids

valid_bearings1 <- valid_bearings %>% mutate(gridx = ((x+5)%/%10)*10,
                                             gridy = ((y+5)%/%10)*10,
                                             id = paste(gridx,gridy,sep = ','))

### the problem with this is that there are lots of points that fall outside the 
# actual grid we've created and this doesn't recognise that

####2. Instead we assign the grid ids using a spatial operator. 


valid_bearings2 <- valid_bearings


valid_bearings2 %>% mutate(x = jitter(x,.0001),y = jitter(y,.0001))

valid_bearings2 <- valid_bearings2 %>% filter(!is.na(x),!is.na(y)) %>% 
  mutate(X = x, Y = y)%>%
  st_as_sf(coords = c('x','y'))

ints <- st_intersects(valid_bearings2, gr2)

is.na(ints) <- lengths(ints) == 0

valid_bearings2$id <- gr2$id[unlist(lapply(ints,'[',1))]
#valid_bearings2$id <- gr2$id[unlist(ints)]

grid_obs <- ggplot()+
  geom_sf(data = polygon)+
  geom_sf(data = gr2)+
  geom_sf(data = valid_bearings2[!is.na(valid_bearings2$id),])+
  theme_classic() +
  labs(x = "x position (metres)", y = "y position (metres)", title = "Raw Observation Data")+
theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10), 
  plot.title = element_text(hjust = 0.5, size = 14),
  legend.key.size = unit(5, "mm"),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12))

### matched grid ids, now getting counts and PA 
valid_bearings2$session_id <- factor(valid_bearings2$session_id)

###We make the ids a factor too making all grids within the survey area into a 
#factor. This helps with the tabulate absences where there were no detections at all
valid_bearings2$id <- factor(valid_bearings2$id, levels = unique(gr2$id))

###Note that this does not take the number of individuals into account. 
# Just the number of detections in a grid per session
final_data <- data.frame(table(valid_bearings2$session_id,valid_bearings2$id))

colnames(final_data) <- c('session_id','id','counts')

#add presence absence
final_data$pa <- ifelse(final_data$counts > 0,1,0)

### Now I have session id to match session/temporal covariates and grid id to match spatial covariates, 

#create spattial covar df
spatcovar <- gr2
spatcovar <- as.data.frame(st_set_geometry(spatcovar, "geometry"))
spatcovar <- select(spatcovar, -geometry)
#rename some columns
spatcovar$dist2sq <- spatcovar$distance
spatcovar$rl_the_sq <-  spatcovar$RL_theoretical
spatcovar$rl_emp_sq <- spatcovar$RL_emperical
spatcovar <- select(spatcovar, -distance, -RL_theoretical, -RL_emperical)
class(spatcovar)

# add unique ID key to tie together session id and id
final_data$key <- paste(final_data$session_id, final_data$id, sep = "_")
valid_bearings2$key <- paste(valid_bearings2$session_id, valid_bearings2$id, sep = "_")

# great now join with final any other spatial covar?
final_data <- left_join(final_data, spatcovar, by = "id")
final_data <- distinct(final_data)
head(final_data)

#convert back to df
valid_bearings2 <- data.frame(valid_bearings2)
names(valid_bearings2)
names(final_data)

#  make new df with session covariates
sessioncovar <- valid_bearings2 %>%
  select(session_id, jul_day, start_hr, dur_hr, treatment, obs, fish) %>% 
  group_by(session_id) %>% 
  filter(row_number() == 1)

final_data <- left_join(final_data, sessioncovar, by = "session_id")

#seperate id into x, y

final_data <- final_data %>% separate(id, sep = ",", c("x", "y"), remove = F)

# write csv 
write.csv(final_data, "~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv") 

