#Learn GIS

library(sf)
library(GISTools) # a wrapper for sp, rgeos etc
# load some data
data(georgia)
class(georgia) #sp
# convert to sf
georgia_sf <- st_as_sf(georgia)
class(georgia_sf)
# convert back to sp
georgia_v2 <- as(georgia_sf, "Spatial")
class(georgia_v2)

##spatial intersection
data(tornados)

library(tmap)
library(sf)
# convert to sf objects
torn_sf <- st_as_sf(torn)
us_states_sf <- st_as_sf(us_states)
# plot extent and grey background
tm_shape(us_states_sf) +
  tm_polygons("grey90") +
  # add the torn points 
  tm_shape(torn_sf) +
  tm_dots(col = "#FB6A4A", size = 0.04, shape = 1, alpha = 0.5) +
  # map the state borders
  tm_shape(us_states_sf) +
  tm_borders(col = "black") +
  tm_layout(frame = F) 

#point in polygon counts
torn.count <- poly.counts(torn, us_states) 
head(torn.count)

#area calculations
proj4string(us_states)
st_crs(us_states_sf) #Retrieve coordinate reference system from object
poly.areas(us_states2) #Given a set of polygons, returns the area of each polygon
st_area(us_states_sf) #computes the area for a set of geometries
# hectares
poly.areas(us_states2) / (100 * 100)
st_area(us_states2_sf) / (100 * 100)
# square kilometres
poly.areas(us_states2) / (1000 * 1000)
st_area(us_states2_sf) / (1000 * 1000)


#######

# load packages
library(sf)
library(dplyr)
library(ggplot2)
library(scico)
library(rnaturalearth)
library(purrr)
library(smoothr)
library(rgbif)


# world map
worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")

# country subset
CRpoly <- worldMap %>% filter(sovereignt == "Costa Rica")

# random points tables as a named list
sp_occ <- rerun(12, st_sample(CRpoly, sample(3:20, 1)))
names(sp_occ) <- paste0("sp_", letters[1:length(sp_occ)])

# to sf object
sflisss <-
  purrr::map(sp_occ, st_sf) %>%  #map a list to create a st_sf
  map2(., names(.), ~ mutate(.x, id = .y))

sp_occ_sf <- sflisss %>% reduce(rbind) #combine into one dataframe

# to multipoint
sp_occ_sf <- sp_occ_sf %>%
  group_by(id) %>%
  summarise()

# trim to study area
limsCR <- st_buffer(CRpoly, dist = 0.7) %>% st_bbox()  #st_bbox() bounds box

sf_use_s2(FALSE)
# neighboring countries
adjacentPolys <- st_touches(CRpoly, worldMap)
neighbours <- worldMap %>% slice(pluck(adjacentPolys, 1))

# countries
divpolPlot <-
  ggplot() +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
divpolPlot

# plot points
spPointsPlot <-
  ggplot() +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly) +
  geom_sf(data = sp_occ_sf, aes(fill = id), pch = 21) +
  scale_fill_scico_d(palette = "davos", direction = -1, end = 0.9, guide = FALSE) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
spPointsPlot

# convex hulls
spEOOs <- st_convex_hull(sp_occ_sf) %>% smooth() #creates the convex hull of a set of points

# plot hulls
hullsPlot <-
  ggplot() +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly) +
  geom_sf(data = spEOOs, aes(fill = id), alpha = 0.7) +
  scale_fill_scico_d(palette = "davos", direction = -1, end = 0.9, guide = FALSE) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
hullsPlot

#### create a grid to calculate species richness
# grid
CRGrid <- CRpoly %>%
  st_make_grid(cellsize = 0.2) %>%
  st_intersection(CRpoly) %>%
  st_cast("MULTIPOLYGON") %>%
  st_sf() %>%
  mutate(cellid = row_number())

# cell richness
richness_grid <- CRGrid %>%
  st_join(sp_occ_sf) %>%
  mutate(overlap = ifelse(!is.na(id), 1, 0)) %>%  #if its not NA return 1, if else return 0
  group_by(cellid) %>%
  summarize(num_species = sum(overlap))

# richness for convex hulls
richness_gridEOO <- CRGrid %>%
  st_join(spEOOs) %>%
  mutate(overlap = ifelse(!is.na(id), 1, 0)) %>%
  group_by(cellid) %>%
  summarize(num_species = sum(overlap))

# empty grid
gridPlot <-
  ggplot() +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly) +
  geom_sf(data = CRGrid) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
gridPlot

# richness
gridRichCR <-
  ggplot(richness_grid) +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly, fill = "grey", size = 0.1) +
  geom_sf(aes(fill = num_species), color = NA) +
  scale_fill_scico(palette = "davos", direction = -1, end = 0.9) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  ) + labs(fill = "richness")
gridRichCR

# richness for convex hulls
gridRichCR_eoo <-
  ggplot(richness_gridEOO) +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly, fill = "grey", size = 0.1) +
  geom_sf(aes(fill = num_species), color = NA) +
  scale_fill_scico(palette = "davos", direction = -1, end = 0.9) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  ) + labs(fill = "richness")
gridRichCR_eoo

###################################################
# bat data
name_suggest(q = "chiroptera")
CRbatsout <- occ_search(
  orderKey = 734, country = "CR",
  basisOfRecord = "PRESERVED_SPECIMEN", limit = 3000
)$data
# species and lat long data
CRbatsXY <- CRbatsout %>%
  dplyr::select(species, decimalLongitude, decimalLatitude) %>%
  na.omit()
# to sf object, specifying variables with coordinates and projection
CRbatsXYsf <- st_as_sf(CRbatsXY, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  group_by(species) %>%
  summarize()

# plot points
batPointsPlot <-
  ggplot() +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly) +
  geom_sf(data = CRbatsXYsf, pch = 21) +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
batPointsPlot

# bat richness
bat_richness_grid <- CRGrid %>%
  st_join(CRbatsXYsf) %>%
  mutate(overlap = ifelse(!is.na(species), 1, 0)) %>%
  group_by(cellid) %>%
  summarize(num_species = sum(overlap))

# plot
batRichCR <-
  ggplot(bat_richness_grid) +
  geom_sf(data = neighbours, color = "white") +
  geom_sf(data = CRpoly, fill = "grey", size = 0.1) +
  geom_sf(aes(fill = num_species), color = NA) +
  scale_fill_scico(palette = "davos", direction = -1, end = 0.9, name = "Bat species richness") +
  coord_sf(
    xlim = c(limsCR["xmin"], limsCR["xmax"]),
    ylim = c(limsCR["ymin"], limsCR["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#f1f2f3"),
    panel.background = element_rect(fill = "#2F4051"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    line = element_blank(),
    rect = element_blank()
  ) + labs(fill = "richness")
batRichCR
