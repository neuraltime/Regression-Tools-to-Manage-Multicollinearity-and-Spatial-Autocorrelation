#### Regression-Tools-to-Manage-Multicollinearity-and-Spatial-Autocorrelation
#### This study looks at the challenges of using linear regression to analyze geospatial data.
### R Code fore ACS Data for Sonoma, Marin, Mendocino, Lake, and Napa Counties
####  Taken from Kyle Walker
####  https://walker-data.com/census-r/modeling-us-census-data.html
### 1 Install Packages
install.packages("tidyverse")
install.packages("stats")
install.packages("plotly")
install.packages("tidycensus")
install.packages("tidyBLS")
install.packages("bea.R")
install.packages("segregation")
install.packages("patchwork")
install.packages("ggfortify")
install.packages("corrr")
install.packages("car")
install.packages("spdep")
install.packages("spatialreg")
install.packages("GWmodel")
### 2 Loading Packages
library(tidyverse)
library(tidycensus)
library(tmap)
library(segregation)
library(tigris)
library(sf)
library(patchwork)
library(stats)
library(ggfortify)
library(car)
library(spdep)
library(units)
library(corrr)
library(spatialreg)
library(GWmodel)

##########################
######library(bea.R)
##########################
#### census_api
##########################

nbay_counties <- c("Mendocino", "Lake", "Sonoma", 
                  "Marin", "Napa")
                  
variables_to_get_nbay <- c(
  Agg_Income_19313 = "B19313_001", 
  Agg_HH_Income_19025 = "B19025_001",
  median_value = "B25077_001",
  median_rooms = "B25018_001",
  median_income = "DP03_0062",
  total_population = "B01003_001",
  median_age = "B01002_001",
  pct_college = "DP02_0068P",
  pct_foreign_born = "DP02_0094P",
  pct_white = "DP05_0077P",
  median_year_built = "B25037_001",
  Agg_Earnings_for_HH = "B19061_001", 
  Agg_Wage_Salary_Inc = "B19062_001", 
  Agg_Social_Sec_Inc = "B19065_001", 
  TotPopulation = "B01003_001",
  TotPop_for_groups = "B03002_001",
  PopFill = "B03002_002",
  White = "B03002_003",
  Black = "B03002_004",
  NativeAmer = "B03002_005",
  Asian = "B03002_006",
  PacificIsland = "B03002_007",
  OtherRace = "B03002_008",
  Hispanic = "B03002_012",
  Housing_Tot = "B07013_001",
  Owner_occupied = "B07013_002",
  Renter_occupied = "B07013_003",
  Same_House_1yr_Past = "B07013_004",
  Owner_Same_House_1yr_Past = "B07013_005",
  percent_ooh = "DP04_0046P"
)

nbay_data_science <- get_acs(
  geography = "tract",
  variables = variables_to_get_nbay,
  state = "CA",
  county = nbay_counties,
  geometry = TRUE,
  output = "wide",
  survey = "acs5",
  year = 2021
) %>%
  select(-NAME) %>%
  st_transform(2226) %>% st_transform(2227) # NAD83 / Northern California

library(tidyverse)
library(patchwork)

###############################################################################
### Clean Data
###############################################################################

options(max.print = 100000)
glimpse(nbay_data_science)

summary(nbay_data_science)

nbay_nacheck <- is.na(nbay_data_science)
##################################################
########### Write table code
##################################################

write.table(nbay_nacheck,
            file = "nbay_nacheck2.csv",
            sep = "\t",
            row.names = FALSE,
)
#####################################################################

##########################################################################################
#########  Map Data Median Value ######################
library(sf)
library(units)

nbay_data_science_for_model <- nbay_data_science %>%
  mutate(pop_density = as.numeric(set_units(total_populationE / st_area(.), "1/km2")),
         median_structure_age = 2018 - median_year_builtE,
         Agg_Income_19313_person = Agg_Income_19313E / total_populationE,
         Wage_to_Social = Agg_Wage_Salary_IncE / Agg_Social_Sec_IncE) %>%
  select(!ends_with("M")) %>% 
  rename_with(.fn = ~str_remove(.x, "E$")) %>%
  na.omit()

formula <- "log(median_value) ~ median_rooms + median_income + pct_college + pct_foreign_born + pct_white + median_age + median_structure_age + percent_ooh + pop_density + total_population +Wage_to_Social+ Hispanic"

model1 <- lm(formula = formula, data = nbay_data_science_for_model)

summary(model1)

library(corrr)

nbay_estimates_science <- nbay_data_science_for_model %>%
  select(-GEOID, -median_year_built) %>%
  st_drop_geometry()

correlations <- correlate(nbay_estimates_science, method = "pearson")

network_plot(correlations)
print(correlations)

library(car)

vif(model1)

##################################################################
### Model 10
#################################################################


formula10 <- "log(median_value) ~ median_rooms+ median_income + pct_college + median_age + percent_ooh"

model10 <- lm(formula = formula10, data = nbay_data_science_for_model)

summary(model10)



##################################################################
### Model 10 validity check
#################################################################
vif(model10)

##################################################################
### Model 2
#################################################################


formula2 <- "log(median_value) ~ median_income + pct_college + median_age + percent_ooh"

model2 <- lm(formula = formula2, data = nbay_data_science_for_model)

summary(model2)



##################################################################
### Model 2 validity check
#################################################################
vif(model2)
#########################################################################
#### A core assumption of the linear model is that the errors are independent 
#### of one another and normally distributed. we can check this by adding the 
#### residuals for model2 to our dataset and drawing a histogram
#########################################################################
nbay_data_science_for_model$residuals <- residuals(model2)

ggplot(nbay_data_science_for_model, aes(x = residuals)) + 
  geom_histogram(bins = 100, alpha = 0.5, color = "navy",
                 fill = "navy") + 
  theme_minimal()

library(spdep)

#########################################################################
#### Moran's I - test residuals
#########################################################################

wts <- nbay_data_science_for_model %>%
  poly2nb() %>%
  nb2listw()

moran.test(nbay_data_science_for_model$residuals, wts)

######################################################################
### The Moran's Itest statistic is modest and positive (0.21) but is statistically significant.
###This can be visualized with a Moran scatterplot:
####################################################################
nbay_data_science_for_model$lagged_residuals <- lag.listw(wts, nbay_data_science_for_model$residuals)

ggplot(nbay_data_science_for_model, aes(x = residuals, y = lagged_residuals)) + 
  theme_minimal() + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", color = "red")

###################################################################
### The plot illustrates the positive spatial autocorrelation in the residuals,
### suggesting that the assumption of independence in the model error term is violated. 
### To resolve this issue, we can turn to spatial regression methods.
###################################################################
library(spatialreg)

lag_model <- lagsarlm(
  formula = formula2, 
  data = nbay_data_science_for_model, 
  listw = wts
)
summary(lag_model, Nagelkerke = TRUE)

###################################################################
### Spatial Error Models
#################################################################
error_model <- errorsarlm(
  formula = formula2, 
  data = nbay_data_science_for_model, 
  listw = wts
)

summary(error_model, Nagelkerke = TRUE)

###################################################################
### Moran test of Spatial lag error model
#################################################################
moran.test(lag_model$residuals, wts)

###################################################################
### Moran test of Spatial error error model
#################################################################
moran.test(error_model$residuals, wts)

###################################################################
### The lm.LMtests() function can be used with an input linear model to compute these tests. 
### We'll use model2, the home value model with median household income omitted, 
### to compute the tests.
###################################################################
lm.LMtests(
  model2, 
  wts, 
  test = c("LMerr", "LMlag", "RLMerr", "RLMlag")
)

###################################################################
### Geographically weighted regression
### Test  It is possible that a relationship between a predictor and
### the outcome variable that is observed for the entire region on average
### may vary significantly from neighborhood to neighborhood. 
### This type of phenomenon is called spatial non-stationarity, and can be explored
### with geographically weighted regression, or GWR (Brunsdon, Fotheringham,
### and Charlton 1996).
###################################################################
libary(GWmodel)

nbay_data_science_sp <- nbay_data_science_for_model %>%
  as_Spatial()

bw <- bw.gwr(
  formula = formula2, 
  data = nbay_data_science_sp, 
  kernel = "bisquare",
  adaptive = TRUE
)
###################################################################
### Fitting and evaluating the GWR model
###################################################################
formula2gw <- "log(median_value) ~ median_rooms + pct_college + pct_foreign_born + pct_white + median_age + median_structure_age + percent_ooh + pop_density + total_population"

gw_model <- gwr.basic(
  formula = formula2gw, 
  data = nbay_data_science_sp, 
  bw = bw,
  kernel = "bisquare",
  adaptive = TRUE
)

names(gw_model)

gw_model_results <- gw_model$SDF %>%
  st_as_sf() 

names(gw_model_results)

###################################################################
### Below, we use ggplot2 and geom_sf() to map the local R-squared,
#####################################################################
ggplot(gw_model_results, aes(fill = Local_R2)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c() + 
  theme_void()


#######################################################################
###  Recall from the global model that this coefficient was negative 
###  and statistically significant.
#######################################################################
ggplot(gw_model_results, aes(fill = percent_ooh)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c() + 
  theme_void() + 
  labs(fill = "Local Percent of Occupants > 1 year")

#######################################################################
### The dark purple areas on the map are those areas where the global relationship
###  in the model reflects the local relationship, as local parameter estimates are negative.
#################################################################################
ggplot(gw_model_results, aes(fill = pop_density)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c() + 
  theme_void() + 
  labs(fill = "Local population density")

