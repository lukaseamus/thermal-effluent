# Impacts of thermal effluent on *Posidonia oceanica* and associated macrofauna
This repository contains data and annotated R code accompanying article 10.3354/meps14261 in *Marine Ecology Progress Series*. The repository containes the five data and code folders **Physchem**, **Seagrass sub-plant**, **Seagrass plant**, **Seagrass habitat** and **Fauna** alongside `map.R` and all final **Figures** in vector format. Below is a description of each column within each dataset and the input and output of each R script.

`map.R`: Code to build and visualise the map of Malta and the study site (the Mediterranean map was built in QGIS).
- **Input** = `land_polygons.shp`, `WCMC013014-Seagrasses-Py-v7_1.shp`, which can be downloaded from https://osmdata.openstreetmap.de/data/land-polygons.html and https://data.unep-wcmc.org/datasets/7 respectively
- **Output** = Figure 1

**Physchem**
1. `physchem.csv`: Temperature, pH and salinity data for the study sites.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Month** = categorical variable with levels February, June, July and September
    - **Position** = categorical variable with levels Surface and Meadow, indicating the position of measurment and enabling detection of clines
    - **Depth** = continuous alternative to **Position** given in metres below the sea surface 
    - **Distance** = distance from the thermal outlet given in metres
    - **Salinity** = salinity (‰)
    - **pH** = pH
    - **Temperature** = temperature (°C)
2. `physchem.R`: Code to analyse and visualise temperature, pH and salinity data.
    - **Input** = `physchem.csv`
    - **Output** = Figure 2 and section 3.1 in `Manuscript.pdf`

**Seagrass sub-plant**
1. `seagrass.subplant.csv`: *Posidonia oceanica* metrics at the sub-plant level, where plant is defined as a single shoot.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Plot** = categorical variable of sampling plots nested within sites with levels A, B, C and D in order of decreasing distance
    - **Distance** = distance from the thermal outlet given in metres, a continuous reexpression of **Plot**
    - **Stage** = categorical variable of within-shoot leaf growth stage with levels Juvenile, Intermediate and Adult
    - **Leaves** = number of leaves of the respective growth stage per shoot
    - **Area** = leaf area of the respective growth stage given in square centimetres per shoot
    - **Biomass** = leaf dry mass of the respective growth stage given in grams per shoot
    - **Density** = shoot density (all growth stages combined) per plot given as number of shoots per square metre
2. `seagrass.subplant.R`: Code to analyse and visualise *Posidonia oceanica* metrics at the sub-plant level.
    - **Input** = `seagrass.subplant.csv`
    - **Output** = Figures S1 and S2 in `Supplement.pdf` and parts of section 3.2 in `Manuscript.pdf`

**Seagrass plant**
1. `seagrass.plant.csv`: *Posidonia oceanica* metrics at the plant level, where plant is defined as a single shoot.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Plot** = categorical variable of sampling plots nested within sites with levels A, B, C and D in order of decreasing distance
    - **Distance** = distance from the thermal outlet given in metres, a continuous reexpression of **Plot**
    - **Plant** = categorical variable with levels Seagrass (*Posidonia oceanica*) and Epiphytes (algae growing on *Posidonia oceanica*)
    - **Leaves** = total number of leaves per shoot
    - **Area** = total leaf area given in square centimetres per shoot
    - **Biomass** = total leaf dry mass given in grams per shoot
    - **Shoots** = shoot density per plot given as number of shoots per square metre (same as **Density** in `seagrass.subplant.csv`)
2. `seagrass.plant.R`: Code to analyse and visualise *Posidonia oceanica* metrics at the plant level.
    - **Input** = `seagrass.plant.csv`
    - **Output** = Figure 3 and most of section 3.2 in `Manuscript.pdf`
    
**Seagrass habitat**
1. `seagrass.habitat.csv`: Areal *Posidonia oceanica* metrics at the habitat level, calculated by multiplying mean (n = 3) shoot density per square metre of meadow for a given plot (see `seagrass.subplant.csv` or `seagrass.plant.csv`) by the other variables per shoot in the same plot.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Plot** = categorical variable of sampling plots nested within sites with levels A, B, C and D in order of decreasing distance
    - **Distance** = distance from the thermal outlet given in metres, a continuous reexpression of **Plot**
    - **Plant** = categorical variable with levels Seagrass (*Posidonia oceanica*) and Epiphytes (algae growing on *Posidonia oceanica*)
    - **Leaves** = total number of leaves per square metre
    - **Area** = total leaf area given in square metres per square metre of meadow
    - **Biomass** = total leaf dry mass given in grams per square metre of meadow
2. `seagrass.habitat.R`: Code to analyse and visualise areal *Posidonia oceanica* metrics at the habitat level.
    - **Input** = `seagrass.habitat.csv`
    - **Output** = Figure S3 in `Supplement.pdf` and parts of section 3.2 in `Manuscript.pdf`

**Fauna**
1. `fauna.csv`: Species-abundance matrix of macrofauna associated with *Posidonia oceanica*.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Plot** = categorical variable of sampling plots nested within sites with levels A, B, C and D in order of decreasing distance
    - **Distance** = distance from the thermal outlet given in metres, a continuous reexpression of **Plot**
    - **Label** = unique sample ID fusing **Site** and **Plot** with the replicate number
    - the remaining columns of the dataframe are species (or the lowest identifiable unique taxon above the species level) and each row is the number of individuals of the given species in the sample (captured by a net dragged through 24 square metres of meadow)
2. `fauna.R`: Code to analyse and visualise the species-abundance matrix of macrofauna associated with *Posidonia oceanica*.
    - **Input** = `fauna.csv`
    - **Output** = Tables S1 and S2 in `Supplement.pdf` and Figures 4 and 5 and section 3.3 in `Manuscript.pdf`

Luka Seamus Wright, 19 March 2023
