# Impacts of thermal effluent on *Posidonia oceanica* and associated macrofauna
This repository contains data and annotated R code accompanying article 10.3354/meps14261 in *Marine Ecology Progress Series*. The repository containes the five data and code folders **Physchem**, **Seagrass sub-plant**, **Seagrass plant**, **Seagrass habitat** and **Fauna** alongside `map.R` and all final **Figures** in vector format. Below is a description of each column within each dataset and the input and output of each R script.

`map.R`: Code to build the map of Malta and the study site (the Mediterranean map was built in QGIS).
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
2. `physchem.R`: Code to analyse temperature, pH and salinity data.
    - **Input** = `physchem.csv`
    - **Output** = Figure 2 and section 3.1 in `Manuscript.pdf`

**Seagrass sub-plant**
1. `seagrass.subplant.csv`: *Posidnia oceanica* metrics at the sub-plant level, where plant is defined as a single shoot.
    - **Site** = categorical variable with levels Kbira (short for Il-Ħofra l-Kbira, control site) and Zghira (short for Il-Ħofra ż-Żgħira, treatment site)
    - **Plot** = categorical variable of sampling plots nested within sites with levels A, B, C and D in order of decreasing distance
    - **Distance** = distance from the thermal outlet given in metres, a continuous reexpression of **Plot**
    - **Stage** = categorical variable of within-shoot leaf growth stage with levels Juvenile, Intermediate and Adult
    - **Leaves** = number of leaves of the respective growth stage per shoot
    - **Area** = leaf area of the respective growth stage per shoot
2. `seagrass.subplant.R`: Code to analyse and visualise *Phorcus turbinatus* density and distribution data.
    - **Input** = `density.csv`, `distribution.csv` from **Predation**
    - **Output** = Figure 1, *Phorcus turbinatus* site-specific density and distribution results

**Seagrass plant**
1. `distribution.csv`: Vertical distribution data for all studied species.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras
    - **date** = date given as DD.MM.YY
    - **time** = time given as HH:MM
    - **species** = categorical variable with levels *Phorcus turbinatus*, *Stramonita haemastoma*, *Thalassoma pavo*, *Hermodice carunculata* and *Hexaplex trunculus*
    - **position** = position in metres in relation to sea level at the time of observation
    - **tide** = tidal level in metres above lowest astronomical tide at the time of observation provided by the Lampedusa tide station 
    - **lat.position** = position in relation to lowest astronomical tide given in metres
2. `predation.csv`: *In situ* predation data.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras
    - **date** = date given as DD.MM.YY
    - **time** = time given as HH:MM
    - **position** = position in metres in relation to sea level at the time of deployment
    - **tide** = tidal level in metres above lowest astronomical tide at the time of deployment provided by the Lampedusa tide station
    - **lat.position** = position in relation to lowest astronomical tide given in metres
    - **length** = *Phorcus turbinatus* shell length given in millimetres
    - **mass** = *Phorcus turbinatus* mass given in grams
    - **predation** = predation event recorded as 1
    - **notes** = notes on what predator ate *Phorcus turbinatus* or what organism occupied its empty shell
3. `muricids.csv`: *In vitro* predation data.
    - **date** = date given as DD.MM.YY
    - **level** = aquarium water level given in metres
    - **volume** = aquarium water volume given in litres
    - **space** = vertical air space between water level and mesh given in metres
    - **fullness** = aquarium water volume expressed as a proportion of total aquarium volume
    - **predator** = categorical variable with levels *Stramonita haemastoma* and *Hexaplex trunculus*
    - **prey.length** = *Phorcus turbinatus* shell length given in millimetres
    - **prey.mass** = *Phorcus turbinatus* mass given in grams
    - **predator.length** = muricid shell length given in millimetres
    - **predator.mass** = muricid mass given in grams
    - **temperature** = average aquarium water temperature over the 24-hour experiment (°C)
    - **sst** = average sea surface temperature for the day of the experiment (°C) provided by the Maltese weather station
    - **difference** = temperature difference between aquarium and sea surface (°C)
    - **predation** = predation event recorded as 1
4. `predation.R`: Code to analyse and visualise *in situ* predation and distribution data.
    - **Input** = `predation.csv`, `distribution.csv` 
    - **Output** = Figure 2, *in situ* predation and distribution results across sites
5. `muricids.R`: Code to analyse and visualise *in vitro* predation data.
    - **Input** = `muricids.csv` 
    - **Output** = Figure S4, *in vitro* predation results
    
**Seagrass habitat**
1. `length.mass.csv`: Paired measurements of length and mass for all studied gastropods.
    - **species** = categorical variable with levels *Phorcus turbinatus*, *Stramonita haemastoma* and *Hexaplex trunculus*
    - **length** = shell length given in millimetres
    - **mass** = gastropod mass given in grams
2. `physiology.csv`: Immersed, predator-free survival data for *Phorcus turbinatus*.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras
    - **date** = date given as DD.MM.YY
    - **time** = time given as HH:MM
    - **position** = position in metres in relation to sea level at the time of deployment
    - **tide** = tidal level in metres above lowest astronomical tide at the time of deployment provided by the Lampedusa tide station
    - **lat.position** = position in relation to lowest astronomical tide given in metres
    - **length** = *Phorcus turbinatus* shell length given in millimetres
    - **mass** = *Phorcus turbinatus* mass given in grams (NA when shell was empty)
    - **day** = time since deployment (initiation of full immersion) given in days
    - **survival** = survival recorded as 1
2. `Irradiance.R`: Code to analyse and visualise length-mass relationship and survival data.
    - **Input** = `length.mass.csv`, `physiology.csv`
    - **Output** = Figure 3, Figure S2, physiology and allometry results

**Fauna**


Luka Seamus Wright, 28 February 2023
