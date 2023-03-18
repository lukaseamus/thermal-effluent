# Impacts of thermal effluent on *Posidonia oceanica* and associated macrofauna
This repository contains data and annotated R code accompanying article 10.3354/meps14261 in *Marine Ecology Progress Series*. The repository containes the four data and code folders **Introduction**, **Density**, **Predation** and **Physiology** alongside all final **Figures** in vector format. Below is a description of each column within each dataset and the input and output of each R script.

**Introduction**
1. `Steromphala.csv`: Data on the lower and upper limits of the vertical distributions of *Steromphala cineraria* and *S. umbilicalis*.
    - **reference** = author name and year of referenced data
    - **doi** = digital object identifier or uniform resource locator of referenced data
    - **species** = categorical variable with levels *Steromphala cineraria* and *Steromphala umbilicalis*
    - **lower** = lower distributional limit in relation to lowest astronomical tide given in metres
    - **upper** = upper distributional limit in relation to lowest astronomical tide given in metres
2. `Steromphala.R`: Code to calculate introductory descriptive statistics.
    - **Input** = `Steromphala.csv`
    - **Output** = *Steromphala* spp. distribution descriptive statistics in paragraph two of the introduction

**Density**
1. `density.csv`: Density data for *Phorcus turbinatus*, *Stramonita haemastoma* and *Thalassoma pavo*.
    - **site** = categorical variable with levels Xwejni, Dwejra and Ras 
    - **date** = date given as DD.MM.YY
    - **species** = categorical variable with levels *Phorcus turbinatus*, *Stramonita haemastoma* and *Thalassoma pavo*
    - **original** = count of individuals in defined area (40-cm quadrat for *Phorcus turbinatus* and *Stramonita haemastoma*, 1×10-m transect for *Thalassoma pavo*)
    - **adjusted** = density given per square metre
2. `density.R`: Code to analyse and visualise *Phorcus turbinatus* density and distribution data.
    - **Input** = `density.csv`, `distribution.csv` from **Predation**
    - **Output** = Figure 1, *Phorcus turbinatus* site-specific density and distribution results

**Predation**
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
    
**Physiology**
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

Luka Seamus Wright, 28 February 2023
