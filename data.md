# Data

This page describes our data sources and the code used to clean and merge the data.

## Overview
Our data are organized into input and output subfolders. The input folder contains our underlying data sources, and the output folder contains the intermediate data files ultimately called in the analysis. Details of our input sources and the code to create our data outputs are described below.

## Data Sources
We take as inputs the data from two seperate repos listed below. 

- [AHA Data](https://github.com/imccart/aha-data), which houses the code to manage the American Hospital Association's Annual Survey data. This repo creates the **aha_data** referenced in the **buid-data.R** code file.
    - We also collected data on hospital mergers, closures, etc. from the AHA summary of changes files. These data were manually collected from our RA at Emory, Jonathan Kim. The manual data are available [here](https://docs.google.com/spreadsheets/d/19CyPL0wYFJJfPxTjr89kvKIiGIqHZGC69AXm4Iz6rrY/edit#gid=1195668795), restricted to Emory users.
    - The OIG also publishes reports on historic hospital closures. We have collected these data but have not digitized them for inclusion in the analysis.

- [Critical Access Hospital Data](https://github.com/imccart/cah), which houses the code to manage the critical access hospital data. This repo creates the **cah_data** referenced in the **build-data.R** code file. Note that these data are from the Flex Monitoring Team and are separate from the critical access information in the AHA data. These data sources are relatively consistent; however, the flex monitoring program identifies CAH designations much earlier than the AHA data, as they include information on demonstration programs that predate the federal CAH program.


## Code
The [data-code](/data-code/) folder contains the code used to clean and merge the data used in the analysis. It also calls a **functions.R** script that contains functions used in the data cleaning and merging process (namely the haversine distance calculation), and an **api-keys.R** script that contains the API keys used to access some geographic crosswalk data. The API key is not under version control but can be easily requested from the Census Bureau [here](https://api.census.gov/data/key_signup.html).
