# Hospital Closures

Analysis of the effects of critical access hospital designation on hospital closures in the United States. The analysis is organized around four subfolders:

1. **Background**: Contains background information on the critical access hospital program and the data used in the analysis. This file is not yet under git version control.
2. **Data**: Contains the data used in the analysis. This folder is not under git version control.
3. [Data Code](/data_code/): Contains the code used to clean and merge the data used in the analysis. It takes as inputs the data from two seperate repos listed below. It also calls a **functions.R** script that contains functions used in the data cleaning and merging process (namely the haversine distance calculation).
    - [AHA Data](https://github.com/imccart/aha-data), which houses the code to manage the American Hospital Association's Annual Survey data. This repo creates the **aha_data** referenced in the **buid-data.R** code file.
    - [Critical Access Hospital Data](https://github.com/imccart/cah), which houses the code to manage the critical access hospital data. This repo creates the **cah_data** referenced in the **build-data.R** code file. Note that these data are from the Flex Monitoring Team and are separate from the critical access information in the AHA data. These data sources are relatively consistent; however, the flex monitoring program identifies CAH designations much earlier than the AHA data, as they include information on demonstration programs that predate the federal CAH program.
4. [Analysis](/analysis/): Contains the code used to analyze the data. The code is organized into three files:
        - **_run-analysis.R**: This file loads the data from the data folder and cleans and merges the data. It then calls the **sum-stats.R** code and the **dd-estimates.R** code. It takes as input the 'aha_data' and 'cah_data' created by the data code.
