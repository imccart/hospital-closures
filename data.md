# Data

This page describes our data sources and the code used to clean and merge the data.

## Overview
Our data are organized into input and output subfolders. The input folder contains our underlying data sources, and the output folder contains the intermediate data files ultimately called in the analysis. Details of our input sources and the code to create our data outputs are described below.

## Data Sources
We take as inputs the data from three seperate repos listed below. 

- [AHA Data](https://github.com/imccart/aha-data), which houses the code to manage the American Hospital Association's Annual Survey data. This repo creates the **aha_data** referenced in the **buid-data.R** code file.
    - We also collected data on hospital mergers, closures, etc. from the AHA summary of changes files. These data were manually collected from our RA at Emory, Jonathan Kim. The manual data are available [here](https://docs.google.com/spreadsheets/d/19CyPL0wYFJJfPxTjr89kvKIiGIqHZGC69AXm4Iz6rrY/edit#gid=1195668795), restricted to Emory users.
    - The OIG also publishes reports on historic hospital closures. We have collected these data but have not digitized them for inclusion in the analysis.

- [Critical Access Hospital Data](https://github.com/imccart/cah), which houses the code to manage the critical access hospital data. This repo creates the **cah_data** referenced in the **build-data.R** code file. Note that these data are from the Flex Monitoring Team and are separate from the critical access information in the AHA data. These data sources are relatively consistent; however, the flex monitoring program identifies CAH designations much earlier than the AHA data, as they include information on demonstration programs that predate the federal CAH program.

- [HCRIS Data](https://github.com/imccart/HCRIS), which houses the code to manage the Healthcare Cost Report Information System data. This repo creates the **hcris_data** referenced in the **build-data.R** code file. Some of the specific variables we pull from the cost reports are listed below, along with their location in the v1996 cost report data.
  - *Total Charges*: The total amount the hospital billed for services provided to patients (**not** the amount the hospital was paid for those services). Located in worksheet G-3, line 1.
  - *Total Discounts*: The sum of all discounts given to patients, including charity care, bad debt, and contractual allowances. Located in worksheet G-3, line 2.
  - *Net Patient Revenue*: The total amount the hospital was paid for services provided to patients, minus the total discounts given to patients. Located in worksheet G-3, line 3.
  - *Total Operating Expenses*: The total amount the hospital spent on operating expenses. Located in worksheet G-3, line 4.
  - *Medicare Discharges*: The number of discharges for patients with Medicare insurance. Located in worksheet S-3, line 1.
  - *Medicaid Discharges*: The number of discharges for patients with Medicaid insurance. Located in worksheet S-3, line 1.
  - *Uncompensated Care*: The sum of charity care and bad debt. Located in worksheet S-10, line 30.
  - *Cost to Charge*: The ratio of the hospital's total operating expenses to its total charges. Located in worksheet S-10, line 24.
  - *New Capital Assets*: The amount the hospital spent on new capital assets. Located in worksheet A-7, line 9.
  - *Cash on Hand*: The amount of cash the hospital had on hand at the end of the fiscal year. Located in worksheet G, line 1.

- [Form 990 data](https://github.com/imccart/form-990s), which houses the code to manage the Form 990 data. This repo creates the **form990_ahaid** data referenced in the **build-data.R** code file. The Form 990 data are used to collect information on hospital financial outcomes, such as total revenue, total expenses, and net income, particularly for earier years when HCRIS data are not available.


## Code
The [data-code](/data-code/) folder contains the code used to clean and merge the data used in the analysis. It also calls a **functions.R** script that contains functions used in the data cleaning and merging process (namely the haversine distance calculation), and an **api-keys.R** script that contains the API keys used to access some geographic crosswalk data. The API key is not under version control but can be easily requested from the Census Bureau [here](https://api.census.gov/data/key_signup.html).

The **build-data.R** script is the main script that calls the other scripts in the data-code folder to clean and merge the data. The script outputs three datasets for futher analysis: 1) hospital data formed from the combination of AHA surveys, HCRIS data, CAH data, and Form 990 data; 2) a dataset of hospital mergers and closures; and 3) a dataset summarizing the distance between a given hospital and its nearest neighbor (useful in assessing eligibility for CAH designation).
