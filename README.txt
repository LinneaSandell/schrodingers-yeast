This directory contains:

Raw data
- bioscreens: Optical density measurements for four different growth rate experiments.
Codes for what line was run in what well.
- FACS_august: flow cytometry files and codes for FACS done for DS1
- FACS_november: flow cytometry files and codes for FACS done for DS2


Scripts
- bioscreen_functions.R Makes spline fits to growth rate data and extracts the maximum growth rate.
- 01_cleandata.R     Takes optical density (OD) measurements from raw-data directory, uses codes.txt to add information about the sample that was run, and extracts summary statistics for each growth curve: initial OD, maximum growth rate, OD at the maximum growth rate, final OD, final time (how long the growth assay was done). Prints out a summary file for each growth rate experiment in data/bioscreens.
- 02_DayBatchEffects.R      Calculates variance attributable to the random effects of Day, machine run, plate. It also calculates correlations between growth rate measurements of identical lines taken out of the freezer and run in the different growth rate experiments.
- 03_Statistics.R is the script where I do the statistical analyses reported on in the manuscript and produces the figures.
- 04_FACS.R takes files from the flow cytometer (.fcs) and print figures into the figures directory to interpret DNA content.