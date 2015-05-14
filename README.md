# Sorter Data Processing Package Problems

## Reading in data

### Identification variables

We need to define a single set of universal identifying variables and a way to feed them into a function. Currently (RIAILs data) this relies on a specific naming structure for the data files. Currently the file naming and directory structure provide the following information:

+ Date
+ Experiment Type (RIAILs, GWAS, NILs, etc.)
+ Round (numerical value such as 1 or 2 for RIAILs)
+ Assay (an alphabetic identifier that serves as a proxy for the date)
+ Plate Number
+ Drug/Condition
+ Concentration (in some cases only)

This information is all included in the final datasets along with the COPAS-provided well row and column values. We need to decide which set of these identifiers are applicable universally and should thusly be included in the data reading functions. As of right now, in order to process the RIAILs data sets, Assay is needed to regress out interday effects, plate number is needed to math sort and score plates to one another, and Drug and Concentration are needed to summarize some replicated plates and to separate different drug concentrations done on the same day.

### Setup Plates

Do all scored plates have corresponding setup plates in all experiments? If not, how should we handle getting the number scored into read-in data frames?

### Data Structure

Once the data are read in, we need to decide whether to keep the data in wide format or switch it to [long format](http://en.wikipedia.org/wiki/Wide_and_narrow_data). My original code from ExperimentRunner/SimpleDataProcess does all of the regression columnwise with the data in wide format. This is a rough and rather Rube Goldberg-ian method method that is very prone to failure. My vote is to switch regression from wide to long format, which may mean switching "normal" data frames to long format during analyses. Is this ok with everyone?

![My regression code from last year works something like this](http://upload.wikimedia.org/wikipedia/commons/a/a6/Professor_Lucifer_Butts.gif)

### Plate-by-Plate vs Whole Assay Analysis

Originally (circa Summer 2014) all analysis was done by joining all plates from a single assay/experiment together in one giant data frame, creating a data frame each for all setup data, score data, and control data. These data frames could then be filtered to just the data needed to the plates being analyzed. Is this the way that we want to keep things (plate reading function will take a directory) or would it be preferable to read plates in individually and allow the user to concatenate data frames or loop as necessary? Perhaps a function for each?





