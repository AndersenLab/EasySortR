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

### Data Structure

Once the data are read in, we need to decide whether to keep the data in wide format or switch it to [long format](http://en.wikipedia.org/wiki/Wide_and_narrow_data). My original code from ExperimentRunner/SimpleDataProcess does all of the regression columnwise with the data in wide format. This is a rough and rather Rube Goldberg-ian method method that is very prone to failure. My vote is to switch regression from wide to long format, which may mean switching "normal" data frames to long format during analyses. Is this ok with everyone?

![My regression code from last year works something like this](http://upload.wikimedia.org/wikipedia/commons/a/a6/Professor_Lucifer_Butts.gif)



