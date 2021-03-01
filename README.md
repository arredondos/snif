# SNIF: Structured Non-stationary Inference Framework



SNIF infers the optimal pair of demographic values from a given dataset, with the aim of obtaining the best IICR that fits the data, according to its set parameters.


SNIF is currently only available in Linux systems and must be called through the command line. It requires python 3.6.9 (minimum) and LaTex (for plotting)
It's currently composed of 4 different python scripts: snif.py; snif2tex.py; collect.py and clean.py 
These can be run independently, or consecutively as the user prefers.
Two first two scripts relate to the inference and plotting of the data, and the last two relate to the management of the output files.

-------------------------------
Guide Structure

```
1.INFERRING AND PLOTTING THE DATA
1.1.INFERRING THE DATA
1.1.1.USING SNIF.py
1.1.1.1.PARAMETERS FILE DESCRIPTION
1.1.1.2.INFERRING PARAMETERS
1.1.1.3.OPTIMIZATION PARAMETERS
1.1.1.4.SETTINGS
1.1.1.5.LOOPING PARAMETERS
1.1.2.OUTPUT FILE FORMAT
1.2.PLOTTING THE DATA
1.2.1.USING SNIF2LATEX.py
1.2.1.1.PARAMETERS FILE DESCRIPTION
1.2.1.2.INFERRING AND PLOTTING 
2.2.OUTPUT FILES MANAGEMENT
2.2.1.USING COLLECT.py
2.2.2.USING CLEAN.py
3.OTHER ISSUES
```

-------------------------------
1.INFERRING AND PLOTTING THE DATA


1.1.INFERRING THE DATA

SNIF currently accepts 3 different data input formats:
PSMC files: from genomic data or simulated sequences
Ms command output files: Make sure the ms command file starts with: "./ms command_input"
Json file: a descriptive file with all the model parameters to be tested.

SNIF tries to infer the optimal n-island model, based on the input data and the given parameters. It is thus very important to define the parameters of the n-island model that we want to test.
Most of the parameters are defined in terms of lower and upper bounds that restrict the search algorithm n-island model. There is no set of parameters that fit all types of data. It is important to consider whether the set parameters fit the type of data we are analysing.


1.1.1.USING SNIF.py
 
Snif.py contains the base script with the computational base. This script invokes and runs all the necessary parameters for an inference and thus should NOT be modified.

The parameters are distributed into three different groups: 
Inference parameters: relates to the inference of the data
Optimisation parameters: relates to the search algorithm
Settings parameters:relates to the internal SNIF library and output files 

It is recommended to create a parameters file that can be freely modified without affecting the base script.
All of the necessary parameters can be obtained in the Snif.py script. 
Feel free to use the parameters_SNIF file as a reference for your data. Unused parameters need to be present however the values can be substituted with "xxx". Here is a description of all required SNIF parameters:


-----
1.1.1.1.PARAMETERS FILE DESCRIPTION

1.1.1.2.INFERRING PARAMETERS
#Set the pathway to the input file, this file must be one in one of the three accepted input files:
data_source


#Specify the type of input file: PSMC, MSMC, JSON, MSCommand. The MSMC file format is not fully implemented yet.
source_type

#Specify the type of IICR: Exact (such as PSMC from genomic data), T_sim (T2 simulations), Seq_sim (Simulated sequences). This is important to distinguish simulated from genomic data.
IICR_type 


#Refers to the population effective size. Required in a MSCommand input file for scaling purposes. JSON files should contain this information in the n_ref option
ms_reference_size

#This correspondes to the number of T2 samples in a T_sim simulation. This value is only used in a JSON source file. MSCommand source files already contain this information
ms_simulations


#Applied mutation rate of the psmc, make sure to check this value before each run!!
psmc_mutation_rate


#Should the source file be a simulated sequence, this corresponds to the number of chromosomes to be simulated
psmc_number_of_sequences


#The length of simulated sequences to be simulated
psmc_length_of_sequences = int(2e6)


#Whether to scale the results. If yes they will be presented in generations times, otherwise T2 times.
infer_scale = True, 

#Specifies the time period where no inferences will occur. May need to be adjusted, always check the obtained plot.
#If scaled the time scale will appear in generation times
#This parameter requires the user discretion in whether very recent or ancient events should be considered as noise or as important information
data_cutoff_bounds,

#Time intervals, corresponds to the p parameters in the PSMC
data_time_intervals = 64,


#The distance function (which returns the fitting) can be calculated on two different ways:
#ApproximatePDF: the preferred approach, uses a weighted method. Time periods with more coalescence data will be given more weight (skews towards more recent time)
#Visual: similar to a manual approach, it calculates the best fit considering all time periods have the same weight. Will infer both ancient and recent events.
distance_function

#regulates the aproximatePDF method. 1 is the default. Omega<1 will give less weight to the recent past, allowing inferences in the most ancient past
#Omega> 1 will give more weight to the recent past, biasing towards most recent events.
distance_parameter


#Distance/Fitting error. Should the minimum amount of iterations be complete, the search algorithm will stop when the obtained fitting error is small than the maximum allowed.
#If this value is not met, the program will run as many iterations as allowed (upper bound of the number of rounds_per_test_bounds).
distance_max_allowed


#Time period where the distance/fitting will be computed.
#This time period needs to be a subset of the cutoffbounds.
#The distance_computation_interval affects the calculation of the distance. Do not confuse with the cutoff bounds which specify time periods where no inferences shall be made.
distance_computation_interval

#Number of allowed iterations. It will ALWAYS run the minimal set of iterations regardless of the distance_max_allowed parameter.
rounds_per_test_bounds

#Number of repetitions a given file will be inferred. This will rerun the file with the same parameters described.
repetitions_per_test
							 
#A component or time-period, designates a period in time in which no demographic changes occur. 
#For n components, n-1 events/changes in migration/size will occur. 
#This value is data dependent, and requires user discretion. The user should analyse what is the best minimal number of components that can best explain its data
number_of_components


#Number of islands to be tested. n=2 is the minimum allowed
#The number of islands inferred will be contained inside this interval. 
#Should the number of inferred islands be equal to the maximum number allowed then the inference should be rerun with a wider interval.
bounds_islands

#Migration bounds to be tested.
#All inferred migration rates will be contained inside this interval. 
bounds_migrations_rates

#What should be the size of the demes?
#The size of the demes will also depend on the number of inferred islands.
bounds_deme_sizes

#If necessary it is possible to define a time periods where the demographic events should be inferred.
#This time periods is a subset of the distance_computation_interval
#Allows to ignore ancient or recent events that should not be inferred
bounds_event_times

#This controls the inference of the effective size. 
bounds_effective_size = (100, 100000)


##Advanced##
It is possible to specify bounds for each given component for the following parameters: migration_rates/ deme_size and event_time
The given bounds (min, max) should be substituted with the following formula, Cn as the maximum number of components to be tested:
 bounds_migrations_rates = [(bounds C1),(bounds C2),...,(bounds Cn)] 
 bounds_deme_sizes = [(bounds C1),(bounds C2),...,(bounds Cn)]
 bounds_event_times = [[(bounds C1),(bounds C2),...,(bounds Cn)] #starts at T1
All three parameters need to be given in the same format
 The bounds parameters can be also given as: [(bounds)]*Cn

-----

1.1.1.3.OPTIMIZATION PARAMETERS
These parameters affect the search algorithm of SNIF.
For a detailed description of how each of these parameters affect the algorithm please #contact Armando Arredondo 

#affects the search algorithm
optimization_parameters
strategy
maxiter
popsize
tol
mutation 
recombination
workers


-----

1.1.1.4.SETTINGS
The settings parameters relate to the internal aspects of the SNIF program


#Where is the location of the library files (libsnif.so) file
#The library file 
#This file is MANDATORY for the scrip to run.

static_library_location

#Allows a custom name tag to be applied to all output files.
custom_filename_tag

#The directory/folder where the output files should go
output_directory

#Default output directory. If no folder is specified in the previous line
default_output_dirname 

-----
-----
1.1.2.LOOPING PARAMETERS

The descriptive format of the parameters file allows for multiple testing of different parameter values should it be necessary.
This can be achieved by including a for loop at the beginning of the script, passing the variable to the parameter we want to loop through.

Reminder: python language is indentation dependent. Should a for loop be added, all the parameters including the call to infer should be indented.

-----
-----
1.1.3.OUTPUT FILE FORMAT
OUTPUT FILE FORMAT:


The output files are named according to the following convention:

 
1      2      3     4  5    67  8
190919-191919_c03YY_dA_w100_SQ[_TAG]
                 NN  V      JT
		     E      MS
			    X

1: Date (YYMMDD)
2: Time (hhmmss)
3: Inferred components, infer scale ([Y]es, [N]o), infer sizes ([Y]es, [N]o)
4: Distance type ([A]pproximate PDF, [V]isual distance, [E]xact PDF)
5: Omega
6: Source type ([S]MC, [J]SON, [M]s command, [X] Simulation), 
7: IICR type ([Q]-matrix (exact), [T]-sim, [S]eq-sim), 
8: Custom tag (optional)



In case of simulated data, the basename will be according to the following convention:

```
1      2      3      4 5    6     7  8    910 11
190919-191919_sc04YY_D_s500_c03YY_dA_w100_SQ[_TAG]
                  NN C         NN  V      JT
				   E      MS
					  X
```

1: Date (YYMMDD)
2: Time (hhmmss)
3: Simulated components, simulate scale ([Y]es, [N]o), simulate sizes ([Y]es, [N]o)
4: Sampling method ([D]iscrete, [C]ontinuous),
5: Sample size
6: Inferred components, infer scale ([Y]es, [N]o), infer sizes ([Y]es, [N]o)
7: Distance type ([A]pproximate PDF, [V]isual distance, [E]xact PDF)
8: Omega
9: Source type ([S]MC, [J]SON, [M]s command, [X] Simulation), 
10: IICR type ([Q]-matrix (exact), [T]-sim, [S]eq-sim), 
11: Custom tag (optional)


Each successful run of SNIF results in 2 csv files and 1 json file
Should SNIF fail to successfully complete its inference, for any reason, the json file will not be created.
The json file has two functions: 
1)a complete report of all tested parameters for a given run
2)can be used as a SNIF input file.



-----
Given a parameters file named parameters_SNIF.py the script can be called with the following command:

python3 parameters_SNIF.py 



--------------------

1.2.PLOTTING THE DATA

At the moment there are 4 distinct plots.
The IICR plot:
The IICR plot is a representation of the input data (source) and of the obtained model (inferred). 
The better the model the closer the inferred is to the source. Ideally both of them should overlap.

Islands.pdf:
This plot represent all the inferred number of islands across a single run.
The y axis corresponds to the allowed bounds in the inference parameter.

Reference-size (or ref-size.pdf):
This plot represent all the inferred population size across a single run.
The y axis corresponds to the allowed bounds in the inference parameter.

Connectivity Graph:
This plot represent the changes in migration rate across time.
The marginal histograms represent the frequency of inference for a given value of migration or time.
This is useful when comparing multiple runs.

----

1.2.1.USING SNIF2LATEX.py

Snif2tex.py: 
Uses LaTex to create basic graphic reports from the results.
This script can be called independently after SNIF finished running or
Given a parameters file, it is possible to add a call Snif2tex after the call on Snif.py

It is necessary to always check the obtained IICR for your inferences. Especially as you may need to change cutoff bounds values.

It is recommended to create a parameters file containing all the plotting parameters. 
This parameters file allows for a better control of the plotting options. 

Types of plots:
excluded: does not generate a plot
minimum: shows minimalistic version of plot
full: show complete plot

----
1.2.1.1.PARAMETERS FILE DESCRIPTION
The required plot parameters are:

#name and path to the output file you want to analyse.
#You should give the csv file name with no extensions
#ex.: 200219-152354_c03YN_dA_w100_SS
SNIF_basename 

#Width of the plot in cm
plot_width_cm  

#Height of the plot in cm
plot_height_cm

#Type of IICR plot: full, minimum or excluded
IICR_plots_style

#Not implemented
PDF_plots_style

#Not implemented
CDF_plots_style

#Create a plot with with the number of obtained inferred islands. The y axis corresponds to the allowed bounds in the inference parameter
islands_plot_style

#Create a plot with with the inferred effective population size. The y axis corresponds to the allowed bounds in the inference parameter
Nref_plot_style 

#The number of tests to be analysed
test_numbers

#Create one plot file for each test
one_file_per_test

#Plot simulations 
versus_plot_style

#Plot Connectivity graph. X-axis is a time scale, y-axis present the inferred migration values.
#Full version contains marginal histograms, they represent the frequency of the inferred values
CG_style

#Plot the changes in population size
CG_size_history

#Source legend for the connectivity graphs
CG_source_legend

#Number of bins for the histogram of the connectivity plots
histograms_bins 

#What is the scalling unit. If specified in years, the generation_time parameter will be used to calculate it.
scaling_units = TimeScale.Generations

#Generation time to be applied. Always check this value before plotting!
generation_time



Thus, given an external script with all the defined parameters named plotting_parameters this script can be called with:
python3 plotting_parameters.py

Should you require help this script recognises the help [-h] command. 
Call with: python3 snif2tex. -h


--------------------
1.2.1.2.INFERRING AND PLOTTING 

It is possible to merge the SNIF and SNIF2TEX parameters file into a single parameter script file.
Both SNIF and LaTex should be working and installed in the same environment.
Since the plotting is an essential part of the evaluation of the inference, and should always be present,
This removes the need to manually call the plotting script for all output files.

-------------------------------

2.2.OUTPUT FILES MANAGEMENT

As previously mentioned each successful run of SNIF results in 3 output files, otherwise it creates only 2.
This creates a fast way to separate successful from unsuccessful runs. 
The Collect and Clean scripts, declutter the output results folder. These two scripts depend on the file names. So be careful when renaming or modifying file names.

--------------------
2.2.1.USING COLLECT.py

Collect.py: 
Groups into a single folder all 3 output files from a single run.
Should the json file with the corresponding name pattern be missing the csv files will not be moved into a new folder. 
This separates successful from incomplete runs.



Should you require help this script recognises the help [-h] command. 
Call with: python3 collect.py -h

--------------------
2.2.2.USING CLEAN.py

Clean.py: 
Removes all files in the results folder that were not grouped. 
This removes all runs that for some reason do not have all 3 necessary files, created during a successful SNIF run.


Should you require help this script recognises the help [-h] command. 
Call with: python3 clean.py -h

3.OTHER ISSUES



