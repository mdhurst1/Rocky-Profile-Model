# MCMC background and Dakota Implementation

The chosen means of implementing multi-objective optimisation was with Dakota optimisation software (Adams et al., 2019). The ‘black-box’ optimisation software was chosen to work with the model because of the flexibility and ease of testing a variety of methods and available functionality within the software. The Queso Bayesian calibration library (Estacio-Hiroms et al., 2016) was utilised to perform a metropolis-hastings MCMC simulation within Dakota. Details of MCMC algorithm not covered in this documentation.

## Compiling the code for Optimisation 

(check dakota is installed)
`$ dakota -v`

1. Clean up directory 
`$ make clean -f RPM_dakota.make`

2. Compile RPM c++ code 
`$ make -f RPM_dakota.make`

3. Execute dakota command 
`$ dakota -i RPM_dakota_CB.in -o RPM_dakota_CB.out`

## Input data needed for multi-objective optimisation:

Data files found in Rocky-Profile-Model/driver_files/Data/ 

1. RSL history data:
Working example uses RSL history from past 8000 years extracted from a GIA model ([Bradley et al., 2008](https://doi.org/10.1002/jqs.1481)). Example file uses RSL (m), where present day is 0 m, rather than RSL rate.

2. Across-shore transect of measured CRN concentrations:
Compatible with chosen CRN, working example uses 10Be concentrations corrected for chemistry and inheritance background, with 1 sigma errors calculated. Position of samples (m) need to be corrected to distance from cliff base (cliff base = 0 m).

3. Across-shore topographic swath profile:
Extracted from a Digital Surface Model. Swath profile of shoreline-normal transect across location of CRN samples. X positions corrected to distance from cliff base (same as CRN data).


## Dakota input file

(RPM_dakota_CB_1.in)

Dakota requires an input file specifying the functionality required from the software package. A Dakota input file is split into 5 sections: environment, method, variables, interface and responses. The purpose and contents of each Dakota input file subsection will be discussed in separate sections below. Here is an example input file used with RPM_CRN:

```
#  DAKOTA INPUT FILE - RPM_test.in
#  Usage:
#  dakota -i RPM_dakota_w.in -o RPM_dakota_w.out > RPM_dakota.stdout

environment
     tabular_data
     tabular_data_file "RPM_dakota_CB_1.dat"

method,
	 bayes_calibration queso 
     metropolis_hastings
     proposal_covariance
     diagonal values 0.1 0.1 0.1
     scaling
     chain_samples = 10000 seed = 348
	 #burn_in_samples = 100
	 chain_diagnostics confidence_intervals
	 export_chain_points_file  'posterior.dat'
	 output debug

variables,
	 uniform_uncertain = 3
	 lower_bounds  1 -5 -2
     upper_bounds  3 -1 -0.8
	 descriptor "FR" "K" "Y"


interface,
	fork
	  asynchronous
	  evaluation_concurrency 1

	  failure_capture recover 999999 999999
	  
	  parameters_file = "params.in"
	  results_file    = "results.out"
	
	  copy_files "template_dir/*"

	  analysis_driver = "python ./RPM_driver_CB.py"


      work_directory
      named "run"
	  directory_tag
      directory_save
      file_save

responses,
	calibration_terms = 2
    response_descriptors "rmse" "CRNrmse"
    primary_scale_types = "value" "value"
    primary_scales = 0.58 932
    weights = 0.5 0.5
	no_gradients
	no_hessians

```
### Dakota Input Environment
This first keyword block is optional and is used to specify general Dakota settings such as the tabular data file. In the working example we have also specified a name for the output tabular datafile. See Dakota manuals for further options although they are not needed for this example.

### Dakota input method
```
bayes_calibration queso 
metropolis_hastings
```
The method keyword block specifies what iterative method Dakota will perform. In this working example, a metropolis hastings algorithm from the bayes calibration queso library has been specified.

```
samples = 10000 seed = 348 
```
The other compulsory option for the method block is to specify the number of chain samples, this tells Dakota how many Markov Chain Monte Carlo posterior samples to run. For full MCMC simulations we have tested 10,000 samples, this may need to be more for different applications. For debugging or if in testing phases, this number can be set to less but be wary that the full set of final results including best fit results, uncertainties, acceptance rates etc. are not viable until a sufficient number of chain samples has been explored. Search MCMC chain diagnostics for more in-depth information on this subject. 

If seed is specified the use of the same random seed in identical studies will generate identical results, making a stochastic study repeatable. Random numbers are drawn in the metropolis hastings algorithm when making the accept/ reject decision. 

If `burn in samples` is specified (commented out in the example), this number of samples will be discarded from the beginning of the sample chain. If the burn in period is particularly large, this can be a good option, so that the resulting chain and final results are less dependent on the starting point of the chain.

`Proposal covariance` is an optional function that tells Dakota what technique to use to generate the MCMC proposal covariance of the proposal distribution (also known as the jumping distribution). In the working example we have declared specific values to be given to each of the variables using a diagonal matrix format. The values we have provided correspond to the values along the diagonal of the matrix. Think of the matrix input set up like this: 


| |	FR | K | Y |
|--|--|--|--|
| **FR** | 0.5 | 0 | 0 |
| **K** | 0 | 0.5 | 0 |
| **Y** | 0 | 0 | 0.5 |

This technique was chosen as we found it gave the user the greatest control over tailoring the acceptance rates. The covariance of the proposal distribution dictates how far a jump the next sample will be in the accepted chain. As a result, the proposal distribution has a direct impact on acceptance rates of the chain. The proposal distribution values should be varied to achieve acceptance rates of ~23% which ensures optimal chain mixing ([Sherlock & Roberts, 2009](https://doi.org/10.3150/08-BEJ176)). Dakota calculates rejection rates (acceptance rates = 100% - rejection percentage) every 500 runs and this can be seen in the diagnostic outputs within the QuesoDiagnostic directory.

**Acceptance rates too low?** If acceptance rates are significantly less than 23% this often means that the proposal variance is set too high as it is jumping too far away from good, previously acceptance results, so rejection is more likely. Set a smaller proposal covariance to increase acceptance rates.

**Acceptance rates too high?** If acceptance rates are significantly higher than 23%, this often means that the proposal covariance is too small. If the next chain sample is selected very close to the previously accepted sample, the likelihood of rejection is low as the chain is already in a good parameter space. Set a higher proposal covariance to reduce acceptance rates.

Other factors also affect acceptance rates that are not discussed here so background reading on MCMC simulations is advised. 

The proposal covariance also impacts the time taken for the chain to converge. Small proposal covariance will result in slower chain convergence as the number of samples needed to explore the full parameter space will be greater. If the proposal covariance is too high, uncertainties on final best fit variables will be great. 

This may take some time to tailor these proposal covariance values using the metropolis hastings MCMC method. Other MCMC methods in Dakota eg. Delayed rejection, adaptive metropolis and DRAM offer adaptive capabilities to tailor proposal distributions automatically. This is not explored in this study but may be suitable for further studies.

As we are optimising 2 objective functions (topographic profile and 10Be concentrations), both of these response functions need to be scaled to a similar range before combining into a single, composite objective function. The `scaling` key word has to be placed in the methods block as well as further information on the response function scaling in the responses block (see 6.3.5). 

Specify name of output accepted chain file as `export_chain_points_file`. In our working example, our accepted results will be written to the posterior.dat file. To compare accepted sample results to all samples explored in the full chain, posterior.dat file should be compared to the tabular data file RPM_dakota_CB_1.dat in this example.

### Dakota input Variables

`uniform_uncertain`: This first term in the variables block specifies the prior distribution assigned to the chosen free parameters. As we have no prior knowledge and want to make no assumptions on the best fit parameter values, a uniform distribution was selected. Uniform_uncertain is one of many aleatory uncertain variable options within Dakota; with aleatory uncertainty meaning uncertainty that comes from a random process. Once we have a better understanding of the parameter space, further variable uncertainty distributions may be explored, such as normal_uncertain that uses a Gaussian distribution.

`Lower_bounds/ Upper_bounds`: Lower and upper bounds on each free parameter have to be set. These have been carefully selected based on relevant literature, namely [Matsumoto et al. (2016)](https://doi.org/10.1016/j.geomorph.2016.05.017) and [Matsumoto et al. (2018)](https://doi.org/10.1002/esp.4422). 

`Descriptor`: Descriptors are the labels assigned to each of the variables. These need to correspond to the descriptors within the curly brackets used in the input template file

FR: Material resistance
K: Maximum weathering efficacy 
Y: Wave height decay rate   


### Dakota input interface

`Fork`: A fork interface tells Dakota to launch a separate application analysis process: an external script is called (`RPM_driver_CB.py`). `Evaluation_concurrency`: The default setting for an `Asynchronous` execution is to launch all available evaluations simultaneously. Setting the evaluation concurrency to 1 limits the number of concurrent evaluations to 1 at a time. Within a MCMC simulation we want the previous simulation to finish before initiating the next.

`Failure capture`: This function helps Dakota respond to failed Rocky-Profile-Model runs. The model will fail when the modelled topographic shore platform profile does not erode to a length of at least the length of the measured topographic profile of the intertidal shore platform. This results in `inf` values when calculating the residuals between modelled and measured data points. Within the RPM driver file, a fail flag detects failed model runs if `inf` values are calculated. Instead of the rmse scores written to the results file, ‘FAIL’ will be output and received by Dakota. The `response` option tells Dakota to replace the 'FAIL' response received with values specified in the Dakota input file and to continue to the next simulation. Both the topographic rmse and CRN rmse will be assigned values of 999999 so that Dakota understands this to be a bad set of parameter values. This function is powerful as when optimising for topographic profile alone we were able to use this to constrain parameter bounds further to remove ‘fail’ parameter space for the multi-objective optimisation.

`Paramateres_file`: name of file to write parameter values for each simulation run (`params.in`). `Params.in` file is `$1` output from Dakota. 

`Results_file`: name of file to write results for each simulation run (`results.out`). Dakota `Results.out` file is `$2` returned to Dakota. 

`Copy_files`: Files from `template_dir` are copied into each working directory.

`Analysis_driver`: provides the name of the executable script, in this case `RPM_driver_CB.py`. This script is run as an external process, as indicated by the `fork` keyword. 

`Work_directory`: Each simulation will have its own separate working directory, where files are copied from the template directory and `params.in`, `results.out` file for each run can be found. `Directory_tag` will tag each run directory with the simulation number, resulting in a series of directories with corresponding Dakota files within: run.1, run.2, run.3 etc. `directory_save` ensures the working directories are not deleted after the evaluation is completed (the default). It can be useful when debugging errors to see exactly what files were used for a specific simulation. Parameter and results files are deleted after an evalulation completes as default; `file_save` will save these files in corresponding run working directories. 

### Dakota input responses

`Calibration_terms`: keyword used for queso baysian calibration methods. We have two calibration terms or objective functions we are looking to optimise: Topographic profile rmse and the 10Be concentration profile rmse. `Response_descriptors`: names assigned to each calibration terms. We have use RMSE for topographic rmse and CRNrmse for 10Be concentration rmse. 

`Primary_scales`: `Primary_scales_type` keyword indicates the type of scaling assigned to each objective function. Value is a specified scaling value by which each response function will be divided. We need to scale both response values to a similar range. `Primary_scales` assigned to each objective function are currently using the average error of each data set.

`Weights`: Here you assign different weights to each objective function in order to conctruct the ‘pareto -front’ of optimised results. The example shows equal 50-50 weighting given to each objective function. Weight values have to be number between 0 and 1. 

`No_gradients` and `no_hessians` are used in this simulation. 

### Dakota driver script

The analysis driver or user simulation code allows Dakota to talk to the RPM model and vice versa by reading and writing short data files. Dakota will write a parameter file containing current variable values, start the user’s simulation code (`RPM_driver_CB.py`) that will run an iteration of the RPM model and will write response data to a result file. This process is repeated until all declared simulation numbers are complete.

A user simulation code is usually divided into 3 parts: pre-processing, analysis and post-processing. Our example differs from this as the analysis or RMSE calculations are run within the RPM driver code and post-processing is not needed as we have shortcut this by writing RPM model responses directly into a file format that Dakota can read (`results.out` as `$2`). The sole purpose for our user simulation code is to pre-process our input variable files from Dakota to input into the RPM model and to run the RPM model.

In the first stage of the user simulation code, Dakota utilises dprepro template processing tool to extract the current variables from a Dakota-produced parameter file (`params.in` (`$1`)) and combines them with the template input file (`input_template.yml`) to create a new input file: `inputs.yml` that can be read by the RPM model via this user simulation script. The `inputs.yml` file will have the current values that Dakota has assigned to each of the free parameters we have established in the Dakota input file. In this example: WaveAttenuationConst, Resistance and WeatheringRate. 

```
input_template = "input_template_w.yml"
inputs = "inputs.yml"
call(["dprepro", sys.argv[1], input_template, inputs])

```
In the second stage of the user simulation code, an iteration of the RPM model is run by calling a subprocess within the python script. The `RPM_dakota_driver.cpp` code is executed by launching a string of arguments expected by the driver file. The make files should be compiled before starting the Dakota simulation. 

## Input template Dakota file 

`input_template.yml`

The dprepro script will search the `input_template.yml` file for fields marked with curly brackets and then create a new file (`inputs.yml`) by replacing these targets with the corresponding numerical values for the variables. The free parameters chosen to vary in the optimisation simulation as stated in the Dakota input file (`RPM_dakota.in`), need to correspond to the fields marked with the curly brackets in the `input_template.yml` file. 


## Dakota output files 

`RPM_dakota_CB.out` Main output file: contains information of every function evaluation. Copy of input file at the beginning of file and timing and summary statistics at the end. Best fit results for each variable can be found in this file. 

`RPM_dakota_CB.dat` All chain samples tab-delimited text file. Summarises inputs and outputs to the function evaluator (non-scaled or weighted responses).  Names of inputs and outputs will match descriptors in the input file. 

`posterior.dat` Accepted chain samples for the posterior distribution. Non-scaled or weighted responses 

`dakota.rst` Dakota restart file 

`NonDQUESOLogLike.txt` Variable values, scaled and weighted responses and log likelihood for every chain sample

`QuesoDiagnostics/display_sub0.txt` Contains stats for rejection rates every 500 runs

`QuesoDiagnostics/mh_output_sub0.m` Final rejection rate score for full chain 

`QuesoDiagnostics/raw_chain.m` Accepted sample positions (no responses) 

`QuesoDiagnostics/raw_chain_sub0.m` Accepted sample positions (no responses) 

`QuesoDiagnostics/raw_chain_loglikelihood.m` Accepted log likelihood responses 

`QuesoDiagnostics/raw_chain_loglikelihood_sub0.m` Accepted log likelihood responses

`QuesoDiagnostics/raw_chain_logtarget.m` Posterior responses (log-posterior = log-likelihood - log-prior)

`QuesoDiagnostics/raw_chain_logtarget_sub0.m` Posterior responses 


##Plotting Dakota results 

Jupyter notebooks can be found in Rocky-Profile-Model/plotting_functions/Dakota_plotting/



