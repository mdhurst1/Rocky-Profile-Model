# MCMC background and Dakota Implementation

The chosen means of implementing multi-objective optimisation was with Dakota optimisation software (Adams et al., 2019). The ‘black-box’ optimisation software was chosen to work with the model because of the flexibility and ease of testing a variety of methods and available functionality within the software. The Queso Bayesian calibration library (Estacio-Hiroms et al., 2016) was utilised to perform a metropolis-hastings MCMC simulation within Dakota. Details of MCMC algorithm not covered in this documentation.

## Compiling the code for Optimisation

Stuff about specific Dakota driver and make files

## Dakota input file

Dakota requires an input file specifying the functionality required from the software package. A Dakota input file is split into 5 sections: environment, method, variables, interface and responses. The purpose and contents of each Dakota input file subsection will be discussed in separate sections below. Here is an example input file used with RPM_CRN:

```
#  DAKOTA INPUT FILE - RPM_dakota.in
#  Usage:
#  dakota -i RPM_dakota.in -o RPM_dakota.out > RPM_dakota.stdout

environment
     tabular_data
     tabular_data_file "RPM_dakota.dat"

method,
	 bayes_calibration queso 
     dram                           #metropolis_hastings/DRAM
     samples = 1000 seed = 348
     #burn_in_samples = 100
	 output debug

variables,
	 uniform_uncertain = 2 
	 initial_point 0.02 0.000002
	 lower_bounds  0.00005 0.
     upper_bounds  1. 0.001
	 descriptor "FR" "K"

interface,
	fork
	  asynchronous
	  evaluation_concurrency 1

	  failure_capture recover 999999999999999
	  
	  parameters_file = "params.in"
	  results_file    = "results.out"
	
	  copy_files "template_dir/*"

	  analysis_driver = "python ./RPM_driver.py"


      work_directory
      named "run"
	  directory_tag
      directory_save
      file_save

	  
responses,
	num_response_functions = 1
    response_descriptors "rmse"  #/likelihood
	no_gradients
	no_hessians

```
### Dakota Input Environment
This first keyword block is optional and is used to specify general Dakota settings such as the tabular data file. In the working example we have also specified a name for the output tabular datafile. See Dakota manuals for further options although they are not needed for this example.

### Dakota input method
```
bayes_calibration queso 
dram
```
The method keyword block specifies what iterative method Dakota will perform. In this working example, a Delayed Rejection Adaptive Metropolis (DRAM) algorithm from the bayes calibration queso library has been specified [Haario et al., 2006](https://doi.org/10.1007/s11222-006-9438-0). 

```
samples = 1000 seed = 348 
```
The other compulsory option for the method block is to specify the number of chain samples, this tells Dakota how many Markov Chain Monte Carlo posterior samples to run. For full MCMC simulations we have tested 10,000 samples, this may need to be more for different applications. For debugging or if in testing phases, this number can be set to less but be wary that the full set of final results including best fit results, uncertainties, acceptance rates etc. are not viable until a sufficient number of chain samples has been explored. Search MCMC chain diagnostics for more in-depth information on this subject. 

If seed is specified the use of the same random seed in identical studies will generate identical results, making a stochastic study repeatable. Random numbers are drawn in the metropolis hastings algorithm when making the accept/ reject decision. 

If `burn in samples` is specified (commented out in the example), this number of samples will be discarded from the beginning of the sample chain. If the burn in period is particularly large, this can be a good option, so that the resulting chain and final results are less dependent on the starting point of the chain.

`Proposal covariance` is an optional function that tells Dakota what technique to use to generate the MCMC proposal covariance of the proposal distribution (also known as the jumping distribution). In the working example we have declared specific values to be given to each of the variables using a diagonal matrix format. The values we have provided correspond to the values along the diagonal of the matrix. Think of the matrix input set up like this: 


| |	FR | K | Y |
|--|--|--|--|
| FR | 0.5 | 0 | 0 |
| K | 0 | 0.5 | 0 |
| Y | 0 | 0 | 0.5 |

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


###
###
###


## Dakota driver script

Info about RPM_driver.py

etc.