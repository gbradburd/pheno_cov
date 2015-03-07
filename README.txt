################################################
################################################
#	README
################################################
################################################

PREAMBLE
Caveat user: this script was designed specifically for
use in the application to data in this project.  It 
has not been designed as a method for general distribution, 
and has therefore not been optimized for speed or generalized for
use with any dataset. 
Please contact Gideon Bradburd at gbradburd(at)ucdavis(dot)edu with questions.

Also, check the git for updates!

README
This repository contains 2 other documents:

(1) pheno_cov_mcmc.R:
	Contains the R code for the MCMC algorithm and,
	at the end, the code for initiating an analysis

(2) sample_dataset.Robj:
	An R object to be loaded into working memory. Contains 
	two lists, which are described with their contents below:
		1) "data" - the list of data to be passed to the MCMC algorithm:
				"dependent.variable" - vector containing the phenotypic variable measured
					in each sample that is being modeled as multi-variate normal
				"NA.pops" - vector with a '1' at the position of all populations
					in North America and a '0' elsewhere
				"Euro.pops" - vector with a '1' at the position of all populations
					in Europe and a '0' elsewhere
				"Enviro" - environmental covariate measured at each sample location, 
					e.g., temperature or precipitation
				"genetic.covariance.matrix" - the covariance in allele frequencies
					taken across loci between all samples
				"identity.matrix" - an identity matrix of the same dimensions as the
					genetic.covariance.matrix
				"sd.dependent.variable" - the standard deviation of the 
					dependent.variable, used to standardize the observations.  
					not necessary for running the code, but useful to keep in the same 
					list for un-transforming the results after the model has run
				"sd.enviro" - the standard deviation of the Enviro variable, used to 
					standardize the observations.  not necessary for running the code, 
					but useful to keep in the same list for un-transforming the results 
					after the model has run
		2) "mcmc.operators" - the list of parameters that specify the behavior of the
				MCMC algorithm
				"Beta.NA.stp" - the scale of the tuning parameter (the standard deviation
					of the normal distribution from which updates are drawn on a 
					parameter) on the fixed effect of being in North America
				"Beta.Euro.stp" - the scale of the tuning parameter on the fixed effect of
					 being in Europe
				"Beta.Enviro.stp" - the scale of the tuning parameter on the fixed effect 
					of a unit of the environmental covariate
				"alpha.stp" - the scale of the tuning parameter on the scalar multiplier
					on the genetic covariance (a quantity related to the heritability of 
					the trait)
				"sigma.stp" - the scale of the tuning parameter on the nugget effect on 
					the covariance (in this case, a scalar multiplier on the identity matrix
					that captures random effects)
				"number.generations" - how many generations the MCMC algorithm should run
				"sample.frequency" - the frequency with which parameter values from the 
					MCMC are sampled.  a sample.frequency of 100 means that a value for
					each parameter is sampled and stored every 100 generations
				"save.frequency" - the frequency with which the MCMC output object is
					saved.
				"print.frequency" - the frequency with which updates from the MCMC
					are printed to the console
