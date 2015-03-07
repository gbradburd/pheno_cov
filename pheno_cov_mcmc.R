##########################################################################################
##########################################################################################
###	Asclepias Data MCMC
##########################################################################################
##########################################################################################

##########################################################################################
#	Required Packages
##########################################################################################

	library(compiler)
	library(mnormt)

##########################################################################################
#	Function Library
##########################################################################################

	Likelihood <- function(dependent.variable,vector.means,covariance){
		LnL <- dmnorm(dependent.variable,vector.means,covariance,log=TRUE)
		return(LnL)
	}

	prior.Beta.NA <- function(Beta.NA){
		dunif(Beta.NA,min=-100,max=100,log=TRUE)
	}

	prior.Beta.Euro <- function(Beta.Euro){
		dunif(Beta.Euro,min=-100,max=100,log=TRUE)
	}

	prior.Beta.Enviro <- function(Beta.Enviro){
		dunif(Beta.Enviro,min=-100,max=100,log=TRUE)
	}	
	
	prior.alpha <- function(alpha){
		dexp(alpha,log=TRUE)
	}

	prior.sigma <- function(sigma){
		dexp(sigma,log=TRUE)
	} 

	Update.Beta.NA <- function(data,last.params,mcmc.operators){
		Beta.NA.prime <- last.params$Beta.NA + rnorm(1,0,mcmc.operators$Beta.NA.stp)
		prior.prob.Beta.NA.prime <- prior.Beta.NA(Beta.NA.prime)
		new.params <- last.params
			if(prior.prob.Beta.NA.prime != -Inf){
				vector.means.prime <- c(Beta.NA.prime*data$NA.pops+last.params$Beta.Euro*data$Euro.pops+last.params$Beta.Enviro*data$Enviro)
				Posterior.probability.prime <- Likelihood(data$dependent.variable,
														vector.means.prime,
														(last.params$alpha*data$genetic.covariance.matrix + 
														last.params$sigma*data$identity.matrix)) + 
															prior.prob.Beta.NA.prime + 
															prior.Beta.Euro(last.params$Beta.Euro) + 
															prior.Beta.Enviro(last.params$Beta.Enviro) + 
															prior.alpha(last.params$alpha) + 
															prior.sigma(last.params$sigma)
				if(exp( Posterior.probability.prime - last.params$Posterior.probability) >= runif(1)){
						new.params$Beta.NA <- Beta.NA.prime
						new.params$Posterior.probability <- Posterior.probability.prime
						new.params$vector.means <- vector.means.prime
						new.params$Beta.NA.accepted.moves <- last.params$Beta.NA.accepted.moves + 1
				}
			}
			new.params$Beta.NA.proposed.moves <- last.params$Beta.NA.proposed.moves + 1
			return(new.params)
	}

	Update.Beta.Euro <- function(data,last.params,mcmc.operators){
		Beta.Euro.prime <- last.params$Beta.Euro + rnorm(1,0,mcmc.operators$Beta.Euro.stp)
		prior.prob.Beta.Euro.prime <- prior.Beta.Euro(Beta.Euro.prime)
		new.params <- last.params
			if(prior.prob.Beta.Euro.prime != -Inf){
				vector.means.prime <- c(last.params$Beta.NA*data$NA.pops+Beta.Euro.prime*data$Euro.pops+last.params$Beta.Enviro*data$Enviro)
				Posterior.probability.prime <- Likelihood(data$dependent.variable,
														vector.means.prime,
														(last.params$alpha*data$genetic.covariance.matrix + 
														last.params$sigma*data$identity.matrix)) + 
															prior.Beta.NA(last.params$Beta.NA) +
															prior.prob.Beta.Euro.prime + 
															prior.Beta.Enviro(last.params$Beta.Enviro) + 
															prior.alpha(last.params$alpha) + 
															prior.sigma(last.params$sigma)
				if(exp( Posterior.probability.prime - last.params$Posterior.probability) >= runif(1)){
						new.params$Beta.Euro <- Beta.Euro.prime
						new.params$Posterior.probability <- Posterior.probability.prime
						new.params$vector.means <- vector.means.prime
						new.params$Beta.Euro.accepted.moves <- last.params$Beta.Euro.accepted.moves + 1
				}
			}
			new.params$Beta.Euro.proposed.moves <- last.params$Beta.Euro.proposed.moves + 1
			return(new.params)
	}

	Update.Beta.Enviro <- function(data,last.params,mcmc.operators){
		Beta.Enviro.prime <- last.params$Beta.Enviro + rnorm(1,0,mcmc.operators$Beta.Enviro.stp)
		prior.prob.Beta.Enviro.prime <- prior.Beta.Enviro(Beta.Enviro.prime)
		new.params <- last.params
			if(prior.prob.Beta.Enviro.prime != -Inf){
				vector.means.prime <- c(last.params$Beta.NA*data$NA.pops+last.params$Beta.Euro*data$Euro.pops+Beta.Enviro.prime*data$Enviro)
				Posterior.probability.prime <- Likelihood(data$dependent.variable,
														vector.means.prime,
														(last.params$alpha*data$genetic.covariance.matrix + 
														last.params$sigma*data$identity.matrix)) + 
															prior.Beta.NA(last.params$Beta.NA) +
															prior.Beta.Euro(last.params$Beta.Euro) + 
															prior.prob.Beta.Enviro.prime + 
															prior.alpha(last.params$alpha) + 
															prior.sigma(last.params$sigma)
				if(exp( Posterior.probability.prime - last.params$Posterior.probability) >= runif(1)){
						new.params$Beta.Enviro <- Beta.Enviro.prime
						new.params$Posterior.probability <- Posterior.probability.prime
						new.params$vector.means <- vector.means.prime
						new.params$Beta.Enviro.accepted.moves <- last.params$Beta.Enviro.accepted.moves + 1
				}
			}
			new.params$Beta.Enviro.proposed.moves <- last.params$Beta.Enviro.proposed.moves + 1
			return(new.params)
	}

	Update.alpha <- function(data,last.params,mcmc.operators){
		alpha.prime <- last.params$alpha + rnorm(1,0,mcmc.operators$alpha.stp)
		prior.prob.alpha.prime <- prior.alpha(alpha.prime)
		new.params <- last.params		
			if(prior.prob.alpha.prime != -Inf){
				Posterior.probability.prime <- Likelihood(data$dependent.variable,
														last.params$vector.means,
														(alpha.prime*data$genetic.covariance.matrix + 
														last.params$sigma*data$identity.matrix)) + 
															prior.Beta.NA(last.params$Beta.NA) +
															prior.Beta.Euro(last.params$Beta.Euro) +
															prior.Beta.Enviro(last.params$Beta.Enviro) + 
															prior.prob.alpha.prime + 
															prior.sigma(last.params$sigma)
				if(exp(Posterior.probability.prime- last.params$Posterior.probability) >= runif(1)){
					new.params$alpha <- alpha.prime
					new.params$Posterior.probability <- Posterior.probability.prime
					new.params$alpha.accepted.moves <- last.params$alpha.accepted.moves + 1
				}											
			}
		new.params$alpha.proposed.moves <- last.params$alpha.proposed.moves + 1
		return(new.params)
	}

	Update.sigma <- function(data,last.params,mcmc.operators){
		sigma.prime <- last.params$sigma + rnorm(1,0,mcmc.operators$sigma.stp)
		prior.prob.sigma.prime <- prior.sigma(sigma.prime)
		new.params <- last.params		
			if(prior.prob.sigma.prime != -Inf){
				Posterior.probability.prime <- Likelihood(data$dependent.variable,
														last.params$vector.means,
														(last.params$alpha*data$genetic.covariance.matrix + 
														sigma.prime*data$identity.matrix)) + 
															prior.Beta.NA(last.params$Beta.NA) +
															prior.Beta.Euro(last.params$Beta.Euro) +
															prior.Beta.Enviro(last.params$Beta.Enviro) + 
															prior.alpha(last.params$alpha) + 
															prior.prob.sigma.prime
				if(exp(Posterior.probability.prime- last.params$Posterior.probability) >= runif(1)){
					new.params$sigma <- sigma.prime
					new.params$Posterior.probability <- Posterior.probability.prime
					new.params$sigma.accepted.moves <- last.params$sigma.accepted.moves + 1
				}											
			}
		new.params$sigma.proposed.moves <- last.params$sigma.proposed.moves + 1
		return(new.params)	
	}

	MCMC <- function(data,mcmc.operators){
		#Declare variables
			Posterior.probability <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.NA <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Euro <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Enviro <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)			
			alpha <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			sigma <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.NA.proposed.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Euro.proposed.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Enviro.proposed.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			alpha.proposed.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			sigma.proposed.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.NA.accepted.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Euro.accepted.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			Beta.Enviro.accepted.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			alpha.accepted.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
			sigma.accepted.moves <- numeric(mcmc.operators$number.generations/mcmc.operators$sample.frequency)
	
		#Initialize parameter values
			Beta.NA[1] <- runif(1,-100,100)
			Beta.Euro[1] <- runif(1,-100,100)
			Beta.Enviro[1] <- runif(1,-100,100)			
			alpha[1] <- rexp(1)
			sigma[1] <- rexp(1)
			vector.means <- c(Beta.NA[1]*data$NA.pops + Beta.Euro[1]*data$Euro.pops + Beta.Enviro[1]*data$Enviro)
		
		#Initialize chain
			Posterior.probability[1] <- Likelihood(data$dependent.variable,
												vector.means,
												(alpha[1]*data$genetic.covariance.matrix + 
													sigma[1]*data$identity.matrix)) +
												prior.Beta.NA(Beta.NA[1]) + 
												prior.Beta.Euro(Beta.Euro[1]) + 
												prior.Beta.Euro(Beta.Enviro[1]) + 												
												prior.alpha(alpha[1]) + 
												prior.sigma(sigma[1])
		
			last.params <- list(Beta.NA[1],Beta.Euro[1],Beta.Enviro[1],alpha[1],sigma[1],
								vector.means,Posterior.probability[1],
								Beta.NA.proposed.moves[1],
									Beta.Euro.proposed.moves[1],
										Beta.Enviro.proposed.moves[1],									
											alpha.proposed.moves[1],
												sigma.proposed.moves[1],
								Beta.NA.accepted.moves[1],
									Beta.Euro.accepted.moves[1],
										Beta.Enviro.accepted.moves[1],
											alpha.accepted.moves[1],
												sigma.accepted.moves[1])
			names(last.params) <- c("Beta.NA","Beta.Euro","Beta.Enviro","alpha","sigma",
								"vector.means","Posterior.probability",
								"Beta.NA.proposed.moves",
									"Beta.Euro.proposed.moves",
										"Beta.Enviro.proposed.moves",									
											"alpha.proposed.moves",
												"sigma.proposed.moves",
								"Beta.NA.accepted.moves",
									"Beta.Euro.accepted.moves",
										"Beta.Enviro.accepted.moves",
											"alpha.accepted.moves",
												"sigma.accepted.moves")
			
		#Run the MCMC
			Updates <- list(Update.Beta.NA,Update.Beta.Euro,Update.Beta.Enviro,Update.alpha,Update.sigma)				
				for(i in 2:mcmc.operators$number.generations){
					x <- sample(c(1:length(Updates)),1)
					new.params <- Updates[[x]](data,last.params,mcmc.operators)
			
					if(i%%mcmc.operators$sample.frequency == 0){
						j <- i/mcmc.operators$sample.frequency
							Posterior.probability[j] <- new.params$Posterior.probability
							Beta.NA[j] <- new.params$Beta.NA
							Beta.Euro[j] <- new.params$Beta.Euro
							Beta.Enviro[j] <- new.params$Beta.Enviro							
							alpha[j] <- new.params$alpha
							sigma[j] <- new.params$sigma
							Beta.NA.proposed.moves[j] <- new.params$Beta.NA.proposed.moves
							Beta.Euro.proposed.moves[j] <- new.params$Beta.Euro.proposed.moves
							Beta.Enviro.proposed.moves[j] <- new.params$Beta.Enviro.proposed.moves
							alpha.proposed.moves[j] <- new.params$alpha.proposed.moves
							sigma.proposed.moves[j] <- new.params$sigma.proposed.moves
							Beta.NA.accepted.moves[j] <- new.params$Beta.NA.accepted.moves
							Beta.Euro.accepted.moves[j] <- new.params$Beta.Euro.accepted.moves
							Beta.Enviro.accepted.moves[j] <- new.params$Beta.Enviro.accepted.moves
							alpha.accepted.moves[j] <- new.params$alpha.accepted.moves
							sigma.accepted.moves[j] <- new.params$sigma.accepted.moves						
					}
			
				last.params <- new.params

					if(i%%mcmc.operators$print.frequency == 0){
						print(i)
						print(new.params$Posterior.probability)
					}
				
					if(i%%mcmc.operators$save.frequency == 0){
						save(last.params,
							data,
							mcmc.operators,
							Posterior.probability,
							Beta.NA,Beta.Euro,Beta.Enviro,
							alpha,
							sigma,
							Beta.NA.proposed.moves,
							Beta.Euro.proposed.moves,
							Beta.Enviro.proposed.moves,							
							alpha.proposed.moves,
							sigma.proposed.moves,
							Beta.NA.accepted.moves,
							Beta.Euro.accepted.moves,
							Beta.Enviro.accepted.moves,
							alpha.accepted.moves,
							sigma.accepted.moves,
							file="MCMC.output.Robj")
					}
				}			
	}

					
##########################################################################################
#	Run the MCMC
##########################################################################################
	
	load(list.files()[grepl("dataset",list.files())])
	
	MCMC( 	data = data,
			mcmc.operators = mcmc.operators	)
	




# data <- dependent.variable,NA.pops,Euro.pops,genetic.covariance.matrix,identity.matrix

# last.params <- Beta.NA,Beta.Euro,alpha,sigma,vector.means,Posterior.probability,Beta.NA.proposed.moves,Beta.Euro.proposed.moves,alpha.proposed.moves,sigma.proposed.moves,Beta.NA.accepted.moves,Beta.Euro.accepted.moves,alpha.accepted.moves,sigma.accepted.moves

# mcmc.operators <- Beta.NA.stp,Beta.Euro.stp,alpha.stp,sigma.stp,number.generations,sample.frequency,save.frequency,print.frequency