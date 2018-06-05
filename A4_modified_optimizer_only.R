library(ape)
library(maps)
library(phytools)
library(mvMORPH)
library(stats4)
library(bbmle)
library(geiger)



	
 neg_yule_old <-function(birth, the_tree_name) {
	return (-1 * log_birth_death(birth, 0, the_tree_name, yule=TRUE))

	}

neg_log_brownian <- function( params, tree, trait_values, model) {
	beginning_mean = params[1]
	rate = params[2]
	extra = params[3]
	
	
	log_likelihood = log_brownian(beginning_mean, rate, tree, trait_values, model, extra)
	
	cat(paste0("beginning_mean=", beginning_mean, "\nrate=", rate, ", \nextra=", extra, , "\nlog_likeihood" = log_likelihood, "\n\n"))
	return (-1 * log_likelihood)
	}


log_birth_death <- function(birth, death, the_tree, yule=FALSE)
{
	if (death >= birth) { return (-1*10^(8) ) }
	#sources for birthdeath:
	#Nee, S.C. & May, Robert & Harvey, P.H.. (1994). The Reconstructed Evolutionary Process. Philosophical transactions of the Royal Society of London. Series B, Biological sciences. 344. 305-11. 10.1098/rstb.1994.0068. 
	#
	# https://www.researchgate.net/profile/Robert_May5/publication/15261438_The_Reconstructed_Evolutionary_Process/links/00b7d5209f7d36a2e0000000/The-Reconstructed-Evolutionary-Process.pdf 
	# page 5, equation 21
	# ape | birthdeath function
	
	
	
	# birth is birth rate
	# death is death rate
	# the_tree is well, the tree that needs to be computed on.
	num_species <- length(the_tree$tip.label) # number of species = number of species names
	branching_times <- c(NA, branching.times(the_tree)) # add Not Available to beginning

	# num_species is length(phy$tip.label), number of species in tree
	# branching_times is branching.times(phy), waiting times

	if (yule==TRUE) {
		death = 0
		}



	r = birth - death #relative birth rate
	a = death / birth # death/birth rate

#	1 * (lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + 
#				N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a)))
#				
#original function from ape | birthdeath function

	log_likelihood = lfactorial(num_species -1) + (num_species - 2) * log(r) + r  * sum(branching_times[3:num_species]) + num_species*log(1 - a) - 2 * sum( log ( exp (r * branching_times[2:num_species]) - a) ) 
	#actual log likelihood
	
	if (log_likelihood == -1*Inf) {
		log_likelihood = -1* 10^(-8)
		}
 
 return (log_likelihood)
  }


	log_brownian <- function(beginning_mean, rate, tree, trait_values, model = "default", extra = 1) 
	{
		#source: https://lukejharmon.github.io/pcm/chapter4_fitbm/
		#equation 4.5
		
		# note that this function is imperfect. The source for the formula gave the likelihood, not the log likelihood. As such the function should work find for small trees, but on one large tree I tried, it gave the determinant of the covariance matrix being infinite, which results in finite/infinite = 0 likelihood, which is wrong.
		# as such, this should be tested on small trees
		#beginning_mean = mean of z at time 0 (mean of root)
		# rate is sigma squared is rate of change is varaince
		# covariance_matrix is sigma^2 [t1+t2, t1; t1, t1+t3] where t1 is parent length, t2, t3 are child length
		#trait_values is   x is an n*1 vector of trait values for the n tip species in the tree, with species in the same order as C
		
		
		# n is number of species
		
		#x is vector of trait values
		n = length(trait_values)
#		print(rate)
#		print(tree)
		
		#covar_matrix = covariance matrix: phylogenetic variance and covariance given length of tree branches
		
		if (model == "kappa") {
			tree = kappaTree(tree, extra) # modifies the  tree to kappa
		}
		
		covar_matrix = vcv(phy = tree, corr=FALSE)*rate
		
		if (model == "lambda"){
			mult_matrix = matrix(extra, dim(covar_matrix), dim(covar_matrix)) + diag(1-extra, dim(covar_matrix)) 
			covar_matrix = covar_matrix * mult_matrix # multiply all non-diagonal elemenrs by extra+0 aka lambda , multiply diagonal elements by 1 (extra+1-extra = 1)
		}
		else if (model == "delta") {
			power_matrix = matrix(extra, dim(covar_matrix), dim(covar_matrix))
			covar_matrix = covar_matrix ^ power_matrix # get each element of covar to power of delta
			}
		
		
		
		#one_vec = n*1 column matrix made up of 1
		one_vec = (matrix(1:n, nrow=n, ncol=1))*0+1
		#, solve is inverse
		# t is transpose
		part_a = (trait_values - beginning_mean*one_vec)
		part_b = solve(covar_matrix)

		part_c = t(part_a) %*% part_b%*%part_a
		part_c = part_c[1,1] # part c is a 1*n atrix times an n*n matrix times an n*1 matrix is a 1*1 matrix. It should be interpreted as a scalar, not a matrix, so I take the value and `convert' it to a scalar
		final = exp((-1/2)*part_c)/( sqrt( (2*pi)^n * det(covar_matrix) ) ) # is likelihood
		log_final = log(final)
		if (!is.finite(log_final)) {log_final = -1*10^8}
		return(log_final) # converts to log likelihood
		}
