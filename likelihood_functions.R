library(ape)
library(maps)
library(phytools)
library(mvMORPH)


log_birth_death <- function(birth, death, the_tree)
{
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





	r = birth - death #relative birth rate
	a = death / birth # death/birth rate

#	1 * (lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + 
#				N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a)))
#				
#original function from ape | birthdeath function

	log_likelihood = lfactorial(num_species -1) + (num_species - 2) * log(r) + r  * sum(branching_times[3:num_species]) + num_species*log(1 - a) - 2 * sum( log ( exp (r * branching_times[2:num_species]) - a) ) 
	#actual log likelihood
 
 return (log_likelihood)
  }
  
  
# brownian motion
# 

#log_brownian(beginning_mean, rate, tree, covariance_matrix) {
	log_brownian <- function(beginning_mean, rate, tree, trait_values) 
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
		covar_matrix = vcv(phy = tree, corr=FALSE)*rate
		#one_vec = n*1 column matrix made up of 1
		one_vec = (matrix(1:n, nrow=n, ncol=1))*0+1
		#, solve is inverse
		# t is transpose
		part_a = (trait_values - beginning_mean*one_vec)
		part_b = solve(covar_matrix)

		part_c = t(part_a) %*% part_b%*%part_a
		part_c = part_c[1,1] # part c is a 1*n atrix times an n*n matrix times an n*1 matrix is a 1*1 matrix. It should be interpreted as a scalar, not a matrix, so I take the value and `convert' it to a scalar
		final = exp((-1/2)*part_c)/( sqrt( (2*pi)^n * det(covar_matrix) ) ) # is likelihood
		return(log(final)) # converts to log likelihood
		}
