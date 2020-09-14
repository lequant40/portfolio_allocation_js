/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function randomVariances
*
* @summary Returns a random variances vector.
*
* @description  This function computes a random n by 1 variances vector, 
* using the unidimensional standard normal distribution with positive support to generate
* n independent random variances.
* 
* @param {number} n the row length of the vector to construct, natural integer greater than or equal to 1.
* @param {number} sigma, the optional standard deviation of the normal distribution used, positive real number; defaults to 1.
*
* @return {Matrix_} the computed vector.
*
*/
self.randomVariances = function(n, sigma) {
	//
	var s = sigma == undefined ? 1 : sigma;
	
	// Generate the random standard deviations vector
	return Matrix_.fill(n, 1, function(i,j) { return pnormrnd_(0, s); });	
}


/**
* @function perturbedVariances
*
* @description This function computes a randomly perturbed version of an input variances vector,
* using the perturbation algorithm described in the first reference.
*
* @see <a href="https://jpm.pm-research.com/content/19/2/6">Chopra, Vijay K; Ziemba, William T; The effect of errors in means, variances, and covariances on optimal portfolio choice; Journal of Portfolio Management; Winter 1993; 19, 2</a>
*
* @param {Matrix_|Array<number>} an n by 1 Matrix representing the variances vector to perturb.
* @param {object} opt optional parameters for the perturbation algorithm.
* @param {string} opt.method the method to use to perturb the standard deviations vector, a string either equals to:
* - "multiplicative-noise", in order to perturb the standard deviations vector using independent normal random variables, c.f. the first reference
; defaults to "multiplicative-noise"
* @param {number} opt.sigma in case opt.method is equal to "multiplicative-noise", the standard deviation of the normal random variable,
* a real number; defaults to 0.05
*
* @return {Matrix_} a Matrix object representing the perturbed variances vector.
*
* @example
* perturbedVariances([0.05, 0.01, 0.01]);
* // ~= Matrix_([0.04, 0.01, 0.0])
*
*/
self.perturbedVariances = function(variancesVect, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "multiplicative-noise";
	}
	
	// Decode the parameters
	var method = opt.method;
	if (method != "multiplicative-noise") {
			throw new Error('unsupported perturbation method');
	}
	
	var sigma = opt.sigma;
	if (opt.sigma === undefined) {
		sigma = 0.05;
	}
	
	// Convert variancesVect to matrix format
	var variancesVect = new Matrix_(variancesVect);
	
	// The input matrix must be a vector
	if (!variancesVect.isVector()) {
		throw new Error('input must be a vector');
	}
	
	// Generate the perturbed vector
	var perturbedVariancesVect = variancesVect.elemMap(function(i,j,val) { 
	                                                       var pvar = val * (1 + sigma*normrnd_(0,1)); 
	                                                       while (pvar <= 0) { 
														       pvar = val * (1 + sigma*normrnd_(0,1));
														   }
														   return pvar;
													   })
	
	// Return it
	return perturbedVariancesVect;
}


