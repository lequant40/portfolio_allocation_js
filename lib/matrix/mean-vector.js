/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function meanVector
*
* @summary Returns the mean vector of series of values.
*
* @description This function computes the mean vector of series of values.
*
* @see <a href="https://en.wikipedia.org/wiki/Sample_mean_and_covariance">Sample mean and covariance</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1891847">DeMiguel, Victor and Martin-Utrera, Alberto and Nogales, Francisco J., Size Matters: Optimal Calibration of Shrinkage Estimators for Portfolio Selection (July 21, 2011)</a>
*
* @param {Array.<Array.<number>>} arr an array of n arrays of m real numbers, with n and m natural integers
* greater than or equal to 1, with n representing the number of series (or features) and m representing the number of observations per series.
* @param {object} opt optional parameters for the mean vector computation.
* @param {string} opt.method the method to use to compute the mean vector, a string either equals to:
* - "sample-mean", in order to compute the sample mean vector, c.f. the first reference
* - "demiguel-shrinked-mean", in order to compute the optimal convex linear combination of the sample mean vector with a constant shrinkage target vector, c.f. the second reference
; defaults to "sample-mean"
*
* @return {Matrix_} a Matrix object representing the mean vector of the input series of values.
*
* @example
* meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]]);
* // == Matrix_([0.023333333333333334, -0.010000000000000002])
*
*/
self.meanVector = function(arr, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "sample-mean";
	}
	
	// Decode the parameters
	var meanMethod = opt.method;
	if (meanMethod != "sample-mean" && meanMethod != "demiguel-shrinked-mean" ) {
			throw new Error('unsupported mean vector computation method');
	}

	
	//
	var nbSeries = arr.length;
	var nbObservations = arr[0].length;

	
	//
	var means;

	
	// In case the mean vector computation method is "sample-mean",
	// proceed with computing the mean of all series, c.f. the first reference.
	if (meanMethod == "sample-mean") {
		// Compute the sample mean vector
		var meanVector = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]); });
		
		//
		means = meanVector;
	}
	else if (meanMethod == "demiguel-shrinked-mean") {
		// Compute the sample mean vector
		var meanVector = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]); });
		
		// Compute the scaling factor nu, equal to the grand mean
		var nu = mean_(meanVector.toArray());
		
		// Compute the prior, a vector made of nu
		var prior = Matrix_.fill(nbSeries, 1, function(i,j) { return nu; });
		
		// Compute the optimal shrinkage intensity alpha, c.f. formula 6 of the second reference
		var variances = Matrix_.fill(nbSeries, 1, function(i,j) { return covariance_(arr[i-1], arr[i-1]); });
		var sigma_sq_b = variances.sum()/nbSeries;
		var alpha = nbSeries/nbObservations * sigma_sq_b / ( nbSeries/nbObservations * sigma_sq_b + Math.pow(Matrix_.xmy(prior, meanVector).vectorNorm('two'), 2) );

		// Compute the optimally shrinked mean vector,
		// c.f. formula 5 of the second reference.		
		means = Matrix_.axpby(alpha, prior, 1-alpha, meanVector);
	}
	else {
		throw new Error('internal error: unsupported mean vector computation method');
	}

	
	// Return the computed vector
	return means;
}


/**
* @function randomMeanVector
*
* @summary Returns a random mean vector.
*
* @description This function computes a random n by 1 mean vector, 
* using the unidimensional standard normal distribution to generate
* n independent random means.
* 
* @param {number} n the row length of the vector to construct, natural integer greater than or equal to 1.
*
* @return {Matrix_} the computed vector.
*
*/
self.randomMeanVector = function(n) {
	// Generate the random mean vector using the standard normal distribution
	var meanVect = Matrix_.fill(n, 1, function(i,j) { return normrnd_(0, 1); });
	
	// Return it
	return meanVect;
}


/**
* @function perturbedMeanVector
*
* @description This function computes a randomly perturbed version of an input mean vector,
* using the perturbation algorithm described in the first reference.
*
* @see <a href="https://jpm.pm-research.com/content/19/2/6">Chopra, Vijay K; Ziemba, William T; The effect of errors in means, variances, and covariances on optimal portfolio choice; Journal of Portfolio Management; Winter 1993; 19, 2</a>
*
* @param {Matrix_} an n by 1 Matrix representing the mean vector to perturb.
* @param {object} opt optional parameters for the perturbation algorithm.
* @param {string} opt.method the method to use to perturb the mean vector, a string either equals to:
* - "chopra-ziemba", in order to perturb the mean vector using independent standard normal random variables, c.f. the first reference
; defaults to "chopra-ziemba"
* @param {number} opt.k in case opt.method is equal to "chopra-ziemba", the quantity by which to multiply the standard normal random variable,
* a real number; defaults to 0.05
*
* @return {Matrix_} a Matrix object representing the perturbed mean vector.
*
* @example
* perturbedMeanVector([0.05, 0.01, 0.01]);
* // ~= Matrix_([0.04, 0.01, 0.0])
*
*/
self.perturbedMeanVector = function(meanVect, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "chopra-ziemba";
	}
	
	// Decode the parameters
	var method = opt.method;
	if (method != "chopra-ziemba") {
			throw new Error('unsupported perturbation method');
	}
	
	var k = opt.k;
	if (opt.k === undefined) {
		k = 0.05;
	}
	
	// The input matrix must be a vector
	if (!meanVect.isVector()) {
		throw new Error('input must be a vector');
	}
	
	// Generate the perturbed vector
	var perturbedMeanVect = meanVect.elemMap(function(i,j,val) { return val * (1 + k*normrnd_(0,1));})
	
	// Return it
	return perturbedMeanVect;
}



