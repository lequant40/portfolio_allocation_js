/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function globalMinimumVarianceWeights
*
* @summary Compute the weights of the global minimum variance portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only global minimum variance portfolio of n assets.
*
* Optionally, the following constraints can be added:
* - Partial investment constraint, replacing the full investment constraint
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* This portfolio is the portfolio with the lowest volatility among all the feasible portfolios.
*
* This portfolio is unique and is mean-variance efficient when the covariance matrix 
* of the assets is definite positive.
* 
* To be noted that when opt.optimizationMethod is set to "critical-line", an error might be raised in certain cases
* when the covariance matrix is not definite positive.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.mu the returns of the n assets in the considered universe, array of n real numbers; defaults to an array of zeroes.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm from Markowitz, c.f. the reference.
* - 'gsmo': usage of  generalized sequential minimization optimization algorithm
* - 'automatic': automatic selection of the optimization method to use
; defaults to 'automatic'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm in case opt.optimizationMethod is set to "critical-line",
* a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo the convergence tolerance of the GSMO algorithm used to solve the minimum-variance optimization problem, 
* a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo the maximum number of iterations of the GSMO algorithm used to solve the minimum-variance optimization problem, 
* a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array.<number>} the weights corresponding to the global minimum variance portfolio, array of n real numbers.
*
* @example
* globalMinimumVarianceWeights([[0.0400, 0.0100], [0.0100, 0.0100]], optimizationMethodParams: {epsGsmo: 1e-10, maxIterGsmo: 10000});
* // [0, 1]
*/
self.globalMinimumVarianceWeights = function (sigma, opt) {
	// Initialize the options structure
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}

	// Decode the options
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var mu;
	if (opt.mu == undefined) {
		mu = Matrix_.zeros(nbAssets, 1);
	}
	else {
		mu = new Matrix_(opt.mu);
	}

	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'automatic';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo' && optimizationMethod != 'automatic') {
		throw new Error('unsupported optimisation method');
	}
	
	
	// ----
	
	
	// Compute the efficient frontier
	var efficientFrontier;
	if (optimizationMethod == "automatic") {
		efficientFrontier = new MeanVarianceEfficientFrontierWrapper(mu, sigma, opt);
	}
	else if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	
	// Return the portfolio with the lowest volatility on the computed 
	// efficient frontier.
	var lowestRiskTolerancePortfolioWeights = efficientFrontier.getLowestRiskTolerancePortfolio();

	// Return the computed weights
	return lowestRiskTolerancePortfolioWeights.toArray();
}

