/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function maximumSharpeRatioWeights
*
* @summary Compute the weights of the portfolio maximizing the Sharpe ratio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only portfolio of n assets maximizing the Sharpe ratio, which is defined as the ratio of the 
* portfolio excess return over a constant risk-free rate to the portfolio volatility.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* When it exists, this portfolio is mean-variance efficient and is unique if the covariance matrix is positive definite.
*
* To be noted that when opt.optimizationMethod is set to "critical-line", an error might be raised in certain cases
* when the covariance matrix is not definite positive.
*
* @see <a href="https://doi.org/10.1111/j.1540-6261.1976.tb03217.x">Elton, E. J., Gruber, M. J. and Padberg, M. W. (1976), SIMPLE CRITERIA FOR OPTIMAL PORTFOLIO SELECTION. The Journal of Finance, 31: 1341-1357</a>
* @see <a href="http://dx.doi.org/10.1080/13504860701255292">S. V. Stoyanov , S. T. Rachev & F. J. Fabozzi (2007) Optimal Financial Portfolios, Applied Mathematical Finance, 14:5, 401-436</a>
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {number} rf the risk-free rate, a real number.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.epsVolatility the volatility below which portfolios are not taken into account in the maximum Sharpe Ratio computation, a strictly positive real number; defaults to 1e-4.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm from Markowitz, c.f. the reference.
* - 'gsmo': usage of  generalized sequential minimization optimization algorithm
* - 'automatic': automatic selection of the optimization method to use
; defaults to 'automatic'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm, a strictly positive natural integer
* or -1 to force an infinite number of iterations; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo in case opt.optimizationMethod is set to "gsmo", the convergence tolerance of the GSMO algorithm used to 
* solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo in case opt.optimizationMethod is set to "gsmo", the maximum number of iterations of the GSMO algorithm 
* used to solve the mean-variance optimization problem, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights corresponding to the portfolio maximizing the Sharpe ratio, array of n real numbers.
*
* @example
* maximumSharpeRatioWeights([0.1, 0.2], [[1, 0.3], [0.3, 1]], 0)
* // [~0.19, ~0.81]
*/
self.maximumSharpeRatioWeights = function(mu, sigma, rf, opt) {
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}

	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'automatic';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo' && optimizationMethod != 'automatic') {
		throw new Error('unsupported optimisation method');
	}
	
	// Decode the minimum value of the volatility, in case the covariance matrix
	// is semi-definite positive.
	var epsVolatility = opt.epsVolatility;
	if (epsVolatility === undefined) {
		epsVolatility = 1e-4;
	}
	
	
	// ------	

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

	
	// Restrict the efficient frontier to the domain of definition
	// of the Sharpe ratio (the portfolios with a strictly positive volatility).
	//
	// To be noted that:
	// - In case the covariance matrix is sufficiently positive definite, 
	// the efficient frontier is not altered, because  no portfolio has a null variance.
	//
	// - In case the covariance matrix is numerically semi-positive definite, the efficient frontier
	// is altered, by moving the minimum variance portfolio to another efficient portfolio
	// located close to it, but with a strictly positive volatility.
	efficientFrontier.restrict("minVolatility", epsVolatility);

	
	// Further restrict the efficient frontier to the domain of strict positivity
	// of the Sharpe ratio.
	//
	// To be noted that the domain of strict positivity of the Sharpe ratio
	// can be empty in case there is no feasible portfolio on the efficient
	// frontier with a strictly positive excess return.
	efficientFrontier.restrict("minReturn", rf);

	
	// Compute the maximum Sharpe ratio efficient portfolio on the restricted
	// efficient frontier.
	var portfolio = efficientFrontier.computeMaximumSharpeRatioEfficientPortfolio(rf);

	
	// Return the computed portfolio weights
	return portfolio[0].toArray();
}
