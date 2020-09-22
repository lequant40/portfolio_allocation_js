/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function mostDiversifiedWeights
*
* @summary Compute the weights of the most diversified portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only
* most diversified portfolio of n assets, defined as the weights which maximizes the diversification ratio of the portfolio.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio of each asset is the same.
*
* To be noted that when opt.optimizationMethod is set to "critical-line", an error might be raised in certain cases
* when the covariance matrix is not definite positive.
*
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1895459">Y. Choueifaty, T. Froidure, J. Reynier, Properties of the Most Diversified Portfolio, Journal of Investment Strategies, Vol.2(2), Spring 2013, pp.49-70.</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.epsVolatility the volatility below which portfolios are not taken into account in the computation, a strictly positive real number; defaults to 1e-4.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm
* - 'gsmo': usage of generalized sequential minimization optimization algorithm
; defaults to 'gsmo'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm in case opt.optimizationMethod is set to "critical-line",
* a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo the convergence tolerance of the GSMO algorithm used to solve the minimum-variance optimization problem, 
* a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo the maximum number of iterations of the GSMO algorithm used to solve the minimum-variance optimization problem, 
* a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
* @param {number} opt.epsVolatility the volatility below which portfolios are not taken into account in the maximum Sharpe Ratio computation, a strictly positive real number; defaults to 1e-4.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.

* @return {Array.<number>} the weights corresponding to the most diversified portfolio, array of n real numbers.
*
* @example
* mostDiversifiedWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // [0.33, 0.67]
*/
self.mostDiversifiedWeights = function (sigma, opt) {
	// Decode the options.
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	
	// Decode the minimum value of the volatility, in case the covariance matrix
	// is semi-definite positive.
	var epsVolatility = opt.epsVolatility;
	if (epsVolatility === undefined) {
		epsVolatility = 1e-4;
	}
	
	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'gsmo';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo') {
		throw new Error('unsupported optimisation method');
	}
	
	// ------

	// The optimization problem related to the most diversified portfolio 
	// is equivalent to the optimization problem related to the maximum Sharpe ratio 
	// portfolio, with returns replaced by volatilies and a null risk-free rate.
	//
	// So, the most diversified portfolio optimization problem will be solved below
	// using a max Sharpe ratio optimization algorithm.
	var vol = new Matrix_(sigma).toCovarianceMatrix().getStandardDeviations(); // assets volatilities
	
	// Compute the efficient frontier
	var efficientFrontier;
	if (optimizationMethod == "critical-line") {
		efficientFrontier = new MeanVarianceEfficientFrontierCla(vol, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(vol, sigma, opt);
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
	
	// Compute the maximum Sharpe ratio associated to the restricted efficient frontier above, 
	// with a null risk-free rate.
	var rf = 0;
	var portfolio = efficientFrontier.computeMaximumSharpeRatioEfficientPortfolio(rf);

	// Return the computed portfolio weights
	return portfolio[0].toArray();
}

