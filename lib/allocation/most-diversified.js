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
* The algorithm used internally is a sequential minimization optimization algorithm,
* which is similar to a cyclical coordinate descent algorithm updating 2 coordinates at each iteration, 
* c.f. the fourth reference, and whose convergence is guaranteed as long as the covariance matrix
* is positive semi-definite
*
* In a previous version of the code, the algorithm used internally was a coordinate descent algorithm,
* c.f. the third reference, kept for historical reference.
*
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1895459">Y. Choueifaty, T. Froidure, J. Reynier, Properties of the Most Diversified Portfolio, Journal of Investment Strategies, Vol.2(2), Spring 2013, pp.49-70.</a>
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* @see <a href="https://link.springer.com/article/10.1023/A:1012431217818">Keerthi, S. & Gilbert, E. Convergence of a Generalized SMO Algorithm for SVM Classifier Design Machine Learning (2002) 46: 351.</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-04.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the most diversified portfolio, array of n real numbers.
*
* @example
* mostDiversifiedWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // [0.33, 0.67]
*/
self.mostDiversifiedWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-4;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	
	// Initializations
	var nbAssets = sigma.nbRows;
	var zeros = Matrix_.zeros(nbAssets, 1);
	var infinitys = Matrix_.fill(nbAssets, 1, function(i,j) { return Number.MAX_VALUE; });

	
	// ----
	
	// The long-only and fully invested most diversified portfolio can be recast
	// as the solution to a convex quadratic program (e.g., the associated matrix 
	// is positive semi-definite, since this is a covariance matrix), c.f. the first reference.
	
		// Build the matrix and the vector of the quadratic program
	var Q = sigma;
	var p = zeros;
	
		// Build the linear equality constraint in the recast space:
		// - Full investment when multiplied by the volatility
	var b = sigma.diagonal().elemMap(function(i,j,val) { return Math.sqrt(val); }); // volatilities
	var r = 1;
	
		// Build the bound constraints:
		// - All lower bounds are zero
		// - All upper bounds are infinite (i.e., no upper bound)
	var l = zeros;
	var u = infinitys;
	
		// Solve the quadratic program
	var sol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: eps, maxIter: maxIterations});
	
	
	// ----
	
	// Extract the rescaled computed portfolio weights.
	var weights = sol[0].normalize();

	// Return the (rescaled) computed weights
	return weights.toArray();
}

