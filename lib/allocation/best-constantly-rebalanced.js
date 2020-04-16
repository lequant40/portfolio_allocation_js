/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function bestConstantlyRebalancedWeights
*
* @summary Compute the weights of the best constantly rebalanced portfolio.
* 
* @description This function returns the weights w_1,...,w_n associated to the 
* the best constantly rebalanced portfolio, which is the portfolio of n assets 
* rebalanced at each period of time so that it holds the same proportion of each assets 
* and with the best return in hindsight, c.f. the reference.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* @see <a href="https://doi.org/10.1109/18.485708">T. M. Cover, Erik Ordentlich, Universal portfolios with side information, 
* IEEE Transactions on Information Theory, Volume 42, Issue 2 March 1996</a>
* 
* @param {Array.<Array.<number>>} priceRelatives an array of n arrays of T real numbers, with 
* priceRelatives[i-1][j-1] the ratio of the final price to the initial price of the i-th asset
* for the j-th period of time, i = 1..n, j = 1..T.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-04.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array.<number>} the weights corresponding to best constantly rebalanced portfolio, array of n real numbers.
*
* @example
* bestConstantlyRebalancedWeights([[1 - 0.05, 1 - 0.05], [1, 1], [1 + 0.05, 1 + 0.05]]);
* // [0, 0, 1]
*/
self.bestConstantlyRebalancedWeights = function (priceRelatives, opt) {
	// Initialize the options structure
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}

	// Initialize the options default values
	if (opt.eps === undefined) {
		opt.eps = 1e-04;
	}
	if (opt.maxIter === undefined) {
		opt.maxIter = 10000;
	}
	
	// Decode the options
	var eps = opt.eps;
	var maxIterations = opt.maxIter;
	var lowerBounds = opt.constraints.minWeights;
	var upperBounds = opt.constraints.maxWeights;
	
	
	// ------
	
	
	// Initializations
	var nbAssets = priceRelatives.length; // m in the reference
	var nbPeriods = priceRelatives[0].length; // n in the reference

	
	// ----
	
	
	// The best constantly rebalanced portfolio is a solution to the following
	// smooth constrained concave optimization problem, c.f. formula 6 of the reference:
	//
	// argmax S_n(b) = argmax ( <b/x_1> * ... * <b/x_n> )
	// s.t. sum b_i = 1
	//      l <= b <= u
	//      (i.e., b belongs to a restricted unit simplex)
	//
	// This optimization problem will be solved using a first-order method
	// for convex minimization, using the facts that:
	// - argmax S_n(b) = argmin -S_n(b)
	// - The restricted unit simplex is a convex set	
	// - -S_n(b) is a convex function on the restricted unit simplex (for instance, it
	//   is log-convex)

	// Define the function representing -S_n(b)
	function f(b) {
		// Initialize the placeholder for the price relatives of all the assets at period k
		var x_k = Matrix_.zeros(nbAssets, 1);

		// Computation of the cumulative product of all the portfolio relatives 
		// <b/x_k>, k = 1..nbPeriods
		var prodPortfolioRelatives = 1.0;
		for (var k = 0; k < nbPeriods; ++k) {
			// Extract the price relatives x_k for all the assets for the period k
			x_k = Matrix_.fill(nbAssets, 1, function(i,j) { return priceRelatives[i-1][k]; }, x_k);
			
			// Compute the portfolio relative <b/x_k> and add it to the cumulative product
			prodPortfolioRelatives *= Matrix_.vectorDotProduct(b, x_k);
		}
		
		// Return the computed function value
		return -prodPortfolioRelatives;
	}
	
	// Define the function representing the gradient of the function -S_n(b).
	//
	// By the multiplicative rule, we have for j = 1..m:
	//
	// d S_n /d b_j (b) = sum_k ( x_k,j * ( Prod_i <b/x_i>, i = 1..n, i<>k ) ), k = 1..n,
	// which is a function costly to evaluate for big m/big n.
	//
	// In case the portfolio relatives <b/x_i>, i = 1..n, are non null (general case),
	// an optimized formula is used.
	function gradf(b) {		
		// Initialize the placeholder for the price relatives of all the assets at period k
		var x_k = Matrix_.zeros(nbAssets, 1);
		
		// Preliminary computation of all the portfolio relatives <b/x_k>, k = 1..nbPeriods,
		// as well as their cumulative product (optimized formula only).
		var portfolioRelatives = typeof Float64Array === 'function' ? new Float64Array(nbPeriods) : new Array(nbPeriods);
		var prodPortfolioRelatives = 1.0;
		var nullPortfolioRelatives = false;
		for (var k = 0; k < nbPeriods; ++k) {
			// Extract the price relatives x_k for all the assets for the period k
			x_k = Matrix_.fill(nbAssets, 1, function(i,j) { return priceRelatives[i-1][k]; }, x_k);
			
			// Compute the portfolio relative <b/x_k>
			var b_d_x_k = Matrix_.vectorDotProduct(b, x_k);
			
			// Save the portfolio relative for future usage
			portfolioRelatives[k] = b_d_x_k;

			// Optimized formula only
			if (nullPortfolioRelatives === false) {
				// Add the portfolio relative to the cumulative product
				prodPortfolioRelatives *= b_d_x_k;
			
				// Determine if the portfolio relative is numerically close to 0,
				// in which case the optimized computation of the gradient of the
				// function S_n(b) cannot be used.
				if (Math.abs(b_d_x_k) <= 1e-12) {
					nullPortfolioRelatives = true;
				}
			}
		}
		
		// Compute the nbAssets products Prod_i <b/x_i>, i = 1..nbPeriods, i <> k, with k = 1..nbAssets:
		//
		// - Optimized formula: if there is no <b/x_k> such that <b/x_k> ~= 0, then the formula
		// Prod_i <b/x_i>, i = 1..nbPeriods / <b/x_k> == Prod_i <b/x_i>, i = 1..nbPeriods, i <> k
		// is used.
		//
		// - Non-optimized formula: Unsupported for now.
		var partialPortfolioRelatives;
		if (nullPortfolioRelatives === false) {
			partialPortfolioRelatives = Matrix_.fill(nbPeriods, 1, 
			                                         function(i,j) { return prodPortfolioRelatives / portfolioRelatives[i-1]; });
		}
		else {
			throw new Error('null portfolio relative detected, unsuported case');
		}

		// Initialize the placeholder for the all the price relatives of the asset k
		var xx_k = Matrix_.zeros(nbPeriods, 1);
		
		// Compute the gradient
		var res = Matrix_.zeros(nbAssets, 1);
		for (var k = 0; k < nbAssets; ++k) {
			// Extract all the price relatives xx_k of the asset k
			xx_k = Matrix_.fill(nbPeriods, 1, function(i,j) { return priceRelatives[k][i-1]; }, xx_k);
			
			// Compute the k-th coordinate of the gradient
			res.data[k] = -Matrix_.vectorDotProduct(xx_k, partialPortfolioRelatives);
		}

		// Return the computed gradient
		return res;
	}

	// Define the characteristic function of the restricted unit simplex
	function g(b) {
		return simplexCharacteristicFunction_(b.data, lowerBounds, upperBounds);
	}
	
	// Define the proximal function associated to g, which is the orthogonal
	// projection on the restricted simplex.
	function proxg(b) {
		return new Matrix_(simplexEuclidianProjection_(b.data, lowerBounds, upperBounds));
	}
	
	// Define the initial point as the projection of the 1 vector 
	// on the restricted unit simplex.
	var x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(nbAssets, 1).data, lowerBounds, upperBounds));

	// Solve the convex optimization program above
	var sol = ccpsolveFISTA_(f, gradf, g, proxg, x0, {eps: eps, maxIter: maxIterations, maxLine: maxIterations});

	
	// ----
	
	
	// Extract the solution, which is the computed portfolio weights
	var weights = sol[0];

	// Return the computed weights
	return weights.toArray();
}

