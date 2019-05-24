/**
 * @file Functions related to grid search portfolios.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function numericalOptimizationWeights
*
* @summary Compute the weights of a portfolio minimizing an arbitrary objective function.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only portfolio
* of n assets minimizing an arbitrary real-valued objective function fct of n real variables defined on the unit simplex of R^n, 
* which is the n-1 dimensional set of R^n containing the points x = (x_1,...,x_n) statisfying sum x_i = 1 and x_i >= 0, i = 1..n, 
* c.f. the first reference.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset to include in the portfolio
* - Maximum weight of each asset to include in the portfolio
*
* Since such a portfolio might not be unique, all the weights corresponding to the same minimum value of the function fct 
* are provided in output.
*
* The minimization of the function fct is achieved through one of the following numerical optimization methods:
* - Grid search on a grid of rational points belonging to the unit simplex of R^n
*
* To be noted that finding the minimum value(s) of an arbitrary real-valued objective function on the unit simplex
* is an NP-hard problem, c.f. the second reference, so that all exact optimization algorithms for this problem 
* are expected to be non-polynomial in n.
* 
* @see <a href="https://en.wikipedia.org/wiki/Simplex">Simplex</a>
* @see <a href="http://www.sciencedirect.com/science/article/pii/S0377221707004262">De Klerk, E., Den Hertog, D., Elabwabi, G.: On the complexity of optimization over the standard simplex. Eur. J Oper. Res. 191, 773–785 (2008)</a>
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
*
* @param {number} nbAssets the number of assets in the considered universe, natural integer superior or equal to 1.
* @param {function} fct a real-valued objective function of n real variables to minimize on the unit simplex of R^n,
* which must take as first input argument an array of n real numbers corresponding to the weights w1,...,wn of the n assets 
* in the considered universe and which must return as output a real number.
* @param {Object} opt parameters for the numerical optimization algorithm.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @param {string} opt.optimizationMethod the optimization method to use in order to minimize the function fct, a string either equals to:
* - 'grid-search': usage of a grid search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the third reference, where k is defined through the parameter opt.optimizationMethodParams.k
; defaults to 'grid-search'.
* @param {number} opt.optimizationMethodParams.k the indice k of the k-th rational grid of the unit simplex of R^n to use in case opt.optimizationMethod is equal to 'grid-search', a natural integer greater than or equal to 1; defaults to n.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to
* the weights of a portfolio minimizing the function fct.
*
* @example
* numericalOptimizationWeights(3, function(arr) { return arr[0];}, {optimizationMethodParams: {k: 2}});
* // [[0,1,0],[0,0.5,0.5],[0,0,1]]
*/
self.numericalOptimizationWeights = function (nbAssets, fct, opt) {
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

	// Initialize the options default values
	if (opt.optimizationMethod === undefined) {
		opt.optimizationMethod = 'grid-search';
	}
	if (opt.optimizationMethod === 'grid-search') {
		if (opt.optimizationMethodParams.k === undefined) {
			opt.optimizationMethodParams.k = nbAssets;
		}
	}
	
	// Select the proper optimisation method
	if (opt.optimizationMethod === 'grid-search') {
		// Call the rational grid search method
		return simplexGridSearch_(fct, nbAssets, opt.optimizationMethodParams.k, opt.constraints.minWeights, opt.constraints.maxWeights);
	}
	else {
	    throw new Error('unsupported optimisation method');
	}
}


