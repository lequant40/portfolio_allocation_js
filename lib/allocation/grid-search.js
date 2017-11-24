/**
 * @file Functions related to grid search portfolios.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function gridSearchWeights
*
* @summary Compute the weights of a portfolio minimizing an arbitrary objective function.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only portfolio
* of n assets minimizing an arbitrary real-valued objective function fct of n real variables defined on the unit simplex of R^n, 
* which is the n-1 dimensional set of R^n containing the points x = (x_1,...,x_n) statisfying sum x_i = 1 and x_i >= 0, i = 1..n, 
* c.f. the first reference.
*
* Since such a portfolio might not be unique, all the weights corresponding to the same minimum value of the function fct 
* are provided in output.
*
* The minimization of the function fct is achieved through one of the following optimisation methods:
* - Deterministic grid search on a grid of rational points belonging to the unit simplex of R^n
*
* To be noted that finding the minimum value(s) of an arbitrary objective function on the unit simplex
* is an NP-hard problem, c.f. the second reference, so that all the optimisation algorithms above are expected to
* be non-polynomial in n.
* 
* @see <a href="https://en.wikipedia.org/wiki/Simplex">Simplex</a>
* @see <a href="http://www.sciencedirect.com/science/article/pii/S0377221707004262">De Klerk, E., Den Hertog, D., Elabwabi, G.: On the complexity of optimization over the standard simplex. Eur. J Oper. Res. 191, 773–785 (2008)</a>
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
*
* @param {number} nbAssets the number of assets in the considered universe, natural integer superior or equal to 1.
* @param {function} fct a real-valued objective function of n real variables to minimize on the unit simplex of R^n,
* which must take as first input argument an array of n real numbers corresponding to the weights w1,...,wn of the n assets 
* in the considered universe and which must return as output a real number.
* @param {Object} opt the optional parameters for the algorithms used by the function.
* @param {string} opt.optimisationMethod the optimisation method to use in order to minimize the function fct, a string either equals to:
* - 'deterministic': usage of a deterministic grid search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the third reference, where k is defined through the parameter opt.rationalGrid.k
* -
* ; defaults to 'deterministic'.
* @param {number} opt.rationalGrid.k the indice k of the k-th rational grid of the unit simplex of R^n to use in case opt.optimisationMethod is equal to 'deterministic', a natural integer greater than or equal to 1; defaults to n.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to
* the weights of a portfolio minimizing the function fct.
*
* @example
* gridSearchWeights(3, function(arr) { return arr[0];}, {optimisationMethod: 'deterministic', rationalGrid: {k: 2}});
* // [[0,1,0],[0,0.5,0.5],[0,0,1]]
*/
self.gridSearchWeights = function (nbAssets, fct, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	opt.optimisationMethod = opt.optimisationMethod || 'deterministic';
	
	// Select the proper optimisation method
	if (opt.optimisationMethod === 'deterministic') {
		// Decode options for the deterministic grid search method
		opt.rationalGrid = opt.rationalGrid || { k: nbAssets };
		var k = opt.rationalGrid.k;
		
		// Call the rational grid search method
		return simplexRationalGirdSearch_(fct, nbAssets, k);
	}
	else {
	    throw new Error('unsupported optimisation method');
	}
}


