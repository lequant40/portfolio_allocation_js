/**
 * @file Misc. optimisation functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.simplexRationalGirdSearch_ = function(fct, n, k) { return simplexRationalGirdSearch_(fct, n, k); }
/* End Wrapper private methods - Unit tests usage only */
 
 
/**
* @function simplexRationalGirdSearch_
*
* @summary Compute the point(s) minimizing a real-valued arbitrary function of several real variables
* defined on the unit simplex, using an exhaustive search algorithm on a grid made of rational points belonging to the unit simplex.
*
* @description This function returns the list of points x = (x_1,...,x_n) belonging to the unit simplex of R^n which 
* minimize an arbitrary real-valued function fct of n real variables defined on the unit simplex of R^n, 
* using an exhaustive search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the reference.

* To be noted that per lemma 1 of the reference, the number of points on such a rational grid is equal to
* factorial(n + k - 1) / (factorial(k - 1) * factorial(n)), i.e., binomial(n+k-1, n-1), 
* so that this method can be of limited use, even for small n.
*
* For instance, n=5 and k=100 already result in 4598126 points on which to evaluate fct.
*
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
*
* @param {function} fct a real-valued function of n real variables defined on the unit simplex of R^n, 
* which must take as first input argument an array of n real numbers corresponding to a point on the unit simplex of R^n 
* and which must return as output a real number.
* @param {number} n the number of variables of the function fct, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to minimize the function fct, a natural integer superior or equal to 1.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to a point of R^n 
* minimizing the function fct on the k-th rational grid of the unit simplex of R^n.
*
*/
function simplexRationalGirdSearch_(fct, n, k) {
	// Initialize the current minimum value and the current list of associated grid points
	var minValue = Number.MAX_VALUE;
	var minValueGridPoints = [];

	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference,
	// using all the k-compositions of the integer n.
	var nextCompositionIterator = new compositionsIterator_(k, n);
	do {
		// Generate a new composition	
		var nextComposition = nextCompositionIterator.next();
		
		// Compute the current rational grid point by normalizing the generated k-composition
		var weights = nextComposition[1];
		for (var i = 0; i < weights.length; ++i) {
			weights[i] = weights[i] / k;
		}
	  
		// Evaluate the function fct at the current grid point
		var fctValue = fct(weights);
	  
		// If the function value at the current grid point is lower than the current minimum value, this value
		// becomes the new minimum value and the current grid point becomes the new (unique for now) associated grid point.
		if (fctValue < minValue) {
			minValue = fctValue;
			minValueGridPoints = [weights];
		}
		// In case of equality of the function value at the current grid point with the current minimum value, 
		// the current grid point is added to the list of grid points associated to the current minimum value.
		else if (fctValue == minValue) {
			minValueGridPoints.push(weights);
		}
		// Otherwise, nothing needs to be done	  
	}
	while (nextComposition[0]);
	
	// Return the list of grid points associated to the minimum value of fct
	return minValueGridPoints;
}

