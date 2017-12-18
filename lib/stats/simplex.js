/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.simplexRationalRounding_ = function(x, r) { return simplexRationalRounding_(x, r); }
self.simplexSampler = simplexSampler_;
self.simplexRationalSampler = simplexRationalSampler_;
/* End Wrapper private methods - Unit tests usage only */

//* @description This function returns the list of points x = (x_1,...,x_n) belonging to the unit simplex of R^n which 
//* minimize an arbitrary real-valued function fct of n real variables defined on the unit simplex of R^n, 
//* using an exhaustive search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the reference.
// * @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
	
// The algorithm used to generated the rational points is to use all the k-compositions of the integer n
// * @see <a href="*https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference,
	// using all the k-compositions of the integer n.
//	* @param {number} n the number of variables of the function fct, natural integer superior or equal to 1.
//* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to minimize the function fct, a natural integer superior or equal to 1.
/*
* @param {number} n a non-negative integer whose composition are desired.
* @param {number} k a non-negative integer indicating the number of parts of desired composition of n.
* @return {function} a function to be used as a generator through its .nextDeterministicSample() method, until it returns null.
*
* @example
* var mySampler = new simplexRationalSampler_(3, 6);
* mySampler.nextDeterministicSample(); mySampler.nextDeterministicSample();
* // [1,0,0]; ~[0.833,0.167,0];
*/
function simplexRationalSampler_(n, k) {
	// Initialize n and x, the coordinates of the point being sampled.
	this.n = n;
	this.k = k;
	this.compositionIterator = new compositionsIterator_(k, n);
	this.compositionIteratorStatus = true;
	
	/**
	* @function sample
	*
	* @summary Returns .
	*
	* @description This function computes the next k-composition of n.
	*
	* The initial k-composition computed by the first call to this function is n00...0, and each subsequent call to 
	* this function will result in a new k-composition until the final k-composition 00...0n is reached.
	*
	* A subsequent call to this function when the final k-composition has been reached will result in
	* the recomputation of all the k-compositions of n, starting from the initial k-composition.
	*
	* @memberof simplexUniformSampler_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one k-composition of n remains to be computed
	* and arr[1] an array of k elements containing the computed k-composition of n.
	*
	*/
	this.nextDeterministicSample = function() {
		// Return null in case there is no more samples to draw
		if (!this.compositionIteratorStatus) {
			return null;
		}
		
		// Generate a new k-composition
		var nextComposition = this.compositionIterator.next();
		
		// Compute the current rational grid point by normalizing the generated k-composition
		var x = nextComposition[1];
		for (var i = 0; i < x.length; ++i) {
			x[i] = x[i] / k;
		}

		// Update the iterator status
		this.compositionIteratorStatus = nextComposition[0];	
		
		// Return the point being sampled
		return x;
	}
}
 
// R.Y. Rubinstein, Generating random vectors uniformly distributed inside and on the surface of different regions, In European Journal of Operational Research, Volume 10, Issue 2, 1982, Pages 205-209
// https://doi.org/10.1016/0377-2217(82)90161-8
// Algorithm 2 of the reference
//  k returns a random variable being uniformly distributed over a standard simplex of dimension k.
function simplexSampler_(n) {
	// Initialize n and x, the coordinates of the point being sampled.
	this.n = n;
	this.x = new Array(n);
	
	/**
	* @function sample
	*
	* @summary Returns .
	*
	* @description This function computes the next k-composition of n.
	*
	* The initial k-composition computed by the first call to this function is n00...0, and each subsequent call to 
	* this function will result in a new k-composition until the final k-composition 00...0n is reached.
	*
	* A subsequent call to this function when the final k-composition has been reached will result in
	* the recomputation of all the k-compositions of n, starting from the initial k-composition.
	*
	* @memberof simplexUniformSampler_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one k-composition of n remains to be computed
	* and arr[1] an array of k elements containing the computed k-composition of n.
	*
	*/
	this.nextRandomSample = function() {
		// Computation of n random variables from EXP(1), which will form the basis
		// of the coordinates of the point being sampled	
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from EXP(1) using the inverse method, with no need for the minus sign
			// as the negative sign would cancel out at the subsequent normalization step
			var u = 1 - Math.random(); // belongs to ]0,1], so that the logarithm below is properly defined
			var e = Math.log(u);

			// Set the i-th coordinate of the point being sampled
			this.x[i] = e;
			
			// Compute the running sum of the exponential variables, for the subsequant normalization step
			sum = sum + e;	
		}

		// Normalization of the computed coordinates of the point being sampled, so that
		// they all belong to [0,1] and sum to 1
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/sum;
		}
		
		// Return a copy of the point being sampled, so that callers can alter it
		return this.x.slice();
	}
}
 
 

/**
* @function simplexRationalRounding_
*
* @summary Compute a rational point on the unit simplex that is the closest to a point on the unit simplex.
*
* @description Given a point x = (x_1,...,x_n) on the standard simplex of R^n, this function computes a proximal point xr = (xr_1,...,xr_n) on 
* the r-th rational grid of the unit simplex of R^n, 1/r * I_n(r), with I_n(r) the set of N^n containing the points m = (m_1,...,m_n) 
* statisfying sum m_i = r with r a strictly positive natural integer, so that the computed proximal point xr is one of the closest points to x 
* on this grid with respect to any norm in a large class, c.f. the first reference.
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258</a>
* @see <a href="https://arxiv.org/abs/1501.00014">Rama Cont, Massoud Heidari, Optimal rounding under integer constraints</a>
* 
* @param {Array.<number>} x a point belonging to the standard simplex of R^n, array of n real numbers.
* @param {number} r the indice of the rational grid of the unit simplex of R^n, 1/r * I_n(r), on which to compute the closest point to x, natural integer greater than or equal to 1.
* @return {Array.<number>} the computed closest point to x on the r-th rational grid of the unit simplex of R^n, array of n real numbers.
*
* @example
* simplexRationalRounding_([0.5759, 0.0671, 0.3570], 20);
* // [0.6, 0.05, 0.35]
*/
function simplexRationalRounding_(x, r) {
	// TODO: Checks, if enabled

	// Compute the integer and fractional parts of the coordinates of the input point multiplied by r, 
	// as described in paragraph 2 of the first reference.
	// In parallel, compute k as also defined in paragraph 2 of the first reference. 
	var k = r;
	var xPartsWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		//
		var rx = r * x[i];
		var integerPart = Math.floor(rx);
		var fractionalPart = rx - integerPart;
		
		//
		k -= integerPart;
		
		//
		xPartsWithIndexes[i] = [integerPart, fractionalPart, i];
	}

	// Re-order the coordinates according to decreasing values of the fractional parts, c.f. theorem 1 of the first reference.
	// In case the fractional parts are equal, re-order the coordinates according to decreasing values of the integer parts,
	// c.f. paragraph 3 of the second reference.
	xPartsWithIndexes.sort(function(a, b) {
		if (b[1] < a[1]) {
			return -1;
		}
		else if (b[1] > a[1]) {
			return 1;
		}
		else { // coordinates have equal fractional parts
			return b[0] - a[0];
		}
	}); 

	// Rounding rule: round the k largest fractional parts up to one and all other fractional parts down to zero,
	// as described in paragraph 2 of the first reference.
	var xr = new Array(x.length);
	for (var i = 0; i < k; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = (xPartsWithIndexes[i][0] + 1) / r;
	}
	for (var i = k; i < xPartsWithIndexes.length; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = xPartsWithIndexes[i][0] / r;
	}

	// Return the computed point
	return xr;
}
