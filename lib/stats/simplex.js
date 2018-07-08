/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.simplexRationalRounding_ = simplexRationalRounding_;
self.simplexRandomSampler_ = simplexRandomSampler_;
self.simplexDeterministicRationalSampler_ = simplexDeterministicRationalSampler_;
self.simplexRationalGirdSearch_ = simplexRationalGirdSearch_;
self.simplexEuclidianProjection_ = simplexEuclidianProjection_;
self.simplexSparseEuclidianProjection_ = simplexSparseEuclidianProjection_;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function simplexSparseEuclidianProjection_
*
* @summary Returns a closest point on the standard simplex subject to a sparsity constraint.
*
* @description This function computes a closest point (relative to the euclidian distance) 
* with at most k non-zero elements on the standard simplex of R^n to a point x = (x_1,...,x_n) in R^n, 
* using an O(n) implementation of the algorithm 1 of the reference.
*
* In other words, this function computes an at most k-sparse euclidian projection of 
* a point x in R^n onto the standard simplex of R^n.
*
* @see <a href="https://arxiv.org/abs/1206.1529">Anastasios Kyrillidis, Stephen Becker, Volkan Cevher and, Christoph Koch, Sparse projections onto the simplex, arXiv:1206.1529 [cs.LG]</a>
*
* @param {Array.<number>} x, a point belonging to R^n, array of n real numbers.
* @param {number} k, a natural integer strictly greater than one corresponding to the maximum desired sparsity
* (i.e., non-zero elements) of the projected point.
* @return {Array.<number>} the computed closest point to x with at most k non-zero elements, array of n real numbers.
*
* @example
* simplexSparseEuclidianProjection_([0.5, 0.1, 0.2], 1);
* //[[1,0,0]]
*/
function simplexSparseEuclidianProjection_(x, k) {
	// Initializations
	var n = x.length;
	
	
	// Short circuit in case there is no sparsity
	if (k === n) {
		return simplexEuclidianProjection_(x);
	}
	
	
	// Otherwise, compute the support of the projection, i.e., the k largest elements of x.
	
	// Initialize the indexes of the elements of x
	var idx = typeof UInt32Array === 'function' ? new UInt32Array(n) : new Array(n);
	for (var i = 0; i < n; ++i) {
		idx[i] = i;
	}
	
	// Per property of the SELECT algorithm of Floyd and Rivest, the array idx is permuted
	// so that the indexes of the largest k elements of x are at the indexes 0..k-1 of the
	// array idx.
	var compareIndexes = function (a, b) {
	    return x[b] - x[a];
	};
	select_(idx, k, compareIndexes);

	// Extract the k largest elements of x
	var x_k = x.slice(0, k);
	for (var i = 0; i < k; ++i) {
		x_k[i] = x[idx[i]];
	}
	
	
	// Compute the projection on the standard simplex of the k largest elements of x.
	var proj_x_k = simplexEuclidianProjection_(x_k);
	
	
	// Compute the final projection by reconciliating the support of the
	// projection above and its complementary set.
	var y = new Array(n);
	for (var i = 0; i < k;  ++i) {
		y[idx[i]] = proj_x_k[i];
	}
	for (var i = k; i < n;  ++i) {
		y[idx[i]] = 0;
	}
	
	
	// Return the computed projection
	return y;
}

/**
* @function simplexEuclidianProjection_
*
* @summary Returns the closest point on the standard simplex.
*
* @description This function computes the closest point (relative to the euclidian distance)
* lying on the standard simplex of R^n to the input point x = (x_1,...,x_n) in R^n.
*
* In other words, this function computes the euclidian projection of the point x in R^n
* onto the standard simplex of R^n.
*
* Internally, the algorithm used is an O(n) algorithm, c.f. the reference.
*
* @see <a href="https://link.springer.com/article/10.1007/s10107-006-0050-z">Kiwiel, K.C., Breakpoint searching algorithms for the continuous quadratic knapsack problem, Math. Program. (2008) 112: 473</a>
*
* @param {Array.<number>} x a point belonging to R^n, array of n real numbers.
* @return {Array.<number>} the computed closest point to x, array of n real numbers.
*
* @example
* simplexEuclidianProjection_([1, 1, 1]);
* // [~0.33, ~0.33, ~0.33]
*/
function simplexEuclidianProjection_(x) {
	// Initializations
	var n = x.length;
	var zeros = Matrix_.zeros(n, 1);
	var ones = Matrix_.ones(n, 1);
	
	// Convert the problem of the euclidian projection on the standard simplex
	// into the associated instance of the continuous quadratic knapsack problem.	
	var d = ones;
	var a = new Matrix_(x);
	var b = ones;
	var r = 1;
	var l = zeros;
	var u = ones;
		
	// Solve this instance.
	var sol = qksolveBS_(d, a, b, r, l, u);
	var y = sol[0];
	
	// Return the computed projection
	return y.toArray();
}


/**
* @function simplexDeterministicRationalSampler_
*
* @summary Returns a function to generate all the points on a rational grid of the unit simplex of R^n.
*
* @description This function constructs a function to generate all the points on the k-th rational grid
* of the unit simplex of R^n, 1/k * I_n(k), c.f. the first reference.
* 
* The algorithm used internally is based on the enumeration of all the k-compositions of the integer n, c.f. the second reference.
*
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic
*  optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
* @see <a href="https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to generate points, 
* a natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing all 
* the points on the k-th rational grid of the unit simplex of R^n.
*
* @example
* var mySampler = new simplexDeterministicRationalSampler_(3, 10);
* mySampler.sample(); mySampler.sample(); ...; mySampler.sample();
* // [1, 0, 0]; [0.9, 0.1, 0]; ...; null
*/
function simplexDeterministicRationalSampler_(n, k) {
	// Initializations
	this.n = n;
	this.k = k;
	this.compositionIterator = new compositionsIterator_(k, n);
	this.compositionIteratorStatus = true;
	
	/**
	* @function sample
	*
	* @summary Returns a point on the k-th rational grid of the unit simplex of R^n.
	*
	* @description This function generates a point on the k-th rational grid of the unit simplex of R^n.
	*
	* Each call to this function results in the generation of a new point on the k-th rational grid
	* of the unit simplex of R^n, until exhaustion of all such points.
	*
	* @memberof simplexDeterministicRationalSampler_
	* @return {Array.<number>|null} an array of n real numbers corresponding to the coordinates of the generated point in R^n,
	* or null in case all such points have been generated.
	*
	*/
	this.sample = function() {
		// Return null in case there is no more samples to draw
		if (!this.compositionIteratorStatus) {
			return null;
		}
		
		// Generate a new k-composition of n
		var nextComposition = this.compositionIterator.next();
		
		// Compute the current rational grid point by normalizing the generated k-composition
		var x = nextComposition[1];
		var n = x.length;
		for (var i = 0; i < n; ++i) {
			x[i] = x[i] / this.k;
		}

		// Update the internal iterator status
		this.compositionIteratorStatus = nextComposition[0];	
		
		// Return the point being sampled
		return x;
	}
}
 
 
/**
* @function simplexRandomSampler_
*
* @summary Returns a function to compute random points on the unit simplex of R^n.
*
* @description This function constructs a function to compute random points uniformly distributed on
* the unit simplex of R^n, using the algorithm 2 of the reference.
* 
* @see <a href="https://doi.org/10.1016/0377-2217(82)90161-8">R.Y. Rubinstein, Generating random vectors uniformly distributed inside and on 
* the surface of different regions, In European Journal of Operational Research, Volume 10, Issue 2, 1982, Pages 205-209</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing random  
* points on the unit simplex of R^n.
*
* @example
* var mySampler = new simplexRandomSampler_(3);
* mySampler.sample();
* // [0.25, 0, 0.75]
*/
function simplexRandomSampler_(n) {
	// Initializations
	this.n = n;
	this.x = new Array(n); // the coordinates of a point being sampled
	
	/**
	* @function sample
	*
	* @summary Returns a random point on the unit simplex of R^n.
	*
	* @description This function computes a point choosen uniformly at random on the unit simplex of R^n,
	* using the O(n) algorithm 2 of the reference.
	*
	* @memberof simplexRandomSampler_
	* @return {Array.<number>} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample = function() {
		// Computation of n random variables from EXP(1), which will form the basis
		// of the coordinates of the point being sampled	
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from EXP(1) using the inverse method, with no need for the minus sign
			// as the negative sign would cancel out at the subsequent normalization step.
			var u = 1 - Math.random(); // belongs to ]0,1], so that the logarithm below is properly defined
			var e = Math.log(u);

			// Set the i-th coordinate of the point being sampled.
			this.x[i] = e;
			
			// Compute the running sum of the exponential variables, for the subsequent normalization step.
			sum += e;
		}

		// Normalization of the computed coordinates of the point being sampled, so that
		// they all belong to [0,1] and sum to 1.
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/sum;
		}
		
		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice();
	}
}
 
 

/**
* @function simplexRationalRounding_
*
* @summary Compute a closest rational point on the unit simplex.
*
* @description Given a point x = (x_1,...,x_n) on the standard simplex of R^n, this function computes a proximal point xr = (xr_1,...,xr_n) on 
* the r-th rational grid of the unit simplex of R^n, 1/r * I_n(r), with I_n(r) the set of N^n containing the points m = (m_1,...,m_n) 
* statisfying sum m_i = r with r a strictly positive natural integer, so that the computed proximal point xr is one of the closest points to x 
* on this grid with respect to any norm in a large class, c.f. the first reference.
*
* @see <a href="https://doi.org/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258</a>
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

	// Re-order the coordinates according to decreasing values of the fractional parts, 
	// c.f. theorem 1 of the first reference.
	// In case the fractional parts are equal, re-order the coordinates according to 
	// decreasing values of the integer parts, c.f. paragraph 3 of the second reference.
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
* which must take as input argument an array of n real numbers corresponding to a point on the unit simplex of R^n 
* and which must return as output a real number.
* @param {number} n the number of variables of the function fct, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to minimize the function fct, a natural integer superior or equal to 1.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to a point of R^n 
* minimizing the function fct on the k-th rational grid of the unit simplex of R^n.
*
*/
function simplexRationalGirdSearch_(fct, n, k) {
	// Initialize the current minimum value and the current list of associated grid points
	var minValue = Number.POSITIVE_INFINITY;
	var minValueGridPoints = [];

	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference.
	var sampler = new simplexDeterministicRationalSampler_(n, k);
	var weights = sampler.sample();
	while (weights !== null) {  
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
		
		// Generate a new grid point
		weights = sampler.sample();
	}
	
	// Return the list of grid points associated to the minimum value of fct
	return minValueGridPoints;
}
