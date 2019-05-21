/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.simplexRationalRounding_ = simplexRationalRounding_;
self.simplexRandomSampler_ = simplexRandomSampler_;
self.simplexDirectionRandomSampler_ = simplexDirectionRandomSampler_;
self.simplexGridSampler_ = simplexGridSampler_;
self.simplexGridSearch_ = simplexGridSearch_;
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
	var y = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
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
* @see <a href="https://link.springer.com/article/10.1007/s10107-006-0050-z">Kiwiel, K.C., Breakpoint searching algorithms 
* for the continuous quadratic knapsack problem, Math. Program. (2008) 112: 473</a>
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
* @function simplexGridSampler_
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
* @param {boolean} useArrayCopy an optional boolean that can be set to false to re-use the same output array throughout
* all the computations (this improves the performances, but requires the caller to NOT alter the output array); defaults to true.
* @return {function} a function to be used through its .sample() method, computing all 
* the points on the k-th rational grid of the unit simplex of R^n.
*
* @example
* var mySampler = new simplexGridSampler_(3, 10);
* mySampler.sample(); mySampler.sample(); ...; mySampler.sample();
* // [1, 0, 0]; [0.9, 0.1, 0]; ...; -1
*/
function simplexGridSampler_(n, k, useArrayCopy) {
	// Initializations
	this.n = n;
	this.k = k;
	
	this.useArrayCopy = useArrayCopy;
	
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.compositionIterator = new compositionsIterator_(k, n, false); // use no array copy in the compositions generation to improve performances
	
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
	* or -1 in case all such points have been generated.
	*/
	this.sample = function() {
		// Generate a new k-composition of n
		var x = this.compositionIterator.next();

		// Return -1 in case there is no more samples to draw
		if (x == -1) {
			return -1;
		}
		
		// Otherwise, compute the current rational grid point by normalizing the generated k-composition
		var n = x.length;
		for (var i = 0; i < n; ++i) {
			this.x[i] = x[i] / this.k;
		}

		// Return either the point being sampled, or a copy of the point being sampled so that callers can alter it
		if (this.useArrayCopy) {
			return this.x.slice(0);
		}
		else {
			return this.x;
		}
	}
}
 
 

/**
* @function simplexRandomSampler_
*
* @summary Returns a function to compute random points on the unit simplex of R^n,
* possibly subject to additional lower bounds and upper bounds constraints on their coordinates.
*
* @description This function constructs a function to compute random points uniformly distributed on either:
* - the unit simplex of R^n, using the algorithm 2 of the first reference
* - the unit simplex of R^n subject to additional lower bounds and upper bounds constraints, i.e. 
* {(x_1,...,x_n), sum x_i = 1 and 0 <= l_i <= x_i <= u_i <= 1, i = 1..n}, where l_i and u_i, i = 1..n are real numbers, 
* using the theorem 1 of the second reference
* 
* @see <a href="https://doi.org/10.1016/0377-2217(82)90161-8">R.Y. Rubinstein, Generating random vectors uniformly distributed inside and on 
* the surface of different regions, In European Journal of Operational Research, Volume 10, Issue 2, 1982, Pages 205-209</a>
* @see <a href="https://doi.org/10.1016/S0167-7152(99)00095-4">Kai-Tai Fang, Zhen-Hai Yang, On uniform design of experiments with restricted
* mixtures and generation of uniform distribution on some domains, Statistics & Probability Letters, Volume 46, Issue 2, 
* 15 January 2000, Pages 113-120</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @param {Array.<number>} l the optional lower bounds contraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
* @param {Array.<number>} u the optional upper bounds contraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
* @param {function} rnd an optional random vector generator in the unit hypercube of R^n-1, a function taking no input argument and
* returning an array of n-1 points each belonging to the interval (0, 1).
* @return {function} a function computing the random points to be used through its .sample() method.
*
* @example
* var mySampler = new simplexRandomSampler_(3);
* mySampler.sample();
* // [0.25, 0, 0.75]
*
* var mySampler = new simplexRandomSampler_(3, [0.1, 0.2, 0.3], [1,1,1]);
* mySampler.sample();
* // [0.20, 0.20, 0.60]
*/
function simplexRandomSampler_(n, l, u, rnd) {
	// Initializations
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.u = typeof Float64Array === 'function' ? new Float64Array(n-1) : new Array(n-1); // the coordinates of a random point in the unit hypercube of R^n-1

	// Initialization of the random number generator
	//
	// By default, it generates points uniformly at random in the unit hypercube of R^n-1
	this.randomGenerator = function(arr) {
		for (var i = 0; i < this.n; ++i) {
			arr[i] = Math.random();
		}
	}
	if (rnd) {
		if (typeof rnd !== "function") {
			throw new Error('the random number generator must be a function');
		}
		this.randomGenerator = rnd;
	}

	// Misc. checks on the optional lower and upper bounds
	if (l) {
		this.lowerBounds = l; // no input checks
		
		if (!u) {
			this.upperBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			
			// Initialization to an array of ones
			for (var i = 0; i < this.n; ++i) {
				this.upperBounds[i] = 1;
			}
		}
	}
	if (u) {
		this.upperBounds = u; // no input checks
		
		if (!l) {
			this.lowerBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);

			// Initialization to an array of zeros
			for (var i = 0; i < this.n; ++i) {
				this.lowerBounds[i] = 0;
			}
		}
	}
	
	// Misc. computations in case the simplex is being sampled through the theorem 1
	// of the second reference:
	//
	// - Feasibility checks on the restricted simplex:
	// -- Lower bounds and upper bounds must satisfy 0 <= l_i <= u_i <= 1, i = 1..n
	// -- Lower bounds and upper bounds must satisfy sum l_i < 1 < sum u_i, c.f. formula 2.2 of the second reference
	//    (In case sum l_i == 1 or sum u_i == 1, the computations are prematurly stopped, because the restricted
	//    simplex is then equal to a point)
	//
	// - Deletion of possible superflous constraints, as described in formula 2.3 of the second reference
	//
	// - Computation of the upper bounds*, as defined after formula 2.3' of the second reference
	//
	// - In parallel of the three steps above, computation of lower, upper and upper* bounds sums, 
	// as well as the upper* bounds running sum
	if (this.lowerBounds || this.upperBounds) {
		// Feasibility checks
		var sumLowerBounds = 0;
		var sumUpperBounds = 0;
		for (var i = 0; i < this.n; ++i) {
			var lowerBound = this.lowerBounds[i];
			var upperBound = this.upperBounds[i];
			
			// Check on lower and upper bounds l_i and u_i
			if (lowerBound > upperBound) {
				throw new Error('infeasible problem detected: lower bound strictly greater than upper bound');
			}
			if (lowerBound < 0) {
				throw new Error('incorrect problem detected: lower bound strictly lower than 0');
			}
			if (upperBound > 1) {
				throw new Error('incorrect problem detected: upper bound strictly greater than 1');
			}

			// Compute the running sum of lower and upper bounds, for subsequent feasibility check
			sumLowerBounds += lowerBound;
			sumUpperBounds += upperBound;
		}
		this.sumLowerBounds = sumLowerBounds;
		this.sumUpperBounds = sumUpperBounds;
		
		if (this.sumLowerBounds > 1 || this.sumUpperBounds < 1) {
			throw new Error('infeasible problem detected: the restricted simplex is empty');
		}
		else if (this.sumLowerBounds == 1 || this.sumUpperBounds == 1) {
			// Nothing to do
		}
		else {
			// Deletion of possible superflous constraints, replacing the lower and upper bounds
			var updatedSumLowerBounds = 0;
			var updatedSumUpperBounds = 0;
			for (var i = 0; i < this.n; ++i) {
				var lowerBound = this.lowerBounds[i];
				var upperBound = this.upperBounds[i];
				
				var updatedLowerBound = Math.max(lowerBound, upperBound + 1 - this.sumUpperBounds);
				var updatedUpperBound = Math.min(upperBound, lowerBound + 1 - this.sumLowerBounds);
				
				this.lowerBounds[i] = updatedLowerBound;
				this.upperBounds[i] = updatedUpperBound;

				updatedSumLowerBounds += updatedLowerBound;
				updatedSumUpperBounds += updatedUpperBound;
			}
			this.sumLowerBounds = updatedSumLowerBounds;
			this.sumUpperBounds = updatedSumUpperBounds;
			
			// Computation of upper* bounds
			this.upperStarBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			this.runningSumUpperStarBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			
			var sumUpperStarBounds = 0;
			for (var i = 0; i < this.n; ++i) {
				this.upperStarBounds[i] = (this.upperBounds[i] - this.lowerBounds[i]) / (1-this.sumLowerBounds);
				
				sumUpperStarBounds += this.upperStarBounds[i];
				this.runningSumUpperStarBounds[i] = sumUpperStarBounds;
			}
		}
	}


	/**
	* @function sample_no_bounds
	*
	* @summary Returns a random point on the unit simplex of R^n.
	*
	* @description This function computes a point choosen uniformly at random on the unit simplex of R^n,
	* using the O(n) algorithm 2 of the first reference.
	*
	* @memberof simplexRandomSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample_no_bounds = function() {
		// Computation of n independent random variables from EXP(1), which will form the basis
		// of the coordinates of the point being sampled	
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from EXP(1) using the inverse method, with no need for the minus sign
			// as the negative sign would cancel out at the subsequent normalization step.
			var u = 1 - Math.random(); // belongs to ]0,1], so that the logarithm below is properly defined
			var e = Math.log(u); // e ~ EXP(1)

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
		return this.x.slice(0);
	}

	/**
	* @function sample_binding_bounds
	*
	* @summary Returns a random point on the unit simplex of R^n subject to exact bounds
	* on its coordinates.
	*
	* @description This function computes a point choosen uniformly at random on the unit simplex of R^n,
	* subject to an exact bounds constaints on its coordinates, which makes the point unique and non random.
	*
	* @memberof simplexRandomSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample_binding_bounds = function() {
		// Determination of the binding bounds
		var bindingBounds = null;
		if (this.sumLowerBounds == 1) {
			bindingBounds = this.lowerBounds;
		}
		else if (this.sumUpperBounds == 1) {
			bindingBounds = this.upperBounds;
		}
		else {
			throw new Error('internal error');
		}

		// Generation of a point on the restricted simplex, with its coordinates equal to the binding bounds
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = bindingBounds[i];
		}
		
		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice(0);
	}
	
	/**
	* @function sample_with_bounds
	*
	* @summary Returns a random point on the unit simplex of R^n subject to additional 
	* lower bounds and upper bounds constraints on its coordinates.
	*
	* @description This function computes a point choosen uniformly at random on the unit simplex of R^n,
	* subject to additional lower bounds and upper bounds constraints on its coordinates, using the algorithm
    * adapted from the theorem 1 of the second reference.
	*
	* @memberof simplexRandomSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample_with_bounds = function() {
		// Generate a point in the unit hypercube of R^n-1
		this.randomGenerator(this.u);
		var u = this.u;

		// Use the theorem 1 of the second reference to generate a random point on T_n(0,upperStarBounds)
		var delta_k = 1;
		for (var k = n; k >= 2; --k) {
			// In case delta_k is numerically null, it means all the remaining coordinates
			// of the random point being genetared on T_n(0,upperStarBounds) must be set to zero.
			//
			// The main loop can then be prematurly stopped.
			if (Math.abs(delta_k) <= 1e-14) {
				for (var kk = k; kk >= 1; --kk) {
					this.x[kk-1] = 0;
				}
				break;
			}

			// Retrieve the k-1th coordinate of the random vector u
			var u_k = u[k-2];

			// Compute the function G of theorem 1 of the second reference
			var d_k = Math.max(0, 1 - this.runningSumUpperStarBounds[k-2]/delta_k);			
			var phi_k = Math.min(1, this.upperStarBounds[k-1]/delta_k);			
			var y_k = delta_k * (1 - Math.pow(u_k*Math.pow(1-phi_k, k-1) + (1-u_k)*Math.pow(1-d_k, k-1), 1/(k-1)));
			
			// Update the k-th coordinate of the generated random point
			this.x[k-1] = y_k;
			
			// Prepare for the next iteration
			delta_k -= y_k;			
		}
		if (k == 1) { // In case the main loop above exited with delta not numerically null
			this.x[0] = delta_k; // Update the 1st coordinate of the generated random point
		}

		// Use the linear mapping described after formula 2.3' of the second reference to map
		// the random point above from T_n(0,upperStarBounds) to T_n(lowerBounds,upperBounds).
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = (1-this.sumLowerBounds)*this.x[i] + this.lowerBounds[i];
		}

		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice(0);
	}
	

	// Definition of the sampling method
	if (this.lowerBounds || this.upperBounds) {
		if (this.sumLowerBounds == 1 || this.sumUpperBounds == 1) {
			this.sample = this.sample_binding_bounds;
		}
		else {
			this.sample = this.sample_with_bounds;
		}
	}
	else {
		this.sample = this.sample_no_bounds;
	}
}


/**
* @function simplexDirectionRandomSampler_
*
* @summary Returns a function to compute random directions to be used on the unit simplex of R^n.
*
* @description This function constructs a function to compute random unit directions uniformly distributed on
* the intersection of the unit hypersphere of R^n and of the hyperplane defined by the equation <(1,1,...,1)/x> = 0,
* using the algorithm 2 of the reference.
* 
* @see <a href="https://projecteuclid.org/euclid.ba/1488337478">Cong, Yulai; Chen, Bo; Zhou, Mingyuan. Fast Simulation of 
* Hyperplane-Truncated Multivariate Normal Distributions. Bayesian Anal. 12 (2017), no. 4, 1017--1037. doi:10.1214/17-BA1052.</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing random  
* directions to be used on the unit simplex of R^n.
*
* @example
* var mySampler = new simplexDirectionSampler_(3);
* mySampler.sample();
* // [-0.5856783494622358, -0.19984292686015526, 0.7855212763223908]
*/
function simplexDirectionRandomSampler_(n) {
	// Initializations
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	
	/**
	* @function sample
	*
	* @summary Returns a random point belonging to the intersection of the unit hypersphere of R^n
	* and of the hyperplane defined by the equation <(1,1,...,1)/x> = 0.
	*
	* @description This function computes a point choosen uniformly at random in the intersection of 
	* the unit hypersphere of R^n and of the hyperplane defined by the equation <(1,1,...,1)/x> = 0,
	* using the O(n) algorithm 2 of the reference.
	*
	* @memberof simplexDirectionSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample = function() {
		// Computation of n independent random variables from N(0,1), which will form the basis
		// of the coordinates of the point being sampled.
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from N(0,1), using the inverse method
			var u = Math.random(); // u ~ U[0,1[
			while (u === 0.0) {
				u = Math.random();
			} // u ~ U]0,1[
			var r = norminv_(u); // r ~ N(0,1)
			
			// Set the i-th coordinate of the point being sampled
			this.x[i] = r;
			
			// Compute the running sum of the normal variables
			sum += r;
		}
		
		// Normalization of the computed coordinates of the point being sampled, so that
		// the associated vector in R^n belongs to the hyperplane <(1,1,...,1)/x> = 0,
		// i.e. sum x_i = 0.
		// - The associated vector in R^n 
		var sum_d_n = sum / this.n;
		var sum_sq = 0;
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i] - sum_d_n;
			
			// Compute the running sum of the squares of the coordinates, for the subsequent normalization step.
			sum_sq += this.x[i] * this.x[i];
		}
		
		// Final normalization of the computed coordinates of the point being sampled, so that
		// the associated vector in R^n belongs to the n-hypersphere, i.e. its 2-norm 
		// is equal to 1.
		var x_two_norm = Math.sqrt(sum_sq);
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/x_two_norm;
		}
	
		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice(0);
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
* @see <a href="https://doi.org/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: 
* Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258</a>
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
* @function simplexGridSearch_
*
* @summary Compute the point(s) minimizing a real-valued function of several real variables
* defined on the unit simplex using a grid search algorithm.
*
* @description This function returns the list of points x = (x_1,...,x_n) belonging to the unit simplex of R^n which 
* minimize a real-valued function fct of n real variables defined on the unit simplex of R^n, 
* using a grid search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), 
* c.f. the reference.
*
* Optionally, lower bounds and upper bounds constraints can be added to the problem, in which case k
* must be chosen so that some of the k-th rational grid points belong to the restricted simplex.
*
* To be noted that per lemma 1 of the reference, the number of points on such a grid is equal to
* factorial(n + k - 1) / (factorial(k - 1) * factorial(n)), i.e., binomial(n+k-1, n-1), 
* so that this method can be of limited use, even for small n.
*
* For instance, n=5 and k=100 already result in 4598126 points on which to evaluate f.
*
* To also be noted that as the optional lower bounds and upper bounds contraints become tigter,
* the volume of the restricted simplex becomes smaller, so that this method can also be of limited use
* because most of the grid points will fall outside of the restricted simplex. 
*
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and 
* quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
* @see <a href="https://doi.org/10.1016/S0167-7152(99)00095-4">Kai-Tai Fang, Zhen-Hai Yang, On uniform design of experiments with restricted
* mixtures and generation of uniform distribution on some domains, Statistics & Probability Letters, Volume 46, Issue 2, 
* 15 January 2000, Pages 113-120</a>
*
* @param {function} f, a function which must take as input argument
* an array of n real numbers corresponding to a point on the unit simplex of R^n and which must return as output a real number 
* corresponding to f(x).
* @param {number} n the number of variables of the function f, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to minimize the function f, 
* a natural integer superior or equal to 1.
* @param {Array.<number>} l the optional lower bounds contraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
* @param {Array.<number>} u the optional upper bounds contraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers
* corresponding to a point of R^n minimizing the function f on the k-th rational grid of the unit simplex of R^n.
*
* @example
* // Minimize f(x,y) = x on the unit simplex of R^2, using the 10-th rational grid
* simplexGridSearch_(function(arr) { return arr[0]; }, 2, 10); 
* // [[0,1]]
*/
function simplexGridSearch_(f, n, k, l, u) {
	// Misc. checks on the optional lower and upper bounds
	var lowerBounds = null;
	var upperBounds = null;
	if (l) {
		lowerBounds = l; // no input checks
		
		if (!u) {
			upperBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			
			// Initialization to an array of ones
			for (var i = 0; i < n; ++i) {
				upperBounds[i] = 1;
			}
		}
	}
	if (u) {
		upperBounds = u; // no input checks
		
		if (!l) {
			lowerBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);

			// Initialization to an array of zeros
			for (var i = 0; i < n; ++i) {
				lowerBounds[i] = 0;
			}
		}
	}
	
	// - Optional feasibility checks on the restricted simplex:
	// -- Lower bounds and upper bounds must satisfy 0 <= l_i <= u_i <= 1, i = 1..n
	// -- Lower bounds and upper bounds must satisfy sum l_i <= 1 <= sum u_i, c.f. formula 2.2 of the second reference
	//    (In case sum l_i == 1 or sum u_i == 1, the computations are prematurly stopped, because the restricted
	//    simplex is then equal to a point)
	if (lowerBounds || upperBounds) {
		// Feasibility checks
		var sumLowerBounds = 0;
		var sumUpperBounds = 0;
		for (var i = 0; i < n; ++i) {
			var lowerBound = lowerBounds[i];
			var upperBound = upperBounds[i];
			
			// Check on lower and upper bounds l_i and u_i
			if (lowerBound > upperBound) {
				throw new Error('infeasible problem detected: lower bound strictly greater than upper bound');
			}
			if (lowerBound < 0) {
				throw new Error('incorrect problem detected: lower bound strictly lower than 0');
			}
			if (upperBound > 1) {
				throw new Error('incorrect problem detected: upper bound strictly greater than 1');
			}

			// Compute the running sum of lower and upper bounds, for subsequent feasibility check
			sumLowerBounds += lowerBound;
			sumUpperBounds += upperBound;
		}
		
		if (sumLowerBounds > 1 || sumUpperBounds < 1) {
			throw new Error('infeasible problem detected: the restricted simplex is empty');
		}
	}
	
	// Initialize the current minimum value and the current list of associated grid points
	var minValue = Number.POSITIVE_INFINITY;
	var minValueGridPoints = [];

	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference.
	var sampler = new simplexGridSampler_(n, k, false); // use no array copy in the simplex grid sampler to improve performances
	var weights = sampler.sample();
	while (weights !== -1) {  
		// Optionally reject the current grid point if it does not belong to the restricted simplex,
		// and generate a new grid point
		var withinBounds = true;
		if (lowerBounds) {
			for (var i = 0; i < n; ++i) {
				if (lowerBounds[i] > weights[i]) {
					withinBounds = false;
					break;
				}
			}
		}
		if (upperBounds) {
			for (var i = 0; i < n; ++i) {
				if (weights[i] > upperBounds[i]) {
					withinBounds = false;
					break;
				}
			}
		}
		if (!withinBounds) {
			weights = sampler.sample();
			continue;
		}
		
		
		// Evaluate the function f at the current grid point
		var fctValue = f(weights);
	  
	  
		// If the function value at the current grid point is lower than the current minimum value, this value
		// becomes the new minimum value and the current grid point becomes the new (unique for now) associated grid point.
		if (fctValue < minValue) {
			minValue = fctValue;
			minValueGridPoints = [weights.slice(0)];
		}
		// In case of equality of the function value at the current grid point with the current minimum value, 
		// the current grid point is added to the list of grid points associated to the current minimum value.
		else if (fctValue == minValue) {
			minValueGridPoints.push(weights.slice(0));
		}
		// Otherwise, nothing needs to be done
		

		// Generate a new grid point
		weights = sampler.sample();
	}
	
	// Return the list of grid points associated to the minimum value of f
	return minValueGridPoints;
}
