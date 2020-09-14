/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.fullSimplexCharacteristicFunction_ = fullSimplexCharacteristicFunction_;
self.simplexCharacteristicFunction_ = simplexCharacteristicFunction_;
self.simplexEmptinessCheck_ = simplexEmptinessCheck_;
self.fullSimplexEmptinessCheck_ = fullSimplexEmptinessCheck_;
self.simplexRationalRounding_ = simplexRationalRounding_;
self.simplexRandomSampler_ = simplexRandomSampler_;
self.simplexDirectionRandomSampler_ = simplexDirectionRandomSampler_;
self.simplexGridSampler_ = simplexGridSampler_;
self.simplexGridSearch_ = simplexGridSearch_;
self.simplexEuclidianProjection_ = simplexEuclidianProjection_;
self.fullSimplexEuclidianProjection_ = fullSimplexEuclidianProjection_;
self.simplexSparseEuclidianProjection_ = simplexSparseEuclidianProjection_;
self.simplexLpSolve_ = simplexLpSolve_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function simplexLpSolve_
*
* @description This function solves a minimization problem of a linear function over the restricted unit simplex of R^n,
* assumed to be provided in the following format:
*
* min f(x) = <c/x>, x in R^n
*
* s.t. sum x_i = 1
*      l <= x <= u (finite bound constraints)
*
* where:
* - c is a an array of n real numbers
* - l and u are two arrays of n real numbers verifying 0 <= l_i <= u_i <= 1, i = 1..n
*
* The internal method used to compute the solution is a greedy method, inspired from 
* the greedy algorithm used to solve the continuous knapsack problem, c.f. the reference.
*
* Note that there can be several solutions to the problem, but all these solutions have
* the same minimal function value.
*
* @see <a href="https://en.wikipedia.org/wiki/Continuous_knapsack_problem">Continuous_knapsack_problem</a>
*
* @param {Array.<number>|Matrix} c, the coefficients of the linear function, array of n real numbers.
* @param {Array.<number>|Matrix} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>|Matrix} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
*
* @return {Array<Object>} an array arr containing three elements: 
* - arr[0] an array of n elements corresponding to a solution x^* to the problem
* - arr[1] the optimal value of the function f, f(x^*) = <c/x^*>
* - arr[2] a boolean indicating if there is a unique optimal solution (true) or not (false)
*
* @example
* simplexLpSolve_([0.1, 0.9]);
* // == [[1, 0], 0.1]
*/
function simplexLpSolve_(c, l, u) {
	// Initializations
	var eps = 1e-16; // the numerical zero
	
	var c = new Matrix_(c);
	var n = c.nbRows;
	var l = l ? new Matrix_(l) : Matrix_.zeros(n, 1);
	var u = u ? new Matrix_(u) : Matrix_.ones(n, 1);

	// Check that the problem is feasible (i.e., that the restricted unit simplex on which
	// the optimization is taking place is not empty.
	simplexEmptinessCheck_(n, l.toArray(), u.toArray());

	// Order the coefficients of the linear function in ascending order w.r.t. their values
	var idx_c = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
	for (var j = 0; j < n; ++j) {		
		idx_c[j] = j + 1;
	}
	idx_c.sort(function(a, b) { 
		return c.getValue(a, 1) - c.getValue(b, 1);
	});
	
	// Initialize the solution with the imposed lower bounds
	var x = new Matrix_(l);
	
	// Starting from the lowest-value coefficient of the linear function, set the 
	// associated coordinates of x to their highest possible value until the sum 
	// of the coordinates of x is equal to one.
	var delta_sum_x = 1 - x.sum();
	var idx_i = -1;
	for (var i = 0; i < n; ++i) {	
		// In case the new delta sum of the coordinates of x is
		// numerically equal to zero, the loop can be stopped.
		if (Math.abs(delta_sum_x) <= eps) {
			break;
		}
		
		// Extract the index of the coefficient of the linear function and its associated x coordinate
		idx_i = idx_c[i];
		var x_idx_i = x.getValue(idx_i, 1);
					
		// Compute the highest possible value for the increment in the x coordinate
		// and set the new value of the x coordinate.
		//
		// At the same time, compute the new delta sum of the coordinates of x for the next iteration.
		if (delta_sum_x >= u.getValue(idx_i, 1) - x_idx_i) {
			x.setValue(idx_i, 1, u.getValue(idx_i, 1));
			delta_sum_x -= u.getValue(idx_i, 1) - x_idx_i;
		}
		else {
			x.setValue(idx_i, 1, x_idx_i + delta_sum_x);
			delta_sum_x = 0;
		}
	}	
	
	// Return the computed solution as well as the associated function value
	return [x.toArray(), Matrix_.vectorDotProduct(c, x)];
}


/**
* @function simplexCharacteristicFunction_
*
* @summary The characteristic function of the restricted unit simplex of R^n.
*
* @description This function is the characteristic function of the restricted unit simplex of R^n,
* c.f. the reference, so that it returns:
* - 0 <=> x belongs to the restricted unit simplex of R^n
* - +oo <=> x does not belong to the restricted unit simplex of R^n
*
* By definition, x belongs to the restricted simplex of R^n if and only if its coordinates satisfy:
* - sum x_i == 1
* - l_i <= x_i <= u_i <= 1, i = 1..n
*
* @see <a href="https://en.wikipedia.org/wiki/Characteristic_function_(convex_analysis)">Characteristic function (convex analysis)</a>
*
* @param {Array.<number>} x, a point belonging to R^n, array of n real numbers.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @return {number} In case the point x belongs to the restricted unit simplex of R^n, returns 0, otherwise, returns Number.POSITIVE_INFINITY.
*
* @example
* simplexCharacteristicFunction_([0.1, 0.9]);
* // == 0
*
* @example
* simplexCharacteristicFunction_([0.1, 1.9]);
* // == Number.POSITIVE_INFINITY
*/
function simplexCharacteristicFunction_(x, l, u) {
	// Initializations
	var n = x.length;
	var eps_tol = 1e-12; // used to numerically determine some conditions


	// Emptiness check on the restricted simplex.
	//
	// In case the restricted simplex is empty, an exception is thrown, so that
	// the process is (violently) stopped here.
	simplexEmptinessCheck_(n, l, u);
	
	
	// Check if the point x belongs to the restricted unit simplex
	var sum_xi = 0;
	for (var i = 0; i < n; ++i) {	
		// Extract the i-th coordinate of the restricted simplex lower bound
		var lb_i = 0;
		if (l) {
			lb_i = l[i];
		}
		
		// Extract the i-th coordinate of the restricted simplex upper bound
		var up_i = 1;
		if (u) {
			up_i = u[i];
		}
		
		// Extract the i-th coordinate of the point x
		var x_i = x[i];

		// Compare the i-th coordinate of the point x with the lower and the upper bounds
		// of the restricted simplex.
		if (x_i < lb_i || x_i > up_i) {
			return Number.POSITIVE_INFINITY;
		}
		
		// Compute the sum of the coordinates of the point x for later use
		sum_xi += x_i;
	}
	
	// Check if the point x belongs to the restricted unit simplex, second step.
	//
	// Note: Due to limited numerical precision, the test on the sum of the coordinates
	// of the point x (which must be equal to 1) cannot be exact.
	if (Math.abs(sum_xi - 1) > eps_tol) {
		return Number.POSITIVE_INFINITY;
	}
	

	// At this stage, the input point x belongs to the restricted unit simplex
	return 0;
}


/**
* @function fullSimplexCharacteristicFunction_
*
* @summary The characteristic function of the restricted unit full simplex of R^n.
*
* @description This function is the characteristic function of the restricted unit full simplex of R^n,
* c.f. the reference, so that it returns:
* - 0 <=> x belongs to the restricted unit full simplex of R^n
* - +oo <=> x does not belong to the restricted unit full simplex of R^n
*
* By definition, x belongs to the restricted unit full simplex of R^n if and only if its coordinates satisfy:
* - sum x_i <= 1
* - l_i <= x_i <= u_i <= 1, i = 1..n
*
* @see <a href="https://en.wikipedia.org/wiki/Characteristic_function_(convex_analysis)">Characteristic function (convex analysis)</a>
*
* @param {Array.<number>} x, a point belonging to R^n, array of n real numbers.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @return {number} In case the point x belongs to the restricted unit full simplex of R^n, returns 0, otherwise, returns Number.POSITIVE_INFINITY.
*
* @example
* fullSimplexCharacteristicFunction_([0.1, 0.8]);
* // == 0
*
* @example
* fullSimplexCharacteristicFunction_([0.1, 1.9]);
* // == Number.POSITIVE_INFINITY
*/
function fullSimplexCharacteristicFunction_(x, l, u) {
	// Initializations
	var n = x.length;
	var eps_tol = 1e-12; // used to numerically determine some conditions


	// Emptiness check on the restricted full simplex.
	//
	// In case the restricted full simplex is empty, an exception is thrown, so that
	// the process is (violently) stopped here.
	fullSimplexEmptinessCheck_(n, l, u);
	
	
	// Check if the point x belongs to the restricted unit full simplex
	var sum_xi = 0;
	for (var i = 0; i < n; ++i) {	
		// Extract the i-th coordinate of the restricted simplex lower bound
		var lb_i = 0;
		if (l) {
			lb_i = l[i];
		}
		
		// Extract the i-th coordinate of the restricted simplex upper bound
		var up_i = 1;
		if (u) {
			up_i = u[i];
		}
		
		// Extract the i-th coordinate of the point x
		var x_i = x[i];

		// Compare the i-th coordinate of the point x with the lower and the upper bounds
		// of the restricted simplex.
		if (x_i < lb_i || x_i > up_i) {
			return Number.POSITIVE_INFINITY;
		}
		
		// Compute the sum of the coordinates of the point x for later use
		sum_xi += x_i;
	}
	
	// Check if the point x belongs to the restricted unit full simplex, second step.
	//
	// Note: Due to limited numerical precision, the test on the sum of the coordinates
	// of the point x (which must be <= 1) cannot be exact.
	if (sum_xi - 1 > eps_tol) {
		return Number.POSITIVE_INFINITY;
	}
	

	// At this stage, the input point x belongs to the restricted unit full simplex
	return 0;
}


/**
* @function fullSimplexEmptinessCheck_
*
* @summary Checks the emptiness of the restricted unit full simplex.
*
* @description This function checks the emptiness of the restricted unit full simplex
* of R^n, which is defined as the unit full simplex of R^n subject to lower and upper 
* bounds constraints.
*
* In more details, this functions checks that the lower bounds l_i, i = 1..n 
* and the upper bounds u_i, i = 1..n satisfy:
* - 0 <= l_i <= u_i <= 1, i = 1..n
* - sum l_i <= 1
*
* Which is a necessary and sufficient condition for the restricted unit full simplex to not be empty.
*
* @param {number} n the dimension of the unit full simplex of R^n, natural integer superior or equal to 1.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @throws Throws an error in case the restricted unit full simplex is empty.
* @return {Array.<number>} In case the restricted unit full simplex is not empty, returns an array of 2 real numbers:
* - the sum of the lower bounds, sum l_i
* - the sum of the upper bounds, sum u_i
*
* @example
* fullSimplexEmptinessCheck_(3, [0.5, 0.1, 0.2]);
* //[0.8, 3]
*
* @example
* fullSimplexEmptinessCheck_(3, [0.7, 0.3, 0.2]);
* // new Error("infeasible problem detected: the restricted simplex is empty")
*/
function fullSimplexEmptinessCheck_(n, l, u) {
	// Initializations
	var sumLowerBounds = 0;
	var sumUpperBounds = 0;
	
	// - Feasibility checks on the restricted unit full simplex:
	// -- Lower bounds and upper bounds must satisfy 0 <= l_i <= u_i <= 1, i = 1..n
	// -- Lower bounds and upper bounds must satisfy sum l_i <= 1
	//
	// - Let's prove that these are necessary and sufficient conditions for the restricted unit full simplex to not be empty:
	// -- If the restricted unit full simplex is not empty, there exist a point x such that l_i <= x_i <= u_i <= 1, i = 1..n and
	//    such that sum x_i <= 1, per definition of the the restricted unit full simplex.
	//    This implies that sum l_i <= 1.
	//
	// -- If the restricted unit full simplex is such that 0 <= l_i <= u_i <= 1, i=1..n and sum l_i <= 1, then, the point
	//    x = (l_1,...,l_n) is such that sum x_i <= 1 and 0 <= l_i = x_i <= u_i <= 1, so that it is not empty
	for (var i = 0; i < n; ++i) {
		var lowerBound = 0;
		if (l) {
			lowerBound = l[i];
		}

		var upperBound = 1;
		if (u) {
			upperBound = u[i];
		}

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
	
	if (sumLowerBounds > 1) {
		throw new Error('infeasible problem detected: the restricted simplex is empty');
	}
	
	// Return the sum of lower bounds and the sum of upper bounds,
	// for potential usage in the calling function.
	return [sumLowerBounds, sumUpperBounds];
}

/**
* @function simplexEmptinessCheck_
*
* @summary Checks the emptiness of the restricted unit simplex.
*
* @description This function checks the emptiness of the restricted unit simplex
* of R^n, which is defined as the unit simplex of R^n subject to lower and upper 
* bounds constraints.
*
* In more details, this functions checks that the lower bounds l_i, i = 1..n 
* and the upper bounds u_i, i = 1..n satisfy:
* - 0 <= l_i <= u_i <= 1, i = 1..n
* - sum l_i <= 1 <= sum u_i, c.f. formula 2.2 of the reference
*
* @see <a href="https://doi.org/10.1016/S0167-7152(99)00095-4">Kai-Tai Fang, Zhen-Hai Yang, On uniform design of experiments with restricted
* mixtures and generation of uniform distribution on some domains, Statistics & Probability Letters, Volume 46, Issue 2, 
* 15 January 2000, Pages 113-120</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @param {Array.<number>} l the optional lower bounds contraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds contraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @throws Throws an error in case the restricted unit simplex is empty.
* @return {Array.<number>} In case the restricted unit simplex is not empty, returns an array of 2 real numbers:
* - the sum of the lower bounds, sum l_i
* - the sum of the upper bounds, sum u_i
*
* @example
* simplexEmptinessCheck_(3, [0.5, 0.1, 0.2]);
* //[0.8, 3]
*
* @example
* simplexEmptinessCheck_(3, [0.7, 0.3, 0.2]);
* // new Error("infeasible problem detected: the restricted simplex is empty")
*/
function simplexEmptinessCheck_(n, l, u) {
	// Initializations
	var sumLowerBounds = 0;
	var sumUpperBounds = 0;
	
	// - Feasibility checks on the restricted simplex:
	// -- Lower bounds and upper bounds must satisfy 0 <= l_i <= u_i <= 1, i = 1..n
	// -- Lower bounds and upper bounds must satisfy sum l_i <= 1 <= sum u_i, c.f. formula 2.2 of the reference
	for (var i = 0; i < n; ++i) {
		var lowerBound = 0;
		if (l) {
			lowerBound = l[i];
		}

		var upperBound = 1;
		if (u) {
			upperBound = u[i];
		}

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
	
	// Return the sum of lower bounds and the sum of upper bounds,
	// for potential usage in the calling function.
	return [sumLowerBounds, sumUpperBounds];
}


/**
* @function simplexSparseEuclidianProjection_
*
* @summary Returns a closest point on the unit simplex subject to a sparsity constraint.
*
* @description This function computes an at most k-sparse euclidean projection of a point
* x = (x_1,...,x_n) in R^n onto the unit simplex of R^n.
* 
* Optionally, lower bounds and upper bounds constraints can be added to the problem, in which case 
* the unit simplex becomes a restricted unit simplex, and the bounds constraints are then to be understood
* as applying only to the at most k coordinates selected to be part of the computed projection.
*
* In case there are no bounds constraints, the algorithm used is an O(n) implementation of
* the algorithm 1 of the first reference.
*
* In case there are bounds constraints, the problem is NP-complete in general, c.f. the third reference,
* so that the only way to guarantee an exact solution is through an exhaustive computation 
* over all the subsets of size <= k of the set of indexes {1,...,n}, searching for an optimal and feasible
* index set solving the projection problem, c.f. section 5.1 of the second reference.
*
* This is unfortunately only tractable for small n in general (n <= 20), or certain
* (n, k) pairs with k << n or k ~ n.
* 
* @see <a href="https://arxiv.org/abs/1206.1529">Anastasios Kyrillidis, Stephen Becker, Volkan Cevher and, Christoph Koch, Sparse projections onto the simplex, arXiv:1206.1529 [cs.LG]</a>
* @see <a href="https://link.springer.com/article/10.1007/s11425-016-9124-0">Fengmin Xu, Yuhong Dai, Zhihu Zhao and Zongben Xu, Efficient projected gradient methods for cardinality constrained optimization, Science China Mathematics volume 62, pages245–268(2019)</a>
* @see <a href="https://epubs.siam.org/doi/abs/10.1137/140978077">Oleg P. Burdakov, Christian Kanzow, and Alexandra Schwartz, Mathematical Programs with Cardinality Constraints: Reformulation by Complementarity-Type Conditions and a Regularization Method, SIAM J. Optim., 26(1), 397–425.</a>
*
* @param {Array.<number>} x, a point belonging to R^n, array of n real numbers.
* @param {number} k, a natural integer strictly greater than one corresponding to the maximum desired sparsity
* (i.e., non-zero elements) of the projected point.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; 
* defaults to an array made of zeros if u is provided, otherwise to null.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n;
defaults to an array made of ones if l is provided, otherwise to null.
*
* @return {Array.<number>|Float64Array} the computed closest point to x with at most k non-zero elements, array of n real numbers.
*
* @example
* simplexSparseEuclidianProjection_([0.5, 0.1, 0.2], 1);
* //[1,0,0]
*/
function simplexSparseEuclidianProjection_(x, k, l, u) {
	// Initializations
	var n = x.length;
	
	
	// Short circuit in case there is no sparsity
	if (k === n) {
		return simplexEuclidianProjection_(x, l, u);
	}
	
	
	// In case no bounds constraints are provided, use the algorithm 1 of the first reference.
	if (!l && !u) {
		// Compute the support of the projection, i.e., the k largest elements of x.
		
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
		
		
		// Compute the projection on the unit simplex of the k largest elements of x.
		var proj_x_k = simplexEuclidianProjection_(x_k);
		
		
		// Compute the final projection by re-conciliating the support of the
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
	// In case bounds constraints are provided, the computation is unfortunately 
	// NP-complete in general, c.f. the third reference.
	else {
		// Convert x to Matrix format to ease the computations below
		var x = new Matrix_(x);
		
		// Initialize the current minimum euclidean distance between x and the
		// best at most k-sparse projection of x, as well as the current list of 
		// associated indexes.
		var minDistanceValue = Infinity;
		var minDistanceIndexes = [];
		var minDistanceProjX = null;
		
		// An exhaustive enumeration of all the subsets of the set {1,...,n} 
		// of size between 1 and k is done, searching for the best feasible 
		// projection over all these subsets.
		for (var K = 1; K <= k; ++K) {
			var nextKSubsetIterator = new kSubsetsIterator_(n, K, false);
			var nextKSubset = nextKSubsetIterator.next();
			
			while (nextKSubset != -1) {
				// Extract the selected indexes of {1..n}
				var subsetNbIndexes = nextKSubset.length;
				var subsetIndexes = typeof UInt32Array === 'function' ? new UInt32Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				for (var i = 0; i < subsetNbIndexes; ++i) {
					subsetIndexes[i] = nextKSubset[i];
				}

				// Extract the coordinates of the vector x associated to the selected indexes
				var subsetX = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				for (var i = 0; i < subsetNbIndexes; ++i) {
					subsetX[i] = x.data[subsetIndexes[i]-1];
				}

				// Extract the lower and upper bounds constraints associated to the selected indexes
				var subsetL = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				if (l) {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetL[i] = l[subsetIndexes[i]-1];
					}
				}
				else {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetL[i] = 0;
					}
				}
				var subsetU = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				if (u) {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetU[i] = u[subsetIndexes[i]-1];
					}
				}
				else {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetU[i] = 1;
					}
				}
				
				// Compute the projection of the selected indexes of x on the restricted simplex.
				//
				// If the projection is better (in the sense of the euclidean distance to x) than the current
				// best projection, it becomes the new best projection and the current subset of indexes
				// becomes the new optimal subset of indexes.
				//
				// Note: because the restricted simplex associated to the subset of selected indexes
				// might be empty, special care must be taken.
				try {
					// Compute the projection
					var proj_subsetX = simplexEuclidianProjection_(subsetX, subsetL, subsetU);
					
					// Transform the projection into an at most k sparse vector of R^n
					var proj_x = Matrix_.zeros(n, 1);
					for (var i = 0; i < subsetNbIndexes; ++i) {
						proj_x.data[subsetIndexes[i] - 1] = proj_subsetX[i];
					}
					
					// Compute the euclidean distance between the initial x and the 
					// at most k sparse vector above.
					var d_x_proj_x = Matrix_.axpby(1, x, -1, proj_x).vectorNorm('two');
					
					// Determine if the projection above is better than the current best projection
					if (d_x_proj_x < minDistanceValue) {
						minDistanceValue = d_x_proj_x;
						minDistanceIndexes = subsetIndexes;
						minDistanceProjX = proj_x;
					}
					
				}
				catch (e) {
					if (e.message !== "infeasible problem detected: the restricted simplex is empty") {
						throw(e);
					}
				}

				// Generate a new subset	
				var nextKSubset = nextKSubsetIterator.next();
			}
		}
		
		// In case a best projection has been found, return it.
		if (minDistanceValue != Infinity) {
			return minDistanceProjX.toArray();
		}
		else {
			throw new Error('infeasible problem detected');
		}
	}
}


/**
* @function fullSimplexEuclidianProjection_
*
* @summary Returns the closest point on the unit full simplex.
*
* @description This function computes the euclidean projection of a point
* x = (x_1,...,x_n) in R^n onto the unit full simplex of R^n.
*
* Optionally, lower bounds and upper bounds constraints can be added to the problem, in which case 
* the unit simplex becomes a restricted unit full simplex.
*
* Internally, the algorithm used is an O(n) algorithm, c.f. the reference.
*
* @see Amir Beck, First-Order Methods in Optimization, MOS-SIAM Series on Optimization. SIAM. 2017
*
* @param {Array.<number>} x a point belonging to R^n, array of n real numbers.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @return {Array.<number>} the computed closest point to x, array of n real numbers.
*
* @example
* fullSimplexEuclidianProjection_([1, 1, 1]);
* // [~0.33, ~0.33, ~0.33]
*/
function fullSimplexEuclidianProjection_(x, l, u) {
	// Initializations
	var n = x.length;
	
	// Emptiness check on the restricted unit full simplex.
	//
	// In case the restricted simplex is empty, an exception is thrown, so that
	// the process is (violently) stopped here.
	fullSimplexEmptinessCheck_(n, l, u);


	// C.f. example 6.32 of the reference, the problem of the euclidean projection
	// on the restricted unit full simplex can be reduced to the problem of the
	// euclidean projection on the restricted unit simplex.
	
	// Compute the projection of the input point on the box 0 <= l_i <= u_i <= 1
	var p_box = new Matrix_(x).elemMap(function(i,j,val) { 
	                                      var l_i = l ? l[i-1] : 0;
										  var u_i = u ? u[i-1] : 1;
	                                      
										  return Math.max(-l_i, Math.min(u_i, val));
									   });
	
	// Compute the scalar product <p_box / (1,...,1)> to determine if a projection
	// on the restricted unit simplex is necessary.
	var a = Matrix_.ones(n, 1);
	var ps = Matrix_.vectorDotProduct(a, p_box);
	
	var p;
	if (ps <= 1) {
		p = p_box.toArray();
	}
	else {
		p = simplexEuclidianProjection_(x, l, u);
	}
	
	// Return the computed projection
	return p;
}


/**
* @function simplexEuclidianProjection_
*
* @summary Returns the closest point on the unit simplex.
*
* @description This function computes the euclidean projection of a point
* x = (x_1,...,x_n) in R^n onto the unit simplex of R^n.
*
* Optionally, lower bounds and upper bounds constraints can be added to the problem, in which case 
* the unit simplex becomes a restricted unit simplex.
*
* Internally, the algorithm used is an O(n) algorithm, c.f. the reference.
*
* @see <a href="https://link.springer.com/article/10.1007/s10107-006-0050-z">Kiwiel, K.C., Breakpoint searching algorithms 
* for the continuous quadratic knapsack problem, Math. Program. (2008) 112: 473</a>
*
* @param {Array.<number>} x a point belonging to R^n, array of n real numbers.
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
* @return {Array.<number>} the computed closest point to x, array of n real numbers.
*
* @example
* simplexEuclidianProjection_([1, 1, 1]);
* // [~0.33, ~0.33, ~0.33]
*/
function simplexEuclidianProjection_(x, l, u) {
	// Initializations
	var n = x.length;

	
	// Emptiness check on the restricted simplex.
	//
	// In case the restricted simplex is empty, an exception is thrown, so that
	// the process is (violently) stopped here.
	simplexEmptinessCheck_(n, l, u);


	// Convert the problem of the euclidean projection on the restricted unit simplex
	// into the associated instance of the continuous quadratic knapsack problem.	
	var d = Matrix_.ones(n, 1);
	var a = new Matrix_(x);
	var b = Matrix_.ones(n, 1);
	var r = 1;
	var lb = l ? new Matrix_(l) : Matrix_.zeros(n, 1);
	var ub = u ? new Matrix_(u) : Matrix_.ones(n, 1);

	
	// Solve this instance
	//
	// Due to numerical rounding issues when projecting big numbers, several successive projections are done
	// until the computed point belong to the simplex.
	var sol = qksolveBS_(d, a, b, r, lb, ub);
	var y = sol[0].toArray();
	while (simplexCharacteristicFunction_(y, l, u) == Number.POSITIVE_INFINITY) {
		a = new Matrix_(y);
		sol = qksolveBS_(d, a, b, r, lb, ub);
		y = sol[0].toArray();
	}
	
	// Return the computed projection
	return y;
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
* @param {boolean} reuseOutputArray an optional boolean that can be set to true to re-use the same output array throughout
* all the computations (this improves the performances, but requires the caller to NOT alter the output array); defaults to false.
* @return {function} a function to be used through its .sample() method, computing all 
* the points on the k-th rational grid of the unit simplex of R^n.
*
* @example
* var mySampler = new simplexGridSampler_(3, 10);
* mySampler.sample(); mySampler.sample(); ...; mySampler.sample();
* // [1, 0, 0]; [0.9, 0.1, 0]; ...; -1
*/
function simplexGridSampler_(n, k, reuseOutputArray) {
	// Initializations
	this.n = n;
	this.k = k;
	
	this.reuseOutputArray = reuseOutputArray;
	
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.compositionIterator = new compositionsIterator_(k, n, true); // reuse the ouput array for better performances
	
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
	* @return {Array.<number>|Float64Array|-1} an array of n real numbers corresponding to the coordinates of the generated point in R^n,
	* or -1 in case all such points have been generated.
	*/
	this.sample = function() {
		// Generate a new k-composition of n
		var comp = this.compositionIterator.next();

		// Return -1 in case there is no more samples to draw
		if (comp == -1) {
			return -1;
		}
		
		// Otherwise, compute the current rational grid point by normalizing the generated k-composition
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = comp[i] / this.k;
		}

		// Return either the point being sampled, or a copy of the point being sampled so that callers can alter it
		if (this.reuseOutputArray) {
			return this.x;
		}
		else {
			return this.x.slice(0);
		}
	}
}
 
 

/**
* @function simplexRandomSampler_
*
* @summary Returns a function to compute random points on the unit simplex of R^n.
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
* @param {Array.<number>} l the optional lower bounds constraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
* @param {Array.<number>} u the optional upper bounds constraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n.
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
	if (this.lowerBounds || this.upperBounds) {
		// Emptiness check on the restricted simplex.
		//
		// In case the restricted simplex is empty, an exception is thrown, so that
		// the process is (violently) stopped here.
		var sumBounds = simplexEmptinessCheck_(n, this.lowerBounds, this.upperBounds);
		this.sumLowerBounds = sumBounds[0];
		this.sumUpperBounds = sumBounds[1];
		
		// In case lower bounds or upper bounds are binding (sum l_i == 1 or sum u_i == 1), 
		// the computations are prematurely stopped, because the restricted simplex is then equal to a point.
		//
		// Otherwise:
		// - Deletion of possible superfluous constraints, as described in formula 2.3 of the second reference
		//
		// - Computation of the upper bounds*, as defined after formula 2.3' of the second reference
		//
		// - In parallel of the steps above, computation of upper* bounds sums, 
		// as well as the upper* bounds running sum
		if (this.sumLowerBounds == 1 || this.sumUpperBounds == 1) {
			// Nothing to do
		}
		else {
			// Deletion of possible superfluous constraints, replacing the lower and upper bounds
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
	* @description This function computes a point chosen uniformly at random on the unit simplex of R^n,
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
	* @description This function computes a point chosen uniformly at random on the unit simplex of R^n,
	* subject to an exact bounds constraints on its coordinates, which makes the point unique and non random.
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
	* @description This function computes a point chosen uniformly at random on the unit simplex of R^n,
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
			// of the random point being generated on T_n(0,upperStarBounds) must be set to zero.
			//
			// The main loop can then be prematurely stopped.
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
* @see Nicholas J. Higham. 2002. Accuracy and Stability of Numerical Algorithms (2nd ed.). Soc. for Industrial and Applied Math., Philadelphia, PA, USA. 
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
	* @description This function computes a point chosen uniformly at random in the intersection of 
	* the unit hypersphere of R^n and of the hyperplane defined by the equation <(1,1,...,1)/x> = 0,
	* using the O(n) algorithm 2 of the reference.
	*
	* @memberof simplexDirectionSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample = function() {
		// Computation of n independent random variables from N(0,1), which will form
		// the coordinates of the point x being sampled.
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from N(0,1)
			var r = normrnd_(0, 1);
			
			// Set the i-th coordinate of the point being sampled
			this.x[i] = r;
			
			// Compute the running sum of the normal variables
			sum += r;
		}
		
		// Normalization of the computed coordinates of the point being sampled, so that
		// the associated vector in R^n also belongs to the hyperplane <(1,1,...,1)/x> = 0,
		// i.e. sum x_i = 0.
		// 
		// In parallel, compute the 2 norm of the vector, for subsequent 
		// normalization, with an accurate algorithm by S. J. Hammarling,
		// c.f. problem 27.5 of the second reference.
		//
		// Note: The algorithm 2 of the reference stops here.
		var sum_d_n = sum / this.n;
		var t = 0;
		var s = 1;
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i] - sum_d_n;
			
			// Compute the running 2 norm of the associated vector.
			var absX = Math.abs(this.x[i]);
			if (absX != 0) {
				if (absX > t) {
					s = 1 + s * (t/this.x[i]) * (t/this.x[i]);
					t = absX;
				}
				else  {
					s = s + (this.x[i]/t) * (this.x[i]/t);
				}
			}
		}
		
		// Final normalization of the computed coordinates of the point being sampled, so that
		// the associated vector in R^n belongs to the n-hypersphere, i.e. its 2-norm 
		// is equal to 1.
		var x_two_norm = t * Math.sqrt(s);
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
* @description Given a point x = (x_1,...,x_n) on the unit simplex of R^n, this function computes a proximal point xr = (xr_1,...,xr_n) on 
* the r-th rational grid of the unit simplex of R^n, 1/r * I_n(r), with I_n(r) the set of N^n containing the points m = (m_1,...,m_n) 
* satisfying sum m_i = r with r a strictly positive natural integer, so that the computed proximal point xr is one of the closest points to x 
* on this grid with respect to any norm in a large class, c.f. the first reference.
*
* @see <a href="https://doi.org/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the unit simplex: 
* Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258</a>
* @see <a href="https://arxiv.org/abs/1501.00014">Rama Cont, Massoud Heidari, Optimal rounding under integer constraints</a>
* 
* @param {Array.<number>} x a point belonging to the unit simplex of R^n, array of n real numbers.
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
* Optionally, lower bounds and upper bounds constraints can be added to the problem, in which case 
* the unit simplex becomes a restricted unit simplex.
*
* To be noted that per lemma 1 of the reference, the number of points on the k-th rational grid
*  of the unit simplex of R^n is equal to factorial(n + k - 1) / (factorial(k - 1) * factorial(n)), 
* i.e., binomial(n+k-1, n-1), so that this method might be of limited use, even for small n.
*
* For instance, n=5 and k=100 already result in 4598126 points on which to evaluate f.
*
* To also be noted that as the optional lower and upper bounds contraints become tigter,
* the volume of the restricted simplex becomes smaller, so that this method might also be of limited use
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
* @param {Array.<number>} l the optional lower bounds contraints, an array of n real numbers l_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of zeros.
* @param {Array.<number>} u the optional upper bounds contraints, an array of n real numbers u_i which must satisfy 0 <= l_i <= u_i <= 1, i = 1..n; defaults to an array made of ones.
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
	
	// Emptiness check on the restricted simplex.
	//
	// In case the restricted simplex is empty, an exception is thrown, so that
	// the process is (violently) stopped here.
	if (lowerBounds || upperBounds) {
		var sumBounds = simplexEmptinessCheck_(n, lowerBounds, upperBounds);
	}
	
	// Initialize the current minimum value and the current list of associated grid points
	var minValue = Number.POSITIVE_INFINITY;
	var minValueGridPoints = [];

	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference.
	var sampler = new simplexGridSampler_(n, k, true); // use no array copy in the simplex grid sampler to improve performances
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
