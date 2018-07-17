/**
 * @file Misc. computational geometry functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.geometricCenter_ = geometricCenter_;
self.geometricMedian_ = geometricMedian_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function geometricCenter_
*
* @summary Compute the geometric center of a finite set of points belonging to R^n.
*
* @description This function returns the geometric center of m points x_1,...x_m 
* belonging to R^n, which is defined as the component-wise arithmetic mean of the m points.
*
* The geometric center of the m points x_1, ..., x_m is also the point y which 
* minimizes the sum of the squared Euclidean distances between itself and each point:
*
* y = argmin_x in R^n f(x) = sum ||y - x_i||_2^2, i = 1..m
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error
* in the computation of the component wise mean, c.f. the second reference.
*
* @see <a href="https://en.wikipedia.org/wiki/Centroid">Centroid</a>
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
*
* @param {Array.<Matrix_>} x an array of m n by 1 matrices, corresponding to the coordinates 
* of the m points belonging to R^n.
* @return {Matrix_} the geometric center of the m points x_1,...x_m
*
* @example
* geometricCenter_([new Matrix([0,1,2]), new Matrix([1,2,3])]);
* // new Matrix([0.5,1.5,2.5]) 
*/
function geometricCenter_(x) {
	// TODO: Checks
	
	// Initialisations
	var m = x.length;
	var n = x[0].nbRows;

	// Instanciate the geometric center
	var y = Matrix_.zeros(n, 1);
	
	// For each coordinate i of the input points:
	// - Compute the mean over the m points of the coordinate i (first pass)
	// - Compute the correction factor (second pass), c.f. M_3 formula of the 
	// second reference
	// - Set the geometric center coordinate i to the corrected mean over the m points
	// of the coordinate i
	for (var i = 1; i <= n; ++i) {
		// Mean computation
		var sum_i = 0.0;
		for (var k = 0; k < m; ++k) {
			sum_i += x[k].getValue(i, 1);
		}
		var tmpMean_i = sum_i/m;

		// Correction factor computation
		var sumDiff_i = 0.0;
		for (var k = 0; k < m; ++k) {
			sumDiff_i += (x[k].getValue(i, 1) - tmpMean_i);
		}

		// Corrected mean computation
		y.setValue(i, 1,
		           (sum_i + sumDiff_i)/m);
	}
	
	// Return the computed geometric center
	return y;
}

/**
* @function geometricMedian_
*
* @summary Compute the geometric median of a finite set of points belonging to R^n.
*
* @description This function returns the geometric median of m points x_1,...x_m 
* belonging to R^n, which is defined as the point y which minimizes 
* the sum of the Euclidean distances between itself and each point:
*
* y = argmin_x in R^n f(x) = sum ||y - x_i||_2, i = 1..m
*
* The algorithm implemented uses a serie of successive hyperbolic approximations of
* the euclidian norms appearing in the function f above, c.f. the second reference,
* which allows to compute the geometric median using a standard first-order convex
* optimization algorithm.
*
* @see <a href="https://en.wikipedia.org/wiki/Geometric_median">Geometric median</a>
* @see <a href="http://dx.doi.org/10.1287/opre.23.3.581">Robert F. Love, James G. Morris, (1975) Technical Note—Solving Constrained Multi-Facility Location Problems Involving lp Distances Using Convex Programming. Operations Research 23(3):581-587.</a>
* @see <a href="https://arxiv.org/abs/1606.05225">Michael B. Cohen, Yin Tat Lee, Gary Miller, Jakub Pachocki, Aaron Sidford. Geometric Median in Nearly Linear Time. arXiv:1606.05225 [cs.DS]</a>
*
* @param {Array.<Matrix_>} x an array of m n by 1 matrices, corresponding to the coordinates 
* of the m points belonging to R^n.
* @return {Matrix_} the geometric median of the m points x_1,...x_m
*
* @example
* geometricMedian_([new Matrix([0,1,2]), new Matrix([1,2,3])]);
* // new Matrix([0.5, 1.5, 2.5]) 
*/
function geometricMedian_(x) {
    // Internal function to compute the function C_ph, 
	// approximation of the function f, c.f. formula 1 
	// of the second reference.
	function f_eps(y) {
		var sum = 0.0;

		for (var k = 0; k < m; ++k) {
			// Compute Math.SQRT(||y - x_k||_2^2 + eps), using an inlined version
			// of the vectorNorm('two') Matrix function and a stable way to
			// compute the square root of two numbers squared.
			//var y_m_x_k = Matrix_.xmy(y, x[k], tmp_vec_n);
			//var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');
			var t = 0;
			var s = 1;
			for (var i = 1; i <= n; ++i) {
				var val = y.getValue(i, 1) - x[k].getValue(i, 1); // y_i - (x_k)_i
				var absVal = Math.abs(val);
				if (absVal != 0) {
					if (absVal > t) {
						s = 1 + s * (t/val) * (t/val);
						t = absVal;
					}
					else  {
						s = s + (val/t) * (val/t);
					}
				}
			}
			var y_m_x_k_two_norm = t * Math.sqrt(s);
			
			sum += hypot_(y_m_x_k_two_norm, eps_f);
		}
		return sum;
	}
	
    // Internal function to compute the function grad(C_ph),
	// approximation of the gradient of the function f, c.f.
	// formula 1 of the second reference.
	function gradf_eps(y) {
		var res = Matrix_.zeros(n, 1);
		
		for (var k = 0; k < m; ++k) {
			// Compute (y - x_k)/||y - x_k||_2 and add it
			// to the currently computed gradient, using an inlined version
			// of the vectorNorm('two') Matrix function and a stable way to
			// compute the square root of two numbers squared.
			//var y_m_x_k = Matrix_.xmy(y, x[k-1], tmp_vec_n);
			//var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');
			var t = 0;
			var s = 1;
			for (var i = 1; i <= n; ++i) {
				var val = y.getValue(i, 1) - x[k].getValue(i, 1);  // y_i - (x_k)_i
				tmp_vec_n.setValue(i, 1, 
				                   val);
				
				var absVal = Math.abs(val);
				if (absVal != 0) {
					if (absVal > t) {
						s = 1 + s * (t/val) * (t/val);
						t = absVal;
					}
					else  {
						s = s + (val/t) * (val/t);
					}
				}
			}
			var y_m_x_k = tmp_vec_n;
			var y_m_x_k_two_norm = t * Math.sqrt(s);

			res = Matrix_.axpby(1, res, 1/hypot_(y_m_x_k_two_norm, eps_f), y_m_x_k, res);
		}
		return res;
	}
	
	
	// TODO: Checks
	
	
	// Initialisations
	var m = x.length; // the number of points provided in input
	var n = x[0].nbRows; // the dimension of each point provided in input

	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	
	
	// The geometric median is computed using successive hyperbolic approximations of
	// the euclidian norms appearing in its objective function, c.f. the second
	// reference.
	//
	// Each hyperbolic approximation of the objective function of the geometric median
	// problem is a smooth convex(/strictly convex) function, so that the associated
	// minimization problem can be solved using a standard first-order convex optimization
	// algorithm (here, FISTA-like).
	
	
	// Compute a proper starting point for the optimization algorithm.
	//
	// Per lemma 18 of the third reference, the geometric center is a 
	// 2-approximation of the geometric median.
	var x0 = geometricCenter_(x);
	
	
	// Define additional functions used by the optimization algorithm
	
	// The projection on R^n
	var g = function(x) {
		return 0;
	}
		
	// The proximal function associated to g is the orthogonal
	// projection on R^n, i.e., the identity.
	var proxg = function(x, mu) {
		return x;
	}
		
			
	// Compute the minimum of the function f_eps on R^n, for successive decreasing values
	// of epsilon (which is actually squared in the computation of f_eps and gradf_eps).
	//
	// Precision in the early stages is not of paramount importance.
	//
	// Note: the associated loop has been unrolled.
	var eps_f = 1e-3;
	var sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, x0, {eps: 1e-2, maxIter: -1, maxLine: -1});

	eps_f = 1e-4;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-2, maxIter: -1, maxLine: -1});
	
	eps_f = 1e-5;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-3, maxIter: -1, maxLine: -1});
	
	eps_f = 1e-6;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-3, maxIter: -1, maxLine: -1});

	eps_f = 1e-7;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-4, maxIter: -1, maxLine: -1});

	eps_f = 1e-8;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-4, maxIter: -1, maxLine: -1});


	// Return the computed optimal solution to the last hyperbolic approximation of the 
	// geometric median.
	return sol[0];
}
