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
* The algorithm implemented uses a serie of successive converging hyperbolic approximations of
* the euclidian norms appearing in the function f above, c.f. the second and third references,
* which allows to compute the geometric median using a standard first-order convex
* optimization algorithm.
*
* @see <a href="https://en.wikipedia.org/wiki/Geometric_median">Geometric median</a>
* @see <a href="http://dx.doi.org/10.1287/opre.23.3.581">Robert F. Love, James G. Morris, (1975) Technical Note—Solving Constrained Multi-Facility Location Problems Involving lp Distances Using Convex Programming. Operations Research 23(3):581-587.</a>
* @see <a href="https://dx.doi.org/10.4169%2Famer.math.monthly.121.02.095%23sthash.QTTb5Z6T.dpuf">Eric C. Chi and Kenneth Lange, A Look at the Generalized Heron Problem through the Lens of Majorization-Minimization. Am Math Mon. 2014 Feb; 121(2): 95–108.</a>
* @see <a href="https://arxiv.org/abs/1606.05225">Michael B. Cohen, Yin Tat Lee, Gary Miller, Jakub Pachocki, Aaron Sidford. Geometric Median in Nearly Linear Time. arXiv:1606.05225 [cs.DS]</a>
* @see <a href="https://doi.org/10.1007/BFb0083582">Ben-Tal A., Teboulle M. (1989) A smoothing technique for nondifferentiable optimization problems. In: Dolecki S. (eds) Optimization. Lecture Notes in Mathematics, vol 1405. Springer, Berlin, Heidelberg</a>
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
			// Compute Math.SQRT(||y - x_k||_2^2 + eps^2), using a stable way to
			// compute the square root of two numbers squared.
			var y_m_x_k = Matrix_.xmy(y, x[k], tmp_vec_n);
			var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');
			
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
			// Compute (y - x_k)/Math.SQRT(||y - x_k||_2^2 + eps^2) and add it
			// to the currently computed gradient, using a stable way to
			// compute the square root of two numbers squared.
			var y_m_x_k = Matrix_.xmy(y, x[k], tmp_vec_n);
			var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');

			res = Matrix_.axpby(1, res, 1/hypot_(y_m_x_k_two_norm, eps_f), y_m_x_k, res);
		}
		
		return res;
	}

	
	// TODO: Checks
	
	
	// Initialisations
	var m = x.length; // the number of points provided in input
	var n = x[0].nbRows; // the dimension of each point provided in input

	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	
	
	// The geometric median will be computed using successive epsilon-approximations of
	// the euclidian norms appearing in its objective function, c.f. the second and the third
	// references.
	//
	// - Each epsilon-approximation of the objective function of the geometric median
	// problem is a smooth convex(/strictly convex) function, so that the associated
	// minimization problem can be solved using a standard first-order convex optimization
	// algorithm (below, FISTA-like).
	//
	// - As the epsilon parameter defining these epsilon-approximations converges to zero, 
	// the solution found by the convex optimization algorithm is proven to converge 
	// to the true geometric median, c.f. the third reference.
	
	// Compute a proper starting point for the optimization algorithm.
	//
	// Per lemma 18 of the fourth reference, the geometric center is a 
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
		
			
	// TODO: BETTER COMMENT + add 1e-2 as parametrized in input
	// Compute the minimum of the function f_eps on R^n, for successive decreasing values
	// of epsilon (which is actually squared in the computation of f_eps and gradf_eps).
	//
	// 0 <= G(x^*_eps) - G(x^*) <= eps*m  =>  G(x^*_eps) - eps*m <= G(x^*) <= G(x^*_eps)
	// c.f. example 3.3 of the fifth reference 
	var eps_f = 1e-2 / m;
	var sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, x0, {eps: 1e-4, maxIter: -1, maxLine: -1});


	// Return the computed optimal solution to the hyperbolic approximation of the 
	// geometric median.
	return sol[0];
}
