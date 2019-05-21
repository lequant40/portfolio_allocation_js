/**
 * @file Misc. optimisation functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.lpsolvePDHG_ = lpsolvePDHG_;
self.qpsolveGSMO_ = qpsolveGSMO_;
self.qksolveBS_ = qksolveBS_;
self.ccpsolveFISTA_ = ccpsolveFISTA_;
/* End Wrapper private methods - Unit tests usage only */
 

/**
* @function ccpsolveFISTA_
*
* @summary Returns an optimal solution to a compositive convex problem, 
* using a FISTA-like accelerated first-order algorithm.
*
* @description This function computes an optimal solution to a compositive convex
* problem using a FISTA-like accelerated first-order algorithm, c.f. the first reference.
*
* The compositive convex problem to solve is assumed to be provided in the
* following format:
*
* min F(x) = f(x) + g(x), x in R^n
*
* f : R^n -> R is a continuously differentiable convex function with a Lipschitz continuous gradient
* g : R^n -> R u {+oo} is a (proximable) proper closed convex function
* gradf : R^n -> R^n is the Lipschitz continuous gradient of f
* proxg : R^n x R^+* -> R^n is the proximal operator associated to g defined as 
* proxg(x, mu) = argmin u in R^n ( g(u) + 1/(2*mu) * ||u - x||_2^2 )
*
* The problem is assumed to be solvable, i.e., argmin F(x), x in R^n, is
* assumed to be non-empty.
*
* The algorithm used internally is based on the FISTA-BKTR algorithm of the third 
* reference, which is an optimal first-order method for a smooth problem (i.e., 
* it ensures a convergence rate of O(1/k^2)), with the following additions:
* - The usage of a convergence criterion based on the gradient of f and on a subdifferential of g,
* c.f. the fourth reference
* - The usage of both a fixed and of an adaptative restart mechanism, c.f. the fifth reference
* - The usage of a Barzilai and Borwein like step size, c.f. the sixth reference
*
* @see <a href="https://doi.org/10.1137/080716542">Amir Beck and Marc Teboulle, A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems, SIAM Journal on Imaging Sciences 2009 2:1, 183-202</a>
* @see <a href="https://doi.org/10.1109/TIP.2009.2028250">A. Beck, M. Teboulle, "Fast gradient-based algorithms for constrained total variation image denoising and deblurring problems", IEEE Trans. Image Process., vol. 18, no. 11, pp. 2419-2434, 2009</a>
* @see <a href="https://doi.org/10.1007/s10208-014-9189-9">Scheinberg, K., Goldfarb, D. & Bai, X. Fast First-Order Methods for Composite Convex Optimization with Backtracking Found Comput Math (2014) 14: 389.</a>
* @see <a href="https://arxiv.org/abs/1411.3406">T. Goldstein, C. Studer, and R. G. Baraniuk, “A field guide to forward-backward splitting with a FASTA implementation,” Nov. 2014</a>
* @see <a href="https://doi.org/10.1137/16M1055323">Bo Wen, Xiaojun Chen, and Ting Kei Pong. Linear Convergence of Proximal Gradient Algorithm with Extrapolation for a Class of Nonconvex Nonsmooth Minimization Problems. SIAM Journal on Optimization 2017 27:1, 124-145</a>
* @see <a href="https://doi.org/10.1007/s10589-006-6446-0">Gradient Methods with Adaptive Step-Sizes. Zhou, B., Gao, L. & Dai, YH. Comput Optim Applic (2006) 35: 69.</a>
*
* @param {function} f, a function representing the function f above, which must take as input argument
* a n by 1 matrix x corresponding to a point in R^n and which must return as output a real number 
* corresponding to f(x).
* @param {function} gradf, a function representing the gradient of the function f above, 
* which must take as input argument a n by 1 matrix x corresponding to a point in R^n and 
* which must return as output a n by 1 matrix gradf(x) corresponding to gradf(x).
* @param {function} g, a function representing the function g above, which must take as input argument
* a n by 1 matrix x corresponding to a point in R^n and which must return as output a real number 
* or Number.POSITIVE_INFINITY corresponding to g(x).
* @param {function} proxg, a function representing the proximal operator associated to 
* the function g above, which must take as input arguments a n by 1 matrix x corresponding to 
* a point in R^n and a strictly positive real number mu corresponding to a step size and which 
* must return as output a n by 1 matrix corresponding to proxg(x, mu).
* @param {Matrix_} x0 an n by 1 matrix corresponding to the point on which to
* start the FISTA algorithm (usually, the best possible guess of the optimal solution).
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the absolute tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-04.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @param {number} opt.maxLine the maximum number of line searches in one iteration of the algorithm, a strictly positive natural integer or -1 to force an infinite number of line searches; defaults to 100.
* @param {number} opt.beta the step size multiplicative shrinkage factor used in the backtracking procedure, a real number belonging to ]0,1[; defaults to 0.5.
* @param {number} opt.alphaMin the minimum value of the step size, a strictly positive real number; defaults to 1e-10.
* @param {number} opt.alphaMax the maximum value of the step size, a strictly positive real number; defaults to 1e10.
* @param {number} opt.restartPeriod the restart period, expressed in a number of iterations, of the fixed restart mechanism of the algorithm; defaults to 1000 iterations.
* @return {Array<Object>} an array arr containing two elements: 
* - arr[0] an n by 1 matrix containing the optimal solution x^* to the compositive convex problem
* - arr[1] the optimal value of the function f, i.e. F(x^*)
*
* @example
* ccpsolveFISTA_(function(x) { return Math.exp((x.getValue(1, 1) - 0.7)*(x.getValue(1, 1) - 0.7)); }, // f(x) = exp((x - 0.7)^2)
*                function(x) { return new Matrix_([2 * (x.getValue(1, 1) - 0.7) * Math.exp((x.getValue(1, 1) - 0.7)*(x.getValue(1, 1) - 0.7))]); },  // gradf(x) = 2*(x - 0.7)*exp((x - 0.7)^2)
*				 function(x) { if (0 > x.getValue(1, 1) || x.getValue(1, 1) > 1) {
*				                   return Number.POSITIVE_INFINITY;
*			                   }
*			                   else {
*                                  return 0;
*	                           }
*			                 }, // g is the usual indicator function of a convex set, here [0,1]
*                function(x, mu) { return new Matrix_([Math.max(0, Math.min(x.getValue(1, 1), 1))]); }, // proxg(x, mu) = orthogonal projection of x on [0,1]
*                new Matrix_([0]) // the starting point of the algorithm
*               )
* // new Matrix_([~0.7])
*/
function ccpsolveFISTA_(f, gradf, g, proxg, x0, opt) {
	// Internal function to compute F(x) = f(x) + g(x), 
	// c.f. formula 1.1 of the third reference.
	function F(x) {
		return f(x) + g(x);
	}

	// Internal function to compute Q_mu(u,v) = f(v) + <u - v/gradf(v)> + 1/(2*mu) * ||u - v||_2^2 + g(u), 
	// c.f. formula 2.2 of the third reference.
	function Q(mu, u, v, gradf_v) {
		// Compute f(v)
		var f_v = f(v);

		// Compute u - v and ||u - v||_2
		var u_m_v = Matrix_.xmy(u, v);
		var u_m_v_two_norm = u_m_v.vectorNorm('two');
		
		// Compute g(u)
		var g_u = g(u);

		// Compute Q_mu
		var Q_mu_u_v = f_v + Matrix_.vectorDotProduct(u_m_v, gradf_v) + 1/(2 * mu) * u_m_v_two_norm * u_m_v_two_norm + g_u;
		
		// Return the computed value
		return Q_mu_u_v;
	}

	// Internal function to compute p_mu(v) = argmin_u Q_mu(u,v), 
	// c.f. formula 2.3 of the third reference.
	//
	// This function is shown to be equal to proxg(v - mu*gradf(v), mu)
	// in formula 3.13 of the second reference.
	function p(mu, v, gradf_v) {
		// Compute v - mu*gradf(v)
		var v_m_mu_gradf_v = Matrix_.axpby(1, v, -mu, gradf_v);
		
		// Compute p_mu
		var p_mu_v = Matrix_.copy(proxg(v_m_mu_gradf_v, mu));
		
		// Return both values
		return [v_m_mu_gradf_v, p_mu_v];
	}
	
	
    // ------
    
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-04;
	var maxIterations = opt.maxIter || 10000;
	var maxLineSearches = opt.maxLine || 100;
	var beta = opt.beta || 0.5;
	var alphaMin = opt.alphaMin || 1e-10;
	var alphaMax = opt.alphaMax || 1e10;
	var restartPeriod = opt.restartPeriod || 1000;
	
	
	// ------
	
	
	// Misc; initializations
	var n = x0.nbRows;
	var eps_tol = 1e-12; // used to numerically determine some conditions (backtrack, adaptative restart, stepsize)
	
	// Initializations, c.f. line 0 of the Algorithm 2 of the third reference			
	// Prediction parameter
	var t_km; 
	var t_k;

	// Theta parameter
	var theta_km;
	var theta_k;

	// Step size parameters
	var tau_k = 0.5;
	var mu_k_0 = alphaMin;
	
	// x iterates
	var x_k; // placeholder for the x_k vector
	var x_km; // placeholder for the x_k-1 vector
	var x_kmm; // placeholder for the x_k-2 vector
	var x_km_m_x_kmm; //  placeholder for the x_k-1 - x_k-2 vector
	var gradf_x_k = Matrix_.zeros(n, 1); // placeholder for the gradf(x_k) vector
	var gradf_x_km = Matrix_.zeros(n, 1);; // placeholder for the gradf(x_k-1) vector
	var gradf_x_kmm = Matrix_.zeros(n, 1);; // placeholder for the gradf(x_k-2) vector
	var gradf_x_km_m_gradf_x_kmm; // // placeholder for the gradf(x_k-1) - gradf(x_k-2) vector

	// y iterates
	var y_k = Matrix_.zeros(n, 1); // placeholder for the y_k vector
	var gradf_y_k = Matrix_.zeros(n, 1); // placeholder for the gradf(y_k) vector	
	
	
	// Main loop of the Algorithm 2 of the third reference,
	// guaranteed to converge per theorem 3.2 of the third reference.
	var restart = true; // a first initialization is needed at the first iteration 

	var iter = 0;	
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		
		
		// Update the number of iterations
		++iter;
		
		
		// (-) Check the condition for a fixed restart of the algorithm, 
		// c.f. section 3.3 of the fifth reference.
		if (iter % restartPeriod === 0) {
			x0 = x_k;
			restart = true;
		}
		
		
		// (-) Restart of the algorithm as needed
		if (restart === true) {
			// Initialization of the prediction parameter
			t_km = 0; 
			t_k = 1;
			
			// Initialization of the theta parameter
			theta_km = 1;
			theta_k = null;
			
			// Initialization of the x iterates
			x_k = Matrix_.copy(x0);
			x_km = Matrix_.copy(x_k);
			x_kmm = Matrix_.copy(x_km);
			x_km_m_x_kmm = Matrix_.zeros(n, 1);
			gradf_x_k = Matrix_.copy(gradf(x_k), gradf_x_k);
			gradf_x_km = Matrix_.copy(gradf_x_k, gradf_x_km);
			gradf_x_kmm = Matrix_.copy(gradf_x_km, gradf_x_kmm);
			gradf_x_km_m_gradf_x_kmm = Matrix_.zeros(n, 1);

			// Initialization of the y iterates
			y_k = Matrix_.copy(x_k, y_k); 
			gradf_y_k = Matrix_.copy(gradf_x_k, gradf_y_k);
			
			// No update of the step size parameters, as the step size can take any value
			// per algorithm 2 of the third reference.
			
			// The restart is completed
			restart = false;
		}
		
		
		// (1) of the Algorithm 2 of the third reference
		// - Initialization of the initial stepsize for the current iteration
		var mu_k = mu_k_0;
		
		
		// (2) of the Algorithm 2 of the third reference
		// - Optimized backtracking line search
		var p_mu = p(mu_k, y_k, gradf_y_k);
		var y_k_m_mu_k_gradf_y_k = p_mu[0];
		var p_mu_k_y_k = p_mu[1];
		
		var iter_ls = 0;
		while ( F(p_mu_k_y_k) > Q(mu_k, p_mu_k_y_k, y_k, gradf_y_k) + eps_tol ) {
			// Check the number of iterations
			if (maxLineSearches !== -1 && iter_ls > maxLineSearches) {
				throw new Error('maximum number of line searches reached: ' + maxLineSearches + ' at iteration: ' + iter);
			}
			
			// Update the number of line search iterations
			++iter_ls;
			
			// Reduction in step size, and associated update of the theta_k parameter
			mu_k = beta * mu_k;
			theta_km = theta_km/beta;
			
			// Update of the current t_k and y_k iterates, due to the change in the theta_km
			// parameter.
			//
			// This is FistaStep(xk−1, xk−2, tk−1, θk−1)
			t_k = ( 1 + Math.sqrt(1 + 4*theta_km*t_km*t_km) ) / 2;
			y_k = Matrix_.axpby(1, x_km, (t_km - 1)/t_k, x_km_m_x_kmm, y_k);
			
			// Recomputation of gradf(y_k) and p_mu_k(y_k) for the next iteration
			//
			// The naive way to recompute gradf(y_k), i.e., computing the
			// gradient of f at point y_k can be improved by noticing
			// that the line search procedure do not update x_km and x_kmm,
			// c.f. the discussion after formula 3.17 of the third reference.
			gradf_y_k = Matrix_.axpby(1, gradf_x_km, (t_km - 1)/t_k, gradf_x_km_m_gradf_x_kmm, gradf_y_k);
			p_mu = p(mu_k, y_k, gradf_y_k);
			y_k_m_mu_k_gradf_y_k = p_mu[0];
			p_mu_k_y_k = p_mu[1];
		}
		
		
		// (3) of the Algorithm 2 of the third reference
		// - Computation of the x_k and t_k+1, y_k+1 iterates
		// - Computation of the initial stepsize for the next iteration
		// - Update of the theta_k parameter
		
		// Computation of the current x_k iterate
		x_k = p_mu_k_y_k;
		gradf_x_k = Matrix_.copy(gradf(x_k), gradf_x_k); // gradf(x_k)
		
		var x_k_m_x_km = Matrix_.xmy(x_k, x_km); // x_k - x_k-1
		var gradf_x_k_m_gradf_x_km = Matrix_.xmy(gradf_x_k, gradf_x_km); // gradf(x_k) - gradf(x_k-1)
		
		// Computation of the initial stepsize for the next iteration,
		// using a Barzilai and Borwein like stepsize, c.f. algorithm SS of the sixth reference.
		var s_k_d_s_k = Matrix_.vectorDotProduct(x_k_m_x_km, x_k_m_x_km); // <x_k - x_k-1/x_k - x_k-1>
		var z_k_d_z_k = Math.max(Matrix_.vectorDotProduct(gradf_x_k_m_gradf_x_km, gradf_x_k_m_gradf_x_km), // <gradf(x_k) - gradf(x_k-1)/gradf(x_k) - gradf(x_k-1)>
		                         eps_tol); // to avoid numerical issues in the division below
		var s_k_d_z_k = Matrix_.vectorDotProduct(x_k_m_x_km, gradf_x_k_m_gradf_x_km); // <x_k - x_k-1/gradf(x_k) - gradf(x_k-1)>

		if (s_k_d_z_k <= eps_tol) {
			mu_kp_0 = alphaMax;
		}
		else {
			var alpha_k_1 = Math.max(alphaMin, Math.min(s_k_d_s_k/s_k_d_z_k, alphaMax));
			var alpha_k_2 = Math.max(alphaMin, Math.min(s_k_d_z_k/z_k_d_z_k, alphaMax));
			
			if (alpha_k_2/alpha_k_1 <= tau_k) {
				mu_kp_0 = alpha_k_2;
				tau_k = 0.9 * tau_k;
			}
			else {
				mu_kp_0 = alpha_k_1;
				tau_k = 1.1 * tau_k;
			}
		}
		
		// Update of the theta_k parameter
		theta_k = mu_k/mu_kp_0;
		
		// Computation of the current t_k+1 and y_k+1 iterates,
		// plus misc. associated vectors.
		//
		// This is FistaStep(xk , xk−1, tk, θk)
		var t_kp = ( 1 + Math.sqrt(1 + 4*theta_k*t_k*t_k) ) / 2; // t_k+1
		var y_kp = Matrix_.axpby(1, x_k, (t_k - 1)/t_kp, x_k_m_x_km); // y_k+1	
		
		var gradf_y_kp = Matrix_.axpby(1, gradf_x_k, (t_k - 1)/t_kp, gradf_x_k_m_gradf_x_km); // gradf(y_k+1)


		// (-) Check the absolute convergence criteria (not in the third reference), 
		// c.f. paragraph 4.6 of the fourth reference
		var subgradg_x_k = Matrix_.axpby(1/mu_k, y_k_m_mu_k_gradf_y_k, -1/mu_k, x_k);
		var r_k = Matrix_.xpy(gradf_x_k, subgradg_x_k);
		if (r_k.vectorNorm('infinity') <= eps) {
			break;
		}
		
		
		// (-) Check the condition for an adaptative restart of the algorithm, 
		// c.f. section 3.3 of the fifth reference.
		var gs_k = Matrix_.vectorDotProduct(Matrix_.xmy(y_k, x_k), x_k_m_x_km);
		if (gs_k >= eps_tol) {
			x0 = x_k;
			restart = true;
		}

		
		// (-) Preparation of the next iteration:
		// Update of the step size
		mu_k_0 = mu_kp_0;
		
		// Update of the x_k iterates
		x_kmm = x_km;
		x_km = x_k;
		x_km_m_x_kmm = x_k_m_x_km;
		gradf_x_kmm = Matrix_.copy(gradf_x_km, gradf_x_kmm);
		gradf_x_km = Matrix_.copy(gradf_x_k, gradf_x_km);
		gradf_x_km_m_gradf_x_kmm = gradf_x_k_m_gradf_x_km;
		
		// Update of the y_k iterates
		y_k = Matrix_.copy(y_kp, y_k);
		gradf_y_k = Matrix_.copy(gradf_y_kp, gradf_y_k);
		
		// Update of the thera parameter
		theta_km = theta_k;
		
		// Update of the prediction parameters
		t_km = t_k;
		t_k = t_kp;
	}
	
	// Return the computed x_k value, as well as F(x_k)
	return [x_k, F(x_k)];
}


/**
* @function qksolveBS_
*
* @summary Returns an optimal solution to the continuous quadratic knapsack problem, 
* using a breakpoint searching algorithm.
*
* @description This function computes an optimal solution to the continuous quadratic
* knapsack problem using an O(n) breakpoint searching algorithm, c.f. the first reference.
*
* The problem to solve is assumed to be provided in the following format:
*
* min f(x) = 1/2 * <d*x/x> - <a/x>
*
* s.t. <b/x> = r (single linear equality constraint)
*      l <= x <= u (finite bound constraints)
*
* with:
* - d an n by 1 matrix with strictly positive elements, representing a diagonal n by n matrix with strictly positive elements
* - a an n by 1 matrix
* - r a real number
* - b an n by 1 matrix with strictly positive elements
* - l an n by 1 matrix
* - u an n by 1 matrix
* 
* To be noted that the algorithm used internally is able to detect the non-feasibility of the problem
* thanks to modifications based on the second reference, in which case an error is returned.
* 
* @see <a href="https://link.springer.com/article/10.1007/s10107-006-0050-z">Kiwiel, K.C., Breakpoint searching algorithms for the continuous quadratic knapsack problem, Math. Program. (2008) 112: 473</a>
* @see <a href="https://www.sciencedirect.com/science/article/pii/0167637784900105">Peter Brucker, An O(n) algorithm for quadratic knapsack problems, Operations Research Letters, Volume 3, Issue 3, 1984, Pages 163-166</a>
*
* @param {Matrix_} d an n by 1 matrix with strictly positive elements.
* @param {Matrix_} a an n by 1 matrix.
* @param {Matrix_} b an n by 1 matrix with strictly positive elements.
* @param {number} r a real number.
* @param {Matrix_} l an n by 1 matrix.
* @param {Matrix_} u an n by 1 matrix.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance when assessing the numerical equality <b/x> = r , a strictly positive real number; defaults to 1e-16.
* @param {boolean} opt.outputLagrangeMultiplier boolean indicating if the Lagrange multiplier associated to the optimal solution of the problem 
* must be provided in output (true) or not (false); defaults to false.
* @return {Array<Object>} an array arr containing:
* - If opt.outputLagrangeMultiplier is set to false, two elements:
* -- arr[0] an n by 1 matrix containing the optimal solution x^* to the problem
* -- arr[1] the optimal value of the function, f(x^*)
*
* - If opt.outputLagrangeMultiplier is set to true, three elements:
* -- arr[0] an n by 1 matrix containing the optimal solution x^* to the problem
* -- arr[1] the optimal value of the function, f(x^*) 
* -- arr[2] the Lagrange multiplier associated to the equality constraint of the problem, t^*
*
* @example
* qksolveBS_(Matrix_([1, 1]), Matrix_([1, 1]), Matrix_([1, 1]), 1, Matrix_([0, 0]), Matrix_([1, 1])); // Compute the projection of the point [1,1] on the standard simplex of R^2
* // [Matrix_([0.5, 0.5]), -0.75]
*/
function qksolveBS_(d, a, b, r, l, u, opt) {
	// Internal function to resize an array.
	function resizeArray(arr, n) {
		if (arr instanceof Array) { // this restrict the size of the array to the first n elements
			arr.length = n; 
			return arr;
		}
		else if (arr instanceof Float64Array || arr instanceof Int32Array) { // this constructs a view on the first n array elements
			return arr.subarray(0, n);
		}
	}
	
	// Internal function to compute g(t) = <b/x(t)>, 
	// c.f. remark 3.2 b) of the first reference.
	function g(t) {
		// Compute g(t) using the formula 3.5 of the first reference.
		
		// Initialize g(t) with the right part of formula 3.5.
		var g_t = (p - t * q) + s;
		
		// Finalize the computation of g(t) by adding the left part of formula 3.5,
		// sum b_i * x_i(t), i belonging to the set of indices I.		
		var I_it = new I.iterator();
		var el = I_it.next();
		while (el != 0) {
			// Get the index belonging to the set I
			var i = el - 1;

			
			var x_i;
			if (t <= T_u[i]) {
				x_i = u.data[i];
			}
			else if (T_u[i] <= t && t <= T_l[i]) {
				x_i = (a.data[i] - t*b.data[i]) / d.data[i];
			}
			else if (T_l[i] <= t) {
				x_i = l.data[i];
			}
			
			g_t += b.data[i] * x_i;
			
			
			// Prepare the next iteration
			el = I_it.next();
		}
		
		// Return the value of g(t).
		return g_t;
	}
	
	// Internal function to compute the vector x(t).
	function x(t) {
		// Initialize x(t).
		var x_t = Matrix_.zeros(n, 1);
		
		// Compute x(t) componentwise, using formula 2.6 of the first reference.
		for (var i = 0; i < n; ++i) {
			var x_i;
			if (t_star <= T_u[i]) {
				x_i = u.data[i];
			}
			else if (T_u[i] <= t_star && t_star <= T_l[i]) {
				x_i = (a.data[i] - t_star*b.data[i]) / d.data[i];
			}
			else if (T_l[i] <= t_star) {
				x_i = l.data[i];
			}
			x_t.data[i] = x_i;
		}
		
		// Return the value of x(t).
		return x_t;
	}
	
	
    // ------
    
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-16;
	var outputLagrangeMultiplier = false || opt.outputLagrangeMultiplier;
	

	// ------

	// Misc. checks
	if (!(d instanceof Matrix_) || !d.isVector()) {
		throw new Error('first input must be a vector');
	}
	if (!(a instanceof Matrix_) || !a.isVector()) {
		throw new Error('second input must be a vector');
	}
	if (!(b instanceof Matrix_) || !b.isVector()) {
		throw new Error('third input must be a vector');
	}
	if (!(l instanceof Matrix_) || !l.isVector()) {
		throw new Error('fifth input must be a vector');
	}
	if (!(u instanceof Matrix_) || !u.isVector()) {
		throw new Error('sixth input must be a vector');
	}
	
	if (d.nbRows !== a.nbRows) {
		throw new Error('first and second inputs number of rows do not match: ' + d.nbRows + '-' + a.nbRows);
	}
	if (d.nbRows !== b.nbRows) {
		throw new Error('first and third inputs number of rows do not match: ' + d.nbRows + '-' + b.nbRows);
	}
	if (d.nbRows !== l.nbRows) {
		throw new Error('first and fifth inputs number of rows do not match: ' + d.nbRows + '-' + l.nbRows);
	}
	if (d.nbRows !== u.nbRows) {
		throw new Error('first and sixth inputs number of rows do not match: ' + d.nbRows + '-' + u.nbRows);
	}


	// ------
	
	// Initialisations    
	var n = b.nbRows;
	
	var abs_r = Math.abs(r);
	
	var T = typeof Float64Array === 'function' ? new Float64Array(2*n) : new Array(2*n); // the set of breakpoints T
	var T_l = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the list of breakpoints t_l_i, i=1..n
	var T_u = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the list of breakpoints t_u_i, i=1...n
	
	var I = new BitSet_().setRange(1, n); // the set of indices I
	var p = 0;
	var q = 0;
	var s = 0;
	
	
	// ------
	
	// Computation of the breakpoints t_l_i and t_u_i, i = 1..n, c.f. formula 2.5 of the first reference.
	// 
	// In parallel:
	// - Computation of t_1 and t_r, c.f. formula 7 of the second reference.
	// - Basic checks on the problem constraints.
	var t_1 = Infinity;
	var t_r = -Infinity;
	for (var i = 0, j = 0; i < n; ++i) {
		// Check on lower and upper bounds l_i and u_i
		if (l.data[i] > u.data[i]) {
			throw new Error('infeasible problem detected');
		}
		
		// Check the strict positivity of b_i
		if (b.data[i] <= 0) {
			throw new Error('negative element detected in b');
		}

		// Check the strict positivity of d_i
		if (d.data[i] <= 0) {
			throw new Error('negative element detected in d');
		}		
	
		// Computation of t_l_i
		var t_l_i = (a.data[i] - l.data[i]*d.data[i]) / b.data[i];
		T_l[i] = t_l_i;
		T[j++] = t_l_i;
		
		// Computation of t_u_i
		var t_u_i = (a.data[i] - u.data[i]*d.data[i]) / b.data[i];
		T_u[i] = t_u_i;
		T[j++] = t_u_i;

		// Potential update of t_1 and t_r
		//
		// To be noted that as t_u_i <= t_l_i, i=1..n:
		// - t_1 is necessarily found amongst t_u_i, i=1..n
		// - t_r is necessarily found amongst t_l_i, i=1..n
		if (t_l_i > t_r) {
			t_r = t_l_i;
		}
		if (t_u_i < t_1) {
			t_1 = t_u_i;
		}	
	}

	// Check the feasibility of the problem , c.f. line 2 of the 
	// algorithm of the second reference.
	var g_t_1 = g(t_1);
	var g_t_r = g(t_r);
	if (g_t_1 < r || g_t_r > r) {
		throw new Error('infeasible problem detected');
	}

	// If the problem is feasible, it admits a unique solution x(t^*), and the problem
	// is then to compute t^*, c.f. the theorem of the second reference.
	var t_star = null;

	// Check if t_1 or t_r is optimal, in which case the algorithm can be stopped, c.f.
	// line 1 of the algorithm of the second reference.	
	if (Math.abs(g_t_1 - r) <= eps * abs_r) {
		t_star = t_1;
	}
	else if (Math.abs(g_t_1 - r) <= eps * abs_r) {
		t_star = t_r;
	}
	// Otherwise, proceed with the core algorithm 3.1 of the first reference.
	else {
		// Step 0: initialisations, with some initialisations already done.
		var t_l = t_1; // t_1 was already computed to check feasibility, so there is no need to use t_0
		var t_u = t_r; // t_r was already computed to check feasibility, so there is no need to use t_rp

		while (T.length != 0) {
			// Step 1: Breakpoint selection
			// 
			// This step uses the SELECT algorithm of Floyd and Rivest to compute the median,
			// with added benefit for steps 4 and 5 below that this algorithm partition the array T 
			// with the elements lower than the median at the left of the median position and the 
			// elements greater than the median at the right of the median position.
			var median_indice = Math.ceil(T.length/2);
			var t_hat = select_(T, median_indice);

			// Step 2: Computing g(t^)
			var g_t_hat = g(t_hat);

			// Step 3: Optimality check
			if (Math.abs(g_t_hat - r) <= eps * abs_r) {
				t_star = t_hat;
				break;
			}
			
			// Step 4: Lower breakpoint removal
			else if (g_t_hat > r) {
				t_l = t_hat;

				// Update T, with T = {t in T : t^ < t}.
				//
				// Thanks to the SELECT algorithm used to compute the median
				// of T, updating T simply means extracting the elements T.length/2...T.length
				var j = 0;
				for (var i = median_indice; i < T.length; ++i) {
					T[j++] = T[i];
				}
				T = resizeArray(T, j);
			}

			// Step 5: Upper breakpoint removal
			else if (g_t_hat < r) {
				t_u = t_hat;
				
				// Update T, with T = {t in T : t < t^}.
				//
				// Thanks to the SELECT algorithm used to compute the median
				// of T, updating T simply means extracting the elements 0...T.length/2
				T = resizeArray(T, median_indice - 1);

			}
			
			// Step 6: 
			// - Update of I, p, q and s following the formula 3.8 of the first reference
			//
			// The elements of I which are kept after all the iterations below are copied
			// at the beginning of the array I, and the array I is resized to the number
			// of kept elements.		
			var I_it = new I.iterator();
			var el = I_it.next();
			while (el != 0) {
				// Get the index belonging to the set I
				var i = el - 1;
				var remove_i = false;
				
				if (T_l[i] <= t_l) {
					s += b.data[i] * l.data[i];

					remove_i = true;
				}
				if (t_u <= T_u[i]) {
					s += b.data[i] * u.data[i];

					remove_i = true;
				}
				if (T_u[i] <= t_l && t_u <= T_l[i]) {
					var b_d = b.data[i] / d.data[i];
					p += a.data[i] * b_d;
					q += b.data[i] * b_d;

					remove_i = true;
				}
				
				if (remove_i === true) {
					I.unset(i+1);
				}
				
				// Prepare the next iteration
				el = I_it.next();
			}
		}

		// Step 6: 
		// - Stopping criterion: test on the size of T
		if (T.length == 0) {
			t_star = (p + s - r) / q;
		}
	}

	// Now that t^* has been computed, the last step is to compute the optimal
	// solution of the problem, x^* = x(t^*), c.f. remark 3.2 d) of the first reference.
	var x_star = x(t_star);
	
	// Compute the optimal function value f(x^*).
	var fctVal = 1/2 * Matrix_.vectorDotProduct(x_star, Matrix_.elementwiseProduct(x_star, d)) - Matrix_.vectorDotProduct(a, x_star);
	
	// Return the computed solution.
	if (outputLagrangeMultiplier === true) {
		return [x_star, fctVal, t_star];
	}
	else {
		return [x_star, fctVal];
	}
}
 
 
/**
* @function qpsolveGSMO_
*
* @summary Returns an optimal solution to a quadratic program with a single linear constraint
* and finite bound constraints, using a generalized sequential minimization optimization algorithm.
*
* @description This function computes an optimal solution to a quadratic program
* with a single linear constraint and finite bound constraints using a GSMO algorithm, 
* c.f. the first reference.
*
* The quadratic program to solve is assumed to be provided in the following format:
*
* min f(x) = 1/2 * <Q*x/x> + <p/x>
*
* s.t. <b/x> = r (single linear equality constraint)
*      l <= x <= u (bound constraints)
*
* with:
* - Q an n by n square symetric positive semi-definite matrix
* - p an n by 1 matrix
* - r a real number
* - b an n by 1 matrix with strictly positive elements
* - l an n by 1 matrix
* - u an n by 1 matrix
* 
* To be noted that the algorithm used internally requires that the feasible set F of this quadratic program is non empty
* and that f is bounded below on F to converge (i.e., admit a finite optimal solution).
*
* Since the feasible set, if non empty, is bounded by definition, the main assumption 
* is then that the feasible set is non-empty (i.e., that the problem is feasible), 
* and if this is not the case, an error is returned.
*
* @see <a href="https://link.springer.com/article/10.1023/A:1012431217818">Keerthi, S. & Gilbert, E. Convergence of a Generalized SMO Algorithm for SVM Classifier Design Machine Learning (2002) 46: 351.</a>
* @see <a href="http://ieeexplore.ieee.org/document/6789464/">S. S. Keerthi, S. K. Shevade, C. Bhattacharyya and K. R. K. Murthy, "Improvements to Platt's SMO Algorithm for SVM Classifier Design," in Neural Computation, vol. 13, no. 3, pp. 637-649, March 1 2001.</a>
* @see <a href="https://www.ncbi.nlm.nih.gov/pubmed/15941003">Takahashi N, Nishi T., Rigorous proof of termination of SMO algorithm for support vector machines., IEEE Trans Neural Netw. 2005 May;16(3):774-6.</a>
*
* @param {Matrix_} Q a square symetric positive semi-definite n by n matrix.
* @param {Matrix_} p an n by 1 matrix.
* @param {Matrix_} b an n by 1 matrix with strictly positive elements.
* @param {number} r a real number.
* @param {Matrix_} l an n by 1 matrix.
* @param {Matrix_} u an n by 1 matrix.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-04.
* @param {number} opt.maxIter maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @return {Array<Object>} an array arr containing two elements: 
* - arr[0] an n by 1 matrix containing the optimal solution x^* (in case Q is positive definite) 
* or an optimal solution x^* (in case Q is positive semi-definite) to the quadratic program
* - arr[1] the optimal value of the function f, i.e. f(x^*)
*
* @example
* qpsolveGSMO_(Matrix_([[2, 1], [1, 1]]), Matrix_([0, 0]), Matrix_([1, 1]), 1, Matrix_([0, 0]), Matrix_([1, 1])); // Solves min x^2 + xy + y^2/2 on the unit simplex of R^2
* // [Matrix_([0, 1]), 0.5]
*/
 function qpsolveGSMO_(Q, p, b, r, l, u, opt) {
    // ------
    
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-04;
	var maxIterations = opt.maxIter || 10000;
	
	
	// ------

	// Misc. checks
	if (!(Q instanceof Matrix_) || !Q.isSquare()) {
		throw new Error('first input must be a square matrix');
	}
	if (!(p instanceof Matrix_) || !p.isVector()) {
		throw new Error('second input must be a vector');
	}
	if (!(b instanceof Matrix_) || !b.isVector()) {
		throw new Error('third input must be a vector');
	}
	if (!(l instanceof Matrix_) || !l.isVector()) {
		throw new Error('fifth input must be a vector');
	}
	if (!(u instanceof Matrix_) || !u.isVector()) {
		throw new Error('sixth input must be a vector');
	}
	
	if (Q.nbRows !== p.nbRows) {
		throw new Error('first and second inputs number of rows do not match: ' + Q.nbRows + '-' + a.nbRows);
	}
	if (Q.nbRows !== b.nbRows) {
		throw new Error('first and third inputs number of rows do not match: ' + Q.nbRows + '-' + b.nbRows);
	}
	if (Q.nbRows !== l.nbRows) {
		throw new Error('first and fifth inputs number of rows do not match: ' + Q.nbRows + '-' + l.nbRows);
	}
	if (Q.nbRows !== u.nbRows) {
		throw new Error('first and sixth inputs number of rows do not match: ' + Q.nbRows + '-' + u.nbRows);
	}

	
	// ------
	
	// Initializations
	var n = Q.nbRows;

	
	// ------
	
	// Implementation of the algorithm GSMO of the first reference,
	// with some implementation details provided in the second reference.
	
	// Compute a feasible initial point x.
	//
	// This is done below by projecting the "centroid" vector
	// r/n * (1/b_1,...,1/b_n) on the constraints set, which is an O(n)
	// operation.
	var centroid = Matrix_.fill(n, 1, 
								function(i,j) { 
									return r/n* 1/b.data[i-1];
								});
	var p_centroid = qksolveBS_(Matrix_.ones(n, 1), centroid, b, r, l, u);
	var x = p_centroid[0];
	
	// Compute the gradient of the function f at the point x, using formula
	// grad(f)(x) = Q*x + p.
	//
	// This step is the most expansive code portion, since matrix-vector
	// multiplication is O(n^2).
	var grad_f_x = Matrix_.xpy(Matrix_.xy(Q, x), p);
	
	// Main loop of the GSMO algorithm, which convergence is guaranteed
	// by theorem 1 of the first reference and by theorem 1 of the
	// third reference.
	var iter = 0;
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Update the number of iterations
		++iter;
		
		// Choose (i,j) a pair of indices such that (i,j) = (i_up,i_low) is 
		// the most violating pair on the sets I_up and I_low, c.f. formulas
		// 5.1a and 5.1b of the second reference:
		// - i_low is the indice such that F_i_low = max {F_i, i belonging to I_low}
		// - i_up is the indice such that F_i_up = min {F_i, i belonging to I_up}
		//
		// with:
		// - I_low = I_0 u I_3 u I_4
		// - I_up = I_0 u I_1 u I_2
		// - I_0 = {i, l(i) < x(i) < u(i)}
		// - I_1 = {i, b(i) > 0, x(i) = l(i)}
		// - I_2 = {i, b(i) < 0, x(i) = u(i)}, empty here as b > 0
		// - I_3 = {i, b(i) > 0, x(i) = u(i)}
		// - I_4 = {i, b(i) < 0, x(i) = l(i)}, empty here as b > 0
		//
		// This choise corresponds to the method 2 described at the point 5 
		// of the section 5 of the second reference, with no optimisation 
		// to compute (i,j) first on set I_0 and only then on the sets I_up and I_low.
		//
		// To be noted that the first reference does not describe this particular 
		// choice in details, because the GSMO algorithm described is more generic.
		var i_low = -1;
		var F_i_low = -Infinity;
		var i_up = -1;
		var F_i_up = Infinity;
		for (var i = 0; i < n; ++i) {
			// Compute F_i, c.f. formula 1 of the first reference.
			F_i = grad_f_x.data[i] / b.data[i];
			
			// If i belongs to I_0, both i_low and i_up can potentially be updated.
			if (l.data[i] < x.data[i] && x.data[i] < u.data[i]) {
				if (F_i > F_i_low) {
					F_i_low = F_i;
					i_low = i;
				}
				if (F_i < F_i_up) {
					F_i_up = F_i;
					i_up = i;
				}
			}
			// If i belongs to I_1, only i_up can potentially be updated.
			else if (x.data[i] == l.data[i]) {
				if (F_i < F_i_up) {
					F_i_up = F_i;
					i_up = i;
				}		
			}
			// If i belongs to I_3, only i_low can potentially be updated.
			else if (x.data[i] == u.data[i]) {
				if (F_i > F_i_low) {
					F_i_low = F_i;
					i_low = i;
				}		
			}
		}
		
		// Stopping condition: check the formula 6 of the first reference:
		// min {F_i, i belonging to I_up} >= max {F_i, i belonging to I_low} - eps,
		// which is equivalent by definition to checking F_i_low >= F_i_up - eps.
		if (F_i_low - F_i_up <= eps) {
			break;
		}
		
		// Minimization of the function f on the rectangle [l(i), u(i)] x [l(j), u(j)] 
		// while varying only the coordinates (i,j) of the point x, with
		// i = i_up and j = i_low per choice of the (i,j) indices above.
		// 
		// From the section 3 of the first reference, this problem is equivalent 
		// to minimizing a function phi(t) = phi(0) + phi'(0)*t + phi''(0)*t^2/2
		// with:
		// - phi(0) irrelevant for the analysis
		// - phi'(0) = F_i - F_j
		// - phi''(0) = Q(i,i)/b(i)^2 + Q(j,j)/b(j)^2 -2*Q(i,j)/(b(i)*b(j))
		// - t_min <= t <= t_max, with t_min and t_max determined such that
		// l(i) <= x(i) + t/b(i) <= u(i) and l(j) <= x(j) - t/b(j) <= u(j),
		// that is:
		// - t_min = max( (l(i) - x(i))*b(i) , (x(j) - u(j))*b(j) )
		// - t_max = min( (u(i) - x(i))*b(i) , (x(j) - l(j))*b(j) )
		//
		// As the matrix Q is supposed to be positive semi-definite, there are 
		// only two possibilities:
		// - phi''(0) > 0, in which case phi is a second order polynomial with a strictly
		// positive leading coefficient. The unconstrained minimum of phi is then reached 
		// at t^* = -phi'(0)/phi''(0), and the contrained minimum of phi is then reached at
		// max(t_min, min(t^*, t_max)).
		//
		// - phi''(0) = 0, in which case phi is a linear function. Since phi'(0) <> 0 
		// per selection of the pair (i,j) = (i_low,i_up), phi is not constant. The
		// minimum of phi is then reached at t^* = t_min if phi'(0) > 0, or at t^* = t_max
		// if phi'(0) < 0.
		var i = i_up;
		var j = i_low;
		
		// Compute t_min
		var t_min_l_i = (l.data[i] - x.data[i])*b.data[i];
		var t_min_u_j = (x.data[j] - u.data[j])*b.data[j];
		var t_min = t_min_l_i >= t_min_u_j ? t_min_l_i : t_min_u_j;

		// Compute t_max
		var t_max_u_i = (u.data[i] - x.data[i])*b.data[i];
		var t_max_l_j = (x.data[j] - l.data[j])*b.data[j];
		var t_max = t_max_u_i <= t_max_l_j ? t_max_u_i : t_max_l_j;
		
		// Compute t^*
		var dphi_0 = F_i_up - F_i_low;
		var ddphi_0 = Q.data[i*Q.nbColumns + i]/(b.data[i] * b.data[i]) + Q.data[j*Q.nbColumns + j]/(b.data[j] * b.data[j]) - 2*Q.data[i * Q.nbColumns + j]/(b.data[i] * b.data[j]);
		var t_star;
		if (ddphi_0 > 0) { // phi''(0) > 0
			t_star = -dphi_0/ddphi_0;
			if (t_star > t_max) {
				t_star = t_max;
			}
			else if (t_star < t_min) {
				t_star = t_min;
			}
		}
		else { // phi''(0) = 0, as phi''(0) < 0 would imply Q is not positive semi-definite
			if (dphi_0 > 0) {
				t_star = t_min;
			}
			else {
				t_star = t_max;
			}
		}
		
		// Once t^* minimizing phi on [t_min,t_max] has been computed,
		// the point x can be updated using the parametric values of its (i,j) coordinates,
		// c.f. section 3 of the first reference:
		// - x(i)_new = x(i)_old + t^*/b(i) 
		// - x(j)_new  = x(j)_old - t^*/b(j)
		//
		// Nevertheless, to avoid loss of numerical precision when either x(i) or x(j) reaches one
		// of its boundaries, a specific logic is implemented below to make sure both x(i) and x(j)
		// stay feasible.
		var old_x_i = x.data[i];
		var old_x_j = x.data[j];
		var delta_x_i = t_star/b.data[i];
		var delta_x_j = -t_star/b.data[j];		
		x.data[i] += delta_x_i;
		x.data[j] += delta_x_j;
		
		// Specific process in case of boundaries crossing.
		//
		// Note that there is no if/else below because it might happen that t_min == t_max
		// and/or that t_min_l_i == t_min_u_j and/or ...
		if (t_star == t_min) {
			if (t_min == t_min_l_i) {
				x.data[i] = l.data[i];
				delta_x_i = x.data[i] - old_x_i;
				// No numerical update for x(j), because if it has not reached a boundary, it is strictly within [l(j), u(j)]
			}
			if (t_min == t_min_u_j) {
				x.data[j] = u.data[j];
				delta_x_j = x.data[j] - old_x_j;
				// No numerical update for x(i), because if it has not reached a boundary, it is strictly within [l(i), u(i)]
			}
		}
		if (t_star == t_max) {
			if (t_max == t_max_u_i) {
				x.data[i] = u.data[i];
				delta_x_i = x.data[i] - old_x_i;
				// No numerical update for x(j), etc.
			}
			if (t_max == t_max_l_j) {
				x.data[j] = l.data[j];
				delta_x_j = x.data[j] - old_x_j;
				// No numerical update for x(i), etc.
			}		
		}		
		
		// Compute the gradient of the function f at the updated point x.
		//
		// To be noted that since the updated point x is different from the previous point x
		// by only its (i,j) coordinates, this gradient can be computed using the formula:
		// grad_f_x_new = Q*x_new + p 
		//              = Q*(x_old + (0...t^*/b(i)...0...-t^*/b(j)...0)^t) + p
		//              = Q*x_old + p + Q*(0...t^*/b(i)...0...-t^*/b(j)...0)^t
		//              = grad_f_x_old + (Q(1,i)*t^*/b(i) - Q(1,j)*t^*/b(j), ..., Q(n,i)*t^*/b(i) - Q(n,j)*t^*/b(j))
		for (var k = 0; k < n; ++k) {
			grad_f_x.data[k] = grad_f_x.data[k] + Q.data[k*Q.nbColumns + i]*delta_x_i + Q.data[k*Q.nbColumns + j]*delta_x_j;
		}
	}
	
	// Compute the optimal function value f(x^*).
	//
	// This step is also the most expansive code portion, since matrix-vector
	// multiplication is O(n^2).
	var fctVal = 1/2 * Matrix_.vectorDotProduct(x, Matrix_.xy(Q, x)) + Matrix_.vectorDotProduct(p, x);
	
	// Return the computed solution.
	return [x, fctVal];
}

/**
* @function lpsolvePDHG_
*
* @summary Returns an optimal solution to a linear program, using a primal-dual hybrid gradient algorithm.
*
* @description This function computes an optimal solution to a linear program using a 
* preconditioned primal-dual hybrid gradient (PDHG) algorithm, c.f. the first reference.
*
* The linear program to solve is assumed to be provided in the following format:
*
* min f(x) = <c/x>
*
* s.t. Ae*x = be (equality constraints)
*      Ai*x <= bi (inequality constraints)
*      lb <= x <= ub (bound constraints)
*
* with:
* - c an n by 1 matrix
* - Ae an optional me by n matrix
* - be an optional me by 1 matrix
* - Ai an optional mi by n matrix
* - bi an optional mi by 1 matrix
* - lb an optional n by 1 matrix, which can contain negative infinity values (-Infinity) corresponding to unbounded variables on the negative axis
* - ub an optional n by 1 matrix, which can contain positive infinity values (Infinity) corresponding to unbounded variables on the positive axis
*
* and with:
* - lb assumed to be an n by 1 matrix made of zeroes if not provided
* - ub assumed to be an n by 1 matrix made of positive infinity values if not provided
* 
* To be noted that the algorithm used internally requires the linear problem to be feasible and bounded to converge
* (i.e., admit a finite optimal solution).
* 
* @see <a href="http://ieeexplore.ieee.org/document/6126441/">T. Pock and A. Chambolle, "Diagonal preconditioning for first order primal-dual algorithms in convex optimization" 2011 International Conference on Computer Vision, Barcelona, 2011, pp. 1762-1769.</a>
* @see <a href="https://arxiv.org/abs/1305.0546">Tom Goldstein, Min Li, Xiaoming Yuan, Ernie Esser, Richard Baraniuk, "Adaptive Primal-Dual Hybrid Gradient Methods forSaddle-Point Problems", 05/2013, eprint arXiv:1305.0546</a>
* @see <a href="http://www.numerical.rl.ac.uk/reports/drRAL2001034.pdf">D. Ruiz ,A scaling algorithm to equilibrate both rows and column norms in matrices, Tech.Report RT/APO/01/4, ENSEEIHT-IRIT, 2001.</a>
*
* @param {Matrix_} Ae an optional me by n matrix; must be null if not provided.
* @param {Matrix_} be an optional me by 1 matrix; must be null if not provided.
* @param {Matrix_} Ai an optional mi by n matrix; must be null if not provided.
* @param {Matrix_} bi an optional mi by 1 matrix; must be null if not provided.
* @param {Matrix_} c an n by n matrix.
* @param {Matrix_} lb an optional n by 1 matrix; must be null if not provided.
* @param {Matrix_} ub an optional n by 1 matrix; must be null if not provided.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-08.
* @param {number} opt.maxIter maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100000.
* @return {Array<Object>} an array arr containing two elements: 
* - arr[0] an n by 1 matrix containing the optimal solution x^* to the linear program
* - arr[1] the optimal value of the function f, i.e. f(x^*)
*
* @example
* lpsolvePDHG_(Matrix_([[1, 1]]), Matrix_([1]), null, null, Matrix_([1, 2]), null, null); // Solves min x + 2*y on the unit simplex of R^2
* // [Matrix_([~1, ~0]), ~1]
*/
 function lpsolvePDHG_(Ae, be, Ai, bi, c, lb, ub, opt) {
    // ------
    
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-08;
	var maxIterations = opt.maxIter || 100000;
	

	// ------

	// Misc. checks
	var eqContraints = false;
	var ineqContraints = false;
	var boundContraints = false;
	
	if (Ae !== null && be !== null) {
		eqContraints = true;
	}
	else if (Ae !== null && be === null) {
		throw new Error('equality constraints vector is missing');
	}
	else if (Ae === null && be !== null) {
		throw new Error('equality constraints matrix is missing');
	}

	if (Ai !== null && bi !== null) {
		ineqContraints = true;
	}
	else if (Ai !== null && bi === null) {
		throw new Error('inequality constraints vector is missing');
	}
	else if (Ai === null && bi !== null) {
		throw new Error('inequality constraints matrix is missing');
	}
	
	if (lb !== null && ub !== null) {
	    boundContraints = true;
	}
	else if (lb !== null && ub === null ) {
    	throw new Error('upper bounds constraints vector is missing');
	}
	else if (lb === null && ub !== null ) {
    	throw new Error('lower bounds constraints vector is missing');
	}
	
	if (!(c instanceof Matrix_)) {
		throw new Error('fifth input must be a matrix');
	}
	if (c.nbColumns !== 1) {
		throw new Error('fifth input is not a vector: ' + c.nbColumns + '-' + c.nbRows);
	}
	
	if (eqContraints) {
		if (!(Ae instanceof Matrix_)) {
			throw new Error('first input must be a matrix');
		}
		if (!(be instanceof Matrix_)) {
			throw new Error('second input must be a matrix');
		}
		if (Ae.nbRows !== be.nbRows) {
			throw new Error('first and second inputs number of rows do not match: ' + Ae.nbRows + '-' + be.nbRows);
		}
		if (Ae.nbColumns !== c.nbRows) {
			throw new Error('first input number of columns and fifth input number of rows do not match: ' + Ae.nbColumns + '-' + c.nbRows);
		}
		if (be.nbColumns !== 1) {
			throw new Error('second input is not a vector: ' + be.nbColumns + '-' + be.nbRows);
		}
	}
	
	if (ineqContraints) {
		if (!(Ai instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (!(bi instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (Ai.nbRows !== bi.nbRows) {
			throw new Error('third and fourth inputs number of rows do not match: ' + Ai.nbRows + '-' + bi.nbRows);
		}
		if (Ai.nbColumns !== c.nbRows) {
			throw new Error('third input number of columns and fifth input number of rows do not match: ' + Ai.nbColumns + '-' + c.nbRows);
		}
		if (bi.nbColumns !== 1) {
			throw new Error('fourth input is not a vector: ' + bi.nbColumns + '-' + bi.nbRows);
		}
	}
	
	if (boundContraints) {
		if (lb.nbRows !== null) {
			if (!(lb instanceof Matrix_)) {
				throw new Error('sixth input must be a matrix');
			}
			if (lb.nbRows !== c.nbRows) {
				throw new Error('sixth input number of rows and fifth input number of rows do not match: ' + lb.nbRows + '-' + c.nbRows);
			}
		}
        if (ub.nbRows !== null) {
			if (!(ub instanceof Matrix_)) {
				throw new Error('seventh input must be a matrix');
			}
			if (ub.nbRows !== c.nbRows) {
				throw new Error('seventh input number of rows and fifth input number of rows do not match: ' + ub.nbRows + '-' + c.nbRows);
			}
		}
	}

	
	// ------
	
	// Initializations
	// Constraints
	var me = 0;
	var ye_k = null;
	var ye_kp = null;
	var res_ye_kp_ye_k = null;
	if (eqContraints) {
	    me = Ae.nbRows; // the number of equality constaints
    	ye_k = Matrix_.zeros(me, 1); // the dual iterate ye_k
    	ye_kp = Matrix_.zeros(me, 1); // the dual iterate ye_k+1
    	res_ye_kp_ye_k = Matrix_.zeros(me, 1); // the residual ye_kp - ye_k
	}
	var mi = 0;
	var yi_k = null;
	var yi_kp = null;
	var res_yi_kp_yi_k = null;
	if (ineqContraints) {
	    mi = Ai.nbRows; // the number of inequality (<=) constaints
    	yi_k = Matrix_.zeros(mi, 1); // the dual iterate yi_k
    	yi_kp = Matrix_.zeros(mi, 1); // the dual iterate yi_k+1
    	res_yi_kp_yi_k = Matrix_.zeros(mi, 1); // the residual yi_kp - yi_k
	}
	var m = me + mi; // the total number of constraints
	
    // Variables
    var n = c.nbRows;
	var x_k = Matrix_.ones(n, 1); // the primal iterate x_k
	var x_kp = Matrix_.ones(n, 1); // the primal iterate x_k+1
	var z_k = Matrix_.ones(n, 1); // the relaxed iterate z_k = 2*x_k+1 - x_k
	var res_x_kp_x_k = Matrix_.zeros(n, 1); // the residual x_kp - x_k

    // Misc.
	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var ttmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var tmp_vec_me = null;
	if (eqContraints) {
		tmp_vec_me = Matrix_.zeros(me, 1); // a temporary placeholder vector of dimension me
	}
	var tmp_vec_mi = null;
	if (ineqContraints) {
		tmp_vec_mi = Matrix_.zeros(mi, 1); // a temporary placeholder vector of dimension mi
	}
	
	// Computation of the diagonal matrices T and S = [Se Si]^t with alpha = 1
	// and mu = 0.9995, nu = 0.9995 (so that mu * nu < 1), c.f. formula 10 and
	// remark 3 of the first reference.
	//
	// Note: in case a column or a row of the matrix [Ae Ai]^t has all its elements equal to zero,
	// the associated value in T or S is replaced by 1; this is standard practice, c.f. for instance
	// section 2 of the third reference.
	var mu = 0.9995;
	var nu = 0.9995;
	var T = Matrix_.fill(n, 1, 
				function(i,j) { 
						var columnNorm = 0;
						if (eqContraints) {
						    var aeColNorm = Ae.vectorNorm('one', 'column', i); 
							columnNorm += aeColNorm;
						}
						if (ineqContraints) {
						    var aiColNorm = Ai.vectorNorm('one', 'column', i);
							columnNorm += aiColNorm;
						}
						if (columnNorm == 0) {
							columnNorm = 1;
						}
						return mu * 1/columnNorm;
				});
	
	var Se = null;
	if (eqContraints) {
	    Se = Matrix_.fill(me, 1, 
	                      function(i,j) { 
						    var aeRowNorm = Ae.vectorNorm('one', 'row', i);
							if (aeRowNorm == 0) {
								aeRowNorm = 1;
							}
							return nu * 1/aeRowNorm;
				          });
	}
	var Si = null;
	if (ineqContraints) {
	    Si = Matrix_.fill(mi, 1, 
	                      function(i,j) { 
						    var aiRowNorm = Ai.vectorNorm('one', 'row', i);
							if (aiRowNorm == 0) {
								aiRowNorm = 1;
							}
			        	    return nu * 1/aiRowNorm;
				          });
	}
	
	// ------
	
	// Main loop of the algorithm, c.f. formula 18 of the first reference.
	//
	// The convergence is guaranteed by the theorem 1 of the first reference
	// in case the primal linear program admit a finite optimal solution.
	//
	// To be noted that the formulation used below is slightly different from
	// the formulation in formula 17 of the first reference (equality constraints only).
	//
	// To arrive at this form, the formulation below mixes formula 17 of the first reference,
	// formula 26 of the second reference, and additional bound constraints on x.
	var iter = 0;
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Update the number of iterations
		++iter;

		// Primal update: 
		// - x_k+1 = proj(x_k - T*(A^t*y_k + c))_[0, +infinity[, with A = [Ae Ai]^t and y_k = [ye yi]^t (no bound constraints)
		// - x_k+1 = proj(x_k - T*(A^t*y_k + c))_[lb, ub], with A = [Ae Ai]^t and y_k = [ye yi]^t (bound constraints)
		if (eqContraints && ineqContraints) {
			x_kp = Matrix_.xpy(Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), Matrix_.txy(Ai, yi_k, ttmp_vec_n), tmp_vec_n), c, x_kp);
		}
		else if (eqContraints) {
			x_kp = Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), c, x_kp);
		}
		else if (ineqContraints) {
			x_kp = Matrix_.xpy( Matrix_.txy(Ai, yi_k, tmp_vec_n), c, x_kp);
		}
		x_kp = Matrix_.xmy(x_k, Matrix_.elementwiseProduct(x_kp, T, tmp_vec_n), x_kp);
		if (boundContraints) {
			// Projection on the interval [lb, ub]
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < lb.data[i-1]) {
					x_kp.data[i-1] = lb.data[i-1];
				}
				else if (x_kp.data[i-1] > ub.data[i-1]) {
					x_kp.data[i-1] = ub.data[i-1];
				}
			}
		}
		else {			
			// Projection on the non-negative orthant
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < 0) {
					x_kp.data[i-1] = 0;
				}
			}
		}

		// Relaxed iterate update:
		// - z_k = 2*x_k+1 - x_k
		z_k = Matrix_.axpby(2, x_kp, -1, x_k, z_k);
		
		// Dual update:
		// - ye_k+1 = ye_k + Se*(Ae*z_k - be) (equality constraints)
		// - yi_k+1 = proj(yi_k + Si*(Ai*z_k - bi))_[0, +infinity[ (inequality constraints)
		if (eqContraints) {
			ye_kp = Matrix_.xpy(ye_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ae, z_k, tmp_vec_me), be, tmp_vec_me), Se, tmp_vec_me), ye_kp);
		}
		if (ineqContraints) {
			yi_kp = Matrix_.xpy(yi_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ai, z_k, tmp_vec_mi), bi, tmp_vec_mi), Si, tmp_vec_mi), yi_kp);
			
			// Projection on the non-negative orthant
			for (var i = 1; i <= mi; ++i) {
				if (yi_kp.data[i-1] < 0) {
					yi_kp.data[i-1] = 0;
				}
			}
		}
	
		// Convergence conditions for (x_k, y_k = [ye yi]^t) to be a saddle point of the min-max problem:
		// - Convergence of the primal iterates (relative) x_k: ||x_k+1 - x_k||_inf <= eps * ||x_k+1||_inf
		// - Convergence of the dual iterates (relative) y_k: ||y_k+1 - y_k||_inf <= eps * ||y_k+1||_inf
		res_x_kp_x_k = Matrix_.xmy(x_kp, x_k, res_x_kp_x_k);
		var res_x_kp_x_k_inf_norm = res_x_kp_x_k.vectorNorm('infinity');
		var x_kp_inf_norm = x_kp.vectorNorm('infinity');
		
		var res_ye_kp_ye_k_inf_norm = 0;
		var ye_kp_inf_norm = 0;
		if (eqContraints) {
		    res_ye_kp_ye_k = Matrix_.xmy(ye_kp, ye_k, res_ye_kp_ye_k);
		    res_ye_kp_ye_k_inf_norm = res_ye_kp_ye_k.vectorNorm('infinity');
											 
			ye_kp_inf_norm = ye_kp.vectorNorm('infinity');
		}
		
		var res_yi_kp_yi_k_inf_norm = 0;
		var yi_kp_inf_norm = 0;
		if (ineqContraints) {
		    res_yi_kp_yi_k = Matrix_.xmy(yi_kp, yi_k, res_yi_kp_yi_k);
		    res_yi_kp_yi_k_inf_norm = res_yi_kp_yi_k.vectorNorm('infinity');
			
			yi_kp_inf_norm = yi_kp.vectorNorm('infinity');
		}
		
		if (res_x_kp_x_k_inf_norm <= eps * x_kp_inf_norm  && 
		    res_ye_kp_ye_k_inf_norm <= eps * ye_kp_inf_norm &&
			res_yi_kp_yi_k_inf_norm <= eps * yi_kp_inf_norm) {
			break;
		}
		
		// Prepare the next iteration:
		// - x_k = x_k+1
		// - y_k = y_k+1 <=> ye_k = ye_k+1, yi_k = yi_k+1
		x_k = Matrix_.copy(x_kp, x_k);
		if (eqContraints) {
		    ye_k = Matrix_.copy(ye_kp, ye_k);
		}
		if (ineqContraints) {
		    yi_k = Matrix_.copy(yi_kp, yi_k);
		}
	}
	
	// Compute the objective function value.
	var fctVal = Matrix_.vectorDotProduct(c, x_kp);

	// Return the computed primal iterate and the associated objective
	// function value.
	return [x_kp, fctVal];		
 }
 
 
