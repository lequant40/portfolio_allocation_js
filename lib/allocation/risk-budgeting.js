/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function riskBudgetingWeights
*
* @summary Compute the weights of the risk budgeting portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with risk budgeting constraints.
*
* This portfolio has the property that the total contribution 
* of each asset to the risk of the portfolio is equal to a pre-determined budget weight.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* To be noted that in case weights constraints are defined, the concept of risk budget does not make any sense, 
* c.f. the sixth reference, and the associated optimization problem might not have any solution.
*
* The algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is semi-definite positive.
*
* @see <a href="https://ssrn.com/abstract=2009778">Bruder, Benjamin and Roncalli, Thierry, Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012).</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* @see <a href="https://link.springer.com/article/10.1023/A:1017501703105">Tseng P., Convergence of a Block Coordinate Descent Method for Nondifferentiable Minimization, Journal of Optimization Theory and Applications, 109(3), pp. 475-494. (2001)</a>
* @see <a href="https://doi.org/10.1090/mcom/3530">Stephen J. Wright, Ching-pei Lee, Analyzing random permutations for cyclic coordinate descent, Math. Comp.</a>
* @see <a href="http://proceedings.mlr.press/v29/Glasmachers13.html">Tobias Glasmachers, Urun Dogan, Accelerated Coordinate Descent with Adaptive Coordinate Frequencies, Proceedings of the 5th Asian Conference on Machine Learning, PMLR 29:72-86, 2013.</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3331184">RONCALLI Thierry, RICHARD Jean-Charles, CONSTRAINED RISK BUDGETING PORTFOLIOS, </a>
*
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {Array.<number>} rb the risk budgets, array of n real strictly positive numbers summing to one.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-10.
* @param {number} opt.epsSdp the tolerance parameter for testing the semi-definite positiveness of the covariance matrix, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.maxCycles the maximum number of cycles of the algorithm, a strictly positive natural integer or -1 to allow an infinite number of cycles; defaults to 10000.
* @param {number} opt.nbCycles the exact number of cycles of the algorithm, a strictly positive natural integer, in which case the values of opt.eps and opt.maxCycles are discarded,
* or -1; defaults to -1.
* @param {string} opt.coordinatesSampler, the type of coordinates sampler to use, a string either equal to:
* - 'cyclic', in order to use a cyclic coordinates sampler, called "cyclic CD" in the fourth reference
* - 'shuffledCyclic', in order to use a uniformly randomly shuffled cyclic coordinates sampler, called "random-permutations CD" in the fourth reference
* - 'randomized', in order to use a uniformly randomized coordinates sampler, called "fully randomized CD" in the fourth reference
* - 'acf', in order to use the adaptive coordinate frequencies coordinates sampler, as described in the fifth reference
*; defaults to 'cyclic'.
* @param {boolean} opt.outputPortfolioVolatility a boolean indicating whether the portfolio volatility should be provided in output
* (if set to true) or not (if set to false); defaults to false.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
* @return {Array.<number>|Array.<Array.<number>>} if opt.outputPortfolioVolatility is set to false, the weights corresponding to the risk budgeting portfolio, 
* array of n real numbers, and if opt.outputPortfolioVolatility is set to true, an array arr of two elements:
* - arr[0], the weights corresponding to the risk budgeting portfolio, array of n real numbers
* - arr[1], the volatility of the computed risk budgeting portfolio, a real number
*
* @example
* riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75], {eps: 1e-10, maxCycles: 1000});
* // [~0.45, ~0.55]
*/
self.riskBudgetingWeights = function (sigma, rb, opt) {
	// Internal function to compute the volatility of a portfolio
	function computeVolatility(x, sigma_x) {
		// Compute the variance x'*SIGMA*x
		var sigma_x_x = Matrix_.vectorDotProduct(sigma_x, x);		

		// In case the variance is numerically zero, which can occur with 
		// a semi-positive definite covariance matrix, it is replaced with zero.
		//
		// Otherwise, if the variance is negative, stops the algorithm,
		// since the covariance matrix is then not numerically semi-positive definite.
		if (Math.abs(sigma_x_x) <= epsSdp) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -epsSdp) {
			throw new Error('negative volatility, covariance matrix might not be semi-definite positive');
		}
		
		// Compute the volatility SQRT(x'*SIGMA*x)
		var s_x = Math.sqrt(sigma_x_x);
	   
		// Return the computed volatility
		return s_x;
	}

	
	// Internal function to compute the solution x^*(lambda) to the minimization problem
	// described in the second and sixth references, using a coordinate descent
	// method.
	//
	// It returns the computed solution x^*(lambda).
	function coordinatesDescentSolve(lambda, lowerBounds, upperBounds) {
		// Initial point for the algorithm is an equal weight vector
		var x = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
		var log_x = x.elemMap(function(i,j,val) { return Math.log(val);});
		
		// Placeholder for the previous x
		var x_old = Matrix_.copy(x);

		// Preparational computations
		var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
		var s_x = computeVolatility(x, sigma_x); // sigma(x)
		
		var rb_d_log_x = Matrix_.vectorDotProduct(rb, log_x); // Matrix_.vectorDotProduct(rb, log_x)
		var obj_old = 0.5 * s_x - rb_d_log_x; // the objective function to be minimized, c.f. formula 3 of the second reference
		
		// Initialization of the coordinates sampler
		var cs = new coordinatesSampler(nbAssets);
		
		// Main loop until convergence, guaranteed as per hypotheses on sigma and b
		// in case of essentially cyclic coordinates sampling.
		//
		// To be noted that what is guaranteed is the convergence of
		// the objective function values, c.f. the third reference, 
		// so that the convergence criteria will be based on the objective function values.
		//
		// To also be noted that the convergence of the x_k sequence is not guaranteed 
		// (the only guarantee is that all the cluster points of the x_k sequence 
		// are minimizers of the objective function).
		var cycle = 0;
		while (true) {
			// Check the exact number of cycles
			if	(nbCycles != -1 && cycle == nbCycles) {
				break;
			}

			// Check the maximum number of cycles
			if (nbCycles == -1 && maxCycles != -1 && cycle >= maxCycles) {
				throw new Error('maximum number of cycles reached: ' + maxCycles);
			}

			// Update the number of cycles
			++cycle;


			// Generate the coordinates to use in the next cycle of coordinates descent
			cs.generateCoordinates();

		   
			// Perform one full cycle of coordinates descent
			var obj_old_i
			var obj_new_i;
			while (true) {
				// Generate the index of the current coordinate to optimize in 1D,
				// in case such an index is remaining.
				var i = cs.sampleCoordinate();
				if (i == -1) {
					break;
				}

				// In case the coordinate sampler is ACF, 
				// compute the old value of the objective function
				if (opt.coordinatesSampler === 'acf') {
					obj_old_i = 0.5 * s_x - rb_d_log_x;
				}
				
				
				// Save the old asset weight before any update
				var xi_old = x.data[i-1];
				var log_xi_old = log_x.data[i-1];


				// Define the coefficients of the second order polynomial (a*x_i)^2 + b*x_i + c_i, c.f. the second reference
				var a = sigma.data[(i-1)*sigma.nbColumns + (i-1)]; // sigma_i^2, always > 0
				var b = sigma_x.data[i-1] - x.data[i-1] * a; // (SIGMA*x)_i - x_i*sigma_i^2, might be any sign
				var c = -lambda * rb.data[i-1] * s_x; // -lambda * b_i * sigma(x), always <= 0 (== 0 iff the covariance matrix is semi-definite positive and sigma(x) == 0)


				// Extract the strictly positive root x_i^* of the equation (a*x_i)^2 + b*x_i + c = 0, using a stable numerical formula
				var b_p = b/2; // reduced discriminant
				var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
				var disc = b_p*b_p - a*c;
				if (disc < 0) {
				   throw new Error('internal error: negative discriminant detected, the covariance matrix might not be semi-definite positive');
				}
				var q = -(b_p + sign_b_p * Math.sqrt(disc));
				var r1 = q/a;
				var r2 = c/q;
				
				var xi_star;
				if (r1 > 0) {
					xi_star = r1;
				}
				else if (r2 > 0) {
					xi_star = r2;
				}
				else {
					throw new Error('internal error: no strictly positive root detected, the covariance matrix might not be semi-definite positive');
				}
			   

				// Possibly truncate the strictly positive root above in case min/max weights are provided
				if (lowerBounds && xi_star < lowerBounds[i-1]) {
					xi_star = lowerBounds[i-1];
				}
				if (upperBounds && xi_star > upperBounds[i-1]) {
					xi_star = upperBounds[i-1];
				}
				
				
				// Update the asset weight
				x.data[i-1] = xi_star;
				
				// Update the asset weight log value, and <rb / log(x)>
				log_x.data[i-1] = Math.log(xi_star);
				
				rb_d_log_x -= log_xi_old * rb.data[i-1];
				rb_d_log_x += log_x.data[i-1] * rb.data[i-1];


				// Compute the updated SIGMA*x and SQRT(x'*SIGMA*x) elements for next loop evaluation.
				//
				// The update of the vector SIGMA*x uses the efficient update procedure described
				// in the second reference, based on the fact that only one coordinate of the vector x
				// changes per iteration.
				//
				// The update of the value x'*SIGMA*x does not use an "efficient" procedure,
				// because it requires a dot product, which is then equivalent to the full recomputation
				// of the volatility from SIGMA*x.

				// Compute the updated SIGMA*x
				for (var j = 1; j <= nbAssets; ++j) {
					sigma_x.data[j-1] += sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_star;
					sigma_x.data[j-1] -= sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_old;
				}

				// Compute the updated volatility SQRT(x'*SIGMA*x)
				s_x = computeVolatility(x, sigma_x);
				
				
				// In case the coordinate sampler is ACF:
				// - Compute the new value of the objective function
				// - Compute the gain for the current coordinate
				// - Update the scheduling preference for the current coordinate index i,
				//   or update the average gain
				if (opt.coordinatesSampler === 'acf') {
					// Compute the new value of the objective function
					obj_new_i = 0.5 * s_x - rb_d_log_x;
					
					// Compute the gain for the current coordinate, with 
					// gain = f(old) - f(new), from the code implementation of the fifth
					// reference by their authors.
					var gain = obj_old_i - obj_new_i;
					
					// If the cycle is the first one, update the average gain, and otherwise,
					// update the scheduling preference for the current coordinate index i.
					if (cycle == 1) {
						cs.updateAverageGain(gain);
					}
					else {
						cs.updateSchedulingPreference(i, gain);
					}
				}	
			}

		   
			// Compute the updated value of the objective function
			var obj_new = 0.5 * s_x - rb_d_log_x;


			// Check the necessary and sufficient convergence condition: |obj* - obj| <= eps,
			// unless an exact number of cycles is required.
			if (nbCycles == -1 && Math.abs(obj_new - obj_old) <= eps) {
				break;
			}
			
			
			// Prepare the next cycle:
			// - Update the previous objective function value
			obj_old = obj_new;
		}
		
		// Return the computed solution
		return x;
	}
	
	
	// Internal function to sample coordinates in a cyclic way.
	//
	// This is an essentially cyclic coordinates sampler, as
	// defined in the third reference.
	function cyclicCoordinatesSampler(n) {
		// Initializations
		this.n = n;
		this.k = n;
	   
	   
		// The sampling function, returning 1,2,...,n in this order,
		// and -1 to indicate that the sampling is finished.
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				return ++this.k;
			}        
		}
		
		
		// The function generating the coordinates for one cycle
		this.generateCoordinates = function() {
			// Re-set the index
			this.k = 0;
		}
	}

	
	// Internal function to sample coordinates in a randomly shuffled 
	// cyclic way.
	//
	// This is an essentially cyclic coordinates sampler, as
	// defined in the third reference.
	function shuffledCyclicCoordinatesSampler(n) {
		// Initializations
		this.n = n;
		this.k = n;
		
		// The coordinates set
		this.r = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.r[i] = i+1;
		}
	   
	   
		// The sampling function, returning the elements this.r[0],...,this.r[n-1] in this order,
		// and -1 to indicate that the sampling is finished.
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				return this.r[this.k++];
			}        
		}
		
		
		// The function generating the coordinates for one cycle
		this.generateCoordinates = function() {
			// Re-shuffle the coordinates
			this.r = new randomPermutationsIterator_(undefined, this.r, true).next();
			
			// Re-set the index
			this.k = 0;
		}
	}
	

	// Internal function to sample coordinates in a uniform random 
	// way.
	function randomizedCoordinatesSampler(n) {
		// Initializations
		this.n = n;
		this.k = n;
		
		var p = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			p[i] = 1/n;
		}
		this.s = new aliasMethodSampler_(p);
	   
	   
		// The sampling function, returning a coordinate chosen uniformly at random
		// among the set {1,...,n}, and -1 to indicate that the sampling is finished
		// after n samplings.
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				++this.k;
				return this.s.sample() + 1;
			}        
		}
		
		
		// The function generating the coordinates for one cycle
		this.generateCoordinates = function() {
			// Re-set the index
			this.k = 0;
		}
	}

	
	// Internal function to sample coordinates following the
	// adaptive coordinate frequencies of the fifth reference.
	//
	// This is an essentially cyclic coordinates sampler, as
	// defined in the third reference.
	function adaptiveCoordinateFrequenciesSampler(n) {
		// Initializations
		this.n = n;
		this.k = 0;
		
		// Initializations of the ACF strategy constants
		this.change_rate = 1/5; // c in the fifth reference
		this.p_min = 1/20;
		this.p_max = 20;	
		
		// Placeholder for the coordinates set I, allocated with the
		// maximum possible size of I which is 2n, c.f. the fifth reference.
		this.idx_set_tmp = typeof Uint32Array === 'function' ? new Uint32Array(2*n) : new Array(2*n);

		// The coordinates set I, null at initialization, which will contain 
		// between n and 2n coordinates for a given cycle, c.f. fifth reference.
		this.idx_set = null;
		this.idx_set_len = 0;
		
		// The preferences for scheduling p_i, i=1..n
		this.pref = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.pref[i] = 1;
		}
		this.prefsum = n; // sum p_i, i=1..n
		
		// The accumulators a_i, i=1..n
		//
		// Note: contrary to the fifth reference, the accumulators are not initialized to 0,
		// but to 0.5, which is the value the authors of the fifth reference used in their 
		// implementations.
		this.acc = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.acc[i] = 0.5;
		}
			
		// Variables for the self-adaptation of the performances monitoring
		this.gain_learning_rate = 1/n;
		this.average_gain = 0;

		
		// The sampling function, returning the elements this.r[0],...,this.r[n-1] in this order,
		// and -1 to indicate that the sampling is finished.
		this.sampleCoordinate = function() {
			if (this.k == this.idx_set_len) {
				return -1;
			}
			else {
				return this.idx_set[this.k++];
			}        
		}
		
		
		// The function generating the coordinates for one cycle
		this.generateCoordinates = function() {
			// Re-set the index
			this.k = 0;
			
			// Re-compute the index set
			this.idx_set_len = 0;
			var q = this.n / this.prefsum;
			for (var i = 0; i < this.n; ++i) {
				var a_i = this.acc[i] + q * this.pref[i];
				var a_i_f = Math.floor(a_i);
				
				for (var j = 0; j < a_i_f; ++j) {
					this.idx_set_tmp[this.idx_set_len] = i + 1;
					this.idx_set_len++;
				}
				
				this.acc[i] = a_i - a_i_f;
			}
		
			// Shuffle the index set
			this.idx_set = new randomPermutationsIterator_(undefined, this.idx_set_tmp.slice(0, this.idx_set_len), true).next();
		}
		
		
		//
		this.updateSchedulingPreference = function(i, gain) {
			// Compute the new scheduling preference for index i
			var p_new = this.pref[i-1] * Math.exp(this.change_rate * (gain / this.average_gain - 1));
			
			// Truncate the new scheduling preference for index i, if required
			if (p_new < this.p_min) {
				p_new = this.p_min;
			}
			else if (p_new > this.p_max) {
				p_new = this.p_max;
			}
			
			// Update the scheduling preference for index i
			this.prefsum = this.prefsum + p_new - this.pref[i-1];
			this.pref[i-1] = p_new;
			
			// Update the average gain
			this.average_gain = (1 - this.gain_learning_rate) * this.average_gain + this.gain_learning_rate * gain;
		}
		
		
		//
		this.updateAverageGain = function(gain) {
			this.average_gain += gain / this.n;
		}
	}
	
	// ------

	// Initialize the options structure
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	
	// Initialize the options default values
	if (opt.eps === undefined) {
		opt.eps = 1e-10;
	}
	if (opt.epsSdp === undefined) {
		opt.epsSdp = 1e-12;
	}
	if (opt.maxCycles === undefined) {
		opt.maxCycles = 10000;
	}
	if (opt.nbCycles === undefined) {
		opt.nbCycles = -1;
	}
	if (opt.outputPortfolioVolatility === undefined) {
		opt.outputPortfolioVolatility = false;
	}
	if (opt.coordinatesSampler === undefined) {
		opt.coordinatesSampler = 'cyclic';
	}
	
	// Decode the options
	var eps = opt.eps;
	var epsSdp = opt.epsSdp;
	var maxCycles = opt.maxCycles;
	var nbCycles = opt.nbCycles;
	
	var outputPortfolioVolatility = opt.outputPortfolioVolatility; 
	
	var coordinatesSampler;
	if (opt.coordinatesSampler === 'cyclic') {
		coordinatesSampler = cyclicCoordinatesSampler;
	}
	else if (opt.coordinatesSampler === 'shuffledCyclic') {
		coordinatesSampler = shuffledCyclicCoordinatesSampler;
	}	
	else if (opt.coordinatesSampler === 'randomized') {
		coordinatesSampler = randomizedCoordinatesSampler;
	}
	else if (opt.coordinatesSampler === 'acf') {
		coordinatesSampler = adaptiveCoordinateFrequenciesSampler;
	}
	else {
		throw new Error('unsupported coordinates sampler');
	}
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	
	// Convert rb to vector format
	var rb = new Matrix_(rb);
		
	// Decode the lower/upper bounds constraints, if provided
	var lowerBounds = null;
	var upperBounds = null;
	if (opt.constraints.minWeights) {
		lowerBounds = opt.constraints.minWeights; // no input checks
				
		if (!opt.constraints.maxWeights) {
			upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
			
			// Initialization to an array of ones
			for (var i = 0; i < upperBounds.length; ++i) {
				upperBounds[i] = 1;
			}
		}
	}
	if (opt.constraints.maxWeights) {
		upperBounds = opt.constraints.maxWeights; // no input checks
		
		if (!opt.constraints.minWeights) {
			lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);

			// Initialization to an array of zeros
			for (var i = 0; i < lowerBounds.length; ++i) {
				lowerBounds[i] = 0;
			}
		}
	}
	
	// Check that rb contains strictly positive numbers summing to one
	// Check that sigma and rb are rows compatible

	
	// ------
	
	
	// Compute the weights of the unconstrained risk budgeting portfolio, which
	// are needed even in case weights constraints are provided.
	//
	// This is done using any value of the lambda penalty parameter (below, 1)
	// and then normalizing the computed weights, c.f. the second reference.
	var lambda = 1;
	var weightsURb = coordinatesDescentSolve(lambda);
	weightsURb = weightsURb.normalize();
	
	
	// In case weights constraints are provided, the associated restricted simplex is first  
	// checked for emptiness, and the bisection algorithm described in the sixth reference is 
	// then used to compute the weights of the constrained risk budgeting portfolio.
	//
	// Otherwise, there is nothing else to do.
	var x;
	if (lowerBounds || upperBounds) {
		// In case the restricted simplex is empty, an exception is thrown, so that
		// the process is (violently) stopped here.
		var sumBounds = simplexEmptinessCheck_(nbAssets, lowerBounds, upperBounds);
		

		// Compute the initial interval for the bisection algorithm,
		// c.f. remark 7 of the sixth reference for the starting values of
		// a and b.
		//
		// Experimentations showed that the initial interval is sometimes not
		// a bracketing interval, so that the possibility to enlarge it has
		// been implemented below.
		//
		// In addition, experimentations showed that it can be that the sum of
		// the weights is numerically equal to 1, but with sign issues, so that the 
		// bisection below is not usable; this specific case is managed below.
		var volURb = computeVolatility(weightsURb, Matrix_.xy(sigma, weightsURb));
		
		var a = 0.5 * volURb;
		var weights_a = coordinatesDescentSolve(a, lowerBounds, upperBounds);
		var nbIterSearch = 0;
		var maxNbIterSearch = 54; // The default is taken to be 54, because 2^-54 is quite close to an already unreasonable low value
		var aEarlyStop = false;
		while ( weights_a.sum() - 1 > 0 ) {
			// Increment the number of iterations and check that the number of iterations stays reasonable
			++nbIterSearch;
			if (nbIterSearch > maxNbIterSearch) {
				throw new Error('internal error: maximum number of iterations reached when searching for a bracketing interval, the problem might be infeasible');
			}
	
			// It is possible that sum a_i is numerically equal to 1 !
			if ( Math.abs(weights_a.sum() - 1) <= eps ) {
				aEarlyStop = true;
				break;
			}
			
			// Halve the value of a, and compute the associated portfolio weights
			a = 0.5 * a;
			weights_a = coordinatesDescentSolve(a, lowerBounds, upperBounds);
		}
		
		var b = 2 * volURb;
		var weights_b = coordinatesDescentSolve(b, lowerBounds, upperBounds);
		var nbIterSearch = 0;
		var maxNbIterSearch = 54; // The default is taken to be 54, because 2^54 is quite close to an already unreasonable high value
		var bEarlyStop = false;
		while ( weights_b.sum() - 1 < 0 ) {
			// Increment the number of iterations and check that the number of iterations stays reasonable
			++nbIterSearch;
			if (nbIterSearch > maxNbIterSearch) {
				throw new Error('internal error: maximum number of iterations reached when searching for a bracketing interval, the problem might be infeasible');
			}

			// It is possible that sum b_i is numerically equal to 1 !
			if ( Math.abs(weights_b.sum() - 1) <= eps ) {
				bEarlyStop = true;
				break;
			}
			
			// Double the value of b, and compute the associated portfolio weights
			b = 2 * b;
			weights_b = coordinatesDescentSolve(b, lowerBounds, upperBounds);
		}
		
		
		if (aEarlyStop) {
			x = weights_a;
		}
		else if (bEarlyStop) {
			x = weights_b;
		}
		else {	
			// Compute the value of the penalty parameter lambda^* such that
			// sum x(lambda^*)_i = 1, i=1..nbAssets, using the bisection algorithm.
			var lambda_star = bisection_(function (lambda) { 
											 var weights = coordinatesDescentSolve(lambda, lowerBounds, upperBounds); 
											 return weights.sum() - 1; 
										 }, 
										 a, b);
			
			
			// Compute the associated portfolio weights
			x = coordinatesDescentSolve(lambda_star, lowerBounds, upperBounds);
		}
	}
	else {
		// The solution to the unconstrained risk budgeting optimization problem
		// has already been computed.
		x = weightsURb;
	}


	// Depending on what is requested in output, return the computed portfolio weights
	// and possibly the associated portfolio volatility.
	if (outputPortfolioVolatility === true) {
		return [x.toArray(), computeVolatility(x, Matrix_.xy(sigma, x))];
	}
	else {
		return x.toArray();
	}
}