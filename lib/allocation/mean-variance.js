/**
 * @file Functions related to mean variance efficient portfolios.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.efficientFrontier_ = efficientFrontier_;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function meanVarianceOptimizationWeights
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a target return
* or volatility constraint.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only efficient portfolio (i.e. belonging to the mean-variance efficient frontier) subject to
* either a target return constraint (in which case this portfolio has the lowest attainable volatility
* among all the portfolios satisfying the target return constraint) or a target volatility constraint
* (in which case this portfolio has the highest attainable return among all the portfolios satisfying 
* the target volatility constraint).
*
* The algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* To be noted that the portfolio volatility is defined as the standard deviation of the portfolio
* variance.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights corresponding to the computed mean-variance efficient portfolio, array of n real numbers.
*
* @example
* meanVarianceOptimizationWeights([0.1, 0.2], [[1, 0.3], [0.3, 1]], { constraints: {targetReturn: 0.15}})
* // [0.5, 0.5] 
*/
self.meanVarianceOptimizationWeights = function(mu, sigma, opt) {
	// Internal function to compute the (at most) two corner portfolios enclosing the
	// efficient portfolio with a return/volatility equals to the provided return/volatility.
	function computeEnclosingCornerPortfolios_(targetFct, targetFctVal, efficientFrontier) {
		// The numerical accuracy for testing equality
		var eps = 1e-8;
		
		// The efficient frontier portfolios are provided from highest return/variance
		// to lowest return/variance, so that *_min below refers to properties of the portfolio
		// with the lowest return/variance.
		var idx_min = efficientFrontier.length - 1;
		var idx_max = 0

		var weights_min = efficientFrontier[idx_min][0];
		var weights_max = efficientFrontier[idx_max][0];

		var fctVal_min = targetFct(weights_min);
		var fctVal_max = targetFct(weights_max);

		// If the target function value is not reachable within numerical accuracy, 
		// return immediately.
		if (targetFctVal - fctVal_max > eps || -eps > targetFctVal - fctVal_min) {
			return [];
		}
		
		// If the target function value has already been numerically reached on one of the
		// two extremal corner portfolios, return immediately.
		if (Math.abs(targetFctVal - fctVal_min) <= eps) {
			return [[weights_min, fctVal_min]];
		}
		else if (Math.abs(targetFctVal - fctVal_max) <= eps) {
			return [[weights_max, fctVal_max]];
		}
		
		// Otherwise, determine the two adjacent corner portfolios enclosing the portfolio
		// with a target function value numerically equals to the provided target function value,
		// using a binary search algorithm.
		//
		// Using a binary search algorithm is possible because the corner portfolios are
		// provided in decreasing return/variance values on the efficient frontier.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle point
			var idx_middle = Math.floor((idx_min + idx_max)/2);
			var weights_middle = efficientFrontier[idx_middle][0];
			var fctVal_middle = targetFct(weights_middle);
			
			// Determine in which sub-interval ]idx_max, idx_middle[ or ]idx_middle, idx_min[
			// lies the target function value.
			if (fctVal_middle - targetFctVal > eps) {
				idx_max = idx_middle;
				fctVal_max = fctVal_middle;
				weights_max = weights_middle;
			}
			else if (fctVal_middle - targetFctVal < -eps) {
				idx_min = idx_middle;
				fctVal_min = fctVal_middle;
				weights_min = weights_middle;
			}
			else { // the target function value is exactly attained on the idx_middle-th corner portfolio
				return [[weights_middle, fctVal_middle]];
			}
		}

		
		// Return the computed adjacent corner portfolios, as well as
		// the associated function values.
		return [[weights_min, fctVal_min], [weights_max, fctVal_max]];
	}
	
	
	// ------	

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	var targetReturn = opt.constraints.targetReturn;
	var targetVolatility = opt.constraints.targetVolatility;
	
	if (targetReturn === undefined && targetVolatility === undefined) {
		throw new Error('target return or target volatility is mandatory');
	}
	else if (targetReturn !== undefined && targetVolatility !== undefined) {
		throw new Error('target return and target volatility cannot be both provided');
	}
	
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	
	
	// ------
	
	var nbAssets = sigma.nbColumns; // the number of assets in the universe
	
	// Compute the efficient frontier
	var efficientFrontier = efficientFrontier_(mu, sigma, opt);
	
	// Set the target function and function value
	var targetFct;
	var targetFctVal;
	if (targetReturn !== undefined) { // the target function is the portfolio return
		targetFct = function(weights) {
			return Matrix_.vectorDotProduct(mu, weights);
		};
		targetFctVal = targetReturn;
	}
	else { // the target function is the portfolio volatility (i.e., standard deviation), convert it to variance
		targetFct = function(weights) {
			return Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights), weights);
		};
		targetFctVal = targetVolatility*targetVolatility;
	}

	// Compute the (at most) two efficient portfolios on the efficient frontier
	// which enclose the efficient portfolio with a target function value equals to 
	// the desired target function value.
	var efficientEnclosingPortfolios = computeEnclosingCornerPortfolios_(targetFct, targetFctVal, efficientFrontier);

	
	// Then:
	// - In case the desired target value function is not reachable, stop the process 
	//
	// - In case there is a unique computed efficient portfolio with a target function value
	// equals to the desired target function value, return the associated portfolio weights
	//
	// - In case there are two efficient portfolios (strictly) enclosing the efficient portfolio with 
	// a target function value equals to the desired target function value, the weights associated
	// to this efficient portfolio are a (strict) convex combination of the two computed enclosing
	// efficient portfolios (c.f. the reference): w = t*w_min + (1-t)*w_max, t in ]0,1[, 
	// with t now to be determined.
	var weights;
	if (efficientEnclosingPortfolios.length == 0) {
		throw new Error('target return or volatility not reachable');
	}
	else if (efficientEnclosingPortfolios.length == 1) {
		var weights_min = efficientEnclosingPortfolios[0][0];
		weights = weights_min;		
	}
	else {
		// Extract information about the computed efficient corner portfolios
		var weights_min = efficientEnclosingPortfolios[0][0];
		var fctVal_min = efficientEnclosingPortfolios[0][1];
		
		var weights_max = efficientEnclosingPortfolios[1][0];
		var fctVal_max = efficientEnclosingPortfolios[1][1];
		
		// Depending on the desired target function, the procedure to compute t above is different:
		// - If the target is return, E(w) = <mu/w> and by linearity of E, we have
		// E(w) = t*E(w_min) + (1-t)*E(w_max) and E(w) = targetReturn
		// <=>
		// t = (E(w_max) - targetReturn)/(E(w_max) - E(w_min))
		//
		// - If the target is volatility (i.e., standard deviation), let the volatility be V(x) = <Sigma*w/w>.
		// Then, by symmetry and bilinerarity of V, we have V(w) = t^2*V(w_min) + (1-t)^2*V(w_max) + 2*t*(1-t)*<Sigma*w_min/w_max>
		// and V(w) = targetVolatility^2
		//	<=> t is the solution belonging to ]0,1[ of the second order polynomial equation
		// t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) -2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) - targetVolatility^2 = 0
		var t;
		if (targetReturn !== undefined) {
			t = (fctVal_max - targetFctVal)/(fctVal_max - fctVal_min);
		}
		else {
			// Define the coefficients of the second order polynomial at^2 + tx + c
    	    var fctVal_cross = Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
			var a = fctVal_min + fctVal_max - 2 * fctVal_cross; // always >= 0, by semi-definite positivity of the covariance matrix
    	    var b = -2 * (fctVal_max - fctVal_cross); // 
    	    var c = fctVal_max - targetVolatility*targetVolatility; //always > 0
    	    
    	    // Extract the root t of the equation at^2 + bt + c = 0 belonging to ]0,1[, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
			    throw new Error('internal error, the covariance matrix might not be semi-definite positive');
			}
    	    var q = -(b_p + sign_b_p * Math.sqrt(disc));
    	    var r1 = q/a;
    	    var r2 = c/q;
			
			if (r1 > 0 && r1 < 1) {
				t = r1;
			}
			else if (r2 > 0 && r2 < 1) {
				t = r2;
			}
			else {
				throw new Error('internal error, the covariance matrix might not be semi-definite positive');
			}
		}
		
		// Compute the final efficient portfolio weights
		weights = Matrix_.fill(nbAssets, 1, 
								function(i,j) { 
									return t*weights_min.getValue(i, 1) + (1-t)*weights_max.getValue(i, 1); 
								})
	}
	
	// Return the computed portfolio weights
	return weights.toArray();
}


/**
* @function meanVarianceEfficientFrontier
*
* @summary Compute the weights of all the corner portfolios belonging to the mean-variance efficient frontier.
*
* @description This function returns the weights w_i1,...,w_in, i = 1..m,
* associated to the m fully invested and long-only corner portfolios defining the mean-variance
* efficient frontier.
*
* The algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights of all the corner portfolios, an array of arrays of n real numbers
*
* @example
* meanVarianceEfficientFrontier([0.1, 0.2], [[1, 0.3], [0.3, 1]])
* // [[0, 1]], [0.5, 0.5]] 
*/
self.meanVarianceEfficientFrontier = function(mu, sigma, opt) {
	// Compute the efficient frontier
	var efficientFrontier = efficientFrontier_(mu, sigma, opt);
	
	// Convert the output of the internal function efficientFrontier_
	// to a list of portfolios weights.
	var efficientFrontierWeights = new Array(efficientFrontier.length);
	for (var i = 0; i < efficientFrontier.length; ++i) {
		efficientFrontierWeights[i] = efficientFrontier[i][0].toArray();
	}
	
	// Return the computed list of portfolios weights
	return efficientFrontierWeights;
}


/**
* @function efficientFrontier_
*
* @summary Compute all the corner portfolios belonging to the mean-variance efficient frontier.
*
* @description This function returns the weights w_i1,...,w_in as well as the risk aversion parameters lambda_i,
* i = 1..m, associated to the m fully invested and long-only corner portfolios belonging to the mean-variance
* efficient frontier.
*
* The algorithm used internally is the Markowitz critical line algorithm, c.f. the first reference.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
* @see <a href="https://doi.org/10.1007/978-0-387-77439-8_12">Niedermayer A., Niedermayer D. (2010) Applying Markowitz’s Critical Line Algorithm. In: Guerard J.B. (eds) Handbook of Portfolio Construction. Springer, Boston, MA</a>
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<Object>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
*
* @example
* efficientFrontier_([0.1, 0.2], [[1, 0.3], [0.3, 1]])
* // [[new Matrix_([0, 1]), 7], [new Matrix_([0.5, 0.5]), 0]] 
*/
function efficientFrontier_(mu, sigma, opt) {	
	// The numerical tolerance for testing equality
	var eps = 1e-8; 
	
	// Internal object managing the statuses of the asset and lambda variables
	function variablesStatusManager_(nbAssets, nbEqualityConstraints) {
		// Variables statuses constants
		this.STATUS_UNDEF = -1;
		this.STATUS_IN = 0;
		this.STATUS_LOW = 1;
		this.STATUS_UP = 2;
		
		// The structure holding the variables status
		this.nbAssets = nbAssets;
		this.nbEqualityConstraints = nbEqualityConstraints;
		
		this.in = new BitSet_();
		this.in.resize(nbAssets + nbEqualityConstraints);
		this.out = new BitSet_();
		this.out.resize(nbAssets + nbEqualityConstraints);
		this.low = new BitSet_();
		this.low.resize(nbAssets + nbEqualityConstraints);
		this.up = new BitSet_();
		this.up.resize(nbAssets + nbEqualityConstraints);
			
		// Public functions to set the status of variales
		this.setIn = function(idx) {
			this.in.set(idx);
			this.out.unset(idx);
			this.low.unset(idx);
			this.up.unset(idx);
		}
		this.setOnLowerBound = function(idx) {
			this.low.set(idx);
			this.out.set(idx);
			this.in.unset(idx);
			this.up.unset(idx);
		}
		this.setOnUpperBound = function(idx) {
			this.up.set(idx);
			this.out.set(idx);
			this.in.unset(idx);
			this.low.unset(idx);
		}
		this.setLambdasIn = function() {
			for (var i = this.nbAssets + 1; i <= this.nbAssets + this.nbEqualityConstraints; ++i) {
				this.in.set(i);
				this.out.unset(i);
				this.low.unset(i);
				this.up.unset(i);
			}
		}
		this.setAssetsOnLowerBounds = function() {
			for (var i = 1; i <= this.nbAssets; ++i) {
				this.low.set(i);
				this.out.set(i);
				this.in.unset(i);
				this.up.unset(i);
			}
		}
		
		// Public functions to get the status of a variable
		this.isAsset = function(idx) {
			return (idx >= 1) && (idx <= this.nbAssets);
		}
		this.isLambda = function(idx) {
			return (idx >= this.nbAssets + 1) && (idx <= this.nbAssets + this.nbEqualityConstraints);
		}
		this.isIn = function(idx) {
			return this.in.get(idx);
		}
		this.isOnLowerBound = function(idx) {
			return this.low.get(idx);
		}
		this.isOnUpperBound = function(idx) {
			return this.up.get(idx);
		}
		this.isOut = function(idx) {
			return this.out.get(idx);
		}
		
		// Public functions to iterate over the different sets.
		this.getInIndexes = function() {
			return this.in.toArray();
		}
		this.getOutIndexes = function() {
			return this.out.toArray();
		}
	}

	// Internal function to compute the E-maximizing portfolio, 
	// c.f. the method "STARTING-SOLUTION" of the second reference.
	//
	// This function replaces the simplex algorithm described in the
	// chapter 8 of the first reference in case:
	// - The only equality constraint on the assets is that their weights
	// sum to one
	//
	// - The only inequality constraints on the assets are positive lower bounds
	// and upper bounds on their weights
	//
	// - There is a unique optimal solution to the E-maximizing portfolio linear 
	// program
	function computeMaxReturnPortfolio_(mu, lowerBounds, upperBounds) {		
		// Check that the problem is feasible (l_i <= u_i, sum l_i <= 1 and 1 <= sum u_i,
		// c.f. paragraph 12.3.1 of the second reference).
		var sum_lb = 0;
		var sum_ub = 0;
		for (var i = 1; i <= nbAssets; ++i) {
			var lb_i = lowerBounds.getValue(i, 1);
			sum_lb += lb_i;
			
			var ub_i = upperBounds.getValue(i, 1);
			sum_ub += ub_i;
			
			if (lb_i > ub_i) {
				throw new Error('infeasible problem detected');
			}
		}
		if (sum_lb > 1 || sum_ub < 1) {
			throw new Error('infeasible problem detected');
		}
		
		// Order the assets in descending order w.r.t. their returns
		var mu_idx = typeof Uint32Array === 'function' ? new Uint32Array(nbAssets) : new Array(nbAssets);
		for (var j = 0; j < nbAssets; ++j) {		
			mu_idx[j] = j + 1;
		}
		mu_idx.sort(function(a, b) { 
			return mu.getValue(b, 1) - mu.getValue(a, 1);
		});

		// Check that the assets returns are all distinct, which is a sufficient condition
		// for the unicity of the E-maximizing portfolio.
		for (var i = 1; i < nbAssets; ++i) {
			if (mu.getValue(mu_idx[i], 1) == mu.getValue(mu_idx[i-1], 1)) {
				throw new Error('unsupported problem detected');
			}
		}

		// Initialize the E-maximizing portfolio weights with the assets lower bounds
		var maxReturnWeights = new Matrix_(lowerBounds);
		
		// Set the assets statuses to LOW
		variablesStatusManager.setAssetsOnLowerBounds();
	
		// Starting from the asset with the highest return, set each asset weight to its
		// highest possible value until the sum of the weights of all the assets is equal
		// to one.
		var delta_sum_weights = 1 - maxReturnWeights.sum();
		var idx_i = -1;
		for (var i = 0; i < nbAssets; ++i) {	
			// In case the new delta sum of the weights of all the assets is
			// numerically equal to zero, the loop can be stopped.
			if (Math.abs(delta_sum_weights) <= eps) {
				break;
			}
			
			// Extract the asset index and its current weight
			idx_i = mu_idx[i];
			var weight_asset = maxReturnWeights.getValue(idx_i, 1);
						
			// Compute the highest possible value for the increment in the asset weight
			var inc_weight_asset = Math.min(upperBounds.getValue(idx_i, 1) - weight_asset, delta_sum_weights);
			
			// Set the new weight of the asset, together with its status
			var new_weight_asset = weight_asset + inc_weight_asset;
			if (new_weight_asset >= upperBounds.getValue(idx_i, 1)) {			
				// In this case, as the highest possible weight for an asset is its upper bound, 
				// the asset weight must be capped to its upper bound.
				maxReturnWeights.setValue(idx_i, 1, upperBounds.getValue(idx_i, 1));
				
				// Set the asset UP status
				variablesStatusManager.setOnUpperBound(idx_i);
			}
			else {
				 // In this case, the asset lies strictly between its lower and upper bounds,
				 // and the new delta sum below will be zero.
				maxReturnWeights.setValue(idx_i, 1, new_weight_asset);
				
				// Set the asset IN status
				variablesStatusManager.setIn(idx_i);
			}
					
			// Compute the new delta sum of the weights of all the assets for the next iteration.
			//
			// Note: doing the computation this way (i.e. without calling .sum() again) allows
			// for a more efficient algorithm, at the price of a slight loss of numerical 
			// precision.
			delta_sum_weights -= inc_weight_asset;

		}

		// At this stage, there are four possibilities:
		// - The loop above has not started because the sum of the initial weights of all
		// the assets (i.e., their lower bounds) is numerically equal to one
		//
		// In this case, all assets are LOW plus the linear program is degenerate
		//
		//
		// - The loop above has not prematurely stopped, which implies - because the linear
		// program is feasible - that the sum of the upper bounds of all the assets is numerically
		// equal to one
		//
		// In this case, all assets are UP, plus the linear program is degenerate
		//
		// (In both cases above, there is no real issue as the efficient frontier is then made
		// of only one portfolio, the E-maximizing portfolio, c.f. paragraph 12.3.1 of the 
		// second reference, so that the critical line algorithm will not be started.)
		//
		//
		// - The loop above has stopped on an asset because this asset lies strictly between 
		// its lower bound and its upper bound
		//
		// In this case, this asset is IN, all the assets with a return strictly higher than 
		// this asset are UP, and all the assets with a return strictly lower than this asset
		// are LOW, plus the linear program is is not generate
		//
		//
		// - The loop above has stopped on an asset because the sum of the weights of 
		// all the assets is numerically equal to one
		//
		// In this case, all the assets with a return strictly higher than this asset are UP,
		// this asset is UP, and all the assets with a return strictly lower than this asset
		// are LOW, plus the linear program is degenerate
		//
		// To circumvene the degeneracy, this asset is forced to IN thanks to a numerical 
		// perturbation of its upper bound (its weight is already strictly greater than its
		// lower bound, otherwise, the loop would have stopped on the previous asset)
		if (idx_i == -1 || i == nbAssets) {
			// First two cases above, nothing to do
		}
		else {
			// Last two cases above
			if (variablesStatusManager.isIn(idx_i)) {
				// Nothing to do
			}
			else {
				// Force the asset on which the loop has stopped 
				// (i.e., the last UP asset) to IN.
				variablesStatusManager.setIn(idx_i);
				
				// In order to justify the asset change from UP to IN,
				// numerically alter the asset upper bound.
				upperBounds.setValue(idx_i, 1, 
									 upperBounds.getValue(idx_i, 1) + 2*eps);
			}
		}

		// Return the computed portfolio weights
		return maxReturnWeights;
	}
	
	
	// ------
	
	
	// TODO: Checks, if enabled

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	var maxIterations = opt.maxIter || 1000;
	var lowerBounds = opt.constraints.minWeights;
	var upperBounds = opt.constraints.maxWeights;
	
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	

	// ------
	
	// Initializations	
	var nbAssets = sigma.nbColumns; // the number of assets in the universe
	
	// The only equality constraint supported by the algorithm below
	// is that the weights of the assets must sum to one, but variables
	// below are kept generic.
	var A = Matrix_.ones(1, nbAssets); // the matrix holding the equality constraints
	var nbEqualityConstraints = A.nbRows; // the number of equality constraints 
	var b = Matrix_.ones(1, 1); // the vector holding the right member of the equality constraints
	
	var lb = lowerBounds ? new Matrix_(lowerBounds) : Matrix_.zeros(nbAssets, 1);
	var ub = upperBounds ? new Matrix_(upperBounds) : Matrix_.ones(nbAssets, 1);
	
	var cornerPortfoliosWeights = [];
	var currentCornerPortfolioWeights = null;
	
	var variablesStatusManager = new variablesStatusManager_(nbAssets, nbEqualityConstraints);
	
	// ----	
	
	
	// Step 1: compute the rightmost corner portfolio, corresponding to the E-maximizing 
	// portfolio (i.e., the portfolio achieving the maximum return), c.f. chapter 8 of the
	// first reference and paragraph 12.3.1 of the second reference.
	//
	// To be noted that if there is more than one E-maximizing portfolio, the critical line 
	// algorithm requires a specific E-maximizing portfolio to be computed, c.f. chapter 8 
	// of the first reference.
	//
	// As such a computation is not supported by the algorithm below, the efficient frontier
	// computation is limited to the case when all assets have different returns, which is a 
	// sufficient condition to guarantee the unicity of the E-maximizing portfolio.
	//
	// A practical workaround to this issue, suggested in chapter 9 of the first reference, 
	// is to slightly alter the assets returns and to relaunch the algorithm.
	currentCornerPortfolioWeights = computeMaxReturnPortfolio_(mu, lb, ub);
	var Ai = Matrix_.ones(1, 1);

	
	// Step 1 bis: eliminate degenerate cases when the whole efficient frontier 
	// consists of only one portfolio (i.e., sum l_i = 1 or sum u_i = 1).
	if (Math.abs(1 - lb.sum()) <= eps || Math.abs(1 - ub.sum()) <= eps) {
		cornerPortfoliosWeights.push([currentCornerPortfolioWeights, 0]);
		return cornerPortfoliosWeights;
	}

	
	// Step 2: Initialization of the critical line algorithm, c.f. chapter 13 of the 
	// first reference, paragraph "The Critical Line Algorithm, Setting Up for the CLA",
	// and chapter 13 of the first reference, paragraph 
	// "The Critical Line Algorithm, Appendix A, Module CLA, <C1>-<C6>".	

	// Get the new IN/OUT variables indexes
	//
	// At this stage, only assets variables are set
	var assetsInIdx = variablesStatusManager.getInIndexes();
	var assetsOutIdx = variablesStatusManager.getOutIndexes();


	// Initialize of the xi vector
	var xi = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);

	
	// Initialize the OUT elements of alpha and beta vectors,
	// c.f. formula 13.16 of the first reference:
	// - alpha(out) = X(out)
	// - beta(out) = 0	
	var alpha = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
	for (var i = 1; i <= assetsOutIdx.length; ++i) {
		var out_idx_i = assetsOutIdx[i-1];
		
		alpha.setValue(out_idx_i, 1, 
					   currentCornerPortfolioWeights.getValue(out_idx_i, 1));
	}
	var beta = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);

		
	// Construct the matrix M, with M = [[Sigma A^t], [A 0]],
	// c.f. formula 13.8 of the first reference.
	//
	// Note: the memory used for allocating matrix M is suboptimal,
	// as the matrices Sigma and A are already allocated.
	//
	// This "problem" could be solved through defining M as a matrix
	// defined through a function.
	var M = Matrix_.fill(nbAssets + nbEqualityConstraints, nbAssets + nbEqualityConstraints, 
								function(i,j) { 
									if (i <= nbAssets && j <= nbAssets) {
										return sigma.data[(i-1)*sigma.nbColumns + (j-1)]; // Sigma(i, j)
									}
									else if (i >= nbAssets + 1 && j <= nbAssets) {
										return A.data[(i-nbAssets-1)*A.nbColumns + (j-1)]; // A(i-nbAssets, j)
									}
									else if (i <= nbAssets && j >= nbAssets + 1) {
										return A.data[(j-nbAssets-1)*A.nbColumns + (i-1)]; // A(j-nbAssets, i) == A(i-nbAssets, j)^t
									}
									else {
										return 0;
									}
								});	


	// Construct the Mi matrix, c.f. formula 13.19 of the first reference,
	// Mi = [0 Ai], [Ai^t -Ai^t * Sigma(IN, IN) * Ai]].
	//
	// Because the only equality constraint supported by the algorithm below is
	// that the sum of the assets weights must be equal to 1, there is only one
	// asset IN at this step, and the matrices A_in and Ai of the first reference
	// are then both equal to (1).
	//
	// To be noted, though, that the code below is generic, so that A_in and Ai 
	// are supposed to be nbEqualityConstraints by nbEqualityConstraints matrices since
	// there is nbEqualityConstraints assets IN at this step.
	//
	// As a consequence of this genericity, and for ease of subsequent computations, 
	// a full matrix is allocated for storing Mi instead of a 
	// (nbEqualityConstraints+1) by (nbEqualityConstraints+1) matrix.
	var Mi = Matrix_.zeros(nbAssets + nbEqualityConstraints, nbAssets + nbEqualityConstraints);
		
	// Copy the matrix Ai in the upper right portion of the matrix Mi
	// Copy the matrix Ai^t into the lower left portion of the matrix Mi
	// Copy the matrix -Ai^t * Sigma(IN, IN) * Ai into the lower right portion of the matrix Mi
	var T = Matrix_.atxy(-1, Ai, Matrix_.xy(sigma.submatrix(assetsInIdx, assetsInIdx), Ai));
	for (var j = 1; j <= assetsInIdx.length; ++j) {
		for (var k = 1; k <= assetsInIdx.length; ++k) {
			var var_in_idx_j = assetsInIdx[j-1];
			
			var Ai_j_k = Ai.getValue(j, k);
			Mi.setValue(nbAssets + k, var_in_idx_j, 
						Ai_j_k);		
			Mi.setValue(var_in_idx_j, nbAssets + k, 
						Ai_j_k);
			
			Mi.setValue(nbAssets + j, nbAssets + k, 
						T.getValue(j, k)); 
		}
	}

		
	// Add the lambda variables to the IN set
	variablesStatusManager.setLambdasIn();

	
	// Construct the b_bar vector, c.f. formula 13.14 of the first reference,
	// b_bar(in) = [0 b]^t - M(in, out) * X(out)
	
	// Construct the b_bar vector for the assets variables IN
	var b_bar = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
	for (var i = 1; i <= assetsInIdx.length; ++i) {
		var in_idx_i = assetsInIdx[i-1];
		
		// Initialization of b_bar(idx_in(i)) with [0 b]^t(idx_in(i))
		var b_bar_in_idx_i = 0;

		// Computation of the remaining part of b_bar(idx_in(i))
		for (var j = 1; j <= assetsOutIdx.length; ++j) {
			var out_idx_j = assetsOutIdx[j-1];
			
			b_bar_in_idx_i -= M.getValue(in_idx_i, out_idx_j) * currentCornerPortfolioWeights.getValue(out_idx_j, 1);
		}
		b_bar.setValue(in_idx_i, 1, 
					   b_bar_in_idx_i);
	}
	
	// Construct the b_bar vector for the lambda variables IN (with indexes from
	// nbAssets + 1 to nbAssets + nbEqualityConstraints), which have been
	// added to the IN set just before the b_bar vector computation step.
	for (var i = nbAssets + 1; i <= nbAssets + nbEqualityConstraints; ++i) {
		var in_idx_i = i;
		
		// Initialization of b_bar(idx_in(i)) with [0 b]^t(idx_in(i))
		var b_bar_in_idx_i = b.getValue(in_idx_i-nbAssets, 1);
		
		// Computation of the remaining part of b_bar(idx_in(i))
		for (var j = 1; j <= assetsOutIdx.length; ++j) {
			var out_idx_j = assetsOutIdx[j-1];
			
			b_bar_in_idx_i -= M.getValue(in_idx_i, out_idx_j) * currentCornerPortfolioWeights.getValue(out_idx_j, 1);
		}
		b_bar.setValue(in_idx_i, 1, 
					   b_bar_in_idx_i);
	}	
	
	
	// Step 3: Main loop of the critical line algorithm, c.f. chapter 13 of the 
	// first reference, paragraph "The Critical Line Algorithm, CLA Iteration",
	// and chapter 13 of the first reference, paragraph 
	// "The Critical Line Algorithm, Appendix A, Module CLA, <C10>-<C14>".

	// In each iteration (excepted the first one), the asset that was determined 
	// by the previous iteration to become IN or to become OUT is done so.
	//
	// The different lambdas (lambda_out and lambda_in) are then updated in order 
	// to compute the value of lambda_e corresponding to the next corner portfolio.
	//
	// Once lambda_e is known, the weights of the next corner portfolio can be
	// computed, and the process continues until the value of lambda_e becomes
	// null or negative.
	var iter = 0;
	var lambda_e = 0;	
	var idx_out = -1;
	var lambda_out = 0;
	var status_out = variablesStatusManager.STATUS_UNDEF;
	var idx_in = -1;
	var lambda_in = 0;
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		
		// Update the number of iterations
		++iter;
		
		
		// In case this iteration is not the first one, set the new status of the assets
		// determined by the previous iteration.
		if (iter >= 2) {
			if (lambda_out >= lambda_in) { // an asset IN goes OUT
				// Update the vectors alpha and beta for the asset idx_out going OUT,
				// c.f. formula 13.16 of the first reference:
				// - alpha(idx_out) = X(idx_out)
				// - beta(idx_out) = 0
				alpha.setValue(idx_out, 1, 
							   currentCornerPortfolioWeights.getValue(idx_out, 1));
				beta.setValue(idx_out, 1, 
							  0);
				
				
				// Set the asset idx_out to OUT, with the proper LOW or UP status
				if (status_out == variablesStatusManager.STATUS_LOW) {
					variablesStatusManager.setOnLowerBound(idx_out);
				}
				else {
					variablesStatusManager.setOnUpperBound(idx_out);
				}

				
				// Get the new IN variables indexes
				var variablesInIdx = variablesStatusManager.getInIndexes();
				
				
				// Update the matrix Mi for the asset idx_out going OUT,
				// c.f. formula 13.20 of the reference, reversed:
				// Mi(NEW_IN,NEW_IN) -= Mi(NEW_IN, idx_out) * Mi(idx_out, NEW_IN) / Mi(idx_out, idx_out), with NEW_IN = IN \ {idx_out}				
				var Mi_out_idx_out_idx = Mi.getValue(idx_out, idx_out);
				if (Math.abs(Mi_out_idx_out_idx) <= eps) {
					throw new Error('division by a numerical zero detected, the covariance matrix might not be positive semi-definite');
				}
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					for (var j = 1; j <= variablesInIdx.length; ++j) {
						var in_idx_j = variablesInIdx[j-1];
						
						Mi.setValue(in_idx_i, in_idx_j, 
								    Mi.getValue(in_idx_i, in_idx_j) - Mi.getValue(in_idx_i, idx_out) * Mi.getValue(idx_out, in_idx_j) / Mi_out_idx_out_idx);
					}
					
				}
				
				
				// Update the b_bar vector, c.f. formula 13.22 of the 
				// first reference, reversed:
				// - b_bar(NEW_IN) -= M(NEW_IN, idx_out) * X(idx_out), with NEW_IN = IN \ {idx_out}
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					b_bar.setValue(in_idx_i, 1, 
								   b_bar.getValue(in_idx_i, 1) - M.getValue(in_idx_i, idx_out) * currentCornerPortfolioWeights.getValue(idx_out, 1));
				}
			}
			else { // an asset OUT goes IN				
				// Get the new IN variables indexes
				var variablesInIdx = variablesStatusManager.getInIndexes();

				
				// Update the matrix Mi for the asset idx_in going IN,
				// c.f. formula 13.20 of the first reference:
				// - xi = Mi(IN,IN) * M(IN, idx_in)
				// - xi_j = M(idx_in, idx_in) - <M(IN, idx_in)/xi>
				//
				// - Mi(IN, IN) += (xi * xi^t)(IN, IN) / xi_j
				// - Mi(idx_in, IN) = Mi(IN, idx_in) = -xi(IN) / xi_j
				// - Mi(idx_in, idx_in) = 1 / xi_j
								
				// Compute the vector xi
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					var xi_in_idx_i = 0;
					for (var j = 1; j <= variablesInIdx.length; ++j) {
						var in_idx_j = variablesInIdx[j-1];
						
						xi_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * M.getValue(in_idx_j, idx_in);
					}	
					xi.setValue(in_idx_i, 1, 
								xi_in_idx_i);
				}
				
				// Compute the scalar xi_j
				var xi_j = M.getValue(idx_in, idx_in);
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					xi_j -= M.getValue(idx_in, in_idx_i) * xi.getValue(in_idx_i, 1);
				}
				if (Math.abs(xi_j) <= eps) {
					throw new Error('division by a numerical zero detected, the covariance matrix might not be positive semi-definite');
				}
				
				// Update the matrix Mi
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					for (var j = 1; j <= variablesInIdx.length; ++j) {
						var in_idx_j = variablesInIdx[j-1];
						
						Mi.setValue(in_idx_i, in_idx_j, 
									Mi.getValue(in_idx_i, in_idx_j) + xi.getValue(in_idx_i, 1) * xi.getValue(in_idx_j, 1) / xi_j);
					}
					Mi.setValue(in_idx_i, idx_in, 
								-xi.getValue(in_idx_i, 1)/xi_j);
					Mi.setValue(idx_in, in_idx_i, 
								-xi.getValue(in_idx_i, 1)/xi_j);
				}
				Mi.setValue(idx_in, idx_in, 
							1/xi_j);
				
				
				// Update the b_bar vector, c.f. formulas 13.21 and 13.22 of the 
				// first reference:
				// - b_bar(IN) += M(IN, idx_in) * X(idx_in)
				// - b_bar(idx_in) = -M(idx_in, NEW_OUT) * X(NEW_OUT), with NEW_OUT = OUT \ {idx_in}
				
				// Update the b_bar vector for the current IN variables
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					b_bar.setValue(in_idx_i, 1, 
								   b_bar.getValue(in_idx_i, 1) + M.getValue(in_idx_i, idx_in) * currentCornerPortfolioWeights.getValue(idx_in, 1));
				}
								
				// Set the asset idx_in as IN
				variablesStatusManager.setIn(idx_in);
				
				// Get the new OUT variables indexes, which consists of only assets per construction
				var assetsOutIdx = variablesStatusManager.getOutIndexes();
				
				// Update the b_bar vector for the new IN asset
				var b_bar_in_idx_i = 0;
				for (var i = 1; i <= assetsOutIdx.length; ++i) {
					var out_idx_i = assetsOutIdx[i-1];
					
					b_bar_in_idx_i -= M.getValue(idx_in, out_idx_i) * currentCornerPortfolioWeights.getValue(out_idx_i, 1);
				}
				b_bar.setValue(idx_in, 1, 
							   b_bar_in_idx_i);
				
			}
		}
		
		// Get the new IN/OUT variables indexes
		var variablesInIdx = variablesStatusManager.getInIndexes();
		var assetsOutIdx = variablesStatusManager.getOutIndexes(); // only assets indexes per construction
		
		
		// Determine the next asset IN to be set OUT

		// Update the alpha vector, c.f. formula 13.15 of the first reference:
		// - alpha(IN) = Mi(IN,IN) * b_bar(IN)
		//
		// Update the beta vector, c.f. formula 13.16 of the first reference:
		// - beta(IN) = Mi(IN,IN) * [mu(IN) 0]^t
		//
		// Compute lambda_out, c.f. formula 13.17 of the first reference:
		// - lambda_out = max( (L(i) - alpha(i))/beta(i), beta(i) > 0, i in IN; (U(i) - alpha(i))/beta(i), beta(i) < 0, i in IN)
		idx_out = -1;
		lambda_out = 0;
		status_out = variablesStatusManager.STATUS_UNDEF;
		for (var i = 1; i <= variablesInIdx.length; ++i) {
			var in_idx_i = variablesInIdx[i-1];

			// For all variables IN, compute alpha(idx_in(i)) and beta(idx_in(i))
			var alpha_in_idx_i = 0;
			var beta_in_idx_i = 0;
			for (var j = 1; j <= variablesInIdx.length; ++j) {
				var in_idx_j = variablesInIdx[j-1];
				
				var Mi_in_idx_i_in_idx_j = Mi.getValue(in_idx_i, in_idx_j);
				
				alpha_in_idx_i += Mi_in_idx_i_in_idx_j * b_bar.getValue(in_idx_j, 1);
			
				if (in_idx_j <= nbAssets) {
					beta_in_idx_i += Mi_in_idx_i_in_idx_j * mu.getValue(in_idx_j, 1);
				}
			}
			alpha.setValue(in_idx_i, 1, 
						   alpha_in_idx_i);
			beta.setValue(in_idx_i, 1, 
						  beta_in_idx_i);
			
			
			// For assets variables IN, proceed with the formula 13.17
			if (variablesStatusManager.isAsset(in_idx_i)) {
				// Check for asset reaching the lower limit lb
				if (beta_in_idx_i > eps) {
					var lb_idx_in_i = lb.getValue(in_idx_i, 1);
					
					var tmp_lambda_out = (lb_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
					if (tmp_lambda_out >= lambda_out) {
						idx_out = in_idx_i;
						lambda_out = tmp_lambda_out;
						status_out = variablesStatusManager.STATUS_LOW;
					}
				}
				
				// Check for asset reaching the upper limit ub
				else if (beta_in_idx_i < -eps) {
					var ub_idx_in_i = ub.getValue(in_idx_i, 1);
					
					var tmp_lambda_out = (ub_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
					if (tmp_lambda_out >= lambda_out) {
						idx_out = in_idx_i;
						lambda_out = tmp_lambda_out;
						status_out = variablesStatusManager.STATUS_UP;
					}
				}
			}
		}


		// Determine the next asset OUT to be set IN
		
		// Compute the gamma and delta vectors, c.f. formula 7.10b of the first reference:
		// - gamma(OUT) = [C A^t](OUT, ALL) * alpha, gamma(IN) = 0
		// - delta(OUT) = [C A^t](OUT, ALL) * beta - mu(OUT), delta(IN) = 0
		//
		// In parallel, compute lambda_in, c.f. formula 13.18 of the first reference:
		// - lambda_in = max( -gamma(i)/delta(i), delta(i) > 0, i in LO; -gamma(i)/delta(i), delta(i) < 0, i in UP)
		idx_in = -1;
		lambda_in = 0;
		for (var i = 1; i <= assetsOutIdx.length; ++i) {
			var out_idx_i = assetsOutIdx[i-1];
			
			// Compute gamma(idx_out(i)) and delta(idx_out(i))
			var gamma_out_idx_i = 0;
			var delta_out_idx_i = -mu.getValue(out_idx_i, 1);
			for (var j = 1; j <= nbAssets + nbEqualityConstraints; ++j) {
				var M_out_idx_i_j = M.getValue(out_idx_i, j);
				
				gamma_out_idx_i += M_out_idx_i_j * alpha.getValue(j, 1);
				delta_out_idx_i += M_out_idx_i_j * beta.getValue(j, 1);
			}
			
			// Check for eta_i reaching zero
			// Check for asset coming off lower limit
			if (variablesStatusManager.isOnLowerBound(out_idx_i)) {
				if (delta_out_idx_i > eps) {
					var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
				
					if (tmp_lambda_in >= lambda_in) {
						idx_in = out_idx_i;
						lambda_in = tmp_lambda_in;
					}
				}
			}
			// Check for asset coming off upper limit
			else {
				if (delta_out_idx_i < -eps) {
					var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
					
					if (tmp_lambda_in >= lambda_in) {
						idx_in = out_idx_i;
						lambda_in = tmp_lambda_in;
					}
				}
			}
		}

		
		// The value of lambda_e for the next corner portfolio is the maximum of
		// the two lambda_out and lambda_in computed above.
		//
		// In case lambda_e == lambda_out, it means an asset first goes OUT as lambda_e
		// is decreased; otherwise, in case lambda_e == lambda_in, it means an asset
		// first goes IN as lambda_e is decreased.
		lambda_e = Math.max(lambda_out, lambda_in, 0);
		
		
		// Compute the weights of the next corner portfolio
		for (var i = 1; i <= variablesInIdx.length; ++i) {
			var in_idx_i = variablesInIdx[i-1];
			
			// In case the variable IN is an asset variable, update the current corner portfolio
			if (variablesStatusManager.isAsset(in_idx_i)) {
				currentCornerPortfolioWeights.setValue(in_idx_i, 1, 
													   alpha.getValue(in_idx_i, 1) + lambda_e * beta.getValue(in_idx_i, 1));
			}
		}
		
		// Save the current corner portfolio
		cornerPortfoliosWeights.push([new Matrix_(currentCornerPortfolioWeights), lambda_e]);			

		
		// When the value of lambda_e becomes numerically null or negative, the critical
		// line algorithm can be stopped.
		if (lambda_e < eps) {
			break;
		}
	}
	
	
	// Return the computed efficient frontier array, filtered for numerically identical
	// corner portfolios.
	var finalCornerPortfoliosWeights = new Array(cornerPortfoliosWeights.length);
	
	// First, the E-maximizing portfolio is always included on the efficient frontier
	var idx = 0;
	var cornerPortfolio = cornerPortfoliosWeights[idx][0];
	var lambda = cornerPortfoliosWeights[idx][1];
	finalCornerPortfoliosWeights[idx] = [cornerPortfolio, lambda]; 
	
	// Then, for each computed corner portfolio:
	// - If it is numerically identical to the last corner portfolio included in the
	// output efficient frontier, replace this last portfolio with it
	//
	// - Otherwise, add it
	for (var i = 1; i < cornerPortfoliosWeights.length; ++i) {
		var cornerPortfolio_i = cornerPortfoliosWeights[i][0];
		var lambda_i = cornerPortfoliosWeights[i][1];
		
		if (Matrix_.areEqual(cornerPortfolio_i, cornerPortfolio, eps)) {
			finalCornerPortfoliosWeights[idx] = [cornerPortfolio_i, lambda_i];
		}
		else {
			cornerPortfolio = cornerPortfolio_i;
			lambda = lambda_i;
			
			++idx;
			finalCornerPortfoliosWeights[idx] = [cornerPortfolio, lambda];
		}
	}
	
	// Resize the output efficient frontier array as required
	finalCornerPortfoliosWeights.length = idx + 1;
	
	// Return the final efficient frontier array
	return finalCornerPortfoliosWeights;
}
