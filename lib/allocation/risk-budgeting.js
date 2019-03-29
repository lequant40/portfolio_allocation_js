/**
 * @file Functions related to risk budgeting portfolio.
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
* of n assets with risk budgeting contraints.
*
* This portfolio has the property that the total contribution 
* of each asset to the risk of the portfolio is equal to a pre-determined budget weight.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* To be noted that the algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is semi-definite positive.
*
* @see <a href="https://ssrn.com/abstract=2009778">Bruder, Benjamin and Roncalli, Thierry, Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012).</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {Array.<number>} rb the risk budgets, array of n real strictly positive numbers summing to one.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.outputPortfolioVolatility a boolean indicating whether the portfolio volatility should be provided in output
* (if set to true) or not (if set to false); defaults to false.
* @return {Array.<number>|Array.<Array.<number>>} if opt.outputPortfolioVolatility is set to false, the weights corresponding to the risk budgeting portfolio, 
* array of n real numbers, and if opt.outputPortfolioVolatility is set to true, an array arr of two elements:
* - arr[0], the weights corresponding to the risk budgeting portfolio, array of n real numbers
* - arr[1], the volatility of the computed risk budgeting portfolio, a real number
*
* @example
* riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75], {eps: 1e-10, maxIter: 10000});
* // [~0.45, ~0.55]
*/
self.riskBudgetingWeights = function (sigma, rb, opt) {
	// The numerical zero 
	var eps_tol = 1e-12;
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	var outputPortfolioVolatility = false || opt.outputPortfolioVolatility; 
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Convert rb to vector format
	var rb = new Matrix_(rb);

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that rb contains strictly positive numbers summing to one
	// Check that sigma and rb are rows compatible


	// ------

	var nbAssets = sigma.nbRows;
	
	
	// ------
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
	
	// Preparational computations
	var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x)); // sigma(x)
	var obj = 0.5 * s_x - Matrix_.vectorDotProduct(rb, x.elemMap(function(i,j,val) { return Math.log(val);})); // the objective function to be minimized, c.f. formula 3 of the second reference
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma and b
	var iter = 0;
	var converged = false;
	while (!converged) {
        // Convergence condition is false if any of the coordinate-wise convergence condition is false
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {
			// Save the old asset weight i before any update
			var xi_old = x.data[i-1];

			
			// Define the coefficients of the second order polynomial ax_i^2 + b_ix + c_i, c.f. the second reference
    	    var a = sigma.data[(i-1)*sigma.nbColumns + (i-1)]; // sigma_i^2, always > 0
    	    var b = sigma_x.data[i-1] - x.data[i-1] * a; // (SIGMA*x)_i - x_i*sigma_i^2, might be any sign
    	    var c = -rb.data[i-1] * s_x; // -b_i * sigma(x), always <= 0 (== 0 iff the covariance matrix is semi-definite positive and sigma(x) == 0)

			
    	    // Extract the strictly positive root x_i^* of the equation ax_i^2 + bx_i + c = 0, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
			    throw new Error('Negative discriminant during iteration ' + iter + ', covariance matrix might not be semi-definite positive');
			}
    	    var q = -(b_p + sign_b_p * Math.sqrt(disc));
    	    var r1 = q/a;
    	    var r2 = c/q;
    	    var xi_star = r1 > 0 ? r1 : r2;
    	    
			
			// Update the asset weight i
			x.data[i-1] = xi_star;
    	    
			
    	    // Compute the updated SIGMA*x and x'*SIGMA*x elements for convergence condition evaluation,
			// and next loop evaluation.
			//
			// The update of the vector SIGMA*x uses the efficient update procedure described 
			// in the second reference, based on the fact that only one coordinate of the vector x
			// changes per iteration.
			//
			// To be noted that the update of the value x'*SIGMA*x does not use this procedure, 
			// because it requires a dot product, which is then equivalent to the full recomputation
			// of the volatility from SIGMA*x.
			
			// Compute the updated SIGMA*x
			for (var j = 1; j <= nbAssets; ++j) {
				sigma_x.data[j-1] += sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_star;
				sigma_x.data[j-1] -= sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_old;
			}
			
			// Compute the updated x'*SIGMA*x
			var sigma_x_x = Matrix_.vectorDotProduct(sigma_x, x);		
			
			// In case the variance is numerically negative zero, which can occur with 
			// a semi-positive definite covariance matrice, it is replaced with zero
			if (sigma_x_x < 0) {
				if (-sigma_x_x < eps_tol) {
					sigma_x_x = 0;
				}
				else {
					throw new Error('Negative volatility during iteration ' + iter + ', covariance matrix might not be semi-definite positive');
				}
			}
			
			// Compute the updated volatility SQRT(x'*SIGMA*x)
			s_x = Math.sqrt(sigma_x_x);

			
    	    // Update the convergence condition: |RC_i* - b_i| <= eps, i = 1..nbAssets,
			// used in case the risk contributions are properly defined (i.e., a strictly
			// positive portfolio volatility).
			if (s_x != 0) {
				var rci_star = x.data[i-1] * sigma_x.data[i-1] / s_x;

				if (Math.abs(rci_star - rb.data[i-1]) > eps) {
					converged = false;
				}
			}
		}

		// Update the generic convergence condition: |obj - obj*| <= eps,
		// used in all cases
		var obj_star = 0.5 * s_x - Matrix_.vectorDotProduct(rb, x.elemMap(function(i,j,val) { return Math.log(val);}));		
		if (Math.abs(obj_star - obj) > eps) {
			converged = false;
			obj = obj_star;
		}


		// Update the number of iterations
		++iter;


		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}

	// Normalize the computed weights, and compute the normalization constant
	var sum_x = x.sum();
	x = x.normalize(x);

	// Depending on what is requested in output, return the computed normalized weights
	// and possibly the associated portfolio volatility.
	if (outputPortfolioVolatility === true) {
		return [x.toArray(), s_x/sum_x];
	}
	else {
		return x.toArray();
	}
}

