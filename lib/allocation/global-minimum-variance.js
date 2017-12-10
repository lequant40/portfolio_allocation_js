/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function globalMinimumVarianceWeights
*
* @summary Compute the weights of the global minimum variance portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only
* global minimum variance portfolio of n assets, defined as the weights which minimizes the variance of the portfolio.
*
* This portfolio lies on the Markowitz efficient frontier.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* The algorithm used is a cyclical coordinate descent, c.f. the reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the global minimum variance portfolio, array of n real numbers.
*
* @example
* globalMinimumVarianceWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // XX
*/
self.globalMinimumVarianceWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1).normalize();
	
	// Preparational computations
	var sigma_x = Matrix_.product(sigma, x); // SIGMA*x
	var var_x = Matrix_.vectorDotProduct(sigma_x, x); // sigma(x), the portfolio variance
	var obj = 0.5 * var_x - x.sum();
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma
	var iter = 0;
	var converged = false;
	while (!converged) {
        // By default, let's assume the new iteration will lead to the convergence of the algorithm
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {		
			// Define the coefficients of the second order polynomial ax_i^2 + b_ix + c, c.f. the reference
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i) - 1; // (SIGMA*x)_i - x_i*sigma_i^2 - 1, might be any sign
    	    var c = 0;
    	    
    	    // Note: what follows is not detailled in the reference, as lambda_erc is supposed there to be non zero.
    	    //
    	    // The equation ax_i^2 + bx_i + c = 0 has two roots: 0 and another one, possibly strictly negative, 0 or strictly positive:
			// Case #1 - If the other root is 0 or strictly negative, x_i^* is then equal to 0
			// Case #2 - If the other root is strictly positive, x_i^* is equal to the value minimizing the objective function in 1D

			// Extract the "interesting" root of the equation ax_i^2 + bx_i = 0, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
    	    var q = -(b_p + sign_b_p * Math.abs(b_p));
    	    var r1 = q/a;		
			
			// Case #1 always needs to be tested
			var xi_star_1 = 0;
			
			  // Case #1 weights update
			x.setValueAt(i, 1, xi_star_1);

    	      // Case #1 updated SIGMA*x and x'*SIGMA*x products
    	    sigma_x = Matrix_.product(sigma, x)
    	    var_x = Matrix_.vectorDotProduct(sigma_x, x);
			
			  // Case #1 objective function computation
			var obj_star = 0.5 * var_x - x.sum();
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				var sigma_x_2 = Matrix_.product(sigma, x)
				var var_x_2 = Matrix_.vectorDotProduct(sigma_x_2, x);
				
				// Case #2 objective function computation
				var obj_star_2 = 0.5 * var_x_2 - x.sum();
				
				// Selection between the two cases
				// Note: Portfolio sparsity is priviledged in case of equality
				if (obj_star_2 < obj_star) {
					// Case #2 weights update is already done
					
					// Update the products for convergence condition evaluation + next loop evaluation
					sigma_x = sigma_x_2;
					var_x = var_x_2;
					obj_star = obj_star_2;
				}
				else {
					// Case #1 weights re-update, as erased by case #2 test
					x.setValueAt(i, 1, xi_star_1);
					
					// No update on the products, as case #1 is the default case
				}
			}
        }
		
		// Update the convergence condition: |obj - obj*| <= eps
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
	
	// Return the computed weights, after normalization
	x = x.normalize();
	return x.toArray();
}

