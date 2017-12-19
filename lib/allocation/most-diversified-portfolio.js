/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function mostDiversifiedWeights
*
* @summary Compute the weights of the most diversified portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only
* most diversified portfolio of n assets, defined as the weights which maximizes the diversification ratio of the portfolio.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio of each asset is the same.
*
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/JPM.2008.35.1.40">Toward Maximum Diversification by Y. Choueifaty, Y. Coignard, The Journal of Portfolio Management, Fall 2008, Vol. 35, No. 1: pp. 40-51</a>
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the most diversified portfolio, array of n real numbers.
*
* @example
* mostDiversifiedWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // [0.33, 0.67]
*/
self.mostDiversifiedWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Extract the standard deviations of the assets in vector format
	var stddevs = sigma.diagonal().elemMap(function(i,j,val) { return Math.sqrt(val); });

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1);
	x = x.normalize(x);
	
	// Preparational computations
	var sigma_x = Matrix_.product(sigma, x); // SIGMA*x
	var sigma_x_2 = Matrix_.product(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x)); // sigma(x)
	var dr = Matrix_.vectorDotProduct(stddevs, x) / s_x; // diversification ratio DR
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma
	var iter = 0;
	var converged = false;
	while (!converged) {
        // By default, let's assume the new iteration will lead to the convergence of the algorithm
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {
    	    // Define the coefficients of the second order polynomial ax_i^2 + b_ix + c, c.f. the second reference
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i) - stddevs.getValueAt(i,1); // (SIGMA*x)_i - x_i*sigma_i^2 - sigma_i, might be any sign
    	    var c = 0;
    	    
    	    // Note: what follows is not detailled in the second reference, as lambda_erc is supposed there to be non zero.
    	    //
    	    // The equation ax_i^2 + bx_i + c = 0 has two roots: 0 and another one, possibly strictly negative, 0 or strictly positive:
			// Case #1 - If the other root is 0 or strictly negative, x_i^* is then equal to 0
			// Case #2 - If the other root is strictly positive, x_i^* is equal to the value maximizing the diversification ratio in 1D

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
    	    sigma_x = Matrix_.product(sigma, x, sigma_x)
    	    s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x));
			
			  // Case #1 diversification ratio computation
			var dr_star = Matrix_.vectorDotProduct(stddevs, x) / s_x; // Weighted average volatility of the portfolio divided by the volatility of the portfolio
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				sigma_x_2 = Matrix_.product(sigma, x, sigma_x_2);
				var s_x_2 = Math.sqrt(Matrix_.vectorDotProduct(sigma_x_2, x));
				
				// Case #2 diversification ratio computation
				var dr_star_2 = Matrix_.vectorDotProduct(stddevs, x) / s_x_2;
				
				// Selection between the two cases
				// Note: Portfolio sparsity is priviledged in case of equality
				if (dr_star_2 > dr_star) {
					// Case #2 weights update is already done
					
					// Update the products for convergence condition evaluation + next loop evaluation
					sigma_x = Matrix_.copy(sigma_x_2, sigma_x);
					s_x = s_x_2;
					dr_star = dr_star_2;
				}
				else {
					// Case #1 weights re-update, as erased by case #2 test
					x.setValueAt(i, 1, xi_star_1);
					
					// No update on the products, as case #1 is the default case
				}
			}
        }
		
        // Update the convergence condition: |dr - dr*| <= eps
		if (Math.abs(dr_star - dr) > eps) {
			converged = false;
			dr = dr_star;
		}
		
		// Update the number of iterations
		++iter;
		
		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}
	
	// Return the computed weights, after normalization
	x = x.normalize(x);
	return x.toArray();
}

