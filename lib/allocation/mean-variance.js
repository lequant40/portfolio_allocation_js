/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.computeCornerPortfolios_ = computeCornerPortfolios_;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function meanVarianceOptimizationWeights
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a return
* constraint, to misc. volatility constraints, or to a risk tolerance constraint.
*
* @description This function returns the weights w_1,...,w_n associated to the 
* fully invested long-only mean-variance efficient portfolio of n assets subject to either (exclusive):
* - a return constraint, in which case this portfolio, if it exists,
* has the lowest attainable volatility among all the feasible portfolios satisfying the return constraint
* - a volatility constraint, in which case this portfolio, if it exists,
* has the highest attainable return among all the feasible portfolios satisfying the volatility constraint
* - a risk tolerance constraint, in which case this portfolio, which always exists,
* is associated to a risk tolerance satisfying the risk tolerance constraint (this portfolio is a solution
* to the optimization problem min_w <Vw/w>/2 - rt*<w/e>, s.t. <w/e>=1 and l <= x <= u)
* - a maximum volatility constraint, in which case this portfolio, if it exists,
* has a volatility equal to the maximum feasible volatility lower than or equal to the maximum volatility constraint,
* and has the highest attainable return among all the feasible portfolios with the same volatility

* Optionally, the following constraints can be added:
* - Partial investment constraint, replacing the full investment constraint
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The main algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt parameters for the mean-variance optimization algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
* @param {number} opt.constraints.return, the desired value for the return of the portfolio, a real number.
* @param {number} opt.constraints.volatility, the desired value for the volatility of the portfolio, a positive real number.
* @param {number} opt.constraints.riskTolerance, the desired value for the risk tolerance parameter, a positive real number.
* @param {number} opt.constraints.maxVolatility, the maximum desired value for the volatility of the portfolio, a positive real number.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights corresponding to the mean-variance efficient portfolio, array of n real numbers.
*
* @example
* meanVarianceOptimizationWeights([0.1, 0.2], [[1, 0.3], [0.3, 1]], { constraints: {return: 0.15}})
* // [0.5, 0.5] 
*/
self.meanVarianceOptimizationWeights = function(mu, sigma, opt) {	
	var eps = 1e-8; // the numerical zero
	
	// Internal function to compute the return of a portfolio
	function portfolioReturn(x, mu, sigma) {
		return Matrix_.vectorDotProduct(mu, x);
	}
	
	// Internal function to compute the volatility of a portfolio
	function portfolioVolatility(x, mu, sigma) {	
		// Compute the variance x'*SIGMA*x
		var sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(sigma, x), x);		

		// In case the variance is numerically zero, which can occur with 
		// a semi-positive definite covariance matrix, it is replaced with zero.
		//
		// Otherwise, if the variance is negative, stops the algorithm,
		// since the covariance matrix is then not numerically semi-positive definite.
		if (Math.abs(sigma_x_x) <= eps) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -eps) {
			throw new Error('negative volatility, covariance matrix might not be semi-definite positive');
		}
		
		// Compute the volatility SQRT(x'*SIGMA*x)
		var s_x = Math.sqrt(sigma_x_x);
	   
		// Return the computed volatility
		return s_x;
	}
	
	
	// ------	

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	
	// Decode the full investment constraint
	var fullInvestmentContraint = true;
	if (opt.constraints.fullInvestment !== undefined) {
		fullInvestmentContraint = opt.constraints.fullInvestment;
	}
	
	// Decode the mutually exclusive return/volatility constraints
	var returnConstraint = opt.constraints["return"]; // .return does not work within Google Script
	var volatilityConstraint = opt.constraints.volatility;
	var maxVolatilityConstraint = opt.constraints.maxVolatility;
	var riskToleranceConstraint = opt.constraints.riskTolerance;
	if (returnConstraint === undefined && volatilityConstraint === undefined &&
		maxVolatilityConstraint === undefined && riskToleranceConstraint == undefined) {
		throw new Error('missing return, volatility or risk tolerance constraints');
	}
	if ( (returnConstraint !== undefined && volatilityConstraint !== undefined) ||
		 (returnConstraint !== undefined && maxVolatilityConstraint !== undefined) ) {
		throw new Error('simultaneous return and volatility constraints');
	}
	if ( riskToleranceConstraint !== undefined && 
	     (returnConstraint !== undefined || maxVolatilityConstraint !== undefined || volatilityConstraint !== undefined) ) {
		throw new Error('simultaneous risk tolerance and return or volatility constraints');
	}
	
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	
	
	// ------
	
	
	// Initializations
	var nbAssets = sigma.nbColumns;
	
	
	// Compute the corner portfolios defining the efficient frontier
	//
	// In case of partial investment constraint, the input assets are altered by
	// the addition of a virtual risk free asset with (0,0,0) return/volatility/covariance.
	//
	// The corner portfolios defining the altered efficient frontier are then computed
	// instead of the corner portfolios defining the initial efficient frontier.
	var cornerPortfolios;
	if (fullInvestmentContraint) {
		cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	}
	else {
		// Extend the returns/covariances of the input assets with a risk free asset,
		// of index (nbAssets + 1).
		mu = Matrix_.fill(nbAssets + 1, 1, 
							function(i, j) { 
								if (i <= nbAssets) {
									return mu.getValue(i, 1);
								}
								else {
									return 0;
								}
							});
		sigma = Matrix_.fill(nbAssets + 1, nbAssets + 1, 
							function(i, j) { 
								if (i <= nbAssets && j <= nbAssets) {
									return sigma.getValue(i, j);
								}
								else {
									return 0;
								}
							});
										
		// In case minimum/maximum weights constraints are provided, extend them.
		var minWeights = null;
		if (opt.constraints.minWeights) {
			minWeights = typeof Float64Array === 'function' ? new Float64Array(nbAssets + 1) : new Array(nbAssets + 1); 
			for (var i = 0; i < nbAssets; ++i) {
				minWeights[i] = opt.constraints.minWeights[i];
			}
			minWeights[nbAssets] = 0; // risk free asset, no min weight
		}
		var maxWeights = null;
		if (opt.constraints.maxWeights) {
			maxWeights = typeof Float64Array === 'function' ? new Float64Array(nbAssets + 1) : new Array(nbAssets + 1); 
			for (var i = 0; i < nbAssets; ++i) {
				maxWeights[i] = opt.constraints.maxWeights[i]
			}
			maxWeights[nbAssets] = 1; // risk free asset, no max weight
		}
			
		// Create and fill a new options structure
		var newOpt = { maxIter: opt.maxIter, constraints: opt.constraints };
		if (opt.constraints.minWeights) {
			newOpt.constraints.minWeights = minWeights;
		}
		if (opt.constraints.maxWeights) {
			newOpt.constraints.maxWeights = maxWeights;
		}
			
		// Compute the corner portfolios corresponding to the new efficient frontier
		// for the altered input assets.
		cornerPortfolios = computeCornerPortfolios_(mu, sigma, newOpt);			
	}
	
	
	// Depending on the constraint, proceed with a different algorithm
	// to compute the requested efficient portfolio.
	var efficientPortfolioWeights;
	if (returnConstraint !== undefined) {
		var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, {constraint: "return", 
																						              constraintValue: returnConstraint});
		if (efficientPortfolio.length == 0) {
			throw new Error('no matching efficient portfolio with a return equal to ' + returnConstraint);
		}
		else {
			efficientPortfolioWeights = efficientPortfolio[0];
		}
	}
	else if (volatilityConstraint !== undefined) {
		var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, {constraint: "volatility", 
																						              constraintValue: volatilityConstraint});
		if (efficientPortfolio.length == 0) {
			throw new Error('no matching efficient portfolio with a volatility equal to ' + volatilityConstraint);
		}
		else {
			efficientPortfolioWeights = efficientPortfolio[0];
		}
	}
	else if (riskToleranceConstraint !== undefined) {
		// Compute the portfolio with the highest/lowest risk tolerance on the efficient frontier,
		// which corresponds to the rightmost/leftmost portfolios.
		var highestRiskTolerancePortfolioWeights = cornerPortfolios[0][0];
		var highestRiskTolerance = cornerPortfolios[0][1];
		var lowestRiskTolerancePortfolioWeights = cornerPortfolios[cornerPortfolios.length-1][0];
		var lowestRiskTolerance = cornerPortfolios[cornerPortfolios.length-1][1];

		
		// If the desired risk tolerance is greater than the highest attainable
		// risk tolerance on the efficient frontier, the efficient portfolio is computed
		// as the rightmost portfolio on the efficient frontier, because this portfolio 
		// also corresponds to all risk tolerances greater than or equal to the
		// highest attainable risk tolerance on the efficient frontier.
		//
		// If the desired risk tolerance is lower than the lowest attainable
		// risk tolerance on the efficient frontier, the efficient portfolio is computed
		// as the leftmost portfolio on the efficient frontier, because this is a numerical
		// rounding issue, since the lowest attainable risk tolerance on the efficient frontier must be 
		// equal to 0.
		//
		// Otherwise, the efficient portfolio is computed as the efficient portfolio 
		// with a risk tolerance strictly equal to the desired risk tolerance.
		if (riskToleranceConstraint > highestRiskTolerance - eps) {
			efficientPortfolioWeights = highestRiskTolerancePortfolioWeights;
		}
		else if (riskToleranceConstraint < lowestRiskTolerance + eps) {
			efficientPortfolioWeights = lowestRiskTolerancePortfolioWeights;
		}
		else {
			var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, {constraint: "riskTolerance", 
																										  constraintValue: riskToleranceConstraint});
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no matching efficient portfolio with a risk tolerance equal to ' + riskToleranceConstraint);
			}
			else {
				efficientPortfolioWeights = efficientPortfolio[0];
			}
		}																						  
	}
	else if (maxVolatilityConstraint !== undefined) {
		// Compute the portfolio with the highest volatility on the efficient frontier,
		// which corresponds to the rightmost portfolio.
		var highestVolatilityPortfolioWeights = cornerPortfolios[0][0];
		var highestVolatility = portfolioVolatility(highestVolatilityPortfolioWeights, mu, sigma);

		// If the desired maximum volatility is greater than the highest attainable
		// volatility on the efficient frontier, the efficient portfolio is computed
		// as the rightmost portfolio on the efficient frontier.
		//
		// Otherwise, the efficient portfolio is computed as the efficient portfolio 
		// with a volatility strictly equal to the desired maximum volatility.		
		if (maxVolatilityConstraint > highestVolatility - eps) {
			efficientPortfolioWeights = highestVolatilityPortfolioWeights;
		}
		else {
			var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, {constraint: "volatility", 
																						                  constraintValue: maxVolatilityConstraint});
			if (efficientPortfolio.length == 0) {
				throw new Error('no matching efficient portfolio with a volatility lower than or equal to ' + maxVolatilityConstraint);
			}
			else {
				efficientPortfolioWeights = efficientPortfolio[0];
			}
		}
	}	

	
	// In case of partial investment constraint, transform the fully invested efficient
	// portfolio with altered input assets into the associated partially invested efficient 
	// portfolio with the input assets.
	if (!fullInvestmentContraint) {
		efficientPortfolioWeights = Matrix_.fill(nbAssets, 1, 
												function(i, j) { 
													if (i <= nbAssets) {
														return efficientPortfolioWeights.getValue(i, 1);
													}
												});
	}

	
	// Return the computed portfolio weights
	return efficientPortfolioWeights.toArray();
}


/**
* @function computeMeanVarianceEfficientPortfolio_
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to an exact return
* constraint, to an exact volatility constraint or to an exact risk tolerance constraint.
*
* @description This function returns the weights w_1,...,w_n associated to the 
* fully invested long-only mean-variance efficient portfolio of n assets subject to either (exclusive):
* - a return constraint, in which case this portfolio, if it exists,
* has the lowest attainable volatility among all the feasible portfolios satisfying the return constraint
* - a volatility constraint, in which case this portfolio, if it exists,
* has the highest attainable return among all the feasible portfolios satisfying the volatility constraint
* - a risk tolerance constraint, in which case this portfolio, which always exists,
* is associated to a risk tolerance satisfying the risk tolerance constraint
*
* The algorithm used internally is based on a bisection search on the mean-variance efficient frontier, c.f. the reference.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {Array<Array.<Object>>} cornerPortfolios the list of corner portfolios defining the efficient frontier on which to compute the mean-variance efficient portfolio, as computed
* by the internal computeCornerPortfolios_ function
* @param {object} opt parameters for the algorithm.
* @param {string} opt.constraint, the type of constraint, a string either equal to "return" for a return constraint, to "volatility" for a volatility constraint,
* or to "riskTolerance" for a risk tolerance constraint.
* @param {number} opt.constraintValue, the desired constraint value (i.e., the desired return or volatility), a real number.
*
* @return {Array.<Matrix_, number, number>} the weights corresponding to the mean-variance efficient portfolio, matrix of n real numbers.
*
*/
function computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, opt) {
	var eps = 1e-8; // the numerical zero
	
	// Internal function to compute the return of a corner portfolio
	function portfolioReturn(x, mu, sigma, lambda) {
		return Matrix_.vectorDotProduct(mu, x);
	}
	
	// Internal function to compute the volatility of a corner portfolio
	function portfolioVolatility(x, mu, sigma, lambda) {	
		// Compute the variance x'*SIGMA*x
		var sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(sigma, x), x);		

		// In case the variance is numerically zero, which can occur with 
		// a semi-positive definite covariance matrix, it is replaced with zero.
		//
		// Otherwise, if the variance is negative, stops the algorithm,
		// since the covariance matrix is then not numerically semi-positive definite.
		if (Math.abs(sigma_x_x) <= eps) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -eps) {
			throw new Error('negative volatility, covariance matrix might not be semi-definite positive');
		}
		
		// Compute the volatility SQRT(x'*SIGMA*x)
		var s_x = Math.sqrt(sigma_x_x);
	   
		// Return the computed volatility
		return s_x;
	}
	
	// Internal function to compute the risk tolerance of a corner portfolio
	function portfolioRiskTolerance(x, mu, sigma, lambda) {
		return lambda;
	}
	
	// Internal function to compute the (at most) two corner portfolios enclosing the
	// efficient portfolio with a target constraint value (return, volatility or risk tolerance),
	// using a binary search algorithm.
	//
	// The usage of a binary search algorithm is justified because the corner portfolios
	// return, volatility and risk tolerance are decreasing as soon as there are at least two corner
	// portfolios on the efficient frontier.
	function computeEnclosingCornerPortfolios(fct, fctValue, cornerPortfolios) {		
		// The efficient frontier portfolios are provided from highest return/volatility/risk tolerance
		// to lowest return/volatility/risk tolerance, so that *_min below refers to properties of the portfolio
		// with the lowest return/volatility/risk tolerance.
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];
		
		var lambda_min = cornerPortfolios[idx_min][1];
		var lambda_max = cornerPortfolios[idx_max][1];

		var fct_min = fct(weights_min, mu, sigma, lambda_min);
		var fct_max = fct(weights_max, mu, sigma, lambda_max);

		// If the desired target constraint function value is numerically strictly greater than the highest attainable
		// constraint function value on the efficient frontier, or if the desired target constraint function 
		// value is numerically strictly lower than the lowest attainable constraint function value on the efficient frontier, 
		// there is no enclosing corner portfolio on the efficient frontier.
		if (fctValue > fct_max + eps || fctValue < fct_min - eps) {
			return [];
		}
		
		// Otherwise, if the target constraint function value is numerically reached on one of the
		// two extremal corner portfolios, return immediately.
		if (Math.abs(fctValue - fct_min) <= eps) {
			return [[idx_min, weights_min, fct_min]];
		}
		else if (Math.abs(fctValue - fct_max) <= eps) {
			return [[idx_max, weights_max, fct_max]];
		}
		
		// Otherwise, determine the two adjacent corner portfolios enclosing the portfolio
		// with a constraint function value numerically equals to the target constraint function value, 
		// using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle point
			var idx_middle = Math.floor(idx_max + (idx_min - idx_max)/2); // formula avoiding numerical overflow
			var weights_middle = cornerPortfolios[idx_middle][0];
			var lambda_middle = cornerPortfolios[idx_middle][1];
			var fct_middle = fct(weights_middle, mu, sigma, lambda_middle);
			
			// Determine in which sub-interval ]idx_max, idx_middle[ or ]idx_middle, idx_min[
			// lies the portfolio with the target constraint function value.
			if (fct_middle > fctValue) {
				idx_max = idx_middle;
				fct_max = fct_middle;
				weights_max = weights_middle;
			}
			else if (fct_middle < fctValue) {
				idx_min = idx_middle;
				fct_min = fct_middle;
				weights_min = weights_middle;
			}
			else { // the target constraint function value is exactly attained on the idx_middle-th corner portfolio
				return [[idx_middle, weights_middle, fct_middle]];
			}
		}

		
		// Return the computed adjacent corner portfolios, as well as
		// the associated function values.
		return [[idx_min, weights_min, fct_min], [idx_max, weights_max, fct_max]];
	}
	
	
	// ------	

	
	// Decode the input constraints
	var constraint = opt.constraint;
	if (constraint === undefined || constraint === null) {
		throw new Error('internal error: missing constraint');
	}
	
	var constraintFct;
	if (constraint == "return") {
		constraintFct = portfolioReturn;
	}
	else if (constraint == "volatility") {
		constraintFct = portfolioVolatility;
	}
	else if (constraint == "riskTolerance") {
		constraintFct = portfolioRiskTolerance;
	}
	else {
		throw new Error('internal error: unknown constraint');
	}
	
	var constraintValue = opt.constraintValue;
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('internal error: missing constraint value');
	}
	
	
	// Safety checks on the corner portfolios
	if (cornerPortfolios.length === 0) {
		throw new Error('internal error: no corner portfolios');
	}

	
	// Compute the (at most) two corner portfolios enclosing the efficient 
	// portfolio with a constraint function value equals to the target constraint function value.
	var enclosingCornerPortfolios = computeEnclosingCornerPortfolios(constraintFct, 
	                                                                 constraintValue, 
																	 cornerPortfolios);

	
	// Then:
	// - In case the target constraint function value is not reachable, return an empty portfolio 
	//
	// - In case there is a unique computed corner portfolio with a constraint function value
	// equals to the target constraint function value, return the associated portfolio weights
	//
	// - In case there are two corner portfolios enclosing the efficient portfolio with 
	// a constraint function value equals to the target constraint function value, the weights associated 
	// to this efficient portfolio are a convex combination of the weights of the two computed enclosing
	// corner portfolios (c.f. the reference): w = t*w_min + (1-t)*w_max, t in [0,1], with t now to be determined.
	if (enclosingCornerPortfolios.length == 0) {
		return [];
	}
	else if (enclosingCornerPortfolios.length == 1) {
		var idx_min = enclosingCornerPortfolios[0][0];
		var weights = enclosingCornerPortfolios[0][1];
		
		// Return the computed portfolio weights
		return [weights, idx_min];
	}
	else {
		// Extract information about the computed efficient corner portfolios
		var idx_min = enclosingCornerPortfolios[0][0];
		var weights_min = enclosingCornerPortfolios[0][1];
		var fct_min = enclosingCornerPortfolios[0][2];
			
		var idx_max = enclosingCornerPortfolios[1][0];
		var weights_max = enclosingCornerPortfolios[1][1];
		var fct_max = enclosingCornerPortfolios[1][2];
		
		// Compute t above
		var t;
		
		// If the constraint function value to compute is "return", then, the procedure to 
		// compute t above is the following:
		//
		// E(w) = <mu/w> and by linearity of E, we have
		// E(w) = t*E(w_min) + (1-t)*E(w_max) and E(w) = constraintValue
		// <=>
		// t = (E(w_max) - constraintValue)/(E(w_max) - E(w_min))
		//
		//
		// If the constraint function value to compute is "volatility", then, the procedure to 
		// compute t above is the following:
		//
		// Let the volatility be V(w) = <Sigma*w/w>.
		// Then, by symmetry and bi-linearity of V, V(w) = t^2*V(w_min) + (1-t)^2*V(w_max) + 2*t*(1-t)*<Sigma*w_min/w_max>
		// and V(w) = constraintValue^2
		// <=> t is the solution belonging to ]0,1[ of the second order polynomial equation
		// t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) -2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) - constraintValue^2 = 0
		//
		//
		// If the constraint function value to compute is "riskTolerance", then, the procedure to 
		// compute t above is the following:
		//
		// On the efficient segment [w_min, w_max], it exist vectors alpha and beta such that
		// w_min = alpha + lambda_min * beta
		// w_max = alpha + lambda_max * beta
		// w     = alpha + constraintValue * beta
		// c.f. for instance formula 7.10a of the second reference.
		// From the two first equations, it is possible to deduce alpha = w_min - lambda_a * (w_max - w_min)/(lambda_max - lambda_min)
		// and beta = (w_max - w_min)/(lambda_max - lambda_min), so that w = ...
		// <=>
		// t = (lambda_max - constraintValue)/(lambda_max - lambda_min)
		//
		if (constraint === "return") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}		
		else if (constraint === "volatility") {		
			// Define the coefficients of the second order polynomial at^2 + bt + c
			var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
			var a = fct_min * fct_min + fct_max * fct_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
			var b = -2 * (fct_max * fct_max - variance_cross); // 
			var c = fct_max * fct_max - constraintValue*constraintValue; //always > 0
			
			// Extract the root t of the equation at^2 + bt + c = 0 belonging to ]0,1[, using a stable numerical formula
			var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
				throw new Error('internal error; the covariance matrix might not be semi-definite positive');
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
				throw new Error('internal error: the covariance matrix might not be semi-definite positive');
			}
		}
		else if (constraint === "riskTolerance") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}
		else {
			throw new Error('internal error: unknown constraint');
		}

		// Compute the final efficient portfolio weights
		var weights = Matrix_.fill(weights_min.nbRows, 1, 
							   	   function(i,j) { 
									   return t*weights_min.getValue(i, 1) + (1-t)*weights_max.getValue(i, 1); 
								   });
		
		// Return the computed portfolio weights
		return [weights, idx_min, idx_max];
	}
}


/**
* @function meanVarianceEfficientFrontierNearestWeights
*
* @summary Compute the weights of the nearest portfolio located on the mean variance
* efficient frontier.
*
* @description This function returns the weights w_1,...,w_n, associated to the portfolio located on the
* efficient frontier nearest to the input portfolio in a l^2 norm sense.
*
* The algorithm used internally to compute the efficient frontier is the Markowitz critical line algorithm, c.f. the reference.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {Array.<number>} inputWeights the weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, array of n real numbers.
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
* @return {Array<Array.<number>} the weights corresponding to the mean-variance efficient portfolio, array of n real numbers.
*
* @example
* meanVarianceEfficientFrontierNearestWeights()
* // TODO
*/
self.meanVarianceEfficientFrontierNearestWeights = function(inputWeights, mu, sigma, opt) {
	var eps = 1e-8; // the numerical zero
	
	// ------	

	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	
	// Convert the input weights to matrix format
	var inputWeights = new Matrix_(inputWeights);
	
	// ------

	
	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	if (cornerPortfolios.length == 0) {
		throw new Error('internal error: no corner portfolios');
	}
	
	
	// Compute the projection, in a l^2 norm sense, of the input portfolio on each of
	// the efficient segments making up the efficient frontier.
	//
	// The projection of the input portfolio on the efficient frontier, in a l^2 norm sense, 
	// is then by definition the projected portfolio having the minimum distance with 
	// the input portfolio.
	var weights;
	var minDist = Number.POSITIVE_INFINITY;
	for (var i = 0; i < cornerPortfolios.length - 1; ++i) {
		// Extract the end points of the current efficient segment
		var w_e = cornerPortfolios[i][0];
		var w_b = cornerPortfolios[i+1][0];
		
		// Project the input portfolio on the current efficient segment
		var proj = lineSegmentEuclidianProjection_(inputWeights, w_b, w_e);
		
		// Compute the l^2 distance between the input portfolio and the projected portfolio
		var dist = Matrix_.xmy(inputWeights, new Matrix_(proj)).vectorNorm('two');
		
		// Check if the projected portfolio is a candidate to be the projection of the
		// input portfolio on the whole efficient frontier.
		if (dist <= minDist) {
			weights = proj;
			minDist = dist;
		}
	}
	
	// Return the computed weights
	return weights;
}



/**
* @function meanVarianceEfficientFrontierPortfolios
*
* @summary Compute the weights, returns and volatilities of portfolios belonging 
* to the mean-variance efficient frontier.
*
* @description This function returns the weights w_i1,...,w_in, the returns r_i and the volatilities
* std_i, with i = 1..nbPortfolios, associated to nbPortfolios fully invested and long-only portfolios 
* of n assets belonging to the mean-variance efficient frontier, ordered from the lowest return/volatility
* portfolio to the highest return/volatility portfolio.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The main algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* The algorithm used internally generates the portfolios uniformly on the efficient frontier,
* with regard to the risk tolerance value.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.nbPortfolios the number of efficient portfolios to compute, a strictly positive natural integer; defaults to 100.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array.<Array.<Object>>} the weights, returns and volatilities of the computed efficient portfolios, an array of nbPortfolios arrays of three elements:
* - arr[0..nbPortfolios-1][0], the weights corresponding to an efficient portfolio, an array of n real numbers
* - arr[0..nbPortfolios-1][1], the return of the efficient portfolio, a real number
* - arr[0..nbPortfolios-1][2], the volatility of the efficient portfolio, a real number
*
* @example
* meanVarianceEfficientFrontierPortfolios([0.1, 0.2], [[1, 0.3], [0.3, 1]], {nbPortfolios: 5})
* // [[[0.5, 0.5], ~0.15, ~0.806], [[0.375, 0.625], 0.1625, ~0.820], [[0.25, 0.75], ~0.175, ~0.859], [[0.125, 0.875], ~0.1875, ~0.920], [[0, 1], 0.2, 1]]
*/
self.meanVarianceEfficientFrontierPortfolios = function(mu, sigma, opt) {
	var eps = 1e-8; // the numerical zero
	
	// Internal function to compute the return of a portfolio
	function portfolioReturn(x, mu, sigma) {
		return Matrix_.vectorDotProduct(mu, x);
	}
	
	// Internal function to compute the volatility of a portfolio
	function portfolioVolatility(x, mu, sigma) {	
		// Compute the variance x'*SIGMA*x
		var sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(sigma, x), x);		

		// In case the variance is numerically zero, which can occur with 
		// a semi-positive definite covariance matrix, it is replaced with zero.
		//
		// Otherwise, if the variance is negative, stops the algorithm,
		// since the covariance matrix is then not numerically semi-positive definite.
		if (Math.abs(sigma_x_x) <= eps) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -eps) {
			throw new Error('negative volatility, covariance matrix might not be semi-definite positive');
		}
		
		// Compute the volatility SQRT(x'*SIGMA*x)
		var s_x = Math.sqrt(sigma_x_x);
	   
		// Return the computed volatility
		return s_x;
	}

	
	// ------
	
	
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	var nbPortfolios = opt.nbPortfolios || 100;

	
	// ------
	
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);

	
	// Compute the corner portfolios defining the efficient frontier,
	// as well the minimum/maximum values of the risk aversion parameter.
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	
	// Initializations
	var nbAssets = sigma.nbColumns;
	var efficientFrontier = new Array(nbPortfolios);
	
	
	// Limit cases: 
	// - If there is only one corner portfolio on the efficient frontier,
	// return it directly.
	//
	// - Otherwise, the number of portfolios to compute must be greater than two
	// for the algorithm below to be valid
	if (cornerPortfolios.length == 0) {
		throw new Error('efficient frontier made of no corner portfolios: internal error');
	}
	else if (cornerPortfolios.length == 1) {
		if (nbPortfolios != 1) {
			throw new Error('efficient frontier made of only one corner portfolio: only one efficient portfolio can be computed');
		}
		else {
			var weights = cornerPortfolios[0][0];
			var ret = portfolioReturn(weights, mu, sigma);
			var vol = portfolioVolatility(weights, mu, sigma);
			
			return [[weights.toArray(), ret, vol]];
		}
	}
	else { // cornerPortfolios.length >= 2
		if (nbPortfolios <= 1) {
			throw new Error('efficient frontier made of several corner portfolios: at least two efficient portfolios must be computed');
		}
	}
	
	
	// Generate nbPortfolios regularly spaced distinct points lambda_i, i=0..nbPortfolios-1, 
	// belonging to the interval [lambda_min, lambda_max], using the formula
	// lambda_i = lambda_min + i * (lambda_max - 1)/(nbPortfolios - 1).
	//
	// Then, for each of these points, compute the efficient portfolio with a target 
	// risk tolerance parameter equal to lambda_i.

	// Initializations
	var lambda_min = cornerPortfolios[0][1];
	var lambda_max = cornerPortfolios[cornerPortfolios.length-1][1];
	var delta_lambda = (lambda_max - lambda_min)/(nbPortfolios - 1);
	
	for (var i = 0; i < nbPortfolios; ++i) {
		// Generate the current point t_i
		var lambda_i = lambda_min + i * delta_lambda;
		
		// Compute the efficient portfolio with a risk tolerance equal to lambda_i
		var weights = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, {constraint: "riskTolerance", 
																						   constraintValue: lambda_i});																						   
		if (weights.length == 0) {
			throw new Error('internal error: no matching efficient portfolio with a risk tolerance equal to ' + lambda_i);
		}
		else {
			weights = weights[0];
			
			var ret = portfolioReturn(weights, mu, sigma);
			var vol = portfolioVolatility(weights, mu, sigma);
			
			// Save the computed values, with the highest return/volatility portfolios 
			// provided last and the lowest return/volatility portfolios to be provided first.
			efficientFrontier[(nbPortfolios - 1) - i] = [weights.toArray(), ret, vol];
		}
	}
	
	
	// Return the computed list of portfolios weights, returns and volatilities
	return efficientFrontier;
}



/**
* @function meanVarianceEfficientFrontierCornerPortfolios
*
* @summary Compute the weights, returns and volatilities of the corner portfolios defining
* the mean-variance efficient frontier.
*
* @description This function returns the weights w_i1,...,w_in, the returns r_i and the volatilities
* std_i, with i = 1..m, associated to the m fully invested and long-only corner portfolios 
* of n assets defining the mean-variance efficient frontier, ordered from the lowest return/volatility
* portfolio to the highest return/volatility portfolio.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array.<Array.<Object>>} the weights, returns and volatilities of the m corner portfolios, an array of m arrays of three elements:
* - arr[0..m-1][0], the weights corresponding to a corner portfolio, an array of n real numbers
* - arr[0..m-1][1], the return of the corner portfolio, a real number
* - arr[0..m-1][2], the volatility of the corner portfolio, a real number
*
* @example
* meanVarianceEfficientFrontierCornerPortfolios([0.1, 0.2], [[1, 0.3], [0.3, 1]])
* // [[[0.5, 0.5], ~0.15, ~0.806], [[0, 1], 0.2, 1]]
*/
self.meanVarianceEfficientFrontierCornerPortfolios = function(mu, sigma, opt) {
	//
	var eps = 1e-8; // the numerical zero
	
	// Internal function to compute the return of a portfolio
	function portfolioReturn(x, mu, sigma) {
		return Matrix_.vectorDotProduct(mu, x);
	}
	
	// Internal function to compute the volatility of a portfolio
	function portfolioVolatility(x, mu, sigma) {	
		// Compute the variance x'*SIGMA*x
		var sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(sigma, x), x);		

		// In case the variance is numerically zero, which can occur with 
		// a semi-positive definite covariance matrix, it is replaced with zero.
		//
		// Otherwise, if the variance is negative, stops the algorithm,
		// since the covariance matrix is then not numerically semi-positive definite.
		if (Math.abs(sigma_x_x) <= eps) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -eps) {
			throw new Error('negative volatility, covariance matrix might not be semi-definite positive');
		}
		
		// Compute the volatility SQRT(x'*SIGMA*x)
		var s_x = Math.sqrt(sigma_x_x);
	   
		// Return the computed volatility
		return s_x;
	}

	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);

	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	// Convert the output of the internal function above to a list of 
	// portfolios weights, returns and volatilities.
	//
	// Filter the portfolio weights for numerically identical portfolios.

	// Initializations
	var efficientFrontier = [];
	
	// First, the E-maximizing portfolio is always included on the efficient frontier
	var weights = cornerPortfolios[0][0];
	var ret = portfolioReturn(weights, mu, sigma);
	var vol = portfolioVolatility(weights, mu, sigma);
	efficientFrontier.push([weights.toArray(), ret, vol]);
			
	// Then, for each computed corner portfolio:
	// - If it is numerically identical to the last corner portfolio included in the
	// output efficient frontier, do nothing
	//
	// - Otherwise, add it
	var currentWeights = weights;
	for (var i = 1; i < cornerPortfolios.length; ++i) {
		var weights = cornerPortfolios[i][0];
		
		//
		if (Matrix_.areEqual(weights, currentWeights, eps)) {
			continue;
		}
		
		//
		var ret = portfolioReturn(weights, mu, sigma);
		var vol = portfolioVolatility(weights, mu, sigma);
		efficientFrontier.push([weights.toArray(), ret, vol]);
		
		// Prepare for the next iteration
		currentWeights = weights;
	}

	// Return the computed list of portfolios weights, returns and volatilities,
	// sorted from lowest return/volatility to highest return/volatility.
	return efficientFrontier.reverse();
}


/**
* @function computeCornerPortfolios_
*
* @summary Compute all the corner portfolios belonging to the mean-variance efficient frontier.
*
* @description This function returns the weights w_i1,...,w_in as well as the risk aversion parameters lambda_i,
* i = 1..m, associated to the m fully invested and long-only corner portfolios defining the mean-variance
* efficient frontier.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The algorithm used internally is the Markowitz critical line algorithm, c.f. the first reference.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
* @see <a href="https://doi.org/10.1007/978-0-387-77439-8_12">Niedermayer A., Niedermayer D. (2010) Applying Markowitz’s Critical Line Algorithm. In: Guerard J.B. (eds) Handbook of Portfolio Construction. Springer, Boston, MA</a>
*
* @param {Matrix_} mu the returns of the n assets in the considered universe, a n by 1 matrix.
* @param {Matrix_} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square n by n matrix.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<Object>>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
*
* @example
* computeCornerPortfolios_(new Matrix_([0.1, 0.2]), new Matrix_([[1, 0.3], [0.3, 1]]))
* // [[new Matrix_([0, 1]), 7], [new Matrix_([0.5, 0.5]), 0]] 
*/
function computeCornerPortfolios_(mu, sigma, opt) {	
	var eps = 1e-8; // the numerical zero
	
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
		
		this.varIn = new BitSet_();
		this.varIn.resize(nbAssets + nbEqualityConstraints);
		this.varOut = new BitSet_();
		this.varOut.resize(nbAssets + nbEqualityConstraints);
		this.varLow = new BitSet_();
		this.varLow.resize(nbAssets + nbEqualityConstraints);
		this.varUp = new BitSet_();
		this.varUp.resize(nbAssets + nbEqualityConstraints);
			
		// Public functions to set the status of variables
		this.setIn = function(idx) {
			this.varIn.set(idx);
			this.varOut.unset(idx);
			this.varLow.unset(idx);
			this.varUp.unset(idx);
		}
		this.setOnLowerBound = function(idx) {
			this.varLow.set(idx);
			this.varOut.set(idx);
			this.varIn.unset(idx);
			this.varUp.unset(idx);
		}
		this.setOnUpperBound = function(idx) {
			this.varUp.set(idx);
			this.varOut.set(idx);
			this.varIn.unset(idx);
			this.varLow.unset(idx);
		}
		this.setLambdasIn = function() {
			for (var i = this.nbAssets + 1; i <= this.nbAssets + this.nbEqualityConstraints; ++i) {
				this.varIn.set(i);
				this.varOut.unset(i);
				this.varLow.unset(i);
				this.varUp.unset(i);
			}
		}
		this.setAssetsOnLowerBounds = function() {
			for (var i = 1; i <= this.nbAssets; ++i) {
				this.varLow.set(i);
				this.varOut.set(i);
				this.varIn.unset(i);
				this.varUp.unset(i);
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
			return this.varIn.get(idx);
		}
		this.isOnLowerBound = function(idx) {
			return this.varLow.get(idx);
		}
		this.isOnUpperBound = function(idx) {
			return this.varUp.get(idx);
		}
		this.isOut = function(idx) {
			return this.varOut.get(idx);
		}
		
		// Public functions to iterate over the different sets.
		this.getInIndexes = function() {
			return this.varIn.toArray();
		}
		this.getOutIndexes = function() {
			return this.varOut.toArray();
		}
	}

	// Internal function to compute the E-maximizing portfolio, 
	// c.f. the method "STARTING-SOLUTION" of the second reference, which
	// is the greedy algorithm used to solve the continuous knapsack problem,
	// c.f. https://en.wikipedia.org/wiki/Continuous_knapsack_problem.
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
		// Check that the problem is feasible (i.e., that the restricted unit simplex on which
		// the E-maximization is taking place is not empty, c.f. paragraph 12.3.1 of the second
		// reference).
		//
		// Note: if ever needed for performances reasons, the .toArray() calls below could be
		// replaced with another way.
		var sumBounds = simplexEmptinessCheck_(nbAssets, lowerBounds.toArray(), upperBounds.toArray());
		
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
		// are LOW, plus the linear program is not degenerate
		//
		//
		// - The loop above has stopped on an asset because the sum of the weights of 
		// all the assets is numerically equal to one
		//
		// In this case, all the assets with a return strictly higher than this asset are UP,
		// this asset is UP, and all the assets with a return strictly lower than this asset
		// are LOW, plus the linear program is degenerate
		//
		// To circumvent the degeneracy, this asset is forced to IN thanks to a numerical 
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
	

	// ------
	
	// Initializations	
	var nbAssets = sigma.nbColumns;
	
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
	// of the first reference (the E-maximizing portfolio with the lowest volatility).
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

	// In each iteration (except the first one), the asset that was determined 
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
		if (maxIterations !== -1 && iter >= maxIterations) {
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
				//
				// Numerical note: C.f. chapter 13 of the reference, paragraph "Deleting a variable",
				// the element Mi_out_idx_out_idx below is such that 1/Mi_out_idx_out_idx is strictly positive.
				//
				//                 So, in case it becomes numerically negative or zero (covariance matrix numerically
				// not semi-definite positive, ill-conditioned...), it is needed to truncate it.
				var Mi_out_idx_out_idx = Mi.getValue(idx_out, idx_out);				
				if (Mi_out_idx_out_idx <= eps) {
					Mi_out_idx_out_idx = eps;
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
				// Get the current IN variables indexes
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
				//
				// Numerical note: C.f. chapter 13 of the reference, paragraph "Adding a variable",
				// the element xi_j below is such that xi_j is strictly positive.
				//
				//                 So, in case it becomes numerically negative or zero (covariance matrix numerically
				// not semi-definite positive, ill-conditioned...), it is needed to truncate it.
				var xi_j = M.getValue(idx_in, idx_in);
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					xi_j -= M.getValue(idx_in, in_idx_i) * xi.getValue(in_idx_i, 1);
				}
				if (xi_j <= eps) {
					xi_j = eps;
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
				
				alpha_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * b_bar.getValue(in_idx_j, 1);
			
				if (in_idx_j <= nbAssets) {
					beta_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * mu.getValue(in_idx_j, 1);
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
	
	
	// Return the computed efficient frontier array
	return cornerPortfoliosWeights;	
}
