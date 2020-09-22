/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
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
* To be noted that when opt.optimizationMethod is set to "critical-line", an error might be raised in certain cases
* when the covariance matrix is not definite positive.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt parameters for the mean-variance optimization algorithm.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm from Markowitz, c.f. the reference.
* - 'gsmo': usage of generalized sequential minimization optimization algorithm
* - 'automatic': automatic selection of the optimization method to use
; defaults to 'automatic'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm in case opt.optimizationMethod is set to "critical-line",
* a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo in case opt.optimizationMethod is set to "gsmo", the convergence tolerance of the GSMO algorithm used to 
* solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo in case opt.optimizationMethod is set to "gsmo", the maximum number of iterations of the GSMO algorithm 
* used to solve the mean-variance optimization problem, a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
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
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}

	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'automatic';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo' && optimizationMethod != 'automatic') {
		throw new Error('unsupported optimisation method');
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
	
	
	// ------
		
	
	// Compute the efficient frontier
	var efficientFrontier;
	if (optimizationMethod == "automatic") {
		efficientFrontier = new MeanVarianceEfficientFrontierWrapper(mu, sigma, opt);
	}
	else if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}

	
	// Depending on the constraint, proceed with a different algorithm
	// to compute the requested efficient portfolio.
	var efficientPortfolioWeights;
	if (returnConstraint !== undefined) {
		var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("return", returnConstraint);
		if (efficientPortfolio.length == 0) {
			throw new Error('no matching efficient portfolio with a return equal to ' + returnConstraint);
		}
		else {
			efficientPortfolioWeights = efficientPortfolio[0];
		}
	}
	else if (volatilityConstraint !== undefined) {
		var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("volatility", volatilityConstraint);
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
		var highestRiskTolerancePortfolioWeights = efficientFrontier.getHighestRiskTolerancePortfolio();
		var highestRiskTolerance = efficientFrontier.getHighestRiskTolerance();
		var lowestRiskTolerancePortfolioWeights = efficientFrontier.getLowestRiskTolerancePortfolio();
		var lowestRiskTolerance = efficientFrontier.getLowestRiskTolerance();
		
		// If the desired risk tolerance is greater than the highest attainable
		// risk tolerance on the efficient frontier, the efficient portfolio is computed
		// as the rightmost portfolio on the efficient frontier, because this portfolio 
		// also corresponds to all risk tolerances greater than or equal to the
		// highest attainable risk tolerance on the efficient frontier.
		//
		// If the desired risk tolerance is lower than the lowest attainable
		// risk tolerance on the efficient frontier, the efficient portfolio is computed
		// as the leftmost portfolio on the efficient frontier.
		//
		// Otherwise, the efficient portfolio is computed as the efficient portfolio 
		// with a risk tolerance strictly equal to the desired risk tolerance.
		var eps = efficientFrontier.epsEfficientPortfolioComputation;
		if (riskToleranceConstraint > highestRiskTolerance - eps) {
			efficientPortfolioWeights = highestRiskTolerancePortfolioWeights;
		}
		else if (riskToleranceConstraint < lowestRiskTolerance + eps) {
			efficientPortfolioWeights = lowestRiskTolerancePortfolioWeights;
		}
		else {
			var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("riskTolerance", riskToleranceConstraint);
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
		var highestVolatilityPortfolioWeights = efficientFrontier.getHighestVolatilityPortfolio();
		var highestVolatility = efficientFrontier.getHighestVolatility();
		
		// If the desired maximum volatility is greater than the highest attainable
		// volatility on the efficient frontier, the efficient portfolio is computed
		// as the rightmost portfolio on the efficient frontier.
		//
		// Otherwise, the efficient portfolio is computed as the efficient portfolio 
		// with a volatility strictly equal to the desired maximum volatility.		
		var eps = efficientFrontier.epsEfficientPortfolioComputation;
		if (maxVolatilityConstraint > highestVolatility - eps) {
			efficientPortfolioWeights = highestVolatilityPortfolioWeights;
		}
		else {
			var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("volatility", maxVolatilityConstraint);
			if (efficientPortfolio.length == 0) {
				throw new Error('no matching efficient portfolio with a volatility lower than or equal to ' + maxVolatilityConstraint);
			}
			else {
				efficientPortfolioWeights = efficientPortfolio[0];
			}
		}
	}	

	
	// Return the computed portfolio weights
	return efficientPortfolioWeights.toArray();
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
* Optionally, the following constraints can be added:
* - Partial investment constraint, replacing the full investment constraint
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* To be noted that when opt.optimizationMethod is set to "critical-line", it is required that the provided assets
* returns are all different; otherwise, an error will be raised.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {Array.<number>} inputWeights the weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, array of n real numbers.
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm from Markowitz, c.f. the reference.
* - 'gsmo': usage of generalized sequential minimization optimization algorithm
* - 'automatic': automatic selection of the optimization method to use
; defaults to 'automatic'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm in case opt.optimizationMethod is set to "critical-line",
* a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo in case opt.optimizationMethod is set to "gsmo", the convergence tolerance of the GSMO algorithm used to 
* solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo in case opt.optimizationMethod is set to "gsmo", the maximum number of iterations of the GSMO algorithm 
* used to solve the mean-variance optimization problem, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.optimizationMethodParams.nbPortfoliosGsmo in case opt.optimizationMethod is set to "gsmo", the number of efficient portfolios to generate in order approximate the 
* efficient frontier corner portfolios, a positive integer greater than or equal to 2; defaults to 100.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
* @return {Array<Array.<number>} the weights corresponding to the mean-variance efficient portfolio, array of n real numbers.
*
*/
self.meanVarianceEfficientFrontierNearestWeights = function(inputWeights, mu, sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}

	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'automatic';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo' && optimizationMethod != 'automatic') {
		throw new Error('unsupported optimisation method');
	}

	
	// ------

	// Convert the input weights to matrix format
	var inputWeights = new Matrix_(inputWeights);

	
	// Compute the efficient frontier
	var efficientFrontier;
	
	if (optimizationMethod == "automatic") {
		efficientFrontier = new MeanVarianceEfficientFrontierWrapper(mu, sigma, opt);
	}
	else if (optimizationMethod == "critical-line") {
		efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	
	
	// Compute the portfolios representing the efficient frontier.
	//
	// In case the optimization method is the critical line, these representing
	// portfolios are the corner portfolios, and the computations below will be exact.
	//
	// Otherwise, these representing portfolios are equidistant points on the
	// efficient frontier, aim being to approximate the corner portfolios.
	var efficientFrontierPortfolios;
	if (optimizationMethod == "critical-line" || 
	    efficientFrontier.efficientFrontierOptimizationMethod == "critical-line") {
		efficientFrontierPortfolios = efficientFrontier.getCornerPortfolios();
	}
	else if (optimizationMethod == "gsmo" || 
	         efficientFrontier.efficientFrontierOptimizationMethod == "gsmo") {
		// Decode the number of efficient portfolios to generate
		var nbEfficientPortfolios = opt.optimizationMethodParams.nbPortfoliosGsmo;
		if (nbEfficientPortfolios == undefined) {
			nbEfficientPortfolios = 100;
		}

		// Generate the efficient portfolios as regularly spaced portfolios
		// on the efficient frontier with a risk tolerance value of rt_i, i=0..nbEfficientPortfolios-1
		// with rt_i belonging to the interval [rt_min, rt_max], 
		// using the formula rt_i = rt_min + i * (rt_max - 1)/(nbEfficientPortfolios - 1).
		efficientFrontierPortfolios = new Array(nbEfficientPortfolios);

		var rt_min = efficientFrontier.getLowestRiskTolerance();
		var rt_max = efficientFrontier.getHighestRiskTolerance();	
		var delta_rt = (rt_max - rt_min)/(nbEfficientPortfolios - 1);	
		for (var i = 0; i < nbEfficientPortfolios; ++i) {
			// Generate the current point t_i
			var rt_i = rt_min + i * delta_rt;
			
			// Compute the efficient portfolio with a risk tolerance equal to rt_i
			var weights = efficientFrontier.computeEfficientPortfolio("riskTolerance", rt_i);
			if (weights.length == 0) {
				throw new Error('internal error: no matching efficient portfolio with a risk tolerance of ' + rt_i);
			}
			else {
				efficientFrontierPortfolios[i] = weights[0];
			}
		}
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	
	// Compute the projection, in a l^2 norm sense, of the input portfolio on each of
	// the efficient segments making up the efficient frontier.
	//
	// The projection of the input portfolio on the efficient frontier, in a l^2 norm sense, 
	// is then by definition the projected portfolio having the minimum distance with 
	// the input portfolio.
	var weights;
	var minDist = Number.POSITIVE_INFINITY;
	for (var i = 0; i < efficientFrontierPortfolios.length - 1; ++i) {
		// Extract the end points of the current efficient segment
		var w_e = efficientFrontierPortfolios[i];
		var w_b = efficientFrontierPortfolios[i+1];
		
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
* - Partial investment constraint, replacing the full investment constraint
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The algorithm used internally generates the portfolios uniformly on the efficient frontier,
* with regard to either the return or the volatility or the risk tolerance value, depending on the
* opt.discretization parameter value.
*
* To be noted that when opt.optimizationMethod is set to "critical-line", it is required that the provided assets
* returns are all different; otherwise, an error will be raised.
*
* To be noted that when opt.optimizationMethod is set to "gsmo" and opt.discretizationType is set to "riskTolerance",
* identical portfolios can appear in output.
*
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the algorithm.
* @param {string} opt.discretizationType, a string either equal to:
* - "return" to compute portfolios uniformly spaced w.r.t. their returns
* - "volatility" to compute portfolios uniformly spaced w.r.t. their volatility
* - "riskTolerance" to compute portfolios uniformly spaced w.r.t. their risk tolerance parameter
*; defaults to "riskTolerance".
* @param {number} opt.nbPortfolios the number of efficient portfolios to compute, a strictly positive natural integer; defaults to 100.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio, a string either equals to:
* - 'critical-line': usage of the critical line algorithm from Markowitz, c.f. the reference.
* - 'gsmo': usage of  generalized sequential minimization optimization algorithm
* - 'automatic': automatic selection of the optimization method to use
; defaults to 'automatic'.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm in case opt.optimizationMethod is set to "critical-line",
* a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo in case opt.optimizationMethod is set to "gsmo", the convergence tolerance of the GSMO algorithm used to 
* solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo in case opt.optimizationMethod is set to "gsmo", the maximum number of iterations of the GSMO algorithm 
* used to solve the mean-variance optimization problem, a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
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
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}

	// Decode the optimization algorithm to use
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'automatic';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo' && optimizationMethod != 'automatic') {
		throw new Error('unsupported optimisation method');
	}

	var nbPortfolios = opt.nbPortfolios || 100;

	var discretizationType = opt.discretizationType;
	if (discretizationType == undefined) {
		discretizationType = "return";
	}
	if (discretizationType != "return" && discretizationType != "volatility" && discretizationType != "riskTolerance") {
		throw new Error('unsupported discretization type');
	}
	
	// ------
		
	// Compute the efficient frontier
	var efficientFrontier;
	if (optimizationMethod == "automatic") {
		efficientFrontier = new MeanVarianceEfficientFrontierWrapper(mu, sigma, opt);
	}
	else if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}

	
	// Initializations
	var efficientFrontierPortfolios = new Array(nbPortfolios);
	
	
	// Generate nbPortfolios regularly spaced distinct points t_i, i=0..nbPortfolios-1, 
	// belonging to the interval [t_min, t_max], using the formula
	// t_i = t_min + i * (t_max - 1)/(nbPortfolios - 1).
	//
	// Then, for each of these points, compute the efficient portfolio with a target 
	// constraint value equal to t_i.
	//
	// Limit case: if nbPortfolios == 1, (t_max - 1)/(nbPortfolios - 1) is not defined,
	// and the only computed portfolio corresponds to  t_i = (t_min + t_max)/2.

	// Initializations
	var t_min;
	var t_max;
	if (discretizationType == "return") {
		t_min = efficientFrontier.getLowestReturn();
		t_max = efficientFrontier.getHighestReturn();
	}
	else if (discretizationType == "volatility") {
		t_min = efficientFrontier.getLowestVolatility();
		t_max = efficientFrontier.getHighestVolatility();
	}
	else if (discretizationType == "riskTolerance") {
		t_min = efficientFrontier.getLowestRiskTolerance();
		t_max = efficientFrontier.getHighestRiskTolerance();
	}
	else {
		throw new Error('internal error: unsupported discretization type');
	}
	var delta_t = nbPortfolios == 1 ? -1 : (t_max - t_min)/(nbPortfolios - 1);
	var t_start = nbPortfolios == 1 ? (t_min + t_max)/2 : t_min;
	
	// Core loop
	for (var i = 0; i < nbPortfolios; ++i) {
		// Generate the current point t_i
		var t_i = t_start + i * delta_t;
		
		// Compute the efficient portfolio with a constraint value equal to t_i
		var weights = efficientFrontier.computeEfficientPortfolio(discretizationType, t_i);
		if (weights.length == 0) {
			throw new Error('internal error: no matching efficient portfolio with a constraint value ' + discretizationType + ' equal to ' + t_i);
		}
		else {
			//
			weights = weights[0];
			
			// Compute the portfolio return and volatility
			var ret = efficientFrontier.computePortfolioReturn(weights);
			var vol = efficientFrontier.computePortfolioVolatility(weights);
			
			// Save the computed values
			efficientFrontierPortfolios[i] = [weights.toArray(), ret, vol];
		}
	}
	
	
	// Return the computed list of portfolios weights, returns and volatilities
	return efficientFrontierPortfolios;
}

