/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.MeanVarianceEfficientFrontierGsmo = MeanVarianceEfficientFrontierGsmo;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function MeanVarianceEfficientFrontierGsmo
*
* @description Object representing a mean-variance efficient frontier computed using the
* GSMO algorithm, c.f. the reference, and implementing the methods of the parent virtual
* object MeanVarianceEfficientFrontier.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
* @see <a href="https://link.springer.com/article/10.1023/A:1012431217818">Keerthi, S. & Gilbert, E. Convergence of a Generalized SMO Algorithm for SVM Classifier Design Machine Learning (2002) 46: 351.</a>
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the mean-variance optimization problem.
* @param {number} opt.optimizationMethodParams.maximumRiskToleranceValueOnlyGsmo set to true to only compute the maximum risk tolerance associated to the efficient frontier; defaults to false;
* @param {number} opt.optimizationMethodParams.minimumRiskToleranceValueOnlyGsmo set to true to only compute the minimum risk tolerance associated to the efficient frontier; defaults to false;
* @param {number} opt.optimizationMethodParams.epsGsmo the convergence tolerance of the algorithm used to solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo the maximum number of iterations of the algorithm used to solve the mean-variance optimization problem, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
*/
function MeanVarianceEfficientFrontierGsmo(mu, sigma, opt) {
	// Call the parent constructor
	MeanVarianceEfficientFrontier.call(this, mu, sigma, opt);
	
	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	this.epsGsmo = opt.optimizationMethodParams.epsGsmo;
	if (this.epsGsmo == undefined) {
		this.epsGsmo = 1e-6;
	}
	this.maxIterationsGsmo = opt.optimizationMethodParams.maxIterGsmo;
	if (this.maxIterationsGsmo == undefined) {
		this.maxIterationsGsmo = 10000;
	}
	this.antiCyclingGsmo = opt.optimizationMethodParams.antiCyclingGsmo;
	if (this.antiCyclingGsmo == undefined) {
		this.antiCyclingGsmo = false;
	}
	this.maximumRiskToleranceValueOnly = opt.optimizationMethodParams.maximumRiskToleranceValueOnlyGsmo;
	if (this.maximumRiskToleranceValueOnly == undefined) {
		this.maximumRiskToleranceValueOnly = false;
	}
	this.minimumRiskToleranceValueOnly = opt.optimizationMethodParams.minimumRiskToleranceValueOnlyGsmo;
	if (this.minimumRiskToleranceValueOnly == undefined) {
		this.minimumRiskToleranceValueOnly = false;
	}
	if (this.maximumRiskToleranceValueOnly == true && this.minimumRiskToleranceValueOnly == true) {
		throw new Error("internal error: inconsistent minimum/maximum risk tolerance values computation");
	}
	
	// Initialize the cache of risk tolerance parameters
	this.cachedEfficientPortfolios = new Map();
	
	// Compute the minimum and maximum risk tolerance values defining the efficient frontier
	if (this.minimumRiskToleranceValueOnly == false) {
		var h = computeMaximumRiskTolerancePortfolio.call(this);
		this.highestRiskTolerance = h[1];
		this.highestRiskTolerancePortfolio = h[0];
		
		// Cache the computed portfolio
		this.cachedEfficientPortfolios.set(this.highestRiskTolerance, {weights: this.highestRiskTolerancePortfolio,
		                                                               ret: this.getHighestReturn(), 
		                                                               volatility: this.getHighestVolatility()});
	}
	
	if (this.maximumRiskToleranceValueOnly == false) {
		var l = computeMinimumRiskTolerancePortfolio.call(this);
		this.lowestRiskTolerance = l[1];
		this.lowestRiskTolerancePortfolio = l[0];
		
		// Cache the computed portfolio
		this.cachedEfficientPortfolios.set(this.lowestRiskTolerance, {weights: this.lowestRiskTolerancePortfolio,
		                                                              ret: this.getLowestReturn(), 
		                                                              volatility: this.getLowestVolatility()});
	}

	
	// ------

	/**
	* @function computeMaximumRiskTolerancePortfolio
	*
	* @description This function computes a value for the risk tolerance parameter rt in the 
	* mean-variance optimization formulation, such that the associated portfolio is the 
	* maximum return/volatility portfolio.
	*
	* C.f. the reference, as there is an infinity of such values, the computed value can be any of them
	* and is unlikely to be the lowest one.
	*
	* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
	*
	* @memberof MeanVarianceEfficientFrontierGsmo
	*
	* @return {Array<Object>} an array arr containing two elements: 
	* - arr[0] an n by 1 matrix containing the weights of the computed portfolio
	* - arr[1] the associated risk tolerance
	*
	*/
	function computeMaximumRiskTolerancePortfolio() {
		// Initializations
		var nbAssets = this.nbAssets;
		var mu = this.mu;
		var sigma = this.sigma;
		var lowerBounds = this.lowerBounds;
		var upperBounds = this.upperBounds;
		
		// Compute the maximum return attainable on the efficient frontier
		var maxReturnSolution = simplexLpSolve_(mu.elemMap(function(i,j,val) { return -val;}), lowerBounds, upperBounds);
		var maxReturn = -maxReturnSolution[1];
		
		// Determine a risk tolerance parameter rt_E such that min_w <Vw/w> - 2*rt_E*<w/e>, s.t. <w/e>=1 and l <= x <= u
		// corresponds to the unique efficient portfolio with the maximum return attainable above.
		//
		// This is done by increasing more and more the value or rt.

		// Preliminary definitions
		var Q = sigma;
		var b = Matrix_.ones(nbAssets, 1);
		var r = 1;
		var l = lowerBounds;
		var u = upperBounds;

		// Compute an initial feasible point for the GSMO algorithm below
		var centroid = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
		var p_centroid = qksolveBS_(Matrix_.ones(nbAssets, 1), centroid, b, r, l, u);
		var x0 = p_centroid[0];
		
		// Core loop, which must converges as rt is brought to infinity
		var epsSearch =  1e-12;
		var nbIterSearch = 0;
		var maxNbIterSearch = 54; // The default is taken to be 54, because 2^54 * 32 is quite close to an already unreasonable high value
		var rt = 32;
		var efficientPortfolio;
		do {
			// Increment the number of iterations and check that the number of iterations stays reasonable
			++nbIterSearch;
			if (nbIterSearch > maxNbIterSearch) {
				throw new Error('internal error: maximum number of iterations reached when searching for the efficient portfolio with maximum return');
			}
			
			// Double the risk tolerance value and update the associated vector value
			rt = 2 * rt;
			var p = Matrix_.ax(-rt, mu);
			
			// Compute the efficient portfolio associated to the risk tolerance value above,
			// as well as its return.
			var efficientSol = qpsolveGSMO_(Q, p, b, r, l, u, {x0: x0, antiCycling: this.antiCyclingGsmo, eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
			efficientPortfolio = efficientSol[0];
			
			var ret = this.computePortfolioReturn(efficientPortfolio);
			
			// Cache the computed portfolio
			this.cachedEfficientPortfolios.set(rt, {weights: efficientPortfolio,
			                                        ret: ret, 
			                                        volatility: this.computePortfolioVolatility(efficientPortfolio)});
		}
		while ( Math.abs(ret - maxReturn) > epsSearch );

		// Return the computed portfolio, as well as the associated risk tolerance.
		return [efficientPortfolio, rt];
	}
	
	
	/**
	* @function computeMinimumRiskTolerancePortfolio
	*
	* @description This function computes a value for the risk tolerance parameter rt in the 
	* mean-variance optimization formulation, such that the associated portfolio is the 
	* minimum return/volatility portfolio.
	*
	* C.f. the reference, this value is not necessarily 0, because in case the covariance matrix is 
	* semi-positive definite, there is no reason why a portfolio minimizing the volatility would 
	* also be maximizing the return.
	*
	* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
	*
	* @memberof MeanVarianceEfficientFrontierGsmo
	*
	* @return {Array<Object>} an array arr containing two elements: 
	* - arr[0] an n by 1 matrix containing the weights of the computed portfolio
	* - arr[1] the associated risk tolerance
	*
	*/
	function computeMinimumRiskTolerancePortfolio() {
		// Initializations
		var nbAssets = this.nbAssets;
		var mu = this.mu;
		var sigma = this.sigma;
		var lowerBounds = this.lowerBounds;
		var upperBounds = this.upperBounds;
		
		// Compute the minimum volatility attainable on the efficient frontier
		// Preliminary definitions
		var Q = sigma;
		var p = Matrix_.zeros(nbAssets, 1);
		var b = Matrix_.ones(nbAssets, 1);
		var r = 1;
		var l = lowerBounds;
		var u = upperBounds;

		// Compute an initial feasible point for the GSMO algorithms below
		var centroid = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
		var p_centroid = qksolveBS_(Matrix_.ones(nbAssets, 1), centroid, b, r, l, u);
		var x0 = p_centroid[0];
		
		//
		var minVolatilitySolution = qpsolveGSMO_(Q, p, b, r, l, u, {x0: x0, antiCycling: this.antiCyclingGsmo, eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
		var minVolatility = this.computePortfolioVolatility(minVolatilitySolution[0]);
		
		// Determine a risk tolerance parameter rt_E such that min_w <Vw/w> - 2*rt_E*<w/e>, s.t. <w/e>=1 and l <= x <= u
		// corresponds to the unique efficient portfolio with the minimum volatility attainable above.
		//
		// This is done by decreasing more and more the value or rt.

		// Core loop, which must converges as rt is brought to zero
		//
		// There is no need to set a maximum number of iterations here, because
		// rt will quickly become numerically null, so that the computed portfolio will have the exact
		// same volatility as the minimum volatility portfolio above.
		var epsSearch =  1e-12;
		var rt = 2;
		var efficientPortfolio;
		do {
			// Halve the risk tolerance value and update the associated vector value
			//
			// Note that being more aggressive (ex: /10 instead of /2) will lead to
			// the computed portfolio NOT being efficient in case the covariance matrix
			// is semi-definite positive.
			rt = rt / 2;
			var p = Matrix_.ax(-rt, mu);
			
			// Compute the efficient portfolio associated to the risk tolerance value above,
			// as well as its return.
			var efficientSol = qpsolveGSMO_(Q, p, b, r, l, u, {x0: x0, antiCycling: this.antiCyclingGsmo, eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
			efficientPortfolio = efficientSol[0];
			
			var vol = this.computePortfolioVolatility(efficientPortfolio);
			
			// Cache the computed portfolio
			this.cachedEfficientPortfolios.set(rt, {weights: efficientPortfolio,
			                                        ret: this.computePortfolioReturn(efficientPortfolio), 
												    volatility: vol});
		}
		while ( Math.abs(vol - minVolatility) > epsSearch );

		// Return the computed portfolio, as well as the associated risk tolerance.
		return [efficientPortfolio, rt];
	}

}
MeanVarianceEfficientFrontierGsmo.prototype = Object.create(MeanVarianceEfficientFrontier.prototype);
MeanVarianceEfficientFrontierGsmo.prototype.constructor = MeanVarianceEfficientFrontierGsmo;

MeanVarianceEfficientFrontierGsmo.prototype.getHighestRiskTolerancePortfolio = function(x) {
	//
	return this.highestRiskTolerancePortfolio;
};
MeanVarianceEfficientFrontierGsmo.prototype.getHighestRiskTolerance = function(x) {
	//
	return this.highestRiskTolerance;
};
MeanVarianceEfficientFrontierGsmo.prototype.getLowestRiskTolerancePortfolio = function(x) {
	//
	return this.lowestRiskTolerancePortfolio;
};
MeanVarianceEfficientFrontierGsmo.prototype.getLowestRiskTolerance = function(x) {
	//
	return this.lowestRiskTolerance;
};

/**
* @function computeEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and
* long-only efficient portfolio belonging to the mean-variance efficient frontier, subject to 
* a return constraint, a volatility constraints or to a risk tolerance constraint.
*
* @memberof MeanVarianceEfficientFrontierGsmo
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "return", to specify a return constraint
* - "volatility",  to specify a volatility constraint
* - "riskTolerance", to specify a risk tolerance constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
* @return {Array.<Object>} in case no efficient portfolio can be computed, return an empty array;
* in case an efficient portfolio satisfying the input constraint has been computed, an array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontierGsmo.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	// Decode inputs
	if (constraintType === undefined || constraintType === null) {
		throw new Error('internal error: missing constraint type');
	}

	var constraintFct;
	var constraintFctFromCache;
	var that = this;
	if (constraintType == "return") {
		constraintFct = function (portfolio) { return that.computePortfolioReturn(portfolio); };
		constraintFctFromCache = function (cachedPortfolio) { return cachedPortfolio.ret; };
	}
	else if (constraintType == "volatility") {
		constraintFct = function (portfolio) { return that.computePortfolioVolatility(portfolio); };
		constraintFctFromCache = function (cachedPortfolio) { return cachedPortfolio.volatility; };
	}
	else if (constraintType == "riskTolerance") {
		// Nothing to do, as the risk tolerance is directly accessible.
	}
	else {
		throw new Error('internal error: unknown constraint type');
	}
	
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('internal error: missing constraint value');
	}
	
	
	// Initializations
	var nbAssets = this.nbAssets;
	var mu = this.mu;
	var sigma = this.sigma;
	var lowerBounds = this.lowerBounds;
	var upperBounds = this.upperBounds;
	var eps = this.epsEfficientPortfolioComputation;
	
	// Preliminary definitions
	var Q = sigma;
	var b = Matrix_.ones(nbAssets, 1);
	var r = 1;
	var l = lowerBounds;
	var u = upperBounds;
	
	
	// Core algorithm, using the GSMO optimization algorithm to compute the desired 
	// efficient portfolio.
	//
	// In case a return or a volatility constraint is provided in input, a bisection
	// search is done on the efficient frontier (parametrized by the risk tolerance parameter)
	// in order to find the associated efficient portfolio, if it exists.
	//
	// In case a risk tolerance constraint is provided in input, the associated efficient
	// portfolio is obtained directly.
	var portfolioWeights;
	var riskTolerance;
	if (constraintType === "return" || constraintType === "volatility") {
		// Compute the minimum and maximum constraint function value on the efficient frontier
		var weights_min = this.getLowestRiskTolerancePortfolio();
		var constraintFct_min = constraintFct(weights_min);
		
		var weights_max = this.getHighestRiskTolerancePortfolio()
		var constraintFct_max = constraintFct(weights_max);
		
		// If the desired target constraint function value is numerically strictly greater than the highest attainable
		// constraint function value on the efficient frontier, or if the desired target constraint function 
		// value is numerically strictly lower than the lowest attainable constraint function value on the efficient frontier, 
		// there is no matching portfolio on the efficient frontier.
		if (constraintValue > constraintFct_max + eps || constraintValue < constraintFct_min - eps) {
			return [];
		}

		// Otherwise, if the target constraint function value is numerically reached on one of the
		// two extremal efficient portfolios, return immediately.
		if (Math.abs(constraintValue - constraintFct_min) <= eps) {
			portfolioWeights = weights_min;
			riskTolerance = this.getLowestRiskTolerance();
		}
		else if (Math.abs(constraintValue - constraintFct_max) <= eps) {
			portfolioWeights = weights_max;
			riskTolerance = this.getHighestRiskTolerance();
		}
		else {
			// Otherwise, compute the efficient portfolio with a constraint function value equal
			// to the target constraint function value through a bisection search on the risk tolerance parameter.
			
			// Compute an initial feasible point for the GSMO algorithm below
			var centroid = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
			var p_centroid = qksolveBS_(Matrix_.ones(nbAssets, 1), centroid, b, r, l, u);
			portfolioWeights = p_centroid[0];
			
			//
			var that = this;
			riskTolerance = bisection_(function (rt) { 
											// Check if the portfolio associated to the risk tolerance value rt has already been computed
											var cachedPortfolio = that.cachedEfficientPortfolios.get(rt);
											if (cachedPortfolio) {
												portfolioWeights = cachedPortfolio.weights;
												
												return constraintFctFromCache(cachedPortfolio) - constraintValue;
											}
											
											// Otherwise, compute the efficient portfolio solution to the optimization problem with a given risk tolerance parameter
											//
											// Note: The initial feasible point is taken to be the computed portfolio weights from the previous bisection 
											// iteration to warm start the GSMO algorithm.
											var p = Matrix_.ax(-rt, mu);
											var sol = qpsolveGSMO_(Q, p, b, r, l, u, {x0: portfolioWeights, antiCycling: that.antiCyclingGsmo, eps: that.epsGsmo, maxIter: that.maxIterationsGsmo});
											portfolioWeights = sol[0];
											
											// Cache the portfolio associated to the risk tolerance value rt
											that.cachedEfficientPortfolios.set(rt, {weights: portfolioWeights,
											                                        ret: that.computePortfolioReturn(portfolioWeights), 
											                                        volatility: that.computePortfolioVolatility(portfolioWeights)});

											// Return the value of the function
											return constraintFct(portfolioWeights) - constraintValue; 
										}, 
										this.getLowestRiskTolerance(), this.getHighestRiskTolerance());
		}
	}		
	else if (constraintType === "riskTolerance") {
		// Directly compute the efficient portfolio solution to the optimization problem
		var rt = constraintValue;
		var p = Matrix_.ax(-rt, mu);
		var sol = qpsolveGSMO_(Q, p, b, r, l, u, {antiCycling: this.antiCyclingGsmo, eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
		
		portfolioWeights = sol[0];
		riskTolerance = rt;
	}
	else {
		throw new Error('internal error: unknown constraint type');
	}

		
	// Return the computed portfolio, as well as the associated risk tolerance.
	return [portfolioWeights, riskTolerance];
};

/**
* @function restrict
*
* @description This function restricts the efficient frontier to portfolios satisfying 
* a minimal return constraint or a minimal volatility constraint.
*
* @memberof MeanVarianceEfficientFrontierGsmo
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "minReturn", to specify a return constraint
* - "minVolatility",  to specify a volatility constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
*/
MeanVarianceEfficientFrontierGsmo.prototype.restrict = function(constraintType, constraintValue) {
	// Decode the input parameters
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}

	// Restrict the efficient frontier, if possible
	if (constraintType == "minReturn") {	
		var lowestReturn = this.getLowestReturn();
		if (lowestReturn < constraintValue) {
			// The efficient frontier needs to be restricted, and this might not be possible
			var highestReturn = this.getHighestReturn();
			if (highestReturn < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum return constraint is not feasible');
			}
			
			// Compute the first efficient portfolio with a strictly positive enough return
			var efficientPortfolio = this.computeEfficientPortfolio("return", constraintValue + 0.5 * this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a return greater than ' + constraintValue);
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			
			// Rebuild the minimum and maximum risk tolerance values defining the efficient frontier
			this.lowestRiskTolerancePortfolio = efficientPortfolioWeights;
			this.lowestRiskTolerance = efficientPortfolioRiskTolerance;
			if (this.lowestRiskTolerance >= this.highestRiskTolerance) {
				this.highestRiskTolerancePortfolio = efficientPortfolioWeights;
				this.highestRiskTolerance = efficientPortfolioRiskTolerance;
			}			
		}
	}
	else if (constraintType == "minVolatility") {	
		var lowestVolatility = this.getLowestVolatility();
		if (lowestVolatility < constraintValue) {
			// The efficient frontier needs to be restricted, and this might not be possible
			var highestVolatility = this.getHighestVolatility();
			if (highestVolatility < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum volatility constraint is not feasible');
			}
			
			// Compute the first efficient portfolio with a strictly positive enough volatility
			var efficientPortfolio = this.computeEfficientPortfolio("volatility", constraintValue + 0.5 * this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a volatility greater than ' + constraintValue + ': the covariance matrix might not be semi-definite positive');
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];

			// Rebuild the minimum and maximum risk tolerance values defining the efficient frontier
			this.lowestRiskTolerancePortfolio = efficientPortfolioWeights;
			this.lowestRiskTolerance = efficientPortfolioRiskTolerance;
			if (this.lowestRiskTolerance >= this.highestRiskTolerance) {
				this.highestRiskTolerancePortfolio = efficientPortfolioWeights;
				this.highestRiskTolerance = efficientPortfolioRiskTolerance;
			}			
		}
	}
	else {
		throw new Error('unknown constraint type');
	}
	
	// Nothing to do here
};


/**
* @function computeMaximumSharpeRatioEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only portfolio of n assets maximizing the Sharpe ratio, which is defined as the ratio of the 
* portfolio excess return over a constant risk-free rate to the portfolio volatility.
*
* @see <a href="https://doi.org/10.1111/j.1540-6261.1976.tb03217.x">Elton, E. J., Gruber, M. J. and Padberg, M. W. (1976), SIMPLE CRITERIA FOR OPTIMAL PORTFOLIO SELECTION. The Journal of Finance, 31: 1341-1357</a>
* @see <a href="http://dx.doi.org/10.1080/13504860701255292">S. V. Stoyanov , S. T. Rachev & F. J. Fabozzi (2007) Optimal Financial Portfolios, Applied Mathematical Finance, 14:5, 401-436</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1437644">Leonid Kopman, Scott Liu, Maximizing the Sharpe Ratio and Information Ratio in the Barra Optimizer, MSCI Barra Research</a>
*
* @memberof MeanVarianceEfficientFrontierGsmo
*
* @param {number} rf the risk-free rate, a real number.
*
* @return {Array.<Object>} An array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontierGsmo.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	// The efficient frontier must be restricted to both:
	// - The domain of definition of the Sharpe ratio (the portfolios with a strictly positive volatility)
	// - The domain of strict positivity of the Sharpe ratio

	// On the restricted efficient frontier, the Sharpe ratio is a pseudo-concave 
	// function, c.f. especially the third reference, so that it is unimodal.
	//
	// This property allows to search for the portfolio with the maximum Sharpe ratio using
	// the golden section algorithm.
	
	// Initializations
	var nbAssets = this.nbAssets;
	var mu = this.mu;
	var sigma = this.sigma;
	var lowerBounds = this.lowerBounds;
	var upperBounds = this.upperBounds;
	var eps = this.epsEfficientPortfolioComputation;

	// Preliminary definitions
	var Q = sigma;
	var b = Matrix_.ones(nbAssets, 1);
	var r = 1;
	var l = lowerBounds;
	var u = upperBounds;
	
	// Core algorithm
	var portfolioWeights;
	var that = this;
	
	var riskTolerance = goldenSectionSearch_(function(rt) { 
												  // Compute the efficient portfolio solution to the optimization problem with a given risk tolerance parameter
												  var p = Matrix_.ax(-rt, mu);
												  var sol = qpsolveGSMO_(Q, p, b, r, l, u, {antiCycling: that.antiCyclingGsmo, eps: that.epsGsmo, maxIter: that.maxIterationsGsmo});
														
												  portfolioWeights = sol[0];
													
												  // Return the (opposite) value of the portfolio Sharpe ratio
												  return -that.computePortfolioSharpeRatio(portfolioWeights, rf); 
											 }, 
											 this.getLowestRiskTolerance(), this.getHighestRiskTolerance());
	
	// Return the computed portfolio, as well as the associated risk tolerance.
	return [portfolioWeights, riskTolerance[0]];
};