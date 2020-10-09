/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.MeanVarianceEfficientFrontierWrapper = MeanVarianceEfficientFrontierWrapper;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function MeanVarianceEfficientFrontierWrapper
*
* @description Object representing a mean-variance efficient frontier computed using the
* most efficient algorithm among the available ones, and implementing the methods of the parent virtual
* object MeanVarianceEfficientFrontier.
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the mean-variance optimization problem.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm, a strictly positive natural integer,
* or -1 to force an infinite number of iterations; defaults to 1000.
* @param {number} opt.optimizationMethodParams.epsGsmo the convergence tolerance of the algorithm used to solve the mean-variance optimization problem, a strictly positive number; defaults to 1e-6.
* @param {number} opt.optimizationMethodParams.maxIterGsmo the maximum number of iterations of the algorithm used to solve the mean-variance optimization problem, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @param {boolean} opt.optimizationMethodParams.antiCyclingGsmo activate an anti cycling rule in the algorithm used to solve the mean-variance optimization problem, at the expense of execution time and stochasticity of the result; defaults to false.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
*/
function MeanVarianceEfficientFrontierWrapper(mu, sigma, opt) {
	// First try to construct an efficient frontier using the critical line algorithm,
	// as this is the most efficient algorithm.
	//
	// If the construction fails due to the singularity of the first KKT system,
	// switches to using the GSMO algorithm.
	this.efficientFrontier = null;
	try {
		this.efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
		this.efficientFrontierOptimizationMethod = "critical-line";
	}
	catch(e) {
		if (e.message == "internal error: impossible to solve the KKT system") {
			this.efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
			this.efficientFrontierOptimizationMethod = "gsmo";
		}
		else {
			throw (e);
		}
	}
	
	// Initializations
	this.epsPortfolioVolatility = this.efficientFrontier.epsPortfolioVolatility;
	this.epsEfficientPortfolioComputation = this.efficientFrontier.epsEfficientPortfolioComputation;
	
	this.mu = this.efficientFrontier.mu;
	this.sigma = this.efficientFrontier.sigma;
	this.nbAssets = this.efficientFrontier.mu.nbRows;
	
	this.lowerBounds = this.efficientFrontier.lowerBounds;
	this.upperBounds = this.efficientFrontier.upperBounds;
	
	// ------
}
MeanVarianceEfficientFrontierWrapper.prototype = Object.create(MeanVarianceEfficientFrontier.prototype);
MeanVarianceEfficientFrontierWrapper.prototype.constructor = MeanVarianceEfficientFrontierWrapper;

MeanVarianceEfficientFrontierWrapper.prototype.getHighestRiskTolerancePortfolio = function(x) {
	return this.efficientFrontier.getHighestRiskTolerancePortfolio(x);
};
MeanVarianceEfficientFrontierWrapper.prototype.getHighestRiskTolerance = function(x) {
	return this.efficientFrontier.getHighestRiskTolerance();
};
MeanVarianceEfficientFrontierWrapper.prototype.getLowestRiskTolerancePortfolio = function(x) {
	return this.efficientFrontier.getLowestRiskTolerancePortfolio(x);
};
MeanVarianceEfficientFrontierWrapper.prototype.getLowestRiskTolerance = function(x) {
	return this.efficientFrontier.getLowestRiskTolerance();
};


/**
* @function computeEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and
* long-only efficient portfolio belonging to the mean-variance efficient frontier, subject to 
* a return constraint, a volatility constraints or to a risk tolerance constraint.
*
* @memberof MeanVarianceEfficientFrontierWrapper
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
MeanVarianceEfficientFrontierWrapper.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	return this.efficientFrontier.computeEfficientPortfolio(constraintType, constraintValue);
};

/**
* @function getCornerPortfolios
*
* @description This function returns the weights w_i1,...,w_in associated to the m fully invested and
* long-only corner portfolios defining the mean-variance efficient frontier.
*
* @memberof MeanVarianceEfficientFrontierWrapper
*
* @return {Array<Array.<Object>>} the list of all corner portfolios, an array of arrays of n by 1 matrices
*
*/
MeanVarianceEfficientFrontierWrapper.prototype.getCornerPortfolios = function() {
	if (this.efficientFrontierOptimizationMethod == "critical-line") {
		return this.efficientFrontier.getCornerPortfolios();
	}
	else {
		throw new Error("internal error: unsupported method")
	}
};

/**
* @function restrict
*
* @description This function restricts the efficient frontier to portfolios satisfying 
* a minimal return constraint or a minimal volatility constraint.
*
* @memberof MeanVarianceEfficientFrontierWrapper
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "minReturn", to specify a return constraint
* - "minVolatility",  to specify a volatility constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
*/
MeanVarianceEfficientFrontierWrapper.prototype.restrict = function(constraintType, constraintValue) {
	return this.efficientFrontier.restrict(constraintType, constraintValue);
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
* @memberof MeanVarianceEfficientFrontierWrapper
*
* @param {number} rf the risk-free rate, a real number.
*
* @return {Array.<Object>} An array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontierWrapper.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	return this.efficientFrontier.computeMaximumSharpeRatioEfficientPortfolio(rf);
};