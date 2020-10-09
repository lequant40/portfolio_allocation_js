/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function MeanVarianceEfficientFrontier
*
* @description Virtual object representing a mean-variance efficient frontier, defined as
* the set of portfolios w satisfying the following optimization problem, for all values of
* the risk tolerance parameter rt, c.f. the references:
*
* min_w <sigma * w/w>/2 - rt*<w/mu>
* s.t. <w/mu> = 1 (full investment constraint) OR <w/mu> <= 1 (partial investment constraint)
*      l <= x <= u
*
* It holds:
* - Assets returns
* - Assets covariance matrix
* - Misc. constraints on the weights of the associated efficient portfolios (bounds constraints)
*
* It offers:
* - Methods to compute portfolios returns, volatilities, Sharpe ratios.
* - Methods to compute the highest/lowest attainable returns, volatilities and risk tolerances
* - Methods to compute efficient portfolios subject to a return, a volatility, a maximum volatility or a risk tolerance constraint
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the mean-variance optimization problem.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
*/
function MeanVarianceEfficientFrontier(mu, sigma, opt) {
	// Decode the options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	
	// Initializations
	//
	this.epsBounds = 1e-14; // the tolerance for checking bounds constraints for equality 
	this.epsPortfolioVolatility = 1e-14; // the tolerance for numerically zero volatility
	this.epsEfficientPortfolioComputation = 1e-10; // the tolerance for numerically searching for an efficient portfolio with a given constraint value
	
	//
	this.mu = new Matrix_(mu);
	
	//
	this.sigma = new Matrix_(sigma);
	
	//
	this.nbAssets = this.sigma.nbRows;
	
	this.lowerBounds = opt.constraints.minWeights == undefined ? Matrix_.zeros(this.nbAssets, 1) : new Matrix_(opt.constraints.minWeights);
	this.upperBounds = opt.constraints.maxWeights == undefined  ? Matrix_.ones(this.nbAssets, 1) : new Matrix_(opt.constraints.maxWeights);

	
	// ------	
	
	// Numerically alter tight minimum / maximum weights constraints
	// to have a feasible numerical problem.
	for (var i = 1; i <= this.nbAssets; ++i) {
		var lb_i = this.lowerBounds.getValue(i, 1);
		var ub_i = this.upperBounds.getValue(i, 1);
		
		if (lb_i == ub_i) {
			this.lowerBounds.setValueAt(i, 1, Math.max(0, lb_i - this.epsBounds));
			this.upperBounds.setValueAt(i, 1, Math.min(1, ub_i + this.epsBounds))
		}
	}

	// Check that the problem is feasible (i.e., that the restricted unit simplex on which
	// the optimization is taking place is not empty.
	simplexEmptinessCheck_(this.nbAssets, this.lowerBounds.toArray(), this.upperBounds.toArray());	
};

MeanVarianceEfficientFrontier.prototype.getHighestReturnPortfolio = function(x) {
	return this.getHighestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getHighestReturn = function(x) {
	return this.computePortfolioReturn(this.getHighestReturnPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getLowestReturnPortfolio = function(x) {
	return this.getLowestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getLowestReturn = function(x) {
	return this.computePortfolioReturn(this.getLowestReturnPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getHighestVolatilityPortfolio = function(x) {
	return this.getHighestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getHighestVolatility = function(x) {
	return this.computePortfolioVolatility(this.getHighestVolatilityPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getLowestVolatilityPortfolio = function(x) {
	return this.getLowestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getLowestVolatility = function(x) {
	return this.computePortfolioVolatility(this.getLowestVolatilityPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getHighestRiskTolerancePortfolio = function(x) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getHighestRiskTolerance = function(x) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getLowestRiskTolerancePortfolio = function(x) {
	throw new Error('internal error: unction is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getLowestRiskTolerance = function(x) {
	throw new Error('internal error: function is not implemented');
};

/**
* @function computePortfolioReturn
*
* @description This function computes the return of either:
* - A portfolio with weights x_1,...,x_n 
* - A portfolio with weights x_1,...,x_n, x_n+1, with x_n+1 representing the portion of cash in the portfolio
* using the assets returns vector associated to the efficient frontier.
*
* @memberof MeanVarianceEfficientFrontier
*
* @param {Matrix_} the portfolio weights, a n by 1 matrix or a n+1 by 1 matrix
*
* @return {number} the computed portfolio return
*/
MeanVarianceEfficientFrontier.prototype.computePortfolioReturn = function(x) {
	return Matrix_.vectorDotProduct(this.mu, x);
};

/**
* @function computePortfolioVolatility
*
* @description This function computes the volatility of either:
* - A portfolio with weights x_1,...,x_n 
* - A portfolio with weights x_1,...,x_n, x_n+1, with x_n+1 representing the portion of cash in the portfolio
* using the assets covariance matrix associated to the efficient frontier.
*
* @memberof MeanVarianceEfficientFrontier
*
* @param {Matrix_} the portfolio weights, a n by 1 matrix or a n+1 by 1 matrix
*
* @return {number} the computed portfolio volatility
*/
MeanVarianceEfficientFrontier.prototype.computePortfolioVolatility = function(x) {
	// Compute the variance x'*SIGMA*x.
	//
	// Take into account the partial investment constraint
	var sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, x), x);

	// In case the variance is numerically zero, which can occur with 
	// a semi-positive definite covariance matrix, it is replaced with zero.
	//
	// Otherwise, if the variance is negative, stops the algorithm,
	// since the covariance matrix is then not numerically semi-positive definite.
	if (Math.abs(sigma_x_x) <= this.epsPortfolioVolatility) {
		sigma_x_x = 0;
	}
	else if (sigma_x_x < 0 && sigma_x_x < -this.epsPortfolioVolatility) {
		throw new Error('internal error: negative volatility, covariance matrix might not be semi-definite positive');
	}
	
	// Compute the volatility SQRT(x'*SIGMA*x)
	var s_x = Math.sqrt(sigma_x_x);
   
	// Return the computed volatility
	return s_x;
};

/**
* @function computePortfolioSharpeRatio
*
* @description This function computes the Sharpe ratio of either:
* - A portfolio with weights x_1,...,x_n 
* - A portfolio with weights x_1,...,x_n, x_n+1, with x_n+1 representing the portion of cash in the portfolio
* using the assets returns vector and covariance matrix associated to the efficient frontier.
*
* @memberof MeanVarianceEfficientFrontier
*
* @param {Matrix_} the portfolio weights, a n by 1 matrix or a n+1 by 1 matrix
*
* @return {number} the computed portfolio Sharpe ratio
*/
MeanVarianceEfficientFrontier.prototype.computePortfolioSharpeRatio = function(x, rf) {
	// 
	if (rf == undefined || rf == null) {
		throw new Error('internal error: missing risk free rate');
	}
	
	// The numerator: <mu/x> - rf
	var ret = this.computePortfolioReturn(x);
	var excessRet = ret - rf;
	
	// The denominator: Sqrt(<Sigma*x/x>), which cannot be null 
	// with the way it is computed.
	var vol = Math.max(this.epsPortfolioVolatility, this.computePortfolioVolatility(x));
	if (vol == 0) {
		throw new Error('internal error: null volatility when computing the Sharpe ratio');
	}
	
	// Compute the Sharpe ratio
	var sharpeRatio = excessRet/vol;
	
	// Return the computed Sharpe ratio
	return sharpeRatio;
};

/**
* @function computeEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and
* long-only efficient portfolio belonging to the mean-variance efficient frontier, subject to 
* a return constraint, a volatility constraint or to a risk tolerance constraint.
*
* @memberof MeanVarianceEfficientFrontier
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
MeanVarianceEfficientFrontier.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	throw new Error('internal error: function is not implemented');
};

/**
* @function restrict
*
* @description This function restricts the efficient frontier to portfolios satisfying 
* a minimal return constraint, a minimal volatility constraint or a minimal risk tolerance constraint.
*
* @memberof MeanVarianceEfficientFrontier
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "minReturn", to specify a return constraint
* - "minVolatility",  to specify a volatility constraint
* - "minRiskTolerance", to specify a risk tolerance constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
*/
MeanVarianceEfficientFrontier.prototype.restrict = function(constraintType, constraintValue) {
	throw new Error('internal error: function is not implemented');
};

/**
* @function computeMaximumSharpeRatioEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only portfolio of n assets maximizing the Sharpe ratio, which is defined as the ratio of the 
* portfolio excess return over a constant risk-free rate to the portfolio volatility.
*
* @memberof MeanVarianceEfficientFrontier
*
* @param {number} rf the risk-free rate, a real number.
*
* @return {Array.<Object>} An array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontier.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	throw new Error('internal error: function is not implemented');
};