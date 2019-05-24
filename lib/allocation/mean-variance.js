/**
 * @file Functions related to mean variance efficient portfolios.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.computeCornerPortfolios_ = computeCornerPortfolios_;
/* End Wrapper private methods - Unit tests usage only */



/**
* @function meanVarianceOptimizationWeights
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a return
* constraint or to misc. volatility constraints.
*
* @description This function returns the weights w_1,...,w_n associated to the 
* long-only mean-variance efficient portfolio of n assets subject to either:
* - a return constraint (if optimizationMethod is 'targetReturn'), in which case this portfolio, if it exists, is fully invested 
* and has the lowest attainable volatility among all the feasible portfolios satisfying the return constraint
* - a volatility constraint (if optimizationMethod is 'targetVolatility'), in which case this portfolio, if it exists, is fully invested
* and has the highest attainable return among all the feasible portfolios satisfying the volatility constraint
* - a maximum volatility constraint (if optimizationMethod is 'maximumTargetVolatility'), in which case this portfolio is potentially 
* not fully invested and has the highest attainable return among all the feasible portfolios satisfying the maximum volatility constraint
* - a minimum variance constraint (if optimizationMethod is 'minimumVariance'), in which case this portfolio is fully invested and 
* has the lowest attainable volatility among all the feasible portfolios
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset to include in the portfolio
* - Maximum weight of each asset to include in the portfolio
*
* The main algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt parameters for the mean-variance optimization algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.optimizationMethod the mandatory mean-variance optimization algorithm to use, a string either equal to:
* - 'targetReturn', to compute the mean-variance efficient portfolio subject to a return constraint
* - 'targetVolatility', to compute the mean-variance efficient portfolio subject to a volatility constraint
* - 'maximumTargetVolatility', to compute the mean-variance efficient portfolio subject to a maximum volatility constraint
* - 'minimumVariance', to compute the global minimum variance efficient portfolio
* @param {number} opt.constraints.return in case opt.optimizationMethod is equal to 'targetReturn', the target return of the portfolio, a real number.
* @param {number} opt.constraints.volatility in case opt.optimizationMethod is equal to 'targetVolatility', the target volatility of the portfolio, a positive real number.
* @param {number} opt.constraints.maxVolatility in case opt.optimizationMethod is equal to 'maximumTargetVolatility', the maximum target volatility of the portfolio, a positive real number.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights corresponding to the mean-variance efficient portfolio, array of n real numbers.
*
* @example
* meanVarianceOptimizationWeights([0.1, 0.2], [[1, 0.3], [0.3, 1]], { optimizationMethod: 'targetReturn', constraints: {return: 0.15}})
* // [0.5, 0.5] 
*/
self.meanVarianceOptimizationWeights = function(mu, sigma, opt) {	
	// ------	

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	
	// Decode the optimization method
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		throw new Error('missing mean-variance optimization method');
	}
	if (optimizationMethod !== 'targetReturn' &&  
	    optimizationMethod !== 'targetVolatility' &&
		optimizationMethod !== 'minimumVariance' &&
		optimizationMethod !== 'maximumTargetVolatility') {
		throw new Error('unsupported mean-variance optimization method');
	}
	
	// Decode the optimization method parameters
	var targetReturn = opt.constraints["return"]; // .return does not work within Google Script
	if (optimizationMethod === 'targetReturn' && targetReturn === undefined) {
		throw new Error('missing mean-variance optimization method return constraint');
	}
	
	var targetVolatility = opt.constraints.volatility;
	if (optimizationMethod === 'targetVolatility' && targetVolatility === undefined) {
		throw new Error('missing mean-variance optimization method volatility constraint');
	}

	var maxTargetVolatility = opt.constraints.maxVolatility;
	if (optimizationMethod === 'maximumTargetVolatility' &&  maxTargetVolatility === undefined) {
		throw new Error('missing mean-variance optimization method maximum volatility constraint');
	}
	
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	
	
	// ------
	
	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	
	// Depending on the target function, proceed with a different algorithm
	// to compute the requested efficient portfolio.
	var efficientPortfolio;
	if (optimizationMethod == 'targetReturn') {
		efficientPortfolio = computeTargetReturnEfficientPortfolio_(mu, targetReturn, cornerPortfolios);
		if (efficientPortfolio.length == 0) {
			throw new Error('return not reachable');
		}
		else {
			efficientPortfolio = efficientPortfolio[0];
		}
	}
	else if (optimizationMethod == 'targetVolatility') {
		efficientPortfolio = computeTargetVolatilityEfficientPortfolio_(sigma, targetVolatility, cornerPortfolios);
		if (efficientPortfolio.length == 0) {
			throw new Error('volatility not reachable');
		}
		else {
			efficientPortfolio = efficientPortfolio[0];
		}
	}
	else if (optimizationMethod == 'minimumVariance') {
		efficientPortfolio = computeMinimumVarianceEfficientPortfolio_(cornerPortfolios);
	}
	else if (optimizationMethod == 'maximumTargetVolatility') {
		// The options for the potential internal mean variance optimization algorithm
		var opt_mv = { maxIter: opt.maxIter, constraints: opt.constraints };
		
		efficientPortfolio = computeMaximumTargetVolatilityEfficientPortfolio_(mu, sigma, maxTargetVolatility, cornerPortfolios, opt_mv);
	}
	else {
		throw new Error('internal error');
	}
	
	// Return the computed portfolio weights
	var weights = efficientPortfolio;
	return weights.toArray();
}



/**
* @function maximumSharpeRatioWeights
*
* @summary Compute the weights of the portfolio maximizing the Sharpe ratio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only portfolio of n assets maximizing the Sharpe ratio, which is defined as the ratio of the 
* portfolio excess return over a constant risk-free rate to the portfolio volatility.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset to include in the portfolio
* - Maximum weight of each asset to include in the portfolio
*
* When it exists, this portfolio is mean-variance efficient and is unique.
*
* The main algorithm used internally is the Markowitz critical line algorithm, c.f. the second reference.
*
* @see <a href="https://doi.org/10.1111/j.1540-6261.1976.tb03217.x">Elton, E. J., Gruber, M. J. and Padberg, M. W. (1976), SIMPLE CRITERIA FOR OPTIMAL PORTFOLIO SELECTION. The Journal of Finance, 31: 1341-1357</a>
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see <a href="https://www.jstor.org/stable/1924119">Lintner, John. The Valuation of Risk Assets and the Selection of Risky Investments in Stock Portfolios and Capital Budgets. The Review of Economics and Statistics 47, no. 1 (1965): 13-37.</a>
*
* @param {<Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {number} rf the risk-free rate, a real number.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array<Array.<number>} the weights corresponding to the portfolio maximizing the Sharpe ratio, array of n real numbers.
*
* @example
* maximumSharpeRatioWeights([0.1, 0.2], [[1, 0.3], [0.3, 1]], 0)
* // [~0.19, ~0.81]
*/
self.maximumSharpeRatioWeights = function(mu, sigma, rf, opt) {
	// Internal function to compute the return of a portfolio
	function computeReturn_(mu, weights) {
		return Matrix_.vectorDotProduct(mu, weights);
	}

	// Internal function to compute the volatility of a portfolio
	function computeVolatility_(sigma, weights) {
		return Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights), weights));
	}
	
	// Internal function to compute the Sharpe ratio of a portfolio,
	// as well as the intermediate values of return and volatility.
	function computeSharpeRatio_(mu, sigma, rf, weights) {
		var eps = 1e-8; // the numerical zero
		
		// The numerator: <mu/w> - rf
		var ret = computeReturn_(mu, weights);
		var excessRet = ret - rf;
		
		// The denominator: Sqrt(<Sigma*w/w>)
		var vol = computeVolatility_(sigma, weights);
		
		// In case the denominator is numerically null, which can occur with
		// semi-positive definite covariance matrices, replace it with the
		// value of the numerical zero.
		if (Math.abs(vol) < eps) {
			vol = eps;
		}
		
		// Compute the Sharpe ratio
		var sharpeRatio = excessRet/vol;
		
		// Return the computed Sharpe ratio and intermediate values
		return [sharpeRatio, ret, vol];
	}

	// Internal function to compute the corner portfolio which maximizes the
	// Sharpe ratio on the efficient frontier restricted to portfolios with:
	// - A strictly positive volatility
	// - A strictly positive excess return
	//
	// This function uses a binary search algorithm, which is justified because
	// the Sharpe ratio is a strictly unimodular function on the above restricted
	// efficient frontier, c.f. the first reference.
	function computeMaximumSharpeRatioCornerPortfolio_(mu, sigma, rf, cornerPortfolios) {
		var eps = 1e-8; // the numerical zero
		
		// The efficient frontier portfolios are provided from highest return/volatility
		// to lowest return/volatility, so that *_min below refers to properties of the portfolio
		// with the lowest return/volatility.
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0;

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];

		var sharpeRatio_min = computeSharpeRatio_(mu, sigma, rf, weights_min);
		var sharpeRatio_max = computeSharpeRatio_(mu, sigma, rf, weights_max);
		
		// In case there is only one corner portfolio on the efficient frontier,
		// exit immediately.
		if (idx_min == idx_max) {
			return [idx_min, weights_min, sharpeRatio_min];
		}
		
		// Otherwise, determine the corner portfolio with the maximum Sharpe ratio 
		// using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle points
			var idx_middle = Math.floor(idx_max + (idx_min - idx_max)/2); // formula avoiding numerical overflow
			var weights_middle = cornerPortfolios[idx_middle][0];
			var sharpeRatio_middle = computeSharpeRatio_(mu, sigma, rf, weights_middle);

			var idx_middle_p = idx_middle + 1; 
			var weights_middle_p = cornerPortfolios[idx_middle_p][0];
			var sharpeRatio_middle_p = computeSharpeRatio_(mu, sigma, rf, weights_middle_p);

			// Determine in which sub-interval [idx_max, idx_middle+1] or [idx_middle+1, idx_min]
			// lies the corner portfolio with the maximum Sharpe ratio.
			if (sharpeRatio_middle[0] > sharpeRatio_middle_p[0]) {
				idx_min = idx_middle;
				weights_min = weights_middle;
				sharpeRatio_min = sharpeRatio_middle;
			}
			else if (sharpeRatio_middle[0]	< sharpeRatio_middle_p[0]) {
				idx_max = idx_middle;
				weights_max = weights_middle;
				sharpeRatio_max = sharpeRatio_middle;
			}
			else {
				// In case the Sharpe ratio is equal on both corner portfolios, 
				// it means its maximum is attained somewhere between these two portfolios, 
				// due to its strict unimodality.
				//
				// The binary search procedure can then be prematurely stopped, although
				// this case is (numerically) highly improbable.
				idx_min = idx_middle_p;
				weights_min = weights_middle_p;
				sharpeRatio_min = sharpeRatio_middle_p;
				
				idx_max = idx_middle;
				weights_max = weights_middle;
				sharpeRatio_max = sharpeRatio_middle;

				break;
			}
		}

		
		// Return the computed corner portfolio
		return [idx_min, weights_min, sharpeRatio_min];
	}

	// Internal function to compute the efficient portfolio which maximizes the
	// Sharpe ratio on an efficient segment defined by two adjacent corner portfolios
	// with:
	// - A strictly positive volatility
	// - A strictly positive excess return
	//
	// On such an efficient segment, the weights associated this portfolio are a 
	// convex combination of the weights of the two adjacent corner portfolios,
	// so that w = t*w_min + (1-t)*w_max, t in [0,1], with t to be determined,
	// c.f. the second reference.
	//
	// With E(w) = <mu/w> the portfolio return and V(w) = <Sigma*w/w> the portfolio 
	// variance, the Sharpe ratio is defined as SR(w) = (E(w) - rf)/SQRT(V(w)).
	//
	// Because SR(w) > 0 on the efficient segment, maximizing SR(w) is equivalent
	// to maximizing SR(w)^2, which is equal to (E(w) - rf)^2/V(w).
	//
	// By linearity of E(w) and bilinearity/symmetry of V(w), SR(w)^2 is also equal to
	// a rational fraction in t:
	//
	// (E(w) - rf)^2/V(w)
	// =
	// ( E(t*w_min + (1-t)*w_max) - rf )^2 / ( V(t*w_min + (1-t)*w_max) )
	// =
	// ( t*(E(w_min) - E(w_max)) + E(w_max) - rf )^2 / ( t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) - 2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) )
	// = ( t^2*(E(w_min) - E(w_max))^2 + 2*(E(w_min) - E(w_max))*(E(w_max) - rf) + (E(w_max) - rf)^2 ) / ( t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) - 2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) )
	//
	// So, maximizing SR(w) on the efficient segment is equivalent to maximizing
	// SR(t)^2, t in [0,1].
	//
	// Since SR(t)^2 is a differentiable function on [0,1] and since [0,1] is a closed convex set,
	// its maximum is either reached on its boundary (i.e., {0,1}) or on a critical interior point
	// (i.e., a point belonging to ]0,1[ on which the derivative of SR(t)^2 vanishes).
	//
	// Evaluating SR(t) on each of these (at most) four points and selecting t
	// as the value which maximizes SR(t) then allows to compute the weights 
	// of the efficient portfolio which maximizes the Sharpe ratio on the efficient segment.
	function computeLocalMaximumSharpeRatioEfficientPortfolio_(mu, sigma, rf, idx_min, weights_min, idx_max, weights_max) {
		// Compute properties of the two adjacent corner portfolios
		var sharpeRatio_min = computeSharpeRatio_(mu, sigma, rf, weights_min);
		var sr_min = sharpeRatio_min[0];
		var return_min = sharpeRatio_min[1];
		var volatility_min = sharpeRatio_min[2];
		var variance_min = volatility_min * volatility_min;
		
		var sharpeRatio_max = computeSharpeRatio_(mu, sigma, rf, weights_max);
		var sr_max = sharpeRatio_max[0];
		var return_max = sharpeRatio_max[1];
		var volatility_max = sharpeRatio_max[2];
		var variance_max = volatility_max * volatility_max;
		
		// Define the coefficients of the fractional function SR(t)^2 = ( at^2 + bt + c ) / ( dt^2 + et + f )
		var return_min_m_max = return_min - return_max;
		var return_max_m_rf = return_max - rf;
		var a = return_min_m_max * return_min_m_max;
		var b = 2 * return_min_m_max * return_max_m_rf;
		var c = return_max_m_rf * return_max_m_rf;
		
		var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
		var d = variance_min + variance_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
		var e = -2 * (variance_max - variance_cross); // 
		var f = variance_max; //always > 0
		
		// Define the coefficients of the second order polynomial aat^2 + bbt + cc equal to the
		// numerator of the derivative d(SR(t)^2)/dt.
		var aa = a*e - b*d;
		var bb = 2*(a*f - c*d);
		var cc = b*f - c*e;
		
		// Extract the roots t1 and t2 of the equation d(SR(t)^2)/dt = 0, using a stable numerical formula.
		var bb_p = bb/2; // reduced discriminant
		var sign_bb_p = (bb_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for bb_p == 0 this returns 1
		var disc = bb_p*bb_p - aa*cc;
		if (disc < 0) {
			throw new Error('internal error, the covariance matrix might not be semi-definite positive');
		}
		var qq = -(bb_p + sign_bb_p * Math.sqrt(disc));
		var t1 = qq/aa;
		var t2 = cc/qq;
		
		// Compute and order the Sharpe ratios for all the efficient 
		// portfolios with t corresponding to {0, 1, t1, t2}.
		var candidateSharpeRatios = [[weights_min, sr_min], [weights_max, sr_max]]; // t = 0 and t = 1 portfolios are always present
		
		if (t1 > 0 && t1 < 1) { // t1 belongs to ]0,1[
			var weights_t1 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t1*weights_min.getValue(i, 1) + (1-t1)*weights_max.getValue(i, 1); 
										})
			var sharpeRatio_t1 = computeSharpeRatio_(mu, sigma, rf, weights_t1);
			var sr_t1 = sharpeRatio_t1[0];
			
			candidateSharpeRatios.push([weights_t1, sr_t1]);
		}

		if (t2 > 0 && t2 < 1) { // t2 belongs to ]0,1[
			var weights_t2 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t2*weights_min.getValue(i, 1) + (1-t2)*weights_max.getValue(i, 1); 
										})
			var sharpeRatio_t2 = computeSharpeRatio_(mu, sigma, rf, weights_t2);
			var sr_t2 = sharpeRatio_t2[0];

			candidateSharpeRatios.push([weights_t2, sr_t2]);
		}

		// Return the efficient portfolio which maximizes the Sharpe ratio
		// on the efficient segment.
		var compareSharpeRatios = function (a, b) {
			return a[1] - b[1];
		};
		return max_(candidateSharpeRatios, compareSharpeRatios)[0];
	}
	
	
	// ------	

	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	
	
	// ------
	
	// Initializations
	var eps = 1e-8; // the numerical zero

	
	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	
	// Retrict the efficient frontier to the domain of definition
	// of the Sharpe ratio by determining the "first" efficient portfolio 
	// with a strictly positive volatility.
	//
	// To be noted that:
	// - In case the covariance matrix is positive definite, this is the minimum 
	// variance portfolio, so that the efficient frontier is not altered
	//
	// - In case the covariance matrix is semi-positive definite, this is an efficient
	// portfolio located "close" to the minimum variance portfolio, because the 
	// variance of corner portfolios is strictly increasing	
	var idx = cornerPortfolios.length - 1;
	var weights = cornerPortfolios[idx][0];
	var volatility = computeVolatility_(sigma, weights);
	if (volatility < eps) {
		// The domain of definition of the Sharpe ratio is not the
		// whole efficient frontier, so that the efficient frontier
		// needs to be restricted.
		
		// Compute the first efficient portfolio with a strictly positive 
		// volatility.
		var efficientPortfolio = computeTargetVolatilityEfficientPortfolio_(sigma, eps, cornerPortfolios);
		if (efficientPortfolio.length == 0) {
			throw new Error('no corner portfolio with a strictly positive volatility: the covariance matrix might not be semi-definite positive');
		}
		var efficientPortfolioWeights = efficientPortfolio[0];
		var cornerPortfolioIndexMin = efficientPortfolio[1];
			
		// Add it as a replacement of the last corner
		// portfolio with a volatility lower than or equal to it.
		//
		// A risk aversion parameter of -1 is set so that subsequent processes 
		// can distinguish between an original corner portfolio and a 
		// replacement corner portfolio.
		cornerPortfolios[cornerPortfolioIndexMin] = [efficientPortfolioWeights, -1];
		
		// Remove the corner portfolios with a volatility strictly lower than
		// the it.
		cornerPortfolios.length = cornerPortfolioIndexMin + 1;
	}

	
	// Further retrict the efficient frontier to the domain of strict positivity
	// of the Sharpe ratio.
	//
	// To be noted that the domain of strict positivity of the Sharpe ratio
	// can be empty in case there is no feasible portfolio on the efficient
	// frontier with a strictly positive excess return.
	var idx = cornerPortfolios.length - 1;
	var weights = cornerPortfolios[idx][0];
	var ret = computeReturn_(mu, weights);
	if (ret < rf + eps) {
		// The domain of strict positivity of the Sharpe ratio is not the
		// whole efficient frontier, so that the efficient frontier
		// needs to be restricted.
		
		// Compute the first efficient portfolio with a strictly positive 
		// excess return.
		var efficientPortfolio = computeTargetReturnEfficientPortfolio_(mu, rf + eps, cornerPortfolios);
		if (efficientPortfolio.length == 0) {
			throw new Error('no corner portfolio with a strictly positive excess return');
		}
		var efficientPortfolioWeights = efficientPortfolio[0];
		var cornerPortfolioIndexMin = efficientPortfolio[1];
		
		// Add it as a replacement of the last corner
		// portfolio with an excess return lower than or equal to it.
		//
		// A risk aversion parameter of -1 is set so that subsequent processes 
		// can distinguish between an original corner portfolio and a 
		// replacement corner portfolio.
		cornerPortfolios[cornerPortfolioIndexMin] = [efficientPortfolioWeights, -1];
		
		// Remove the corner portfolios with an excess return strictly lower than it.
		cornerPortfolios.length = cornerPortfolioIndexMin + 1;
	}


	// On the retricted efficient frontier, the Sharpe ratio is a pseudo-concave 
	// function, c.f. the first reference, so that it is strictly unimodal.
	//
	// This property allows to search for the corner portfolio with the maximum
	// Sharpe ratio using a binary search algorithm.
	var cornerPortfolio = computeMaximumSharpeRatioCornerPortfolio_(mu, sigma, rf, cornerPortfolios);
	var idx_middle = cornerPortfolio[0];
	var weights_middle = cornerPortfolio[1];
	var sr_middle = cornerPortfolio[2][0];
	

	// The corner portfolio with the maximum Sharpe ratio is adjacent to at most
	// two other corner portfolios, depending on its position on the efficient
	// frontier:
	// - Unique portfolio => zero adjacent corner portfolio
	// - Non unique leftmost or rightmost corner portfolio => one adjacent corner portfolio
	// - Non unique any other corner portfolio => two adjacent corner portfolios
	//
	// In the first case, the efficient portfolio with the maximum Sharpe ratio
	// is the same as the corner portfolio with the maximum Sharpe ratio
	//
	// In the last two cases, because of the strict unimodality of the Sharpe ratio
	// on the restricted efficient frontier, the efficient portfolio with the maximum
	// Sharpe ratio is guaranteed to belong to the efficient segment(s) connecting
	// the corner portfolio with the maximum Sharpe ratio to its adjacent corner
	// portfolio(s).
	//
	// So, computing the efficient portfolio with the maximum Sharpe ratio is equivalent
	// to computing the efficient portfolio with the maximum Sharpe ratio on the efficient
	// segment(s) connecting the corner portfolio with the maximum Sharpe ratio
	// to its adjacent corner portfolio(s).
	
	// Add the corner portfolio with the maximum Sharpe ratio as a candidate 
	// for being the efficient portfolio with the maximum Sharpe ratio.
	var candidateSharpeRatios = [[weights_middle, sr_middle]];
	
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// on the efficient segment [idx_middle+1, idx_middle], if existing.
	var idx_min = idx_middle + 1;
	if (idx_min <= cornerPortfolios.length - 1) {
		var weights_min = cornerPortfolios[idx_min][0];
		
		var weights_min_middle = computeLocalMaximumSharpeRatioEfficientPortfolio_(mu, sigma, rf, 
		                                                                           idx_min, weights_min, 
																			       idx_middle, weights_middle);
		candidateSharpeRatios.push(weights_min_middle);
	}
	
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// on the efficient segment [idx_middle, idx_middle-1], if existing.	
	var idx_max = idx_middle - 1;
	if (idx_max >= 0) {
		var weights_max = cornerPortfolios[idx_max][0];
		
		var weights_middle_max = computeLocalMaximumSharpeRatioEfficientPortfolio_(mu, sigma, rf, 
		                                                                           idx_middle, weights_middle, 
																			       idx_max, weights_max);
        candidateSharpeRatios.push(weights_middle_max);
	}
	
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// by merging the efficient portfolios locally maximizing
	// the Sharpe ratio on each efficient segment.
	var compareSharpeRatios = function (a, b) {
		return a[1] - b[1];
	};
	var maxSharpeRatio = max_(candidateSharpeRatios, compareSharpeRatios)[0];
	var weights = maxSharpeRatio[0];


	// Return the computed portfolio weights
	return weights.toArray();
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
* - Minimum weight of each asset to include in the portfolios
* - Maximum weight of each asset to include in the portfolios
*
* The main algorithm used internally is the Markowitz critical line algorithm, c.f. the reference.
*
* The algorithm used internally generates the portfolios uniformly on the efficient frontier with
* regard to the interval of variation of the risk aversion parameter.
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
			var portfolioWeights = cornerPortfolios[0][0];
			var portfolioReturn = Matrix_.vectorDotProduct(mu, portfolioWeights);
			var portfolioVolatility = Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, portfolioWeights), 
																	     portfolioWeights));
			
			return [[portfolioWeights.toArray(), portfolioReturn, portfolioVolatility]];
		}
	}
	else { // cornerPortfolios.length >= 2
		if (nbPortfolios <= 1) {
			throw new Error('efficient frontier made of several corner portfolios: at least two efficient portfolios must be computed');
		}
	}
	
	// Generate nbPortfolios distinct points lambda_i, i=1..nbPortfolios, corresponding to
	// strictly increasing values of the risk aversion parameter, uniformly
	// spaced on the interval [lambda_min, lambda_max], using the formula
	// lambda_i = lambda_min + i * (lambda_max - lambda_min)/(nbPortfolios - 1).
	//
	// Then, for each of these points, compute the two enclosing corner portfolios 
	// w_i_min, w_i_max satisfying lambda_i_min <= lambda_i < lambda_i_max or
	// lambda_i_min < lambda_i <= lambda_i_max.
	//
	// In this case, the weights corresponding to the associated efficient
	// portfolio are a convex combination of the weights of the two computed enclosing
	// corner portfolios (c.f. the reference): w_i = t*w_i_min + (1-t)*w_i_max, t in [0,1],
	// with t now to be determined.
	//
	// As the relationship between lambda_i and and w_i is the identity, we have
	// lambda_i = (1-t)*lambda_i_min + t*lambda_i_max
	// <=>
	// t = (lambda_i - lambda_i_min)/(lambda_i_max - lambda_i_min)
	
	// Initializations
	var lambda_min = cornerPortfolios[cornerPortfolios.length - 1][1];
	var lambda_max = cornerPortfolios[0][1];
	var delta_lambda = (lambda_max - lambda_min)/(nbPortfolios - 1);
	
	var lambda_i = lambda_min; // the first lambda_i point is lambda_min
	var lambda_i_min_idx = cornerPortfolios.length - 1;
	var lambda_i_max_idx = lambda_i_min_idx - 1;
	var lambda_i_min = cornerPortfolios[lambda_i_min_idx][1];
	var lambda_i_max = cornerPortfolios[lambda_i_max_idx][1];
	var w_i_min = cornerPortfolios[lambda_i_min_idx][0];
	var w_i_max = cornerPortfolios[lambda_i_max_idx][0];
	
	// Specific process for the first efficient portfolio (the
	// minimum variance portfolio).
	{
		var minimumVariancePortfolioWeights = cornerPortfolios[lambda_i_min_idx][0];
		var portfolioReturn = Matrix_.vectorDotProduct(mu, minimumVariancePortfolioWeights);
		var portfolioVolatility = Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, minimumVariancePortfolioWeights), 
																	minimumVariancePortfolioWeights));
		
		efficientFrontier[0] = [minimumVariancePortfolioWeights.toArray(), portfolioReturn, portfolioVolatility];
	}
	
	
	// Core process for the nbPortfolios-2 middle efficient portfolios
	for (var i = 1; i < nbPortfolios - 1; ++i) {
		// Generate the current risk aversion point
		var lambda_i = lambda_min + i * delta_lambda;
		
		// Compute the two enclosing corner portfolios 
		//
		// Note: the associated indexes and values are updated
		// only when the current risk aversion point goes beyond
		// the current [lambda_i_min, lambda_i_max] interval
		while (lambda_i > lambda_i_max) {
			--lambda_i_min_idx;
			--lambda_i_max_idx;
			
			lambda_i_min = cornerPortfolios[lambda_i_min_idx][1];
			lambda_i_max = cornerPortfolios[lambda_i_max_idx][1];
			
			w_i_min = cornerPortfolios[lambda_i_min_idx][0];
			w_i_max = cornerPortfolios[lambda_i_max_idx][0];
		}
				
		// Compute the efficient portfolios weights, returns and volatilities
		var t = (lambda_i - lambda_i_min)/(lambda_i_max - lambda_i_min);
		var portfolioWeights = Matrix_.fill(nbAssets, 1, 
								function(i,j) { 
									return (1-t)*w_i_min.getValue(i, 1) + t*w_i_max.getValue(i, 1); 
								})
		var portfolioReturn = Matrix_.vectorDotProduct(mu, portfolioWeights);
		var portfolioVolatility = Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, portfolioWeights), 
																	   portfolioWeights));
		
		efficientFrontier[i] = [portfolioWeights.toArray(), portfolioReturn, portfolioVolatility];
	}
	
	
	// Specific process for the last efficient portfolio (the
	// maximum return portfolio).
	lambda_i = lambda_max;
	lambda_i_min_idx = 1;
	lambda_i_max_idx = lambda_i_min_idx - 1;
	{
		var maximumReturnPortfolioWeights = cornerPortfolios[lambda_i_max_idx][0];
		var portfolioReturn = Matrix_.vectorDotProduct(mu, maximumReturnPortfolioWeights);
		var portfolioVolatility = Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, maximumReturnPortfolioWeights), 
																	 maximumReturnPortfolioWeights));
		
		efficientFrontier[nbPortfolios - 1] = [maximumReturnPortfolioWeights.toArray(), portfolioReturn, portfolioVolatility];
	}
	
	
	// Return the computed list of portfolios weights, returns and volatilities
	return efficientFrontier;
}


/**
* @function meanVarianceCornerPortfolios
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
* - Minimum weight of each asset to include in the portfolios
* - Maximum weight of each asset to include in the portfolios
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
* meanVarianceCornerPortfolios([0.1, 0.2], [[1, 0.3], [0.3, 1]])
* // [[[0.5, 0.5], ~0.15, ~0.806], [[0, 1], 0.2, 1]]
*/
self.meanVarianceCornerPortfolios = function(mu, sigma, opt) {
	// Convert mu and sigma to matrix format
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);

	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	// Initializations
	var efficientFrontier = new Array(cornerPortfolios.length);
	
	// Convert the output of the internal function above to a list of 
	// portfolios weights, returns and volatilities.
	for (var i = 0; i < cornerPortfolios.length; ++i) {
		var portfolioWeights = cornerPortfolios[i][0];
		var portfolioReturn = Matrix_.vectorDotProduct(mu, portfolioWeights);
		var portfolioVolatility = Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, portfolioWeights), 
																	       portfolioWeights));
		
		efficientFrontier[cornerPortfolios.length - 1 - i] = [portfolioWeights.toArray(), portfolioReturn, portfolioVolatility];
	}
	
	// Return the computed list of portfolios weights, returns and volatilities
	return efficientFrontier;
}



/**
* @function computeTargetReturnEfficientPortfolio_
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a target return
* constraint, as well as the index(es) of its enclosing corner portfolio(s) on the efficient frontier.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only mean-variance efficient portfolio of n assets subject to a target return constraint, 
* as well as the index(es) of its enclosing corner portfolio(s) on the efficient frontier.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Matrix_} mu the returns of the n assets in the considered universe, a n by 1 matrix of real numbers.
* @param {number} targetReturn the desired return of the portfolio, a real number.
* @param {Array<Array.<Object>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
* @return {Array.<Object>} an array arr of either two or three elements:
* - arr[0][0], the weights of the efficient portfolio with the desired return, a n by 1 Matrix_ of n real numbers
* - arr[0][1], the index of the corner portfolio with a return equal to (in case arr is made of two elements) or strictly lower than (in case arr is made of three elements)
* the return of the efficient portfolio, a natural integer
* - arr[0][2], the index of the corner portfolio with a return strictly greater than the return of the efficient portfolio (only in case arr is made of three elements), a natural integer
*/
function computeTargetReturnEfficientPortfolio_(mu, targetReturn, cornerPortfolios) {
	// Internal functon to compute the return of a portfolio
	function computeReturn_(mu, weights) {
		return  Matrix_.vectorDotProduct(mu, weights);
	}

	// Internal function to compute the (at most) two corner portfolios strictly enclosing the
	// efficient portfolio with a target return, using a binary search algorithm.
	//
	// The usage of a binary search algorithm is justified because the corner portfolios
	// return is strictly decreasing as soon as there are at least two corner
	// portfolios on the efficient frontier.
	function computeEnclosingCornerPortfolios_(targetReturn, mu, cornerPortfolios) {
		var eps = 1e-8; // the numerical zero
		
		// The efficient frontier portfolios are provided from highest return
		// to lowest return, so that *_min below refers to properties of the portfolio
		// with the lowest return.
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];

		var return_min = computeReturn_(mu, weights_min);
		var return_max = computeReturn_(mu, weights_max);

		// If the target return is not reachable within numerical accuracy, 
		// return immediately.
		if (targetReturn - return_max > eps || -eps > targetReturn - return_min) {
			return [];
		}
		
		// If the target return is numerically reached on one of the
		// two extremal corner portfolios, return immediately.
		if (Math.abs(targetReturn - return_min) <= eps) {
			return [[idx_min, weights_min, return_min]];
		}
		else if (Math.abs(targetReturn - return_max) <= eps) {
			return [[idx_max, weights_max, return_max]];
		}
		
		// Otherwise, determine the two adjacent corner portfolios strictly enclosing the portfolio
		// with a return numerically equals to the target return, using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle point
			var idx_middle = Math.floor(idx_max + (idx_min - idx_max)/2); // formula avoiding numerical overflow
			var weights_middle = cornerPortfolios[idx_middle][0];
			var return_middle = computeReturn_(mu, weights_middle);
			
			// Determine in which sub-interval ]idx_max, idx_middle[ or ]idx_middle, idx_min[
			// lies the portfolio with the target return.
			if (return_middle - targetReturn > eps) {
				idx_max = idx_middle;
				return_max = return_middle;
				weights_max = weights_middle;
			}
			else if (return_middle - targetReturn < -eps) {
				idx_min = idx_middle;
				return_min = return_middle;
				weights_min = weights_middle;
			}
			else { // the target return is exactly attained on the idx_middle-th corner portfolio
				return [[idx_middle, weights_middle, return_middle]];
			}
		}

		
		// Return the computed adjacent corner portfolios, as well as
		// the associated function values.
		return [[idx_min, weights_min, return_min], [idx_max, weights_max, return_max]];
	}
	
	
	// ------	

	
	// Compute the (at most) two corner portfolios strictly enclosing the efficient 
	// portfolio with a return equals to the target return.
	var enclosingCornerPortfolios = computeEnclosingCornerPortfolios_(targetReturn, mu, cornerPortfolios);

	
	// Then:
	// - In case the desired target return is not reachable, return an empty portfolio 
	//
	// - In case there is a unique computed corner portfolio with a return
	// equals to the target return, return the associated portfolio weights
	//
	// - In case there are two corner portfolios strictly enclosing the efficient portfolio with 
	// a return equals to the target return, the weights associated to this efficient portfolio are
	// a (strict) convex combination of the weights of the two computed enclosing corner portfolios
	// (c.f. the reference): w = t*w_min + (1-t)*w_max, t in ]0,1[, with t now to be determined.
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
		var return_min = enclosingCornerPortfolios[0][2];
		
		var idx_max = enclosingCornerPortfolios[1][0];
		var weights_max = enclosingCornerPortfolios[1][1];
		var return_max = enclosingCornerPortfolios[1][2];
		
		// The procedure to compute t above is the following:
		// E(w) = <mu/w> and by linearity of E, we have
		// E(w) = t*E(w_min) + (1-t)*E(w_max) and E(w) = targetReturn
		// <=>
		// t = (E(w_max) - targetReturn)/(E(w_max) - E(w_min))
		var t = (return_max - targetReturn)/(return_max - return_min);

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
* @function computeTargetVolatilityEfficientPortfolio_
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a target volatility
* constraint, as well as the index(es) of its enclosing corner portfolio(s) on the efficient frontier.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only mean-variance efficient portfolio of n assets subject to a target volatility constraint, 
* as well as the index(es) of its enclosing corner portfolio(s) on the efficient frontier.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Matrix_} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, a n by n matrix of real numbers.
* @param {number} targetVolatility the desired volatility of the portfolio, a positive real number.
* @param {Array<Array.<Object>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
* @return {Array.<Object>} an array arr of either two or three elements:
* - arr[0][0], the weights of the efficient portfolio with the desired volatility, a n by 1 Matrix_ of n real numbers
* - arr[0][1], the index of the corner portfolio with a volatility equal to (in case arr is made of two elements) or strictly lower than (in case arr is made of three elements)
* the volatility of the efficient portfolio, a natural integer
* - arr[0][2], the index of the corner portfolio with a volatility strictly greater than the volatility of the efficient portfolio (only in case arr is made of three elements), a natural integer
*/
function computeTargetVolatilityEfficientPortfolio_(sigma, targetVolatility, cornerPortfolios) {
	// Internal function to compute the volatility of a portfolio
	function computeVolatility_(sigma, weights) {
		return Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights), weights));
	}
	
	// Internal function to compute the (at most) two corner portfolios strictly enclosing the
	// efficient portfolio with a target volatility, using a binary search algorithm.
	//
	// The usage of a binary search algorithm is justified because the corner portfolios
	// volatility is strictly decreasing as soon as there are at least two corner
	// portfolios on the efficient frontier.
	function computeEnclosingCornerPortfolios_(targetVolatility, sigma, cornerPortfolios) {
		// The numerical accuracy for testing equality
		var eps = 1e-8;
		
		// The efficient frontier portfolios are provided from highest volatility
		// to lowest volatility, so that *_min below refers to properties of the portfolio
		// with the lowest volatility.
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];

		var volatility_min = computeVolatility_(sigma, weights_min);
		var volatility_max = computeVolatility_(sigma, weights_max);

		// If the target volatility is not reachable within numerical accuracy, 
		// return immediately.
		if (targetVolatility - volatility_max > eps || -eps > targetVolatility - volatility_min) {
			return [];
		}
		
		// If the target volatility is numerically reached on one of the
		// two extremal corner portfolios, return immediately.
		if (Math.abs(targetVolatility - volatility_min) <= eps) {
			return [[idx_min, weights_min, volatility_min]];
		}
		else if (Math.abs(targetVolatility - volatility_max) <= eps) {
			return [[idx_max, weights_max, volatility_max]];
		}
		
		// Otherwise, determine the two adjacent corner portfolios strictly enclosing the portfolio
		// with a volatility numerically equals to the target volatility, using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle point
			var idx_middle = Math.floor(idx_max + (idx_min - idx_max)/2); // formula avoiding numerical overflow
			var weights_middle = cornerPortfolios[idx_middle][0];
			var volatility_middle = computeVolatility_(sigma, weights_middle);
			
			// Determine in which sub-interval ]idx_max, idx_middle[ or ]idx_middle, idx_min[
			// lies the portfolio with the target volatility.
			if (volatility_middle - targetVolatility > eps) {
				idx_max = idx_middle;
				volatility_max = volatility_middle;
				weights_max = weights_middle;
			}
			else if (volatility_middle - targetVolatility < -eps) {
				idx_min = idx_middle;
				volatility_min = volatility_middle;
				weights_min = weights_middle;
			}
			else { // the target volatility is exactly attained on the idx_middle-th corner portfolio
				return [[idx_middle, weights_middle, volatility_middle]];
			}
		}

		
		// Return the computed adjacent corner portfolios, as well as
		// the associated function values.
		return [[idx_min, weights_min, volatility_min], [idx_max, weights_max, volatility_max]];
	}
	
	
	// ------	

	
	// Compute the (at most) two corner portfolios strictly enclosing the efficient 
	// portfolio with a volatility equals to the target volatility.
	var enclosingCornerPortfolios = computeEnclosingCornerPortfolios_(targetVolatility, sigma, cornerPortfolios);

	
	// Then:
	// - In case the desired target volatility is not reachable, return an empty portfolio 
	//
	// - In case there is a unique computed corner portfolio with a volatility
	// equals to the target volatility, return the associated portfolio weights
	//
	// - In case there are two corner portfolios strictly enclosing the efficient portfolio with 
	// a volatility equals to the target volatility, the weights associated to this efficient portfolio are
	// a strict convex combination of the weights of the two computed enclosing corner portfolios
	// (c.f. the reference): w = t*w_min + (1-t)*w_max, t in ]0,1[, with t now to be determined.
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
		var volatility_min = enclosingCornerPortfolios[0][2];
		var variance_min = volatility_min * volatility_min;
		
		var idx_max = enclosingCornerPortfolios[1][0];
		var weights_max = enclosingCornerPortfolios[1][1];
		var volatility_max = enclosingCornerPortfolios[1][2];
		var variance_max = volatility_max * volatility_max;
		
		// The procedure to compute t above is the following:
		// Let the volatility be V(w) = <Sigma*w/w>.
		// Then, by symmetry and bilinerarity of V, we have V(w) = t^2*V(w_min) + (1-t)^2*V(w_max) + 2*t*(1-t)*<Sigma*w_min/w_max>
		// and V(w) = targetVolatility^2
		// <=> t is the solution belonging to ]0,1[ of the second order polynomial equation
		// t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) -2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) - targetVolatility^2 = 0
		
		// Define the coefficients of the second order polynomial at^2 + bt + c
		var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
		var a = variance_min + variance_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
		var b = -2 * (variance_max - variance_cross); // 
		var c = variance_max - targetVolatility*targetVolatility; //always > 0
		
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
		
		var t;
		if (r1 > 0 && r1 < 1) {
			t = r1;
		}
		else if (r2 > 0 && r2 < 1) {
			t = r2;
		}
		else {
			throw new Error('internal error, the covariance matrix might not be semi-definite positive');
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
* @function computeMinimumVarianceEfficientPortfolio_
*
* @summary Compute the weights of the efficient mean-variance portfolio with the lowest volatility.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only mean-variance efficient portfolio of n assets with the lowest attainable volatility,
* i.e. the global minimum variance portfolio/leftmost portfolio on the efficient frontier.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Matrix_} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, a n by n matrix of real numbers.
* @param {Array<Array.<Object>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
* @return {Matrix_} the weights of the efficient portfolio with the minimum variance, a n by 1 Matrix_ of n real numbers
*/
function computeMinimumVarianceEfficientPortfolio_(cornerPortfolios) {
	if (cornerPortfolios.length === 0) { // this case should never occur
		throw new Error('internal error: empty list of corner portfolios');
	}
					
	// Return the computed portfolio weights
	return cornerPortfolios[cornerPortfolios.length - 1][0];
}


/**
* @function computeMaximumTargetVolatilityEfficientPortfolio_
*
* @summary Compute the weights of an efficient mean-variance portfolio subject to a maximum volatility
* constraint.
*
* @description This function returns the weights w_1,...,w_n associated to the maximally invested 
* and long-only mean-variance efficient portfolio of n assets subject to a maximum volatility constraint.
*
* The computed efficient portfolio is either:
* - Fully invested in case its volatility is equal to or lower than the maximum desired volatility
* - Maximally partially invested so that its volatility is equal to the maximum desired volatility, 
* in case a full investment would result in its volatility being strictly greater than the desired maximum volatility
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
*
* @param {Matrix_} mu the returns of the n assets in the considered universe, a n by 1 matrix of real numbers.
* @param {Matrix_} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, a n by n matrix of real numbers.
* @param {number} maxVolatility the desired maximum volatility of the portfolio, a strictly positive real number.
* @param {Array<Array.<Object>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
* - The corner portfolio risk aversion parameter, a positive real number
* @param {object} opt optional and/or mandatory parameters for the algorithm.
* @param {number} opt.maxIter the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
* @param {number} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {number} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the minimum weights for the assets to include in the portfolios with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Matrix_} the weights of the efficient portfolio with the desired maximum volatility, a n by 1 Matrix_ of n real numbers
*/
function computeMaximumTargetVolatilityEfficientPortfolio_(mu, sigma, maxVolatility, cornerPortfolios, opt) {
	// Internal function to compute the volatility of a portfolio
	function computeVolatility_(sigma, weights) {
		return Math.sqrt(Matrix_.vectorDotProduct(Matrix_.xy(sigma, weights), weights));
	}
	
	// ------	
	
	// Initializations
	var eps = 1e-8; // the numerical accuracy for testing equality	
	
	var nbAssets = sigma.nbColumns;
	
	var weights; // the weights to be computed

	
	// 1 - If the desired maximum volatility is greater than the highest attainable
	// volatility on the efficient frontier, the efficient portfolio is computed
	// as the rightmost portfolio on the efficient frontier.
	var highestVolatilityWeights = cornerPortfolios[0][0];
	var highestVolatility = computeVolatility_(sigma, highestVolatilityWeights);
	
	if (maxVolatility > highestVolatility + eps) {
		weights = highestVolatilityWeights;
	}
	else {
		// 2 - Otherwise, if the desired maximum volatility is lower than the lowest attainable
		// volatility on the efficient frontier, a new efficient frontier is computed
		// including a risk free asset with (0,0,0) return/volatility/covariance, and
		// the efficient portfolio is computed as the efficient portfolio with a target
		// volatility strictly equal to desired maximum volatility using the new efficient 
		// frontier.
		var lowestVolatilityWeights = cornerPortfolios[cornerPortfolios.length - 1][0];
		var lowestVolatility = computeVolatility_(sigma, lowestVolatilityWeights);
		
		if (maxVolatility < lowestVolatility - eps) {
			// Extend the returns/covariances of the input assets with a risk free asset,
			// of index (nbAssets + 1).
			var newMu = Matrix_.fill(nbAssets + 1, 1, 
									function(i, j) { 
										if (i <= nbAssets) {
											return mu.getValue(i, 1);
										}
										else {
											return 0;
										}
									});
			var newSigma = Matrix_.fill(nbAssets + 1, nbAssets + 1, 
										function(i, j) { 
											if (i <= nbAssets && j <= nbAssets) {
												return sigma.getValue(i, j);
											}
											else {
												return 0;
											}
										});

			// Compute the new corner portfolios defining the new efficient frontier
			// for the input assets with a risk free asset.
			var newCornerPortfolios = computeCornerPortfolios_(newMu, newSigma, opt);
			
			// Compute the efficient portfolio with a target volatility strictly equal to
			// desired maximum volatility.
			var newWeights;
			
			// Limit case: in case there is only one portfolio located on the new efficient frontier,
			// (for instance if all assets returns are negative), this portfolio must be the
			// (0, 0) portfolio, with a 0 variance.
			//
			// In this case, the efficient portfolio is unique and is entirely made of the risk free asset,
			// so that, in terms of original assets, this portfolio is not invested.
			//
			// (Otherwise, it means another portfolio with 0 variance possesses a strictly positive return,
			// which is not possible since the lowest attainable volatility on the input
			// efficient frontier was greater than the desired maximum volatility, and so
			// greater than 0).
			if (newCornerPortfolios.length === 1) { 
				newWeights = newCornerPortfolios[0][0];
				
				if (newCornerPortfolios[0][0].getValue(nbAssets + 1, 1) != 1) { // this case should never occur, as per description above
					throw new Error('internal error');
				}
			}
			else {				
				var efficientPortfolio = computeTargetVolatilityEfficientPortfolio_(newSigma, maxVolatility, newCornerPortfolios);
				
				if (efficientPortfolio.length === 0) { // this case should never occur, since the new efficient frontier includes a 0 volatility portfolio and at least 2 corner portfolios
					throw new Error('internal error');
				}
				else {
					newWeights = efficientPortfolio[0];
				}
			}
			
			// Transform the fully invested efficient portfolio with a risk free asset
			// into the associated partially invested efficient portfolio without a
			// risk free asset.
			weights = Matrix_.fill(nbAssets, 1, 
									function(i, j) { 
										if (i <= nbAssets) {
											return newWeights.getValue(i, 1);
										}
									});
		}
		
		// 3 - Otherwise, the efficient portfolio is computed as the efficient portfolio 
		// with a target volatility strictly equal to desired maximum volatility using
		// the input efficient frontier.
		else {
			var efficientPortfolio = computeTargetVolatilityEfficientPortfolio_(sigma, maxVolatility, cornerPortfolios);
			if (efficientPortfolio.length === 0) { // this case should never occur, since maxVolatility belongs to [lowestVolatility, highestVolatility]
				throw new Error('internal error');
			}
			else {
				weights = efficientPortfolio[0];
			}
		}
	}
	
	// Return the computed weights
	return weights;
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
* - Minimum weight of each asset to include in the portfolios
* - Maximum weight of each asset to include in the portfolios
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
			
		// Public functions to set the status of variales
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
		for (var i = 1; i <= nbAssets; ++i) {
			var lb_i = lowerBounds.getValue(i, 1);
			var ub_i = upperBounds.getValue(i, 1);
			
			if (lb_i > ub_i) {
				throw new Error('infeasible problem detected');
			}
		}

		var sum_lb = lowerBounds.sum();
		var sum_ub = upperBounds.sum();
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
		// are LOW, plus the linear program is is not degenerate
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
				var Mi_out_idx_out_idx = Mi.getValue(idx_out, idx_out);
				if (Math.abs(Mi_out_idx_out_idx) <= eps) {
					Mi_out_idx_out_idx = eps;
					//throw new Error('division by a numerical zero detected, the covariance matrix might not be positive semi-definite');
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
				var xi_j = M.getValue(idx_in, idx_in);
				for (var i = 1; i <= variablesInIdx.length; ++i) {
					var in_idx_i = variablesInIdx[i-1];
					
					xi_j -= M.getValue(idx_in, in_idx_i) * xi.getValue(in_idx_i, 1);
				}
				if (Math.abs(xi_j) <= eps) {
					xi_j = eps;
					//throw new Error('division by a numerical zero detected, the covariance matrix might not be positive semi-definite');
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
		lambda_e = max_([lambda_out, lambda_in, 0])[0];
		
		
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
		
		if (!Matrix_.areEqual(cornerPortfolio_i, cornerPortfolio, eps)) {
			++idx;
		}
		
		// Either replace a current portfolio or add a new one
		finalCornerPortfoliosWeights[idx] = [cornerPortfolio_i, lambda_i];
		
		// Prepare for the next iteration
		cornerPortfolio = cornerPortfolio_i;
		lambda = lambda_i;
	}
	
	// Resize the output efficient frontier array as required
	finalCornerPortfoliosWeights.length = idx + 1;
	
	// Return the final efficient frontier array
	return finalCornerPortfoliosWeights;
}
