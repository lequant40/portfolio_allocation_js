/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



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
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
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

	// Internal function to compute the Sharpe ratio of a portfolio,
	// as well as the intermediate values of return and volatility.
	function computeSharpeRatio_(mu, sigma, rf, weights) {
		// The numerator: <mu/w> - rf
		var ret = portfolioReturn(weights, mu, sigma);
		var excessRet = ret - rf;
		
		// The denominator: Sqrt(<Sigma*w/w>)
		var vol = portfolioVolatility(weights, mu, sigma);
		
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

	
	// Compute the corner portfolios defining the efficient frontier
	var cornerPortfolios = computeCornerPortfolios_(mu, sigma, opt);
	
	
	// Restrict the efficient frontier to the domain of definition
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
	var volatility = portfolioVolatility(weights, mu, sigma);
	if (volatility < eps) {
		// The domain of definition of the Sharpe ratio is not the
		// whole efficient frontier, so that the efficient frontier
		// needs to be restricted.
		
		// Compute the first efficient portfolio with a strictly positive 
		// volatility.
		var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, 
		                                                               {constraint: "volatility",
															            constraintValue: eps});		
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

	
	// Further restrict the efficient frontier to the domain of strict positivity
	// of the Sharpe ratio.
	//
	// To be noted that the domain of strict positivity of the Sharpe ratio
	// can be empty in case there is no feasible portfolio on the efficient
	// frontier with a strictly positive excess return.
	var idx = cornerPortfolios.length - 1;
	var weights = cornerPortfolios[idx][0];
	var ret = portfolioReturn(weights, mu, sigma);
	if (ret < rf + eps) {
		// The domain of strict positivity of the Sharpe ratio is not the
		// whole efficient frontier, so that the efficient frontier
		// needs to be restricted.
		
		// Compute the first efficient portfolio with a strictly positive 
		// excess return.
		var efficientPortfolio = computeMeanVarianceEfficientPortfolio_(mu, sigma, cornerPortfolios, 
		                                                                {constraint: "return",
															             constraintValue: rf + eps});
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


	// On the restricted efficient frontier, the Sharpe ratio is a pseudo-concave 
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
