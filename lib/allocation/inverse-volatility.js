/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function inverseVolatilityWeights
*
* @summary Compute the weights of the inverse volatility portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets weighted by their inverse volatility, defined as w_i = 1/sigma_i / (sum(1/sigma_j), j=1..n), i=1..n, with sigma_i the standard deviation 
* of the asset i.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the assets mean returns are proportional to their volatilities and all pair-wise correlations are equal.
* 
* @see <a href="https://doi.org/10.3905/jpm.2012.38.3.056">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* @see <a href="https://doi.org/10.1007/s10479-017-2474-7">Ardia, D., Bolliger, G., Boudt, K. et al. The impact of covariance misspecification in risk-based portfolios. Ann Oper Res 254, 1–16 (2017).</a>
* 
* @param {Matrix_|<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, an n by 1 matrix (i.e., vector) or an array of n real numbers satisfying sigma[i-1] = sigma_i.
* @param {object} opt optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to computed portfolio, array of n real numbers.
*
* @example
* inverseVolatilityWeights([0.1, 0.2]);
* // ~[0.59, 0.41]
*/
self.inverseVolatilityWeights = function (sigma, opt) {
	// TODO: Checks, if enabled
	// Check that the values of sigma are strictly positive
	
	// ------
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Matrix_(sigma).elemMap(function(i,j,val) { return 1/Math.sqrt(val); })
	weights = weights.normalize(weights);
	
	// Return the computed weights
	return weights.toArray();
}

