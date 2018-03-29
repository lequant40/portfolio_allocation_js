/**
 * @file Functions related to equal risk budget portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskBudgetWeights
*
* @summary Compute the weights of the equal risk budget portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk budgets, defined as w_i = 1/sigma_i / (sum(1/sigma_j), j=1..n), i=1..n, with sigma_i the standard deviation 
* of the asset i.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the assets mean returns are proportional to their volatilities and all pair-wise correlations are equal.
* 
* @see <a href="https://doi.org/10.3905/jpm.2012.38.3.056 ">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* 
* @param {Matrix_|<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, an n by 1 matrix (i.e., vector) or an array of n real numbers statisfying sigma[i-1] = sigma_i.
* @param {object} opt optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to the equal risk budget portfolio, array of n real numbers.
*
* @example
* equalRiskBudgetWeights([0.1, 0.2]);
* // ~[0.59, 0.41]
*/
self.equalRiskBudgetWeights = function (sigma, opt) {
	// TODO: Checks, if enabled
	// Check that the values of sigma are strictly positive
	
	// ------
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Matrix_(sigma).elemMap(function(i,j,val) { return 1/Math.sqrt(val); })
	weights = weights.normalize(weights);
	
	// Return the computed weights
	return weights.toArray();
}

