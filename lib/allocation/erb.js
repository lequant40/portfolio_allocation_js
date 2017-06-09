/**
 * @file Functions related to equal risk budget portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskBudgetWeights
*
* @summary Compute the weights of the equally weighted risk budget portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk budgets, defined as w_i = 1/sigma_i / (sum(1/sigma_j), j=1..n), i=1..n, with sigma_i the standard deviation 
* of the asset i.
*
* These weights ensure that the risk budget of each asset, defined as the product of the asset’s weight and its standard deviation, is equal.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio for each stock is the same and all pair-wise correlations are equal.
* 
* @see <a href="https://ssrn.com/abstract=1949003">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* 
* @param {<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, array of n real numbers statisfying sigma[i-1] = sigma_i.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to the equal risk budget portfolio, array of n real numbers.
*
* @example
* equalRiskBudgetWeights([0.1, 0.2]);
* // ~[0.59, 0.41]
*/
self.equalRiskBudgetWeights = function (sigma, opt) {
	// TODO: Checks, if enabled
	// Check that the values of sigma are strictly positive
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Vector_(sigma).elementWisePower(-1/2).normalize();
	
	// Return the computed weights
	return weights.toArray();
}

