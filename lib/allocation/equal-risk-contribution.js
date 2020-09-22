/**
 * @file Functions related to equal risk contributions portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskContributionWeights
*
* @summary Compute the weights of the equal risk contribution portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk contributions.
*
* This portfolio has the property that the contribution of each asset to the risk of the portfolio is equal,
* and is a special case of the more generic risk budgeting portfolio, with all risk budgets
* equal, c.f. the first reference.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* To be noted that in case weights constraints are defined, the concept of "equal risk contribution"
* does not make any sense, c.f. the third reference, and the associated optimization problem might not have any solution.
*
* The algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is semi-definite positive.
*
* @see <a href="https://doi.org/10.3905/jpm.2010.36.4.060">Maillard, S., Roncalli, T., Teiletche, J.: The properties of equally weighted risk contribution portfolios. J. Portf. Manag. 36, 60–70 (2010)</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3331184">RONCALLI Thierry, RICHARD Jean-Charles, CONSTRAINED RISK BUDGETING PORTFOLIOS, </a>
*
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.epsSdp the tolerance parameter for testing the semi-definite positiveness of the covariance matrix, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.maxCycles the maximum number of cycles of the algorithm, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.nbCycles the exact number of cycles of the algorithm, a strictly positive natural integer, in which case the values of opt.eps and opt.maxCycles are discarded,
* or -1; defaults to -1.
* @param {boolean} opt.outputPortfolioVolatility a boolean indicating whether the portfolio volatility should be provided in output
* (if set to true) or not (if set to false); defaults to false.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
* @return {Array.<number>|Array.<Array.<number>>} if opt.outputPortfolioVolatility is set to false, the weights corresponding to the equal risk contribution portfolio, 
* array of n real numbers, and if opt.outputPortfolioVolatility is set to true, an array arr of two elements:
* - arr[0], the weights corresponding to the equal risk contribution portfolio, array of n real numbers
* - arr[1], the volatility of the computed equal risk contribution portfolio, a real number
*
* @example
* equalRiskContributionWeights([[0.1,0], [0,0.2]]);
* // ~[0.59, 0.41]
*/
self.equalRiskContributionWeights = function (sigma, opt) {	
	// The ERC portfolio is a specific case of the more general risk budgeting portfolio, with equal risk budgets.
	//
	// Generate equal risk budgets: rb_i = 1/nbAssets, i=1..nbAssets
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var rb = Matrix_.fill(nbAssets, 1, function (i,j) { return 1/nbAssets; });

	// Compute the associated risk budgeting weights
	return self.riskBudgetingWeights(sigma, rb, opt);
}

