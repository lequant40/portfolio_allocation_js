/**
 * @file Functions related to equal risk bounding portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskBoundingWeights
*
* @summary Compute the weights of the equal risk bounding portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk bounding, defined as the portfolio for which the contribution of each asset 
* included in the portfolio to the risk of the portfolio is equal, c.f. the reference.
*
* A property of the equal risk bounding portfolio with no bounds contraints is that it is either equal to the equal risk contribution portfolio
* (if all the assets are included in the portfolio), or to a portfolio with a smaller variance than the 
* equal risk contribution portfolio and for which the risk contribution of each asset included in the 
* portfolio is strictly smaller than in the equal risk contribution portfolio (if some assets are not included
* in the portfolio).
*
* Hence, the equal risk bounding portfolio is an equal risk contribution portfolio restricted
* to a (possibly strict) subset of the n assets.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* The algorithm used to solve the associated optimisation problem is an exhaustive computation of all the ERC portfolios
* defined on all the 2^n subsets of the n assets.
* This approach is expected to produce a solution within a reasonable amount of time for small n (e.g. n <= 20),
* but due to the combinatorial nature of the problem, the computation for greater values of n will not be tractable.
*
* To be noted that in case minimum/maximum weights are provided, they are to be understood as applying only to
* the assets selected by the optimization algorithm to be included in the portfolio, due to the combinatorial nature of
* the algorithm.
* So, for example, in case a minimum weight constraint is defined for a non-selected asset,
* this minimum weight constraint is discarded.
*
* @see <a href="https://doi.org/10.1007/s10898-016-0477-6">Cesarone, F. & Tardella F., Equal Risk Bounding is better than Risk Parity for portfolio selection, J Glob Optim (2017) 68: 439</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the underlying ERC algorithm, a strictly positive real number; defaults to 1e-10.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
* @return {Array.<number>} the weights corresponding to the equal risk bounding portfolio, array of n real numbers.
*
* @example
* equalRiskBoundingWeights([[1,-9/10, 3/5], [-9/10, 1,-1/5],[ 3/5, -1/5, 4]]);
* // [0.5, 0.5, 0]
*/
self.equalRiskBoundingWeights = function (sigma, opt) {	
	// Create the options structure, if not defined,
	// and request the portfolio volatility to be providede
	// in output of the ERC algorithm.
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Determine the size of the universe
	var nbAssets = sigma.nbRows;
	
	// Initialize the current minimum value of the risk contribution and the current list of associated assets/assets weights
	var minRCValue = Infinity;
	var minRCAssetsIndexes = [];
	var minRCAssetsWeights = [];

	// Proceed to an exhaustive enumeration of all the subsets of the set {1,...,nbAssets},
	// in order to find the x^ERB as detailled in section 3.3.1, formula 22, of the reference.
	//
	// The empty set is skipped.
	var nextSubsetIterator = new subsetsIterator_(nbAssets);
	var nextSubset = nextSubsetIterator.next(); // empty set
	var nextSubset = nextSubsetIterator.next(); // "true" first set
	while (nextSubset != -1) {	
		// Extract the selected assets indexes
		var subsetAssetsIdx = nextSubset;
		var sizeSubset = subsetAssetsIdx.length;
		
		// Extract the covariance matrix of the selected assets
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);

		// Extract the lower/upper bounds constraints of the selected assets, if applicable
		var subsetMinWeights;
		if (opt.constraints.minWeights) {
			subsetMinWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubset) : new Array(sizeSubset);
			for (var i = 0; i < sizeSubset; ++i) {
				subsetMinWeights[i] = opt.constraints.minWeights[subsetAssetsIdx[i]-1];
			}
		}
		
		var subsetMaxWeights;
		if (opt.constraints.maxWeights) {
			subsetMaxWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubset) : new Array(sizeSubset);
			for (var i = 0; i < sizeSubset; ++i) {
				subsetMaxWeights[i] = opt.constraints.maxWeights[subsetAssetsIdx[i]-1];
			}
		}

		// Compute ERC weights for the selected assets, taking into account possible
		// subset non feasibility.
		try {
			var sol = self.equalRiskContributionWeights(subsetSigma, {eps: opt.eps,
			                                                          maxCycles: -1, 
																	  outputPortfolioVolatility: true,
			                                                          constraints: {minWeights:subsetMinWeights, maxWeights:subsetMaxWeights}});
			var subsetAssetsWeights = sol[0];
			var subsetPortfolioVolatility = sol[1];

			// Compute lambda_erc, c.f. the formula following the formula 3 of the reference.
			var rcValue = subsetPortfolioVolatility * subsetPortfolioVolatility / sizeSubset;
			
			// If the risk contribution of the current subset is lower than the current minimum risk contribution, it
			// becomes the new minimum risk contribution and the current subset becomes the new list of selected assets.
			//
			// Otherwise, nothing needs to be done.
			if (rcValue < minRCValue) {
				minRCValue = rcValue;
				minRCAssetsIndexes = subsetAssetsIdx;
				minRCAssetsWeights = subsetAssetsWeights;
			}
		}
		catch (e) {
			if (e.message !== 'infeasible problem detected: the restricted simplex is empty' && 
			    e.message !== 'internal error: maximum number of iterations reached when searching for a bracketing interval, the problem might be infeasible') {
					throw (e);
				}
		}
	
		// Generate a new subset	
		var nextSubset = nextSubsetIterator.next();
	}
	
	// In case no feasible portfolio has been generated, throw an error
	if (minRCValue == Infinity) {
		throw new Error('no feasible portfolio generated');
	}
	
	// Compute the original assets weights, following formula 22 of the reference
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < minRCAssetsIndexes.length; ++i) {
		weights.setValueAt(minRCAssetsIndexes[i], 1, 
		                   minRCAssetsWeights[i]);
	}
	
	// Return the computed weights (already normalized)
	return weights.toArray();
}

