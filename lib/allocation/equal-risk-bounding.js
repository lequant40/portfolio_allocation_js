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
* of n assets with equal risk bounding, defined as the weights with the property that the contribution 
* of each asset to the risk of the portfolio is equal.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* This portfolio is a trade-off between the MV portfolio and the EW portfolio.
* It can be viewed as an ERB portfolio tilted toward the assets less correlated with other assets.
*
*When short selling is not allowed, there is a unique RP portfolio
and it contains all assets in the market. In this case, the ERB approach might lead to the RP
portfolio or it might lead to portfolios with smaller variance that do not contain all assets,
and where the risk contributions of each asset included in the portfolio is strictly smaller
than in the RP portfolio

We show here that the Risk Parity approach
is theoretically dominated by an alternative similar approach that does not actually require
equally weighted risk contribution of all assets but only an equal upper bound on all such
risks.
This result implies that for every solution to the ERB model (6) there exists a subset P of
the set N = {1, . . . , n} of indices such that only the assets with indices in P have nonzero
weights, and such weights are those of the RP portfolio of the assets with indices in P. More
formally, we have the following:
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* Generate all subsets of the portfolios 2^n, so that for small n in the order of 15, will work in a reasonable time, but not for higher ns. (several mnutes for 20 assets)
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-016-0477-6">Cesarone, F. & Tardella F., Equal Risk Bounding is better than Risk Parity for portfolio selection, J Glob Optim (2017) 68: 439</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the underlying ERC algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the underlying ERC algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the equal risk bounding portfolio, array of n real numbers.
*
* @example
* equalRiskBoundingWeights([[0.1,0], [0,0.2]]);
* // XXX
*/
self.equalRiskBoundingWeights = function (sigma, opt) {	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Determine the size of the universe
	var nbAssets = sigma.nbRows;
	
	// Initialize the current minimum value of the risk contribution and the current list of associated assets/assets weights
	var minRCValue = Number.MAX_VALUE; // 
	var minRCAssetsIndexes = [];
	var minRCAssetsWeights = [];

	// Proceed to an exhaustive enumeration of all the subsets of the set {1,...,nbAssets},
	// in order to find the x^ERB as detailled in section 3.3.1, formula 22, of the reference.
	//
	// Skip the empty set, for obvious reasons.
	var nextSubsetIterator = new subsetsIterator_(nbAssets);
	var nextSubset = nextSubsetIterator.next();
	do {
		// Generate a new subset	
		var nextSubset = nextSubsetIterator.next();
		
		// Extract the associated assets indexes
		var assetsIndexes = nextSubset[1];
		
		// Extract the covariance matrix of the associated assets
		var subsetSigma = sigma.getSubmatrix(assetsIndexes, assetsIndexes);
		
		// Compute ERC weights for these assets
		var assetsWeights = self.equalRiskContributionWeights(subsetSigma, opt);
		
		// Compute the risk contribution of the first asset
		//
		// Note: all risk contributions are equal, by definition of the ERC portfolio
		assetsWeightsMatrix = new Matrix_(assetsWeights);
		var rcValue =  assetsWeightsMatrix.getValueAt(1,1) * Matrix_.product(subsetSigma, assetsWeightsMatrix).getValueAt(1,1) // = x_1 * (Sigma*x)_1
		
		// If the risk contribution of the current subset is lower than the current minimum risk contribution, it
		// becomes the new minimum risk contribution and the current subset becomes the new list of associated assets.
		if (rcValue < minRCValue) {
			minRCValue = rcValue;
			minRCAssetsIndexes = assetsIndexes;
			minRCAssetsWeights = assetsWeights;
		}
		// Otherwise, nothing needs to be done
	}
	while (nextSubset[0]);
	
	// At the end, map selected assets weight to full weights (all zero by default, those selected are the weights computed
	// Compute original assets weights, following formula 22 of the reference
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < minRCAssetsIndexes.length; ++i) {
		weights.setValueAt(minRCAssetsIndexes[i], 1, minRCAssetsWeights[i]);
	}
	
	// Return the computed weights (already normalized)
	return weights.toArray();
}

