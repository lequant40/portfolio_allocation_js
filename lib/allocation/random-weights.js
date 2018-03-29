/**
 * @file Functions related to random weights portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function randomWeights
*
* @summary Compute the weights of a randomly generated portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only 
* random portfolio of n assets.
*
* Optionally, the following constraints can be added:
* - Minimum number of assets to include in the portfolio
* - Maximum number of assets to include in the portfolio
*
* This portfolio is not unique.
*
* Random portfolios have several applications in asset allocation, c.f. the first reference, as
* well as in trading strategies evaluation, c.f. the second reference.
*
* To be noted that the algorithms used internally allow the random portfolios to be generated uniformly
* among all the feasible portfolios.
*
* @see <a href="https://arxiv.org/abs/1008.3718">William T. Shaw, Monte Carlo Portfolio Optimization for General Investor Risk-Return Objectives and Arbitrary Return Distributions: a Solution for Long-only Portfolios</a>
* @see <a href="https://doi.org/10.1007/3-540-36626-1_11">Burns P. (2007) Random Portfolios for Performance Measurement. In: Kontoghiorghes E.J., Gatu C. (eds) Optimisation, Econometric and Financial Analysis. Advances in Computational Management Science, vol 9. Springer, Berlin, Heidelberg</a>
*
* @param {number} nbAssets the number of assets in the universe, natural integer superior or equal to 1.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.constraints.minAssets the minimum number of assets to include in the portfolio, an integer i satisfying 1 <= i <= nbAssets; defaults to 1.
* @param {number} opt.constraints.maxAssets the maximum number of assets to include in the portfolio, an integer j satisfying i <= j <= nbAssets; defaults to nbAssets.
* @return {Array.<number>} the weights corresponding to a random portfolio, array of real numbers of length nbAssets.
*
* @example
* randomWeights(5);
* // ~[0, 0 0.33, 0.33, 0.33]
*/
self.randomWeights = function (nbAssets, opt) {
	// TODO: Checks, if enabled

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	var nbMinAssets = opt.constraints.minAssets || 1;
	var nbMaxAssets = opt.constraints.maxAssets || nbAssets;
	
	// 1 - Generate the number of assets to include in the portfolio (uniform generation)
	var nbSelectedAssets = Math.floor(Math.random() * (nbMaxAssets - nbMinAssets +1)) + nbMinAssets;
	
	// 2 - Generate the indices of the assets to include in the portfolio (uniform generation)
	var selectedAssetsIdx = new randomKSubsetIterator_(nbAssets, nbSelectedAssets).next();
	
	// 3 - Generate the weights of the assets to include in the portfolio (uniform generation)
	//
	// Extra caution needs to be taken in case one of the weights is zero,
	// as exactly nbSelectedAssets must be included in the portfolio.
	var selectedAssetsWeights;
	var simplexSampler = new simplexRandomSampler_(nbSelectedAssets);
	while (true) {
		// Generate a sample of assets weights
		selectedAssetsWeights = simplexSampler.sample();
		
		// Reject the sample if there is one asset with a zero weight
		var rejectSample = false;
		for (var i = 0; i < nbSelectedAssets; ++i) {
			if (selectedAssetsWeights[i] == 0) {
				rejectSample = true;
				break;
			}
		}
		if (!rejectSample) {
			break;
		}
	}
	
	// Compute the final weights vector:
	// - The weights associated to assets not included in the portfolio at step 2 are set to zero
	// - The weights associated to assets included in the portfolio at step 2 are set to their values generated at step 3
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < nbSelectedAssets; ++i) {
		// Extract included assets information
		var assetIdx = selectedAssetsIdx[i];
		var assetWeight = selectedAssetsWeights[i];
		
		// Update the weights vector
		weights.setValueAt(assetIdx, 1, assetWeight);
	}

	// Return the computed weights
	return weights.toArray();
}

