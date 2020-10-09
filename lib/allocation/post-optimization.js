/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function postOptimizationWeights
*
* @summary Compute investable (integer, rational) portfolio weights from numerical (floating point) portfolio weights. 
*
* @description Given n (floating point) weights w_1,...,w_n associated to a long-only portfolio of n assets at most fully invested, 
this function computes n integer weights wi_1,...,wi_n associated to a long-only portfolio of total value opt.roundingMethodParams.portfolioValue, 
* composed of n assets with prices opt.roundingMethodParams.assetsPrices and which must be bought by multiples of opt.roundingMethodParams.sizeLots quantities, 
* satisfying wi_j * opt.roundingMethodParams.assetsPrices[j-1] * opt.roundingMethodParams.sizeLots[j-1], j=1..n is the quantity of asset j which must be possessed
* so that the proportion of the total value of asset j in the portfolio is close to the numerical weight w_j in the sense defined in the second reference.
*
*   Note: Due to the integer-valued nature of the underlying optimization problem, computing an optimal round lot solution is an NP-hard problem,
*   c.f. the second and third references, so that the heuristic optimization algorithm described in the fourth reference, called "Threshold Accepting", is used.
*         While the Threshold Accepting algorithm is not guaranteed to provide an optimal solution, it is reasonably guaranteed to provide a "good enough" solution.
*         One caveat though, because the Threshold Accepting algorithm is stochastic, different executions of this algorithm might return different weights.
*
* @see <a href="https://www.msci.com/documents/10199/b166dc05-b842-48fe-812d-3e310756c21c">Rong Xu and Scott Liu, Managing Odd Lot Trades with the Barra Optimizer, Research Insight, September 2013</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2261131">Steiner, Andreas, Accuracy and Rounding in Portfolio Construction (April 30, 2013)</a>
* @see <a href="https://link.springer.com/article/10.1007/s10732-010-9138-y">Gilli, M., Schumann, E. Optimal enough?. J Heuristics 17, 373–387 (2011).</a>
* @see <a href="https://pubsonline.informs.org/doi/abs/10.1287/ijoc.7.1.109">Kurt M. Bretthauer, Bala Shetty, Siddhartha Syam, (1995) A Branch and Bound Algorithm for Integer Quadratic Knapsack Problems. ORSA Journal on Computing 7(1):109-116.</a>
* 
* @param {Array.<number>} originalWeights the weights w_1,...,w_n associated to a long-only portfolio of n assets, array of n real numbers; in case the weights sum to < 1, the
* portfolio cash position is automatically computed; in case the weights sum to > 1, the closest weights (in the euclidean sense) summing to 1 are computed instead to have
* a feasible problem.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.portfolioValue the total value of the portfolio whose numerical weights w_1,...,w_n are provided
* @param {Array.<number>} opt.assetsPrices the prices of the n assets whose numerical weights w_1,...,w_n are provided
* @param {Array.<number>} opt.sizeLots the size of the lots to buy for the n assets whose numerical weights w_1,...,w_n are provided; defaults to an array of ones (i.e., minimal quantity to buy for all the assets is 1)
* @param {number} opt.cashWeight the optional cash weight of the portfolio whose numerical assets weights w_1,...,w_n are provided
*
* @return {Array.<Array.<number>, Array.<number>, number>}:
* - the integer weights wi_1,...,wi_n, array of n integers
* - the rational weights wrr_1,...,wrr_n corresponding to the integer weights wi_1,...,wi_n, array of n real numbers, computed with 
* wrr_j = wi_j * opt.roundingMethodParams.assetsPrices[j-1] * opt.roundingMethodParams.sizeLots[j-1] / opt.roundingMethodParams.portfolioValue, j=1..n
* - the quantity of remaining cash in the portfolio, a positive real number
*
*/
self.postOptimizationWeights = function (originalWeights, opt) {
	// Initialize the options structure
	if (opt === undefined) {
		opt = { roundingOptimizationMethodParams: {} };
	}
	if (opt.roundingOptimizationMethodParams === undefined) {
		opt.roundingOptimizationMethodParams = {};
	}

	// In case of roundlotting rounding method:
	// - Initialize the current portfolio value, as well as the current
	// assets prices and associated minimum round lots.
	//
	// - Initialize the default parameters for the Threshold Accepting
	// heuristic.
	var portfolioValue = opt.portfolioValue;
	var assetsPrices = opt.assetsPrices;
	var sizeLots = opt.sizeLots;
	if (portfolioValue == undefined || portfolioValue == 0) {
		throw new Error('missing portfolio value');
	}
	if (assetsPrices == undefined) {
		throw new Error('missing assets prices');
	}
	if (sizeLots == undefined) {
		sizeLots = Matrix_.ones(originalWeights.length, 1).data;
	}
	
	var nRounds = opt.roundingOptimizationMethodParams.nRounds;
	if (nRounds === undefined) {
		nRounds = 3;
	}
	var nSteps = opt.roundingOptimizationMethodParams.nSteps;	
	if (nSteps === undefined) {
		nSteps = 5000;
	}
	var nDeltas = opt.roundingOptimizationMethodParams.nDeltas;
	if (nDeltas === undefined) {
		nDeltas = nSteps;
	}
	
	
	// ------
	
	// Initialize the dimension
	var nbAssets = originalWeights.length;

	// Compute the closest weights of the provided weights, belonging to the unit full simplex (summing to <= 1)
	var alteredAssetsWeights = fullSimplexEuclidianProjection_(originalWeights);
		
	// Convert the altered weights to Matrix format, taking into account the cash
	var targetWeights = Matrix_.fill(nbAssets + 1, 1, function(i,j) { if (i <= nbAssets) { return alteredAssetsWeights[i-1]; } else { return 0;} });
	targetWeights.data[nbAssets] = 1 - targetWeights.sum();
	
	// The numerical zero
	var eps = 1e-12;

	// Initialize the price of the lots, taking into account the cash (lot size = 1, price = 1)
	var priceLots = Matrix_.fill(nbAssets + 1, 1, function(i,j) { if (i <= nbAssets) { return sizeLots[i-1] * assetsPrices[i-1]; } else { return 1;} });
	
	// Define the function representing the weights RMSE error to minimize,
	// taking in input q = [q_1,...q_nbAssets, q_cash], the quantity of lots associated
	// to each asset, with q_cash being the quantity of cash of the portfolio.
	function f(q, opt) {
		//
		if (opt) {
			var q_curr = opt.x_c;
			var q_curr_updated_idx = opt.x_c_updatedIndexes;
			var f_curr = opt.f_x_c;
		}
		
		// Compute the RMSE
		//
		// Note that the value below is not exactly the RMSE, because for speed-up
		// purposes, the full RMSE is not recomputed at every function call, and dropping 
		// the factor 1/2, as well as the square root, is then easier for incremental calculations.
		//
		// For speed-up purposes, if the indexes of the updated coordinates are provided,
		// the new RMSE is computed as follows: RMSE_new = RMSE_old + sum_i (i in updated coordinates) ( (w_orig[i] - q_new[i])^2 - (w_orig[i] - q_old[i])^2 )
		var rmse = 0;
		var invPortfolioValue = 1/portfolioValue;
		if (q_curr_updated_idx) {
			rmse = f_curr;
			for (var j = 0; j < q_curr_updated_idx.length; ++j) {
				var i = q_curr_updated_idx[j];
				
				// Cash (i == nbAssets) and no cash indexes are handled the same
				var w_orig_i = targetWeights.data[i];
				var w_i = q[i] * priceLots.data[i] * invPortfolioValue;
				var w_curr_i = q_curr[i] * priceLots.data[i] * invPortfolioValue;

				rmse += Math.pow(w_orig_i - w_i, 2) - Math.pow(w_orig_i - w_curr_i, 2);
			}
		}		
		else {
			// Compute the "RMSE" with the non-incremental formula as sum_i (w_orig[i] - q_new[i])^2
			for (var i = 0; i < nbAssets + 1; ++i) { // assets positions and cash position (i == nbAssets)
				w_orig_i = targetWeights.data[i];
				w_i = q[i] * priceLots.data[i] * invPortfolioValue;
				
				rmse += Math.pow(w_orig_i - w_i, 2);
			}
		}
		
		//
		return {f_x: rmse, f_x_context: {}};
	}

	// Define the function which generates a feasible solution "close" to another
	// feasible solution, used by the threshold accepting algorithm.
	//
	// The input q represents q = [q_1,...q_nbAssets, q_cash], the quantity of lots associated
	// to each asset, with q_cash being the quantity of cash of the portfolio.
	function neighbourGenerator(q, neighbourGeneratorParameters) {					
		// Internal function to compute a random integer
		// C.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
		function randomIndex(min, max) { // 0, n --> 0, n-1
			return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
		}
		
		// Randomly shuffle the sellable assets (i.e., those with potentially non-null weights), including cash
		var sellableAssets = neighbourGeneratorParameters.workingAssetsRandomIterator.next();
		
		// Determine an asset to sell among the randomly shuffled working assets.
		// 
		// As soon as the cash index is reached, the loop ends.
		var toSellIndex;
		for (var j = 0; j < sellableAssets.length; ++j) {
			//
			toSellIndex = sellableAssets[j];
			
			//
			if (q[toSellIndex] != 0 || toSellIndex == nbAssets) {
				break;
			}
		}
		
		// 
		var updatedAssetsIdx = [];
		
		// If the selected asset to sell is not cash, sell a random number of lots of a random asset
		if (toSellIndex != nbAssets) {
			// Determine a random number of lots of the selected asset to sell (minimum 1, maximum full quantity)
			var toSellNbLots = randomIndex(1, q[toSellIndex] + 1);
			
			// Sell them
			q[toSellIndex] -= toSellNbLots;
			var cashAddition = priceLots.data[toSellIndex] * toSellNbLots;
			q[nbAssets] += cashAddition;
			
			// 
			updatedAssetsIdx.push(toSellIndex);
			updatedAssetsIdx.push(nbAssets);
		}
		
		// Randomly shuffle the assets which can be bought (i.e., these with non null original weights), including cash
		var buyableAssets = neighbourGeneratorParameters.workingAssetsRandomIterator.next();

		// Determine an asset to buy among the randomly shuffled buyable assets.
		// 
		// As soon as the cash index is reached, the loop ends.
		var toBuyIndex;
		var toBuyMaxNbLots;
		for (var j = 0; j < buyableAssets.length; ++j) {
			//
			toBuyIndex = buyableAssets[j];
			
			//
			if (toBuyIndex == nbAssets) {
				break;
			}
			
			// Determine the maximum number of lots of the selected asset to buy (0 if not enough cash)
			toBuyMaxNbLots = Math.floor( q[nbAssets] / priceLots.data[toBuyIndex] );

			// If the maximum number of lots of the selected asset to buy is not 0, ends the loop
			if (toBuyMaxNbLots != 0) {
				break;
			}
		}
				
		// If the selected asset to buy is not cash, buy a random number of lots of the selected asset
		if (toBuyIndex != nbAssets) {
			// Determine a random number of lots of the selected asset to buy (between 1 and toBuyMaxNbLots)
			var toBuyNbLots = randomIndex(1, toBuyMaxNbLots + 1);

			// Buy them
			q[toBuyIndex] += toBuyNbLots;
			var cashRemoval = priceLots.data[toBuyIndex] * toBuyNbLots;
			q[nbAssets] -= cashRemoval;
			
			// 
			updatedAssetsIdx.push(toBuyIndex);
			updatedAssetsIdx.push(nbAssets);
		}	
		
		// Return the updated current portfolio lots
		var uniqueUpdatedAssetsIdx = updatedAssetsIdx.filter(function (x, i, a) { return a.indexOf(x) == i; })
		return {xx: q, x_updated_indexes: uniqueUpdatedAssetsIdx};
	}
	
	// Determine once for all the assets which can be bought (i.e., these with non null original weights), which
	// are also the same assets that can potentially be bought, including cash.
	var nbWorkingAssets = 1; // cash
	for (var i = 0; i < nbAssets; ++i) {
		if (targetWeights.data[i] > eps) {
			++nbWorkingAssets;
		}
	}
	var workingAssets = typeof UInt32Array === 'function' ? new UInt32Array(nbWorkingAssets) : new Array(nbWorkingAssets);
	for (var i = 0, j = 0; i < nbAssets; ++i) {
		if (targetWeights.data[i] > eps) {
			workingAssets[j] = i;
			++j;
		}
	}
	workingAssets[nbWorkingAssets - 1] = nbAssets; // cash index, at the end
	
	// Generate an initial feasible point for the threshold accepting algorithm,
	// c.f. "Implementing a Target Allocation" section of the third reference.
	//
	// The procedure belows guarantees that the cash is always in excess.
	var q0Vector = Matrix_.zeros(nbAssets + 1, 1);
	var q0 = q0Vector.data;
	q0[nbAssets] = portfolioValue; // cash position
	for (var j = 0; j < nbWorkingAssets - 1; ++j) { // excludes cash position
		var i = workingAssets[j];

		q0[i] = Math.floor( targetWeights.data[i] * portfolioValue / priceLots.data[i] );

		q0[nbAssets] -= q0[i] * priceLots.data[i];
	}
	if (q0[nbAssets] < 0) {
		throw new Error('internal error: negative cash during intitial feasible point generation');
	}
	
	// Now, in addition of having an initial feasible solution, iterate over all 
	// the assets i for which buying an additional lot is both feasible (i.e., enough cash)
	// and helps decreasing the value of the objective function.
	//
	// If several such assets are found, buy the asset which helps decreases the
	// objection function the most.
	//
	// Otherwise, stops the process.
	//
	// This procedure is similar in spirit to the heuristic described in the section 
	// "Feasible Integer Solution Heuristic" of the fifth reference, except that
	// the chosen iteration order is not dependent of the assets indexes.
	while (true) {
		var assetIdx = -1;
		var fMin = f(q0).f_x;
		for (var j = 0; j < nbWorkingAssets - 1; ++j) { // excludes cash position
			var i = workingAssets[j];
			
			// Buying an additional lot of asset i is feasible
			if (priceLots.data[i] <= q0[nbAssets]) {
				// Update the initial feasible solution with a potential update
				// to compute the change in the objective function value.
				q0[i] += 1;
				q0[nbAssets] -= priceLots.data[i];
				
				var fNew = f(q0).f_x;
				if (fNew < fMin) {
					assetIdx = i;
					fMin = fNew;
				}
				
				// Undo the update above
				q0[i] -= 1;
				q0[nbAssets] += priceLots.data[i];
			}
		}
		 
		// In case improving the initial feasible solution is not possible anymore,
		// stops the process.
		if (assetIdx == -1) {
			break;
		}
		
		// Buy one additional lot of the determined asset
		q0[assetIdx] += 1;
		q0[nbAssets] -= priceLots.data[assetIdx];
	}
	
	// The best initial feasible solution has been found
	if (q0[nbAssets] < 0) {
		throw new Error('internal error: negative cash during intitial feasible point generation');
	}
			
	//
	var workingAssetsRandomIterator = new randomPermutationsIterator_(null, workingAssets.slice(), true);
			
	// Solve the NP-hard optimization problem with an heuristic optimization method,
	// proven in the fourth reference to be effective for index tracking problems.
	var q = thresholdAcceptingSolve_(f, q0,
									 {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
									  neighbourGenerator: neighbourGenerator,
									  neighbourGeneratorParameters: {workingAssetsRandomIterator: workingAssetsRandomIterator}})[0];

	// Format the computed weights
	var roundedWeights = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1] * priceLots.data[i-1] / portfolioValue; }).toArray();
	var qLots = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1]; }).toArray();
	var cash = q[nbAssets];
	
	// Return the formatted computed weights
	return [qLots, roundedWeights, cash];
}