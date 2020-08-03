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
* @description Given n (floating point) weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, 
this function either:
* - If opt.roundingMethod equals to "rationalizing", computes n (rational) weights wr_1,...,wr_n associated to a fully invested 
* and long-only portfolio of n assets satisfying:
* -- opt.roundingMethodParams.k * wr_i is a natural integer, i=1..n
* -- wr_1,...,wr_n are the closest weights to w_1,...,w_n, in the sense defined in the first reference
*
*   Note: Typical values of opt.roundingMethodParams.k are 10 (rounding to 10%), 20 (rounding to 5%) and 100 (rounding to 1%).
*
*
* - If opt.roundingMethod equals to "roundlotting", computes n integer weights wi_1,...,wi_n associated to a fully invested 
* and long-only portfolio of total value opt.roundingMethodParams.portfolioValue, composed of n assets with prices 
* opt.roundingMethodParams.assetsPrices and which must be bought by multiples of opt.roundingMethodParams.sizeLots quantities, satisfying
* wi_j * opt.roundingMethodParams.assetsPrices[j-1] * opt.roundingMethodParams.sizeLots[j-1], j=1..n is the quantity of asset j which must be possessed
* so that the proportion of the total value of asset j in the portfolio is close to the numerical weight w_j in the sense defined in the third reference.
*
*   Note: Due to the integer-valued nature of the underlying optimization problem, computing an optimal round lot solution is an NP-hard problem,
*   c.f. the third and fifth references, so that the heuristic optimization algorithm described in the fourth reference, called "Threshold Accepting", is used.
*         While the Threshold Accepting algorithm is not guaranteed to provide an optimal solution, it is reasonably guaranteed to provide a "good enough" solution.
*         One caveat though, because the Threshold Accepting algorithm is stochastic, different executions of this algorithm might return different weights.
*
* @see <a href="https://doi.org/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258.</a>
* @see <a href="https://www.msci.com/documents/10199/b166dc05-b842-48fe-812d-3e310756c21c">Rong Xu and Scott Liu, Managing Odd Lot Trades with the Barra Optimizer, Research Insight, September 2013</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2261131">Steiner, Andreas, Accuracy and Rounding in Portfolio Construction (April 30, 2013)</a>
* @see <a href="https://link.springer.com/article/10.1007/s10732-010-9138-y">Gilli, M., Schumann, E. Optimal enough?. J Heuristics 17, 373–387 (2011).</a>
* @see <a href="https://pubsonline.informs.org/doi/abs/10.1287/ijoc.7.1.109">Kurt M. Bretthauer, Bala Shetty, Siddhartha Syam, (1995) A Branch and Bound Algorithm for Integer Quadratic Knapsack Problems. ORSA Journal on Computing 7(1):109-116.</a>
* 
* @param {Array.<number>} originalWeights the weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, array of n real numbers.
* @param {object} opt optional parameters for the algorithm.
* @param {string} opt.roundingMethod the rounding method to use in order to compute the investable portfolio weights, a string either equals to:
* - 'rationalizing': usage of the rational rounding method described in the first reference
* - 'roundlotting': usage of the round lot optimization method described in the third reference, using the "Threshold Accepting" heuristic described 
* in the fourth reference to approximately compute "good enough" investable portfolio weights; parameters opt.roundingMethod.portfolioValue and 
* opt.roundingMethod.assetsPrices are mandatory in this case
; defaults to 'rationalizing'.
* @param {number} opt.roundingMethodParams.k the value to which the rounded weights will be a multiple of the inverse in case opt.roundingMethod is equal to "rationalizing",
* natural integer greater than or equal to 1; defaults to 20.
* @param {number} opt.roundingMethodParams.portfolioValue the total value of the portfolio whose numerical weights w_1,...,w_n are provided in case opt.roundingMethod is equal to "roundlotting"
* @param {Array.<number>} opt.roundingMethodParams.assetsPrices the prices of the n assets whose numerical weights w_1,...,w_n are provided in case opt.roundingMethod is equal to "roundlotting"
* @param {Array.<number>} opt.roundingMethodParams.sizeLots the size of the lots to buy for the n assets whose numerical weights w_1,...,w_n are provided in case opt.roundingMethod is equal to 
* "roundlotting"; defaults to an array of ones (i.e., minimal quantity to buy for all the assets is 1)
*
* @return {Array.<number>} in case opt.roundingMethod is equal to "rationalizing", the rational rounded weights wr_1,...,wr_n, array of n real numbers.
* @return {Array.<Array.<number>, Array.<number>, number>} in case opt.roundingMethod is equal to "roundlotting":
* - the integer weights wi_1,...,wi_n, array of n integers
* - the rational weights wrr_1,...,wrr_n corresponding to the integer weights wi_1,...,wi_n, array of n real numbers, computed with 
* wrr_j = wi_j * opt.roundingMethodParams.assetsPrices[j-1] * opt.roundingMethodParams.sizeLots[j-1] / opt.roundingMethodParams.portfolioValue, j=1..n
* - the quantity of remaining cash in the portfolio, a positive real number
*
* @example
* postOptimizationWeights([0.5759, 0.0671, 0.3570], {roundingMethodParams: 10});
* // [0.6, 0.1, 0.3]
* postOptimizationWeights([0.5759, 0.0671, 0.3570], {roundingMethodParams: 20});
* // [0.6, 0.05, 0.35]
* postOptimizationWeights([0.5759, 0.0671, 0.3570], {roundingMethodParams: 100});
* // [0.57, 0.07, 0.36]
*/
self.postOptimizationWeights = function (originalWeights, opt) {
	// Initialize the options structure
	if (opt === undefined) {
		opt = { roundingMethodParams: {}, 
		        roundingOptimizationMethodParams: {} };
	}
	if (opt.roundingOptimizationMethodParams === undefined) {
		opt.roundingOptimizationMethodParams = {};
	}
	if (opt.roundingMethodParams === undefined) {
		opt.roundingMethodParams = {};
	}
	
	// Initialize the options default values
	if (opt.roundingMethod === undefined) {
		opt.roundingMethod = "rationalizing";
	}
	if (opt.roundingMethod == "rationalizing" && opt.roundingMethodParams.k == undefined) {
		opt.roundingMethodParams.k = 20;
	}
	
	// Decode the options
	var roundingMethod = opt.roundingMethod;
	if (roundingMethod != 'rationalizing' && roundingMethod != 'roundlotting') {
		throw new Error('unsupported rounding method');
	}
	
	// In case of rationalizing rounding method
	var k;
	if (roundingMethod == "rationalizing") {
		k = opt.roundingMethodParams.k;
	}

	// In case of roundlotting rounding method:
	// - Initialize the current portfolio value, as well as the current
	// assets prices and associated minimum rount lots, adding an additional
	// cash as last asset.
	//
	// - Initialize the default parameters for the Threshold Accepting
	// heuristic.
	var portfolioValue = opt.roundingMethodParams.portfolioValue;
	var assetsPrices = opt.roundingMethodParams.assetsPrices;
	var sizeLots = opt.roundingMethodParams.sizeLots;
	if (roundingMethod == 'roundlotting') {
		if (portfolioValue == undefined || portfolioValue == 0) {
			throw new Error('missing portfolio value');
		}
		if (assetsPrices == undefined) {
			throw new Error('missing assets prices');
		}
		if (sizeLots === undefined) {
			sizeLots = typeof Float64Array === 'function' ? new Float64Array(originalWeights.length) : new Array(originalWeights.length);		
			for (var i = 0; i < sizeLots.length; ++i) {
				sizeLots[i] = 1;
			}
		}
	}
	
	var nRounds;
	var nSteps;
	var nDeltas;
	if (roundingMethod == 'roundlotting') {
		nRounds = opt.roundingOptimizationMethodParams.nRounds;
		if (nRounds === undefined) {
			nRounds = 3;
		}
		
		nSteps = opt.roundingOptimizationMethodParams.nSteps;
		if (nSteps === undefined) {
			nSteps = 5000;
		}

		nDeltas = opt.roundingOptimizationMethodParams.nDeltas;
		if (nDeltas === undefined) {
			nDeltas = nSteps;
		}
	}
	
	
	// ------
	
	//
	if (roundingMethod == "rationalizing") {
		// Call to the simplex rational rounding method
		var roundedWeights = simplexRationalRounding_(originalWeights, k);
		
		// Return the computed weights
		return roundedWeights;
	}
	else if  (roundingMethod == "roundlotting") {
		// The numerical zero
		var eps = 1e-12;
				
		// Initialize the dimension
		var nbAssets = originalWeights.length;

		// Convert the original weights to matrix format,
		// adding a cash position as an additional last asset.
		var originalWeights = new Matrix_.fill(nbAssets + 1, 1, function(i,j) { if (i <= nbAssets) { return originalWeights[i-1]; } else { return 0;} });
		
		// Define the function representing the weights RMSE error to minimize,
		// taking in input q = [q_1,...q_nbAssets, q_cash], the quantity of lots associated
		// to each asset, with q_cash being the quantity of cash of the portfolio.
		function f(q) {
			//
			var w = new Matrix_.fill(nbAssets + 1, 1, 
			                         function(i,j) { 
									    if (i <= nbAssets) { // true asset
											return q[i-1] * sizeLots[i-1] * assetsPrices[i-1] / portfolioValue;
										}
                                        else { // cash
											return q[nbAssets] / portfolioValue;
									    }
									  });

			// Note that the value below is not exactly the RMSE, but the RMSE up
			// to a constant.
			var rmse = Matrix_.xmy(originalWeights, w).vectorNorm('two');
			
			//
			return rmse;
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
			
			// Determine the assets which can be sold (excluding cash)
			var nbSellableAssets = 0;
			for (var i = 0; i < nbAssets; ++i) {
				if (q[i] > 0) {
					++nbSellableAssets;
				}
			}
			var sellableAssets = typeof UInt32Array === 'function' ? new UInt32Array(nbSellableAssets) : new Array(nbSellableAssets);
			for (var i = 0, j = 0; i < nbAssets; ++i) {
				if (q[i] > 0) {
					sellableAssets[j] = i;
					++j;
				}
			}
			
			// If the portfolio is not 100% in cash, sell a random number of lots of a random asset
			if (sellableAssets.length != 0) {
				// Determine a random asset to sell
				var toSellIndex = sellableAssets[randomIndex(0, nbSellableAssets)];

				// Determine a random number of lots of the selected asset to sell (minimum 1, maximum full quantity)
				var toSellNbLots = randomIndex(1, q[toSellIndex] + 1);
				
				// Sell them
				q[toSellIndex] -= toSellNbLots;
				var cashAddition = assetsPrices[toSellIndex] * sizeLots[toSellIndex] * toSellNbLots;
				q[nbAssets] += cashAddition;
			}
			
			// Determine the assets which can be bought (i.e., these with non null original weights), including cash
			var nbBuyableAssets = 1; // cash
			for (var i = 0; i < nbAssets; ++i) {
				if (originalWeights.getValue(i + 1, 1) > eps) {
					++nbBuyableAssets;
				}
			}
			var buyableAssets = typeof UInt32Array === 'function' ? new UInt32Array(nbBuyableAssets) : new Array(nbBuyableAssets);
			for (var i = 0, j = 0; i < nbAssets; ++i) {
				if (originalWeights.getValue(i + 1, 1) > eps) {
					buyableAssets[j] = i;
					++j;
				}
			}
			buyableAssets[nbBuyableAssets - 1] = nbAssets;
			
			
			// Determine a random asset to buy
			var toBuyIndex = buyableAssets[randomIndex(0, nbBuyableAssets)];
			
			// If the selected asset to buy is not cash, buy a random number of lots of the selected asset
			if (toBuyIndex != nbAssets) {
				// Determine the maximum number of lots of the selected asset to buy
				var toBuyMaxNbLots = Math.floor( q[nbAssets] / (assetsPrices[toBuyIndex] * sizeLots[toBuyIndex]) );
				
				// Determine a random number of lots of the selected asset to buy (minimum 0 if not enough cash, maximum maximum quantity above)
				var toBuyNbLots = randomIndex(0, toBuyMaxNbLots + 1);

				// Buy them
				q[toBuyIndex] += toBuyNbLots;
				var cashRemoval = assetsPrices[toBuyIndex] * sizeLots[toBuyIndex] * toBuyNbLots;
				q[nbAssets] -= cashRemoval;
			}	

			// Return the updated current portfolio lots
			return q;
		}
		
		// Generate an initial feasible point for the threshold accepting algorithm,
		// c.f. "Implementing a Target Allocation" section of the third reference.
		//
		// The procedure belows guarantees that the cash is always in excess.
		var q0 = typeof Float64Array === 'function' ? new Float64Array(nbAssets + 1) : new Array(nbAssets + 1);
		var val = 0;
		for (var i = 0; i < nbAssets; ++i) {
			var lotPrice = assetsPrices[i] * sizeLots[i];
			q0[i] = Math.floor( originalWeights.getValue(i + 1, 1) * portfolioValue / lotPrice );

			val += q0[i] * lotPrice;
		}
		var cash = portfolioValue - val;
		if (cash < 0) {
			throw new Error('internal error: negative cash during intitial feasible point generation');
		}
		q0[nbAssets] = cash;

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
			var fMin = f(q0);
			for (var i = 0; i < nbAssets; ++i) {
				var lotPrice = assetsPrices[i] * sizeLots[i];
				
				// Buying an additional lot of asset i is feasible
				if (lotPrice < cash) {
					// Update the initial feasible solution with a potential update
					// to compute the change in the objective function value.
					q0[i] += 1;
					q0[nbAssets] = cash - lotPrice;
					
					var fNew = f(q0);
					if (fNew <= fMin) {
						assetIdx = i;
						fMin = fNew;
					}
					
					// Undo the update above
					q0[i] -= 1;
					q0[nbAssets] = cash;
				}
			}
			 
			// In case improving the initial feasible solution is not possible anymore,
			// stops the process.
			if (assetIdx == -1) {
				break;
			}
			
			// Buy one additional lot of the determined asset
			q0[assetIdx] += 1;
			cash -= assetsPrices[assetIdx] * sizeLots[assetIdx];
		}
		
		// The best initial feasible solution has been found
		if (cash < 0) {
			throw new Error('internal error: negative cash during intitial feasible point generation');
		}
		q0[nbAssets] = cash;

		// Solve the NP-hard optimization problem with an heuristic optimization method,
		// proven in the fourth reference to be effective for index tracking problems.
		var q = thresholdAcceptingSolve_(f, q0,
					  					 {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
										  neighbourGenerator: neighbourGenerator})[0];

		// Format the computed weights
		var roundedWeights = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1] * sizeLots[i-1] * assetsPrices[i-1] / portfolioValue; }).toArray();
		var qLots = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1]; }).toArray();
		var cash = q[nbAssets];
		
		// Return the formatted computed weights
		return [qLots, roundedWeights, cash];
	}
	else {
		throw new Error('internal error: unsupported rounding method');
	}
}