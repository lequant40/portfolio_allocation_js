/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function minimumTrackingErrorWeights
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested
* and long-only portfolio of n assets which minimizes the tracking error with regard to
* a provided benchmark.
*
* The definition of the tracking error taken is described in the first and second references, and is
* the sum of the squared deviations of returns between the portfolio and the benchmark,
* i.e. the tracking error volatility.
*
* Optionally, the following constraints can be added:
* - Minimum number of assets to include in the portfolio
* - Maximum number of assets to include in the portfolio
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* The algorithm used internally to solve the associated optimization problem is a FISTA-like 
* convex composite optimization algorithm.
*
* In case cardinality constraints are provided:
* - The associated optimization problem becomes strongly NP-hard, c.f. the second reference,
* so that an exhaustive computation of all the portfolios minimizing the tracking error 
* for each possible subset of assets is the only possible way to find an exact solution.
*
* Due to the combinatorial nature of the problem, it is not possible to exactly solve the problem
* for any real life value of n, so that the heuristic heuristic optimization algorithm 
* described in the third reference, called "Threshold Accepting", is used. This heuristic is not guaranteed
* to provide an optimal solution, but is reasonably guaranteed to provide a "good enough" 
* solution.
*
* One caveat though, because the Threshold Accepting algorithm is stochastic, different executions
* of this algorithm might return different weights.
*
* - The minimum/maximum weight of each asset is then to be understood as applying only to
* the assets selected by the optimization algorithm to be included in the portfolio.
* So, for example, in case a minimum weight constraint is defined for a non-selected asset,
* this minimum weight constraint is discarded.
*
* @see <a href="https://doi.org/10.1016/S0378-4266%2898%2900076-4">Markus Rudolf and Hans-jurgen Wolter and Heinz Zimmermann. A linear model for tracking error minimization. Journal of Banking and Finance. 1998</a>
* @see <a href="https://doi.org/10.1016/j.cor.2017.09.002">Purity Mutunge and Dag Haugland. Minimizing the tracking error of cardinality constrained portfolios. Computers & Operations Research Volume 90, February 2018, Pages 33-41</a>
* @see <a href="https://link.springer.com/chapter/10.1007/978-1-4757-5226-7_1">Gilli M., Këllezi E. (2002) The Threshold Accepting Heuristic for Index Tracking. In: Pardalos P.M., Tsitsiringos V.K. (eds) Financial Engineering, E-commerce and Supply Chain. Applied Optimization, vol 70. Springer, Boston, MA</a>
*
* @param {Array.<Array.<number>>} assetsReturns an array of n arrays of T real numbers, 
* with assetsReturns[i-1][j-1] the return of the i-th asset for the j-th period of time,
* i = 1..n, j = 1..T.
* @param {<Array.<number>} benchmarkReturns an array of T real numbers, 
* with benchmarkReturns[j-1] the return of the benchmark for the j-th period of time,
* j = 1..T..
* @param {object} opt optional parameters for the algorithm.
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @param {number} opt.constraints.minNbAssets the minimum number of assets to include in the portfolio, an integer i satisfying 1 <= i <= nbAssets; defaults to 1 if opt.constraints.maxNbAssets is set.
* @param {number} opt.constraints.maxNbAssets the maximum number of assets to include in the portfolio, an integer j satisfying i <= j <= nbAssets; defaults to nbAssets if opt.constraints.minNbAssets is set.
* @param {number} opt.optimizationMethodParams.maxIter the optional maximum number of iterations of the optimization algorithm in case of no cardinality constraints; defaults to 10000.
* @param {number} opt.optimizationMethodParams.maxIterationsInitPoint the optional maximum number of iterations of the algorithm to generate an initial feasible portfolio in case of cardinality constraints; defaults to 10000.
* @param {number} opt.optimizationMethodParams.eps the optional tolerance parameter for the convergence of the optimization algorithm in case of no cardinality constraints; defaults to 1e-04.
* @param {number} opt.optimizationMethodParams.nSteps the optional number of steps per threshold to use in case of cardinality constraints; defaults to 5000
* @param {number} opt.optimizationMethodParams.nDeltas the optional number of steps per threshold to use in case of cardinality constraints; defaults to opt.optimizationMethodParams.nSteps
* @param {number} opt.optimizationMethodParams.nRounds the optional number of random steps to generate the thresholds in case of cardinality constraints; defaults to 3
*
* @return {Array.<number>} the weights corresponding to the computed portfolio, array of n real numbers.
*
*/
self.minimumTrackingErrorWeights = function (assetsReturns, benchmarkReturns, opt) {
	// Initialize the options structure
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	
	// Initialize the options default values
	if (opt.optimizationMethodParams.eps === undefined) {
		opt.optimizationMethodParams.eps = 1e-04;
	}
	if (opt.optimizationMethodParams.maxIter === undefined) {
		opt.optimizationMethodParams.maxIter = 10000;
	}
	if (opt.constraints.minNbAssets === undefined && opt.constraints.maxNbAssets) {
		opt.constraints.minNbAssets = 1;
	}
	if (opt.constraints.maxNbAssets === undefined && opt.constraints.minNbAssets) {
		opt.constraints.maxNbAssets = assetsReturns.length;
	}

	// Decode the options
	var eps = opt.optimizationMethodParams.eps;
	var maxIterations = opt.optimizationMethodParams.maxIter;

	var fullInvestmentConstraint = true;
	if (opt.constraints.fullInvestment !== undefined) {
		fullInvestmentConstraint = opt.constraints.fullInvestment;
	}
	
	var minNbAssets = opt.constraints.minNbAssets;
	var maxNbAssets = opt.constraints.maxNbAssets;
	var cardinalityConstraints = minNbAssets || maxNbAssets;
	
	
	// Initialize the default parameters for the Threshold Accepting
	// heuristic in case of cardinality constraints.
	var nRounds = opt.optimizationMethodParams.nRounds;
	if (nRounds === undefined) {
		nRounds = 3;
	}
	var nSteps = opt.optimizationMethodParams.nSteps;
	if (nSteps === undefined) {
		nSteps = 5000;
	}
	var nDeltas = opt.optimizationMethodParams.nDeltas;
	if (nDeltas === undefined) {
		nDeltas = nSteps;
	}
	var maxIterationsInitPoint = opt.optimizationMethodParams.maxIterationsInitPoint;
	if (maxIterationsInitPoint === undefined) {
		maxIterationsInitPoint = 10000;
	}

	
	// ------
	
	
	// Initializations
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;

	// Convert the benchmark returns to matrix format
	var benchmarkReturns = new Matrix_(benchmarkReturns);

	// Extract the assets returns
	var assetsReturns = Matrix_.fill(nbPeriods, nbAssets, function(i,j) { return assetsReturns[j-1][i-1]; });
	
	// Define the lower/upper bounds
	var lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
	var upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		lowerBounds[i] = opt.constraints.minWeights ? opt.constraints.minWeights[i] : 0;
		upperBounds[i] = opt.constraints.maxWeights ? opt.constraints.maxWeights[i] : 1;
	}

	
	// ----
	
	
	// Define the weights of the portfolio
	var weights;
	
	
	// In case no cardinality constraints are imposed, the portfolio minimizing the 
	// tracking error is the solution of a convex program.
	//
	// In case cardinality constraints are imposed, a proven heuristic optimization algorithm is used in order to find
	// a "best" approximation of the optimal portfolio, c.f. the third referenc.
	if (!cardinalityConstraints) {
		// The portfolio minimizing the tracking error volatility v.s. a provided benchmark
		// is a solution to the following smooth constrained convex optimization problem, 
		// c.f. formula 2 of the first reference:
		//
		// argmax f(w) = 1/2 * ||X*w - Y||_2^2, with:
		// - X the nbPeriods by nbAssets matrix of returns of the assets
		// - Y the nbPeriods vector of benchmark returns
		// s.t. sum w_i = 1
		//      l <= w <= u
		//      (i.e., b belongs to a restricted unit simplex)
		//
		// This optimization problem is solved using a first-order method
		// for convex minimization.
		//
		// To be noted that in case the problem is not feasible, this method throws an exception.
		
		// Define the function representing f(w), the tracking error
		function f(w) {
			var te = Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns).vectorNorm('two');
			
			return 0.5 * te * te;
		}

		// Define the function representing the gradient of the function f(w),
		// which is equal to X^t * (X*w - Y), c.f. proposition 2.1 of the second reference.
		function gradf(w) {
			var gte = Matrix_.txy(assetsReturns, Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns));
			
			return gte;
		}

		// Define the characteristic function of the restricted unit simplex
		function g(w) {
			return simplexCharacteristicFunction_(w.data, lowerBounds, upperBounds);
		}

		// Define the proximal function associated to g, which is the orthogonal
		// projection on the restricted unit simplex.
		function proxg(w) {
			return new Matrix_(simplexEuclidianProjection_(w.data, lowerBounds, upperBounds));
		}

		// Define the initial point as the projection of the upper bounds vector 
		// on the restricted unit simplex.
		var x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(assetsReturns.nbColumns, 1).data, lowerBounds, upperBounds));

		// Solve the convex optimization problem
		var sol = ccpsolveFISTA_(f, gradf, g, proxg, x0, {eps: eps, maxIter: maxIterations, maxLine: maxIterations});
		
		// Return the solution, whose first element is the computed portfolio weights
		weights = sol[0];
	}
	else {
		// Define the function representing the tracking error to minimize
		function f(w, opt) {
			//
			var w = new Matrix_(w);
			
			// Compute the RMSE
			//
			// For speed-up purposes, if the indexes of the updated coordinates are provided,
			// the new RMSE is computed incrementally, using the fact that the matrix - vector product
			// assetsReturns * w_new = assetsReturns * w_old + v, with v depending only of the
			// updated w coordinates and the associated assetsReturns matrix coordinates.
			var assetsReturnsTw_m_b;
			if (opt) {
				//
				var w_curr = opt.x_c;
				var w_curr_updated_idx = opt.x_c_updatedIndexes;
				var f_curr = opt.f_x_c;
				
				var f_curr_context = opt.f_x_c_context;
				var assetsReturnsTx_c_m_b = f_curr_context.assetsReturnsTx_m_b;
				
				// Optimized formula to compute assetsReturns * w - benchmarkReturns 
				assetsReturnsTw_m_b = new Matrix_(assetsReturnsTx_c_m_b);
				for (var k = 0; k < w_curr_updated_idx.length; ++k) {
					var i = w_curr_updated_idx[k];
					
					var w_i_delta = w.data[i] - w_curr[i];
					for (var j = 0; j < assetsReturns.nbRows; ++j) {
						assetsReturnsTw_m_b.data[j] += assetsReturns.data[j*assetsReturns.nbColumns + i] * w_i_delta;
					}
				}
			}
			else {
				// Non-optimized formulas to compute assetsReturns * w - benchmarkReturns
				assetsReturnsTw_m_b = Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns);
			}
			
			// 
			var te = assetsReturnsTw_m_b.vectorNorm('two');
							
			//
			var rmse = 0.5 * te * te;
			
			//
			return {f_x: rmse, f_x_context: {assetsReturnsTx_m_b: assetsReturnsTw_m_b}};
		}

		
		// Define the function which generates a feasible solution "close" to another
		// feasible solution, used by the threshold accepting algorithm.
		//
		// Initially adapted from Remarks on 'A comparison of some heuristic optimization methods'
		// http://enricoschumann.net/R/remarks.htm, but diverged after for performances 
		// and numerical stability reasons.
		function neighbourGenerator(x, neighbourGeneratorParameters) {		
			// Internal function to compute a random integer
			// C.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
			function randomIndex(min, max) { // 0, n --> 0, n-1
				return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
			}
			
			// Internal function to compute a random number between two values
			// C.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
			function random(min, max) { // [min, max[
				return Math.random() * (max - min) + min;
			}

			
			// Decode the input parameters
			var l = neighbourGeneratorParameters.lowerBounds;
			var u = neighbourGeneratorParameters.upperBounds;
			
			
			// For performances optimization reason, do not copy the current portfolio weights,
			// but alter it instead !
			var xx = x;
		
			// Compute the cardinality of the current portfolio weights, and proceed with sanity checks
			var nbNonNullAssets = 0;
			var cashPosition = 1;
			for (var i = 0; i < nbAssets; ++i) {				
				if ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) {
					++nbNonNullAssets;
					cashPosition -= xx[i];
					
					if (xx[i] < l[i]) {
						throw new Error('internal error: asset present in the portfolio, but strictly below its lower bound');
					}
					
					if (xx[i] > u[i]) {
						throw new Error('internal error: asset present in the portfolio, but strictly above its upper bound');
					}
				}
			}
			cashPosition = Math.max(0, cashPosition);
			
			if (nbNonNullAssets > maxNbAssets) {
				throw new Error('internal error: number of assets strictly greater than the allowed maximum number of assets');
			}
			else if (nbNonNullAssets < minNbAssets) {
				throw new Error('internal error: number of assets strictly lower than the allowed minimum number of assets');
			}				

			//
			var assetsSwitches = [];
			
			// Generate a random order to visit the assets in order to compute feasible switches
			var assetsIndexes = assetsRandomIterator.next();
			
			for (var k = 0; k < nbAssets; ++k) {
				var i = assetsIndexes[k] - 1;

				// Determine if the asset can be sold, with its associated assets that can be bought.
				//
				// An asset which can be sold is:
				// - An asset that can be sold partially, with proceedings fully re-injected in the portfolio
				//   in order to increase the position size of an asset already existing (no impact on the number of assets) OR
				//   to buy a new asset (+1 asset)
				//
				// - An asset that can be fully sold, with proceeding fully re-injected in the portfolio
				//   in order to increase the position size of an asset already existing (-1 asset) OR
				//   to buy a new asset (no impact on the number of assets)
				//
				// To be noted that in case partial investment constraint is imposed, selling an asset either partially
				// (no impact on the number of assets) or fully (-1 asset) to increase the cash position is an additional possibility.
				if ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) { // The asset i is existing in the portfolio and can be partially of fully sold
					// Compute the minimum and maximum position sizes to re-allocate in case the asset i is partially  sold
					var alphaPartialSellMax = xx[i] - l[i];
					var alphaPartialSellMin = Math.min(0, alphaPartialSellMax);
					
					// Compute the position size to re-allocate in case the asset is fully sold
					var alphaFullSell = xx[i];
					
					// In case of partial investment constraint:
					// - Partially selling asset i for cash is always feasible
					// - Fully selling asset i for cash is feasible if the current number of non null assets is 
					//   strictly greater than the minimum number of assets
					if (!fullInvestmentConstraint) {
						// Push the switch in the list of feasible assets switches
						assetsSwitches.push([i, -1, alphaPartialSellMin, alphaPartialSellMax]);
						
						if (nbNonNullAssets > minNbAssets) {
							// Push the switch in the list of feasible assets switches
							assetsSwitches.push([i, -1, alphaFullSell]);
						}
					}						
					
					// 
					for (var j = 0; j < nbAssets; ++j) {
						// The asset j is NOT existing in the portfolio and it might be possible to buy it.
						if ((xx[j] < l[j] || xx[j] == 0) && xx[j] < u[j]) {
							// Compute the minimum and maximum positions size that would be possible to
							// to buy for asset j.
							var alphaBuyMin = l[j];
							var alphaBuyMax = Math.max(alphaBuyMin, u[j] - xx[j]);
							
							// If the current number of non null assets is strictly lower than the maximum number of assets,
							// it might be possible to buy asset j in case the asset i is partially sold.
							if (nbNonNullAssets < maxNbAssets) {
								// If the minimum position size to buy for asset j is greater than the 
								// minimum position size to re-allocate in case the asset i is partially sold AND 
								// lower than the maximum position size to re-allocate in case the asset i is partially sold,
								// the partial move asset i -> asset j is feasible.
								if (alphaPartialSellMin <= alphaBuyMin && alphaBuyMin <= alphaPartialSellMax) {
									// Compute the minimum and maximum position sizes for the partial move asset i -> asset j
									var alphaMin = alphaBuyMin;
									var alphaMax = Math.min(alphaBuyMax, alphaPartialSellMax);
									
									// Push the switch in the list of feasible assets switches
									assetsSwitches.push([i, j, alphaMin, alphaMax]);
								}
							}
							
							// Whatever the current number of assets, it might be possible to buy asset j in case
							// asset i is fully sold.

							// If the minimum position size to buy for asset j is lower than the 
							// position size to re-allocate in case the asset i is fully sold AND if
							// the maximum position size to buy for asset j is greater than the
							// position size to re-allocate in case the asset i is fully sold,
							// the full move asset i -> asset j is feasible.
							if (alphaBuyMin <= alphaFullSell && alphaFullSell <= alphaBuyMax) {
								// Push the switch in the list of feasible assets switches
								assetsSwitches.push([i, j, alphaFullSell]);
							}
						}
						
						// The asset j is ALREADY existing in the portfolio and it might be possible to buy more of it.
						if (((xx[j] > 0 && l[j] == 0) || (xx[j] >= l[j] && l[j] > 0)) && xx[j] < u[j]) {
							// Compute the minimum and maximum positions size that would be possible to
							// to buy for asset j.
							var alphaBuyMax = u[j] - xx[j];
							var alphaBuyMin = Math.min(0, alphaBuyMax);
							
							// If the current number of non null assets is strictly greater than the minimum number of assets,
							// it might be possible to buy more asset j in case the asset i is fully sold.
							if (nbNonNullAssets > minNbAssets) {
								// If the minimum position size to buy for asset j is lower than the 
								// position size to re-allocate in case the asset i is fully sold AND if
								// the maximum position size to buy for asset j is greater than the position 
								// size to re-allocate in case the asset i is fully sold,
								// the full move asset i -> asset j is feasible.
								if (alphaBuyMin <= alphaFullSell && alphaFullSell <= alphaBuyMax) {
									// Push the switch in the list of feasible assets switches
									assetsSwitches.push([i, j, alphaFullSell]);
								}
							}
							
							// Whatever the current number of assets, it might be possible to buy more asset j in case
							// asset i is partially sold.
							
							// If the minimum position size to buy for asset j is greater than the 
							// minimum position size to re-allocate in case the asset i is partially sold AND 
							// lower than the maximum position size to re-allocate in case the asset i is partially sold,
							// the partial move asset i -> asset j is feasible.
							if (alphaPartialSellMin <= alphaBuyMin && alphaBuyMin <= alphaPartialSellMax) {
								// Compute the minimum and maximum position sizes for the partial move asset i -> asset j
								var alphaMin = alphaBuyMin;
								var alphaMax = Math.min(alphaBuyMax, alphaPartialSellMax);
																
								// Push the switch in the list of feasible assets switches
								assetsSwitches.push([i, j, alphaMin, alphaMax]);
							}
						}
					}
				}
				
				// Break further below after the cash investment test, as soon as first assets feasible switches have been found, 
				// since the assets loop is randomized
				var assetsSwitchesFound = false;
				if (assetsSwitches.length != 0) {
					assetsSwitchesFound = true;
				}
				
				// In case of partial investment constraint, buying an asset already existing in the portfolio (no impact on the number of assets)
				// is always feasible as long as the cash position is sufficient to do it.
				//
				// Adding a new asset (+1 asset) is feasible as long as the cash position is sufficient to do it AND if the current number
				// of assets allows it.
				if (!fullInvestmentConstraint) {
					if ( ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) && xx[i] < u[i] ){ // The asset i is existing in the portfolio and is not at its upper bound
						// Compute min and max position that it would be possible to buy more for asset i 
						var alphaBuyMax = u[i] - xx[i];
						var alphaBuyMin = Math.min(0, alphaBuyMax);
						
						// If the cash position is sufficient, it is possible to buy the asset i
						if (cashPosition >= alphaBuyMin) {
							// Compute the minimum and maximum position sizes for the move cash -> asset i
							var alphaMin = alphaBuyMin;
							var alphaMax = Math.min(alphaBuyMax, cashPosition);
										
							// Push the switch in the list of feasible assets switches
							assetsSwitches.push([-1, i, alphaMin, alphaMax]);
						}
					}

					if ((xx[i] < l[i] || xx[i] == 0) && xx[i] < u[i]) { // The asset i is NOT existing in the portfolio and is not at its upper bound
						if (nbNonNullAssets < maxNbAssets) {
							// Compute the minimum and maximum positions size that would be possible to
							// to buy for asset i.
							var alphaBuyMin = l[i];
							var alphaBuyMax = Math.max(alphaBuyMin, u[i] - xx[i]);
							
							// If the cash position is sufficient, it is possible to buy the asset i
							if (cashPosition >= alphaBuyMin) {
								// Compute the minimum and maximum position sizes for the move cash -> asset i
								var alphaMin = alphaBuyMin;
								var alphaMax = Math.min(alphaBuyMax, cashPosition);

								// Push the switch in the list of feasible assets switches
								assetsSwitches.push([-1, i, alphaMin, alphaMax]);
							}
						}
					}
				}

				// Break when the first feasible switches have been found, since the assets loop is randomized
				if (assetsSwitchesFound) {
					break;
				}
			}
			

			// Compute a random index, corresponding to a random feasible switch asset i -> asset j 
			if (assetsSwitches.length == 0) {
				throw new Error('internal error: no feasible switch');
			}
			var assetSwitch = assetsSwitches[randomIndex(0, assetsSwitches.length)];
			
			// 
			var updatedAssetsIdx = [];
			
			// Extract the asset to sell, the asset to buy, and compute, if applicable, a random
			// quantity to sell/buy.
			//
			// Update the current portfolio weights
			var sellAssetIdx = assetSwitch[0];
			var buyAssetIdx = assetSwitch[1];
			var quantity;
			if (assetSwitch.length == 3) {
				quantity = assetSwitch[2];
				if (sellAssetIdx != -1) {
					xx[sellAssetIdx] = 0; // full sell case
					updatedAssetsIdx.push(sellAssetIdx);
				}
				if (buyAssetIdx != -1) {
					xx[buyAssetIdx] += quantity;
					updatedAssetsIdx.push(buyAssetIdx);
				}
			}
			else if (assetSwitch.length == 4) {
				quantity = random(assetSwitch[2], assetSwitch[3]);
				if (sellAssetIdx != -1) {
					xx[sellAssetIdx] -= quantity;
					updatedAssetsIdx.push(sellAssetIdx);
				}
				if (buyAssetIdx != -1) {
					xx[buyAssetIdx] += quantity;
					updatedAssetsIdx.push(buyAssetIdx);
				}
			}
			else {
				throw new Error('internal error: unexpected lenght of the switches structure');
			}
			
			
			// Return the updated current portfolio weights
			var uniqueUpdatedAssetsIdx = updatedAssetsIdx.filter(function (x, i, a) { return a.indexOf(x) == i; })
			return {xx: xx, x_updated_indexes: uniqueUpdatedAssetsIdx};
		}

		// Generate an initial feasible point for the threshold accepting algorithm.
		//
		// In case no such feasible point is generated, abort the computation.
		var x0;
		try {
			// In case of partial investment constraint, the initial feasible point
			// does not need to have its assets positions summing to one.
			var minExposure = 1;
			if (!fullInvestmentConstraint) {
				minExposure = Math.random();
			}
			
			x0 = self.randomWeights(nbAssets, { maxIter: maxIterationsInitPoint,
												constraints: {
													minExposure: minExposure, maxExposure: 1,
													minNbAssets: minNbAssets, maxNbAssets: maxNbAssets, 
													minWeights: lowerBounds, maxWeights: upperBounds
												}
											  });
		}
		catch (e) {
			if (e.message === "maximum number of iterations reached") {
				throw new Error("infeasible problem detected");
			}
			else {
				throw(e);
			}
		}
		
		//
		var assetsRandomIterator = new randomPermutationsIterator_(nbAssets, undefined, true);
		
		// Solve the NP-hard optimization problem with an heuristic optimization method,
		// proven in the third reference to be effective for index tracking problems.
		var sol = thresholdAcceptingSolve_(f, x0,
										   {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
											neighbourGenerator: neighbourGenerator, 
											neighbourGeneratorParameters: { assetsRandomIterator: assetsRandomIterator,
																			lowerBounds: lowerBounds, 
																			upperBounds: upperBounds }});

		//
		var weights = new Matrix_(sol[0]);
	}


	// Return the computed weights
	return weights.toArray();
}

