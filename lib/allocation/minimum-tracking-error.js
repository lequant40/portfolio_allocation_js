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
* This approach is expected to produce an exact solution within a reasonable amount of time
* for small n (e.g. n <= 15), but due to the combinatorial nature of the problem, 
* the computation for greater values of n will not be tractable, unless the maximum 
* number of assets is very small or very big compared to n.
*
* For these intractable cases, it is possible to use the heuristic optimization algorithm 
* described in the third reference, called "Threshold Accepting", which is not guaranteed
* to provide an optimal solution, but which is reasonably guaranteed to provide a "good enough" 
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
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-04.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.constraints.minNbAssets the minimum number of assets to include in the portfolio, an integer i satisfying 1 <= i <= nbAssets; defaults to 1 if opt.constraints.maxNbAssets is set.
* @param {number} opt.constraints.maxNbAssets the maximum number of assets to include in the portfolio, an integer j satisfying i <= j <= nbAssets; defaults to nbAssets if opt.constraints.minNbAssets is set.
* @param {string} opt.optimizationMethod the optimization method to use in order to compute the portfolio in case opt.constraints.minNbAssets or opt.constraints.maxNbAssets is set, a string either equals to:
* - 'combinatorial': usage of an exhaustive combinatorial search to exactly compute the optimal portfolio
* - 'thresholdAccepting': usage of the "Threshold Accepting" heuristic described in the third reference to approximately compute a "good enough" portfolio
; defaults to 'combinatorial'.
* @param {number} opt.optimizationMethodParams.nSteps the optional number of steps per threshold to use in case opt.optimizationMethod is equal to 'thresholdAccepting', defaults to 5000
* @param {number} opt.optimizationMethodParams.nDeltas the optional number of steps per threshold to use in case opt.optimizationMethod is equal to 'thresholdAccepting', defaults to opt.optimizationMethodParams.nSteps
* @param {number} opt.optimizationMethodParams.nRounds the optional number of random steps to generate the thresholds in case opt.optimizationMethod is equal to 'thresholdAccepting', defaults to 3
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
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
	if (opt.eps === undefined) {
		opt.eps = 1e-04;
	}
	if (opt.maxIter === undefined) {
		opt.maxIter = 10000;
	}
	if (opt.constraints.minNbAssets === undefined && opt.constraints.maxNbAssets) {
		opt.constraints.minNbAssets = 1;
	}
	if (opt.constraints.maxNbAssets === undefined && opt.constraints.minNbAssets) {
		opt.constraints.maxNbAssets = assetsReturns.length;
	}
	
	// Decode the options
	var eps = opt.eps;
	var maxIterations = opt.maxIter;
	
	var minNbAssets = opt.constraints.minNbAssets;
	var maxNbAssets = opt.constraints.maxNbAssets;
	var cardinalityConstraints = minNbAssets || maxNbAssets;
	
	var lowerBounds = opt.constraints.minWeights;
	var upperBounds = opt.constraints.maxWeights;
	
	
	// In case cardinality constraints are provided, the
	// default optimization method is a combinatorial search
	// in order to compute the optimal solution.
	var optimizationMethod;
	if (cardinalityConstraints) {
		optimizationMethod = opt.optimizationMethod;
		
		if (optimizationMethod === undefined) {
			optimizationMethod = 'combinatorial';
		}
		
		if (optimizationMethod != 'combinatorial' && optimizationMethod != 'thresholdAccepting') {
			throw new Error('unsupported optimisation method');
		}
	}

	
	// Initialize the default parameters for the Threshold Accepting
	// heuristic.
	var nRounds;
	var nSteps;
	var nDeltas;
	if (optimizationMethod == 'thresholdAccepting') {
		nRounds = opt.optimizationMethodParams.nRounds;
		if (nRounds === undefined) {
			nRounds = 3;
		}
		
		nSteps = opt.optimizationMethodParams.nSteps;
		if (nSteps === undefined) {
			nSteps = 5000;
		}

		nDeltas = opt.optimizationMethodParams.nDeltas;
		if (nDeltas === undefined) {
			nDeltas = nSteps;
		}
	}
	
	// ------
	
	
	// Initializations
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;
	
	// Convert the benchmark returns to matrix format
	var benchmarkReturns = new Matrix_(benchmarkReturns);
	
	
	// ----
	

	// Internal function to compute a portfolio minimizing the tracking error
	// v.s. a provided benchmark for a given number of assets nbAssets.
	//
	// The portfolio minimizing the tracking error volatility v.s. a provided benchmark
	// is a solution to the following smooth constrained convex optimization problem, 
	// c.f. formula 2 of the first reference, and formula 4 of the second reference
	// for the formulation of the soft inequality constraints:
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
	function computeMinimumTrackingErrorVolatilityPortfolio(assetsReturns, benchmarkReturns, 
															lowerBounds, upperBounds) {
		//
		var nbAssets = assetsReturns.nbColumns;
		
		// Define the matrix X
		var X = assetsReturns;

		// Define the vector Y
		var Y = benchmarkReturns;
		
		// Define the function representing f(w), the tracking error
		function f(w) {
			var te = Matrix_.xmy(Matrix_.xy(X, w), Y).vectorNorm('two');
			
			return 0.5 * te * te;
		}

		// Define the function representing the gradient of the function f(w),
		// which is equal to X^t * (X*w - Y), c.f.
		// proposition 2.1 of the second reference.
		function gradf(w) {
			var gte = Matrix_.txy(X, Matrix_.xmy(Matrix_.xy(X, w), Y));
			
			return gte;
		}

		// Define the characteristic function of the restricted unit simplex
		function g(w) {
			return simplexCharacteristicFunction_(w.data, lowerBounds, upperBounds);
		}

		// Define the proximal function associated to g, which is the orthogonal
		// projection on the restricted simplex.
		function proxg(w) {
			return new Matrix_(simplexEuclidianProjection_(w.data, lowerBounds, upperBounds));
		}

		// Define the initial point as the projection of the 1 vector 
		// on the restricted unit simplex.
		var x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(nbAssets, 1).data, lowerBounds, upperBounds));

		// Solve the convex optimization problem
		var sol = ccpsolveFISTA_(f, gradf, g, proxg, x0, {eps: eps, maxIter: maxIterations, maxLine: maxIterations});
		
		// Return the solution, whose first element is the computed portfolio weights
		return sol;
	} 
	
	
	// ----
	
	
	// Define the weights of the portfolio
	var weights;
	
	
	// In case no cardinality constraints are imposed, the portfolio minimizing the 
	// tracking error is the solution of a convex program.
	//
	// In case cardinality constraints are imposed:
	// - Either an exhaustive enumeration of all the subsets of the set {1,...,nbAssets} 
	// of size between minNbAssets and maxNbAssets is done, searching for the feasible 
	// portfolio minimizing the tracking error over all these subsets - "combinatorialSearch" algorithm
	//
	// - Or a proven heuristic optimization algorithm is used in order to find a "best" approximation
	// of the optimal portfolio, c.f. the third reference - "heuristicSearch" algorithm
	if (!cardinalityConstraints) {
		// Extract the assets returns
		var assetsReturns = Matrix_.fill(nbPeriods, nbAssets, function(i,j) { return assetsReturns[j-1][i-1]; });
		
		// Compute the solution of the convex program associated to the
		// portfolio minimizing the tracking error.
		weights = computeMinimumTrackingErrorVolatilityPortfolio(assetsReturns, benchmarkReturns, 
																 lowerBounds, upperBounds)[0];		
	}
	else {
		if (optimizationMethod == 'thresholdAccepting') {
			// Extract the assets returns
			var assetsReturns = Matrix_.fill(nbPeriods, nbAssets, function(i,j) { return assetsReturns[j-1][i-1]; });
					
			// Define default lower and upper bounds constraints in case they are not provided
			if (!lowerBounds) {
				lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
				for (var i = 0; i < nbAssets; ++i) {
					lowerBounds[i] = 0;
				}
			}
			if (!upperBounds) {
				upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
				for (var i = 0; i < nbAssets; ++i) {
					upperBounds[i] = 1;
				}
			}
			
			// Define the function representing the tracking error to minimize
			function f(w) {
				var w = new Matrix_(w);
				
				var te = Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns).vectorNorm('two');
				
				return 0.5 * te * te;
			}

		    
			// Define the function which generates a feasible solution "close" to another
			// feasible solution, used by the threshold accepting algorithm.
			//
			// Adapted from Remarks on 'A comparison of some heuristic optimization methods'
			// http://enricoschumann.net/R/remarks.htm.
			function neighbourGenerator(x, neighbourGeneratorParameters) {		
				// Internal function to compute a random integer
				// C.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
				function randomIndex(min, max) { // 0, n --> 0, n-1
					return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
				}

				// Decode the input parameters
				var l = neighbourGeneratorParameters.lowerBounds;
				var u = neighbourGeneratorParameters.upperBounds;
		
				// The maximum default value of the position size to buy/sell
				var alpha = 5/100;
				
				// The numerical zero
				var eps = 1e-12;
				
				// Initialize the dimension
				var nbAssets = x.length;
				
				// For performances optimization reason, do not copy the current portfolio weights,
				// but alter it instead !
				var xx = x;
				
				// Compute the cardinality of the current portfolio weights
				var nbNonNullAssets = 0;
				for (var i = 0; i < nbAssets; ++i) {
					if (xx[i] > eps) {
						++nbNonNullAssets;
					}
				}
				
				// Determine the possible list of assets to buy
				//
				// - In case the number of non null assets is equal to the maximum
				// number of allowed assets:
				// -- An asset already present in the portfolio can always have its position
				// size increased (provided its upper bound is not reached) because this will not change
				// the number of non null assets present in the portfolio
				// -- An asset NOT already present in the portfolio can be added if 
				// another asset already present in the portfolio can be removed with the same position size
				// as the added asset
				//
				// - Otherwise, there is no specific restriction on the list of assets to buy
				var toBuyToSellIndexes = new Array(0);
				if (nbNonNullAssets == maxNbAssets) {
					// Case #1: Assets already present in the portfolio
					for (var i = 0; i < nbAssets; ++i) {
						if (xx[i] > eps && xx[i] < u[i] - eps) { // increasing the position size of asset i is feasible
							toBuyToSellIndexes.push([i, -1]);
						}
					}
					
					// Case #2: Assets NOT already present in the portfolio
					for (var i = 0; i < nbAssets; ++i) {
						// Determine possible assets j already present in the portfolio 
						// that can be swapped with asset i not already present in the portfolio.
						var toSellIndexes = new Array(0);
						if (xx[i] <= eps && xx[i] < u[i] - eps) { // asset i is not already in the portfolio and adding it to the portfolio is feasible
							// Compute the maximum potential position size to buy,
							// which is constrained by both the asset upper bound and
							// the asset lower bound
							var alphab = Math.max(l[i], Math.min(u[i] - xx[i], alpha));
							
							//
							for (var j = 0; j < nbAssets; ++j) {
								if (xx[j] > eps && (xx[j] - alphab <= eps)) { // asset j is already in the portfolio and its position size can be fully swapped with asset i
									toSellIndexes.push([j, xx[j]]);
								}
							}
						}
						
						// In case such possible assets j exist, it is feasible
						// to replace any one of them with asset i.
						if (toSellIndexes.length != 0) {
							toBuyToSellIndexes.push([i, toSellIndexes]);
						}
					}
				}
				else if (nbNonNullAssets < maxNbAssets) {
					for (var i = 0; i < nbAssets; ++i) {
						if (xx[i] < u[i] - eps) { // increasing the position size of asset i, or adding asset i is feasible
							toBuyToSellIndexes.push([i, -1]);
						}
					}
				}
				else {
					throw new Error('internal error: number of assets strictly greater than the allowed maximum number of assets');
				}
				
				// Generate a random asset to buy
				var toBuyToSell = toBuyToSellIndexes[randomIndex(0, toBuyToSellIndexes.length)];
				
				// 
				var toBuyIndex = toBuyToSell[0];

				// Compute the maximum position size that is possible to buy:
				// - In all cases, this position size is constrained by the asset upper bound
				// - Additionally, if the asset to buy is NOT already in the portfolio,
				// this position size must at minimum be equal to the asset lower bound
				alpha = Math.min(u[toBuyIndex] - xx[toBuyIndex], alpha);
				if (xx[toBuyIndex] <= eps) {
					alpha = Math.max(l[toBuyIndex], alpha);
				}
				
				// In case the possible assets to sell with their position quantity to sell
				// have not been determined yet, proceed with their determination.
				var toSellIndexes;
				if (toBuyToSell[1] == -1) {
					// Determine the possible list of assets to sell
					//
					// - In case the number of non null assets is equal to the minimum
					// number of allowed assets:
					// -- An asset already present in the portfolio can have its position size
					// reduced provided the resulting reduced position size will not force to sell the asset, 
					// because this will not change the number of non null assets present in the portfolio
					// (-- An asset already present in the portfolio whose size would reduced below its lower bound
					// but above 0 would not feasible, since this asset would then need to be fully sold, leaving cash
					// available after buying the asset to buy)
					//
					// - Otherwise, an asset already present in the portfolio can have its position size
					// reduced provided the resulting reduced position size will not result in having a position size
					// below its lower bound but above 0
					//
					toSellIndexes = new Array(0);
					if (nbNonNullAssets == minNbAssets) {
						for (var i = 0; i < nbAssets; ++i) {
							if (xx[i] > eps && xx[i] > eps + l[i]) { // decreasing the position size of asset i is feasible because the resulting weight will be > l[i], and so, will not force a full sell of the position
								// Compute the maximum position size that is possible to sell
								var alphas = Math.min(alpha, xx[i] - l[i] - 2*eps); // -2*eps is in order to prevent null positions in case there is no lower bound for the asset i, otherwise, the cardinality constraint would be violated
								
								// If the asset to buy is NOT already in the portfolio,
								// its minimum position size must be equal to its lower bound,
								// so that the position size of the asset to sell must be
								// sufficient.
								//
								// Otherwise, there is no additional constraint on the position size
								// if the asset to sell.
								if (xx[toBuyIndex] <= eps) {
									if (alphas >= l[toBuyIndex]) {
										toSellIndexes.push([i, alphas]);
									}
								}
								else {
									toSellIndexes.push([i, alphas]);
								}
								
							}
						}
					}
					else if (nbNonNullAssets > minNbAssets) {
						for (var i = 0; i < nbAssets; ++i) {
							if (xx[i] > eps && (xx[i] - alpha <= eps)) { // full sell of the asset i is feasible
								toSellIndexes.push([i, xx[i]]);
							}
							else if (xx[i] > eps && (xx[i] > eps + l[i])) { // decreasing the position size of asset i is feasible
								toSellIndexes.push([i, Math.min(alpha, xx[i] - l[i])]);
							}
						}
					}
					else {					
						throw new Error('internal error: number of assets strictly lower than the allowed minimum number of assets');
					}
				}
				else {
					// The possible list of assets to sell has already been determined
					toSellIndexes = toBuyToSell[1];
				}
				
				// Generate a random asset to sell
				var toSellIndex = toSellIndexes[randomIndex(0, toSellIndexes.length)];
					
				// Compute the maximum position size to sell,
				// which has already been done at asset selection time above.
				alpha = toSellIndex[1];
				
				// Update the current portfolio weights
				xx[toSellIndex[0]] -= alpha;			
				xx[toBuyIndex] += alpha;
				
				// Return the updated current portfolio weights
				return xx;
			}
			
			// Generate an initial feasible point for the threshold accepting algorithm
			var x0 = self.randomWeights(nbAssets, { constraints: {
											            minNbAssets: minNbAssets, maxNbAssets: maxNbAssets, 
											            minWeights: lowerBounds, maxWeights: upperBounds 
											        }
												  });
												  
			// Solve the NP-hard optimization problem with an heuristic optimization method,
			// proven in the third reference to be effective for index tracking problems.
			weights = thresholdAcceptingSolve_(f, x0,
											   {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
											    neighbourGenerator: neighbourGenerator, 
												neighbourGeneratorParameters: { lowerBounds: lowerBounds, 
		                                                                        upperBounds: upperBounds }})[0];
			weights = new Matrix_(weights);
		}
		else if (optimizationMethod == 'combinatorial') {
			// Initialize the current minimum value of the tracking error
			// and the current list of associated assets/assets weights.
			var minTEValue = Infinity;
			var minTEAssetsIndexes = [];
			var minTEAssetsWeights = [];

			
			// Proceed to an exhaustive enumeration of all the subsets of the set {1,...,nbAssets}
			// satisfying the cardinality constraints in order to find 
			// the portfolio minimizing the tracking error over all these subsets.
			//
			// The feasibility of the enumerated subsets is guaranteed per their construction.
			for (var K = minNbAssets; K <= maxNbAssets; ++K) {
				var nextKSubsetIterator = new kSubsetsIterator_(nbAssets, K, false);
				var nextKSubset = nextKSubsetIterator.next();
				
				while (nextKSubset != -1) {
					// Extract the selected assets indexes
					var subsetNbAssets = nextKSubset.length;
					var subsetAssetsIdx = typeof UInt32Array === 'function' ? new UInt32Array(subsetNbAssets) : new Array(subsetNbAssets);
					for (var i = 0; i < nextKSubset.length; ++i) {
						subsetAssetsIdx[i] = nextKSubset[i];
					}

					// Extract the returns of the selected assets
					var subsetAssetsReturns = Matrix_.fill(nbPeriods, subsetNbAssets, 
														   function(i,j) { 
															   return assetsReturns[subsetAssetsIdx[j-1]-1][i-1]; 
														   });
					
					// Extract the lower and upper bounds constraints, if applicable
					var subsetLowerBounds;
					if (lowerBounds) {
						subsetLowerBounds = typeof Float64Array === 'function' ? new Float64Array(subsetNbAssets) : new Array(subsetNbAssets);
						for (var i = 0; i < subsetNbAssets; ++i) {
							subsetLowerBounds[i] = lowerBounds[subsetAssetsIdx[i]-1];
						}
					}
					var subsetUpperBounds;
					if (upperBounds) {
						subsetUpperBounds = typeof Float64Array === 'function' ? new Float64Array(subsetNbAssets) : new Array(subsetNbAssets);
						for (var i = 0; i < subsetNbAssets; ++i) {
							subsetUpperBounds[i] = upperBounds[subsetAssetsIdx[i]-1];
						}
					}
									
					// Compute the weights of the minimum tracking error portfolio for the selected assets
					//
					// Note: because the restricted simplex associated to the subset of selected assets
					// might be empty, special care must be taken.
					var subsetAssetsWeights;
					var subsetPortfolioTrackingError = Infinity;
					try {
						var subsetSol = computeMinimumTrackingErrorVolatilityPortfolio(subsetAssetsReturns, benchmarkReturns, 
																					   subsetLowerBounds, subsetUpperBounds);
																					   
						subsetAssetsWeights = subsetSol[0];
						subsetPortfolioTrackingError = subsetSol[1];
					}
					catch (e) {
						if (e.message !== "infeasible problem detected: the restricted simplex is empty") {
							throw(e);
						}
					}

					// If the tracking error of the current subset is lower than the current tracking error,
					// it becomes the new minimum tracking error and the current subset becomes the new 
					// list of selected assets.
					if (subsetPortfolioTrackingError < minTEValue) {
						minTEValue = subsetPortfolioTrackingError;
						minTEAssetsIndexes = subsetAssetsIdx;
						minTEAssetsWeights = subsetAssetsWeights;
					}
					// Otherwise, nothing needs to be done

					// Generate a new subset	
					var nextKSubset = nextKSubsetIterator.next();
				}
			}

			
			// Compute the original assets weights, in case a feasible minimum tracking error
			// portfolio has been found.
			if (minTEValue != Infinity) {
				weights = Matrix_.zeros(nbAssets, 1);
				for (var i = 0; i < minTEAssetsIndexes.length; ++i) {
					weights.data[minTEAssetsIndexes[i] - 1] = minTEAssetsWeights.data[i];
				}
			}
			else {
				throw new Error('infeasible problem detected');
			}
		}
		else {
			throw new Error('internal error: unsupported optimisation method');
		}
	}


	// Return the computed weights
	return weights.toArray();
}

