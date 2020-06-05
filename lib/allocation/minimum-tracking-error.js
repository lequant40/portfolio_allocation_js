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
* The definition of the tracking error taken is described in the first and third references, and is
* the sum of the squared deviations of returns between the portfolio and the benchmark,
* i.e. the tracking error volatility.
*
* Optionally, the following constraints can be added:
* - Minimum number of assets to include in the portfolio
* - Maximum number of assets to include in the portfolio
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* Optionally, the following soft constraints can be added (these constraints will be taken into
* account as best as possible, but their violation is possible):
* - Any linear inequality constraints on the portfolio weights, in the form Ai * w <= bi,
* with Ai an IE by n matrix and bi an IE by 1 matrix
*
* The algorithm used internally to solve the associated optimization problem is a FISTA-like 
* convex composite optimization algorithm.
*
* In case cardinality constraints are provided:
* - The associated optimization problem becomes strongly NP-hard, c.f. the third reference,
* so that an exhaustive computation of all the portfolios minimizing the tracking error 
* for each possible subset of assets is the only possible way to find an exact solution.
*
* This approach is expected to produce an exact solution within a reasonable amount of time
* for small n (e.g. n <= 15), but due to the combinatorial nature of the problem, 
* the computation for greater values of n will not be tractable, unless the maximum 
* number of assets is very small or very big compared to n.
*
* - The minimum/maximum weight of each asset is then to be understood as applying only to
* the assets selected by the optimization algorithm to be included in the portfolio.
* So, for example, in case a minimum weight constraint is defined for a non-selected asset,
* this minimum weight constraint is discarded.
*
* @see <a href="https://doi.org/10.1016/S0378-4266%2898%2900076-4">Markus Rudolf and Hans-jurgen Wolter and Heinz Zimmermann. A linear model for tracking error minimization. Journal of Banking and Finance. 1998</a>
* @see <a href="https://epubs.siam.org/doi/10.1137/0917020">R. Bramley and B. Winnicka. Solving Linear Inequalities in a Least Squares Sense. SIAM J. Sci. Comput., 17(1), 275–286.</a>
* @see <a href="https://doi.org/10.1016/j.cor.2017.09.002">Purity Mutunge and Dag Haugland. Minimizing the tracking error of cardinality constrained portfolios. Computers & Operations Research Volume 90, February 2018, Pages 33-41</a>
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
* @param {Array.<number>} opt.constraints.minWeights an array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @param {Array.<Array.<number>>} opt.softConstraints.Ai an array of IE arrays of nbAssets real numbers, 
* with opt.softConstraints.Ai[i-1] the i-th lhs of a soft linear inequality constraint on the weights of the 
* assets to include in the portfolio, i = 1..IE.
* @param {<Array.<number>} opt.softConstraints.bi an array of IE real numbers, 
* with opt.softConstraints.bi[i-1] the i-th rhs of a soft linear inequality constraint on the weights of the
* assets to include in the portfolio, i = 1..IE.
* @param {number} opt.softConstraints.lambdai the penalty parameter associated to the soft linear inequality constraints, 
* a strictly positive real number; defaults to 1.
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
	if (opt.softConstraints  === undefined) {
		opt.softConstraints = {};
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
	if (opt.softConstraints.lambdai === undefined && (opt.softConstraints.Ai && opt.softConstraints.bi)) {
		opt.softConstraints.lambdai = 1;
	}
	
	// Decode the options
	var eps = opt.eps;
	var maxIterations = opt.maxIter;
	
	var minNbAssets = opt.constraints.minNbAssets;
	var maxNbAssets = opt.constraints.maxNbAssets;
	var cardinalityConstraints = (minNbAssets || maxNbAssets) ? true : false;
	
	var lowerBounds = opt.constraints.minWeights;
	var upperBounds = opt.constraints.maxWeights;
	
	var Ai = opt.softConstraints.Ai;
	var bi = opt.softConstraints.bi;
	var lambdai = opt.softConstraints.lambdai;
	var softInequalityConstraints = (Ai && bi) ? true : false;
	
	
	// ------
	
	
	// Initializations
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;
	var nbSoftInequalityConstraints = softInequalityConstraints ? bi.length : 0;
	
	// Convert the benchmark returns to matrix format
	var benchmarkReturns = new Matrix_(benchmarkReturns);
	
	// Convert the rhs of the soft inequality constraints to matrix format, if applicable
	if (softInequalityConstraints) {
		var bi = new Matrix_(bi);
	}		
	
	// ----
	

	// Internal function to compute a portfolio minimizing the tracking error
	// v.s. a provided benchmark for a given number of assets subsetNbAssets.
	//
	// The portfolio minimizing the tracking error volatility v.s. a provided benchmark
	// is a solution to the following smooth constrained convex optimization problem, 
	// c.f. formula 2 of the first reference, and formula 4 of the second reference
	// for the formulation of the soft inequality constraints:
	//
	// argmax f(w) = 1/2 * ||X*w - Y||_2^2 + lambdai/2 * ||(Ai*w - bi)_+||_2^2, with:
	// - X the nbPeriods by subsetNbAssets matrix of returns of the assets
	// - Y the nbPeriods vector of benchmark returns
	// - Ai the nbInequalityConstraints by subsetNbAssets matrix of the soft inequality constraints
	// - bi the nbInequalityConstraints vector of the soft inequality constraints
	// s.t. sum w_i = 1
	//      l <= w <= u
	//      (i.e., b belongs to a restricted unit simplex)
	//
	// This optimization problem is solved using a first-order method
	// for convex minimization.
	//
	// To be noted that in case the problem is not feasible, this method throws an exception.
	function computeMinimumTrackingErrorVolatilityPortfolio(subsetNbAssets, 
	                                                        subsetAssetsReturns, benchmarkReturns, 
															subsetLowerBounds, subsetUpperBounds, 
															subsetAi, subsetBi) {
		// Define the matrix X
		var X = subsetAssetsReturns;

		// Define the vector Y
		var Y = benchmarkReturns;
		
		// Define the function representing f(w)
		function f(w) {
			var te = Matrix_.xmy(Matrix_.xy(X, w), Y).vectorNorm('two');
			
			if (softInequalityConstraints) {
				var ineqe = Matrix_.xmy(Matrix_.xy(subsetAi, w), subsetBi).elemMap(function(i,j,val) { return Math.max(0, val); }).vectorNorm('two');
				
				return 0.5 * te * te + 0.5 * lambdai * ineqe * ineqe;
			}
			else {			
				return 0.5 * te * te;
			}
		}

		// Define the function representing the gradient of the function f(w),
		// which is equal to X^t * (X*w - Y) + lambdai * Ai^t * (Ai*w - bi)_+, c.f.
		// proposition 2.1 of the second reference.
		function gradf(w) {
			var gte = Matrix_.txy(X, Matrix_.xmy(Matrix_.xy(X, w), Y));
			
			if (softInequalityConstraints) {
				var gineqe = Matrix_.txy(subsetAi, Matrix_.xmy(Matrix_.xy(subsetAi, w), subsetBi).elemMap(function(i,j,val) { return Math.max(0, val); }));
				
				return Matrix_.axpby(1, gte, lambdai, gineqe);
			}
			else {
				return gte;
			}
		}

		// Define the characteristic function of the restricted unit simplex
		function g(w) {
			return simplexCharacteristicFunction_(w.data, subsetLowerBounds, subsetUpperBounds);
		}

		// Define the proximal function associated to g, which is the orthogonal
		// projection on the restricted simplex.
		function proxg(w) {
			return new Matrix_(simplexEuclidianProjection_(w.data, subsetLowerBounds, subsetUpperBounds));
		}

		// Define the initial point as the projection of the 1 vector 
		// on the restricted unit simplex.
		var x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(subsetNbAssets, 1).data, subsetLowerBounds, subsetUpperBounds));

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
	// In case cardinality constraints are imposed, an exhaustive enumeration of
	// all the subsets of the set {1,...,nbAssets} of size between minNbAssets and
	// maxNbAssets is done, searching for the feasible portfolio minimizing the tracking error
	// over all these subsets.
	if (!cardinalityConstraints) {
		// Extract the assets returns
		var assetsReturns = Matrix_.fill(nbPeriods, nbAssets, function(i,j) { return assetsReturns[j-1][i-1]; });
		
		// Extract the soft inequality constraints lhs, if applicable
		var Ai;
		if (softInequalityConstraints) {
			Ai = new Matrix_(Ai);
		}
		
		// Compute the solution of the convex program associated to the
		// portfolio minimizing the tracking error.
		weights = computeMinimumTrackingErrorVolatilityPortfolio(nbAssets, 
																 assetsReturns, benchmarkReturns, 
																 lowerBounds, upperBounds,
																 Ai, bi)[0];		
	}
	else {
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
				
				// Extract the soft inequality constraints lhs, if applicable
				var subsetAi;
				if (softInequalityConstraints) {
					subsetAi = Matrix_.fill(nbSoftInequalityConstraints, subsetNbAssets, 
											function(i,j) { 
												return Ai[i-1][subsetAssetsIdx[j-1]-1]; 
											});
				}			
				
				// Compute the weights of the minimum tracking error portfolio for the selected assets
				//
				// Note: because the restricted simplex associated to the subset of selected assets
				// might be empty, special care must be taken.
				var subsetAssetsWeights;
				var subsetPortfolioTrackingError = Infinity;
				try {
					var subsetSol = computeMinimumTrackingErrorVolatilityPortfolio(subsetNbAssets, 
																				   subsetAssetsReturns, benchmarkReturns, 
																				   subsetLowerBounds, subsetUpperBounds,
																				   subsetAi, bi);
																				   
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


	// Return the computed weights
	return weights.toArray();
}

