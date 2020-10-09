/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.MeanVarianceEfficientFrontierCla = MeanVarianceEfficientFrontierCla;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function MeanVarianceEfficientFrontierCla
*
* @description Object representing a mean-variance efficient frontier computed using the
* critical line algorithm, c.f. the references, and implementing the methods of the parent virtual
* object MeanVarianceEfficientFrontier.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt optional parameters for the mean-variance optimization problem.
* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm, a strictly positive natural integer 
* or -1 to force an infinite number of iterations; defaults to 1000.
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
*
*/
function MeanVarianceEfficientFrontierCla(mu, sigma, opt) {
	// Call the parent constructor
	MeanVarianceEfficientFrontier.call(this, mu, sigma, opt);
	
	// Compute the corner portfolios defining the efficient frontier
	this.cornerPortfolios = computeCornerPortfolios(this.mu, this.sigma, this.lowerBounds, this.upperBounds, this.epsBounds, opt);
	
	// Safety check
	if (this.cornerPortfolios.length == 0) {
		throw new Error('internal error: no corner portfolio could be computed');
	}
	
	
	// ------


	/**
	* @function computeCornerPortfolios
	*
	* @summary Compute all the corner portfolios belonging to the mean-variance efficient frontier.
	*
	* @description This function returns the weights w_i1,...,w_in as well as the risk aversion parameters lambda_i,
	* i = 1..m, associated to the m fully invested and long-only corner portfolios defining the mean-variance
	* efficient frontier.
	*
	* The algorithm used internally is the Markowitz critical line algorithm, c.f. the first reference.
	*
	* @see Harry M. Markowitz, Mean-Variance Analysis in Portfolio Choice and Capital Markets, Revised issue (2000), McGraw-Hill Professional;
	* @see <a href="https://www.hudsonbaycapital.com/documents/FG/hudsonbay/research/599440_paper.pdf">Harry Markowitz, David Starer, Harvey Fram, Sander Gerber, Avoiding the Downside: A Practical Review of the Critical Line Algorithm for Mean-Semivariance Portfolio Optimization</a>
	*
	* @param {Matrix_} mu the returns of the n assets in the considered universe, a n by 1 matrix.
	* @param {Matrix_} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square n by n matrix.
	* @param {Matrix_} lowerBounds the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n, a n by 1 matrix
	* @param {Matrix_} upperBounds the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n, a n by 1 matrix
	* @param {object} opt optional parameters for the algorithm.
	* @param {number} opt.optimizationMethodParams.maxIterCriticalLine the maximum number of iterations of the critical line algorithm, a strictly positive natural integer; defaults to 1000.
	* @return {Array<Array.<Object>>} the list of all corner portfolios as well as their associated risk aversion parameter, an array made of arrays of two elements:
	* - The corner portfolio weights, a n by 1 Matrix_ of n real numbers
	* - The corner portfolio risk aversion parameter, a positive real number
	*/	
	function computeCornerPortfolios(mu, sigma, lowerBounds, upperBounds, epsBounds, opt) {	
		var eps = 1e-8; // the numerical zero
		
		// Internal object managing the statuses of the asset and lambda variables
		function variablesStatusManager_(nbAssets) {
			// Variables statuses constants
			this.STATUS_UNDEF = -1;
			this.STATUS_IN = 0;
			this.STATUS_LOW = 1;
			this.STATUS_UP = 2;
			
			// The structure holding the variables status
			this.nbAssets = nbAssets;
			
			this.varIn = new BitSet_();
			this.varIn.resize(nbAssets);
			this.varLow = new BitSet_();
			this.varLow.resize(nbAssets);
			this.varUp = new BitSet_();
			this.varUp.resize(nbAssets);
				
			// Public functions to set the status of variables
			this.setIn = function(idx) {
				this.varIn.set(idx);
				this.varLow.unset(idx);
				this.varUp.unset(idx);
			}
			this.setOnLowerBound = function(idx) {
				this.varLow.set(idx);
				this.varIn.unset(idx);
				this.varUp.unset(idx);
			}
			this.setOnUpperBound = function(idx) {
				this.varUp.set(idx);
				this.varIn.unset(idx);
				this.varLow.unset(idx);
			}

			// Public functions to get the status of a variable
			this.isIn = function(idx) {
				return this.varIn.get(idx);
			}
			this.isOnLowerBound = function(idx) {
				return this.varLow.get(idx);
			}
			this.isOnUpperBound = function(idx) {
				return this.varUp.get(idx);
			}
			this.isOut = function(idx) {
				return this.varLow.get(idx) || this.varUp.get(idx);
			}
			
			// Public functions to iterate over the different sets.
			this.getInIndexes = function() {
				return this.varIn.toArray();
			}
			this.getOutIndexes = function() {
				return this.varLow.toArray().concat(this.varUp.toArray());
			}
		}
		
		// Internal function to compute the maximum return efficient portfolio, as well as the status
		// status LOW, IN, UP of the different assets constituting it.
		//
		// To be noted that if there is more than one E-maximizing portfolio, the critical line 
		// algorithm does not support this case "by the book", c.f. chapter 8 of the first reference.
		//
		// A practical workaround to this issue, suggested in chapter 9 of the first reference, 
		// is to slightly alter the assets returns and to relaunch the algorithm.
		//
		// Nevertheless, in function below, the E-maximizing portfolio with minimum variance is computed,
		// so that the starting portfolio is always efficient.
		function computeMaxReturnEfficientPortfolio(mu, sigma, lowerBounds, upperBounds, epsBounds) {
			//
			var nbAssets = sigma.nbRows;
			
			
			// Determine if all assets have different returns, in order to continue 
			// with the determination of the unique efficient max return portfolio with a linear solver,
			// or with a quadratic solver.
			var mu_idx = typeof Uint32Array === 'function' ? new Uint32Array(nbAssets) : new Array(nbAssets);
			for (var j = 0; j < nbAssets; ++j) {		
				mu_idx[j] = j + 1;
			}
			mu_idx.sort(function(a, b) {  // Order the assets in descending order w.r.t. their returns
				return mu.getValue(b, 1) - mu.getValue(a, 1);
			});

			var returnsDifferent = true;
			for (var i = 1; i < nbAssets; ++i) {
				if (mu.getValue(mu_idx[i], 1) == mu.getValue(mu_idx[i-1], 1)) {
					returnsDifferent = false;
					break;
				}
			}

			
			// Compute the maximum return portfolio
			var maxReturnPortfolioWeights;
			if (returnsDifferent == true) {
				var maxReturnSolution = simplexLpSolve_(mu.elemMap(function(i,j,val) { return -val;}), lowerBounds, upperBounds);
				maxReturnPortfolioWeights = new Matrix_(maxReturnSolution[0]);
			}
			else {
				var efficientFrontierGsmo = new MeanVarianceEfficientFrontierGsmo(mu, sigma, {optimizationMethodParams: {antiCyclingGsmo: true, 
				                                                                                                         maximumRiskToleranceValueOnlyGsmo: true,
																							                             maxIterGsmo: -1}, 
				                                                                              constraints: {minWeights: lowerBounds,  maxWeights: upperBounds}});
				maxReturnPortfolioWeights = efficientFrontierGsmo.getHighestReturnPortfolio();
			}
			
			
			// Note: Tests for equality of lower/upper bounds below can be done numerically with a high
			// precision due to the near exact bounds computation of the maxReturnPortfolioWeights above.
					
			
			//
			var nbAssetsIn = 0;
			var nbAssetsLow = 0;
			var nbAssetsUp = 0;
			var variablesStatusManager = new variablesStatusManager_(nbAssets);
			for (var i = 1; i <= nbAssets; ++i) {
				if (Math.abs(maxReturnPortfolioWeights.getValue(i, 1) - lowerBounds.getValue(i, 1)) <= epsBounds) {
					variablesStatusManager.setOnLowerBound(i);
					++nbAssetsLow;
				}
				else if (Math.abs(maxReturnPortfolioWeights.getValue(i, 1) - upperBounds.getValue(i, 1)) <= epsBounds) {
					if (upperBounds.getValue(i, 1) < 1 - epsBounds) { // True UP
						variablesStatusManager.setOnUpperBound(i);
						++nbAssetsUp;
					}
					else { // False UP: upper bound numerically equal to 1
						variablesStatusManager.setIn(i);
						++nbAssetsIn;
					}
				}
				else {
					variablesStatusManager.setIn(i);
					++nbAssetsIn;
				}
			}
			
			
			// Return the computed portfolio weights, as well as associated data
			return {
					weights: maxReturnPortfolioWeights,
					variablesStatusManager: variablesStatusManager,
					nbAssetsLow: nbAssetsLow,
					nbAssetsIn: nbAssetsIn,
					nbAssetsUp: nbAssetsUp,
					};
		}

		
		// ------
		
		
		// TODO: Checks, if enabled

		// Decode options
		if (opt === undefined) {
			opt = { constraints: {} };
		}
		if (opt.optimizationMethodParams === undefined) {
			opt.optimizationMethodParams = {};
		}
		
		// Initialize the options default values
		var maxIterations = opt.optimizationMethodParams.maxIterCriticalLine;
		if (maxIterations == undefined) {
			maxIterations = 1000;
		}
		
		
		// ------

		// Initializations	
		var nbAssets = sigma.nbColumns;
		var lb = lowerBounds;
		var ub = upperBounds;

		var cornerPortfoliosWeights = [];


		// ------
		
		// The only equality constraint supported by the algorithm below
		// is that the weights of the assets must sum to one, but variables
		// below are kept generic.
		var A = Matrix_.ones(1, nbAssets); // the matrix holding the equality constraints
		var nbEqualityConstraints = A.nbRows; // the number of equality constraints 
		var b = Matrix_.ones(1, 1); // the vector holding the right member of the equality constraints
		
		var P = Matrix_.fill(nbAssets, nbAssets + nbEqualityConstraints, 
									function(i,j) { 
										if (j <= nbAssets) {
											return sigma.data[(i-1)*sigma.nbColumns + (j-1)]; // Sigma(i, j)
										}
										else {
											return A.data[(j-nbAssets-1)*A.nbColumns + (i-1)]; // A(j-nbAssets, i) == A(i-nbAssets, j)^t
										}
									});
				
		// ----	

		
		// Step 1: - Compute the rightmost corner portfolio, corresponding to the E-maximizing 
		// portfolio (i.e., the efficient portfolio achieving the maximum return), c.f. chapter 8
		// of the first reference and paragraph 12.3.1 of the second reference.
		//
		//         - Compute as well the status of the different assets constituting it.
		var maxReturnPortfolio = computeMaxReturnEfficientPortfolio(mu, sigma, lb, ub, epsBounds);
		var currentCornerPortfolioWeights = maxReturnPortfolio.weights;
		var variablesStatusManager = maxReturnPortfolio.variablesStatusManager;
		

		// Step 1 bis: manage degeneracy (no asset IN, that is, no asset strictly between its bounds)
		// in the max return efficient portfolio computation.
		//
		// Three cases can occur:
		// - All assets are on LOW (on their lower bounds)
		// - All assets are on UP (on their upper bounds)
		// - All assets are either LOW or UP
		//
		// In the first two cases, the whole efficient frontier 
		// consists of only one portfolio, and it means lower bounds or upper bounds
		// constraints are tight (i.e., they sum to one).
		//
		// In the last case, it means that one or several UP assets are actually IN, but appear UP 
		// due to the "sum to 1" constraint.
		//
		// In order to determine which is (are) this (these) asset(s), the max return portfolio is computed
		// again, with slightly increased upper bounds, so that the "sum to 1" constraint is not binding anymore.
		if (maxReturnPortfolio.nbAssetsLow == nbAssets || maxReturnPortfolio.nbAssetsUp == nbAssets) {
			var weights = new Matrix_(currentCornerPortfolioWeights);
			cornerPortfoliosWeights.push([weights, 0]);

			return cornerPortfoliosWeights;
		}
		else if (maxReturnPortfolio.nbAssetsIn == 0) {
			// To be noted that upper bounds which are already at 1 cannot be increased.
			//
			// Nevertheless, this case has already been managed before, because an asset
			// on its upper bound equal to 1 is actually IN.
			var epsUpperBounds = 1e-8;
			var ubb = ub.elemMap(function(i,j,val) { return Math.min(val + epsUpperBounds, 1);});
			
			// Compute a slightly relaxed max return efficient portfolio, and update the 
			// statuses of the assets accordingly (but not the weights, since the weights
			// computed initially are perfectly fine).
			//
			// In case there is still no asset IN, this is an internal error, since
			// the "sum to 1" constraint cannot be binding anymore due to the updated upper
			// bounds constraints.
			var maxReturnPortfolioB = computeMaxReturnEfficientPortfolio(mu, sigma, lb, ubb, epsBounds);
			variablesStatusManager = maxReturnPortfolioB.variablesStatusManager;
			if (maxReturnPortfolioB.nbAssetsIn == 0) {
				throw new Error("internal error: impossible to determine the IN variables associated to the maximum return portfolio");
			}
		}

		var iter = 0;
		while (true) {
			// Check the number of iterations
			if (maxIterations !== -1 && iter >= maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}

			
			// Update the number of iterations
			++iter;
			
			
			// Step 2: (A13 to A15) of the second reference
			// - Get the current IN/LOW/UP sets
			// - Construct the KKT matrix and related vectors associated to the current IN/LOW/UP sets
			
			// Get the current IN/LOW/UP sets
			var assetsInIdx = variablesStatusManager.getInIndexes();
			var assetsOutIdx = variablesStatusManager.getOutIndexes();
			
			// Construct the matrix Mbar, with Mbar = [[Cbar Abar^t], [Abar 0]],
			// c.f. formula 12 of the second reference, which is incorret, c.f. the accompanying 
			// code and text for the proper formula.
			var Mbar = Matrix_.fill(nbAssets + nbEqualityConstraints, nbAssets + nbEqualityConstraints, 
										function(i,j) { //Cbar
											if (i <= nbAssets && j <= nbAssets) {
												/* 
												This is the incorrect formula, kept for reference
												if (variablesStatusManager.isIn(i) && variablesStatusManager.isIn(j)) {
													return sigma.data[(i-1)*sigma.nbColumns + (j-1)]; // Sigma(i, j)
												}
												else if (!variablesStatusManager.isIn(i) && i == j) {
													return 1;
												}
												else {
													return 0;
												}
												*/
												if (!variablesStatusManager.isIn(i)) {
													if (i == j) {
														return 1;
													}
													else {
														return 0;
													}
												}
												else {
													return sigma.data[(i-1)*sigma.nbColumns + (j-1)]; // Sigma(i, j)
												}									
											}
											else if (i >= nbAssets + 1 && j <= nbAssets) { //Abar
												if (!variablesStatusManager.isIn(j)) {
													return 0;
												}
												else {
													return A.data[(i-nbAssets-1)*A.nbColumns + (j-1)]; // A(i-nbAssets, j)
												}
											}
											else if (i <= nbAssets && j >= nbAssets + 1) { //Abar^t
												if (!variablesStatusManager.isIn(i)) {
													return 0;
												}
												else {
													return A.data[(j-nbAssets-1)*A.nbColumns + (i-1)]; // A(j-nbAssets, i) == A(i-nbAssets, j)^t
												}
											}
											else { // 0
												return 0;
											}
										});	

			// Construct the right hand side vectors associated to alpha and beta vectors,
			// which requires computing the vectors mubar and k, c.f. formulas 11, 12 and 13 of the second reference.
			//
			// To be noted that formula 11 is incorrect, c.f. the associated code and text for the proper formula.
			
			// Right hand side vector associated to alpha
			var k = Matrix_.fill(nbAssets, 1, 
									  function(i,j) {
										  if (variablesStatusManager.isOnUpperBound(i)) {
											  return ub.data[i-1];
										  }
										  else if (variablesStatusManager.isOnLowerBound(i)) {
											  return lb.data[i-1];
										  }
										  else {
											  return 0;
										  }
									  });
			var Ak = Matrix_.xy(A, k);
			var rhsalpha = Matrix_.fill(nbAssets + nbEqualityConstraints, 1, 
											function(i,j) { 
												if (i <= nbAssets) {  
													return k.data[i-1]; /* Incorrect formula is return 0; */
												}
												else {
													return b.data[i-nbAssets-1] - Ak.data[i-nbAssets-1];
												}
											});
			
			// Right hand side vector associated to beta
			var mubar = mu.elemMap(function(i,j,val) { if (!variablesStatusManager.isIn(i)) { return 0; } else { return val; }})
			var rhsbeta = Matrix_.fill(nbAssets + nbEqualityConstraints, 1, 
											function(i,j) { 
												if (i <= nbAssets) {
													return mubar.data[i-1];
												}
												else {
													return 0;
												}
											});
			
			
			// Step 3: (A16) of the second reference
			// - Solve the KKT linear equations in order to compute alpha and beta vectors
			// - Compute gamma and delta vectors
			
			// Compute an LU decomposition of the Mbar matrix
			// 
			// It is guaranteed to be invertible per the critical line algorithm,
			// PROVIDED the first such matrix (i.e., for the IN/LOW/UP set associated to
			// the maximum return efficient portfolio) is invertible, which might not
			// be the case if more than one asset is IN.
			//
			// If Mbar is not invertible, there will be an error at this step.
			try {
				var lu = Matrix_.luDecomposition(Mbar);
				var l = lu.lowerTriangular;
				var u = lu.upperTriangular;
				var p = lu.rowPermutation;
				var q = lu.columnPermutation;

				// Compute alpha
				var z = Matrix_.linsolveForwardSubstitution(lu.lowerTriangular, Matrix_.xy(lu.rowPermutation, rhsalpha));
				var y = Matrix_.linsolveBackSubstitution(lu.upperTriangular, z);
				var alpha = Matrix_.xy(lu.columnPermutation, y);
				
				// Compute beta
				z = Matrix_.linsolveForwardSubstitution(lu.lowerTriangular, Matrix_.xy(lu.rowPermutation, rhsbeta));
				y = Matrix_.linsolveBackSubstitution(lu.upperTriangular, z);
				var beta = Matrix_.xy(lu.columnPermutation, y);
			}
			catch (e) {
				throw new Error('internal error: impossible to solve the KKT system');
			}

			// Compute gamma
			var gamma = Matrix_.xy(P, alpha);
			
			// Compute delta
			var delta = Matrix_.xmy(Matrix_.xy(P, beta), mu);
			
			
			// Step 4: (A17 to A21, A23, A24, A25) of the second reference
			// - Determine the next asset IN to be set OUT
			// - Determine the next asset OUT to be set IN
			// - Compute the current value of lambda_e
			
			// - Determine the next asset IN to be set OUT
			// - Compute lambda_out, c.f. formula 13.17 of the first reference, OR formulas 21 and 20 of the second reference:
			// - lambda_out = max( (L(i) - alpha(i))/beta(i), beta(i) > 0, i in IN; (U(i) - alpha(i))/beta(i), beta(i) < 0, i in IN)
			var idx_out = -1;
			var lambda_out = 0;
			var status_out = variablesStatusManager.STATUS_UNDEF;
			for (var i = 1; i <= assetsInIdx.length; ++i) {
				//
				var in_idx_i = assetsInIdx[i-1];
				var alpha_in_idx_i = alpha.data[in_idx_i-1];
				var beta_in_idx_i = beta.data[in_idx_i-1];

				// Check for asset reaching the lower limit lb
				if (beta_in_idx_i > eps) {
					var lb_idx_in_i = lb.data[in_idx_i-1];
					
					var tmp_lambda_out = (lb_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
					if (tmp_lambda_out >= lambda_out) {
						idx_out = in_idx_i;
						lambda_out = tmp_lambda_out;
						status_out = variablesStatusManager.STATUS_LOW;
					}
				}
				
				// Check for asset reaching the upper limit ub
				else if (beta_in_idx_i < -eps) {
					var ub_idx_in_i = ub.data[in_idx_i-1];
					
					var tmp_lambda_out = (ub_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
					if (tmp_lambda_out >= lambda_out) {
						idx_out = in_idx_i;
						lambda_out = tmp_lambda_out;
						status_out = variablesStatusManager.STATUS_UP;
					}
				}

			}

			// Determine the next asset OUT to be set IN
			// - Compute lambda_in, c.f. formula 13.18 of the first reference, OR formulas 23 and 22 of the second reference:
			// - lambda_in = max( -gamma(i)/delta(i), delta(i) > 0, i in LO; -gamma(i)/delta(i), delta(i) < 0, i in UP)
			var idx_in = -1;
			var lambda_in = 0;
			for (var i = 1; i <= assetsOutIdx.length; ++i) {
				//
				var out_idx_i = assetsOutIdx[i-1];
				var gamma_out_idx_i = gamma.data[out_idx_i-1];
				var delta_out_idx_i = delta.data[out_idx_i-1];
			
				// Check for asset LOW going IN
				if (variablesStatusManager.isOnLowerBound(out_idx_i)) { 
					if (delta_out_idx_i > eps) {
						var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
					
						if (tmp_lambda_in >= lambda_in) {
							idx_in = out_idx_i;
							lambda_in = tmp_lambda_in;
						}
					}
				}
				// Check for asset UP going IN
				else if (variablesStatusManager.isOnUpperBound(out_idx_i)) { 
					if (delta_out_idx_i < -eps) {
						var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
						
						if (tmp_lambda_in >= lambda_in) {
							idx_in = out_idx_i;
							lambda_in = tmp_lambda_in;
						}
					}
				}
				else {
					throw new Error("internal error: OUT variable neither on its lower of upper bound");
				}
			}
				
			// The value of lambda_e for the next corner portfolio is the maximum of
			// the two lambda_out and lambda_in computed above.
			//
			// In case lambda_e == lambda_out, it means an asset first goes OUT as lambda_e
			// is decreased; otherwise, in case lambda_e == lambda_in, it means an asset
			// first goes IN as lambda_e is decreased.
			lambda_e = Math.max(lambda_out, lambda_in, 0);

			
			// Step 5: (A27, A28) of the second reference
			// - Compute the weights of the corner portfolio associated to the current IN/LOW/UP set
			// - Save these weights
			var weights = Matrix_.fill(nbAssets, 1, 
											function(i,j) { 
												if (i <= nbAssets) {
													// Take into account the lower and upper bounds, in
													// case of numerical errors.
													var x_i = alpha.data[i-1] + lambda_e * beta.data[i-1];
													var l_i = lb.data[i-1];
													var u_i = ub.data[i-1];
													
													return Math.max(Math.min(x_i, u_i), l_i);
												}
											});
			cornerPortfoliosWeights.push([weights, lambda_e]);

			
			// Step 6: (A22, A26) of the second reference
			// - Test for the termination criteria of the critical line algorithm
			// - Prepare for the next iteration
			
			// When the value of lambda_e becomes numerically null or negative, the critical
			// line algorithm can be stopped.
			if (lambda_e < eps) {
				break;
			}
			
			// Update the assets statuses with the new IN/OUT asset
			if (lambda_out >= lambda_in) { // an asset IN goes OUT
				// Set the asset idx_out to OUT, with the proper LOW or UP status
				if (status_out == variablesStatusManager.STATUS_LOW) {
					variablesStatusManager.setOnLowerBound(idx_out);
				}
				else {
					variablesStatusManager.setOnUpperBound(idx_out);
				}
			}
			else { // an asset OUT goes IN	
				// Set the asset idx_in as IN
				variablesStatusManager.setIn(idx_in);
			}		
		}
		
		// Return the computed efficient frontier array
		return cornerPortfoliosWeights;	
	}

};
MeanVarianceEfficientFrontierCla.prototype = Object.create(MeanVarianceEfficientFrontier.prototype);
MeanVarianceEfficientFrontierCla.prototype.constructor = MeanVarianceEfficientFrontierCla;

MeanVarianceEfficientFrontierCla.prototype.getHighestRiskTolerancePortfolio = function(x) {
	//
	return this.cornerPortfolios[0][0];
};
MeanVarianceEfficientFrontierCla.prototype.getHighestRiskTolerance = function(x) {
	//
	return this.cornerPortfolios[0][1];
};
MeanVarianceEfficientFrontierCla.prototype.getLowestRiskTolerancePortfolio = function(x) {
	//
	return this.cornerPortfolios[this.cornerPortfolios.length-1][0];
};
MeanVarianceEfficientFrontierCla.prototype.getLowestRiskTolerance = function(x) {
	//
	return this.cornerPortfolios[this.cornerPortfolios.length-1][1];
};


/**
* @function computeEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and
* long-only efficient portfolio belonging to the mean-variance efficient frontier, subject to 
* a return constraint, a volatility constraints or to a risk tolerance constraint.
*
* @memberof MeanVarianceEfficientFrontierCla
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "return", to specify a return constraint
* - "volatility",  to specify a volatility constraint
* - "riskTolerance", to specify a risk tolerance constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
* @return {Array.<Object>} in case no efficient portfolio can be computed, return an empty array;
* in case an efficient portfolio satisfying the input constraint has been computed, an array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontierCla.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	//
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}

	var cornerPortfolioConstraintFct;
	var that = this;
	if (constraintType == "return") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return that.computePortfolioReturn(cornerPortfolio[0]); };
	}
	else if (constraintType == "volatility") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return that.computePortfolioVolatility(cornerPortfolio[0]); };
	}
	else if (constraintType == "riskTolerance") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return cornerPortfolio[1]; };
	}
	else {
		throw new Error('unknown constraint type');
	}
	
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}
	
	// Compute the (at most) two corner portfolios enclosing the efficient 
	// portfolio with a constraint function value equals to the target constraint function value.
	var enclosingCornerPortfolios = computeEnclosingCornerPortfolios(cornerPortfolioConstraintFct, 
	                                                                 constraintValue, 
																	 this.cornerPortfolios,
																	 this.epsEfficientPortfolioComputation);

	
	// Then:
	// - In case the target constraint function value is not reachable, return an empty portfolio 
	//
	// - In case there is a unique computed corner portfolio with a constraint function value
	// equals to the target constraint function value, return the associated portfolio weights
	//
	// - In case there are two corner portfolios enclosing the efficient portfolio with 
	// a constraint function value equals to the target constraint function value, the weights associated 
	// to this efficient portfolio are a convex combination of the weights of the two computed enclosing
	// corner portfolios (c.f. the reference): w = t*w_min + (1-t)*w_max, t in [0,1], with t now to be determined.
	if (enclosingCornerPortfolios.length == 0) {
		return [];
	}
	else if (enclosingCornerPortfolios.length == 1) {
		var idx_min = enclosingCornerPortfolios[0];
		var weights = this.cornerPortfolios[idx_min][0];
		var lambda = this.cornerPortfolios[idx_min][1];
		
		// Return the computed portfolio weights
		return [weights, lambda];
	}
	else {
		// Extract information about the computed efficient corner portfolios
		var idx_min = enclosingCornerPortfolios[0];
		var weights_min = this.cornerPortfolios[idx_min][0];
		var fct_min = cornerPortfolioConstraintFct(this.cornerPortfolios[idx_min]);

		var idx_max = enclosingCornerPortfolios[1];
		var weights_max = this.cornerPortfolios[idx_max][0];
		var fct_max = cornerPortfolioConstraintFct(this.cornerPortfolios[idx_max]);
		
		// Compute t above
		var t;
		
		// If the constraint function value to compute is "return", then, the procedure to 
		// compute t above is the following:
		//
		// E(w) = <mu/w> and by linearity of E, we have
		// E(w) = t*E(w_min) + (1-t)*E(w_max) and E(w) = constraintValue
		// <=>
		// t = (E(w_max) - constraintValue)/(E(w_max) - E(w_min))
		//
		//
		// If the constraint function value to compute is "volatility", then, the procedure to 
		// compute t above is the following:
		//
		// Let the volatility be V(w) = <Sigma*w/w>.
		// Then, by symmetry and bi-linearity of V, V(w) = t^2*V(w_min) + (1-t)^2*V(w_max) + 2*t*(1-t)*<Sigma*w_min/w_max>
		// and V(w) = constraintValue^2
		// <=> t is the solution belonging to ]0,1[ of the second order polynomial equation
		// t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) -2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) - constraintValue^2 = 0
		//
		//
		// If the constraint function value to compute is "riskTolerance", then, the procedure to 
		// compute t above is the following:
		//
		// On the efficient segment [w_min, w_max], it exist vectors alpha and beta such that
		// w_min = alpha + lambda_min * beta
		// w_max = alpha + lambda_max * beta
		// w     = alpha + constraintValue * beta
		// c.f. for instance formula 7.10a of the second reference.
		// From the two first equations, it is possible to deduce alpha = w_min - lambda_a * (w_max - w_min)/(lambda_max - lambda_min)
		// and beta = (w_max - w_min)/(lambda_max - lambda_min), so that w = ...
		// <=>
		// t = (lambda_max - constraintValue)/(lambda_max - lambda_min)
		//
		if (constraintType === "return") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}		
		else if (constraintType === "volatility") {			
			// Define the coefficients of the second order polynomial at^2 + bt + c
			var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
			var a = fct_min * fct_min + fct_max * fct_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
			var b = -2 * (fct_max * fct_max - variance_cross); // 
			var c = fct_max * fct_max - constraintValue*constraintValue; //always > 0
			
			// Extract the root t of the equation at^2 + bt + c = 0 belonging to ]0,1[, using a stable numerical formula
			var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
				throw new Error('internal error; the covariance matrix might not be semi-definite positive');
			}
			var q = -(b_p + sign_b_p * Math.sqrt(disc));
			var r1 = q/a;
			var r2 = c/q;
			
			if (r1 > 0 && r1 < 1) {
				t = r1;
			}
			else if (r2 > 0 && r2 < 1) {
				t = r2;
			}
			else {
				throw new Error('internal error: the covariance matrix might not be semi-definite positive');
			}
		}
		else if (constraintType === "riskTolerance") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}
		else {
			throw new Error('internal error: unknown constraint type');
		}

		// Compute the final efficient portfolio weights
		var that = this;
		var weights = Matrix_.fill(weights_min.nbRows, 1,
							   	   function(i,j) { 
										// Take into account the lower and upper bounds, in
										// case of numerical errors.
										var x_i = t*weights_min.getValue(i, 1) + (1-t)*weights_max.getValue(i, 1);
										var l_i = that.lowerBounds.data[i-1];
										var u_i = that.upperBounds.data[i-1];
										
										return Math.max(Math.min(x_i, u_i), l_i);
								   });
		
		// Compute the associated risk tolerance parameter
		var lambda = t*this.cornerPortfolios[idx_min][1] + (1-t)*this.cornerPortfolios[idx_max][1];
		
		// Return the computed portfolio weights
		return [weights, lambda];	
	}
	
	
	// ------
	
	
	// Internal function to compute the (at most) two corner portfolios enclosing the
	// efficient portfolio with a target constraint value (return, volatility or risk tolerance),
	// using a binary search algorithm.
	//
	// The usage of a binary search algorithm is justified because the corner portfolios
	// return, volatility and risk tolerance are decreasing as soon as there are at least two corner
	// portfolios on the efficient frontier.
	function computeEnclosingCornerPortfolios(fct, fctValue, cornerPortfolios, eps) {		
		// The efficient frontier portfolios are provided from highest return/volatility/risk tolerance
		// to lowest return/volatility/risk tolerance, so that *_min below refers to properties of the portfolio
		// with the lowest return/volatility/risk tolerance.
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];
		
		var fct_min = fct(cornerPortfolios[idx_min]);
		var fct_max = fct(cornerPortfolios[idx_max]);

		// If the desired target constraint function value is numerically strictly greater than the highest attainable
		// constraint function value on the efficient frontier, or if the desired target constraint function 
		// value is numerically strictly lower than the lowest attainable constraint function value on the efficient frontier, 
		// there is no enclosing corner portfolio on the efficient frontier.
		if (fctValue > fct_max + eps || fctValue < fct_min - eps) {
			return [];
		}

		// Otherwise, if the target constraint function value is numerically reached on one of the
		// two extremal corner portfolios, return immediately.
		if (Math.abs(fctValue - fct_min) <= eps) {
			return [idx_min];
		}
		else if (Math.abs(fctValue - fct_max) <= eps) {
			return [idx_max];
		}
		
		// Otherwise, determine the two adjacent corner portfolios enclosing the portfolio
		// with a constraint function value numerically equals to the target constraint function value, 
		// using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle point
			var idx_middle = Math.floor((idx_min + idx_max)/2);
			var fct_middle = fct(cornerPortfolios[idx_middle]);
			
			// Determine in which sub-interval ]idx_max, idx_middle[ or ]idx_middle, idx_min[
			// lies the portfolio with the target constraint function value.
			if (fct_middle > fctValue) {
				idx_max = idx_middle;
			}
			else if (fct_middle < fctValue) {
				idx_min = idx_middle;
			}
			else { // the target constraint function value is exactly attained on the idx_middle-th corner portfolio
				return [idx_middle];
			}
		}

		// Return the computed adjacent corner portfolios.
		return [idx_min, idx_max];
	}
};


/**
* @function getCornerPortfolios
*
* @description This function returns the weights w_i1,...,w_in associated to the m fully invested and
* long-only corner portfolios defining the mean-variance efficient frontier.
*
* @memberof MeanVarianceEfficientFrontierCla
*
* @return {Array<Array.<Object>>} the list of all corner portfolios, an array of arrays of n by 1 matrices
*
*/
MeanVarianceEfficientFrontierCla.prototype.getCornerPortfolios = function() {
	//
	var portfolios = new Array(this.cornerPortfolios.length);
	
	//
	for (var i = 0; i < portfolios.length; ++i) {
		portfolios[i] = new Matrix_(this.cornerPortfolios[i][0]);
	}
	
	//
	return portfolios;
};


/**
* @function restrict
*
* @description This function restricts the efficient frontier to portfolios satisfying 
* a minimal return constraint or a minimal volatility constraint.
*
* @memberof MeanVarianceEfficientFrontierCla
*
* @param {string} constraintType, the type of constraint, a string either equal to:
* - "minReturn", to specify a return constraint
* - "minVolatility",  to specify a volatility constraint
* @param {number} constraintValue, the value of the constraint, a real number.
*
*/
MeanVarianceEfficientFrontierCla.prototype.restrict = function(constraintType, constraintValue) {
	// Decode the input parameters
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}

	// Restrict the efficient frontier, if possible
	if (constraintType == "minReturn") {	
		var lowestReturn = this.getLowestReturn();
		if (lowestReturn < constraintValue) {
			// The efficient frontier needs to be restricted, and this might not be possible
			var highestReturn = this.getHighestReturn();
			if (highestReturn < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum return constraint is not feasible');
			}
			
			// Compute the first efficient portfolio with a strictly positive enough return
			var efficientPortfolio = this.computeEfficientPortfolio("return", constraintValue + 0.5 * this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a return greater than ' + constraintValue);
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			
			// Rebuild the list of corner portfolios, removing all the corner portfolios 
			// with a return lower than the desired minimum return, and adding
			// the first efficient portfolio with a strictly positive enough return computed
			// above.
			var updatedCornerPortfolios = new Array();
			for (i = 0; i < this.cornerPortfolios.length - 1; ++i) {	
				if (this.cornerPortfolios[i][1] > efficientPortfolioRiskTolerance) {
					updatedCornerPortfolios.push(this.cornerPortfolios[i]);
				}
				else {
					break;
				}
			}
			updatedCornerPortfolios.push([efficientPortfolioWeights, efficientPortfolioRiskTolerance]);
			this.cornerPortfolios = updatedCornerPortfolios;
		}
		
	}
	else if (constraintType == "minVolatility") {
		var lowestVolatility = this.getLowestVolatility();
		if (lowestVolatility < constraintValue) {
			// The efficient frontier needs to be restricted, and this might not be possible
			var highestVolatility = this.getHighestVolatility();
			if (highestVolatility < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum volatility constraint is not feasible');
			}
			
			// Compute the first efficient portfolio with a strictly positive enough volatility
			var efficientPortfolio = this.computeEfficientPortfolio("volatility", constraintValue + 0.5 * this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a volatility greater than ' + constraintValue + ': the covariance matrix might not be semi-definite positive');
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];

			// Rebuild the list of corner portfolios, removing all the corner portfolios 
			// with a volatility lower than the desired minimum volatility, and adding
			// the first efficient portfolio with a strictly positive enough volatility computed
			// above.
			var updatedCornerPortfolios = new Array();
			for (i = 0; i < this.cornerPortfolios.length - 1; ++i) {	
				if (this.cornerPortfolios[i][1] > efficientPortfolioRiskTolerance) {
					updatedCornerPortfolios.push(this.cornerPortfolios[i]);
				}
				else {
					break;
				}
			}
			updatedCornerPortfolios.push([efficientPortfolioWeights, efficientPortfolioRiskTolerance]);
			this.cornerPortfolios = updatedCornerPortfolios;
		}
	}
	else {
		throw new Error('unknown constraint type');
	}
	
	// Nothing to do here
};


/**
* @function computeMaximumSharpeRatioEfficientPortfolio
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and 
* long-only portfolio of n assets maximizing the Sharpe ratio, which is defined as the ratio of the 
* portfolio excess return over a constant risk-free rate to the portfolio volatility.
*
* @see <a href="https://doi.org/10.1111/j.1540-6261.1976.tb03217.x">Elton, E. J., Gruber, M. J. and Padberg, M. W. (1976), SIMPLE CRITERIA FOR OPTIMAL PORTFOLIO SELECTION. The Journal of Finance, 31: 1341-1357</a>
* @see <a href="http://dx.doi.org/10.1080/13504860701255292">S. V. Stoyanov , S. T. Rachev & F. J. Fabozzi (2007) Optimal Financial Portfolios, Applied Mathematical Finance, 14:5, 401-436</a>
*
* @memberof MeanVarianceEfficientFrontierCla
*
* @param {number} rf the risk-free rate, a real number.
*
* @return {Array.<Object>} An array of 2 elements: 
* -- The computed efficient portfolio weights
* -- The risk tolerance associated to the computed portfolio
*/
MeanVarianceEfficientFrontierCla.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	// The efficient frontier must be restricted to both:
	// - The domain of definition of the Sharpe ratio (the portfolios with a strictly positive volatility)
	// - The domain of strict positivity of the Sharpe ratio

	// On the restricted efficient frontier, the Sharpe ratio is a pseudo-concave 
	// function, c.f. the first or second reference, so that it is unimodal.
	//
	// This property allows:
	// - To first search for the corner portfolio with the maximum
	// Sharpe ratio using a binary search algorithm.
	//
	// - Then, because the corner portfolio with the maximum Sharpe ratio is adjacent to at most
	// two other corner portfolios, to search for the efficient portfolio with the maximum
	// Sharpe ratio.
	//
	//  Indeed, depending on the position of the corner portfolio with the maximum Sharpe
	// ratio on the efficient frontier:
	// - Unique corner portfolio => zero adjacent corner portfolio
	// - Non unique leftmost or rightmost corner portfolio => one adjacent corner portfolio
	// - Non unique any other corner portfolio => two adjacent corner portfolios
	//
	// In the first case, the efficient portfolio with the maximum Sharpe ratio
	// is the same as the corner portfolio with the maximum Sharpe ratio
	//
	// In the last two cases, because of the pseudo concavity of the Sharpe ratio
	// on the restricted efficient frontier, the efficient portfolio with the maximum
	// Sharpe ratio is guaranteed to belong to the efficient segment(s) connecting
	// the corner portfolio with the maximum Sharpe ratio to its adjacent corner
	// portfolio(s).
	//
	// So, computing the efficient portfolio with the maximum Sharpe ratio is equivalent
	// to computing the efficient portfolio with the maximum Sharpe ratio on the efficient
	// segment(s) connecting the corner portfolio with the maximum Sharpe ratio
	// to its adjacent corner portfolio(s).
	var idx = computeMaximumSharpeRatioCornerPortfolio.call(this, rf);

	// Add the corner portfolio with the maximum Sharpe ratio as a candidate 
	// for being the efficient portfolio with the maximum Sharpe ratio.
	var weights_idx = this.cornerPortfolios[idx][0];
	var sr_idx = this.computePortfolioSharpeRatio(weights_idx, rf);
	var candidatePortfolios = [[weights_idx, sr_idx, this.cornerPortfolios[idx][1]]];
		
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// on the efficient segment [idx+1, idx], if existing.
	if (idx <= this.cornerPortfolios.length - 2) {
		candidatePortfolios.push( computeMaximumSharpeRatioEfficientSegmentPortfolio.call(this, rf, idx + 1, idx) );
	}
	
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// on the efficient segment [idx, idx-1], if existing.	
	if (idx >= 1) {
        candidatePortfolios.push( computeMaximumSharpeRatioEfficientSegmentPortfolio.call(this, rf, idx, idx - 1) );
	}
	
	// Compute the efficient portfolio maximizing the Sharpe ratio
	// by merging the efficient portfolios locally maximizing
	// the Sharpe ratio on each efficient segment.
	var compareSharpeRatios = function (a, b) {
		return a[1] - b[1];
	};
	var maxSharpeRatioPortfolio = max_(candidatePortfolios, compareSharpeRatios)[0];
	
	// Return this portfolio, as well as its risk tolerance parameter
	return [maxSharpeRatioPortfolio[0], maxSharpeRatioPortfolio[2]];
	
	
	// ------
	
	// Internal function to compute the corner portfolio which maximizes the
	// Sharpe ratio on the efficient frontier restricted to portfolios with:
	// - A strictly positive volatility
	// - A strictly positive excess return
	//
	// This function uses a binary search algorithm, which is justified because
	// the Sharpe ratio is a pseudo concave function on its domain of strict positivity,
	// so that it is a unimodular function, c.f. the first or second reference.
	function computeMaximumSharpeRatioCornerPortfolio(rf) {
		// The efficient frontier portfolios are provided from highest return/volatility
		// to lowest return/volatility, so that *_min below refers to properties of the portfolio
		// with the lowest return/volatility.
		var idx_min = this.cornerPortfolios.length - 1;
		var idx_max = 0;
		
		// In case there is only one corner portfolio on the efficient frontier,
		// exit immediately.
		if (idx_min == idx_max) {
			return idx_min;
		}
		
		// Otherwise, determine the corner portfolio with the maximum Sharpe ratio 
		// using a binary search algorithm.
		while (idx_min - idx_max != 1) { 
			// Compute properties on the middle points
			var idx_middle = Math.floor((idx_min + idx_max)/2);
			var weights_middle = this.cornerPortfolios[idx_middle][0];
			var sharpeRatio_middle = this.computePortfolioSharpeRatio(weights_middle, rf);

			var idx_middle_p = idx_middle + 1; 
			var weights_middle_p = this.cornerPortfolios[idx_middle_p][0];
			var sharpeRatio_middle_p = this.computePortfolioSharpeRatio(weights_middle_p, rf);

			// Determine in which sub-interval [idx_max, idx_middle] or [idx_middle, idx_min]
			// lies the corner portfolio with the maximum Sharpe ratio.
			if (sharpeRatio_middle > sharpeRatio_middle_p) {
				idx_min = idx_middle;
			}
			else if (sharpeRatio_middle	< sharpeRatio_middle_p) {
				idx_max = idx_middle;
			}
			else {
				// In case the Sharpe ratio is equal on both corner portfolios, 
				// it means its maximum is attained somewhere between these two portfolios, 
				// due to its strict unimodality.
				//
				// The binary search procedure can then be prematurely stopped, although
				// this case is (numerically) highly improbable.
				idx_min = idx_middle_p;		
				idx_max = idx_middle;

				break;
			}
		}
		
		// Return the computed corner portfolio index
		return idx_min;
	}

	// Internal function to compute the efficient portfolio which maximizes the
	// Sharpe ratio on an efficient segment defined by two adjacent corner portfolios
	// with:
	// - A strictly positive volatility
	// - A strictly positive excess return
	//
	// On such an efficient segment, the weights associated this portfolio are a 
	// convex combination of the weights of the two adjacent corner portfolios,
	// so that w = t*w_min + (1-t)*w_max, t in [0,1], with t to be determined,
	// c.f. the third reference.
	//
	// With E(w) = <mu/w> the portfolio return and V(w) = <Sigma*w/w> the portfolio 
	// variance, the Sharpe ratio is defined as SR(w) = (E(w) - rf)/SQRT(V(w)).
	//
	// Because SR(w) > 0 on the efficient segment, maximizing SR(w) is equivalent
	// to maximizing SR(w)^2, which is equal to (E(w) - rf)^2/V(w).
	//
	// By linearity of E(w) and bilinearity/symmetry of V(w), SR(w)^2 is also equal to
	// a rational fraction in t:
	//
	// (E(w) - rf)^2/V(w)
	// =
	// ( E(t*w_min + (1-t)*w_max) - rf )^2 / ( V(t*w_min + (1-t)*w_max) )
	// =
	// ( t*(E(w_min) - E(w_max)) + E(w_max) - rf )^2 / ( t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) - 2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) )
	// = ( t^2*(E(w_min) - E(w_max))^2 + 2*(E(w_min) - E(w_max))*(E(w_max) - rf) + (E(w_max) - rf)^2 ) / ( t^2*(V(w_min) + V(w_max) - 2*<Sigma*w_min/w_max>) - 2*t*(V(w_max) - <Sigma*w_min/w_max>) + V(w_max) )
	//
	// So, maximizing SR(w) on the efficient segment is equivalent to maximizing
	// SR(t)^2, t in [0,1].
	//
	// Since SR(t)^2 is a differentiable function on [0,1] and since [0,1] is a closed convex set,
	// its maximum is either reached on its boundary (i.e., {0,1}) or on a critical interior point
	// (i.e., a point belonging to ]0,1[ on which the derivative of SR(t)^2 vanishes).
	//
	// Evaluating SR(t) on each of these (at most) four points and selecting t
	// as the value which maximizes SR(t) then allows to compute the weights 
	// of the efficient portfolio which maximizes the Sharpe ratio on the efficient segment.
	function computeMaximumSharpeRatioEfficientSegmentPortfolio(rf, idx_min, idx_max) {
		//
		var cornerPortfolios = this.cornerPortfolios;
		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];
		
		// Compute properties of the two adjacent corner portfolios
		var sr_min = this.computePortfolioSharpeRatio(weights_min, rf);
		var return_min = this.computePortfolioReturn(weights_min);
		var volatility_min = this.computePortfolioVolatility(weights_min);
		var variance_min = volatility_min * volatility_min;
		
		var sr_max = this.computePortfolioSharpeRatio(weights_max, rf);
		var return_max = this.computePortfolioReturn(weights_max);
		var volatility_max = this.computePortfolioVolatility(weights_max);
		var variance_max = volatility_max * volatility_max;
		
		// Define the coefficients of the fractional function SR(t)^2 = ( at^2 + bt + c ) / ( dt^2 + et + f )
		var return_min_m_max = return_min - return_max;
		var return_max_m_rf = return_max - rf;
		var a = return_min_m_max * return_min_m_max;
		var b = 2 * return_min_m_max * return_max_m_rf;
		var c = return_max_m_rf * return_max_m_rf;
		
		var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
		var d = variance_min + variance_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
		var e = -2 * (variance_max - variance_cross); // 
		var f = variance_max; //always > 0
		
		// Define the coefficients of the second order polynomial aat^2 + bbt + cc equal to the
		// numerator of the derivative d(SR(t)^2)/dt.
		var aa = a*e - b*d;
		var bb = 2*(a*f - c*d);
		var cc = b*f - c*e;
		
		// Extract the roots t1 and t2 of the equation d(SR(t)^2)/dt = 0, using a stable numerical formula.
		var bb_p = bb/2; // reduced discriminant
		var sign_bb_p = (bb_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for bb_p == 0 this returns 1
		var disc = bb_p*bb_p - aa*cc;
		if (Math.abs(disc) <= 1e-16) { // In case of numerically semi-positive definite covariance matrix
			disc = 0;
		}
		if (disc < 0) {
			throw new Error('internal error, the covariance matrix might not be semi-definite positive');
		}
		var qq = -(bb_p + sign_bb_p * Math.sqrt(disc));
		var t1 = qq/aa;
		var t2 = cc/qq;
		
		// Compute and order the Sharpe ratios for all the efficient 
		// portfolios with t corresponding to {0, 1, t1, t2}.
		var candidateSharpeRatios = [[weights_min, sr_min, cornerPortfolios[idx_min][1]], [weights_max, sr_max, cornerPortfolios[idx_max][1]]]; // t = 0 and t = 1 portfolios are always present
		
		if (t1 > 0 && t1 < 1) { // t1 belongs to ]0,1[
			var weights_t1 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t1*weights_min.getValue(i, 1) + (1-t1)*weights_max.getValue(i, 1);
										})
			var sr_t1 = this.computePortfolioSharpeRatio(weights_t1, rf);
			var lambda_t1 = t1*cornerPortfolios[idx_min][1] + (1-t1)*cornerPortfolios[idx_max][1];

			candidateSharpeRatios.push([weights_t1, sr_t1, lambda_t1]);
		}

		if (t2 > 0 && t2 < 1) { // t2 belongs to ]0,1[
			var weights_t2 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t2*weights_min.getValue(i, 1) + (1-t2)*weights_max.getValue(i, 1);
										})
			var sr_t2 = this.computePortfolioSharpeRatio(weights_t2, rf);
			var lambda_t2 = t2*cornerPortfolios[idx_min][1] + (1-t2)*cornerPortfolios[idx_max][1];
			
			candidateSharpeRatios.push([weights_t2, sr_t2, lambda_t2]);
		}

		
		// Return the efficient portfolio which maximizes the Sharpe ratio
		// on the efficient segment.
		var compareSharpeRatios = function (a, b) {
			return a[1] - b[1];
		};
		return max_(candidateSharpeRatios, compareSharpeRatios)[0];
	}
};