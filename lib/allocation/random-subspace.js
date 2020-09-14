/**
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function randomSubspaceOptimizationWeights
*
* @summary Compute the weights of a portfolio using the random subspace optimization method 
* applied to an arbitrary portfolio optimization algorithm.
*
* @description This function returns the weights w_1,...,w_n associated to a long-only fully or partially 
* invested portfolio of n assets, as computed by the random subspace optimization method described informally 
* in the first reference and more formally in the second and third references (in the specific case of 
* the mean-variance portfolio optimization algorithm).
*
* Optionally, the following constraints can be added:
* - Minimum weight of each asset that is included in the portfolio
* - Maximum weight of each asset that is included in the portfolio
*
* This algorithm combines the usage of a random subspace optimization method with an arbitrary portfolio
* optimization algorithm the following way:
* - If subsetsGenerationMethod is equal to 'random', repeat nbRandomSubsets times:
* -- Generate uniformly at random a subset of sizeSubsets assets from the n assets 
* -- Using an arbitrary portfolio optimization algorithm over the generated subset of assets,
* compute the weights associated to an optimal portfolio
*
* - Else if subsetsGenerationMethod is equal to 'deterministic', repeat Binomial(nbAssets, sizeSubsets) times:
* -- Generate a subset of sizeSubsets assets from the n assets, without replacement
* -- Using an arbitrary portfolio optimization algorithm over the generated subset of assets,
* compute the weights associated to an optimal portfolio
*
* Note: in both cases, it can happen that the portfolio optimization algorithm is not able to
* compute an optimal portfolio over some generated subsets of assets, for instance because
* there is no feasible portfolio over these subsets of assets. The parameter maxInfeasibleSubsetsRatio
* allows to define a maximum allowed proportion of such infeasible subsets over the total number of generated subsets,
* beyond which the random subspace optimization method is aborted.
*
* - In both cases, compute the weights of a final portfolio as 
* either the arithmetic average of all the previously computed portfolios weights
* (if subsetsAggregationMethod is equal to 'average'), which is ex-ante optimal 
* c.f. the third reference, or as the geometric median of all the previously computed 
* portfolios weights (if subsetsAggregationMethod is equal to 'median'), which is more robust to
* a bad luck of the draw than the average.
*
* Note: The algorithm used internally for the uniform selection at random is the method D
* of J.S. Vitter, c.f. the sixth reference.
*
* @see <a href="https://cssanalytics.wordpress.com/2013/10/10/rso-mvo-vs-standard-mvo-backtest-comparison/">RSO MVO vs Standard MVO Backtest Comparison, OCTOBER 10, 2013</a>
* @see <a href="https://aaai.org/ocs/index.php/AAAI/AAAI17/paper/view/14443">SHEN, W.; WANG, J.. Portfolio Selection via Subset Resampling. AAAI Conference on Artificial Intelligence, North America, feb. 2017.</a>
* @see <a href="http://www.hss.caltech.edu/content/subset-optimization-asset-allocation">Benjamin J. Gillen, Subset Optimization for Asset Allocation, SOCIAL SCIENCE WORKING PAPER 1421, June 1, 2016</a>
* @see <a href="https://doi.org/10.1007/978-3-642-31537-4_13">Oshiro T.M., Perez P.S., Baranauskas J.A. (2012) How Many Trees in a Random Forest?. In: Perner P. (eds) Machine Learning and Data Mining in Pattern Recognition. MLDM 2012. Lecture Notes in Computer Science, vol 7376. Springer, Berlin, Heidelberg</a>
* @see <a href="https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf">Breiman, L (2002), Manual On Setting Up, Using, And Understanding Random Forests V3.1</a>
* @see J.S. Vitter. An efficient algorithm for sequential random sampling. RR-0624, INRIA. 1987. <inria-00075929>
*
* @param {number} nbAssets the number of assets in the considered universe, natural integer superior or equal to 1.
* @param {function} subsetOptimizationFct, a function representing a portfolio optimization method to apply on each generated subset of assets, which must take:
* - As a first input argument, an array subsetAssetsIdx of sizeSubsets different integers belonging to [1..nbAssets], 
* in increasing order, representing the indexes of the original assets belonging to the generated subset of assets
* - As a second input argument, a JavaScript object subsetOptimizationFctOpt, representing optional parameters
* to be provided to the function subsetOptimizationFct, with subsetOptimizationFctOpt.constraints.minWeights and 
* subsetOptimizationFctOpt.constraints.maxWeights automatically computed from
* opt.constraints.minWeights and opt.constraints.maxWeights if provided
* and which must return an array of sizeSubsets real numbers belonging to the unit simplex of R^sizeSubsets, representing
* the weights w_1,...,w_sizeSubsets of the portfolio computed by the portfolio optimization method applied to the generated 
* subset of assets. In case these weights cannot be computed, the funtion subsetOptimizationFct is expected to throw an instance of
* the JavaScript Error class with the exact error message "infeasible portfolio optimization problem".
* @param {object} opt optional parameters for the random subspace optimization method.
* @param {number} opt.sizeSubsets the number of assets to include in the generated subsets of assets, a positive natural integer satisfying 2 <= sizeSubsets < n; 
* defaults to the floored value of SQRT(nbAssets).
* @param {string} opt.subsetsGenerationMethod the method used to generate the subset of assets, a string either equal to:
* - 'random' in order to generate the subsets of assets uniformly at random
* - 'deterministic' in order to generate the subsets of assets deterministically, through the enumeration of all the Binomial(nbAssets, sizeSubsets) subsets of assets
*; defaults to 'random'
* @param {number} opt.nbRandomSubsets the number of subsets of assets to generate in case opt.subsetsGenerationMethod is set to 'random', a strictly positive natural integer; defaults to 128.
* @param {number} opt.maxInfeasibleSubsetsRatio the maximum allowed proportion of infeasible subsets of assets over all the generated subsets of assets,
* a positive real number satisfying 0 <= maxInfeasibleSubsetsRatio < 1; defaults to 0.
* @param {string} opt.subsetsAggregationMethod the method used to compute the final portfolio weights from all the generated portfolios weights,
* a string equal to:
* - 'average' in order to compute the final portfolio weights as the arithmetic average of all the computed portfolios weights
* - 'median' in order to compute the final portfolio weights as the geometric median of all the computed portfolios weights
; defaults to 'average'.
* @param {object} opt.subsetOptimizationFctOpt optional parameters to be provided as is to the portfolio optimization method represented by the function subsetOptimizationFct, a JavaScript object 
* @param {Array.<number>} opt.constraints.minWeights an optional array of size n (l_i),i=1..n containing the minimum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of zeros.
* @param {Array.<number>} opt.constraints.maxWeights an optional array of size n (u_i),i=1..n containing the maximum weights of the assets that are included in the portfolio with 0 <= l_i <= u_i <= 1, i=1..n; defaults to an array made of ones.
* @return {Array.<number>} the weights corresponding to the computed portfolio, array of n real numbers.
*
* @example
* randomSubspaceOptimizationWeights(3,
									function subsetOptimizationFct(subsetAssetsIdx, subsetOptimizationFctOpt) {
										// Return a constant value, independently of the selected assets
										return [0, 1];
									},
									{sizeSubsets: 2})
*/
self.randomSubspaceOptimizationWeights = function(nbAssets, subsetOptimizationFct, opt) {			
	// Initialize the options structure
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetOptimizationFctOpt === undefined) {
		opt.subsetOptimizationFctOpt = {};
	}
	if (opt.subsetOptimizationFctOpt.constraints === undefined) {
		opt.subsetOptimizationFctOpt.constraints = {};
	}

	// ------
	
	// Limit cases: 
	// - If the number of assets is lower than or equal to 1, return immediately
	if (nbAssets <= 1) {
		return [1];
	}
	
	
	// ------	
	
	// Decode options

	// The default size of the subsets to generate is obtained following 
	// the fifth reference, i.e., the square root of the number of assets.
	var sizeSubsets = opt.sizeSubsets;
	if (sizeSubsets === undefined) {
		sizeSubsets = Math.max(2, Math.floor(Math.sqrt(nbAssets)));
	}
	if (sizeSubsets <= 1 || sizeSubsets >= nbAssets + 1) {
		throw new Error('the size of the subsets of assets must be between 2 and ' + nbAssets.toString());
	}

	
	// The subsets optimization method is equal to the function subsetOptimizationFct

	// The default method of subsets generation is uniform at random 
	var subsetsGenerationMethod = opt.subsetsGenerationMethod;
	if (subsetsGenerationMethod === undefined) {
		subsetsGenerationMethod = 'random';
	}
	
	// The default number of random subsets to generate is 128,
	// c.f. the fourth reference.
	//
	// Limit cases:
	// - If the size of the subsets to generate is equal to the number of assets,
	// there would be nbSubsets identical computations made in the core process below,
	// to which their average/median would also be identical.
	//
	// Thus, this case is explicitly managed for performances reasons.
	var nbRandomSubsets = opt.nbRandomSubsets;
	if (nbRandomSubsets === undefined) {
		nbRandomSubsets = 128;		
	}
	if (sizeSubsets === nbAssets) {
		nbRandomSubsets = 1;
	}

	
	// The default number of subsets to generate depends on the method of
	// subsets generation.
	var nbSubsets;
	if (subsetsGenerationMethod === 'random') {
		nbSubsets = nbRandomSubsets;
	}
	else if (subsetsGenerationMethod === 'deterministic') {
		nbSubsets = binomial_(nbAssets, sizeSubsets);
	}
	else {
		throw new Error('unsupported subsets of assets generation method');
	}
	
	// The default method of aggregation of the random portfolios into the
	// final portfolio.
	var subsetsAggregationMethod = opt.subsetsAggregationMethod;
	if (subsetsAggregationMethod === undefined) {
		subsetsAggregationMethod =  'average';
	}
	if (subsetsAggregationMethod !== 'average' && 
	    subsetsAggregationMethod !== 'median') {
		throw new Error('unsupported aggregation method');
	}
	
	// The default proportion of allowed infeasible generated portfolios over 
	// all the generated portfolios is zero.
	var maxInfeasibleSubsetsRatio = opt.maxInfeasibleSubsetsRatio;
	if (maxInfeasibleSubsetsRatio === undefined) {
		maxInfeasibleSubsetsRatio =  0;
	}

	
	// ------
	
	
	// Core process
	
	// Initializations
	// The assets subsets generator
	var subsetAssetsIdxIterator; 
	if (subsetsGenerationMethod === 'random') {
		subsetAssetsIdxIterator = new randomKSubsetsIterator_(nbAssets, sizeSubsets, false); // use no array copy in the subsets generation to improve performances
	}
	else if (subsetsGenerationMethod === 'deterministic') {
		subsetAssetsIdxIterator = new kSubsetsIterator_(nbAssets, sizeSubsets, false); // use no array copy in the subsets generation to improve performances
	}
	else {
		throw new Error('unsupported subsets generation method');
	}
	
	// The options to provide to the portfolio optimization algorithm used
	// to compute the weights associated to the selected assets.
	var subsetOptimizationFctOpt = opt.subsetOptimizationFctOpt;

    // The number of generated feasible portfolios
	var nbFeasibleGeneratedPortfolios = 0;

	// The storage space for the generated portfolios weights
	var generatedPortfoliosWeights = new Array(nbSubsets); 

	// The storage space for the optional subsets minWeights and maxWeights constraints.
	if (opt.constraints.minWeights) {
		subsetOptimizationFctOpt.constraints.minWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubsets) : new Array(sizeSubsets);
	}
	if (opt.constraints.maxWeights) {
		subsetOptimizationFctOpt.constraints.maxWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubsets) : new Array(sizeSubsets);
	}
	
	// Generation of the nbSubsets portfolios weights
	for (var k = 0; k < nbSubsets; ++k) {
		// Select either uniformly at random or deterministically sizeSubsets assets
		// from the nbAssets assets.
		//
		// Note: the ex-ante optimality of the final portfolio (if the aggregation mode is
		// 'average') relies on the fact that the subsets of assets are generated uniformly,
		// c.f. the third reference.
		//
		// The "uniformness" of the random algorithm used is then very important.
		var subsetAssetsIdx = subsetAssetsIdxIterator.next();
		
		// In case minimum/maximum weights constraints are provided, automatically map
		// these constraints to the selected assets.
		if (opt.constraints.minWeights) {
			for (var i = 0; i < sizeSubsets; ++i) {
				subsetOptimizationFctOpt.constraints.minWeights[i] = opt.constraints.minWeights[subsetAssetsIdx[i]-1];
			}
		}
		if (opt.constraints.maxWeights) {
			for (var i = 0; i < sizeSubsets; ++i) {
				subsetOptimizationFctOpt.constraints.maxWeights[i] = opt.constraints.maxWeights[subsetAssetsIdx[i]-1]
			}
		}
		
		// Compute the optimal portfolio for the selected assets using the 
		// provided portfolio optimization algorithm.
		//
		// In case the weights cannot be computed due to the infeasibility of the problem
		// (for instance, because of lower/upper bounds), another portfolio is generated.
		try {
			var subsetWeights = subsetOptimizationFct(subsetAssetsIdx, subsetOptimizationFctOpt);
			
			if (subsetWeights.length != sizeSubsets) {
				throw new Error('internal error: the portfolio optimization method did not return the expected number of weights');
			}
		}
		catch (e) {
			if (e.message === "infeasible portfolio optimization problem") {
				continue;
			}
			else {
				throw(e);
			}
		}

		// Transform the computed weights for the selected assets into their equivalent 
		// computed weights for the original assets (e.g. adding zero weights on non-selected
		// assets).
		var weights = Matrix_.zeros(nbAssets, 1);
		for (var i = 0; i < sizeSubsets; ++i) {
			weights.setValueAt(subsetAssetsIdx[i], 1, 
							   subsetWeights[i]);
		}

		// Save the resulting original assets weights
		generatedPortfoliosWeights[nbFeasibleGeneratedPortfolios++] = weights;
	}
	
	// Computation of the final portfolio weights

	// Resize of the storage space for the generated portfolios weights
	generatedPortfoliosWeights.length = nbFeasibleGeneratedPortfolios;

	// In case there was no generated portfolios, because they were all
	// infeasible, abort the process.
	//
	// Otherwise, depending on the proportion of infeasible portfolios generated,
	// possibly also abort the process.
	if (nbFeasibleGeneratedPortfolios === 0) {
		throw new Error('no feasible portfolio generated');
	}
	else {
		var nbInfeasibleGeneratedPortfolios = nbSubsets - nbFeasibleGeneratedPortfolios;
		
		if (nbInfeasibleGeneratedPortfolios >= 1 && nbInfeasibleGeneratedPortfolios/nbSubsets >= maxInfeasibleSubsetsRatio) {
			throw new Error('too many infeasible portfolios generated');
		}
	}
	
	// Note: the geometric center and the geometric median of m points in R^n
	// both lie within the convex hull of these m points.
	//
	// As a consequence:
	// - Long-only or long-short constraints
	// - Partial or full investment constraints
	// - Minimum/maximum weights constraints
	// imposed on the subset portfolios are automatically satisfied by the 
	// final portfolio.
	//
	// This is also the case for any convex and/or linear constraints (e.g.
	// 	maximum volatility constraint).
	var weights = null;
	if (subsetsAggregationMethod == 'average') {
		weights = geometricCenter_(generatedPortfoliosWeights);
	}
	else if (subsetsAggregationMethod == 'median') {
		weights = geometricMedian_(generatedPortfoliosWeights);
	}
	else  {
		throw new Error('internal error');
	}
	
	
	// ------
	
	// Return the computed portfolio weights
	return weights.toArray();
}



/**
* @function randomSubspaceMeanVarianceOptimizationWeights
*
* @summary Compute the weights of a portfolio using the random subspace optimization method 
* applied to the Markowitz mean-variance optimization algorithm.
*
* @description This function returns the weights w_1,...,w_n associated to a long-only fully or partially 
* invested portfolio of n assets, as computed by the random subspace optimization method applied 
* to the Markowitz mean-variance optimization algorithm of the first reference.
*
* C.f. the method randomSubspaceOptimizationWeights for more details.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see <a href="https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf">Breiman, L (2002), Manual On Setting Up, Using, And Understanding Random Forests V3.1</a>
*
* @param {Array.<number>} mu the returns of the n assets in the considered universe, array of n real numbers.
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt parameters for the random subspace optimization method, described in the method randomSubspaceOptimizationWeights, with 
* opt.sizeSubsets defaulting to the floored positive solution of the equation x^2 + 3x - SQRT(2*n*(n+3)) = 0.
* @param {number} opt.subsetsOpt parameters for the mean-variance optimization algorithm, 
* described in the method meanVarianceOptimizationWeights
* @return {Array.<number>} the weights corresponding to the computed portfolio, array of n real numbers.
*
* @example
* randomSubspaceMeanVarianceOptimizationWeights([0.1, 0.2, 0.15], [[1, 0.3, -0.2], [0.3, 1, 0.6], [-0.2, 0.6, 1]], 
*                                   			{ subsetsGenerationMethod: 'deterministic', 
*		  										  subsetsOpt: {
*													  constraints: {
*														volatility: Math.sqrt(0.10)
*													  }
*												  }
*												})
* // ~[0.09, 0.19, 0.12]
*/
self.randomSubspaceMeanVarianceOptimizationWeights = function(mu, sigma, opt) {	
	// Initialize the options structure
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetsOpt === undefined) {
		opt.subsetsOpt = {};
	}
	if (opt.subsetsOpt.constraints === undefined) {
		opt.subsetsOpt.constraints = {};
	}
	
	// ------
	
	
	// Initializations
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbColumns;


	// ------	
	
	
	// Decode options

	// The default size of the subsets to generate is obtained following 
	// the sixth reference: for a random forest based on m features, the default 
	// value of the number of sub features to use in each random tree is SQRT(m).
	//
	// Applied to the case of mean-variance optimization:
	// - The number m of features is nbAssets (the number of returns to estimate) +
	// nbAssets*(nbAssets+1)/2 (the number of variances/covariances to estimate), 
	// which is equal to nbAssets*(nbAssets+3)/2
	//
	// - The number of assets X corresponding to a number of features SQRT(m)
	// is the solution of the equation X*(X+3)/2 = SQRT(m),
	// i.e. the solution of the equation X^2 + 3X - 2*SQRT(m) = 0.
	if (opt.sizeSubsets === undefined) {
		// Define the coefficients of the second order polynomial x^2 + 3x - 2*SQRT(m)
		var a = 1;
		var b = 3;
		var c = -2*Math.sqrt(nbAssets*(nbAssets+3)/2);
				
		// Extract the strictly positive root x^* of the equation x^2 + 3x - 2*SQRT(m) = 0, using a stable numerical formula
		var b_p = b/2; // reduced discriminant
		var sign_b_p = 1;
		var disc = b_p*b_p - a*c; // > 0 because a,b are positive and c is negative
		var q = -(b_p + Math.sqrt(disc));
		var r2 = c/q; // positive because c and q are negative
		
		opt.sizeSubsets = Math.max(2, Math.floor(r2));
	}

	
	// ------
	
	
	// Core process
	
	// Define the options for the mean-variance optimization algorithm
	opt.subsetOptimizationFctOpt = opt.subsetsOpt;
	

	// Definition of the portfolio optimization algorithm to use on the subsets
	function subsetMeanVarianceOptimization(subsetAssetsIdx, subsetOptimizationFctOpt) {
		// Extract the returns of the selected assets
		var subsetMu = mu.submatrix(subsetAssetsIdx, [1]);

		// Extract the covariance matrix of the selected assets
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);
		
		// Return the weights of the mean-variance optimal portfolio of the selected assets
		//
		// Catches non reachable return/volatility constraint, as well as infeasible problem, which all
		// can be raised due to the subsetting of assets.
		try {
			return self.meanVarianceOptimizationWeights(subsetMu, subsetSigma, subsetOptimizationFctOpt);
		}
		catch (e) {
			if (e.message.includes('no matching efficient portfolio') ||
				e.message === 'infeasible problem detected: the restricted simplex is empty') {
				throw new Error('infeasible portfolio optimization problem');
			}
			else {
				throw(e);
			}
		}
		
	}

	// Return the computed portfolio weights using the generic random subspace optimization method
	return self.randomSubspaceOptimizationWeights(nbAssets, subsetMeanVarianceOptimization, opt);
}


/**
* @function randomSubspaceGlobalMinimumVarianceOptimizationWeights
*
* @summary Compute the weights of a portfolio using the random subspace optimization method 
* applied to the Markowitz minimum variance optimization algorithm.
*
* @description This function returns the weights w_1,...,w_n associated to a long-only fully or partially 
* invested portfolio of n assets, as computed by the random subspace optimization method applied 
* to the Markowitz minimum variance optimization algorithm of the first reference.
*
* C.f. the method randomSubspaceOptimizationWeights for more details.
*
* @see Harry M. Markowitz, Portfolio Selection, Efficient Diversification of Investments, Second edition, Blackwell Publishers Inc.
* @see <a href="https://blog.thinknewfound.com/2018/07/machine-learning-subset-resampling-and-portfolio-optimization/">Machine Learning, Subset Resampling, and Portfolio Optimization</a>
* @see <a href="https://aaai.org/ocs/index.php/AAAI/AAAI17/paper/view/14443">SHEN, W.; WANG, J.. Portfolio Selection via Subset Resampling. AAAI Conference on Artificial Intelligence, North America, feb. 2017.</a>
* @see <a href="https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf">Breiman, L (2002), Manual On Setting Up, Using, And Understanding Random Forests V3.1</a>
*
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt parameters for the random subspace optimization method, described in the method randomSubspaceOptimizationWeights, with 
* opt.sizeSubsets defaulting to the floored positive solution of the equation x^2 + x - SQRT(n*(n+1)/2) = 0 if opt.mu is not provided
* or to the floored positive solution of the equation x^2 + 3x - SQRT(n*(n+3)/2) = 0 if opt.mu is provided.
* @param {number} opt.mu the optional returns of the n assets in the considered universe, array of n real numbers; defaults to an array of zeroes.
* @param {number} opt.subsetsOpt optional parameters for the minimum-variance optimization algorithm used in the subsets
* @param {number} opt.subsetsOpt.optimizationMethodParams.epsGsmo the convergence tolerance of the GSMO algorithm used to solve the subsets minimum-variance optimization problem, 
* a strictly positive number; defaults to 1e-6.
* @param {number} opt.subsetsOpt.optimizationMethodParams.maxIterGsmo the maximum number of iterations of the GSMO algorithm used to solve the subsets minimum-variance optimization problem, 
* a strictly positive natural integer; defaults to 10000.
* @param {boolean} opt.subsetsOpt.constraints.fullInvestment parameter set to false in case the full investment constraint of the subsets portfolio must be replaced
* by a partial investment constraint; defaults to true.
*
* @return {Array.<number>} the weights corresponding to the computed portfolio, array of n real numbers.
*
* @example
* randomSubspaceGlobalMinimumVarianceOptimizationWeights([0.1, 0.2, 0.15], [[1, 0.3, -0.2], [0.3, 1, 0.6], [-0.2, 0.6, 1]], 
*                                   				     { subsetsGenerationMethod: 'deterministic' })
* // ~[0.09, 0.19, 0.12] // TODO
*/
self.randomSubspaceGlobalMinimumVarianceOptimizationWeights = function(sigma, opt) {	
	// Initialize the options structure
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetsOpt === undefined) {
		opt.subsetsOpt = {};
	}
	if (opt.subsetsOpt.optimizationMethodParams === undefined) {
		opt.subsetsOpt.optimizationMethodParams = {};
	}
	if (opt.subsetsOpt.constraints === undefined) {
		opt.subsetsOpt.constraints = {};
	}
	
	// ------
	
	
	// Initializations
	var mu = opt.mu; // the optional assets returns
	if (mu != undefined) {
		mu = new Matrix_(mu);
	}
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbColumns;


	// ------	
	
	
	// Decode options

	// The default size of the subsets to generate is obtained following 
	// the fourth reference: for a random forest based on m features, the default 
	// value of the number of sub features to use in each random tree is SQRT(m).
	//
	// Applied to the case of minimum variance optimization:
	// - If the assets returns are provided, the number m of features is nbAssets (the number of returns to estimate) +
	// nbAssets*(nbAssets+1)/2 (the number of variances/covariances to estimate), 
	// which is equal to nbAssets*(nbAssets+3)/2
	//
	//   In this case, the number of assets X corresponding to a number of features SQRT(m)
	// is the solution of the equation X*(X+3)/2 = SQRT(m),
	// i.e. the solution of the equation X^2 + 3X - 2*SQRT(m) = 0.
	//
	//
	// - Otherwise, the number m of features is nbAssets*(nbAssets+1)/2 (the number of variances/covariances to estimate)
	//
	//   In this case, the number of assets X corresponding to a number of features SQRT(m)
	// is the solution of the equation X*(X+1)/2 = SQRT(m),
	// i.e. the solution of the equation X^2 + X - 2*SQRT(m) = 0.
	if (opt.sizeSubsets === undefined) {
		// Define the coefficients of the second order polynomial x^2 + x - 2*SQRT(m) OR x^2 + 3x - 2*SQRT(m)
		var a = 1;
		var b = mu == undefined ? 1 : 3;
		var c = mu == undefined ? -2*Math.sqrt(nbAssets*(nbAssets+1)/2) : -2*Math.sqrt(nbAssets*(nbAssets+3)/2);
				
		// Extract the strictly positive root x^* of the equation x^2 + x - 2*SQRT(m) = 0 OR x^2 + 3x - 2*SQRT(m) = 0, using a stable numerical formula
		var b_p = b/2; // reduced discriminant
		var sign_b_p = 1;
		var disc = b_p*b_p - a*c; // > 0 because a,b are positive and c is negative
		var q = -(b_p + Math.sqrt(disc));
		var r2 = c/q; // positive because c and q are negative
		
		opt.sizeSubsets = Math.max(2, Math.floor(r2));
	}

	
	// ------
	
	
	// Core process
	
	// Define the options for the minimum variance optimization algorithm
	opt.subsetOptimizationFctOpt = opt.subsetsOpt;
	
	// Definition of the portfolio optimization algorithm to use on the subsets
	function subsetGlobalMinimumVarianceOptimization(subsetAssetsIdx, subsetOptimizationFctOpt) {
		// Extract the covariance matrix of the selected assets
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);

		// If provided, extract the returns of the selected assets
		var subsetMu;
		if (mu != undefined) {
			subsetMu = mu.submatrix(subsetAssetsIdx, [1]);
			subsetOptimizationFctOpt.mu = subsetMu;
		}
		
		// Return the weights of the global minimum variance portfolio of the selected assets
		//
		// Catches infeasible problem, which can be raised due to the subsetting of assets.
		try {
			return self.globalMinimumVarianceWeights(subsetSigma, subsetOptimizationFctOpt);
		}
		catch (e) {
			if (e.message === 'infeasible problem detected: the restricted simplex is empty') {
				throw new Error('infeasible portfolio optimization problem');
			}
			else {
				throw(e);
			}
		}
		
	}

	// Return the computed portfolio weights using the generic random subspace optimization method
	return self.randomSubspaceOptimizationWeights(nbAssets, subsetGlobalMinimumVarianceOptimization, opt);
}
