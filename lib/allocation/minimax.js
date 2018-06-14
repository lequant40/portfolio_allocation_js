/**
 * @file Functions related to minimax weights portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function minimaxWeights
*
* @summary Compute the weights of a minimax portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only 
* minimax portfolio of n assets.
*
* Optionally, the following constraint can be added:
* - Partial investment contraint, replacing the full investment contraint
*
* A minimax portfolio has the property that it maximizes the minimum possible return over the period on which it is computed,
* c.f. the first reference.
*
* This portfolio might not be unique.
*
* @see <a href="http://www.jstor.org/stable/2634472">Young, M. (1998). A Minimax Portfolio Selection Rule with Linear Programming Solution. Management Science, 44(5), 673-683.</a>
* @see <a href="https://link.springer.com/article/10.1007/s11135-005-1054-0">Yuanyao Ding, (2006), Portfolio Selection under Maximum Minimum Criterion, Quality & Quantity: International Journal of Methodology, 40, (3), 457-468</a>
*
* @param {Array.<Array.<number>>} assetsReturns an array of n arrays of T real numbers representing the returns of n assets over T periods of time.
* @param {object} opt optional parameters for the algorithm.
* @param {boolean} opt.constraints.fullInvestment parameter set to false in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to true.
* @param {boolean} opt.outputMinimumPortfolioReturn a boolean indicating whether the computed minimum portfolio return should be provided in output
* (if set to true) or not (if set to false); defaults to false.
* @return {Array.<Object>|Array.<number>} if opt.outputMinimumPortfolioReturn is set to true, an array arr of two elements:
* - arr[0], the weights corresponding to a minimax portfolio, array of n real numbers
* - arr[1], the computed minimum portfolio return, a real number
*
* otherwise, the weights corresponding to a minimax portfolio, an array of n real numbers.
*
* @example
* minimaxWeights([[0.01, -0.02, 0.01], [-0.05, 0.03, 0.01]]);
* // ~[ 0.727, 0.273 ]
*/
self.minimaxWeights = function (assetsReturns, opt) {
	// TODO: Checks, if enabled

	// Decode options
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	var fullInvestmentContraint = true;
	if (opt.constraints.fullInvestment !== undefined) {
		fullInvestmentContraint = opt.constraints.fullInvestment;
	}
	var outputMinimumPortfolioReturn = false || opt.outputMinimumPortfolioReturn; 
	
	// Initializations
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;
	
	
	// ----
	
	// The minimax portfolio is the solution to a linear programming problem,
	// choosen to be the one of the formula 21 of the second reference instead of
	// the one of the section 1.1 of the first reference.
	//
	// In other zords, no minimum return is imposed on the portfolio, 
	// so that the linear program is always feasible.

		// Build the objective function (c.f. formula 1a of the first reference):
		// - Maximize the minimum portfolio return
	var c = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 0 : -1; }); // c = [0,...,0,-1]
	
		// Build the equality constraints (c.f. formula 1d of the first reference
		// for the equivalent inequality constraint):
		// - Full investment (optional)
	var Ae = null;
	var be = null;
	if (fullInvestmentContraint) {
		Ae = Matrix_.fill(1, nbAssets + 1, function(i,j) { return j <= nbAssets ? 1 : 0; }); // Ae = [1,...,1,0]
		be = Matrix_.ones(1,1); // be = [1]
	}
	
		// Build the inequality constraints (c.f. formula 1b of the first reference):
		// - Portfolio return greater than or equal to the minimum portfolio return, for each period
		// - Partial investment (optional)
	var Ai = Matrix_.fill(nbPeriods + (fullInvestmentContraint ? 0 : 1), nbAssets + 1, 
						  function(i,j) { 
								if (i <= nbPeriods) { return j <= nbAssets ? -assetsReturns[j-1][i-1] : 1; }
								else if (i == nbPeriods + 1) { return j <= nbAssets ? 1 : 0; }
						  }); // Ai = [[-ret11,...,-retN1,1], ..., [-ret1T,...,-retNT,1] (optional: , [1,...,1,0])]
	var bi = Matrix_.fill(nbPeriods + (fullInvestmentContraint ? 0 : 1), 1, 
						  function(i,j) { 
							  if (i <= nbPeriods) { return 0; }
							  else if (i == nbPeriods + 1) { return 1; }
						  }); // bi = [0, ..., 0 (optional: , 1)]
	
		// Build the bound constraints (c.f. formula 1e of the first reference
		// as well as the definition of M_p):
		// - No short sales
		// - Absence of leverage
		// - Unbounded minimum portfolio return
	var lb = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 0 : -Infinity; }); // lb = [0,...,0, -Infinity]
	var ub = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 1 : Infinity; });  // ub = [1,...,1, Infinity]
	
		// Solve the constructed linear program, which is:
		// - Bounded: the portfolio weights belong to the unit simplex
		//            the minimum portfolio return is bounded below by the minimum portfolio return over all the periods
		// - Feasible: the portfolio with the minimum return over all the periods is a solution to the linear program
		//
		// Note: given the assumptions above, the convergence of the primal-dual hybrid gradient algorithm is guaranteed.
	var lpSolution = lpsolvePDHG_(Ae, be, Ai, bi, c, lb, ub, {maxIter: -1});
	
	
	// ----
	
	// Extract the computed portfolio weights, plus the minimum portfolio return.
	var sol = lpSolution[0];
	var weights = sol.toArray(function(i, j, val) { return i != nbAssets + 1; });
	var minPortfolioReturn = sol.data[nbAssets];

	// Return the computed weights, and potentially return the minimum
	// portfolio return.
	if (outputMinimumPortfolioReturn === true) {
		return [weights, minPortfolioReturn];
	}
	else {
		return weights;
	}
}
