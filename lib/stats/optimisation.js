/**
 * @file Misc. optimisation functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.lpsolveAffineScaling_ = lpsolveAffineScaling_;
/* End Wrapper private methods - Unit tests usage only */
 
 
 // c,x n cost vector; x vector unknows
// b m; RHS of constraints
// A m n, m constraints matrix;  m constraints, n variables  
// Because the optimal solutions to linear programs reside on the boundary of the feasible region, interior point algorithms cannot move to the exact optimal solution
// * @see <a href="https://link.springer.com/article/10.1007/BF01840454">Robert J. Vanderbei; Meketon, Marc; Freedman, Barry (1986). "A Modification of Karmarkar's Linear Programming Algorithm" (PDF). Algorithmica. 1: 395–407.</a>
// * @see <a href="https://link.springer.com/article/10.1007/BF02592024">Barnes, E.R. A variation on Karmarkar’s algorithm for solving linear programming problems. Mathematical Programming (1986) 36: 174.</a>
// * @see <a href="https://link.springer.com/article/10.1007/BF01582276">Vanderbei, R.J. Affine-scaling for linear programs with free variables. Mathematical Programming (1989) 43: 31.</a>
// * @see <a href="https://link.springer.com/article/10.1007%2FBF02592206">Monteiro, R.D.C.; Tsuchiya, T.; Wang, Y.; A simplified global convergence proof of the affine scaling algorithm. Mathematical Programming (1996) 75: 77.</a>
function lpsolveAffineScaling_(A, b, c, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var epsRelOptimality = opt.eps || 1e-12;
	var epsAbsInfeasibility = 1e-14; // Note: must be agressive with this one, to avoid false positives
	var epsAbsFeasibility = 1e-10; // Note: beware with these, as problem can easily be marked as infeasible if too agressive
	var epsRelFeasibility = 1e-08; 
	var t = opt.t || 0.666666666666;
	var maxIterations = opt.maxIter || 10000;
	
	// ------
	
	// Misc. checks
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(c instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	if (A.nbRows > A.nbColumns) {
		throw new Error('first input matrix has more rows than columns: ' + '(' + A.nbColumns + ') v.s. ' + '(' + A.nbRows + ')');
	}
	if (A.nbRows !== b.nbRows) {
		throw new Error('first and second inputs number of rows do not match: ' + A.nbRows + '-' + b.nbRows);
	}
	if (A.nbColumns !== c.nbRows) {
		throw new Error('first input number of columns and third input number of rows do not match: ' + A.nbColumns + '-' + c.nbRows);
	}
	
	// ------

	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	var fctVal = Number.MAX_VALUE;
	var b_inf_norm = b.vectorNorm('infinity');
	var x_k = Matrix_.ones(A.nbColumns, 1); // the solution vector; initial value such that x_0 > 0
	var tx_k = x_k.transpose(); // the transpose of the solution vector
	var ad_k = Matrix_.zeros(m, n); // the rescaled matrix
	var p_k = Matrix_.zeros(m, 1); // the degree of primal feasibility
	var r_k = Matrix_.zeros(n, 1); // the reduced costs vector
	var d_k_r_k = Matrix_.zeros(n, 1); // ??
	var z_k = Matrix_.zeros(n, 1); // the affine scaling direction/step direction
	
	// ------
		
	// Main loop of the algorithm.
	//
	// This algorithm is proved to converge under the following assumptions for 0 < t <= 2/3, 
	// c.f. the fourth reference assumptions 1-4 and theorem 3.3:
	// - Rank A = m (if A is rank-deficient, either there are redundant constraints or the problem is infeasible)
	// - The linear function <x/c> is not constant on the feasible region (TODO: isn't this catched by the algorithm ?)
	// - The feasible region has a non-empty interior
	// - The problem has an optimal solution
	var iter = 0;
	while (true) {
		// Update the number of iterations
		++iter;

		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Compute the objective function value (not in the references)
		fctVal = Matrix_.vectorDotProduct(c, x_k);

		// The scaling matrix D_k will never be computed, and
		// operations using it will be transformed into element wise
		// matrix products.

		// Preliminary computation of the rescaled matrix A*D_k
		tx_k = x_k.transpose(tx_k);
		ad_k = Matrix_.elementwiseProduct(A, tx_k, ad_k);
		
		// ------
		
		// Phase 1 - Primal interior feasibility checks/repair on x_k, c.f. the algorithm
		// of the section 3 of the third reference:
		// - Primal feasibility of x_k: A*x_k == b <=> ||A*x_k - b||_inf == 0
		// - Interiority of x_k: x_k > 0
		//
		// To be noted that the primal affine scaling algorithm theoretically 
		// guarantees that the sequence (x_k)_k=1.. is always primal feasible if x_0 is,
		// c.f. for instance proposition 2.3 of the fourth reference.
		//
		// Unfortunately, due to round off errors, this might not be the case
		// in practice, so that this condition must be tested on each iteration
		// for the proper convergence of the algorithm.
		//
		// To also be noted that there should not be such an issue with the interiority,
		// thanks to the update formula for x_k.
		//
		// The algorithm below is described in the section 3 of the third reference,
		// with sign differences due to the sign of p_k above, which is the opposite of 
		// the sign of rho_k.
		p_k = Matrix_.axpby(1, Matrix_.product(A, x_k), -1, b, p_k); // degree of primal feasibility
		if (p_k.vectorNorm('infinity') > epsAbsFeasibility || p_k.vectorNorm('infinity') > epsRelFeasibility * (1 + b_inf_norm)) {
			// Compute w_k = (A*D_k^2*A^t)^-1*p_k
		    //
		    // This is equivalent to solving the linear system (A*D_k)*(A*D_k)^t w_k = p_k (1)
		    // ---
		    // Compute the lhs of the linear system (1)
		    var lhs = Matrix_.axty(1, ad_k, ad_k);

		    // Compute the rhs of the linear system (1)
		    var rhs = p_k;
		        
            // Solve the linear system (1)
            var w_k = Matrix_.linsolveKaczmarz(lhs, rhs, {maxIter: -1});			
            // ---
          
		    // Compute r_k = A^t*w_k
		    r_k = Matrix_.atxy(1, A, w_k, r_k);
		    
		    // Compute ??
			d_k_r_k = Matrix_.elementwiseProduct(r_k, x_k, d_k_r_k);
		 	
			// Compute the affine scaling direction
			z_k = Matrix_.elementwiseProduct(d_k_r_k, x_k, z_k);
		 	
			// Stopping/infeasibility conditions, c.f. section 3 of the third reference.
			var gamma = d_k_r_k.vectorNorm('infinity'); // degree of complementary slackness
			var delta = -r_k.min(); // degree of dual feasibility
			var beta = Matrix_.vectorDotProduct(p_k, w_k);
			if (gamma + delta * x_k.max() < epsAbsInfeasibility / n) { // the problem is infeasible
				return [null, null, 1];
			}
			else if (gamma < t) { // the problem is feasible, the next iteration will be a phase 2 iteration
				// The update of the current solution computes a primal feasible and strictly interior point.
				x_k = Matrix_.axpby(1, x_k, -1, z_k, x_k);
			}
			else { // progress is made, the next iteration will be a phase 1 iteration
				// Prepare the next iteration: update of the current solution x_k+1 = x_k - t/gamma * z_k
				x_k = Matrix_.axpby(1, x_k, -t/gamma, z_k, x_k);
			}
				
			// In all cases, continue to the next iteration
			continue;
		}
		else if (!x_k.isPositive()) {
			throw new Error('interiority lost for the sequence of iterates at iteration: ' + iter + ', with value: ' + x_k.toString());
		}
		

		// ------
		
		// Phase 2 - "Real" affine scaling algorithm
		//
		// The algorithm below is described in formula 12 of the fourth reference.
		
		// Compute the dual variables vector w_k = (A*D_k^2*A^t)^-1 * A*D_k^2*c
		//
		// This is equivalent to solving the linear system (A*D_k*(A*D_k)^t) w_k = A*D_k*D_k*c (2)
		// ---		
		// Compute the lhs of the linear system (2)
		var lhs = Matrix_.axty(1, ad_k, ad_k);

		// Compute the rhs of the linear system (2)
		var rhs = Matrix_.product(ad_k, Matrix_.elementwiseProduct(c, x_k)); 
			
		// Solve the linear system (1)
		var w_k = Matrix_.linsolveKaczmarz(lhs, rhs, {maxIter: -1});
		// ---
		
		// Compute the reduced costs vector r_k = c - A^t*w_k
		r_k = Matrix_.axpby(1, c, -1, Matrix_.atxy(1, A, w_k), r_k);
		
		// Compute ??
		d_k_r_k = Matrix_.elementwiseProduct(r_k, x_k, d_k_r_k);
		
		// Optimality condition, c.f. proposition 6 of the first reference
		// and step B.6 of the third reference.
		var gamma = d_k_r_k.vectorNorm('infinity'); // degree of complementary slackness
		var delta = -r_k.min(); // degree of dual feasibility
		if (gamma + delta * x_k.max() < epsRelOptimality / n * (1 + Math.abs(fctVal))) { // an approximate optimal interior point has been found
			return [x_k, fctVal, 0];
		}

		// Other conditions, c.f. proposition 1 of the first reference:
		// - If D_k*r_k == 0, every feasible point is optimal (TODO: constant objective value catch ???)
		// - If D_k*r_k <= 0 and D_k*r_k <> 0, the problem is unbounded (TODO: test, OR  is it  < 0 ??)
		if (gamma <= Number.EPSILON) { // every feasible point is optimal
			return [x_k, fctVal, 2];
		}
		if (d_k_r_k.isNonPositive()) { // the problem is unbounded
			return [null, Number.NEGATIVE_INFINITY, 3];
		}
		
		// Prepare the next iteration:
		// - Compute the affine scaling direction/step direction
		// - Update the current solution x_k+1 = x_k - t * D_k^2*r_k / ||D_k * r_k||_inf
		z_k = Matrix_.elementwiseProduct(d_k_r_k, x_k, z_k);
		x_k = Matrix_.axpby(1, x_k, -t/gamma, z_k, x_k);
	}

	
	
	
	
	
	
	
	// * See Vander book for phase 2
	function phaseOne(A, b, c, eps, t, maxIterations) {
		// Initializations
		var m = A.nbRows;
		var n = A.nbColumns;
		var x_k = Matrix_.ones(n, 1); // the solution vector, with initial value such that x_0 > 0
		var d_k = Matrix_.zeros(n, n); // the (diagonal) scaling matrix
		var p_k = Matrix_.zeros(m, 1); // the degree of primal feasibility
		var r_k = Matrix_.zeros(n, 1); // the reduced costs vector
		var d_k_r_k = Matrix_.zeros(n, 1); // ??
		var z_k = Matrix_.zeros(n, 1); // the affine scaling direction/step direction
		
		// Main loop of the algorithm, implemented following the description of the second reference.
		//
		// This algorithm is proved to converge under the following assumptions for 0 < t <= 2/3, 
		// c.f. the fourth reference assumptions 1-4 and theorem 3.3:
		// - Rank A = m (if A is rank-deficient, either there are redundant constraints or the problem is infeasible)
		// - The linear function <x/c> is not constant on the feasible region (TODO: isn't this catched by the algorithm ?)
		// - The feasible region has a non-empty interior
		// - The problem has an optimal solution
		var iter = 0;
		var converged = false;
		var convergenceCode = -1;
		while (!converged) {
			// Update the number of iterations
			++iter;

			// Check the number of iterations
			if (iter > maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			
			// Compute the degree of primal feasibility
			// Note, the sign of p_k is the opposite of the sign of rho_k in the algorithm
			// of the section 3 of the third reference, which changes the computation of z_k below
			p_k = Matrix_.axpby(1, Matrix_.product(A, x_k), -1, b, p_k); 
			
			// Compute the scaling matrix D_k
		    d_k = Matrix_.diagonal(x_k, d_k);

		    // Compute w_k = (A*D_k^2*A^t)^-1*p_k
		    //
		    // This is equivalent to solving the linear system (A*D_k^2*A^t) w_k = p_k (1)
		    // ---
		    // Preliminary computation of the rescaled matrix A*D_k
		    var ad_k = Matrix_.product(A, d_k);
		        
		    // Compute the lhs of the linear system (1)
		    var lhs = Matrix_.product(ad_k, ad_k.transpose());

		    // Compute the rhs of the linear system (1)
		    var rhs = p_k;
		        
            // Solve the linear system (1)
            var w_k = Matrix_.linsolveKaczmarz(lhs, rhs);			
            // ---
            
		    // Compute r_k = A^t*w_k
		    var r_k = Matrix_.product(A.transpose(), w_k, r_k);
		    
		    // Compute ??
			d_k_r_k = Matrix_.product(d_k, r_k, d_k_r_k);
			
			// Compute the affine scaling direction
			z_k = Matrix_.product(d_k, d_k_r_k, z_k);
			
			// Stopping/infeasibility conditions, c.f. section 3 of the third reference.
			var gamma = d_k_r_k.vectorNorm('infinity'); // degree of complementary slackness
			var delta = -r_k.min(); // degree of dual feasibility
			var beta = Matrix_.vectorDotProduct(p_k, w_k);
			if (gamma + delta * x_k.max() < eps / n) {
				converged = true;
				convergenceCode = 1;
				x_k = null;
		        break;
			}
			else if (gamma < t) {
				converged = true;
				convergenceCode = 0;
				x_k = Matrix_.axpby(1, x_k, -1, z_k, x_k);
				break;
			}

		    // Prepare the next iteration:
			// - Update the current solution x_k+1 = x_k - t/gamma * z_k
			x_k = Matrix_.axpby(1, x_k, -t/gamma, z_k, x_k);
		}

        // Return the computed feasible vector
        return [x_k, convergenceCode];
	}
	

	// TODO: more epsilons
	function phaseTwo(A, b, c, x_0, eps, t, maxIterations) {
		// TODO: Check that x_0 is feasible: A*x_0 = b and x_0 > 0

		// Initializations
		var m = A.nbRows;
		var n = A.nbColumns;
		var x_k = Matrix_.copy(x_0); // the solution vector
		var fctVal = Number.MAX_VALUE;
		var d_k = Matrix_.zeros(n, n); // the (diagonal) scaling matrix
		//var p_k = Matrix_.zeros(m, 1); // the degree of primal feasibility
		var r_k = Matrix_.zeros(n, 1); // the reduced costs vector
		var d_k_r_k = Matrix_.zeros(n, 1); // ??
		var z_k = Matrix_.zeros(n, 1); // the affine scaling direction/step direction
		
		// Main loop of the algorithm, implemented following the description of the second reference.
		//
		// This algorithm is proved to converge under the following assumptions for 0 < t <= 2/3, 
		// c.f. the fourth reference assumptions 1-4 and theorem 3.3:
		// - Rank A = m (if A is rank-deficient, either there are redundant constraints or the problem is infeasible)
		// - The linear function <x/c> is not constant on the feasible region (TODO: isn't this catched by the algorithm ?)
		// - The feasible region has a non-empty interior
		// - The problem has an optimal solution
		var iter = 0;
		var converged = false;
		var convergenceCode = -1;
		while (!converged) {
			// Initializations of the loop variables
			var primalFeasible = false;

			// Update the number of iterations
			++iter;

			// Check the number of iterations
			if (iter > maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			
		    // Primal feasibility conditions (not in the references):
			// - A*x_k == b <=> ||A*x_k - b||_inf == 0
			// - x_k >= 0
			//
			// The primal affine scaling algorithm theoretically guarantees that the sequence x_k is primal feasible,
			// c.f. for instance proposition 2.3 of the fourth reference.
			//
			// Nevertheless, due to round off errors, this might not be the case, so that this condition must be tested
			// for the convergence of the algorithm.
			// p_k = Matrix_.axpby(1, Matrix_.product(A, x_k), -1, b, p_k);
			//if (p_k.vectorNorm('infinity') <= eps && x_k.isNonNegative()) {
				primalFeasible = true;
			//}
			
			// Compute the objective function value (not in the references)
			fctVal = Matrix_.vectorDotProduct(c, x_k);
		
			// Compute the scaling matrix D_k
		    d_k = Matrix_.diagonal(x_k, d_k);

		    // Compute the dual variables vector w_k = (A*D_k^2*A^t)^-1 * A*D_k^2*c
		    //
		    // This is equivalent to solving the linear system (A*D_k^2*A^t) w_k = A*D_k^2*c,
		    // equivalent to (A*D_k*(A*D_k)^t) w_k = A*D_k*D_k*c (1)
		    // ---
		    // Preliminary computation of the rescaled matrix A*D_k
		    var ad_k = Matrix_.product(A, d_k);
		        
		    // Compute the lhs of the linear system (1)
		    var lhs = Matrix_.product(ad_k, ad_k.transpose());

		    // Compute the rhs of the linear system (1)
		    var rhs = Matrix_.product(ad_k, Matrix_.product(d_k, c));
		        
            // Solve the linear system (1)
            var w_k = Matrix_.linsolveKaczmarz(lhs, rhs);
            // ---
            
		    // Compute the reduced costs vector r_k = c - A^t*w_k
		    var r_k = Matrix_.axpby(1, c, -1, Matrix_.product(A.transpose(), w_k), r_k);
		    
		    // Compute ??
			d_k_r_k = Matrix_.product(d_k, r_k, d_k_r_k);

			// Optimality condition, c.f. proposition 6 of the first reference
			// and step B.6 of the third reference.
			var gamma = d_k_r_k.vectorNorm('infinity'); // degree of complementary slackness
			var delta = -r_k.min(); // degree of dual feasibility
			if (primalFeasible && gamma + delta * x_k.max() < eps / n * (1 + Math.abs(fctVal))) {
				converged = true;
				convergenceCode = 0;
		        break;
			}

			// Other conditions, c.f. proposition 1 of the first reference:
			// - If D_k*r_k == 0, every feasible point is optimal (TODO: constant objective value catch ???)
			// - If D_k*r_k <= 0 and D_k*r_k <> 0, the problem is unbounded (TODO: test, OR  is it  < 0 ??)
		    if (primalFeasible && gamma <= Number.EPSILON) {
				converged = true;
				convergenceCode = 1;
				break;
			}
			else if (primalFeasible && d_k_r_k.isNonPositive()) {
				converged = true;
				convergenceCode = 2;
				fctVal = Number.NEGATIVE_INFINITY;
				break;
			}
		    
		    // Prepare the next iteration:
			// - Compute the affine scaling direction/step direction
			// - Update the current solution x_k+1 = x_k - t * D_k^2*r_k / ||D_k * r_k||_inf
			z_k = Matrix_.product(d_k, d_k_r_k, z_k);
			x_k = Matrix_.axpby(1, x_k, -t/gamma, z_k, x_k);
		}

        // Return the computed variables vector, plus the computed minimal objective function value
        return [x_k, fctVal, convergenceCode];
	}

	
}


