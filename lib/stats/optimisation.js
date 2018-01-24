/**
 * @file Misc. optimisation functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.lpsolveChambollePock_ = lpsolveChambollePock_;
self.nonNegativeOrthantEuclideanProjection_ = nonNegativeOrthantEuclideanProjection_;
/* End Wrapper private methods - Unit tests usage only */
 
 
function nonNegativeOrthantEuclideanProjection_(x, out) {
	// Misc. checks
	if (!(x instanceof Matrix_)) {
		throw new Error('input must be a matrix');
	}
	if (!x.isVector()) {
		throw new Error('input must be a vector');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(x.nbRows, 1, out);
	
	// Computation of proj(x_i)_[0, +infinity[, i = 1..x.nbRows
	for (var i = 1; i <= x.nbRows; ++i) {
		obj.data[(i-1)*x.nbColumns] = Math.max(0, x.data[(i-1)*x.nbColumns]);
	}
	
	// Return the computed matrix
	return obj;
}

 
 function lpsolveChambollePock_(Ae, be, Ai, bi, c, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var optimalityEps = opt.optimalityEps || 1e-08;
	var maxIterations = opt.maxIter || 100000;
	
	var eqContraints = false;
	var ineqContraints = false;
	if (Ae !== null && be !== null) {
		eqContraints = true;
	}
	if (Ai !== null && bi !== null) {
		ineqContraints = true;
	}

	
	// ------

	// Misc. checks
	if (!(c instanceof Matrix_)) {
		throw new Error('fifth input must be a matrix');
	}
	if (c.nbColumns !== 1) {
		throw new Error('fifth input is not a vector: ' + c.nbColumns + '-' + c.nbRows);
	}
	
	if (eqContraints) {
		if (!(Ae instanceof Matrix_)) {
			throw new Error('first input must be a matrix');
		}
		if (!(be instanceof Matrix_)) {
			throw new Error('second input must be a matrix');
		}
		if (Ae.nbRows !== be.nbRows) {
			throw new Error('first and second inputs number of rows do not match: ' + Ae.nbRows + '-' + be.nbRows);
		}
		if (Ae.nbColumns !== c.nbRows) {
			throw new Error('first input number of columns and fifth input number of rows do not match: ' + Ae.nbColumns + '-' + c.nbRows);
		}
		if (be.nbColumns !== 1) {
			throw new Error('second input is not a vector: ' + be.nbColumns + '-' + be.nbRows);
		}
	}
	
	if (ineqContraints) {
		if (!(Ai instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (!(bi instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (Ai.nbRows !== bi.nbRows) {
			throw new Error('third and fourth inputs number of rows do not match: ' + Ai.nbRows + '-' + bi.nbRows);
		}
		if (Ai.nbColumns !== c.nbRows) {
			throw new Error('third input number of columns and fifth input number of rows do not match: ' + Ai.nbColumns + '-' + c.nbRows);
		}
		if (bi.nbColumns !== 1) {
			throw new Error('fourth input is not a vector: ' + bi.nbColumns + '-' + bi.nbRows);
		}
	}

	
	// ------
	
	// Initializations
	// Constraints
	var me = 0;
	var ye_k = null;
	var ye_kp = null;
	var res_ye_kp_ye_k = null;
	var pe_k = null;
	var be_inf_norm = 0;
	if (eqContraints) {
	    me = Ae.nbRows; // the number of equality constaints
    	ye_k = Matrix_.zeros(me, 1); // the dual iterate ye_k
    	ye_kp = Matrix_.zeros(me, 1); // the dual iterate ye_k+1
    	res_ye_kp_ye_k = Matrix_.zeros(me, 1); // the residual ye_kp - ye_k
    	pe_k = Matrix_.zeros(me, 1); // the degree of primal feasibility on equality constraints
		be_inf_norm = be.vectorNorm('infinity');
	}
	var mi = 0;
	var yi_k = null;
	var yi_kp = null;
	var res_yi_kp_yi_k = null;
	var pi_k = null;
	var bi_inf_norm = 0;
	if (ineqContraints) {
	    mi = Ai.nbRows; // the number of inequality (<=) constaints
    	yi_k = Matrix_.zeros(mi, 1); // the dual iterate yi_k
    	yi_kp = Matrix_.zeros(mi, 1); // the dual iterate yi_k+1
    	res_yi_kp_yi_k = Matrix_.zeros(mi, 1); // the residual yi_kp - yi_k
    	pi_k = Matrix_.zeros(mi, 1); // the degree of primal feasibility on inequality constraints
		bi_inf_norm = bi.vectorNorm('infinity');
	}
	var m = me + mi; // the total number of constraints
	
    // Variables
    var n = c.nbRows;
	var x_k = Matrix_.ones(n, 1); // the primal iterate x_k
	var x_kp = Matrix_.ones(n, 1); // the primal iterate x_k+1
	var z_k = Matrix_.ones(n, 1); // the relaxed iterate z_k = 2*x_k+1 - x_k
	var res_x_kp_x_k = Matrix_.zeros(n, 1); // the residual x_kp - x_k

    // Misc.
	var fctVal = Number.MAX_VALUE; // the value of the objective function

	// Computation of the diagonal matrices T and S = [Se Si]^t with alpha = 1
	// and mu = 0.95, nu = 0.95 (so that mu * nu < 1), c.f. formula 10 and
	// remark 3 of the reference.
	var mu = 0.95;
	var nu = 0.95;
	var T = Matrix_.fill(n, 1, 
				function(i,j) { 
						var columnNorm = 0;
						if (eqContraints) {
						    columnNorm += Ae.vectorNorm('one', 'column', i); 
						}
						if (ineqContraints) {
						    columnNorm += Ai.vectorNorm('one', 'column', i);
						}
						return mu * 1/columnNorm; // Ae and Ai do not contain any null column;
				});
	var Se = null;
	if (eqContraints) {
	    Se = Matrix_.fill(me, 1, 
	                      function(i,j) { 
						    return nu * 1/Ae.vectorNorm('one', 'row', i); // Ae does not contain any null row;
				          });
	}
	var Si = null;
	if (ineqContraints) {
	    Si = Matrix_.fill(mi, 1, 
	                      function(i,j) { 
			        	    return nu * 1/Ai.vectorNorm('one', 'row', i); // Ai does not contain any null row;
				          });
	}

	
	// ------
	
	// Main loop of the algorithm
	// The convergence is guaranteed by the theorem 1 of the reference
	// in case both the primal and the dual problems admit an optimum.
	var iter = 0;
	while (true) {
		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Update the number of iterations
		++iter;

		// Primal update: 
		// - x_k+1 = proj(x_k - T*(A^t*y_k + c))_[0, +infinity[, with A = [Ae Ai]^t and y_k = [ye yi]^t
		var tmp = Matrix_.copy(c);
		if (eqContraints) {
		    tmp = Matrix_.axpby(1, Matrix_.atxy(1, Ae, ye_k), 1, tmp);
		}
		if (ineqContraints) {
		    tmp = Matrix_.axpby(1, Matrix_.atxy(1, Ai, yi_k), 1, tmp);
		}
		var val = Matrix_.axpby(1, x_k, -1, Matrix_.elementwiseProduct(tmp, T));
		x_kp = nonNegativeOrthantEuclideanProjection_(val, x_kp);

		// Relaxed iterate update:
		// - z_k = 2*x_k+1 - x_k
		z_k = Matrix_.axpby(2, x_kp, -1, x_k, z_k);
		
		// Dual update:
		// - ye_k+1 = ye_k + Se*(Ae*z_k - be) (equality constraints)
		// - yi_k+1 = proj(yi_k + Si*(Ai*z_k - bi))_[0, +infinity[ (inequality constraints)
		if (eqContraints) {
            ye_kp = Matrix_.axpby(1, ye_k, 1, Matrix_.elementwiseProduct(Matrix_.axpby(1, Matrix_.product(Ae, z_k), -1, be), Se));
		}
		if (ineqContraints) {
		    var val_2 = Matrix_.axpby(1, yi_k, 1, Matrix_.elementwiseProduct(Matrix_.axpby(1, Matrix_.product(Ai, z_k), -1, bi), Si));
			yi_kp = nonNegativeOrthantEuclideanProjection_(val_2, yi_kp);
		}

		// Compute the new objective function value (not in the reference)
		var newFctVal = Matrix_.vectorDotProduct(c, x_kp);
		
		// Convergence conditions for (x_k, y_k = [ye yi]^t) to be a saddle point of the min-max problem:
		// - Convergence of the primal iterates (relative) x_k: ||x_k+1 - x_k||_inf <= eps * ||x_k+1||_inf
		// - Convergence of the dual iterates (relative) y_k: ||y_k+1 - y_k||_inf <= eps * ||y_k+1||_inf
		// - Convergence of the objective function values (relative): |<c/x_k> - <c/x_k+1>| <= eps * |<c/x_k+1>|
		//
		// Then, as the problem is a linear programming problem, it can be additionally checked:
		// - Primal feasibility (absolute) Ae*x_k = be (equality constraints) and Ai*x_k <= bi (inequality constraints): 
		//  ||Ae*x_k - be||_inf <= eps and ||(Ai*x_k - bi)^+||_inf <= eps
		// - Primal feasibility x_k >= 0: guaranteed per construction of x_k
		res_x_kp_x_k = Matrix_.axpby(1, x_kp, -1, x_k, res_x_kp_x_k);
		var res_x_kp_x_k_inf_norm = res_x_kp_x_k.vectorNorm('infinity');
		var x_kp_inf_norm = x_kp.vectorNorm('infinity');
		
		var res_y_kp_y_k_inf_norm = 0;
		var y_kp_inf_norm = 0;
		if (eqContraints) {
		    res_ye_kp_ye_k = Matrix_.axpby(1, ye_kp, -1, ye_k, res_ye_kp_ye_k);
		    res_y_kp_y_k_inf_norm = Math.max(res_ye_kp_ye_k.vectorNorm('infinity'), 
		                                     res_y_kp_y_k_inf_norm);
											 
			y_kp_inf_norm = Math.max(ye_kp.vectorNorm('infinity'), y_kp_inf_norm);
		}
		if (ineqContraints) {
		    res_yi_kp_yi_k = Matrix_.axpby(1, yi_kp, -1, yi_k, res_yi_kp_yi_k);
		    res_y_kp_y_k_inf_norm = Math.max(res_yi_kp_yi_k.vectorNorm('infinity'), 
		                                     res_y_kp_y_k_inf_norm);
			
			y_kp_inf_norm = Math.max(yi_kp.vectorNorm('infinity'), y_kp_inf_norm);
		}
		
		var res_f_kp_f_k_inf_norm = Math.abs(newFctVal - fctVal);
        var f_kp_inf_norm = Math.abs(newFctVal);
        
		if (res_x_kp_x_k_inf_norm <= optimalityEps * x_kp_inf_norm  && 
		    res_y_kp_y_k_inf_norm <= optimalityEps * y_kp_inf_norm &&
		    res_f_kp_f_k_inf_norm <= optimalityEps * f_kp_inf_norm) {
			// Compute the degree of primal feasibility 
			var p_k_eq = true;
			var p_k_ineq = true;
			
			// Compute the degree of primal feasibility on equality constraints
			// pe_k = Ae*x_k+1 - be
			if (eqContraints) {
			    pe_k = Matrix_.axpby(1, Matrix_.product(Ae, x_kp), -1, be, pe_k);
			    p_k_eq = pe_k.vectorNorm('infinity') <= optimalityEps * be_inf_norm;
			}

			// Compute the degree of primal feasibility on inequality constraints
			// pi_k = (Ai*x_k+1 - bi)^+
			if (ineqContraints) {
			    pi_k = Matrix_.axpby(1, Matrix_.product(Ai, x_kp), -1, bi, pi_k);
				for (var i = 1; i <= mi; ++i) {
					pi_k.setValueAt(i, 1, Math.max(0, pi_k.getValueAt(i, 1)));
				}
				p_k_ineq = pi_k.vectorNorm('infinity') <= optimalityEps * bi_inf_norm;
			}
			
			// Convergence condition
			if (p_k_eq && p_k_ineq) {
				break;
			}
		}
		
		// Prepare the next iteration:
		// - x_k = x_k+1
		// - y_k = y_k+1 <=> ye_k = ye_k+1, yi_k = yi_k+1
		x_k = Matrix_.copy(x_kp, x_k);
		if (eqContraints) {
		    ye_k = Matrix_.copy(ye_kp, ye_k);
		}
		if (ineqContraints) {
		    yi_k = Matrix_.copy(yi_kp, yi_k);
		}
		fctVal = newFctVal;
	}
	
	// Return the computed primal iterate and the associated objective
	// function value.
	return [x_kp, newFctVal];		
 }
 
 
 

