/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function randomCorrelationMatrix
*
* @summary Returns a random correlation matrix.
*
* @description This function computes a random n by n correlation matrix, 
* using the algorithm described in the reference
* 
* @see <a href="https://link.springer.com/article/10.1023/A:1022384216930">Philip I. Davies and Nicholas J. Higham, Numerically Stable Generation of Correlation Matrices and Their Factors,BIT Numerical Mathematics volume 40, pages 640–651 (2000)</a>
*
* @param {number} n the row/column length of the matrix to construct, natural integer greater than or equal to 1.
* @param {object} opt optional parameters for the random correlation matrix generation algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-14.
* @param {Array.<number>} opt.lambda the desired eigenvalues lambda_i, i=1..n of the generated correlation matrix, an array of n real numbers lambda_i,
* which must satisfy 0 <= lambda_i, i = 1..n and sum_i lambda_i = n; defaults to an array of random numbers belonging to [0,1] and summing to one.
* @param {number} opt.epsLambda tolerance for the condition that the provided eigenvalues must sum to one; defaults to 1e-8.
*
* @return {Matrix_} the computed matrix.
*
*/
self.randomCorrelationMatrix = function(n, opt) {
	// Generate the random correlation matrix
	return Matrix_.randomCorrelation(n, opt);
}



/**
* @function nearestCorrelationMatrix
*
* @summary Returns the nearest correlation matrix to a symmetric matrix.
*
* @description This function computes the nearest correlation matrix to a n by n symmetric matrix, 
* in Frobenius norm, using the algorithm described in the first references, with extensions described
* in the second reference.
* 
* @see <a href="https://www.maths.manchester.ac.uk/~higham/narep/narep369.pdf">Nicholas J Higham. Computing the nearest correlation matrix - a problem from finance. IMA Journal of Numerical Analysis, 22(3):329-343, 2002.</a>
* @see <a href="https://link.springer.com/article/10.1007/s11075-015-0078-3">Nicholas J Higham. Natasa Strabic. Anderson acceleration of the alternating projections method for computing the nearest correlation matrix. Numer Algor (2016) 72:1021–1042</a>
* @see <a href="https://link.springer.com/article/10.1007/s11075-015-0078-3">N. Higham, N. Strabic, Anderson acceleration of the alternating projections method for computing the nearest correlation matrix, Numer Algor (2016) 72:1021–1042</a>
*
* @param {Matrix_|Array.<Array.<number>>} mat a n by n matrix, or an array of n arrays of n real numbers.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-6, and in case it is set to 0, force convergence to full precision, which can be time consuming.
* @param {number} opt.maxIter the maximum number of iterations of the alternating projections algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100.
* @param {number} opt.minEigenvalue lower bound on the smallest eigenvalue(s) of the nearest correlation matrix, a positive real number belonging to [0,1]; defaults to 1e-12.
*
* Note: The usual minimum value of opt.minEigenvalue is 1e-8, which does not seem to impact the number of iterations of the algorithm (c.f. the second reference); values
* greater than 0.1 are more impacting. 
*
* @return {Matrix_} the computed matrix.
*
*/
self.nearestCorrelationMatrix = function(mat, opt) {
	// Internal function to project a symmetric matrix on the set U, c.f. formula 3.2
	// of the first reference.
	function unitDiagonalMatricesFrobeniusProjector(mat) {
		// Safety check
		if (!mat.isSymmetric(epsSymmetric)) {
			throw new Error('internal error: input matrix must be symmetric');
		}

		// Set 1's on the input matrix diagonal
		var p = mat.unitDiagonalize();
		
		// Return the computed projection
		return p;
	}

	// Internal function to project a symmetric matrix on the set S_delta, c.f. theorem 3.4
	// of the second reference.
	//
	// For delta = 0, S_0 corresponds to the set S of the first reference, c.f. formula 3.3
	// of the first reference.
	function pdMatricesFrobeniusProjector(mat, delta) {
		// Safety check
		if (!mat.isSymmetric(epsSymmetric)) {
			throw new Error('internal error: input matrix must be symmetric');
		}
	
		// Compute the spectral decomposition of mat
		var jacobi = Matrix_.eig(mat, {epsSymmetric: epsSymmetric});
		var q = jacobi[0];
		var lambdas = jacobi[1];
		
		// Truncate the eigenvalues to minimum delta
		var lambdas_p = lambdas.elemMap(function(i,j,val) { return Math.max(delta, val); });
		
		// Reconstruct the matrix
		var p = Matrix_.axty(1, Matrix_.elementwiseProduct(q, lambdas_p.transpose()), q);

		// Polish the non-diagonal elements, to ensure symmetry
		p = p.symmetrize();
		
		// Return the computed projection
		return p;
	}

		
	// ------

	// Convert mat to matrix format
	var mat = new Matrix_(mat);
	
	
	// ------
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var eps = opt.eps
	if (eps == undefined) {
		eps = 1e-6;
	}
	var maxIterations = opt.maxIter
	if (maxIterations == undefined) {
		maxIterations = 100;
	}
	var minEigenvalue = opt.minEigenvalue
	if (minEigenvalue == undefined) {
		minEigenvalue = 1e-12;
	}
	
	
	// ------
	
	// Checks
	if (!mat.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	
	
	// ------
	
	// Initializations
	var n = mat.nbRows;
	var x_k = Matrix_.copy(mat);
	var y_k = Matrix_.copy(mat);
	var delta_s_k = Matrix_.zeros(n, n);

	var diff_x_rel = Number.POSITIVE_INFINITY;
	var diff_y_rel = Number.POSITIVE_INFINITY;
	var diff_yx_rel = Number.POSITIVE_INFINITY;
	
	// Core loop corresponding to algorithm 3.3 of the first reference, 
	// guaranteed to converge.
	var iter = 0;
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Update the number of iterations
		++iter;
		
		// Compute the modified y vector
		var r = Matrix_.xmy(y_k, delta_s_k);
		
		// Compute the next x_k (project r on the set S or on the set S_delta)
		var x_kp = pdMatricesFrobeniusProjector(r, minEigenvalue);

		// Compute next Dykstra's correction
		var delta_s_kp = Matrix_.xmy(x_kp, r);
		
		// Compute the next y_k (project x on the set U)
		var y_kp = unitDiagonalMatricesFrobeniusProjector(x_kp);
		
		// Update the relative differences in norms, c.f. formula 4.1 of the first reference
		diff_x_rel = Matrix_.xmy(x_kp, x_k).vectorNorm('infinity') / x_k.vectorNorm('infinity');
		diff_y_rel = Matrix_.xmy(y_kp, y_k).vectorNorm('infinity') / y_k.vectorNorm('infinity');
		diff_yx_rel = Matrix_.xmy(y_kp, x_kp).vectorNorm('infinity') / y_kp.vectorNorm('infinity');
		
		// Prepare the next iteration
		x_k = x_kp;
		y_k = y_kp;
		delta_s_k = delta_s_kp;
	
		// Check convergence, c.f. formula 4.1 of the first reference
		//
		// To be noted that checking only the formula 4.1 can lead to non unit diagonal x_k,
		// so that in case eps is set to zero, check is done on y_k being a correlation matrix
		// instead, but this is a costly requirement (i.e., convergence to full precision).
		if (eps == 0.0) {
			if (y_k.isCorrelationMatrix()) {
				break;
			}
		}
		else {
			if (Math.max(diff_x_rel, diff_y_rel, diff_yx_rel) <= eps) {
				break;
			}
		}
	}
	
	 
	// Return the computed matrix, which is numerically unit diagonal, but
	// which might not be strictly unit diagonal in case eps <> 0.0.
	if (eps == 0.0) {
		return y_k;
	}
	else {
		return x_k;
	}
}

/**
* @function repairCorrelationMatrix
*
* @description This function transforms a real symmetric indefinite (or positive semi definite) correlation matrix 
* into a positive semi definite (or positive definite) correlation matrix.
*
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1969689">Riccardo Rebonato, Peter Jaeckel, The Most General Methodology to Create a Valid Correlation Matrix for Risk Management and Option Pricing Purposes</a>
* @see <a href="https://epubs.siam.org/doi/abs/10.1137/140996112">N. Higham, N.Strabic, V. Sego Restoring Definiteness via Shrinking, with an Application to Correlation Matrices with a Fixed Block, SIAM Rev., 58(2), 245–263</a>
* @see <a href="https://link.springer.com/article/10.1007/s11075-015-0078-3">N. Higham, N. Strabic, Anderson acceleration of the alternating projections method for computing the nearest correlation matrix, Numer Algor (2016) 72:1021–1042</a>
*
* @param {Matrix_|Array.<Array.<number>>} corrMat the possibly indefinite correlation matrix to repair (corr_ij),i,j=1..n, array of n arrays of n real numbers satisfying corrMat[i-1][j-1] = corr_ij, or a n by n Matrix_,
* with corrMat symmetric and with unit diagonal.
* @param {object} opt optional parameters for the algorithm.
* @param {string} opt.method the method to use to repair the correlation matrix, a string either equals to:
* - "nearest-correlation-matrix", in order to repair the correlation matrix by computing its nearest correlation matrix as described in the third reference
* - "spectral", in order to repair the correlation matrix using the spectral decomposition method described in the first reference
* - "linear-shrinkage", in order to repair the correlation matrix using the shrinkage method towards a target matrix as described in the second reference
; defaults to "nearest-correlation-matrix"
* @param {number} opt.minEigenvalue in case opt.method is either equal to "linear-shrinkage" or to "nearest-correlation-matrix", a lower bound on the smallest eigenvalue(s) of the repaired correlation matrix, a positive real number belonging to [0,1]; defaults to 0.
* @param {Matrix_|Array.<Array.<number>>} opt.shrinkageTarget in case opt.method is equal to "linear-shrinkage", the target correlation matrix towards which to shrink 
* the input correlation matrix, whose smallest eigenvalue(s) must satisfy the opt.minEigenvalue lower bound; default to the identity matrix.
* @param {number} opt.epsNcm in case opt.method is equal to "nearest-correlation-matrix", tolerance for the convergence of the associated algorithm, a strictly positive real number; defaults to 1e-6.
* @param {number} opt.epsSymmetric tolerance for the numerical symmetry of the input matrix corrMat, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.epsUnitDiagonal tolerance for the numerical unit diagonal of the input matrix corrMat, a strictly positive real number; defaults to 1e-12.
*
* @return {Matrix_} the repaired correlation matrix based on corrMat.
*
*/
self.repairCorrelationMatrix = function(corrMat, opt) {
	// ------
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var method = opt.method;
	if (method === undefined) {
		method = "nearest-correlation-matrix";
	}
	var minEigenvalue = opt.minEigenvalue
	if (minEigenvalue == undefined) {
		minEigenvalue = 0;
	}
	var epsNcm = opt.epsNcm
	if (epsNcm == undefined) {
		epsNcm = 1e-6;
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var epsUnitDiagonal = opt.epsUnitDiagonal
	if (epsUnitDiagonal == undefined) {
		epsUnitDiagonal = 1e-12;
	}
	var shrinkageTarget = opt.shrinkageTarget;
	
	
	// Decode the parameters
	if (method != "spectral" && method != "linear-shrinkage" && method != "nearest-correlation-matrix") {
		throw new Error('unsupported correlation matrix repair method');
	}
	

	// ------
	
	// Convert corrMat to matrix format
	var corrMat = new Matrix_(corrMat);
	
	// Checks
	if (!corrMat.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	if (!corrMat.isUnitDiagonal(epsUnitDiagonal)) {
		throw new Error('input matrix must be unit diagonal');
	}
	// C.f. the third reference for this test
	if (minEigenvalue < 0 || minEigenvalue > 1) {
		throw new Error('lower bound on the smallest eigenvalue(s) of the correlation matrix must belong to interval [0,1]');
	}
	
	
	// ------
	
	// Initializations
	var n = corrMat.nbRows;
	
	
	// Compute the eigenvalues of the input correlation matrix
	// to ensure that it is really broken.
	//
	// If not, simply return it.
	var jacobi = Matrix_.eig(corrMat, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
	if (jacobi[1].data[n-1] >= minEigenvalue) {
		return corrMat;
	}
	
	
	// Several algorithms are implement below to repair a non positive definite correlation matrix:
	// - Repair the eigenvalues of the correlation matrix (spectral method), c.f. the first reference
	// - Shrink the correlation matrix towards a target correlation matrix (linear shrinkage method), c.f. the second reference
	// - Compute its nearest correlation matrix, c.f. the third reference
	var c;
	if (method == "spectral") {
		// Compute the eigenvectors and eigenvalues of the input broken correlation matrix
		var s = jacobi[0];
		var lambda = jacobi[1];
		
		// Construct the matrix SQRT(lambda'), by truncating all the eigenvalues strictly inferior to 0,
		// c.f. formula 9 of the first reference, and then taking their square root.
		var lambdap_sqrt = lambda.elemMap(function(i,j,val) { return Math.sqrt(Math.max(0, val)); });
		
		// Construct the matrix B', c.f. formula 11 of the first reference.
		var bp = Matrix_.elementwiseProduct(s, lambdap_sqrt.transpose());
		
		// Construct the matrix B, c.f. formula 12 of the first reference, by normalizing 
		// the row vectors of B' to unit length.
		var t_sqrt = Matrix_.fill(n, 1, function(i,j) { return 1 / bp.vectorNorm('two', 'row', i); });
		var b = Matrix_.elementwiseProduct(bp, t_sqrt);
		
		// Construct the matrix C, which is the repaired correlation matrix, c.f. formula 13
		// of the first reference.
		c = Matrix_.axty(1, b, b);
		
		// Per construction, the matrix c is strictly symmetric, but its unit diagonal
		// is only numeric, so, fix this.
		c = c.unitDiagonalize();
	}
	else if (method == "linear-shrinkage") {
		// Define the shrinkage target matrix as the identity matrix, c.f. the second reference,
		// or convert the provided shrinkage target matrix to matrix format.
		//
		// In case a shrinkage target matrix is provided, ensure that it is a correlation matrix 
		// and that its smallest eigenvalue is greater than or equal to the desired lower bound 
		// on the smallest eigenvalue(s) of the repaired correlation matrix.
		if (shrinkageTarget == undefined) {
			shrinkageTarget = Matrix_.identity(n);
		}
		else {
			shrinkageTarget = new Matrix_(shrinkageTarget);
			
			//
			if (!shrinkageTarget.isSymmetric(epsSymmetric)) {
				throw new Error('shrinkage target matrix must be symmetric');
			}
			if (!shrinkageTarget.isUnitDiagonal(epsUnitDiagonal)) {
				throw new Error('shrinkage target matrix must be unit diagonal');
			}
			
			//
			var jacobi = Matrix_.eig(shrinkageTarget, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
			if (jacobi[1].data[n-1] < minEigenvalue) {
				throw new Error('smallest eigenvalue of the shrinkage target matrix strictly lower than the desired lower bound on the smallest eigenvalue(s) of the correlation matrix');
			}
		}
		
		// C.f. section 3.1 of the second reference, the optimal alpha belonging to [0,1] is searched
		// using the bisection method, so that the minimum eigenvalue of the matrix S(alpha) =
		// alpha*Target + (1-alpha)*corrMat is greater than or equal to minEigenvalue.
		//
		// This is possible because the function f = lambda_min (S(alpha)) - minEigenvalue is continuous 
		// and concave on [0,1].
		var m0 = corrMat;
		var m1 = shrinkageTarget;
		
		var alpha_star_interval = bisection_(function (alpha) { 
												// Compute the matrix S(alpha)
												var s_alpha = Matrix_.axpby(alpha, m1, 1-alpha, m0);
												
												// Compute the minimal eigenvalue of the matrix S(alpha)
												var jacobi = Matrix_.eig(s_alpha, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
												var lambdas = jacobi[1];
												var lambda_min = lambdas.data[n-1];

												// Return the value of the function f
												return lambda_min - minEigenvalue; 
											}, 
											0, 1, {outputInterval: true});
		
		// The optimal value of alpha, alpha^*, is the upper bound of the interval
		// found by bisection, to ensure that S(alpha^*) has its minimum eigenvalue
		// greater than or equal to minEigenvalue, c.f. the second reference.
		var alpha_star = alpha_star_interval[1];
		
		// Compute the matrix S(alpha^*)
		c = Matrix_.axpby(alpha_star, m1, 1-alpha_star, m0);
		
		// Per construction, the matrix is only numerically symmetric and unit diagonal.
		//
		// Ensure the matrix C is symmetric with unit diagonal
		c = c.symmetrize().unitDiagonalize();
 	}
	else if (method == "nearest-correlation-matrix") {
		// Compute the nearest correlation matrix 
		c = self.nearestCorrelationMatrix(corrMat, {maxIter: -1, eps: epsNcm, minEigenvalue: minEigenvalue});
		
		// Per construction, the matrix C is symmetric with unit diagonal
	}
	else {
		throw new Error('internal error: unsupported repair method');
	}

	 
	// Return the computed matrix
	return c;
}



/**
* @function perturbedCorrelationMatrix
*
* @description This function computes a randomly perturbed version of an input correlation matrix.
*
* @see <a href="https://onlinelibrary.wiley.com/doi/10.1002/9781118467381.ch6">Roza Galeeva, Jiri Hoogland, and Alexander Eydeland, Measuring Correlation Risk for Energy Derivatives</a>
* @see <a href="https://jpm.pm-research.com/content/19/2/6">Chopra, Vijay K; Ziemba, William T; The effect of errors in means, variances, and covariances on optimal portfolio choice; Journal of Portfolio Management; Winter 1993; 19, 2</a>
* @see <a href="https://projecteuclid.org/euclid.aoas/1380804814">Hardin, Johanna; Garcia, Stephan Ramon; Golan, David. A method for generating realistic correlation matrices. Ann. Appl. Stat. 7 (2013), no. 3, 1733--1762.</a>
*
* @param {Matrix_} an n by 1 Matrix representing the correlation matrix to perturb.
* @param {object} opt optional parameters for the perturbation algorithm.
* @param {string} opt.method the method to use to perturb the correlation matrix, a string either equals to:
* - "additive-noise", in order to perturb the correlation matrix off-diagonal entries additively using independent normal random variables, c.f. the third reference
* - "multiplicative-noise", in order to perturb the correlation matrix off-diagonal entries multiplicatively using independent normal random variables, c.f. the second reference
* - "spectral-noise", in order to perturb the correlation matrix eigenvalues using independent normal random variables, c.f. the first reference
; defaults to "additive-noise"
* @param {number} opt.sigma in case opt.method is equal to "spectral-noise" or "multiplicative-noise", the common standard deviation of the independent normal random variables,
* a positive real number; defaults to 0.01
* @param {number} opt.noiseLevelMax in case opt.method is equal to "additive-noise", the maximum noise level, a positive real number; defaults to 0.01
* @param {number} opt.noiseSpaceDimension in case opt.method is equal to "additive-noise", the dimension of the noise space, a positive integer; 
* defaults to 3 to perturb all the elements uniformly, c.f. the third reference.
* @param {number} opt.epsNcm in case pt.method is equal to "multiplicative-noise" or to "additive-noise", the generated correlation matrix might not be semidefinite positive, and
* requires the computation of the nearest correlation matrix; if so, tolerance for the convergence of the associated algorithm, a strictly positive real number; defaults to 1e-6.
* @param {number} opt.maxIterNcm in case pt.method is equal to "multiplicative-noise" or to "additive-noise", the generated correlation matrix might not be semidefinite positive, and
* requires the computation of the nearest correlation matrix; if so, maximum number of iterations for the convergence of the associated algorithm, a strictly positive real number; defaults to 100.

* @param {number} opt.epsCorrMat tolerance for the numerical symmetry and unit diagonal of the input matrix corrMat, a strictly positive real number; defaults to 1e-12.
*
* @return {Matrix_} a Matrix object representing the perturbed correlation matrix.
*
*/
self.perturbedCorrelationMatrix = function(corrMat, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "additive-noise";
	}
	var epsCorrMat = opt.epsCorrMat
	if (epsCorrMat == undefined) {
		epsCorrMat = 1e-12;
	}
	var epsNcm = opt.epsNcm
	if (epsNcm == undefined) {
		epsNcm = 1e-6;
	}
	var maxIterNcm = opt.maxIterNcm
	if (maxIterNcm == undefined) {
		maxIterNcm = 100;
	}
	
	var sigma = opt.sigma;
	if (sigma == undefined) {
		sigma = 0.05;
	}
	
	var noiseLevelMax = opt.noiseLevelMax;
	if (noiseLevelMax == undefined) {
		noiseLevelMax = 0.01;
	}
	var noiseSpaceDimension = opt.noiseSpaceDimension;
	if (noiseSpaceDimension == undefined) {
		noiseSpaceDimension = 3;
	}
	
	// Decode the parameters
	var method = opt.method;
	if (method != "spectral-noise" && method != "multiplicative-noise" && method != "additive-noise") {
			throw new Error('unsupported perturbation method');
	}
	

	// ------
	
	// Convert corrMat to matrix format
	var corrMat = new Matrix_(corrMat);
	
	// Checks
	if (!corrMat.isCorrelationMatrix(epsCorrMat)) {
		throw new Error('input matrix must be symmetric, unit diagonal and positive semidefinite');
	}

	
	// ------
		
	// Initializations
	var n = corrMat.nbRows;
	
	
	// Several algorithms are implement below to perturb a correlation matrix:
	// - Perturb the eigenvalues of the correlation matrix, c.f. section 6.3.4 of the first reference
	// - Perturb the non-diagonal elements of the correlation matrix with multiplicative, c.f. the second reference
	// - Perturb the non-diagonal elements of the correlation matrix with additive noise, c.f. the third reference
	var c;
	if (method == "spectral-noise"){		
		// Compute the eigenvectors and eigenvalues of the input correlation matrix
		var jacobi = Matrix_.eig(corrMat, {epsSymmetric: epsCorrMat, sortedEigenvalues: true});
		var s = jacobi[0];
		var lambdas = jacobi[1];
		
		// Compute the perturbed eigenvalues
		var lambdas_p = lambdas.elemMap(function(i,j,val) { return val * Math.exp(sigma * normrnd_(0, 1)); });
			
		// Construct the perturbed correlation matrix, not yet renormalized
		var cp = Matrix_.axty(1, Matrix_.elementwiseProduct(s, lambdas_p.transpose()), s);
		
		// Renormalize the above perturbed correlation matrix
		c = cp.elemMap(function(i,j,val) { return val / Math.sqrt(cp.data[(i-1) * cp.nbColumns + (i-1)] * cp.data[(j-1) * cp.nbColumns + (j-1)]); });
		
		// Ensure the matrix C is symmetric with unit diagonal
		c = c.symmetrize().unitDiagonalize();
	}
	else if (method == "multiplicative-noise") {
		// Perturb the correlation matrix non-diagonal elements,
		// which must stay bounded in [-1, 1].
		c = Matrix_.fillSymmetric(n, function(i,j) { 
										if (i == j) { 
											return 1; 
										} 
										else { 
											return Math.max(-1, Math.min(1, corrMat.data[(i-1) * corrMat.nbColumns + (j-1)] * (1 + sigma*normrnd_(0,1))));
										}
									});
		
		// Ensure the perturbed correlation matrix is still semidefinite positive
		if (!c.isCorrelationMatrix()) {
			c = self.nearestCorrelationMatrix(c, {maxIter: maxIterNcm, eps: epsNcm});
		}
	}
	else if (method == "additive-noise") {
		// C.f. algorithm 4 of the third reference
		
		// Compute the n unit vectors from R^noiseSpaceDimension using random uniform vectors 
		// on the R^noiseSpaceDimension unit sphere, c.f. paragraph 3.1 of the third reference
		var hypersphereRandomSampler = new hypersphereRandomSampler_(noiseSpaceDimension, true);
		var us = new Array(n);
		for (var i = 0; i < n; ++i) {
			us[i] = new Matrix_(hypersphereRandomSampler.sample());
		}
		
		// Perturb the correlation matrix non-diagonal elements using the 
		// matrix noiseLevelMax * (U^t * U - I), taking care that they 
		// must stay bounded in [-1, 1].
		c = Matrix_.fillSymmetric(n, function(i,j) { 
										if (i == j) { 
											return 1; 
										} 
										else { 
											var corrMat_ij = corrMat.data[(i-1) * corrMat.nbColumns + (j-1)];
											var noise_ij = noiseLevelMax * Matrix_.vectorDotProduct(us[i-1], us[j-1]);
											
											return Math.max(-1, Math.min(1, corrMat_ij + noise_ij));
										}
									});
									
		// Ensure the perturbed correlation matrix is still semidefinite positive
		if (!c.isCorrelationMatrix()) {
			c = self.nearestCorrelationMatrix(c, {maxIter: maxIterNcm, eps: epsNcm});
		}									
	}
	else {
		throw new Error('internal error: unsupported perturbation method');
	}

	 
	// Return it
	return c;
}
