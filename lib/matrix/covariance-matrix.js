/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function covarianceMatrixFromCorrelationMatrix
*
* @summary Returns the covariance matrix associated to a correlation matrix and to a
* variances vector or to a standard deviations vector.
*
* @description Returns the covariance matrix associated to a correlation matrix and to a
* variances vector or to a standard deviations vector.
*
* @see <a href="https://en.wikipedia.org/wiki/Covariance_matrix#Relation_to_the_correlation_matrix">Covariance matrix</a>
*
* @param {Matrix_|Array.<Array.<number>>} corrMat the correlation matrix (corr_ij),i,j=1..n, array of n arrays of n real numbers satisfying corrMat[i-1][j-1] = corr_ij, or a n by n Matrix_,
* with corrMat symmetric, with unit diagonal and with off-diagonal elements belonging to [-1,1].
* @param {Matrix_|<Array.<number>} diagonalVec the variance vector (if opt.diagonalVectorType is equal to "variances") or the standard deviations vector 
* (if opt.diagonalVectorType is equal to "standard-deviations") of the n assets in the considered universe, 
* an n by 1 matrix (i.e., vector) or an array of n real numbers.
* @param {string} opt.diagonalVectorType the type of diagonal vector provided in diagonalVec, a string either equals to:
* - "variances", in order to specify that diagonalVec is a vector of variances
* - "standard-deviations", in order to specify that diagonalVec is a vector of standard deviations
; defaults to "variances"
* @param {number} opt.epsSymmetric tolerance for the numerical symmetry of the input matrix corrMat, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.epsUnitDiagonal tolerance for the numerical unit diagonal of the input matrix corrMat, a strictly positive real number; defaults to 1e-12.
*
* @return {Matrix_} a Matrix object representing the covariance matrix.
*/
self.covarianceMatrixFromCorrelationMatrix = function(corrMat, diagonalVec, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.diagonalVectorType === undefined) {
		opt.diagonalVectorType = "variances";
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var epsUnitDiagonal = opt.epsUnitDiagonal
	if (epsUnitDiagonal == undefined) {
		epsUnitDiagonal = 1e-12;
	}
	
	
	// Decode the parameters
	var diagonalVectorType = opt.diagonalVectorType;
	if (diagonalVectorType != "variances" && diagonalVectorType != "standard-deviations") {
		throw new Error('unsupported diagonal vector type');
	}
	
	
	// Convert corrMat and diagonalVec to matrix format
	var corrMat = new Matrix_(corrMat);
	var diagonalVec = new Matrix_(diagonalVec);
		
	// Compatibility checks
	if (!corrMat.isSymmetric(epsSymmetric)) {
		throw new Error('input correlation matrix not symmetric');
	}
	if (corrMat.nbRows != diagonalVec.nbRows) {
		throw new Error('input correlation matrix dimensions not compatible with input variances vector dimension');
	}
	if ( Matrix_.xmy(corrMat.diagonal(), Matrix_.ones(corrMat.nbRows, 1)).vectorNorm('infinity') > epsUnitDiagonal ) {
		throw new Error('input correlation matrix not unit diagonal');
	}
	
	
	// Compute the covariance matrix associated to the input correlation matrix and
	// to the input variances, using the formula Cov = Diag(stddev) * Corr * Diag(stddev).
	var n = corrMat.nbRows;
	
	var stddev;
	if (diagonalVectorType == "variances") {
		stddev = diagonalVec.elemMap(function(i,j,val) { return Math.sqrt(val); });
	}
	else if (diagonalVectorType == "standard-deviations") {
		stddev = diagonalVec;
	}
	else {
		throw new Error('internal error: unsupported diagonal vector type');
	}
	var c = Matrix_.elementwiseProduct(Matrix_.elementwiseProduct(corrMat, stddev), stddev.transpose());
	
	
	// Polish the covariance matrix to ensure it is symmetric
	c = c.symmetrize();
	
	
	// Return it
	return c;
}

/**
* @function covarianceMatrix
*
* @summary Returns the covariance matrix of series of values.
*
* @description This function computes the covariance matrix of series of values.
*
* @see <a href="https://en.wikipedia.org/wiki/Covariance_matrix">Covariance matrix</a>
* @see <a href="https://en.wikipedia.org/wiki/Sample_mean_and_covariance">Sample mean and covariance</a>
* @see <a href="https://link.springer.com/article/10.1023/A:1022384216930">O. Ledoit, M. Wolf, A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices, Journal of Multivariate Analysis, Volume 88, Issue 2, February 2004, pages 365-411.</a>
* @see <a href="https://jpm.pm-research.com/content/30/4/110">O. Ledoit, M. Wolf, Honey, I Shrunk the Sample Covariance Matrix, The Journal of Portfolio Management Summer 2004, 30 (4) 110-119</a>
* @see <a href="https://pubmed.ncbi.nlm.nih.gov/16646851/">J. Schafer, K. Strimmer, A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics, Statistical Applications in Genetics and Molecular Biology, Volume 4, Issue 1, 2005</a>
* @see <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3400062">G. De Nard, Oops! I Shrunk the Sample Covariance Matrix Again: Blockbuster Meets Shrinkage, Journal of Financial Econometrics, 2020</a>
*
* @param {Array.<Array.<number>>} arr an array of n arrays of m real numbers, with n and m natural integers
* greater than or equal to 1, with n representing the number of series (or features) and m representing the number of observations per series.
* @param {object} opt optional parameters for the covariance matrix computation.
* @param {string} opt.method the method to use to compute the covariance matrix, a string either equals to:
* - "covariance", in order to compute the covariance matrix, c.f. the first reference
* - "sample-covariance", in order to compute the sample covariance matrix, c.f. the second reference
* - "linear-shrinkage", in order to compute the asymptotically optimal convex linear combination of the covariance matrix with a shrinkage target matrix, c.f. 
* the fourth reference
; defaults to "covariance"
* @param {string} opt.shrinkageTarget in case opt.method is equal to "linear-shrinkage", the shrinkage target matrix to use, a string either equals to:
* - "constant-variance-null-correlation", in order to use a multiple of the identity matrix, c.f. the third reference
* - "constant-variance-correlation", in order to use a multiple of the identity matrix plus a multiple of the matrix with 1s everywhere except on the diagonal with 0s, c.f. the sixth reference
* - "null-correlation", in order to use a covariance matrix made of different variances and a null correlation coefficient, c.f. the fifth reference (target D), which is computed using 
* the generic formula in appendix B of fourth reference applied to the used target shrinkage matrix
* - "constant-correlation", in order to use a covariance matrix made of different variances and a constant correlation coefficient, c.f. the fourth reference
*
* @return {Matrix_} a Matrix object representing the covariance matrix of the input series of values.
*
* @example
* covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]]);
* // == Matrix_([[0.00036,  -0.00053], [-0.00053, 0.00107]])
*
* covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {method: "sample-covariance"});
* // == Matrix_([[0.00053, -0.0008], [-0.0008, 0.0016]])
*/
self.covarianceMatrix = function(arr, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "covariance";
	}
	
	
	// Decode the parameters
	var covarianceMethod = opt.method;
	if (covarianceMethod != "covariance" &&
	    covarianceMethod != "sample-covariance" &&
		covarianceMethod != "linear-shrinkage") {
			throw new Error('unsupported covariance matrix computation method');
	}

	var shrinkageTarget = opt.shrinkageTarget;
	if ( covarianceMethod == "linear-shrinkage" &&
	    ["constant-variance-null-correlation", 
		 "constant-variance-correlation", 
		 "null-correlation", 
		 "constant-correlation"].indexOf(shrinkageTarget) == -1 ) {
		 throw new Error('unsupported covariance matrix shrinkage target');
	}
	
	
	//
	var nbSeries = arr.length;
	var nbObservations = arr[0].length;
	
	
	//
	var obj;
	
	// In case the covariance matrix computation method is either "covariance" or "sample-covariance",
	// proceed with a similar computation formula.
	if (covarianceMethod == "covariance" || covarianceMethod == "sample-covariance") {
		// Define the covariance function to use
		var covarianceFunction = covariance_;
		if (covarianceMethod == "sample-covariance") {
			covarianceFunction = sampleCovariance_;
		}

		// Construct the covariance matrix
		obj = allocateMatrix_(nbSeries, nbSeries);
		for (var i = 0; i < obj.nbRows; ++i) {
			// Copy from upper triangular part
			for (var j = 0; j < i; ++j) {
				obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
			}
			
			// Computation part
			for (var j = i; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = covarianceFunction(arr[i], arr[j]);
			}
		}
	}
	else if (covarianceMethod == "linear-shrinkage") {
		// All "Ledoitt-Wolf" like linear shrinkage estimators share
		// a similar structure.
		
		// De-mean the series
		var means = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]) });
		var x = Matrix_.fill(nbSeries, nbObservations, function(i,j) { return arr[i-1][j-1] - means.data[i-1] });
		
		// Compute the sample covariance matrix as defined in the third reference
		//
		// Note: contrary to the usual sample covariance matrix definition, there is no division per n-1, but by n !
		var covMat = Matrix_.axty(1/nbObservations, x, x).toCovarianceMatrix();
		
		// Extract the standard deviations and the correlation matrix
		var stdvarVec = covMat.getStandardDeviations();	
		var corrMat = covMat.getCorrelationMatrix();
		
		// Compute the prior, c.f. appendix A of the fourth reference
		// for the general formula.
		var prior;
		var rBar;
		if (shrinkageTarget == "constant-variance-null-correlation") {
			// C.f. also lemma 3.2 of the third reference.
			var mu = covMat.trace() / nbSeries;
			prior = Matrix_.fill(nbSeries, nbSeries, function(i,j) { if (i == j) { return mu; } else { return 0; }  });
		}
		else if (shrinkageTarget == "constant-variance-correlation") {
			// C.f. also formulas 2.15a and 2.15b of the sixth reference
			var phi = covMat.trace() / nbSeries;
			var nu = 0;
			for (var i = 0; i < nbSeries - 1; ++i) {
				for (var j = i+1; j < nbSeries; ++j) {
					nu += covMat.data[i * covMat.nbColumns + j];
				}
			}
			nu *= 2/(nbSeries*(nbSeries-1)); 
			prior = Matrix_.fill(nbSeries, nbSeries, function(i,j) { if (i == j) { return phi; } else { return nu; }  });
		}	
		else if (shrinkageTarget == "null-correlation") {
			// C.f. also "Target D" in table 2 of the fifth reference
			prior = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									 if (i == j) { 
										 return stdvarVec.data[i-1]*stdvarVec.data[i-1]; 
									 } 
									 else { 
										 return 0;
									 } 
								 });
		}
		else if (shrinkageTarget == "constant-correlation") {
			rBar = 0;
			for (var i = 0; i < nbSeries - 1; ++i) {
				for (var j = i+1; j < nbSeries; ++j) {
					rBar += corrMat.data[i * corrMat.nbColumns + j];
				}
			}
			rBar *= 2/(nbSeries*(nbSeries-1)); 
			prior = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									 if (i == j) { 
										 return stdvarVec.data[i-1]*stdvarVec.data[i-1]; 
									 } 
									 else { 
										return rBar*stdvarVec.data[i-1]*stdvarVec.data[j-1] ;
									} 
								 });
		}

		// Compute pi hat, c.f. appendix B of the fourth reference
		// for the general formula.
		//
		// To be noted that in the case of constant variance/null covariance,
		// pi hat corresponds (up to a constant) to (b_n)^2,
		// c.f. lemma 3.4 of the third reference.
		var piMat = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									//
									var pi_ij = 0;
									for (var k = 0; k < nbObservations; ++k) {
										 pi_ij += Math.pow((arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)], 2);
									}
									pi_ij /= nbObservations;
								
									//
									return pi_ij;
								 });
		var pi = piMat.sum();

		// Compute rho hat, c.f. appendix B of the fourth reference
		// for the general formula.
		var rho;
		if (shrinkageTarget == "constant-variance-null-correlation") {
			// Not needed from the third reference, c.f. also code cov1para.m from the authors
			// of the third reference.
			rho = 0;
		}
		else if (shrinkageTarget == "constant-variance-correlation") {
			// C.f. remark 2.4 of the sixth reference, this coefficient can
			// be neglected for all practical purposes.
			//
			// C.f. appendix A.1 of the sixth reference for the computation
			// formula, if ever needed.
			rho = 0;
		}	
		else if (shrinkageTarget == "null-correlation") {
			// In the specific case of null correlation, all the terms
			// AsyCov are null, so that only the terms AsyVar remain.
			rho = piMat.trace();
		}
		else if (shrinkageTarget == "constant-correlation") {
			rho = 0;
			for (var i = 1; i <= nbSeries; ++i) {
				for (var j = 1; j <= nbSeries; ++j) {
					// The sum defining rho skips i == j
					if (i == j) {
						continue;
					}
					
					// Compute theta_ii__ij and theta_jj__ij
					var theta_ii__ij = 0;
					for (var k = 0; k < nbObservations; ++k) {				
						 theta_ii__ij += ( Math.pow(arr[i-1][k] - means.data[i-1], 2) - covMat.data[(i-1) * covMat.nbColumns + (i-1)] )
										  *
										  ( (arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)] );
					}
					theta_ii__ij /= nbObservations;
					
					var theta_jj__ij = 0;
					for (var k = 0; k < nbObservations; ++k) {
						 theta_jj__ij += ( Math.pow(arr[j-1][k] - means.data[j-1], 2) - covMat.data[(j-1) * covMat.nbColumns + (j-1)] )
										  *
										  ( (arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)] );
					}
					theta_jj__ij /= nbObservations;
			
					// Update the running sum for rho
					rho += (stdvarVec.data[j-1]/stdvarVec.data[i-1]) * theta_ii__ij  
						   +
						   (stdvarVec.data[i-1]/stdvarVec.data[j-1]) * theta_jj__ij;
				}
			}
			rho *= rBar/2;
			rho += piMat.trace();
		}
		
		// Compute gamma hat, described appendix B of the fourth reference.
		//
		// To be noted that in the case of constant variance/null covariance,
		// gamma hat corresponds (up to a constant) to (d_n)^2,
		// c.f. lemma 3.3 of the third reference. 
		var gamma = Math.pow(Matrix_.xmy(covMat, prior).matrixNorm('frobenius'), 2);

		// Compute the optimal shrinkage factor, c.f.  appendix B of the fourth reference.
		var kappa = (pi - rho)/gamma;
		var shrinkage = Math.max(0, Math.min(1, kappa/nbObservations));
		
		// Compute the optimally shrinked covariance matrix,
		// c.f. formula 2 of the fourth reference.
		obj = Matrix_.axpby(shrinkage, prior, 1-shrinkage, covMat);
	}
	else {
		throw new Error('internal error: unsupported covariance matrix computation method');
	}

	// Add covariance matrix methods
	addCovarianceMatrixMethods_(obj);
	
	// Return it
	return obj;
}




/**
* @function addCovarianceMatrixMethods_
*
* @summary Add methods related to a covariance matrix to a Matrix object.
*
* @description This function adds methods related to a covariance matrix to a Matrix object.
*
* @param {Matrix_} a, a matrix.
* @return {void}
*
* @example
* addCovarianceMatrixMethods_(Matrix_([[1,0.1], [0.1,1]]));
* // 
*/
function addCovarianceMatrixMethods_(matrix) {
	var methods = {
    	/**
    	* @function getCorrelationMatrix
    	*
    	* @summary Returns the correlation matrix associated to a covariance matrix.
    	*
    	* @description This function computes a correlation matrix (c_ij),i=1..n,j=1..n from the original 
    	* matrix (a_ij),i=1..n,j=1..n, with coefficients satisfying c_ij = a_ij/(SQRT(a_ii * a_jj)).
    	*
    	* @memberof Matrix_
		* @param {Matrix_} out an optional n by n matrix.
    	* @return {Matrix_} a n by n matrix containing the correlation matrix associated to the covariance matrix,
		* either stored in the matrix out or in a new matrix.
    	*
    	* @example
    	* Matrix_([[1, 1, 8.1], [1, 16, 18], [8.1, 18, 81]]).getCorrelationMatrix();
    	* // Matrix_([[1, 0.25, 0.9], [0.25, 1, 0.5], [0.9, 0.5, 1]])
    	*/
		'getCorrelationMatrix': function(out) { 
			// Result matrix allocation
			var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
			
			// Computation of the correlation matrix
			for (var i = 0; i < obj.nbRows; ++i) {
				// Standard deviation of a_ii
				var stdDevI = Math.sqrt(this.data[i * this.nbColumns + i]);
				
				// Copy from upper triangular part
				for (var j = 0; j < i; ++j) {
					obj.data[i * obj.nbColumns + j] = obj.data[j * this.nbColumns + i];
				}
				
				// Computation part
				for (var j = i; j < obj.nbColumns; ++j) {
					// Standard deviation of a_jj
					var stdDevJ = Math.sqrt(this.data[j * this.nbColumns + j]);
				
					obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j] / ( stdDevI * stdDevJ );
				}
			}

			// Return
			return obj;
		},
		
    	/**
    	* @function getVariances
    	*
    	* @summary Returns the variances associated to a covariance matrix.
    	*
    	* @description This function returns, as a n by 1 matrix, the diagonal elements (a_ii), i=1..n from the original
    	* square matrix (a_ij),i=1..n,j=1..n.
    	*
    	* @memberof Matrix_
		* @param {Matrix_} out an optional n by 1 matrix.
    	* @return {Matrix_} a n by 1 column matrix containing the variances vector associated to the covariance matrix,
		* either stored in the matrix out or in a new matrix.
    	*
    	* @example
    	* Matrix_([[1, 1, 8.1], [1, 16, 18], [8.1, 18, 81]]).getVariances();
    	* // Matrix_([1, 16, 81])
    	*/
		'getVariances': function(out) { 
			return this.diagonal(out);
		},
		
    	/**
    	* @function getStandardDeviations
    	*
    	* @summary Returns the standard deviations associated to a covariance matrix.
    	*
    	* @description This function returns, as a n by 1 matrix, the square of the diagonal elements SQRT(a_ii), i=1..n from the original
    	* square matrix (a_ij),i=1..n,j=1..n.
    	*
    	* @memberof Matrix_
		* @param {Matrix_} out an optional n by 1 matrix.
    	* @return {Matrix_} a n by 1 column matrix containing the standard deviations vector associated to the covariance matrix,
		* either stored in the matrix out or in a new matrix.
    	*
    	* @example
    	* Matrix_([[1, 1, 8.1], [1, 16, 18], [8.1, 18, 81]]).getStandardDeviations();
    	* //  Matrix_([1, 4, 9])
    	*/
		'getStandardDeviations': function(out) { 
			return this.diagonal(out).elemMap(function(i,j,val) { return Math.sqrt(val);}, out);
		},
		
    	/**
    	* @function standardizedGeneralizedVariance
    	*
    	* @summary Returns the standardized generalized variance of a covariance matrix.
    	*
    	* @description This function computes the standardized generalized variance of a covariance matrix,
    	* as described in the reference.
    	*
    	* @see <a href="http://dx.doi.org/10.1016/0047-259X(87)90153-9">Ashis SenGupta, Tests for standardized generalized variances of multivariate normal populations of possibly different dimensions, Journal of Multivariate Analysis, Volume 23, Issue 2, 1987, Pages 209-219</a>
    	* 
    	* @memberof Matrix_
    	* @return {number} the standardized generalized variance of the covariance matrix.
    	*
    	*/
    	'standardizedGeneralizedVariance': function () {
        	// Compute the determinant of the matrix
        	var det = this.determinant();
    		
    		// Check for positivity
    		if (det < 0) {
    		    throw new Error('covariance matrix is not positive semi-definite');
    		}
    		
    		// Compute the SGV if the determinant is positive, as per the formula of the reference.
    		var sgv = Math.pow(det, 1/this.nbRows);
    		
    		// Return it
    		return sgv;
    	},
	};

  // Addition of the methods to the input matrix
  for (var name in methods) {
    matrix[name] = methods[name];
  }
};


