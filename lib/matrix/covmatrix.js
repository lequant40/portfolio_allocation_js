/**
 * @file Functions related to covariance matrix computations for portfolio allocation.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function covarianceMatrix
*
* @summary Returns the covariance matrix of series of values.
*
* @description This function computes the covariance matrix of series of values, provided as 
* a variable number of arrays of real numbers of the same length.
*
* @param {...Array.<number>} var_args, arrays of real numbers of the same length.
* @return {Matrix_} a Matrix object representing the covariance matrix of the input series of values.
*
* @example
* covarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
* // == Matrix_([[0.00036,  -0.00053], [-0.00053, 0.00107]])
*/
self.covarianceMatrix = function(varg_args) {
	// Construct the covariance matrix
	var arrs = arguments;
	
	// Result matrix allocation
	var obj = allocateMatrix_(arrs.length, arrs.length);

	// Computation of the correlation matrix
	for (var i = 0; i < obj.nbRows; ++i) {
		// Copy from upper triangular part
		for (var j = 0; j < i; ++j) {
			obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
		}
		
		// Computation part
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = covariance_(arrs[i], arrs[j]);
		}
	}

	// Add covariance matrix methods
	addCovarianceMatrixMethods_(obj);
	
	// Return it
	return obj;
}

/**
* @function sampleCovarianceMatrix
*
* @summary Returns the sample covariance matrix of series of values.
*
* @description This function computes the sample covariance matrix of series of values, provided as 
* a variable number of arrays of real numbers of the same length.
*
* @param {...Array.<number>} var_args, arrays of real numbers of the same length.
* @return {Matrix_} a Matrix object representing the sample covariance matrix of the input series of values.
*
* @example
* sampleCovarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
* // == Matrix_([[0.00053, -0.0008], [-0.0008, 0.0016]])
*/
self.sampleCovarianceMatrix = function(varg_args) {
	// Construct the covariance matrix
	var arrs = arguments;
	
	// Result matrix allocation
	var obj = allocateMatrix_(arrs.length, arrs.length);

	// Computation of the correlation matrix
	for (var i = 0; i < obj.nbRows; ++i) {
		// Copy from upper triangular part
		for (var j = 0; j < i; ++j) {
			obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
		}
		
		// Computation part
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = sampleCovariance_(arrs[i], arrs[j]);
		}
	}

	// Add covariance matrix methods
	addCovarianceMatrixMethods_(obj);
	
	// Return it
	return obj;
}

/**
* @function toCovarianceMatrix
*
* @summary Returns a copy of the original matrix to which is added the methods of a covariance matrix.
*
*@description This function computes a copy of the original matrix and adds the methods of a covariance matrix to it.
*
* @memberof Matrix_
* @return {Matrix_} a copy of the original matrix to which is added the methods of a covariance matrix.
*
* @example
* Matrix_([[1,0.1], [0.1,1]]).toCovarianceMatrix();
* // == Matrix_([[1,0.1], [0.1,1]]) with covariance matrix methods
*/
Matrix_.prototype.toCovarianceMatrix = function() {
	// Construct the covariance matrix
	var cov = new Matrix_(this);
	
	// No checks: square, symetric, semidefinite positive, etc
	
	// Add covariance matrix methods
	addCovarianceMatrixMethods_(this);
	
	// Return it
	return this;
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
    	* matrix (a_ij),i=1..n,j=1..n, with coefficients satisfying c_ij = a_ij/(a_ii * a_jj).
    	*
    	* @memberof Matrix_
		* @param {Matrix_} out an optional n by n matrix.
    	* @return {Matrix_} a n by n matrix containing the correlation matrix associated to the covariance matrix,
		* either stored in the matrix out or in a new matrix.
    	*
    	* @example
    	* fromDoubleArray([[1,0.1], [0.1,1]]).getCorrelationMatrix();
    	* // Matrix_([[1,0.1], [0.1,1]])
    	*/
		'getCorrelationMatrix': function(out) { 
			// Result matrix allocation
			var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
			
			// Computation of the correlation matrix
			for (var i = 0; i < obj.nbRows; ++i) {
				// Standard devation of a_ii
				var stdDevI = Math.sqrt(this.data[i * this.nbColumns + i]);
				
				// Copy from upper triangular part
				for (var j = 0; j < i; ++j) {
					obj.data[i * obj.nbColumns + j] = obj.data[j * this.nbColumns + i];
				}
				
				// Computation part
				for (var j = i; j < obj.nbColumns; ++j) {
					// Standard devation of a_jj
					var stdDevJ = Math.sqrt(this.data[j * this.nbColumns + j]);
				
					obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j] / ( stdDevI * stdDevJ );
				}
			}

			// Return
			return obj;
		},
		
    	/**
    	* @function getVariancesVector
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
    	* fromDoubleArray([[1,0.1], [0.1,1]]).getVariancesVector();
    	* // Vector_([1,1])
    	*/
		'getVariancesVector': function(out) { 
			return this.diagonal(out);
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
    	* @example
    	* Matrix_([[1,2], [2,1]]).standardizedGeneralizedVariance();
    	* // XX
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


