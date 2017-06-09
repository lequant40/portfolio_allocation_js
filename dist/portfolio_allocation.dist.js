/**
 * @file Header
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Not to be used as is in Google Sheets */
var PortfolioAllocation = PortfolioAllocation || {};

PortfolioAllocation = (function(self) {

/* End Not to be used as is in Google Sheets */
/**
 * @file Functions related to matrix object.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function Matrix_
*
* @summary Construct a matrix object from an array of arrays of numbers.
*
* @description This function constructs an n by m matrix (a_ij),i=1..n,j=1..m from an array arr of n arrays of m real numbers, with coefficients a_ij 
* satisfying a_ij = arr[i-1][j-1].
* 
* @param {Array.<Array.<number>>} arr an array of n arrays of m real numbers with n and m natural integers greater than or equal to 1.
* @return {this} the constructed matrix.
*
* @example
* Matrix_([[1,2,3], [4,5,6]]);
*/
function Matrix_(arr) {
	// Checks
	if (!(arr instanceof Array)) {
		throw new Error('arr must be an array');
	}
	if (!(arr[0] instanceof Array)) {
		throw new Error('arr must be an array of arrays, arr[' + 0 + '] is not an array');
	}
	
	// Initialise the underlying variables
	this.nbRows = arr.length;
	this.nbColumns = arr[0].length;
	this.data = typeof Float64Array === 'function' ? new Float64Array(this.nbRows * this.nbColumns) : new Array(this.nbRows * this.nbColumns);
	
    // Fill the matrix
	for (var i = 0; i < this.nbRows; ++i) {
		if (!(arr[i] instanceof Array)) {
		throw new Error('arr must be an array of arrays, arr[' + i + '] is not an array');
		}
		if (arr[i].length !== this.nbColumns) {
			throw new Error('arr contains arrays of irregular length: index 0 - ' + this.nbColmuns + ', index ' + i + ' - ' + arr[i].length);
		}
		for (var j = 0; j < this.nbColumns; ++j) {
			this.data[i * this.nbColumns + j] = arr[i][j];
		}
	}
	
	// Return
	return this;
}


Matrix_.prototype = {
    constructor: Matrix_,
  
	/**
	* @function toString
	*
	* @summary Returns a string representation of the matrix.
	*
	* @description This function builds a readable string representation of the matrix, with automatically adjusted padding.
	* 
	* @memberof Matrix_
	* @return {string} a string representation of the matrix' contents.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).toString();
	* //
	* [  1  2  3 ]
    * [  4  5 10 ]
	*/
	toString: function() {
		// Calculate a correct padding for the display of the matrix
		var maxLen = 0;
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j=0; j < this.nbColumns; ++j) {
				var len = String(this.data[i * this.nbColumns + j]).length;
				if (len > maxLen) {
					maxLen = len;
				}
			}
		}
		
		// Build the string representation of the matrix
		var strMat = [];
		for (var i = 0; i < this.nbRows; ++i) {
			strMat.push('[ ');
			
			for (var j=0; j < this.nbColumns; ++j) {
				var strVal = String(this.data[i * this.nbColumns + j]);
				strMat.push(new Array(maxLen - strVal.length + 1).join(' ') + strVal + ' ');
			}
			
			strMat.push(']\n');
		}
		
		// Return it
		return strMat.join('');
	},
	
	
	/**
	* @function toArray
	*
	* @summary Returns an array representation of the matrix.
	*
	* @description This function builds an array representation of the matrix arr, satisfying arr[i-1][j-1] =  a_ij, i=1..nbRows, j=1..nbColumns.
	* 
	* @memberof Matrix_
	* @return {Array.<Array.<number>>} an array of array of real numbers, representing the matrix' contents.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).toArray();
	* //
	* [[1,2,3], [4,5,10]]
	*/
	toArray: function() {
		// Allocate the outer array
		var arr = new Array(this.nbRows);
		
		// Fill it
		for (var i = 0; i < this.nbRows; ++i) {
			// Allocate the inner array
			arr[i] = new Array(this.nbColumns);
			
			for (var j=0; j < this.nbColumns; ++j) {
				arr[i][j] = this.data[i * this.nbColumns + j];
			}
		}
		
		// Return it
		return arr;
	},
	
	
	/**
	* @function setValueAt
	*
	* @summary Sets the matrix coefficient at row i, column j.
	*
	* @description This function sets the matrix coefficient at row i, column j to the value val, with i belonging to 1..number of matrix rows,
	* and j belonging to 1..number of matrix columns.
	* 
	* @memberof Matrix_
	* @param {number} i the row index, natural integer belonging to 1..number of matrix rows
	* @param {number} j the column index, natural integer belonging to 1..number of matrix columns.
	* @param {number} val the value to set at row index i and column index j.
	* @return {this} the updated matrix.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).setValueAt(1,1,9).toString();
	* //
	* [  9  2  3 ]
    * [  4  5 10 ]
	*/
	setValueAt: function(i, j, val) {
		// Bounds check
		if (i < 1 || j < 1 || i > this.nbRows || j > this.nbColumns) {
			throw Error(
			'index out of bounds when setting matrix value, (' + i + ',' + j +
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		
		// Value setting
		this.data[(i-1) * this.nbColumns + (j-1)] = val;
		return this;
	},
	
	
	/**
	* @function getValueAt
	*
	* @summary Gets the matrix coefficient at row i, column j.
	*
	* @description This function gets the matrix coefficient at row i, column j, with i belonging to 1..number of matrix rows,
	* and j belonging to 1..number of matrix columns.
	* 
	* @memberof Matrix_
	* @param {number} i the row index, natural integer belonging to 1..number of matrix rows
	* @param {number} j the column index, natural integer belonging to 1..number of matrix columns.
	* @return {number} the value of the matrix coefficient at row i, column j.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).getValueAt(1,1);
	* // 1
	*/
	getValueAt: function(i, j) {
		// Bounds check
		if (i < 1 || j < 1 || i > this.nbRows || j > this.nbColumns) {
			throw Error(
			'index out of bounds when getting matrix value, (' + i + ',' + j +
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		
		// Value getting
		return this.data[(i-1) * this.nbColumns + (j-1)];
	},	

	
	/**
	* @function isSquare
	*
	* @summary Determines if the matrix is square.
	*
	* @description This function determines if the number of rows of the matrix is equal to the number of columns of the matrix.
	* 
	* @memberof Matrix_
	* @return {boolean} true if the matrix is square, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isSquare();
	* // false
	*
	* @example
	* Matrix_([[1,2], [4,5]]).isSquare();
	* // true
	*/
	isSquare: function() {
	    return this.nbRows === this.nbColumns;
	},
	

	/**
	* @function elementWisePower
	*
	* @summary Returns a matrix made of the elements of the original matrix raised to a power.
	*
	* @description This function computes a matrix (c_ij),i=1..n,j=1..m from the original matrix (a_ij),i=1..n,j=1..m, with coefficients
	* satisfying c_ij = a_ij^p.
	*
	* @memberof Matrix_
	* @param {number} p the exposant to which elevate the coefficients of the original matrix.
	* @return {Matrix_} the matrix with the original matrix coefficients elevated to the power p.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).elementWisePower(2);
	* // Matrix_([[1,4,9], [16,25,36]])
	*/
	elementWisePower: function (p) {
		// Result matrix instantiation
		var obj = Object.create(this);
		obj.nbRows = this.nbRows;
		obj.nbColumns = this.nbColumns;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
		
		// Computation of the power p of the coefficients of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = Math.pow(this.data[i * obj.nbColumns + j], p);
			}
		}
		
		// Return
		return obj;
	},
};


/**
* @function matrixIdentical_
*
* @summary Determines if two matrices are identical.
*
* @description This function determines if two matrices a and b are identical, i.e. if a_ij == b_ij, i=1..n,j=1..m where
* n is their common number of rows and m is their common number of rows.
*
* In case an optional tolerance parameter eps is provided, the strict component by component equality condition a_ij == b_ij
* is relaxed into |a_ij - b_ij| <= eps.
* 
* @param {Matrix_} a, a matrix.
* @param {Matrix_} b, a matrix.
* @param {number} eps, an optional real number; by default, eps is equal to zero.
* @return {boolean} true if a and b are identical, false otherwise.
*
* @example
* matrixIdentical_(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,2,3], [4,5,6]]));
* // true
*/
function matrixIdentical_(a, b, eps) {
	// Ensure a,b are matrices
	if (!(a instanceof Matrix_)) {
		throw new Error('a must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('b must be a matrix');
	}
	
	// Easy checks
	if (a.nbRows !== b.nbRows) {
		return false;
	}
	if (a.nbColumns !== b.nbColumns) {
		return false;
	}

	// Real checks
	var nbRows = a.nbRows;
	var nbColumns = a.nbColumns;
	var tol = eps || 0;
	for (var i = 0; i < nbRows; ++i) {
		for (var j = 0; j < nbColumns; ++j) {
			if (Math.abs(a.data[i * nbRows + j] - b.data[i * nbRows + j]) > tol) {
				return false;
			}
		}
	}

	//
	return true;
};


/**
* @function matrixMatrixProduct_
*
* @summary Returns the matrix-matrix product.
*
* @description This function computes the matrix-matrix product a*b of matrices a and b, where a is a n by m matrix and b is a m by p matrix.
*
* The algorithm implemented uses an IKJ form, cache aware.
* 
* @param {Matrix_} a an n by m matrix.
* @param {Matrix_} a an m by p matrix.
* @return {Vector_} the matrix product a*b, an m by p matrix.
*
* @example
* matrixMatrixProduct_(Matrix_([[1,2,3], [4,5,6]]), .Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
function matrixMatrixProduct_(a, b) {
	// Ensure a,b are matrices
	if (!(a instanceof Matrix_)) {
		throw new Error('a must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('b must be a matrix');
	}
	
	// Checks
	if (a.nbColumns !== b.nbRows) {
		throw new Error('Matrices sizes do not match: ' + '(' + a.nbRows + ',' + a.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns === 1) {
		throw new Error('Matrix vector multiplication unsupported');
	}
	
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = a.nbRows;
	obj.nbColumns = b.nbColumns;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Computation of A*B product in IKJ format, cache aware
	for (var i = 0; i < a.nbRows; ++i) {
		for (var j = 0; j < b.nbColumns; ++j) {
			obj.data[i * b.nbColumns + j] = 0;
		}
	
		for (var k = 0; k < a.nbColumns; ++k) {
			var aik = a.data[i * a.nbColumns + k];
			
			for (var j = 0; j < b.nbColumns; ++j) {
				obj.data[i * b.nbColumns + j] += aik * b.data[k * b.nbColumns + j];
			}
		}
	}
	
	// Return
    return obj;
};



/**
 * @file Functions related to vector object.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function Vector_
*
* @summary Construct a vector object from an array of numbers.
*
* @description This function constructs a vector (v_i),i=1..n of size n from an array arr of n real numbers, with coefficients v_i 
* satisfying v_i = arr[i-1].
* 
* @param {Array.<number>} arr an array of n real numbers with n natural integer greater than or equal to 1.
* @return {this} the constructed vector.
*
* @example
* Vector_([[1,2,3]);
*/
function Vector_(arr) {
	// Checks
	if (!(arr instanceof Array)) {
		throw new Error('arr must be an array');
	}

	// Initialise the underlying variables
	this.nbRows = arr.length;
	this.nbColumns = 1;
	this.data = typeof Float64Array === 'function' ? new Float64Array(this.nbRows * this.nbColumns) : new Array(this.nbRows * this.nbColumns);
	
    // Fill the vector
	for (var i = 0; i < this.nbRows; ++i) {
		this.data[i * this.nbColumns] = arr[i];
	}
	
	// Return
	return this;
}
Vector_.prototype = Object.create(Matrix_.prototype); // A vector inherit from all Matrix operations by default


/**
* @function toArray
*
* @summary Returns an array representation of the vector.
*
* @description This function builds an array representation of the vector arr, satisfying arr[i-1] = x_i, i=1..nbRows.
* 
* @memberof Vector_
* @return {<Array.<number>} an array of real numbers, representing the vector' contents.
*
* @example
* Vector_([1,2,3]).toArray();
* //
* [1,2,3]
*/
Vector_.prototype.toArray = function() { // Override Matrix_.toArray
	// Allocate the array
	var arr = new Array(this.nbRows);
	
	// Fill it
	for (var i = 0; i < this.nbRows; ++i) {
		arr[i] = this.data[i];
	}
	
	// Return it
	return arr;
};


/**
* @function normalize
*
* @summary Returns a vector made of the elements of the original vector normalized by their sum.
*
* @description This function computes a vector (v_i),i=1..n with elements v_i satisfying v_i = x_i/sum(x_j),j=1..n, 
* with n the length of the original vector x.
* 
* @memberof Vector_
* @return {Vector_} the normalized vector.
*
* @example
* Vector_([1,2,3]).normalize();
* // Vector_([1/6,1/3,1/2])
*/
Vector_.prototype.normalize = function() {
	// Result vector instantiation
	var obj = Object.create(Vector_.prototype);
	obj.nbRows = this.nbRows;
	obj.nbColumns = this.nbColumns;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);

	// Computation of sum x_i, i=1..nbRows
	var sum = this.sum();
	if (sum === 0.0) {
	    throw new Error('sum of coefficients of vector is null');
	}
	
	// Normalization of the vector
	for (var i = 0; i < this.nbRows; ++i) {
		obj.data[i] = this.data[i] / sum; 
	}
	
	// Return
	return obj;
};


/**
* @function sum
*
* @summary Returns the sum of the elements of a vector.
*
* @description This function computes the sum of the elements of a vector.
* 
* @memberof Matrix_
* @return {number} the sum of the elements of the vector x.
*
* @example
* Vector_([1,2,3]).sum();
* // 6
*/
Vector_.prototype.sum = function() {
	// Computation of sum x_i, i=1..nbRows
	var sum = 0;
	for (var i = 0; i < this.nbRows; ++i) {
		sum += this.data[i]; 
	}
	
	// Return
	return sum;
};


/**
* @function matrixVectorProduct_
*
* @summary Returns the matrix-vector product.
*
* @description This function computes the matrix-vector product a*x, where a is a n by m matrix and x is a vector of size m.
* 
* @param {Matrix_} a an n by m matrix.
* @param {Vector_} x a vector of size m.
* @return {Vector_} the matix-vector product a*x.
*
* @example
* matrixVectorProduct_(Matrix_([[1,2,3], [4,5,6]]), Vector_([1,2,3]));
* // Vector_([14,32])
*/
function matrixVectorProduct_(a, x) {
	// Ensure a is a matrix, x is a vector
	if (!(a instanceof Matrix_)) {
		throw new Error('a must be a matrix');
	}
	if (!(x instanceof Vector_)) {
		throw new Error('x must be a vector');
	}
	
	// Checks
	if (a.nbColumns !== x.nbRows) {
		throw new Error('Matrice size does not match vector size: ' + '(' + a.nbRows + ',' + a.nbColumns + 
		') - ' + '(' + x.nbRows + ',' + x.nbColumns + ')');
	}
	
	// Result vector instantiation
	var obj = Object.create(Vector_.prototype);
	obj.nbRows = a.nbRows;
	obj.nbColumns = 1;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Computation of A*x product
	for (var i = 0; i < a.nbRows; ++i) {
		obj.data[i] = 0;
		
		for (var j = 0; j < a.nbColumns; ++j) {
			obj.data[i] += a.data[i * a.nbColumns + j] * x.data[j] ;
		}
	}
	
	// Return the vector
    return obj;
}


/**
* @function vectorDotProduct_
*
* @summary Returns the vector dot product.
*
* @description This function computes the vector dot product <x/y>, where x and y are vectors of the same size.
* 
* @param {Vector_} x a vector.
* @param {Vector_} y a vector of same size as x.
* @return {number} the product <x/y>.
*
* @example
* vectorDotProduct_(Vector_([1,2,3]), Vector_([1,2,3]));
* // 14
*/
function vectorDotProduct_(x, y) {
	// Ensure x,y are vectors
	if (!(x instanceof Vector_)) {
		throw new Error('x must be a vector');
	}
	if (!(y instanceof Vector_)) {
		throw new Error('y must be a vector');
	}

	// Checks
	if (x.nbRows !== y.nbRows) {
		throw new Error('Vectors sizes do not match: ' + '(' + x.nbRows + ') - ' + '(' + y.nbRows + ')');
	}
	
	// Computation of <x/y>
	var dotProd = 0;
	for (var i = 0; i < x.nbRows; ++i) {
		dotProd += x.data[i] * y.data[i]; 
	}
	
	// Return it
	return dotProd;
}


/**
* @function vectorOnes_
*
* @summary Returns a vector made of ones.
*
* @description This function builds a vector (v_i),i=1..n of size n and satisfying v_i = 1, i=1..n.
* 
* @param {number} n the row length of the vector to construct.
* @return {Vector_} the constructed vector.
*
* @example
* vectorOnes_(3);
* // Vector_([1,1,1])
*/
function vectorOnes_(n) {
	// Result vector instantiation
	var obj = Object.create(Vector_.prototype);
	obj.nbRows = n;
	obj.nbColumns = 1;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Vector filling
	for (var i = 0; i < obj.nbRows; ++i) {
		obj.data[i] = 1;
	}
	
	// Return
	return obj;
}



/**
 * @file Functions related to equal risk budget portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function equalRiskBudgetWeights
*
* @summary Compute the weights of the equally weighted risk budget portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk budgets, defined as w_i = 1/sigma_i / (sum(1/sigma_j), j=1..n), i=1..n, with sigma_i the standard deviation 
* of the asset i.
*
* These weights ensure that the risk budget of each asset, defined as the product of the asset’s weight and its standard deviation, is equal.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio for each stock is the same and all pair-wise correlations are equal.
* 
* @see <a href="https://ssrn.com/abstract=1949003">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* 
* @param {<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, array of n real numbers statisfying sigma[i-1] = sigma_i.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to the equal risk budget portfolio, array of n real numbers.
*
* @example
* equalRiskBudgetWeights([0.1, 0.2]);
* // ~[0.59, 0.41]
*/
self.equalRiskBudgetWeights = function (sigma, opt) {
	// TODO: Checks, if enabled
	// Check that the values of sigma are strictly positive
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Vector_(sigma).elementWisePower(-1/2).normalize();
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to equal risk contributions portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function equalRiskContributionWeights
*
* @summary Compute the weights of the equally weighted risk contribution portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk contributions, defined as the weights with the property that the contribution 
* of each asset to the risk of the portfolio is equal.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* This portfolio is a trade-off between the MV portfolio and the EW portfolio.
* It can be viewed as an ERB portfolio tilted toward the assets less correlated with other assets.
*
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/jpm.2010.36.4.060">Maillard, S., Roncalli, T., Teiletche, J.: The properties of equally weighted risk contribution portfolios. J. Portf. Manag. 36, 60–70 (2010)</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* 
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 1000.
* @return {Array.<number>} the weights corresponding to the equal risk contributions portfolio, array of n real numbers.
*
* @example
* equalRiskContributionWeights([[0.1,0], [0,0.2]]);
* // ~[0.59, 0.41]
*/
self.equalRiskContributionWeights = function (sigma, opt) {
	// The ERC portfolio is a specific case of the more general risk budgeting portfolio, with equal risk budgets.
	//
	// Generate equal risk budgets: rb_i = 1/nbAssets, i=1..nbAssets
	var nbAssets = sigma.length;
	var rb = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
	    rb[i] = 1/nbAssets;
    }

	// Compute the associated risk budgeting weights
	return self.riskBudgetingWeights(sigma, rb, opt);
}


/**
 * @file Functions related to equal weights portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function equalWeights
*
* @summary Compute the weights of the equally weighted portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal weights, defined as w_i = 1/n, i=1..n.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the returns and volatility are the same for all stocks and all pair-wise correlations are equal.
* 
* @see <a href="https://academic.oup.com/rfs/article-abstract/22/5/1915/1592901/Optimal-Versus-Naive-Diversification-How">Victor DeMiguel, Lorenzo Garlappi, Raman Uppal; Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?. Rev Financ Stud 2009; 22 (5): 1915-1953. doi: 10.1093/rfs/hhm075</a>
* 
* @param {number} nbAssets the number of assets in the universe, natural integer superior or equal to 1.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to the equally weighted portfolio, array of real numbers of length nbAssets.
*
* @example
* equalWeights(5);
* // [0.2, 0.2, 0.2, 0.2, 0.2]
*/
self.equalWeights = function (nbAssets, opt) {
	// TODO: Checks, if enabled
	// Check that nbAssets is a strictly positive natural integer

	// Allocate the weights vector: all the weights are equal to 1/nbAssets
	var weights = vectorOnes_(nbAssets).normalize();

	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to risk budgeting portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function riskBudgetingWeights
*
* @summary Compute the weights of the risk budgeting portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with risk budgeting contraints, defined as the weights with the property that the contribution 
* of each asset to the risk of the portfolio is equal to a pre-determined budget weight.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="https://ssrn.com/abstract=2009778">Bruder, Benjamin and Roncalli, Thierry, Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012).</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* 
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {Array.<number>} rb the risk budgets, array of n real strictly positive numbers summing to one.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 1000.
* @return {Array.<number>} the weights corresponding to the risk budgeting portfolio, array of n real numbers.
*
* @example
* riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75], {eps: 1e-10, maxIter: 10000});
* // [~0.45, ~0.55]
*/
self.riskBudgetingWeights = function (sigma, rb, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 1000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Convert rb to vector format
	var rb = new Vector_(rb);

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that rb contains strictly positive numbers summing to one
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = vectorOnes_(nbAssets).normalize();
	
	// Preparational computations
	var sigma_x = matrixVectorProduct_(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(vectorDotProduct_(sigma_x, x)); // sigma(x)
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma and b
	var iter = 0;
	var converged = false;
	while (!converged) {
        // Convergence condition is false if any of the coordinate-wise convergence condition is false
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {
    	    // Define the coefficients of the second order polynomial ax_i^2 + b_ix + c
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i); // (SIGMA*x)_i - x_i*sigma_i^2, might be any sign
    	    var c = -rb.getValueAt(i,1) * s_x; // -b_i * sigma(x), always < 0
    	    
    	    // Extract the strictly positive root x_i^* of the equation ax_i^2 + bx_i + c = 0, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
			    throw new Error('discriminant not positive during iteration ' + iter + ', covariance matrix might not be psd');
			}
    	    q = -(b_p + sign_b_p * Math.sqrt(disc));
    	    r1 = q/a;
    	    r2 = c/q;
    	    xi_star = r1 > 0 ? r1 : r2;
    	    
    	    // Update the weights
    	    x.setValueAt(i,1,xi_star);
    	    
    	    // Compute the updated SIGMA*x and x'*SIGMA*x products for convergence condition evaluation + next loop evaluation
    	    sigma_x = matrixVectorProduct_(sigma, x)
    	    s_x = Math.sqrt(vectorDotProduct_(sigma_x, x));
    	    
    	    // Update the convergence condition: |RC_i* - b_i| <= eps, i = 1..nbAssets
    	    var rci_star = x.getValueAt(i,1) * sigma_x.getValueAt(i,1) / s_x;
    	    if (Math.abs(rci_star - rb.getValueAt(i,1)) > eps) {
    	        converged = false;
    	    }
    	    
    	    // Update the number of iterations
    	    ++iter;
    	    
    	    // Check the number of iterations
    	    if (iter >= maxIterations) {
    	        throw new Error('maximum number of iterations reached: ' + maxIterations);
    	    }
        }
	}
	
	// Return the computed weights, after normalization
	x = x.normalize();
	return x.toArray();
}


/**
 * @file Footer
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Not to be used as is in Google Sheets */
   return self;
  
})(PortfolioAllocation || {});

 
if (typeof module !== 'undefined') {
  module.exports = PortfolioAllocation;
}

/* End Not to be used as is in Google Sheets */
