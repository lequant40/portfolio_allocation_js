/**
 * @file Functions related to matrix object.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.Matrix_ = function(arr) { return new Matrix_(arr); }
self.matrixIdentical_ = function(a,b,eps) { return matrixIdentical_(a,b,eps); }
self.matrixMatrixProduct_ = function(a,b) { return matrixMatrixProduct_(a,b); }
/* End Wrapper private methods - Unit tests usage only */


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
	this.data = FLOAT64_ARRAY(this.nbRows * this.nbColumns);
	
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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


