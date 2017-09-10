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
* @summary Construct a matrix object from several possible inputs: an array of arrays of numbers, an array of numbers or a matrix object.
*
* @description This function constructs an n by m matrix (a_ij),i=1..n,j=1..m from either:
* - An array dblarr of n arrays of m real numbers, with coefficients a_ij satisfying a_ij = dblarr[i-1][j-1]
* - An array arr of n real numbers, with coefficients a_ij satisfying a_ij = arr[i-1] with m=1 (i.e., the constructed matrix is a column matrix, that is, a vector)
* - An n by m matrix (b_ij),i=1..n,j=1..m, with coefficients a_ij satisfying a_ij = b_ij
* 
* @param {Array.<Array.<number>>|<Array.<number>|Matrix_} input either an array of n arrays of m real numbers, or an array of n real numbers or a n by m matrix object,
* with n and m natural integers greater than or equal to 1.
* @return {this} the constructed matrix.
*
* @example
* Matrix_([[1,2,3], [4,5,6]]);
*/
function Matrix_(input) {
	var that = this;
	
	function fromDoubleArray(dblarr) {
		// Initialise the underlying variables
		that.nbRows = dblarr.length;
		that.nbColumns = dblarr[0].length;
		that.data = typeof Float64Array === 'function' ? new Float64Array(that.nbRows * that.nbColumns) : new Array(that.nbRows * that.nbColumns);
		
		// Fill the matrix
		for (var i = 0; i < that.nbRows; ++i) {
			if (!(dblarr[i] instanceof Array)) {
				throw new Error('unsupported input type');
			}
			if (dblarr[i].length !== that.nbColumns) {
				throw new Error('unsupported input type');
			}
			for (var j = 0; j < that.nbColumns; ++j) {
				that.data[i * that.nbColumns + j] = dblarr[i][j];
			}
		}
		
		// Return
		return that;
	}
	
	function fromArray(arr) {
		// Initialise the underlying variables
		that.nbRows = arr.length;
		that.nbColumns = 1;
		that.data = typeof Float64Array === 'function' ? new Float64Array(that.nbRows * that.nbColumns) : new Array(that.nbRows * that.nbColumns);
		
		// Fill the vector
		for (var i = 0; i < that.nbRows; ++i) {
			that.data[i * that.nbColumns] = arr[i];
		}
		
		// Return
		return that;
	}
	
	function fromMatrix(mat) {
		// Initialise the underlying variables
		that.nbRows = mat.nbRows;
		that.nbColumns = mat.nbColumns;
		that.data = typeof Float64Array === 'function' ? new Float64Array(that.nbRows * that.nbColumns) : new Array(that.nbRows * that.nbColumns);
		
		// Fill the matrix
		for (var i = 0; i < that.nbRows; ++i) {
			for (var j = 0; j < that.nbColumns; ++j) {
				that.data[i * that.nbColumns + j] = mat.data[i * mat.nbColumns + j];
			}
		}
		
		// Return
		return that;
	}
	
	// Checks
	if (input instanceof Array && input[0] instanceof Array) { // Standard matrix
		return fromDoubleArray(input);
	}
	else if (input instanceof Array) { // Simplified constructor for a column matrix (i.e., a vector)
		return fromArray(input);
	}
	else if (input instanceof Matrix_) { // Avoid conversion from Matrix to arrays to Matrix
		return fromMatrix(input);
	}
	else {
		throw new Error('unsupported input type');
	}
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
	* @function toDoubleArray
	*
	* @summary Returns a double array containing all the elements of the matrix satisfying a predicate function.
	*
	* @description This function builds a double array arr from the matrix (a_ij),i=1..nbRows,j=1..nbColumns with arr[i-1] 
	* containing all the elements a_ij,j=1..nbColumns of the i-th row of the matrix satisfying fct(i, j, a_ij) == true
	* where fct is an optional predicate function.
	*
	* In case fct is not provided, it defaults to always true, i.e., all the elements of the matrix are selected.
	* 
	* @memberof Matrix_
	* @param {function(number, number, number): number} fct the optional function to call on each element of the original matrix,
	* which should take 3 arguments: row index i=1..n, column index j=1..m, matrix element a_ij and which
	* should return a boolean, true in order to include the element a_ij in the output double array, false to skip it.
	* @return {Array.<Array.<number>>} an array of array of real numbers, representing the matrix' elements selected by the function fct.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).toDoubleArray();
	* //
	* [[1,2,3], [4,5,10]]
	*/
	toDoubleArray: function(fct) {
		// In case the function is not provided
		fct = fct || function(i, j, val) { return true; };
		
		// Allocate the outer array
		var arr = new Array(this.nbRows);
		
		// Fill it with matching matrix elements
		for (var i = 0; i < this.nbRows; ++i) {
			// Allocate the inner array
			arr[i] = new Array(this.nbColumns); // potential size, to be restrained after
			
			var k =0;
			for (var j=0; j < this.nbColumns; ++j) {
				var val = this.data[i * this.nbColumns + j];
				if (fct(i+1, j+1, val)) {
					arr[i][k++] = val;
				}
			}			
			arr[i].length = k; // Restrain the inner array size to effectively matching elements
		}
		
		// Return it
		return arr;
	},
	
	
	/**
	* @function toArray
	*
	* @summary Returns an array containing all the elements of the matrix satisfying a predicate function.
	*
	* @description This function builds an array arr from the matrix (a_ij),i=1..nbRows,j=1..nbColumns
	* containing all the elements a_ij,i=1..nbRows,j=1..nbColumns of the matrix satisfying fct(i, j, a_ij) == true
	* where fct is an optional predicate function.
	*
	* In case fct is not provided, it defaults to always true, i.e., all the elements of the matrix are selected.
	*
	* The order of the elements inside the array arr follows a row major organisation.
	* 
	* @memberof Matrix_
	* @param {function(number, number, number): number} fct the optional function to call on each element of the original matrix,
	* which should take 3 arguments: row index i=1..n, column index j=1..m, matrix element a_ij and which
	* should return a boolean, true in order to include the element a_ij in the output array, false to skip it.
	* @return {<Array.<number>} an array of real numbers, representing the matrix' elements selected by the function fct.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).toArray();
	* //
	* [1,2,3,4,5,10]
	*/
	toArray: function(fct) {
		// In case the function is not provided
		fct = fct || function(i, j, val) { return true; };
		
		// Allocate the array
		var arr = new Array(this.nbRows * this.nbColumns);
		
		// Fill it with matching matrix elements
		var k =0;
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j=0; j < this.nbColumns; ++j) {
				var val = this.data[i * this.nbColumns + j];
				if (fct(i+1, j+1, val)) {
					arr[k++] = val;
				}
			}			
		}
		arr.length = k; // Restrain the array size to effectively matching elements
		
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
	* @function isVector
	*
	* @summary Determines if the matrix is a (column) vector.
	*
	* @description This function determines if the number of columns of the matrix is equal to 1.
	* 
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a column vector, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isVector();
	* // false
	*
	* @example
	* Matrix_([[1], [4]]).isVector();
	* // true
	*/
	isVector: function() {
	    return this.nbColumns === 1;
	},

	/**
	* @function sum
	*
	* @summary Returns the sum of the elements of the matrix.
	*
	* @description This function computes the sum of the elements of the matrix.
	* 
	* @memberof Matrix_
	* @return {number} the sum of the elements of the matrix.
	*
	* @example
	* Matrix_([[1,2,3]]).sum();
	* // 6
	*/
	sum: function() {
		// Computation of sum a_ij, i=1..nbRows, j=1..nbColumns
		var sum = 0;
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				sum += this.data[i * this.nbColumns + j];
			}
		}
		
		// Return
		return sum;
	},
	
	/**
	* @function normalize
	*
	* @summary Returns a matrix made of the elements of the original matrix divided by their common sum.
	*
	* @description This function computes a matrix (c_ij),i=1..n,j=1..m with elements c_ij satisfying c_ij = a_ij/(sum(a_kl),k=1..n,l=1..m),i=1..n,j=1..m, 
	* with n the number of rows of the original matrix and m the number of columns of the original matrix.
	* 
	* @memberof Matrix_
	* @return {Matrix_} the normalized matrix.
	*
	* @example
	* Matrix_([[1,2,3]]).normalize();
	* // Matrix_([[1/6,1/3,1/2]])
	*/
	normalize: function() {
		// Result matrix instantiation
		var obj = Object.create(Matrix_.prototype);
		obj.nbRows = this.nbRows;
		obj.nbColumns = this.nbColumns;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);

		// Computation of the sum of the matrix elements
		var sum = this.sum();
		if (sum === 0.0) {
			throw new Error('sum of coefficients of matrix is null');
		}
		
		// Normalization of the matrix
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j] / sum; 
			}
		}
		
		// Return
		return obj;
	},

	/**
	* @function elemPower
	*
	* @summary Returns a matrix made of the elements of the original matrix raised to a power.
	*
	* @description This function computes a matrix (c_ij),i=1..n,j=1..m from the original matrix (a_ij),i=1..n,j=1..m, with coefficients
	* satisfying c_ij = a_ij^p.
	*
	* @memberof Matrix_
	* @param {number} p the exposant to which raise the coefficients of the original matrix.
	* @return {Matrix_} a matrix with the original matrix coefficients raised to the power p.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).elemPower(2);
	* // Matrix_([[1,4,9], [16,25,36]])
	*/
	elemPower: function (p) {
		// Result matrix instantiation
		var obj = Object.create(Matrix_.prototype);
		obj.nbRows = this.nbRows;
		obj.nbColumns = this.nbColumns;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
		
		// Computation of the power p of the coefficients of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = Math.pow(this.data[i * this.nbColumns + j], p);
			}
		}
		
		// Return
		return obj;
	},
	
	/**
	* @function getDiagonal
	*
	* @summary Returns the diagonal elements of a square matrix.
	*
	* @description This function returns, as a vector of n rows, the diagonal elements (a_ii), i=1..n from the original square matrix (a_ij),i=1..n,j=1..n.
	*
	* @memberof Matrix_
	* @return {Matrix_} a column matrix (i.e., vector) containing the diagonal coefficients of the original matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).getDiagonal();
	* // Matrix_([1,5])
	*/
    getDiagonal: function () {
    	// Checks
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
    	
    	// Result vector instantiation
    	var obj = Object.create(Matrix_.prototype);
    	obj.nbRows = this.nbRows;
    	obj.nbColumns = 1;
    	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
    	
    	// Extraction of diagonal elements of A
    	for (var i = 0; i < this.nbRows; ++i) {
    		obj.data[i] = this.data[i * (this.nbColumns + 1)] ;
    	}
    	
    	// Return
        return obj;
    },
    
    /**
	* @function elemMap
	*
	* @summary Returns a matrix made of the elements of the original matrix transformed by a function.
	*
	* @description This function computes a matrix (c_ij),i=1..n,j=1..m from the original matrix (a_ij),i=1..n,j=1..m, with coefficients
	* satisfying c_ij = fct(a_ij).
	*
	* @memberof Matrix_
	* @param {function(number, number, number): number} fct the function to call on each element of the original matrix,
	* which should take 3 arguments: row index i=1..n, column index j=1..m, original matrix element a_ij and which
	* should return a number, which will be inserted into the new matrix.
	* @return {Matrix_} a matrix with the original matrix coefficients transformed by the function fct.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).elemMap(function(i,j,val) { return i+j;});
	* // Matrix_([[2,3,4], [3,4,5]])
	*/
	elemMap: function (fct) {
		// Result matrix instantiation
		var obj = Object.create(Matrix_.prototype);
		obj.nbRows = this.nbRows;
		obj.nbColumns = this.nbColumns;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
		
		// Computation of the fonction fct applied to the coefficients of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = fct(i+1, j+1, this.data[i * this.nbColumns + j]);
			}
		}
		
		// Return
		return obj;
	},
	
    /**
	* @function getSubmatrix
	*
	* @summary Returns a (possibly non contiguous) submatrix from the original matrix, keeping the elements whose row and column indexes are specified.
	*
	* @description This function computes a matrix (c_ij),i=1..p,j=1..q from the original matrix (a_ij),i=1..n,j=1..m and from the lists of row/column indexes to keep
	* rindexes/cindexes, where p is equal to the length of rindexes and q is equal to the length of cindexes, with coefficients satisfying c_ij = a_rindexes[i]cindexes[j].
	*
	* @memberof Matrix_
	* @param {Array.<number>} rindexes the row indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..number of matrix rows
    * @param {Array.<number>} cindexes the column indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..number of matrix columns
	* @return {Matrix_} a matrix whose elements correspond to the elements of the original matrix whose row/column indexes belong to the input lists of row/column indexes to keep.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).getSubmatrix([1], [2, 3];
	* // Matrix_([[2,3]])
	*/
    getSubmatrix : function(rindexes, cindexes) {
    	// Check that indexes are arrays
    	if (!(rindexes instanceof Array) || rindexes.length == 0) {
    		throw new Error('first parameter must be a non empty array');
    	}
    	if (!(cindexes instanceof Array) || cindexes.length == 0) {
    		throw new Error('second parameter must be a non empty array');
    	}
	
    	// Check that the indexes are sorted arrays of distinct elements
    	for (var i = 1; i < rindexes.length; ++i) {
    	    if (rindexes[i-1] >= rindexes[i]) {
    	        throw new Error('first parameter must be a sorted array');
    	    }
    	    
    	}
        for (var j = 1; j < cindexes.length; ++j) {
    	    if (cindexes[j-1] >= rindexes[j]) {
    	        throw new Error('second parameter must be a sorted array');
    	    }
    	    
    	}
    	
    	// Result matrix instantiation
    	var obj = Object.create(this);
    	obj.nbRows = rindexes.length;
    	obj.nbColumns = cindexes.length;
    	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
    	
    	// Computation of the elements of A
    	for (var i = 0; i < obj.nbRows; ++i) {
    		var rindex = rindexes[i] - 1; // row index of the original matrix
    		
    		for (var j = 0; j < obj.nbColumns; ++j) {
    		    var cindex = cindexes[j] - 1; // column index of the original matrix
    		    
    			obj.data[i * obj.nbColumns + j] = this.data[rindex * this.nbColumns + cindex];
    		}
    	}
    	
    	// Return
        return obj;
    },
	
    /**
	* @function transpose
	*
	* @summary Returns the transpose of the original matrix.
	*
	* @description This function computes the transpose matrix (c_ij),i=1..n,j=1..m of the original matrix (a_ij),i=1..n,j=1..m, with coefficients
	* satisfying c_ij = a_ji.
	*
	* @memberof Matrix_
	* @return {Matrix_} the transpose of the original matrix.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).transpose();
	* // Matrix_([[1,4], [2,5], [3,6]])
	*/
	transpose: function () {
		// Result matrix instantiation
		var obj = Object.create(this);
		obj.nbRows = this.nbColumns;
		obj.nbColumns = this.nbRows;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
		
		// Computation of the transpose of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[j * this.nbColumns + i];
			}
		}
		
		// Return
		return obj;
	},
};


/**
* @function areEqual
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
* areEqual(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,2,3], [4,5,6]]));
* // true
*/
Matrix_.areEqual = function (a, b, eps) {
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
* @function product
*
* @summary Returns the matrix-matrix product.
*
* @description This function computes the matrix-matrix product a*b of matrices a and b, where a is a n by m matrix and b is a m by p matrix.
*
* The algorithm implemented uses an IKJ form, cache aware.
* 
* @param {Matrix_} a an n by m matrix.
* @param {Matrix_} a an m by p matrix.
* @return {Matrix_} the matrix product a*b, an m by p matrix.
*
* @example
* product(Matrix_([[1,2,3], [4,5,6]]), .Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.product = function(a, b) {
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
* @function diagonal
*
* @summary Returns a diagonal matrix constructed from an array of numbers.
*
* @description This function computes a diagonal square matrix (a_ij),i=1..n,j=1..n from an array arr of n real numbers, 
* with coefficients a_ij satisfying a_ij = 0,i <>j and a_ij = arr[i-1],i==j.
*
* @param {Vector_} x a vector of size n.
* @return {Matrix_} a n by n diagonal matrix with its diagonal elements corresponding to the elements of x.
*
* @example
* diagonal([1,2,3]);
* // == Matrix_([[1,0,0], [0,2,0], [0,0,3]])
*/
Matrix_.diagonal = function(arr) {
	// Checks
	if (!(arr instanceof Array)) {
		throw new Error('arr must be an array');
	}

	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = arr.length;
	obj.nbColumns = arr.length;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Computation of diag(arr)
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}	
		obj.data[i * (obj.nbColumns + 1)] = arr[i];
	}
	
	// Return
    return obj;
};


/**
* @function fillSymetric
*
* @summary Returns a symetric matrix constructed from a function.
*
* @description This function computes a diagonal square matrix (a_ij),i=1..n,j=1..n from a function, 
* with coefficients a_ij satisfying a_ij = a_ji = fct(i,j).
*
* The function fct is only called for indices i such that j >= i.
*
* @param {number} n, the order of the matrix to create, a strictly positive integer.
* @param {function(number, number): number} fct the function to call on each (i,j) pair of indexes with j >= i,
* which should take 2 arguments: row index i=1..n, column index j=i..n and which
* should return a number, which will be inserted into the matrix as its a_ij coefficient.
* @return {Matrix_} a n by n symetric matrix with its elements computed by the fonction fct.
*
* @example
* fillSymetric(2, function(i,j) { return 0; });
* // == Matrix_([[0,0], [0,0]])
*/
Matrix_.fillSymetric = function(n, fct) {
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = n;
	obj.nbColumns = n;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Computation of the elements of A
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < i; ++j) {
			obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
		}
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = fct(i+1 ,j+1);
		}	
	}
	
	// Return
    return obj;
};



/**
* @function zeros
*
* @summary Returns a matrix made of zeros.
*
* @description This function builds an n by m matrix (a_ij),i=1..n,j=1..m satisfying a_ij = 0,i=1..n,j=1..m.
* 
* @param {number} n the row length of the matrix to construct, natural integer greater than or equal to 1.
* @param {number} m the column length of the matrix to construct, natural integer greater than or equal to 1.
* @return {Matrix_} the constructed matrix.
*
* @example
* zeros(3, 2);
* // Matrix_([[0,0], [0,0], [0,0]])
*/
Matrix_.zeros = function(n, m) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows out of bounds: ' + n);
	}
	if (m < 1) {
		throw new Error('input number of columns out of bounds: ' + m);
	}
		
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = n;
	obj.nbColumns = m;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Initialisation of the elements of A to 0
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
	}
	
	// Return
    return obj;
}


/**
* @function ones
*
* @summary Returns a matrix made of ones.
*
* @description This function builds an n by m matrix (a_ij),i=1..n,j=1..m satisfying a_ij = 1,i=1..n,j=1..m.
* 
* @param {number} n the row length of the matrix to construct, natural integer greater than or equal to 1.
* @param {number} m the column length of the matrix to construct, natural integer greater than or equal to 1.
* @return {Matrix_} the constructed matrix.
*
* @example
* ones(3, 2);
* // Matrix_([[1,1], [1,1], [1,1]])
*/
Matrix_.ones = function(n, m) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows out of bounds: ' + n);
	}
	if (m < 1) {
		throw new Error('input number of columns out of bounds: ' + m);
	}
		
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = n;
	obj.nbColumns = m;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Initialisation of the elements of A to 0
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 1;
		}
	}
	
	// Return
    return obj;
}

/**
* @function vectorHadamardProduct
*
* @summary Returns the Hadamard product of two vectors.
*
* @description This function computes the Hadamard product x*y of two vectors x and y of the same size.
* 
* @param {Matrix_} x a column matrix.
* @param {Matrix_} y a column matrix of same size as x.
* @return {Matrix_} the Hadamard product x*y.
*
* @example
* vectorHadamardProduct(Vector_([1,2,3]), Vector_([1,2,3]));
* // Matrix_([[1],[4],[9]])
*/
Matrix_.vectorHadamardProduct = function(x, y) {
	// Ensure x,y are vectors
	if (!x.isVector()) {
		throw new Error('first argument must be a vector');
	}
	if (!y.isVector()) {
		throw new Error('second argument must be a vector');
	}

	// Checks
	if (x.nbRows !== y.nbRows) {
		throw new Error('Vectors sizes do not match: ' + '(' + x.nbRows + ') - ' + '(' + y.nbRows + ')');
	}
	
	// Result vector instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = x.nbRows;
	obj.nbColumns = 1;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	
	// Computation of x*y product
	for (var i = 0; i < x.nbRows; ++i) {
		for (var j = 0; j < y.nbRows; ++j) {
			obj.data[i] = x.data[i] * y.data[i];
		}
	}
	
	// Return the vector
    return obj;
}


/**
* @function vectorDotProduct
*
* @summary Returns the dot product of two vectors.
*
* @description This function computes the vector dot product <x/y>, where x and y are vectors of the same size.
* 
* @param {Matrix_} x a column matrix.
* @param {Matrix_} y a column matrix of same size as x.
* @return {number} the dot product <x/y>.
*
* @example
* vectorDotProduct(Matrix_([[1],[2],[3]]), Matrix_([[1],[2],[3]]));
* // 14
*/
Matrix_.vectorDotProduct = function(x, y) {
	// Ensure x,y are vectors
	if (!x.isVector()) {
		throw new Error('first argument must be a vector');
	}
	if (!y.isVector()) {
		throw new Error('second argument must be a vector');
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
 * @file Functions related to covariance matrix computations for portfolio allocation.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 



/**
* @function covarianceMatrix
*
* @summary Returns the covariance matrix of a series of values.
*
* @description This function computes the covariance matrix of a series of values, provided as 
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
	
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = arrs.length;
	obj.nbColumns =  arrs.length;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);

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
* @summary Returns the sample covariance matrix of a series of values.
*
* @description This function computes the sample covariance matrix of a series of values, provided as 
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
	
	// Result matrix instantiation
	var obj = Object.create(Matrix_.prototype);
	obj.nbRows = arrs.length;
	obj.nbColumns =  arrs.length;
	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);

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
    	* @return {Matrix_} the correlation matrix associated to the covariance matrix.
    	*
    	* @example
    	* fromDoubleArray([[1,0.1], [0.1,1]]).getCorrelationMatrix();
    	* // Matrix_([[1,0.1], [0.1,1]])
    	*/
		'getCorrelationMatrix': function() { 
			// Result matrix instantiation
			var obj = Object.create(Matrix_.prototype);
			obj.nbRows = this.nbRows;
			obj.nbColumns =  this.nbColumns;
			obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);

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
    	* @description This function returns, as a vector of n rows, the diagonal elements (a_ii), i=1..n from the original
    	* square matrix (a_ij),i=1..n,j=1..n.
    	*
    	* @memberof Matrix_
    	* @return {Vector_} the variances vector associated to the covariance matrix.
    	*
    	* @example
    	* fromDoubleArray([[1,0.1], [0.1,1]]).getVariancesVector();
    	* // Vector_([1,1])
    	*/
		'getVariancesVector': function() { 
			return this.getDiagonal();
		},
	};

  // Addition of the methods to the input matrix
  for (var name in methods) {
    matrix[name] = methods[name];
  }
};



/**
 * @file Misc. statistical functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
 


/**
* @function rank_
*
* @summary Returns the rank of each value in a serie of values.
*
* @description This function computes the rank of each value in a serie of values, which is computed 
* by first sorting the serie of values, either by ascending or descending order, and then by computing 
* the position of each value in the sorted serie.
*
* Duplicate values in the serie of values all have the same rank, defined as the bottom rank of these duplicate values.
*
* This function mimics the Excel function RANK.EQ.
*
* @param {Array.<number>} x an array of real numbers.
* @param {number} order an integer equals to 0 to sort the serie of values in descending order, or equals to 1 to sort the serie of values in ascending order.
* @return {Array.<number>} an array of real numbers of the same size as x, containing the rank of each value of x.
*
* @example
* rank_([12, 13, 15, 10, 12], 1);
* // [2, 4, 5, 1, 2]
*/
function rank_(x, order) {
	// Transform the input array into an array with indexes
	var xWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		xWithIndexes[i] = [x[i], i];
	}
	
	// Sort the transformed array
	if (order == 0) {
		xWithIndexes.sort(function(a, b) {
			return a[0] > b[0] ? -1 : 1;
		}); 
	}
	else if (order == 1) {
		xWithIndexes.sort(function(a, b) {
			return a[0] < b[0] ? -1 : 1;
		}); 
	}
	
	// Compute the ranks of the values, setting an equal rank for all identical values
	// and skipping the next ranks values
	var xRanks = new Array(x.length);
	xRanks[xWithIndexes[0][1]] = 1; // first rank is always 1
	for (var i = 1; i < x.length; ++i) {
	    if (xWithIndexes[i][0] == xWithIndexes[i-1][0]) {
	  	  xRanks[xWithIndexes[i][1]] = xRanks[xWithIndexes[i-1][1]];
	    }
	    else {
		  xRanks[xWithIndexes[i][1]] = i + 1;
	    }
	}
	
	// Returnt the computed ranks
	return xRanks;
}


/**
* @function ftca_
*
* @summary Compute a clustering representation of a universe of elements using the Fast Threshold Clustering Algorithm (FTCA).
*
* @description This function returns the clustering representation of a universe of n elements based on their correlation
* and a correlation threshold, as computed by the FTCA algorithm of David Varadi described in the reference, initially created
* to deal with financial assets.
*
* This algorithm has many desirable properties that traditional clustering algorithms do not:
* 1) it produces fairly stable clusters
* 2) it is fast and deterministic 
* 3) it is easy to understand. 
*
* By default, David Varadi used a correlation threshold of 0.5 (approximately the level of statistical significance) to separate similar from dissimilar elements (assets).
* The choice of the threshold will change the number and stability of the clusters, with higher thresholds showing more clusters and a greater change in membership than lower thresholds. 
*
* To be noted that the FTCA works similar to the Minimum Correlation Algorithm from the same author in that it uses the average correlation of each element (asset)
* to all other elements (asset)s as a means of determining how closely or distantly related an element (asset) is to the universe of elements (assets) chosen.
*
* @see <a href="https://cssanalytics.wordpress.com/2013/11/26/fast-threshold-clustering-algorithm-ftca/">Fast Threshold Clustering Algorithm (FTCA)</a>
* 
* @param {Matrix_} correlationMatrix the correlation matrix (rho_ij),i,j=1..n of n elements, square n by n Matrix where n is a strictly positive natural integer.
* @param {number} threshold the correlation threshold to use in the FTCA algorithm, a real number typically belonging to interval [-1, 1].
* @return {Array.<Array.<number>>} the list of clusters as computed by the FTCA algorithm, array of m arrays of strictly positive integers representing the indexes of the elements in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
*
* @example
* ftca_(new Matrix_([[1, 0], [0,1]]), 0.5);
*  // [[2],[1]]
*/
function ftca_(correlationMatrix, threshold) {
	// Decode the optional threshold
	var threshold = threshold;
	if (threshold === undefined) {
		threshold = 0.5;
	}
	
	// The list of output clusters, to be populated
	var clusters = [];

	// The list of elements indexes not assigned to any cluster, initialized with all elements indexes (initially, no clusters are existing)
	var nbElements = correlationMatrix.nbRows;
	var unassignedElementsIdx = new Array(nbElements);
	for (var i = 0; i < unassignedElementsIdx.length; ++i) {
		unassignedElementsIdx[i] = i + 1;
	}

	// While there are elements that have not been assigned to a cluster
	while (unassignedElementsIdx.length != 0) {
		// If only one element remaining then
		if (unassignedElementsIdx.length === 1) {
			// Add a new cluster
			// Only member is the remaining element, set as not unassigned anymore
			var newCluster = [unassignedElementsIdx[0]];
			unassignedElementsIdx[0] = null;
			
			// Effectively add the new cluster into the list of clusters
			clusters.push(newCluster);
		}	
		else {
			// Get the (sub)correlation matix of the unassigned elements
			var subCorrMat = correlationMatrix.getSubmatrix(unassignedElementsIdx, unassignedElementsIdx);
			
			// Compute the average correlation of each unassigned element to all the other unassigned elements
			// Computation is done for each row
			var subCorrMatRows = subCorrMat.toDoubleArray(function(i, j, val) {
				return i != j;
			});
			var avgCorrelation = new Array(unassignedElementsIdx);
			for (var i = 0; i < unassignedElementsIdx.length; ++i) {
				avgCorrelation[i] = mean_(subCorrMatRows[i]);
			}
				
			// Find the element with the Highest Average Correlation (HC) to all elements not yet been assigned to a Cluster
			// Find the element with the Lowest Average Correlation (LC) to all elements not yet assigned to a Cluster
			// Note: When only 2 elements are remaining, HC will be equal to LC
			var hc = 0;
			var hcIdx = -1;
			var highestAvgCorrelation = -1;
			var lc = 0;
			var lcIdx = -1;
			var lowestAvgCorrelation = 1;		
			for (var i = 0; i < unassignedElementsIdx.length; ++i) {
				if (avgCorrelation[i] >= highestAvgCorrelation) {
					hc = unassignedElementsIdx[i];
					hcIdx = i;
					highestAvgCorrelation = avgCorrelation[i];
				}
				if (avgCorrelation[i] <= lowestAvgCorrelation) {
					lc = unassignedElementsIdx[i];
					lcIdx = i;
					lowestAvgCorrelation = avgCorrelation[i];
				}
			}
			
			// If Correlation between HC and LC > Threshold
			if (correlationMatrix.getValueAt(hc, lc) > threshold) {
				// Add a new Cluster made of HC and LC and set these two elements as not unassigned anymore
				// (Unless HC == LC, which can happen, for instance when there are only two elements remaining)
				var newClusterHcLc = (hc === lc ? [hc] : [hc, lc]);
				unassignedElementsIdx[hcIdx] = null;
				unassignedElementsIdx[lcIdx] = null;
				
				// Add to Cluster all other elements that have yet been assigned to a Cluster and have an Average Correlation to HC and LC > Threshold
				// Note: In Systematic Investor R code, all remaining elements are put inthe HcLc cluster, disregarding the condition on the correlation above.
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned elements (HC and LC)
						var avgHcLcAssetCorrelation = (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) + correlationMatrix.getValueAt(unassignedElementsIdx[i], lc)) / 2;
						if (avgHcLcAssetCorrelation  > threshold) {
							newClusterHcLc.push(unassignedElementsIdx[i]);
							
							// Set the element as not unassigned anymore
							unassignedElementsIdx[i] = null;
						}
					}
				}
			   
				// Effectively add the new cluster into the list of clusters				
				clusters.push(newClusterHcLc);
			}
			// Else
			else {
				// Add a Cluster made of HC and set this element as not unassigned anymore
				var newClusterHc = [hc];
				unassignedElementsIdx[hcIdx] = null;
				
				// Add to Cluster all other assets that have yet been assigned to a Cluster and have a Correlation to HC > Threshold
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC)
						if (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) > threshold) {
							newClusterHc.push(unassignedElementsIdx[i]);
							
							// Set the element as not unassigned anymore
							unassignedElementsIdx[i] = null;
						}
					}
				}
				
				// Effectively add the new cluster into the list of clusters				
				clusters.push(newClusterHc);

				// Add a Cluster made of LC and set this element as not unassigned anymore
				// (Unless HC == LC, which can happen, for instance when there are only two elements remaining)
				if (hc !== lc) {
					// Note: At this stage, to be noted that the LC element cannot have been assigned to the Hc cluster above if LC <> HC, since
					// otherwise, it would mean corr(lc, hc) > threshold, which is incompatible with the "else" branch in which
					// the code currently is; Lc cluster is thus always non empty.
					var newClusterLc = [lc];
					unassignedElementsIdx[lcIdx] = null;
					
					// Add to Cluster all other assets that have yet been assigned to a Cluster and have Correlation to LC > Threshold
					for (var i = 0; i < unassignedElementsIdx.length; ++i) {
						if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC with its correlated assets, and LC)
							if (correlationMatrix.getValueAt(unassignedElementsIdx[i], lc) > threshold) {
								newClusterLc.push(unassignedElementsIdx[i]);
								
								// Set the element as not unassigned anymore
								unassignedElementsIdx[i] = null;
							}
						}
					}

					// Effectively add the new cluster into the list of clusters				
					clusters.push(newClusterLc);
				}
				
				// Note: In Systematic Investor R code, it is possible for an element to belong to the two clusters Hc and Lc, in which case Hc is the final cluster,
				// which is conform to the description of the reference.
			}
		}
		
		// Effectively remove the assigned elements indexes (now pointing to null)  from the list of unassigned elements
		var newUnassignedElementsIdx = [];
		for (var i = 0; i < unassignedElementsIdx.length; ++i) {
			if (unassignedElementsIdx[i] !== null) {
				newUnassignedElementsIdx.push(unassignedElementsIdx[i]);
			}
		}
		unassignedElementsIdx = newUnassignedElementsIdx;
	}

	// Return the computed list of clusters
	return clusters;
}



/**
* @function mean_
*
* @summary Compute the arithmetic mean of a serie of values.
*
* @description This function returns the arithmetic mean of a serie of values [x_1,...,x_p], 
* which is defined as the sum of the p values x_1,...,x_p, divided by p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
* 
* @param {Array.<number>} x an array of real numbers.
* @return {number} the arithmetic mean of the values of the array x.
*
* @example
* mean_([2,4]); 
* // 3
*/
function mean_(x) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the values of the input numeric array (first pass)
	var tmpMean = 0.0;
	var sum = 0.0;
	for (var i=0; i<nn; ++i) {
		sum += x[i];
	}
	tmpMean = sum/nn;

	// Compute the correction factor (second pass)
	// C.f. M_3 formula of the reference
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		sumDiff += (x[i] - tmpMean);
	}

	// Return the corrected mean
	return (sum + sumDiff)/nn;
}


/**
* @function variance_
*
* @summary Compute the variance of a serie of values.
*
* @description This function returns the variance of a serie of values [x_1,...,x_p], 
* which is defined as the arithmetic mean of the p values (x_1-m)^2,...,(x_p-m)^2, where m is the arithmetic mean
* of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the variance of the values of the array x.
*
* @example
* variance_([4, 7, 13, 16]); 
* // 22.5
*/
function variance_(x) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the input numeric array (first pass)
	var meanX = mean_(x);

	// Compute the squared deviations plus the correction factor (second pass)
	// C.f. S_4 formula of the reference
	var sumSquareDiff = 0.0;
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		var diff = (x[i] - meanX);
		sumSquareDiff += diff * diff;
		sumDiff += diff;
	}

	// Compute the corrected sum of squares of the deviations from the mean
	var S = sumSquareDiff - ((sumDiff * sumDiff) / nn);

	// Return the corrected variance
	return S/nn;
}


/**
* @function sampleVariance_
*
* @summary Compute the sample variance of a serie of values.
*
* @description This function returns the sample variance of a serie of values [x_1,...,x_p], 
* which is defined as the variance of the p values x_1,...,x_p multiplied by p/(p-1).
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function variance_.
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the variance of the values of the array x.
*
* @example
* sampleVariance_([4, 7, 13, 16]); 
* // 30
*/
function sampleVariance_(x) {
	var nn = x.length;
	return variance_(x) * nn/(nn - 1);
}


/**
* @function stddev_
*
* @description Compute the standard deviation of a serie of values.
*
* @description This function returns the standard deviation of a serie of values [x_1,...,x_p], 
* which is defined as the square root of the variance of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function variance_.
*
* @see <a href="https://en.wikipedia.org/wiki/Standard_deviation">https://en.wikipedia.org/wiki/Standard_deviation</a>
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the standard deviation of the values of the array x.
*
* @example
* stddev_([1, 2, 3, 4]); 
* // ~1.12
*/
function stddev_(x) {
	return Math.sqrt(variance_(x));
}


/**
* @function sampleStddev_
*
* @description Compute the sample standard deviation of a serie of values.
*
* @description This function returns the sample standard deviation of a serie of values [x_1,...,x_p], 
* which is defined as the square root of the sample variance of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function sampleVariance_.
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the standard deviation of the values of the array x.
*
* @example
* sampleStddev_([1, 2, 3, 4]); 
* // ~1.29
*/
function sampleStddev_(x) {
	return Math.sqrt(sampleVariance_(x));
}

/**
* @function normcdf_
*
* @summary Compute the the standard normal cumulative distribution function.
*
* @description This function returns an approximation of the standard normal cumulative distribution function, i.e.
* given x a real number, it returns an approximation to p = Pr{Z <= x} where Z is a
* random variable following a standard normal distribution law.
*
* This function is also called Phi in the statistical litterature.
*
* The algorithm uses a Taylor expansion around 0 of a well chosen function of Phi.
* The algorithm has an absolute error of less than 8e−16.
*
* @author George Marsaglia
*
* @see <a href="https://www.jstatsoft.org/article/view/v011i04/v11i04.pdf"> G. Marsaglia. Evaluating the normal distribution. Journal of Statistical Software, 11(4):1–11, 2004.</a>
* 
* @param {number} x a real number.
* @return {number} an approximation to the p value satisfying p = Pr{Z <= x} where Z is a random variable following a standard normal distribution law.
*
* @example
* normcdf_(0);
* // 0.5
*/
function normcdf_(x) {
	// Initialisations
	var s=x;
	var t=0;
	var b=x;
	var q=x*x;
	var i=1;

	// The main loop corresponds to the computation of the Taylor serie of the function B around 0, c.f. page 5 of the reference.
	while (s != t) {
		s = (t = s) + (b *= q/(i += 2));
	}

	// The formula linking Phi and the Taylor expansion above if Phi = 1/2 + normal density * B, c.f. page 5 of the reference.
	return 0.5 + s * Math.exp(-0.5 * q - 0.91893853320467274178)
}


/**
* @function covariance_
*
* @summary Compute the covariance of two serie of values.
*
* @description This function returns the covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p], 
* which is defined as the arithmetic mean of the p values (x_1-m_x)*(y_1-m_y),...,(x_p-m_x)*(y_p-m_y), 
* where m_x is the arithmetic mean of the p values x_1,...,x_p and m_y is the arithmetic mean of the p values y_1,...,y_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
*
* @param {Array.<number>} x an array of real numbers.
* @param {Array.<number>} y an array of real numbers of the same length as x.
* @return {number} the covariance of the values of the arrays x and y.
*
* @example
* covariance_([4, 7, 13, 16], [4, 7, 13, 16]); 
* // 22.5
*/
function covariance_(x, y) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the input numeric arrays (first pass)
	var meanX = mean_(x);
	var meanY = mean_(y);

	// Compute the sum of the product of the deviations plus the correction factor (second pass)
	// C.f. P_4 formula of the reference
	var sumProdDiff = 0.0;
	var sumDiffX = 0.0;
	var sumDiffY = 0.0;
	for (var i=0; i<nn; ++i) {
		var diffX = (x[i] - meanX);
		var diffY = (y[i] - meanY);
		sumProdDiff += diffX * diffY;
		sumDiffX += diffX;
		sumDiffY += diffY;
	}

	// Compute the corrected sum of the product of the deviations from the means
	var C = sumProdDiff - ((sumDiffX * sumDiffY) / nn);

	// Return the corrected covariance
	return C/nn;
}


/**
* @function sampleCovariance_
*
* @summary Compute the sample covariance of two serie of values.
*
* @description This function returns the sample covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p], 
* which is defined as the covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p] multiplied by p/(p-1).
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function covariance_.
*
* @param {Array.<number>} x an array of real numbers.
* @param {Array.<number>} y an array of real numbers of the same length as x.
* @return {number} the covariance of the values of the arrays x and y.
*
* @example
* sampleCovariance_([4, 7, 13, 16], [4, 7, 13, 16]); 
* // 30
*/
function sampleCovariance_(x, y) {
	var nn = x.length;
	return covariance_(x,y) * nn/(nn - 1);
}


/**
 * @file Functions related to cluster risk parity portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function clusterRiskParityWeights
*
* @summary Compute the weights of the cluster risk parity portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only 
* cluster risk parity portfolio of n assets, as computed by the cluster risk parity algorithm described in 
* the reference.
*
* This algorithm combines the use of a clustering method to isolate groups of assets and then allocate 
* both within and across these groups using equal risk contribution (ERC) weights.
*
* To be noted that the choice of the clustering method is not detailled in the reference, so that two possibilites
* are offered by this function:
* - (Default) Using automatically the Fast Threshold Clustering Algorithm (FTCA) from David Varadi, c.f. the ftca_ function.
* - Using a list of clusters provided by the user, typically constructed with the clustering algorithm of her choice
* 
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* @see <a href="https://cssanalytics.wordpress.com/2013/01/03/cluster-risk-parity/">Cluster Risk Parity</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square n by n Matrix or array of n arrays of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the ERC algorithms, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the ERC algorithms, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.clusteringMode the method to use for the clusters computation, a string either equals to 'ftca' (usage of the list of clusters automatically constructed by the FTCA algorithm from David Varadi) or to 'manual' (usage of a list of clusters provided in input in the opt.clusters option); defaults to 'ftca'.
* @param {number} opt.ftcaThreshold the correlation threshold to use in the FTCA algorithm in case opt.clusteringMode is equal to 'ftca', a real number beloning to [0,1]; defaults to 0.5.
* @param {Array.<Array.<number>>} opt.clusters the list of clusters to use in the algorithm in case opt.clusteringMode is equal to 'manual', an array of m arrays of strictly positive integers representing the indexes of the assets in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
* @return {Array.<number>} the weights corresponding to the cluster risk parity portfolio, array of n real numbers.
*
* @example
* clusterRiskParityWeights([[0.1,0], [0,0.2]], {clusteringMode: 'none', clusters: [[1], [2]]});
*  
*/
self.clusterRiskParityWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	
	// Clustering options
	var clusteringMode = opt.clusteringMode || 'ftca';
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma).toCovarianceMatrix();
	var nbAssets = sigma.nbRows;

	
	// ------
	// Clustering logic
	var clusters = [];
	
	// In case clusters are provided, check that they form a partition of the set [1..nbAssets]
	if (clusteringMode === 'manual') {
		// Decode the none clustering options
		clusters = opt.clusters || [];
		
		// Prepare a partition of the set of integers [1..nbAssets]
		var partition = new Array(nbAssets);
		for (var i = 0; i < partition.length; ++i) {
		    partition[i] = 0;
		}
		
		// Count the number of times each asset index appears in the list of clusters
	 	for (var j = 0; j < clusters.length; ++j) {
	        // Detect empty clusters
			if (clusters[j].length === 0) {
				throw new Error('empty cluster at index: ' + j);
			}
			
			// Detect out of bounds asset indexes
			// Count the number of times each asset index appears in the list of clusters
			for (var k = 0; k < clusters[j].length; ++k) {
	            var assetIdx = clusters[j][k];
	            if (assetIdx < 1 || assetIdx > nbAssets) {
	                throw new Error('asset index out of bounds: ' + assetIdx);
	            }
	            else {
	                partition[assetIdx-1]++;
	            }
	        }
	    }
		    
		// Check that each integer in the set [1..nbAssets] is appears once and only once 
		// in the list of clusters
		for (var i = 0; i < partition.length; ++i) {
		    if (partition[i] !== 1) {
		        if (partition[i] === 0) {
		            throw new Error('missing asset index: ' + (i+1));
		        }
		        else if (partition[i] > 1) {
		            throw new Error('duplicate asset index: ' + (i+1));
		        }
		        else {
		            throw new Error('unknown error');
		        }
		    }
		}

		// All checks passed
	}
	// Otherwise, compute the clusters using the Fast Threshold Clustering Algorithm from David Varadi
	else if (clusteringMode === 'ftca') {
		// Extract the correlation matrix from the covariance matrix
		var corrMat = sigma.getCorrelationMatrix();
		
		// Compute the clusters using the FTCA algorithm
		clusters = ftca_(corrMat, opt.ftcaThreshold);
	}
	else {
		//
		throw new Error('unsupported clustering method');
	}
	

	// ------
	// The cluster risk parity portfolio is constructed in three steps:
	// 1 - Clusters creation (done above)
	// 2 - Assets allocation within each cluster using ERC weights
	// 3 - Portfolio allocation across all clusters using ERC weights and considering each cluster as a (synthetic) asset
	
	// 1 - N/A
	var nbClusters = clusters.length;
	
	// Instantiate the matrix mapping the initial assets to the clusters space (e.g., a change of base matrix),
	// to be populated with proper weights in step 2 below.
	var assetsToClustersWeights = Matrix_.zeros(nbClusters, nbAssets); 
	
	// 2 - For each cluster:
	for (var i = 0; i < nbClusters; ++i) {
		// Extract the covariance matrix associated to the assets it contains 
		var assetsIndexes = clusters[i].slice().sort(function(a, b) { return a - b; });
		
		var clusterSigma = sigma.getSubmatrix(assetsIndexes, assetsIndexes);
		
		// Compute ERC weights for these assets
		var assetsWeights = self.equalRiskContributionWeights(clusterSigma, opt);
		
		// - Populate the change of base matrix for the current cluster with the computed ERC weights
		for (var j = 0; j < assetsWeights.length; ++j) {
			assetsToClustersWeights.setValueAt(i+1, assetsIndexes[j], assetsWeights[j]);
		}
	}
	
	// 3a - Compute the the covariance matrix associated to the weighted clusters space, using formula Var(A*X) = A * Var(X) * A'.
	var clustersSigma = Matrix_.product(assetsToClustersWeights, Matrix_.product(sigma, assetsToClustersWeights.transpose()));
	
	// 3b - Compute ERC weights in the clusters space
	var clustersWeights = self.equalRiskContributionWeights(clustersSigma, opt);
	
	// 3c - Compute original assets weights, using formula A' * Y
	var weights = Matrix_.product(assetsToClustersWeights.transpose(), new Matrix_(clustersWeights));
	
	// Return them (already normalized)
	return weights.toArray();
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
* This portfolio maximizes the Sharpe ratio if the assets mean returns are proportional to their volatilities and all pair-wise correlations are equal.
* 
* @see <a href="https://ssrn.com/abstract=1949003">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* 
* @param {Matrix_|<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, column Matrix (i.e., vector) or array of n real numbers statisfying sigma[i-1] = sigma_i.
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
	
	// ------
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Matrix_(sigma).elemPower(-1/2).normalize();
	
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
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio for each stock is the same and all pair-wise correlations are equal.
*
* This portfolio is a special case of the more generic risk budgeting portfolio, with all risk budgets
* equal.
*
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/jpm.2010.36.4.060">Maillard, S., Roncalli, T., Teiletche, J.: The properties of equally weighted risk contribution portfolios. J. Portf. Manag. 36, 60–70 (2010)</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
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
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;	
	var rb = Matrix_.ones(nbAssets, 1).normalize();

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
* This portfolio maximizes the Sharpe ratio if the returns and volatilities are the same for all assets and all pair-wise correlations are equal.
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
	var weights = Matrix_.ones(nbAssets, 1).normalize();

	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function mostDiversifiedWeights
*
* @summary Compute the weights of the most diversified portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only
* most diversified portfolio of n assets, defined as the weights which maximizes the diversification ratio of the portfolio.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio of each asset is the same.
*
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/JPM.2008.35.1.40">Toward Maximum Diversification by Y. Choueifaty, Y. Coignard, The Journal of Portfolio Management, Fall 2008, Vol. 35, No. 1: pp. 40-51</a>
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the most diversified portfolio, array of n real numbers.
*
* @example
* mostDiversifiedWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // [0.33, 0.67]
*/
self.mostDiversifiedWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Extract the standard deviations of the assets in vector format
	var stddevs = sigma.getDiagonal().elemPower(1/2);

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1).normalize();
	
	// Preparational computations
	var sigma_x = Matrix_.product(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x)); // sigma(x)
	var dr = Matrix_.vectorDotProduct(stddevs, x) / s_x; // diversification ratio DR
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma
	var iter = 0;
	var converged = false;
	while (!converged) {
        // By default, let's assume the new iteration will lead to the convergence of the algorithm
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {
    	    // Define the coefficients of the second order polynomial ax_i^2 + b_ix + c, c.f. the second reference
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i) - stddevs.getValueAt(i,1); // (SIGMA*x)_i - x_i*sigma_i^2 - sigma_i, might be any sign
    	    var c = 0;
    	    
    	    // Note: what follows is not detailled in the second reference, as lambda_erc is supposed there to be non zero.
    	    //
    	    // The equation ax_i^2 + bx_i + c = 0 has two roots: 0 and another one, possibly strictly negative, 0 or strictly positive:
			// Case #1 - If the other root is 0 or strictly negative, x_i^* is then equal to 0
			// Case #2 - If the other root is strictly positive, x_i^* is equal to the value maximization the diversification ratio in 1D

			// Extract the "interesting" root of the equation ax_i^2 + bx_i = 0, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
    	    var q = -(b_p + sign_b_p * Math.abs(b_p));
    	    var r1 = q/a;
    	    
			// Case #1 always needs to be tested
			var xi_star_1 = 0;
			
			  // Case #1 weights update
			x.setValueAt(i, 1, xi_star_1);

    	      // Case #1 updated SIGMA*x and x'*SIGMA*x products
    	    sigma_x = Matrix_.product(sigma, x)
    	    s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x));
			
			  // Case #1 diversification ratio computation
			var dr_star = Matrix_.vectorDotProduct(stddevs, x) / s_x; // Weighted average volatility of the portfolio divided by the volatility of the portfolio
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				var sigma_x_2 = Matrix_.product(sigma, x)
				var s_x_2 = Math.sqrt(Matrix_.vectorDotProduct(sigma_x_2, x));
				
				// Case #2 diversification ratio computation
				var dr_star_2 = Matrix_.vectorDotProduct(stddevs, x) / s_x_2;
				
				// Selection between the two cases
				// Note: Portfolio sparsity is priviledged in case of equality
				if (dr_star_2 > dr_star) {
					// Case #2 weights update is already done
					
					// Update the products for convergence condition evaluation + next loop evaluation
					sigma_x = sigma_x_2;
					s_x = s_x_2;
					dr_star = dr_star_2;
				}
				else {
					// Case #1 weights re-update, as erased by case #2 test
					x.setValueAt(i, 1, xi_star_1);
					
					// No update on the products, as case #1 is the default case
				}
			}
    	    
    	    // Update the convergence condition: |dr - dr*| <= eps
    	    if (Math.abs(dr_star - dr) > eps) {
    	        converged = false;
    	        dr = dr_star;
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
 * @file Functions related to minimum correlation (heuristic) portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function minCorrWeights
*
* @summary Compute the weights of the minimum correlation (heuristic) portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only 
* minimum correlation (heuristic) portfolio of n assets, as computed by the minimum correlation algorithm of the reference.
*
* This portfolio is unique.
* 
* This portfolio is meant to approximate the most diversified portfolio.
*
* The algorithm used is the minimum correlation algorithm (MCA), with primary benefits versus the conventional optimization methods:
* - speed of computation and ease of calculation 
* - greater robustness to estimation error
* - superior risk dispersion
*
* @see <a href="www.systematicportfolio.com/mincorr_paper.pdf">David Varadi, Michael Kapler, Henry Bee, Corey Rittenhouse, Sep. 2012, The Minimum Correlation Algorithm: A Practical Diversification Tool</a>
*
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to theminimum correlation (heuristic) portfolio, array of n real numbers.
*
* @example
* minCorrWeights([[0.1, 0], [0, 0.2]]);
* // ~[0.59, 0.41]
*/
self.minCorrWeights = function (sigma, opt) {
	// Convert sigma to covariance matrix format
	var sigma = new Matrix_(sigma).toCovarianceMatrix();
	
	// Extract the correlation matrix and the misc. variances-related vectors
	var variances = sigma.getVariancesVector();
	var invStddevs = variances.elemPower(-1/2);
	var rho = sigma.getCorrelationMatrix();

	// TODO: Checks, if enabled	
	
	// ------
	var nbAssets = rho.nbRows;
	
	// Specific case to be filtered out: number of assets is 2 (or 1, but this case is less useful)
	if (nbAssets <= 2) {
		return self.equalRiskBudgetWeights(variances, opt);
	}
	
	// Step 2: Compute mean and sample standard deviation of the correlation matrix elements
	// Note: only strict upper triangular part is considered, c.f. the example of the reference
	var elemRho = rho.toArray(function(i, j, val) {
		return j > i;
	});
	var elementsMean = mean_(elemRho);
	var elementsStddev = sampleStddev_(elemRho);

	// Step 3: Create Adjusted Correlation Matrix
	var adjustedRho = rho.elemMap(function(i, j, val) { 
			if (i == j) {
				return 0; // Necessary for step 6
			}
			else {
				return 1 - normcdf_((val - elementsMean)/elementsStddev);
			}
		});

	// Step 4: Compute average value for each row (= column, as matrix is symetric) of the Adjusted Correlation Matrix
	// Note: only non diagonal elements are considered, c.f. the example of the reference
	var rowsElements = adjustedRho.toDoubleArray(function(i, j, val) {
		return i != j;
	});
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	
	// Step 5: Compute rank portfolio weight estimates
	// Note: ranks are computed in descending order, c.f. the example of the reference
	var ranks = rank_(rowsAverages, 0);
	var weights = new Matrix_(ranks, 1).normalize();
	
	// Step 6: Combine rank portfolio weight estimates with Adjusted Correlation Matrix
	weights = Matrix_.product(adjustedRho, weights).normalize();
	
	// Step 7: Scale portfolio weights by assets standard deviations
	weights = Matrix_.vectorHadamardProduct(weights, invStddevs).normalize();
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to proportional minimum variance (heuristic) portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 


/**
* @function minVarWeights
*
* @summary Compute the weights of the proportional minimum variance (heuristic) portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only 
* proportional minimum variance (heuristic) portfolio of n assets, as computed by the minimum variance algorithm of the reference.
*
* This portfolio is unique.
* 
* This portfolio is meant to approximate the global minimum variance portfolio.
*
* The algorithm used is the minimum variance algorithm (MVA), with primary benefits versus the conventional optimization methods:
* - speed of computation and ease of calculation 
* - greater robustness to estimation error
* - superior risk dispersion
*
* @see <a href="https://cssanalytics.wordpress.com/2013/04/05/minimum-variance-algorithm-mva-excel-sheet/">Minimum Variance Algorithm (MVA) Excel Sheet</a>
*
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to theminimum variance (heuristic) portfolio, array of n real numbers.
*
* @example
* minVarWeights([[0.1, 0], [0, 0.2]]);
* // ~[0.86, 0.14] 
*/
self.minVarWeights = function (sigma, opt) {
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// TODO: Checks, if enabled	
	
	// ------
	var nbAssets = sigma.nbRows;
	
	// Step 1: Average pairwise covariance, and associated mean/standard deviation
	var rowsElements = sigma.toDoubleArray();
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	var elementsMean = mean_(rowsAverages);
	var elementsStddev = sampleStddev_(rowsAverages);

	// Step 2: Gaussian convertion, and proportional average covar weigth
	var weights = new Matrix_(rowsAverages, 1).elemMap(function(i, j, val) { 
		return 1 - normcdf_((val - elementsMean)/elementsStddev);
	}).normalize();
	
	// Step 3: Scale portfolio weights by assets variances
	var invVariancesWeights = sigma.getDiagonal().elemPower(-1).normalize();
	weights = Matrix_.vectorHadamardProduct(weights, invVariancesWeights).normalize();
	
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
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {Array.<number>} rb the risk budgets, array of n real strictly positive numbers summing to one.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
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
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Convert rb to vector format
	var rb = new Matrix_(rb);

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that rb contains strictly positive numbers summing to one
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1).normalize();
	
	// Preparational computations
	var sigma_x = Matrix_.product(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x)); // sigma(x)
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma and b
	var iter = 0;
	var converged = false;
	while (!converged) {
        // Convergence condition is false if any of the coordinate-wise convergence condition is false
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {
    	    // Define the coefficients of the second order polynomial ax_i^2 + b_ix + c, c.f. the second reference
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i); // (SIGMA*x)_i - x_i*sigma_i^2, might be any sign
    	    var c = -rb.getValueAt(i,1) * s_x; // -b_i * sigma(x), always < 0
    	    
    	    // Extract the strictly positive root x_i^* of the equation ax_i^2 + bx_i + c = 0, using a stable numerical formula
    	    var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
			    throw new Error('Discriminant not positive during iteration ' + iter + ', covariance matrix might not be definite positive');
			}
    	    var q = -(b_p + sign_b_p * Math.sqrt(disc));
    	    var r1 = q/a;
    	    var r2 = c/q;
    	    var xi_star = r1 > 0 ? r1 : r2;
    	    
    	    // Update the weights
    	    x.setValueAt(i,1,xi_star);
    	    
    	    // Compute the updated SIGMA*x and x'*SIGMA*x products for convergence condition evaluation + next loop evaluation
    	    sigma_x = Matrix_.product(sigma, x)
    	    s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x));
    	    
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
