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

 
/* Start Wrapper private methods - Unit tests usage only */
self.Matrix = Matrix_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function Matrix_
*
* @summary Construct a matrix object from several possible inputs: an array of arrays of numbers, an array of numbers or a matrix object.
*
* @description This function constructs an n by m matrix (a_ij),i=1..n,j=1..m from either:
* - An array dblarr of n arrays of m real numbers, with coefficients a_ij satisfying a_ij = dblarr[i-1][j-1]
* - An array arr of n real numbers, with coefficients a_ij satisfying a_ij = arr[i-1] with m=1 (i.e., the constructed matrix is a column matrix, e.g. a vector)
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
		// Result matrix allocation
		that = allocateMatrix_(dblarr.length, dblarr[0].length);
	
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
		// Result matrix allocation
		that = allocateMatrix_(arr.length, 1);
		
		// Fill the vector
		for (var i = 0; i < that.nbRows; ++i) {
			that.data[i * that.nbColumns] = arr[i];
		}
		
		// Return
		return that;
	}
	
	function fromMatrix(mat) {
		// Both steps below are delegated to the matrix copy function:
		// - Result matrix allocation
		// - Matrix copy
		that = Matrix_.copy(mat);

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
	else if (input instanceof Matrix_) { // Equivalent to a "clone" operation on the Matrix
		return fromMatrix(input);
	}
	else {
		throw new Error('unsupported input type');
	}
}

/**
* @function allocateMatrix_
*
* @summary Returns a matrix to be used for subsequent computations.
*
* @description This internal function centralizes the allocation of a n by m matrix to be used for
* subsequent computations.
*
* @param {number} n the number of rows of the matrix to allocate, natural integer greater than or equal to 1.
* @param {number} m the number of columns of the matrix to allocate, natural integer greater than or equal to 1.
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} either a newly allocated matrix if out is not provided or out if out is provided, a n by m matrix.
*
* @example
* allocateMatrix_(2, 3);
*/
function allocateMatrix_(n, m, out) {
	// The logic of the matrix allocation is the following:
	// - If no output matrix is provided, a new n by m matrix is created
	// - If an output is provided, it is checked to be a matrix with proper dimensions, and if, it is re-used
	var obj;
	if (out === undefined) {
    	obj = Object.create(Matrix_.prototype);
    	obj.nbRows = n;
    	obj.nbColumns = m;
    	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	}
	else if (!(out instanceof Matrix_)){
		throw new Error('provided output must be a matrix');
	}
	else if (out.nbRows !== n || out.nbColumns !== m) {
		throw new Error('provided output matrix size does not match expected matrix size: ' + '(' + out.nbRows + ',' + out.nbColumns + 
		') v.s. ' + '(' + n + ',' + m + ')');
	}
	else {
		obj = out;
	}
	
	// Return the allocated matrix
	return obj;
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
	* @function toRowArray
	*
	* @summary Returns a double array containing all the elements of the matrix satisfying a predicate function,
	* organized by row.
	*
	* @description This function builds a double array arr from the matrix (a_ij),i=1..nbRows,j=1..nbColumns with arr[i-1] 
	* containing all the elements a_ij,j=1..nbColumns of the i-th row of the matrix satisfying fct(i, j, a_ij) == true
	* where fct is an optional predicate function.
	*
	* In case fct is not provided, it defaults to true, i.e., all the elements of the matrix are selected.
	* 
	* @memberof Matrix_
	* @param {function(number, number, number): number} fct the optional function to call on each element of the original matrix,
	* which should take 3 arguments: row index i=1..n, column index j=1..m, matrix element a_ij and which
	* should return a boolean, true in order to include the element a_ij in the output double array, false to skip it.
	* @return {Array.<Array.<number>>} an array of array of real numbers, representing the matrix' elements selected by the function fct.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).toRowArray();
	* //
	* [[1,2,3], [4,5,10]]
	*/
	toRowArray: function(fct) {
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
	* In case fct is not provided, it defaults to true, i.e., all the elements of the matrix are selected.
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
			throw new Error(
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
			throw new Error(
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
	* @description This function determines if the number of rows of the matrix
	* is equal to the number of columns of the matrix.
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
	* @summary Determines if the matrix is a vector.
	*
	* @description This function determines if the number of columns of the matrix is equal to 1.
	* 
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a vector, false otherwise.
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
	* @function isNonNegative
	*
	* @summary Determines if the matrix is a nonnegative matrix.
	*
	* @description This function determines if the matrix is a nonnegative matrix,
    * i.e. if all its elements are equal to or greater than zero.
	* 
	* @see <a href="https://en.wikipedia.org/wiki/Nonnegative_matrix">Nonnegative matrix</a>
	*
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a nonnegative matrix, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isNonNegative();
	* // true
	*
	* @example
	* Matrix_([[-1], [4]]).isNonNegative();
	* // false
	*/
	isNonNegative: function() {
		// Check the nonnegativity condition element by element
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] < 0) {
					return false;
				}
			}
		}
		
		// At this stage, the matrix is nonnegative.
		return true;
	},
	
	/**
	* @function isPositive
	*
	* @summary Determines if the matrix is a positive matrix.
	*
	* @description This function determines if the matrix is a positive matrix,
    * i.e. if all its elements are strictly greater than zero.
	* 
	* @see <a href="https://en.wikipedia.org/wiki/Nonnegative_matrix">Nonnegative matrix</a>
	*
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a positive matrix, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isPositive();
	* // true
	*
	* @example
	* Matrix_([[0], [4]]).isPositive();
	* // false
	*/
	isPositive: function() {
		// Check the positivity condition element by element
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] <= 0) {
					return false;
				}
			}
		}
		
		// At this stage, the matrix is positive.
		return true;
	},

	/**
	* @function isNonPositive
	*
	* @summary Determines if the matrix is a nonpositive matrix.
	*
	* @description This function determines if the matrix is a nonpositive matrix,
    * i.e. if all its elements are equal to or lower than zero.
	* 
	* @see <a href="http://mathworld.wolfram.com/NonpositiveMatrix.html"> Weisstein, Eric W. "Nonpositive Matrix." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/NonpositiveMatrix.html </a>
	*
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a nonpositive matrix, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isNonPositive();
	* // false
	*
	* @example
	* Matrix_([[-1], [0]]).isNonPositive();
	* // true
	*/
	isNonPositive: function() {	
		// Check the nonpositivity condition element by element
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] > 0) {
					return false;
				}
			}
		}
		
		// At this stage, the matrix is nonpositive.
		return true;
	},
	
	/**
	* @function isNegative
	*
	* @summary Determines if the matrix is a negative matrix.
	*
	* @description This function determines if the matrix is a negative matrix,
    * i.e. if all its elements are strictly lower than zero.
	* 
	* @see <a href="http://mathworld.wolfram.com/NegativeMatrix.html">Weisstein, Eric W. "Negative Matrix." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/NegativeMatrix.html </a>
	*
	* @memberof Matrix_
	* @return {boolean} true if the matrix is a negative matrix, false otherwise.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,10]]).isNegative();
	* // false
	*
	* @example
	* Matrix_([[-1], [-2]]).isNegative();
	* // true
	*/
	isNegative: function() {
		// Check the negativity condition element by element
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] >= 0) {
					return false;
				}
			}
		}
		
		// At this stage, the matrix is negative.
		return true;
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
		for (var k = 0; k < this.data.length; ++k) {
			sum += this.data[k];
		}
		
		// Return the computed value
		return sum;
	},
	
	/**
	* @function min
	*
	* @summary Returns the minimum element of the matrix.
	*
	* @description This function computes the minimum element of the matrix.
	* 
	* @memberof Matrix_
	* @return {number} the minimum element of the matrix.
	*
	* @example
	* Matrix_([[1,2,3]]).min();
	* // 1
	*/
	min: function() {
		// Computation of the minimum of a_ij, i=1..nbRows, j=1..nbColumns
		var minVal = Number.POSITIVE_INFINITY;
		for (var k = 0; k < this.data.length; ++k) {
		  if (this.data[k] < minVal) {
			minVal = this.data[k];
		  }
		}

		// Return the computed value
		return minVal;
	},
	
	/**
	* @function max
	*
	* @summary Returns the maximum element of the matrix.
	*
	* @description This function computes the maximum element of the matrix.
	* 
	* @memberof Matrix_
	* @return {number} the maximum element of the matrix.
	*
	* @example
	* Matrix_([[1,2,3]]).max();
	* // 3
	*/
	max: function() {
		// Computation of the maximum of a_ij, i=1..nbRows, j=1..nbColumns
		var maxVal = Number.NEGATIVE_INFINITY;
		for (var k = 0; k < this.data.length; ++k) {
		  if (this.data[k] > maxVal) {
			maxVal = this.data[k];
		  }
		}

		// Return the computed value
		return maxVal;
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
	* @param {Matrix_} out an optional n by m matrix.
	* @return {Matrix_} the normalized n by m matrix matrix, either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2,3]]).normalize();
	* // Matrix_([[1/6,1/3,1/2]])
	*/
	normalize: function(out) {
		// Result matrix allocation
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		
		// Computation of the sum of the matrix elements
		var sum = this.sum();
		if (sum === 0.0) {
			throw new Error('sum of coefficients of matrix is null');
		}
		
		// Normalization of the matrix
		for (var k = 0; k < this.data.length; ++k) {
			obj.data[k] = this.data[k] / sum; 
		}
		
		// Return the computed matrix
		return obj;
	},
	
	/**
	* @function diagonal
	*
	* @summary Returns the diagonal elements of a square matrix.
	*
	* @description This function returns, as a column matrix (i.e., vector) of n elements, 
	* the diagonal elements (a_ii), i=1..n from the original square matrix (a_ij),i=1..n,j=1..n.
	*
	* @memberof Matrix_
	* @param {Matrix_} out an optional n by 1 matrix.
	* @return {Matrix_} a n by 1 matrix (i.e., vector) containing the diagonal elements of the original matrix, 
	* either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).diagonal();
	* // Matrix_([1,5])
	*/
    diagonal: function (out) {
    	// Checks
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
    	
		// Result matrix allocation
		var obj = allocateMatrix_(this.nbRows, 1, out);
    	
    	// Extraction of diagonal elements of A
    	for (var i = 0; i < this.nbRows; ++i) {
    		obj.data[i] = this.data[i * (this.nbColumns + 1)] ;
    	}
    	
    	// Return the computed vector
        return obj;
    },
	
	/**
	* @function row
	*
	* @summary Returns the row elements of a matrix.
	*
	* @description This function returns, as a column matrix (i.e., vector) of m elements, 
	* the row elements (a_ij), j=1..m from the original matrix (a_ij),i=1..n,j=1..m.
	*
	* @memberof Matrix_
	* @param {number} i the row index of the matrix for wich to return the elements, a natural integer belonging to 1..number of matrix rows.
	* @param {Matrix_} out an optional n by 1 matrix.
	* @return {Matrix_} a n by 1 matrix (i.e., vector) containing the elements of the i-th row of the original matrix, 
	* either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).row(1);
	* // Matrix_([1,2])
	*/
    row: function (i, out) {
		// Bounds check
		if (i < 1 || i > this.nbRows) {
			throw new Error(
			'index out of bounds when getting matrix row, (' + i + 
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
    	
    	// Result matrix allocation
		var obj = allocateMatrix_(this.nbColumns, 1, out);
		
    	// Extraction of the elements of the i-th row of A
    	for (var j = 0; j < this.nbColumns; ++j) {
    		obj.data[j] = this.data[(i-1) * this.nbColumns + j] ;
    	}
    	
    	// Return the computed vector
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
	* @param {Matrix_} out an optional n by m matrix.
	* @return {Matrix_} a n by m matrix with the original matrix coefficients transformed by the function fct,
	* either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).elemMap(function(i,j,val) { return i+j;});
	* // Matrix_([[2,3,4], [3,4,5]])
	*/
	elemMap: function (fct, out) {
		// Result matrix allocation
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		
		// Computation of the fonction fct applied to the coefficients of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = fct(i+1, j+1, this.data[i * this.nbColumns + j]);
			}
		}
		
		// Return the computed matrix
		return obj;
	},
	
    /**
	* @function submatrix
	*
	* @summary Returns a (possibly non contiguous) submatrix from the original matrix, keeping the elements whose row and column indexes are specified.
	*
	* @description This function computes a matrix (c_ij),i=1..p,j=1..q from the original matrix (a_ij),i=1..n,j=1..m and from the lists of row/column indexes to keep
	* rindexes/cindexes, where p is equal to the length of rindexes and q is equal to the length of cindexes, with coefficients satisfying c_ij = a_rindexes[i]cindexes[j].
	*
	* @memberof Matrix_
	* @param {Array.<number>} rindexes the row indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..n.
    * @param {Array.<number>} cindexes the column indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..m.
	* @param {Matrix_} out an optional rindexes by cindexes matrix.
	* @return {Matrix_} a rindexes by cindexes matrix whose elements correspond to the elements of the original matrix
	* whose row/column indexes belong to the input lists of row/column indexes to keep, either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).submatrix([1], [2, 3];
	* // Matrix_([[2,3]])
	*/
    submatrix : function(rindexes, cindexes, out) {
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
    	
		// Result matrix allocation
		var obj = allocateMatrix_(rindexes.length, cindexes.length, out);
    	
    	// Computation of the elements of A
    	for (var i = 0; i < obj.nbRows; ++i) {
    		var rindex = rindexes[i] - 1; // row index of the original matrix
    		
    		for (var j = 0; j < obj.nbColumns; ++j) {
    		    var cindex = cindexes[j] - 1; // column index of the original matrix
    			obj.data[i * obj.nbColumns + j] = this.data[rindex * this.nbColumns + cindex];
    		}
    	}
    	
    	// Return the computed submatrix
        return obj;
    },
	
    /**
	* @function transpose
	*
	* @summary Returns the transpose of the original matrix.
	*
	* @description This function computes the transpose matrix (c_ij),i=1..n,j=1..m of the original matrix (a_ij),i=1..n,j=1..m,
	* with coefficients satisfying c_ij = a_ji.
	*
	* @memberof Matrix_
	* @param {Matrix_} out an optional m by n matrix.
	* @return {Matrix_} a m by n matrix, transposed of the original matrix, either stored in the matrix out or in a new matrix.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).transpose();
	* // Matrix_([[1,4], [2,5], [3,6]])
	*/
	transpose: function (out) {
		// Result matrix allocation
		var obj = allocateMatrix_(this.nbColumns, this.nbRows, out);
		
		// Computation of the transpose of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[j * this.nbColumns + i];
			}
		}
		
		// Return the computed transpose of A
		return obj;
	},
	
	/**
	* @function determinant
	*
	* @summary Returns the determinant of the square matrix.
	*
	* @description This function computes the determinant of the square matrix (a_ij),i=1..n,j=1..n,
	* using a QR decomposition.
	*
	* @memberof Matrix_
	* @return {number} the determinant of the matrix.
	*
	* @example
	* Matrix_([[1,2], [3,4]]).determinant();
	* // -2
	*/
	determinant: function () {
    	// Checks
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		
		// Compute the (Givens) QR decomposition of the matrix, using only R.
		var R = Matrix_.qrDecomposition(this, {qLess: true});
		
		// By property of the Givens QR decomposition, the determinant of the matrix
		// is then the product of the diagonal elements of the R matrix.
		var det = 1.0;
		for (var k = 0; k < R.data.length; k += R.nbColumns + 1) {
		        det *= R.data[k];
		}

		// Return the computed determinant
		return det;
	},
	
	/**
	* @function matrixNorm
	*
	* @summary Returns a matrix norm of the matrix.
	*
	* @description This function computes a matrix norm of the matrix, the exact norm depending on the value of the parameter p:
	* - 'one', for the matrix norm induced by the vector 1-norm (largest column sum of absolute elements)
	* - 'infinity', for the matrix norm induced by the vector infinity-norm (largest row sum of absolute elements)
	* - 'frobenius', for the matrix Frobenius norm (sqrt(trace(transpose(A) * A)))
	*
	* @see <a href="https://en.wikipedia.org/wiki/Matrix_norm">Matrix norm</a>
	*
	* @memberof Matrix_
	* @param {string} p the matrix norm to compute as detailled in the description, a string either equals to 'one', to 'infinity' or to 'frobenius'.
	* @return {number} the computed matrix norm.
	*
	* @example
	* Matrix_([[1,2], [3,4]]).matrixNorm('one');
	* // 7
	*/
	matrixNorm: function(p) {
		// The supported matrix norms are 1-norm, infinity-norm and Frobenius norm
		if (p == 'one') {
            // Compute the largest column sum of the absolute elements of the matrix
			var maxColSum = 0;
            for (var j = 0; j < this.nbColumns; ++j) {
                var colSum = 0;
			    for (var i = 0; i < this.nbRows; ++i) {
				    colSum += Math.abs(this.data[i * this.nbColumns + j]);
			     }
			     maxColSum = Math.max(maxColSum, colSum);
            }
            return maxColSum;
		}
		else if (p == 'infinity') {
			// Compute the largest row sum of the absolute elements of the matrix
		    var maxRowSum = 0;
		    for (var i = 0; i < this.nbRows; ++i) {
			    var rowSum = 0;
			    for (var j = 0; j < this.nbColumns; ++j) {
				    rowSum += Math.abs(this.data[i * this.nbColumns + j]);
			    }
			    maxRowSum = Math.max(maxRowSum, rowSum);
		    }
		    return maxRowSum;
		}
		else if (p == 'frobenius') {
			// Compute the Frobenius norm, which is equal to the l^2 vector norm of the matrix
			return this.vectorNorm('two');
		}
		else {
			throw new Error('unsupported matrix norm: ' + p);
		}
	},
	
	/**
	* @function vectorNorm
	*
	* @summary Returns a vector norm of the matrix or of a subset of the matrix.
	*
	* @description This function computes a vector norm of the matrix or of a subset of the matrix, 
	* the exact norm depending on the value of the parameter p:
	* - 'one', for the vector l^1 norm (a.k.a. Manhattan norm)
	* - 'two' for the vector l^2 norm (a.k.a. Euclidean norm)
	* - 'infinity' for the vector l^infinity norm (a.k.a. Maximum norm)
	*
	* The value of the optional parameter subset determines the subset of the matrix to consider:
	* - 'matrix', in order to compute a vector norm for the whole matrix.
	* - 'row', in order to compute a vector norm for a specific row of the matrix.
	* - 'column', in order to compute a vector norm for a specific column of the matrix.
	*
	* If the parameter subset is provided and is equal to 'row'/'column', the parameter idx corresponds
	* to the row/column index for which to compute the vector norm.
	*
	* @see <a href="https://en.wikipedia.org/wiki/Norm_(mathematics)">Norm (mathematics)</a>
	* @see Nicholas J. Higham. 2002. Accuracy and Stability of Numerical Algorithms (2nd ed.). Soc. for Industrial and Applied Math., Philadelphia, PA, USA. 
	*
	* @memberof Matrix_
	* @param {string} p the vector norm to compute, a string either equals to 'one', to 'two' or to 'infinity'.
	* @param {string} subset the subset of the matrix for which to compute the vector norm, an optional string either equals to 'matrix', 'row' or 'column'; defaults to 'matrix'.
	* @param {number} idx if subset is equal to 'row'/'column', the row/column index for which to compute the vector norm, a natural integer belonging to 1..nbRows/nbColumns.
	* @return {number} the computed vector norm.
	*
	* @example
	* Matrix_([[1,2], [3,4]]).vectorNorm('one');
	* // 10
	*
	* @example
	* Matrix_([[1,2], [3,4]]).vectorNorm('one', 'row', 1);
	* // 3
	*/
	vectorNorm: function(p, subset, idx) {
		// Initializations
		var idxStartRow;
		var idxEndRow;
		var idxStartColumn;
		var idxEndColumn;
		
		// Supported subsets are 'matrix', 'row', 'column'
		if (subset === undefined || subset == 'matrix') {
			// Indexes assignement
			idxStartRow = 0;
			idxEndRow = this.nbRows;
			idxStartColumn = 0;
			idxEndColumn = this.nbColumns;
		}
		else if (subset == 'row') {
			if (idx !== undefined) {
				// Checks
				if (idx < 1 || idx > this.nbRows) { 
					throw new Error('row index out of bounds: ' + idx);
				}

				// Indexes assignement
				idxStartRow = idx-1;
				idxEndRow = idx;
				idxStartColumn = 0;
				idxEndColumn = this.nbColumns;
			}
			else {
				throw new Error('undefined row index');
			}
		}
		else if (subset == 'column') {
			if (idx !== undefined) {
				// Checks
				if (idx < 1 || idx > this.nbColumns) { 
					throw new Error('column index out of bounds: ' + idx);
				}

				// Indexes assignement
				idxStartRow = 0;
				idxEndRow = this.nbRows;
				idxStartColumn = idx-1;
				idxEndColumn = idx;
			}
			else {
				throw new Error('undefined row index');
			}
		}
		else if (subset !== undefined) {
			throw new Error('unsupported matrix subset: ' + subset);
		}
		
		// The supported vector norms are l^1, l^2 and l^infinity norms
		if (p == 'one') {
            // Compute the sum of the absolute values of the elements of the matrix
			var absSum = 0;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					absSum += Math.abs(this.data[ij_idx]);
					ij_idx++;
				}
				i_m += this.nbColumns;
			}
			return absSum;
		}
		else if (p == 'two') {
			// Compute the l^2 vector norm using an accurate algorithm by S. J. Hammarling
			// C.f. problem 27.5 of the second reference
			var t = 0;
			var s = 1;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
				    var val = this.data[ij_idx];
					var absVal = Math.abs(val);
					if (absVal != 0) {
						if (absVal > t) {
							s = 1 + s * (t/val) * (t/val);
							t = absVal;
						}
						else  {
							s = s + (val/t) * (val/t);
						}
					}
					ij_idx++;
				}
				i_m += this.nbColumns;
			}
			return t * Math.sqrt(s);
		}
		else if (p == 'infinity') {
			// Compute the largest absolute value of the elements of the matrix
		    var maxAbsVal = 0;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
    		    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					 maxAbsVal = Math.max(maxAbsVal, Math.abs(this.data[ij_idx]));
					 ij_idx++;
				}
				i_m += this.nbColumns;
			}
		    return maxAbsVal;
		}
		else {
			throw new Error('unsupported vector norm: ' + p);
		}
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
	var nbElements = nbRows * nbColumns;
	var tol = eps || 0;
	for (var k = 0; k < nbElements; ++k) {
		if (Math.abs(a.data[k] - b.data[k]) > tol) {
			return false;
		}
	}
	
	//
	return true;
};

/**
* @function xpy
*
* @summary Returns the addition of a matrix with another matrix.
*
* @description This function computes X + Y, the addition of a n by m matrix X 
* with another n by m matrix Y.
*
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a n by m matrix.
* @param {Matrix_} out an optional n by m matrix, possibly either the matrix X or the matrix Y.
* @return {Matrix_} the matrix X + Y, either stored in the matrix out or in a new matrix, an n by m matrix.
*
* @example
* xpy(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[7,8,9], [10,11,12]]));
* // Matrix_([[8,10,12], [14,16,18]])
*/
Matrix_.xpy = function(X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	// Ensure X,Y are compatible
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	
	// Computation of X + Y
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = X.data[k] + Y.data[k];
	}
	
	// Return the computed matrix
    return obj;
};

/**
* @function xmy
*
* @summary Returns the difference of a matrix with another matrix.
*
* @description This function computes X - Y, the difference of a n by m matrix X 
* with another n by m matrix Y.
*
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a n by m matrix.
* @param {Matrix_} out an optional n by m matrix, possibly either the matrix X or the matrix Y.
* @return {Matrix_} the matrix X - Y, either stored in the matrix out or in a new matrix, an n by m matrix.
*
* @example
* xmy(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[7,8,9], [10,11,12]]));
* // Matrix_([[-6,-6,-6], [-6,-6,-6]])
*/
Matrix_.xmy = function(X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	// Ensure X,Y are compatible
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	
	// Computation of X - Y
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = X.data[k] - Y.data[k];
	}
	
	// Return the computed matrix
    return obj;
};

/**
* @function axpby
*
* @summary Returns the product of a matrix with a real number plus the product of another matrix with another real number.
*
* @description This function computes a*X + b*Y, the product of a n by m matrix X with a real number a
* added to the product of another n by m matrix Y with another real number b.
*
* @param {number} a a real number.
* @param {Matrix_} X a n by m matrix.
* @param {number} b a real number.
* @param {Matrix_} Y a n by m matrix.
* @param {Matrix_} out an optional n by m matrix, possibly either the matrix X or the matrix Y.
* @return {Matrix_} the matrix a*X + b*Y, either stored in the matrix out or in a new matrix, an n by m matrix.
*
* @example
* axpy(-1, Matrix_([[1,2,3], [4,5,6]]), 1, Matrix_([[7,8,9], [10,11,12]]));
* // Matrix_([[6,6,6], [6,6,6]])
*/
Matrix_.axpby = function(a, X, b, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	// Ensure X,Y are compatible
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	
	// Computation of aX + bY
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = a*X.data[k] + b*Y.data[k];
	}
	
	// Return the computed matrix
    return obj;
};


/**
* @function copy
*
* @summary Returns a copy of a matrix.
*
* @description This function builds a copy of the n by m matrix A.
*
* @param {Matrix_} A a n by m matrix.
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} a copy of the matrix A, either stored in the matrix out or in a new matrix, an n by m matrix
*
* @example
* copy(Matrix_([[1,2,3], [4,5,6]]));
* // Matrix_([[1,2,3], [4,5,6]])
*/
Matrix_.copy = function(A, out) {
	// Ensure A is a matrix
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}

	// Result matrix allocation
	var obj = allocateMatrix_(A.nbRows, A.nbColumns, out);

	// Copy of A into the provided out matrix
	var nbElements = A.nbRows * A.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = A.data[k];
	}
	
	// Return the computed matrix
	return obj;
};



/**
* @function elementwiseProduct
*
* @summary Returns the elementwise product of a matrix with another matrix.
*
* @description This function computes the elementwise product Z = X.*Y of a matrix X with another matrix Y,
* where A is a n by m matrix and B is either a n by m matrix (full matrix elementwise product), or
* a n by 1 matrix (row matrix elementwise product) or a 1 by m matrix (column matrix elementwise product).
*
* When used with a n by 1 matrix, this function mimics the behavior of a left multiplication of X with
* a diagonal matrix made of the n by 1 matrix elements: Z = Diag(Y) * X.
*
* When used with a 1 by m matrix, this function mimics the behavior of a right multiplication of X with
* a diagonal matrix made of the 1 by m matrix elements: Z = X * Diag(Y).
*
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a m by p matrix, a n by 1 matrix or a 1 by m matrix.
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} the matrix elementwise product A.*B, either stored in the matrix out or in a new matrix, a n by m matrix.
*
* @example
* elementwiseProduct(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,2,3]]));
* // Matrix_([[1,4,9], [4,10,18]])
*/
Matrix_.elementwiseProduct = function(X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	// Ensure X,Y are compatible
	if (!(X.nbRows === Y.nbRows && X.nbColumns === Y.nbColumns || 
		  Y.nbRows === X.nbRows && Y.nbColumns === 1 ||
		  Y.nbRows === 1 && Y.nbColumns === X.nbColumns)) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	
	// Computation of X.*Y
	var nbRows = X.nbRows;
	var nbColumns = X.nbColumns;
	if (X.nbRows === Y.nbRows && X.nbColumns === Y.nbColumns) { // full matrix elementwise product	
		var nbElements = nbRows * nbColumns;
		for (var k = 0; k < nbElements; ++k) {
			obj.data[k] = X.data[k] * Y.data[k];
		}
	}
	else if (Y.nbRows === X.nbRows && Y.nbColumns === 1) { // row matrix elementwise product
		var x_ij_idx = 0;
		for (var i = 0; i < nbRows; ++i) {
			var y_i = Y.data[i];
			for (var j = 0; j < nbColumns; ++j) {
				obj.data[x_ij_idx] = X.data[x_ij_idx] * y_i;
				x_ij_idx++;
			}
		}
	}
	else if (Y.nbRows === 1 && Y.nbColumns === X.nbColumns) { // column matrix elementwise product
		var x_ij_idx = 0;
		for (var i = 0; i < nbRows; ++i) {
			for (var j = 0; j < nbColumns; ++j) {
				obj.data[x_ij_idx] = X.data[x_ij_idx] * Y.data[j];
				x_ij_idx++
			}
		}
	}
	
	// Return the computed matrix
    return obj;
};




/**
* @function xy
*
* @summary Returns the product of a matrix with another matrix.
*
* @description This function computes X*Y, the product of a n by m matrix X with another m by p matrix Y.
*
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a m by p matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix X*Y, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* xy(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.xy = function(X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	// Checks
	if (X.nbColumns !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, Y.nbColumns, out);

	// Computation of X*Y product in IKJ format
	var n = X.nbRows;
	var m = X.nbColumns;
	var p = Y.nbColumns;
	
	var x_ik_idx = 0;
	var i_p = 0; // == i*p
	for (var i = 0; i < n; ++i) {
		var obj_ij_idx = i_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = 0; // obj(i,j) = 0
			obj_ij_idx++;
		}
	
		var y_kj_idx = 0;
		for (var k = 0; k < m; ++k) {
			var x_ik = X.data[x_ik_idx]; // x(i,k)
			
			obj_ij_idx = i_p;
			for (var j = 0; j < p; ++j) {
				obj.data[obj_ij_idx] += x_ik * Y.data[y_kj_idx]; // obj(i,j) += x(i,k) * y(k,j);
				y_kj_idx++;
				obj_ij_idx++;
			}
			x_ik_idx++;
		}
		i_p += p;
	}
	
	// Return the computed matrix
    return obj;
};


/**
* @function txy
*
* @summary Returns the product of the transpose of a matrix with another matrix.
*
* @description This function computes X^t*Y, the product of the transpose of 
* a m by n matrix X with another m by p matrix Y.
*
* @param {Matrix_} X a m by n matrix.
* @param {Matrix_} Y a m by p matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix X^t*Y, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* txy(Matrix_([[1,4], [2,5], [3,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.txy = function(X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	
	// Checks
	if (X.nbRows !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}

	// Result matrix allocation
	var obj = allocateMatrix_(X.nbColumns, Y.nbColumns, out);

	// Computation of X*Y product in KIJ format due to X being used in
	// transposed form
	var n = X.nbColumns;
	var m = X.nbRows;
	var p = Y.nbColumns;
	
	var x_ki_idx = 0;
	var obj_ij_idx = 0;
	
	// First k loop unrolled to initialize the matrix with zeros
	var k = 0;
	for (var i = 0; i < n; ++i) {
		var x_ki = X.data[x_ki_idx]; // x(k,i)
		var y_kj_idx = k_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = x_ki * Y.data[j]; // obj(i,j) = x(k,i) *y(k,j)
			obj_ij_idx++;
		}
		x_ki_idx++;
	}
	
	var k_p = p; // == k*p
    for (var k = 1; k < m; ++k) {
        obj_ij_idx = 0;
		for (var i = 0; i < n; ++i) {
            var x_ki = X.data[x_ki_idx];  // x(k,i)
			var y_kj_idx = k_p;
			for (var j = 0; j < p; ++j) {
                obj.data[obj_ij_idx] += x_ki * Y.data[y_kj_idx]; // obj(i,j) = x(k,i) *y(k,j)
				obj_ij_idx++;
				y_kj_idx++;
            }
			x_ki_idx++;
        }
		k_p += p;
    }

	// Return the computed matrix
    return obj;
};

/**
* @function axy
*
* @summary Returns the product of a matrix with a real number and with another matrix.
*
* @description This function computes a*X*Y, the product of a n by m matrix X with a real number a
* and with another m by p matrix Y.
*
* @param {number} a a real number.
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a m by p matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix a*X*Y, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* axy(Matrix_(1, [[1,2,3], [4,5,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.axy = function(a, X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	
	// Checks
	if (X.nbColumns !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, Y.nbColumns, out);

	// Computation of a*X*Y product in IKJ format
	var n = X.nbRows;
	var m = X.nbColumns;
	var p = Y.nbColumns;
	
	var x_ik_idx = 0;
	var i_p = 0; // == i*p
	for (var i = 0; i < n; ++i) {
		var obj_ij_idx = i_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = 0; // obj(i,j) = 0
			obj_ij_idx++;
		}
	
		var y_kj_idx = 0;
		for (var k = 0; k < m; ++k) {
			var x_ik = a * X.data[x_ik_idx]; // a * x(i,k)
			
			obj_ij_idx = i_p;
			for (var j = 0; j < p; ++j) {
				obj.data[obj_ij_idx] += x_ik * Y.data[y_kj_idx]; // obj(i,j) += x(i,k) * y(k,j);
				y_kj_idx++;
				obj_ij_idx++;
			}
			x_ik_idx++;
		}
		i_p += p;
	}
	
	// Return the computed matrix
    return obj;
};

/**
* @function axty
*
* @summary Returns the product of a matrix with a real number and with the
* transpose of another matrix.
*
* @description This function computes a*X*Y^t, the product of a n by m matrix X
* with a real number a and with the transpose of another p by m matrix Y.
*
* @param {number} a a real number.
* @param {Matrix_} X a n by m matrix.
* @param {Matrix_} Y a p by m matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix a*X*Y^t, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* axty(Matrix_(1, [[1,2,3], [4,5,6]]), Matrix_([[1,2,3], [1,2,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.axty = function(a, X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	
	// Checks
	if (X.nbColumns !== Y.nbColumns) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(X.nbRows, Y.nbRows, out);

	// Computation of a*X*Y product in IJK format due to Y being used in
	// transposed form
    for (var i = 0; i < X.nbRows; ++i) {
        for (var j = 0; j < Y.nbRows; ++j) {
            obj.data[i * obj.nbColumns + j] = 0;
            
            for (var k = 0; k < X.nbColumns; ++k) {
                obj.data[i * obj.nbColumns + j] += X.data[i * X.nbColumns + k] * Y.data[j * Y.nbColumns + k];
            }
            
            obj.data[i * obj.nbColumns + j] = a * obj.data[i * obj.nbColumns + j];
        }
    }
	
	// Return the computed matrix
    return obj;
};


/**
* @function atxy
*
* @summary Returns the product of the transpose of a matrix with a real number 
* and with another matrix.
*
* @description This function computes a*X^t*Y, the product of the transpose of 
* a m by n matrix X with a real number a and with another m by p matrix Y.
*
* @param {number} a a real number.
* @param {Matrix_} X a m by n matrix.
* @param {Matrix_} Y a m by p matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix a*X^t*Y, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* atxy(1, Matrix_([[1,4], [2,5], [3,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.atxy = function(a, X, Y, out) {
	// Ensure X, Y are matrices
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	
	// Checks
	if (X.nbRows !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}

	// Result matrix allocation
	var obj = allocateMatrix_(X.nbColumns, Y.nbColumns, out);

	// Computation of a*X*Y product in KIJ format due to X being used in
	// transposed form
    // First k loop unrolled to initialize the matrix with zeros
	var k = 0;
	for (var i = 0; i < X.nbColumns; ++i) {
		var a_ik = a * X.data[k * X.nbColumns + i];
		
		for (var j = 0; j < Y.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}

		for (var j = 0; j < Y.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] += a_ik * Y.data[k * Y.nbColumns + j];
		} 
	}
    for (var k = 1; k < X.nbRows; ++k) {
        for (var i = 0; i < X.nbColumns; ++i) {
            var a_ik = a * X.data[k * X.nbColumns + i];
			
            for (var j = 0; j < Y.nbColumns; ++j) {
                obj.data[i * obj.nbColumns + j] += a_ik * Y.data[k * Y.nbColumns + j];
            } 
        }
    }

	// Return the computed matrix
    return obj;
};

/**
* @function diagonal
*
* @summary Returns a diagonal matrix constructed from a vector.
*
* @description This function computes a diagonal square matrix (a_ij),i=1..n,j=1..n from 
* a vector vec of n elements, with coefficients a_ij satisfying a_ij = 0, i <> j and a_ij = vec_i1, i==j.
*
* @param {Matrix_} x a n by 1 matrix.
* @param {Matrix_} out an optional n by n matrix.
* @return {Matrix_} a n by n diagonal matrix with its diagonal elements corresponding to the elements of x,
* either stored in the matrix out or in a new matrix.
*
* @example
* diagonal(Matrix_([1,2,3]));
* // == Matrix_([[1,0,0], [0,2,0], [0,0,3]])
*/
Matrix_.diagonal = function(x, out) {
	// Checks
	if (!(x instanceof Matrix_) || !x.isVector()) {
		throw new Error('first input must be a vector');
	}

	// Result matrix allocation
	var n = x.nbRows;
	var obj = allocateMatrix_(n, n, out);
	
	// Result matrix computation
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < n; ++j) {
			obj.data[i * n + j] = 0;
		}
		obj.data[i * n + i] = x.data[i];
	}
	
	// Return the computed matrix
    return obj;
};


/**
* @function fillSymetric
*
* @summary Returns a symetric matrix constructed from a function.
*
* @description This function computes a square matrix (a_ij),i=1..n,j=1..n from a function, 
* with coefficients a_ij satisfying a_ij = a_ji = fct(i,j).
*
* The function fct is only called for indices i such that j >= i.
*
* @param {number} n, the order of the matrix to create, a strictly positive integer.
* @param {function(number, number): number} fct the function to call on each (i,j) pair of indexes with j >= i,
* which should take 2 arguments: row index i=1..n, column index j=i..n and which
* should return a number, which will be inserted into the matrix as its a_ij coefficient.
* @param {Matrix_} out an optional n by n matrix.
* @return {Matrix_} a n by n symetric matrix with its elements computed by the fonction fct,
* either stored in the matrix out or in a new matrix..
*
* @example
* fillSymetric(2, function(i,j) { return 0; });
* // == Matrix_([[0,0], [0,0]])
*/
Matrix_.fillSymetric = function(n, fct, out) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows and columns out of bounds: ' + n);
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(n, n, out);
	
	// Computation of the elements of the matrix
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < i; ++j) {
			obj.data[i * n + j] = obj.data[j * n + i];
		}
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * n + j] = fct(i + 1, j + 1);
		}	
	}
	
	// Return the computed matrix
    return obj;
};


/**
* @function fill
*
* @summary Returns a matrix constructed from a function.
*
* @description This function computes a n by m matrix (a_ij),i=1..n,j=1..m from a function, 
* with coefficients a_ij satisfying a_ij = fct(i,j).
*
* @param {number} n the number of rows of the matrix to construct, natural integer greater than or equal to 1.
* @param {number} m the number of columns of the matrix to construct, natural integer greater than or equal to 1.
* @param {function(number, number): number} fct the function to call on each (i,j) pair of indexes,
* which should take 2 arguments: row index i=1..n, column index j=i..n and which
* should return a number, which will be inserted into the matrix as its a_ij coefficient.
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} a n by m matrix with its elements computed by the fonction fct,
* either stored in the matrix out or in a new matrix.
*
* @example
* fill(2, 2, function(i,j) { return 0; });
* // == Matrix_([[0,0], [0,0]])
*/
Matrix_.fill = function(n, m, fct, out) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows out of bounds: ' + n);
	}
	if (m < 1) {
		throw new Error('input number of columns out of bounds: ' + m);
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(n, m, out);
	
	// Computation of the elements of the matrix
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < m; ++j) {
			obj.data[i * m + j] = fct(i + 1 , j + 1);
		}
	}
	
	// Return the computed matrix
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
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} the constructed matrix, either stored in the matrix out or in a new matrix.
*
* @example
* zeros(3, 2);
* // Matrix_([[0,0], [0,0], [0,0]])
*/
Matrix_.zeros = function(n, m, out) {
	// Result matrix allocation
	var obj = allocateMatrix_(n, m, out);
	
	// Result matrix computation
	var nbElements = n * m;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 0;
	}
	
	// Return the computed matrix
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
	// Result matrix allocation
	var obj = allocateMatrix_(n, m);
	
	// Result matrix computation
	var nbElements = n * m;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 1;
	}
	
	// Return the computed matrix
    return obj;
}


/**
* @function identity
*
* @summary Returns a matrix made of ones on the main diagonal and zeros elsewhere, i.e., an identity matrix.
*
* @description This function builds an n by n matrix (a_ij),i=1..n,j=1..n satisfying a_ii = 1, i=1..n and
* a_ij = 0, i=1..n,j=1..n,i<>j.
* 
* @param {number} n the row/column length of the matrix to construct, natural integer greater than or equal to 1.
* @return {Matrix_} the constructed matrix.
*
* @example
* identity(3);
* // Matrix_([[1,0,0], [0,1,0], [0,0,1]])
*/
Matrix_.identity = function(n) {
	// Result matrix allocation
	var obj = allocateMatrix_(n, n);
	
	// Result matrix computation
	var nbElements = n * n;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 0;
		if (k % (n + 1) == 0) {
			obj.data[k] = 1;
		}
	}
	
	// Return the computed matrix
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
	// Delegate the computations to the internal function.
	return Matrix_.elementwiseProduct(x, y);
}


/**
* @function vectorDotProduct
*
* @summary Returns the dot product of two column matrices (e.g. vectors).
*
* @description This function computes the dot product <x/y>, 
* where x and y are two n by 1 column matrices (e.g. vectors).
* 
* @param {Matrix_} x a n by 1 column matrix.
* @param {Matrix_} y a n by 1 column matrix.
* @return {number} the dot product <x/y>, a real number.
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
	var nbElements = x.nbRows;
	for (var i = 0; i < nbElements; ++i) {
		dotProd += x.data[i] * y.data[i]; 
	}
	
	// Return it
	return dotProd;
}


/**
* @function qrDecomposition
*
* @summary Returns a QR decomposition of a matrix, using Givens rotations.
*
* @description This function computes a QR decomposition of an m by n matrix A with m >= n, 
* using Givens rotations as described in the algorithm 5.2.4 (and above discussion therein) of the first reference.
* 
* To be noted that the 'givens' internal function used in this function is not the same
* as the one described in the first reference, but is the continuous one described
* in the algorithm 4 of the second reference.
* 
* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
* @see <a href="http://www.netlib.org/lapack/lawnspdf/lawn150.pdf">Anderson, Edward (4 December 2000). Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem. LAPACK Working Note. University of Tennessee at Knoxville and Oak Ridge National Laboratory.</a>
*
* @param {Matrix_} A an m by n matrix, with m >= n.
* @param {object} opt optional parameters for the algorithm.
* @param {boolean} opt.qLess boolean parameter to be set to true to discard the computation of the Q matrix; defaults to false.
* @return {<Array.<Matrix_>|Matrix_} either an R matrix if opt.qLess is set to true or an array of two matrices [Q, R] otherwise, 
* with Q and R satisfying the following properties:
* - Q is an orthogonal m by m matrix, with a determinant equals to 1 (i.e., a rotation)
* - R is an m by n upper triangular matrix, with its bottom (m−n) rows consisting entirely of zeroes
* - A = Q*R
*
* @example
* qrDecomposition(Matrix_([[2,3], [0,4]]));
* // [Matrix_([[1,0], [0,1]]), Matrix_([[2,3], [0,4]])]
*/
Matrix_.qrDecomposition = function(A, opt) {
	/**
    * @function givens
    *
    * @description Given real numbers a and b, this function computes c = cos(theta), s = sin(theta) and r >= 0 
    * satisfying [[c, s], [-s, c]]^T * [[a], [b]] = [[r], [0]].
    * 
    * C.f. the second refence for details.
    * 
    * @param {number} a, a real number
    * @param {number} b, a real number
    * @return {<Array.<number>} an array of three real numbers - c at array index 0, s at array index 1 and r at array index 2 -, satisfying the matrix equation above.
    */
	function givens(a, b) {
		var c;
		var s;
		var r;
		
		if (b == 0) {
			c = (a >= 0) ? 1 : -1; // emulates sign(a)
			s = 0;
			r = Math.abs(a);
		}
		else if (a == 0) {
			c = 0;
			s = ((b >= 0) ? 1 : -1); // emulates sign(b)
			r = Math.abs(b);
		}
		else if (Math.abs(a) > Math.abs(b)) {
			var t = b/a;
			var u = ((a >= 0) ? 1 : -1) * Math.sqrt(1 + t*t);
			c = 1/u;
			s = t*c;
			r = a*u;
		}
		else {
			var t = a/b;
			var u = ((b >= 0) ? 1 : -1) * Math.sqrt(1 + t*t);
			s = 1/u;
			c = t*s;
			r = b*u;
		}
		
		return [c, -s, r]; // Compared to the second reference, the sign of s is altered so that it is made compatible with the notations of the first reference.
	}

	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var isQLessDecomposition = opt.qLess || false;
	
	// Checks
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (A.nbRows < A.nbColumns) {
		throw new Error('matrix has more columns than rows: ' + '(' + A.nbRows + ') v.s. ' + '(' + A.nbColumns + ')');
	}

	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;

	// Create a copy of A so that it is not overwritten
	var R = new Matrix_(A); // represents R
	
	// Create the matrix that will hold Q
	var Q = null;
	if (!isQLessDecomposition) {
		Q = Matrix_.identity(m); // represents Q
	}
	
	// Core of the algorithm
	for (var j = 1; j <= n; ++j) {
		for (var i = m; i >= j+1; --i) {
			// This loop iteration will set R(i,j) to 0 using an appropriate Givens rotation
			
			// Builds G(i-1, i, theta)	
			var R_im_j_idx = (i-2) * R.nbColumns + (j-1);
			var R_i_j_idx = (i-1) * R.nbColumns + (j-1);
			
			var a = R.data[R_im_j_idx]; // R(i-1, j)
			var b = R.data[R_i_j_idx]; // R(i, j)
			var givensArr = givens(a, b);
			var c = givensArr[0];
			var s = givensArr[1];
			var r = givensArr[2];
			
			// Update R (left multiply R with the transpose of the Givens rotation, c.f. 5.1.9 of the first reference)		
				// Loop below unrolled for k = j, since R(i,j) must be made zero
			R.data[R_im_j_idx] = r; // R(i-1, j) = r
			R.data[R_i_j_idx] = 0; // R(i, j) = 0

				// Loop resumed at k = j+1
			for (var k = j+1; k <= n; ++k) {
				var R_im_k_idx = (i-2) * R.nbColumns + (k-1);
				var R_i_k_idx = (i-1) * R.nbColumns + (k-1);
			
				var t1 = R.data[R_im_k_idx]; // t1 = R(i-1, k)
				var t2 = R.data[R_i_k_idx]; // t2 = R(i, k)
				R.data[R_im_k_idx] = c*t1 - s*t2; // R(i-1, k) = ...
				R.data[R_i_k_idx] = s*t1 + c*t2; // R(i, k) = ...
			}
			
			// Update Q (right multiply Q with the Givens rotation, c.f. 5.1.9 of the first reference)
			if (!isQLessDecomposition) {
				for (var k = 1; k <= m; ++k) {
					var Q_k_im_idx = (k-1) * Q.nbColumns + (i-2);
					var Q_k_i_idx = (k-1) * Q.nbColumns + (i-1);
				
					var t1 = Q.data[Q_k_im_idx] // t1 = Q(k,i-1)
					var t2 = Q.data[Q_k_i_idx] // t2 = Q(k,i)
					Q.data[Q_k_im_idx] = c*t1 - s*t2; // Q(k,i-1) = ...
					Q.data[Q_k_i_idx] = s*t1 + c*t2; // Q(k,i) = ...
				}
			}
		}
	}
	
	// Return either the computed [Q, R] pair or the matrix R
	if (!isQLessDecomposition) {
		return [Q, R];
	}
	else {
		return R;
	}
}

/**
* @function svdDecomposition
*
* @summary Returns a singular value decomposition of a matrix, using a one-sided Jacobi algorithm.
*
* @description This function computes a singular value decomposition of an m by n matrix A with m >= n, 
* using the one-sided Jacobi algorithm described in the algorithm 4.1 of the first reference, 
* together with some computational improvements described in the second reference.
* 
* The computed decomposition can either be the thin or the full singular value decomposition,
* c.f. the third reference for the associated definitions.
*
* To be noted that a singular value decomposition of a rank-deficient matrix might not produce exact zero singular values 
* due to finite numerical precision.
* 
* @see <a href="https://doi.org/10.1137/0613074">James Demmel and Kresimir Veselic, Jacobi’s Method is More Accurate than QR, SIAM Journal on Matrix Analysis and Applications, 1992, Vol. 13, No. 4 : pp. 1204-1245</a>
* @see <a href="https://epubs.siam.org/doi/abs/10.1137/0910023">P. P. M. de Rijk, A One-Sided Jacobi Algorithm for Computing the Singular Value Decomposition on a Vector Computer, SIAM Journal on Scientific and Statistical Computing, 1989, Vol. 10, No. 2 : pp. 359-371</a>
* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
*
* @param {Matrix_} A an m by n matrix, with m >= n.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-16.
* @param {number} opt.maxIter maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100.
* @param {string} opt.svdForm string either equal to 'full' to compute the full singular value decomposition of the matrix A, or equal to 'thin' to compute
* the thin singular value decomposition of the matrix A; defaults to 'thin'.
* @return {<Array.<Matrix_>} an array of three matrices [U, S, V], with U, S and V satisfying the following properties:
* - If opt.svdForm is equal to 'thin':
* -- U is an m by n orthogonal matrix
* -- S is an n by n diagonal matrix, with its diagonal elements S_ii, i=1..n satisfying S_11 >= ... >= S_nn >= 0 
* - If opt.svdForm is equal to 'full':
* -- U is an m by m orthogonal matrix
* -- S is an m by n matrix, with its elements S_ii, i=1..n satisfying S_11 >= ... >= S_nn >= 0 and S_ij = 0, i=1..m, j=1..n, j <> i
* - In all cases:
* - V is an n by n orthogonal matrix 
* - A = U*S*V^t
*
* @example
* svdDecomposition(Matrix_([[1,0], [0,1]]));
* // [Matrix_([[1,0], [0,1]]), Matrix_([[1,0], [0,1]]), Matrix_([[1,0], [0,1]])]
*/
Matrix_.svdDecomposition = function(A, opt) {
	// ------
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-16;
	var maxIterations = opt.maxIter || 100;
	var svdForm = opt.svdForm || 'thin';
	
	
	// ------
	
	// Checks
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (A.nbRows < A.nbColumns) {
		throw new Error('matrix has more columns than rows: ' + A.nbColumns + ' v.s. ' + A.nbRows);
	}
	
	
	// ------
	
	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	
	// Create a copy of A so that it is not overwritten
	var uu = new Matrix_(A); // represents U 
	var u_frob_norm = uu.matrixNorm('frobenius'); 	
	
	// Create the matrix that will hold V
	var vv = Matrix_.identity(n); // represents V
	
	
	// ------
	
	// Core of the algorithm, guaranteed to converge per theorem 4.9 of the first reference
	var iter = 0;
	var u_columns_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	while (true) {
		// Update the number of iterations
		++iter;

		// Check the number of iterations (number of sweeps)
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		
		// Compute/update the squares of the euclidean norms of the columns
		// of U in order to avoid the accumulation of numerical 
		// round-off errors in the Jacobi sweeps
		//
		// To be noted that as the Frobenius norm of U is not
		// altered by Jacobi rotations, no update of this norm is necessary
    	for (var j = 1; j <= n; ++j) {
    		var u_j_two_norm = uu.vectorNorm('two', 'column', j);
			u_columns_two_norm_sq[j-1] = u_j_two_norm * u_j_two_norm;
    	}
		
		// For all pairs (i, j) with i < j, proceed with one Jacobi sweep,
		// that is, n*(n-1)/2 Jacobi rotations
		var converged = true; // the convergence status is updated during the sweep
		for (var j = 2; j <= n; ++j) {
			for (var i = 1; i <= j-1; ++i) {
				// Compute the (i, j) 2x2 submatrix [[aa, cc], [cc, bb]] of U^t*U
				var aa = u_columns_two_norm_sq[i-1]; // aa = sum U(k,i)^2, k=1..m = ||U(:,i)||_2^2
				var bb = u_columns_two_norm_sq[j-1]; // bb = sum U(k,j)^2, k=1..m = ||U(:,j)||_2^2
				var cc = 0; // cc = sum U(k,i)*U(k,j), k=1..m = <U(:,i)/U(:.j)>
			    for (var k = 1; k <= m; ++k) {
					cc += uu.data[(k-1) * uu.nbColumns + (i-1)] * uu.data[(k-1) * uu.nbColumns + (j-1)]; 
				}

				// Test on convergence conditions before applying a Jacobi rotation
				// on columns i and j of U, c.f. formula 4.2.1 of the second reference:
				// - |<U(:,i)/U(:.j)>| > Math.sqrt(m) * eps * ||U(:,i)||_2 * ||U(:.j)||_2
				// - ||U(:,i)||_2 > eps * ||U||_f
				// - ||U(:,j)||_2 > eps * ||U||_f
				if (Math.abs(cc) <= eps * Math.sqrt(m * aa * bb)) { // this condition also covers the case cc = 0, which would make zeta below indefinite
					continue;
				}
				if (Math.sqrt(aa) <= eps * u_frob_norm || Math.sqrt(bb) <= eps * u_frob_norm) {
				    continue;
				}
				
				// The convergence conditions are not satisfied yet, so, a Jacobi
				// rotation is needed
				converged = false;
				
				// Compute the Jacobi rotation which diagonalizes the 
				// 2x2 submatrix [[aa, cc], [cc, bb]]
				var zeta = (aa - bb)/(2 * cc);
				var t = ((zeta >= 0) ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta)); // first part emulates sign(zeta)
				var cs = 1 / Math.sqrt(1 + t*t);
				var sn = cs * t;

				// Update columns i and j of U (right multiply U with the Jacobi rotation)
				for (var k = 1; k <= m; ++k) {
					var t1 = uu.data[(k-1) * uu.nbColumns + (i-1)] // t1 = U(k,i)
					var t2 = uu.data[(k-1) * uu.nbColumns + (j-1)] // t2 = U(k,j)
					uu.data[(k-1) * uu.nbColumns + (i-1)] = cs*t1 + sn*t2 // U(k,i) = ...
					uu.data[(k-1) * uu.nbColumns + (j-1)] = -sn*t1 + cs*t2 // U(k,j) = ...
				}
				
				// Update the squares of the euclidean norms of columns i and j of U,
                // c.f. formula 2.2.10 of the second reference
                //
                // To be noted that this update creates numerical round-off errors, mitigated
                // by the full recomputation of these norms at the beginning of each Jacobi sweep
                //
                // To also be noted that as the Frobenius norm of U is not altered
                // by the Jacobi rotations, no update of this norm is necessary
                u_columns_two_norm_sq[i-1] += t * cc;
                u_columns_two_norm_sq[j-1] -= t * cc;

                // Update columns i and j of V (right multiply V with the Jacobi rotation)
				for (var k = 1; k <= n; ++k) {
					var t1 = vv.data[(k-1) * vv.nbColumns + (i-1)] // t1 = V(k,i)
					var t2 = vv.data[(k-1) * vv.nbColumns + (j-1)] // t2 = V(k,j)
					vv.data[(k-1) * vv.nbColumns + (i-1)] = cs*t1 + sn*t2 // V(k,i) = ...
					vv.data[(k-1) * vv.nbColumns + (j-1)] = -sn*t1 + cs*t2 // V(k,j) = ...
                }
            }
		}
		
		// In case the convergence status is true at the end of the sweep, 
		// the algorithm can be stopped
		if (converged == true) {
			break;
		}
	}

	
	// ------
	
	// At this stage:
	// - U is a rectangular m by n matrix with its non null columns made of orthogonal vectors
	// - V is an orthogonal square n by n matrix
	// - The right singular vectors are the columns of V
	// - The singular values are the norms of the columns of U
	// - The left singular vectors are the normalized columns of U
	//
	// Additional work is needed in order to achieve the thin singular value decomposition of A with:
    // - A = U*S*V^t
	// - U is an orthogonal rectangular m by n matrix
	// - S a diagonal n by n matrix with the n singular values of A on its diagonal 
	//   verifying S(i,i) = sigma_i and sigma_1 >= ... >= sigma_n >= 0
	// - V is an orthogonal square n by n matrix
	
	// Compute and sort the singular values of A in increasing order 
	// together with their V matrix column indexes
	var sigmas = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	var sigmas_idx = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
	for (var j = 1; j <= n; ++j) {
		// Compute sigma_j
		var singVal_j = Math.sqrt(u_columns_two_norm_sq[j-1]); // == uu.vectorNorm('two', 'column', j);
		sigmas[j-1] = singVal_j; 
		
		// TODO: Replace "numerically zero" singular values with 0
		
		sigmas_idx[j-1] = j;
	}
	sigmas_idx.sort(function(a, b) { return sigmas[b-1] - sigmas[a-1]; });
	
	// Compute the thin U,S and V matrices
	var uuu = Matrix_.zeros(m, n); 
	var sss = Matrix_.zeros(n, n);
	var vvv = Matrix_.zeros(n, n);
	for (var j = 1; j <= n; ++j) {
		// Extract singluar values information
		var sigma_j_col_idx = sigmas_idx[j-1];
		var sigma_j = sigmas[sigma_j_col_idx-1];

		// Thin S(j,j) must be equal to sigma_j
		sss.data[(j-1) * sss.nbColumns + (j-1)] = sigma_j;
		
		// Thin U(:,j) must be equal to the normalized column of the intermediate U
		// corresponding to sigma_j
		if (sigma_j != 0) { 
			for (var i = 1; i <= m; ++i) {
				uuu.data[(i-1) * uuu.nbColumns + (j-1)] = uu.data[(i-1) * uu.nbColumns + (sigma_j_col_idx-1)] / sigma_j;
			}
		}
		
		// Thin V(:,j) must be equal to the column of the intermediate V corresponding to sigma_j
		for (var i = 1; i <= n; ++i) {
			vvv.data[(i-1) * vvv.nbColumns + (j-1)] = vv.data[(i-1) * vv.nbColumns + (sigma_j_col_idx-1)];
		}		
	}
		
	// If the thin SVD decomposition of A is requested, return the 
	// computed U,S,V triple
	if (svdForm == 'thin') {
		return [uuu, sss, vvv];
	}
	
	// ------
	
	// Additional work is needed in order to achieve the full singular value decomposition of A with:
    // - A = U*S*V^t
	// - U is an orthogonal square m by m matrix
	// - S a rectangular m by n matrix with the n singular values of A on its diagonal 
	//   verifying S(i,i) = sigma_i and sigma_1 >= ... >= sigma_n >= 0
	// - V is an orthogonal square n by n matrix
	
	// Full U computation is made through a QR decomposition of thin U and completing the missing
	// columns of thin U by those of Q
	//
	// Full S computation is made by extending the thin S with zeroes
	var qr = Matrix_.qrDecomposition(uuu);
	var q = qr[0];
	
	var uuuu = Matrix_.zeros(m, m);
	var ssss = Matrix_.zeros(m, n);
	for (var j = 1; j <= n; ++j) {
		// Extract singluar values information
		var sigma_j_col_idx = sigmas_idx[j-1];
		var sigma_j = sigmas[sigma_j_col_idx-1];
		
		// Full S(j,j) must be equal to sigma_j
		ssss.data[(j-1) * ssss.nbColumns + (j-1)] = sigma_j;
		
		// Full U(:,j) must be equal to the thin U(:,j) if the associated
		// singular value sigma_j is not null, and otherwise to Q(:,j)
		if (sigma_j != 0) {
			for (var i = 1; i <= m; ++i) {
				uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = uuu.data[(i-1) * uuu.nbColumns + (j-1)];
			}
		}
		else {
			for (var i = 1; i <= m; ++i) {
				uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = q.data[(i-1) * q.nbColumns + (j-1)];
			}
		}
	}
	for (var j = n+1; j <= m; ++j) {
		// Full U(:,j) must be equal to Q(:,j)
		for (var i = 1; i <= m; ++i) {
			uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = q.data[(i-1) * q.nbColumns + (j-1)];
		}
	}
	
	// If the full SVD decomposition of A is requested, return the 
	// computed U,S,V triple
	if (svdForm == 'full') {
		return [uuuu, ssss, vvv];
	}
}


/**
* @function nullSpace
*
* @summary Returns an orthonormal basis of the null space of a matrix.
*
* @description This function computes an orthonormal basis of the null space of an m by n matrix A, 
* obtained from a singular value decomposition of this matrix.
* 
* As a singular value decomposition of a rank-deficient matrix might not produce exact zero singular values
* due to finite numerical precision, the dimension of the null space of the matrix A is numerically defined as
* the number of singular values strictly lower than an optional tolerance eps.
*
* To be noted that in case the matrix A is of full-rank, a zero vector is returned.
* 
* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
*
* @param {Matrix_} A an m by n matrix.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps numerical tolerance under which a singular value is considered as zero, a strictly positive real number; defaults to max(m,n) * eps(s_1), where eps(s_1) is the distance from the greatest singular value s_1 of A to the next larger floating-point number.
* @return {Matrix_} either an m by 1 matrix made of zeroes in case the matrix A is of full-rank, or an m by p orthogonal matrix
* where p is equal to the numerical dimension of the null space of the matrix A.
*
* @example
* nullSpace(Matrix_([[2,1],[-4,-2]]));
* // Matrix_([-0.4472135954999578, 0.8944271909999159])
*/
Matrix_.nullSpace = function(A, opt) {
	// ------
	
	// Internal function used to automatically adapt the precision
	// of the numerical rank computation
	function nextUp(x) {
		var EPSILON = Math.pow(2, -52)
		var MAX_VALUE = Number.MAX_VALUE
		var MIN_VALUE = Math.pow(2, -1022)

		if (x !== x) return x
		if (x === -1 / 0) return -MAX_VALUE
		if (x === 1 / 0) return +1 / 0
		if (x === MAX_VALUE) return +1 / 0
		var y = x * (x < 0 ? 1 - EPSILON / 2 : 1 + EPSILON)
		if (y === x) y = MIN_VALUE * EPSILON > 0 ? x + MIN_VALUE * EPSILON : x + MIN_VALUE
		if (y === +1 / 0) y = +MAX_VALUE
		var b = x + (y - x) / 2
		if (x < b && b < y) y = b
		var c = (y + x) / 2
		if (x < c && c < y) y = c
		return y === 0 ? -0 : y
	}
	
    // ------
   
    // Decode options
    if (opt === undefined) {
        opt = {};
    }
    var eps = opt.eps || undefined;

   
    // ------
   
    // Misc. checks
    if (!(A instanceof Matrix_)) {
        throw new Error('first input must be a matrix');
    }   


    // ------
   
    // Initializations
    var m = A.nbRows;
    var n = A.nbColumns;
    var a_ns = null;
   
    // ------
   
    // Depending on the shape of the matrix A, a different algorithm is used:
    // - m >= n: an orthogonal basis of Ker(A) is directly read from a thin SVD
    // decomposition of the matrix A, c.f. corollary 2.4.6 of the reference
    //
    // - m < n: an orthogonal basis of Ker(A) is obtained by computing the orthogonal
    // complement of an orthogonal basis of Range(A^t), this orthogonal basis being directly
    // read from a thin SVD decomposition of A^t, c.f. corollary 2.4.6 of the reference
    //
    // The underlying properties which justifies this approach is that
    // R^m = Range(A^t) _|_ Ker(A) for any m by n real matrix A
    if (m >= n) {
        // Compute a thin SVD decomposition of A: A = U*S*V^t
		var svd = Matrix_.svdDecomposition(A, {maxIter: -1});		
        var u = svd[0];
        var s = svd[1];
        var v = svd[2];
       
		// Compute the default tolerance, if required
		if (eps === undefined) {
			eps = m * (nextUp(s.data[0]) - s.data[0]);
		}
	   
        // Determine the numerical rank of A
        var r;
        for (r = n; r >= 1; --r) {
            // Comparison of the r-th greatest singular value of A with the tolerance parameter to determine
            // the first "zero" singular value of A
            if (s.data[(r-1)*s.nbColumns + (r-1)] > eps) {
                break;
            }
        }
        var p = r;
       
        // If the matrix A is of full rank n, its null space is reduced to the zero vector
        //
        // Otherwise, extract from the matrix V the n-p right singular vectors associated to the
        // n-p zero singular values, which then form an orthogonal basis of Ker(A)
        if (p == n) {
            a_ns = Matrix_.zeros(n, 1);
        }
        else {
            // Extract from V an orthogonal basis of Ker(A) of dimension n-p
            a_ns = new Matrix_.zeros(n, n - p);
            for (var j = p + 1; j <= n; ++j) {
                for (var i = 1; i <= n; ++i) {
                    a_ns.data[(i-1)*a_ns.nbColumns + ((j-(p+1)+1)-1)] = v.data[(i-1)*v.nbColumns + (j-1)];
                }
            }
        }       
    }
    else {
        // Compute a thin SVD decomposition of A^t: A^t = U*S*V^t
        var svd = Matrix_.svdDecomposition(A.transpose(), {maxIter: -1, svdForm: 'full'});
        var u = svd[0];
        var s = svd[1];
        var v = svd[2];
       
		// Compute the default tolerance, if required
		if (eps === undefined) {
			eps = n * (nextUp(s.data[0]) - s.data[0]);
		}
	   
        // Determine the numerical rank of A^t
        var r;
        for (r = m; r >= 1; --r) {
            // Comparison of the r-th greatest singular value of A^t with the tolerance parameter to determine
            // the first "zero" singular value of A^t
            if (s.data[(r-1)*s.nbColumns + (r-1)] > eps) {
                break;
            }
        }
        var p = r;
       
        // If the matrix A^t is of full rank, the null space of the matrix A
        // is reduced to the zero vector
        //
        // Otherwise:
        // - Extract from the matrix U the p left singular vectors associated to the
        // p non zero singular values, which then form an orthogonal basis of Range(A^t)
        //
        // - Compute the orthogonal complement of this orthogonal basis, which then
        // form an orthogonal basis of Ker(A)
        if (p == m) {
            a_ns = Matrix_.zeros(n, 1);
        }
        else {
            // Extract from U an orthogonal basis of Range(A^t) of dimension p,
            // with n > m > p
            var ta_rg = new Matrix_.zeros(n, p);
            for (var j = 1; j <= p; ++j) {
                for (var i = 1; i <= n; ++i) {
                    ta_rg.data[(i-1)*ta_rg.nbColumns + (j-1)] = u.data[(i-1)*u.nbColumns + (j-1)];
                }
            }
           
            // Compute the orthogonal complement, of dimension n-p,
            // of the above orthogonal basis
            var qr = Matrix_.qrDecomposition(ta_rg);
            var q = qr[0];
           
            // The orthogonal complement above (the last n-p column vectors of the matrix Q)
            // is then an orthogonal basis of Ker(A)
            var a_ns = new Matrix_.zeros(n, n-p);
            for (var j = p+1; j <= n; ++j) {
                for (var i = 1; i <= n; ++i) {
                    a_ns.data[(i-1)*a_ns.nbColumns + ((j-(p+1)+1)-1)] = q.data[(i-1)*q.nbColumns + (j-1)];
                }
            }           
        }
    }
   
    // Return the computed null space basis
    return a_ns;
}


/**
* @function linsolveBackwardSubstitution
*
* @summary Returns the solution of an upper triangular square system of linear equations.
*
* @description This function computes the solution of an upper triangular square system of linear equations Ax = b,
* as described in the algorithm 3.1.2 of the reference.
* 
* To be noted that the system of linear equations must be solvable (i.e., no zero elements must be present on the matrix A diagonal).
*
* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
*
* @param {Matrix_} A a n by n matrix.
* @param {Matrix_} b a n by 1 matrix.
* @return {<Matrix_} an n by 1 matrix x^* satisfying Ax^* = b.
*
* @example
* linsolveBackSubstitution(Matrix_([[1,2], [0,1]]), Matrix_([3,1]));
* // Matrix_([1,1])
*/
Matrix_.linsolveBackSubstitution = function(A, b, out) {
	// ------
	
	// Misc. checks
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	if (!A.isSquare()) {
		throw new Error('matrix is not square: ' + '(' + A.nbRows + ',' + A.nbColumns + ')');
	}
	if (A.nbRows !== b.nbRows) {
		throw new Error('matrix and second member sizes do not match: ' + '(' + A.nbRows + ',' + A.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns !== 1) {
		throw new Error('b is not a vector: ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}

	// ------
	
	// Initializations
	var n = A.nbColumns;

	// Result matrix allocation
	var x = allocateMatrix_(n, 1, out); // the solution vector

	// ------
	
	// Compute the solution to Ax = b, with a a square invertible upper triangular matrix,
	// using the row-oriented back substitution algorithm described in section 3.1.2 of the reference.
	for (var i = n; i >= 1; --i) {
		x.data[i-1] = b.data[i-1];
		for (var j = i + 1; j <= n; ++j) {
			x.data[i-1] = x.data[i-1] - A.data[(i-1) * A.nbColumns + (j-1)] * x.data[j-1];
		}
		var a_ii = A.data[(i-1) * A.nbColumns + (i-1)];
		if (a_ii == 0) {
			throw new Error('input matrix is not invertible: zero diagonal coefficient at index ' + i);
		}
		else {
			x.data[i-1] = x.data[i-1] / a_ii;
		}
	}

	// Return the computed solution
	return x;
}


/**
* @function linsolveExtendedKaczmarz
*
* @summary Returns the solution of a system of linear equations, using the extended Kaczmarz method.
*
* @description This function computes the solution of a system of linear equations Ax = b, using the 
* extended Kaczmarz iterative method in either its deterministic form (c.f. the first reference) 
* or its randomized form (c.f. the second reference).
* 
* To be noted that the extended Kaczmarz method allows to compute the solution of a system of 
* linear equations Ax = b even if this system is not solvable; in this case, the computed solution is
* the solution of minimum euclidian norm of the linear least squares problem argmin_x ||Ax - b||_2.
*
* @see <a href="https://www.tandfonline.com/doi/abs/10.1080/00207169508804364">Popa, Constantin. (1995). Least-squares solution of overdetermined inconsistent linear systems using Kaczmarz's relaxation. International Journal of Computer Mathematics. 55. 79-89.</a>
* @see <a href="https://epubs.siam.org/doi/abs/10.1137/120889897">Anastasios Zouzias and Nikolaos M. Freris, Randomized Extended Kaczmarz for Solving Least Squares, SIAM Journal on Matrix Analysis and Applications, 2013, Vol. 34, No. 2 : pp. 773-793</a>
*
* @param {Matrix_} A a m by n matrix.
* @param {Matrix_} b a m by 1 matrix.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.maxIter maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100000.
* @param {boolean} opt.randomized whether to use the deterministic extended Kaczmarz method (true) or the randomized extended Kaczmarz method (false); defaults to false.
* @return {<Matrix_} an n by 1 matrix x^* satisfying the following properties:
* - If the linear system of equations Ax = b is square with A invertible, then x^* is the unique solution of this system
* - If the linear system of equations Ax = b is not square (i.e., overdetermined or underdetermined) or is square with A non invertible, then x^* is the minimum euclidian norm solution
* of the linear least squares problem argmin_x ||Ax - b||_2
*
* @example
* linsolveExtendedKaczmarz(Matrix_([[3,2,-1], [2,-2,4], [-1, 0.5, -1]]), Matrix_([1,-2,0]));
* // Matrix_([1, -2, -2])
*/
Matrix_.linsolveExtendedKaczmarz = function(A, b, opt) {
	// ------
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-12;
	var maxIterations = opt.maxIter || 100000;
	var randomized = false;
	if (opt.randomized !== undefined) {
		randomized = opt.randomized;
	}
	
	
	// ------
	
	// Misc. checks
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	if (A.nbRows !== b.nbRows) {
		throw new Error('matrix and second member sizes do not match: ' + '(' + A.nbRows + ',' + A.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns !== 1) {
		throw new Error('b is not a vector: ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}


	// ------
	
	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	var x_k = Matrix_.zeros(n, 1); // the current solution
	var z_k = Matrix_.copy(b); // the current "adjusted" b
	
	var x_res = Matrix_.zeros(m, 1); // the x residuals vector
	var b_res = Matrix_.zeros(m, 1); // the b residuals vector
	var a_x_k = Matrix_.zeros(m, 1); // A*x_k

	// Preliminary computation of the Frobenius norm of A
	var a_frob_norm = A.matrixNorm('frobenius');
	var a_frob_norm_sq = a_frob_norm * a_frob_norm;	
	
	// Limit case: null matrix
	if (a_frob_norm === 0) {
		return x_k;
	}
	
	// Preliminary computation of the squares of the 2-norms of the rows of A.
	var a_rows_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(m) : new Array(m);
	for (var i = 1; i <= m; ++i) {
		var a_i_two_norm = A.vectorNorm('two', 'row', i);
		a_rows_two_norm_sq[i-1] = a_i_two_norm * a_i_two_norm;
	}
	
	// Preliminary computation of the squares of the 2-norms of the columns of A.
	var a_columns_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	for (var j = 1; j <= n; ++j) {
		var alpha_j_two_norm = A.vectorNorm('two', 'column', j);
		a_columns_two_norm_sq[j-1] = alpha_j_two_norm * alpha_j_two_norm;
	}
	
	// ------
	
	// Deterministic Extended Kaczmarz, c.f. algorithm R of the first reference.
	if (randomized == false) {
		// Main loop until convergence, guaranteed as per theorem 3.1 of the first reference.	
		var iter = 0;
		while (true) {
			// Update the number of iterations
			++iter;

			// Check the number of iterations
			if (maxIterations !== -1 && iter > maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			
			// Orthogonally project the current iterate z_k onto the hyperplane generated 
			// by the columns A(:,j), j=1..n
			for (var j = 1; j <= n; ++j) {
				if (a_columns_two_norm_sq[j-1] == 0) {
					continue;
				}
				
				// Compute <A(:,j)/z_k>
				var a_j_z_k = 0;
				for (var i = 1; i <= m; ++i) {
					a_j_z_k += A.data[(i-1) * A.nbColumns + (j-1)] * z_k.data[(i-1) * z_k.nbColumns + 0];
				}
				
				// Update z_k: z_k+1 = z_k - <A(:,j)/z_k>/||A(:,j)||_2^2 * A(:,j)
				for (var i = 1; i <= m; ++i) {
					z_k.data[(i-1) * z_k.nbColumns + 0] -= a_j_z_k / a_columns_two_norm_sq[j-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
				}
			}
			
			// Orthogonally project the current iterate x_k onto the solution hyperplane 
			// of <A(i,:)/x_k> = b(i) - z_k(i), i=1..m
			for (var i = 1; i <= m; ++i) {
				if (a_rows_two_norm_sq[i-1] == 0) {
					continue;
				}
				
				// Compute r_k = <A(i,:)/x_k> - (b(i) - z_k(i))
				var a_i_x_k = 0;
				for (var j = 1; j <= n; ++j) {
					a_i_x_k += A.data[(i-1) * A.nbColumns + (j-1)] * x_k.data[(j-1) * x_k.nbColumns];
				}
				var r_k = a_i_x_k - (b.data[(i-1) * b.nbColumns + 0] - z_k.data[(i-1) * z_k.nbColumns + 0]);
				
				// Update x_k: x_k+1 = x_k - r_k/||A(i,:)||_2^2 * A(i,:)
				for (var j = 1; j <= n; ++j) {
					x_k.data[(j-1) * x_k.nbColumns] -= r_k / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
				}
			}

			// Convergence condition (adapted from formula 4.3 of the first reference):
			// - ||Ax_k - (b - z_k)||_2 <= eps * ||A||_f * ||x_k||_2
			
			// Compute ||x_k||_2
			var x_k_two_norm = x_k.vectorNorm('two');
			
			// Compute ||Ax_k - (b - z_k)||_2
			a_x_k =  Matrix_.xy(A, x_k, a_x_k);
			b_res = Matrix_.xmy(b, z_k, b_res);
			x_res = Matrix_.xmy(a_x_k, b_res, x_res);
			var x_res_two_norm = x_res.vectorNorm('two');
			
			// Test the convergence condition		
			if (x_res_two_norm <= eps * a_frob_norm * x_k_two_norm) {
				break;
			}
		}
	}
	
	// Randomized Extended Kaczmarz, c.f. algorithm 3 of the second reference.
	else {
		// Initializations
		var ta_z_k = Matrix_.zeros(n, 1); // A^t*z_k
		
		// ------

		// Preliminary computation of the probabilities q_i, with their associated sampler.
		var q = typeof Float64Array === 'function' ? new Float64Array(m) : new Array(m);
		for (var i = 1; i <= m; ++i) {
			q[i-1] = a_rows_two_norm_sq[i-1]/a_frob_norm_sq;
		}
		var qSampler = new aliasMethodSampler_(q);
	
		// Preliminary computation of theprobabilities p_j with their associated sampler.
		var p = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var j = 1; j <= n; ++j) {
			p[j-1] = a_columns_two_norm_sq[j-1]/a_frob_norm_sq;
		}
		var pSampler = new aliasMethodSampler_(p);
	
		
		// ------
		
		// Main loop until convergence, guaranteed as per theorem 4.1 of the second reference.	
		var iter = 0;
		while (true) {
			// Update the number of iterations
			++iter;

			// Check the number of iterations
			if (maxIterations !== -1 && iter > maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			
			// Pick a column index j with probability p_j
			var j = pSampler.sample() + 1;
			
			// Orthogonally project the current iterate z_k onto the hyperplane generated 
			// by the column A(:,j)
				// Compute <A(:,j)/z_k>
			var a_j_z_k = 0;
			for (var i = 1; i <= m; ++i) {
				a_j_z_k += A.data[(i-1) * A.nbColumns + (j-1)] * z_k.data[(i-1) * z_k.nbColumns + 0];
			}
			
				// Update z_k: z_k+1 = z_k - <A(:,j)/z_k>/||A(:,j)||_2^2 * A(:,j)
			for (var i = 1; i <= m; ++i) {
				z_k.data[(i-1) * z_k.nbColumns + 0] -= a_j_z_k / a_columns_two_norm_sq[j-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}

			// Pick a row index i with probability q_i
			var i = qSampler.sample() + 1;
			
			// Orthogonally project the current iterate x_k onto the solution hyperplane 
			// of <A(i,:)/x_k> = b(i) - z_k(i)
				// Compute r_k = <A(i,:)/x_k> - (b(i) - z_k(i))
			var a_i_x_k = 0;
			for (var j = 1; j <= n; ++j) {
				a_i_x_k += A.data[(i-1) * A.nbColumns + (j-1)] * x_k.data[(j-1) * x_k.nbColumns];
			}
			var r_k = a_i_x_k - (b.data[(i-1) * b.nbColumns + 0] - z_k.data[(i-1) * z_k.nbColumns + 0]);

				// Update x_k: x_k+1 = x_k - r_k/||A(i,:)||_2^2 * A(i,:)
			for (var j = 1; j <= n; ++j) {
				x_k.data[(j-1) * x_k.nbColumns] -= r_k / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
			
			// Convergence conditions every 8 min(m, n) iterations: 
			// - ||Ax_k - (b - z_k)||_2 <= eps * ||A||_f * ||x_k||_2
			// - ||A^tz_k||_2 <= eps * ||A||_f^2 * ||x_k||_2
			if (iter % 8 * Math.min(m, n) === 0) {
				// Compute ||x_k||_2
				var x_k_two_norm = x_k.vectorNorm('two');
				
				// Compute ||Ax_k - (b - z_k)||_2
				a_x_k =  Matrix_.xy(A, x_k, a_x_k);
				b_res = Matrix_.xmy(b, z_k, b_res);
				x_res = Matrix_.xmy(a_x_k, b_res, x_res);
				var x_res_two_norm = x_res.vectorNorm('two');
				
				// Compute ||A^tz_k||_2
				ta_z_k = Matrix_.txy(A, z_k, ta_z_k);
				var ta_z_k_two_norm = ta_z_k.vectorNorm('two');
				
				if (x_res_two_norm <= eps * a_frob_norm * x_k_two_norm && 
					ta_z_k_two_norm <= eps * a_frob_norm_sq * x_k_two_norm) {
					break;
				}
			}
		}
	}

	// ------
	
	// Return the computed solution
	return x_k;
}

/**
 * @file Functions related to covariance matrix computations for portfolio allocation.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



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



/**
 * @file Misc. combinatorics functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.aliasMethodSampler_ = aliasMethodSampler_;
self.compositionsIterator_ = compositionsIterator_;
self.subsetsIterator_ = subsetsIterator_;
self.randomKSubsetIterator_ = randomKSubsetIterator_;
self.binomial_ = binomial_;
/* End Wrapper private methods - Unit tests usage only */
 

/**
* @function aliasMethodSampler_
*
* @summary Returns a function to generate random values sampled from a discrete
* finite probability distribution.
*
* @description This function constructs a function to generate random values from the set
* {0,...,n-1} sampled according to the provided discrete finite probability distribution
* {p_0,...,p_n-1}.
* 
* The algorithm used is the Vose's algorithm, which is a variation of the alias method
* allowing to sample random values from a finite discrete probability distribution in O(1) time
* after a O(n) time preprocessing step, c.f. the reference.
* 
* @see <a href="https://doi.org/10.1109/32.92917">M. D. Vose, A linear algorithm for generating random numbers 
* with a given distribution, IEEE Transactions on Software Engineering, vol. 17, no. 9, pp. 972-975, Sep 1991.</a>
*
* @param {Array.<number>} p, an array of n positive real numbers p_0,...,p_n-1 with sum_i p_i = 1.
* @return {function} a function to be used through its .sample() method, generating an integer i from the set
* {0,...,n-1} with probability p_i.
*
* @example
* var mySampler = new aliasMethodSampler_([0, 0.1, 0.4, 0.5]);
* mySampler.sample();
* // 3;
*/
function aliasMethodSampler_(p) {
	// ----
	// init function, c.f. paragraph B of section III of the reference.
	// ----

	// Initializations.
	this.prob = typeof Float64Array === 'function' ? new Float64Array(p.length) : new Array(p.length);
	this.alias = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);

	// TODO: Checks on probabilities (positive, sum to one)

	// Computation of the average probability.
    var avgProb = 1 / p.length;
		 
	// Initializations of the small and large stacks, together with their associated indexes.
	var small = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var s = 0;
	var large = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var l = 0;
		
	// Population of the small and large stacks with the probabilities indexes.
	for (var j = 0; j < p.length; ++j) {
		if (p[j] > avgProb) {
			large[l] = j;
			++l;
		}
		else {
			small[s] = j;
			++s;
		}
	}
		
	// Main loop of the algorithm, populating the prob and alias arrays.
	var p = p.slice(0); // local copy of the probabilities, as they are updated below
	while (s > 0 && l > 0) {
		// Get the index of the small and the large probabilities.
		--s;
		var j = small[s];
		
		--l;
		var k = large[l];
		
		// Update the prob and alias arrays.
		this.prob[j] = p.length * p[j];
		this.alias[j] = k;
		
		// Update the probabilities.
		p[k] = p[k] + (p[j] - avgProb);
		
		// Update the large and small stacks.
		if (p[k] > avgProb) {
			large[l] = k;
			++l;
		}
		else {
			small[s]= k;
			++s;
		}
	}
		
	// Process the remaining elements of the small stack.
	while (s > 0) {
		--s;
		this.prob[small[s]] = 1;
	}
	
	// Process the remaining elements of the large stack.
	//
	// Theoretically not needed, but due to round off errors, practicaly needed.
	while (l > 0) {
		--l;
		this.prob[large[l]] = 1;
	}
	
	
	// ----
	// rand function, c.f. paragraph A of section III of the reference.
	// ----
	
	/**
	* @function sample
	*
	* @summary Returns a random value sampled from the underlying probability distribution.
	*
	* @description This function computes a random value sampled from the underlying 
	* probability distribution using the method described in the reference.
	*
	* @memberof aliasMethodSampler_
	* @return {number} an integer i belonging to the set {0,...,n-1} with probability p_i.
	*/
    this.sample = function() {
		var u = Math.random() * this.prob.length; // Uniform real number belonging to [0..n[
		var j = Math.floor(u);
		if (u - j <= this.prob[j]) {
			return j;
		}
		else {
			return this.alias[j];
		}
    }
}

	
 /**
* @function compositionsIterator_
*
* @summary Returns an iterator to compute all the compositions of a non negative integer.
*
* @description This function constructs an iterator to compute all the k-compositions of a non-negative integer n, 
* using the algorithm NEXCOM described in section 5 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n a non-negative integer whose composition are desired.
* @param {number} k a non-negative integer indicating the number of parts of desired composition of n.
* @return {function} a function to be used as an iterator through its .next() method, computing all 
* the k-compositions of n.
*
* @example
* var myIterator = new compositionsIterator_(6, 3);
* myIterator.next(); myIterator.next();
* // [true, [6,0,0]]; [true, [5,1,0]];
*/
function compositionsIterator_(n, k) {
	// Initialize n and k
	this.n = n;
	this.k = k;
	
	// Variables required for NEXTCOM internal computations,
	// initialized so as to generate the first composition upon
	// the first call to .next() function.
	this.mtc = false;
	this.r = new Array(k);
	this.t = this.n;
	this.h = 0;

	/**
	* @function next
	*
	* @summary Returns the next composition of a non negative integer.
	*
	* @description This function computes the next k-composition of a non negative integer n.
	*
	* The initial k-composition computed by the first call to this function is n00...0, and each subsequent call to 
	* this function will result in a new k-composition until the final k-composition 00...0n is reached.
	*
	* A subsequent call to this function when the final k-composition has been reached will result in
	* the recomputation of all the k-compositions of n, starting from the initial k-composition.
	*
	* @memberof compositionsIterator_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one k-composition of n remains to be computed
	* and arr[1] an array of k elements containing the computed k-composition of n.
	*
	*/
	this.next = function() {
		if (this.mtc) { // There is still a composition to generate
			if (this.t > 1) {
				this.h = 0;
			}
			this.h++;
			this.t = this.r[this.h - 1];
			this.r[this.h - 1] = 0;
			this.r[0] = this.t - 1;
			++this.r[this.h];
		}
		else  { 
		    // No more composition to generate, so, (re) generation of the first composition, equals to n00...0
			this.r[0] = this.n;
			for (var i = 1; i <= this.k - 1; ++i) {
				this.r[i] = 0;
			}
		}
		
		// End logic
		this.mtc = (this.r[this.k - 1] != this.n);
		
		// Return a copy of the r array, so that callers can alter it
		return [this.mtc, this.r.slice()];
	}
}


/**
* @function randomKSubsetIterator_
*
* @summary Returns an infinite iterator to compute random k-subsets of a n-set.
*
* @description This function constructs an iterator to compute random k-subsets of the n-set {1,...,n}, 
* using both the algorithms RANKSB and RKS2 described in section 4 of the reference.
*
* From the discussion following the examples in the reference, the random k-subsets are probably generated
* uniformly, but this is not written in the reference.
*
* The algorithm used to compute the random k-subsets is either RANKSB when k < n/2,
* or RKS2 when k >= n/2, so that performances are in O(k).
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
*
* @param {number} n the number of elements of the n-set {1,...,n} whose k-subsets are desired, a non-negative integer.
* @param {number} k a non-negative integer, with 0 <= k <= n.
* @return {function} a function to be used as an iterator through its .next() method, computing random 
* k-subsets of the n-set {1,...,n}.
*
* @example
* var myIterator = new randomKSubsetIterator_(6, 3);
* myIterator.next();
* // [1, 2, 5];
*/
function randomKSubsetIterator_(n, k) {
	// Initialize n and k
	this.n = n;
	this.k = k;
	
	// Initialize the array to hold the k-subsets
	this.a = new Array(k);

	/**
	* @function next_ranksb
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the algorithm RANKSB of the reference.
	*
	* @memberof randomKSubsetIterator_
	* @return {Array.<number>} a random k-subset of the n-set {1,...,n}, a sorted array of k increasing strictly positive integers.
	*
	*/
	this.next_ranksb = function() {
		// Step A - Initialization of a
		for (var i = 1; i <= this.k; ++i) {
			this.a[i-1] = Math.floor((i - 1) * this.n / this.k);
		}
		
		// Step B
		// Note: in the reference, the c variable is initialized to k and is decremented until 0
		// each time a generated x is accepted: this is a reverse for loop in disguise.
		var x;
		var l;
		for (var c = this.k; c > 0; --c) {
			do {
				var u = Math.random();
				x = 1 + Math.floor(u * this.n);
				l = 1 + Math.floor((x * this.k - 1) / this.n);
			} while (x <= this.a[l-1]);
			this.a[l-1] = this.a[l-1] + 1;
		}
		var p = 0;
		var s = this.k;
		
		// Step C
		// Note: in the reference, the i variable is initialized to 0 and is incremented
		// until k each time: this is a for loop in disguise.
		for (var i = 1; i <= this.k; ++i) {
			if (this.a[i-1] == Math.floor((i - 1) * this.n / this.k)) {
				this.a[i-1] = 0;
			}
			else {
				p = p + 1;
				var m = this.a[i-1];
				this.a[i-1] = 0;
				this.a[p-1] = m;
			}
		}
		
		// Step D
		// Note: in the reference, the p variable is initialized to whatever value it has, and is decremented
		// until 0 each time: this is a reverse for loop in disguise.
		for (; p > 0; --p) {
			l = 1 + Math.floor((this.a[p-1] * this.k - 1) / this.n);
			var delta_s = this.a[p-1] - Math.floor((l - 1) * this.n / this.k);
			this.a[p-1] = 0;
			this.a[s-1] = l;
			s = s - delta_s;			
		}
		l = k;
		
		// Steps E to H
		// Note: in the reference, the l variable is initialized at this step to k, and is decremented
		// until 0 each time: this is a reverse for loop in disguise.
		var r;
		for (; l > 0; --l) {
			// Step E
			var m_0;
			if (this.a[l-1] != 0) {
				r = l;
				m_0 = 1 + Math.floor((this.a[l-1] - 1) * this.n / this.k);
				m = Math.floor(this.a[l-1] * this.n / this.k) - m_0 + 1;
			}

			// Step F
			var u = Math.random();
			x = m_0 + Math.floor(u * m);
			i = l;
			
			// Step G
			++i;
			while (i <= r && x >= this.a[i-1]) {
				this.a[i-2] = this.a[i-1];
				x = x + 1;
				++i;
			}
			
			// Step H
			this.a[i-2] = x;
			m = m - 1;
		}
		
		// Return a copy of the computed array
		return this.a.slice();
	}
	
	/**
	* @function next_rks2
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the algorithm RKS2 of the reference.
	*
	* @memberof randomKSubsetIterator_
	* @return {Array.<number>} a random k-subset of the n-set {1,...,n}, a sorted array of k increasing strictly positive integers.
	*
	*/
	this.next_rks2 = function() {
		// Initializations
		var c_1 = this.k;
		var c_2 = this.n;
		var k_0 = 0;
		var i = 0;

		// Main loop of the RKS2 algorithm
		while (c_1 > 0) {
			++i;
			var u = Math.random();
			if (u <= c_1/c_2) {
				c_1 = c_1 - 1;
				this.a[k_0] = i;
				k_0 = k_0 + 1; // this line is inversed compared to the reference because of JavaScript arrays starting at index 0
			}
			c_2 = c_2 - 1;
		}
		
		// Return a copy of the computed array
		return this.a.slice();
	}
	
	// Initialize the appropriate iterator to keep the required labor to O(k) uniformly for 1 <= k <= n
	if (k < n/2) {
		this.next = this.next_ranksb;
	}
	else {
		this.next = this.next_rks2;
	}
}


/**
* @function subsetsIterator_
*
* @summary Returns an iterator to compute all the subsets of a n-set.
*
* @description This function constructs an iterator to compute all the subsets of the n-set {1,...,n}, 
* using the algorithm NEXSUB described in section 1 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="*https://en.wikipedia.org/wiki/Power_set">Power set</a>
*
* @param {number} n the number of elements of the n-set {1,...,n} whose subsets are desired, a non-negative integer.
* @return {function} a function to be used as an iterator through its .next() method, computing all 
* the subsets of the n-set {1,...,n}.
*
* @example
* var myIterator = new subsetsIterator_(5);
* myIterator.next(); myIterator.next();
* // [true, []]; [true, [1]];
*/
function subsetsIterator_(n) {
	// Initialize n
	this.n = n;
	
	// Variables required for NEXSUB internal computations,
	// initialized so as to generate the first subset upon
	// the first call to .next() function.
	this.mtc = false;
	this.iin = new Array(n);
	this.ncard = 0;
	
	/**
	* @function next
	*
	* @summary Returns the next subset of a set.
	*
	* @description This function computes the next subset of the n-set {1,...,n}.
	*
	* The initial subset computed by the first call to this function is the empty subset {}, and each subsequent call to 
	* this function will result in a new subset until the final subset {1,...,n} is reached.
	*
	* A subsequent call to this function when the final subset {1,...,n} has been reached will result in
	* the recomputation of all the subsets, re-starting from the initial subset.
	*
	* @memberof subsetsIterator_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one subset
	* of the n-set {1,...,n} remains to be computed and arr[1] an array containing the computed sorted subset
	* of the n-set {1,...,n}.
	*
	*/
	this.next = function() {
		// The output array containing the computed subset
		var nextSubset = new Array(0);
		
		if (this.mtc) { // There is still a subset to generate
			var j = 0;
			if (this.ncard % 2 != 0) {
				++j;
				while (this.iin[j - 1] == 0) {
					++j;
				}
			}
			this.iin[j] = 1 - this.iin[j];
			this.ncard = this.ncard + 2*this.iin[j] - 1;

			// Build the output array
			nextSubset = new Array(this.ncard);
			var idx = 0;
			for (var i = 0; i <= this.n - 1; ++i) {
				if (this.iin[i] == 1) {
					nextSubset[idx++] = i + 1;
				}
			}
			
			// End logic
			this.mtc = (this.ncard != this.iin[this.n -1]);
		}
		else  { 
		    // No more subset to generate, so, (re) generation of the first subset, equals to {}
			for (var i = 0; i <= this.n - 1; ++i) {
				this.iin[i] = 0;
			}
			
			// The output array is already built in this case (empty)
			
			// Specific end logic
			this.mtc = true;
		}

		// Return the computed array, not used anymore by this function
		return [this.mtc, nextSubset];
	}
}


/**
* @function binomial_
*
* @summary Returns a binomial coefficient.
*
* @description This function computes the k-th binomial coefficient of order n, which is
* the coefficient of the x^k term in the polynomial expansion of the binomial power (1 + x)^n.
*
* This coefficient is also the number of ways to choose a subset of k elements,
* disregarding their order, from a set of n elements.
*
* The algorithm used is a multiplicative formula, c.f. the reference.
*
* @see <a href="*https://en.wikipedia.org/wiki/Binomial_coefficient">Binomial coefficient</a>
*
* @param {number} n a non-negative integer.
* @param {number} k a non-negative integer, with 0 <= k <= n.
* @return {number} the computed binomial coefficient.
*
* @example
* binomial_(7, 5);
* // 21
*/
function binomial_(n, k) {
    // Checks
    if (n < 0) {
        throw new Error('n must be a positive integer');
    }
    if (k < 0) {
        throw new Error('k must be a positive integer');
    }
    if (k > n) {
        throw new Error('k must be less than or equal to n');
    }
	  
	// Compute the binomial coefficient using the multiplicative formula of the reference.
	var val = 1;
    for (var i = 1; i <= k; ++i) {
        val *= (n + 1 - i);
		val /= i; // Note: separating the computation of val in two steps guarantees (unless compiler optimisations) that val is an integer.
    }
	
	// Return it
    return val;
}
     
/**
 * @file Misc. statistical functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.hypot_ = function(x, y) { return hypot_(x, y); }
self.rank_ = function(x, order) { return rank_(x, order); }
self.ftca_ = function(correlationMatrix, threshold) { return ftca_(correlationMatrix, threshold); }
/* End Wrapper private methods - Unit tests usage only */
 

/**
* @function hypot_
*
* @summary Returns the square root of the sum of the squares of two numbers (i.e., the hypotenuse).
*
* @description This function computes the value of sqrt(abs(x)^2 + abs(y)^2) in a way
* to avoid as much as possible underflow and overflow.
*
* @see <a href="https://en.wikipedia.org/wiki/Hypot#Implementation">Hypot</a>
*
* @param {number} x a real number.
* @param {number} y a real number.
* @return {number} the value of sqrt(abs(x)^2 + abs(y)^2), a real number.
*
* @example
* hypot_(3, 4);
* // 5
*/
function hypot_(x, y) {
    // Initialization
	var r = 0;
    
	// Main algorithm
	var absX = Math.abs(x);
	var absY = Math.abs(y);
	if (absX > absY) {
	   r = y/x;
	   r = absX * Math.sqrt(1 + r*r);
    } 
	else if (y != 0) {
	   r = x/y;
	   r = absY * Math.sqrt(1 + r*r);
    }
	else {
	   r = 0;
    }
    
	// Return the computed value
	return r;
}

 
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
* @param {Array.<Array.<number>>} correlationMatrix the correlation matrix (rho_ij),i,j=1..n, an array arr of n arrays of n real numbers 
* satisfying arr[i-1][j-1] = rho_ij, i,j=1..n, where n is a strictly positive natural integer.
* @param {number} threshold the correlation threshold to use in the FTCA algorithm, a real number typically belonging to interval [-1, 1].
* @return {Array.<Array.<number>>} the list of clusters as computed by the FTCA algorithm, array of m arrays of strictly positive integers representing the indexes of the elements in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
*
* @example
* ftca_([[1, 0], [0,1]]), 0.5);
*  // [[2],[1]]
*/
function ftca_(correlationMatrix, threshold) {
	// Decode the optional threshold
	var threshold = threshold;
	if (threshold === undefined) {
		threshold = 0.5;
	}
	
	// Convert the correlation matrix to matrix format
	var correlationMatrix = new Matrix_(correlationMatrix);
	
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
			var subCorrMat = correlationMatrix.submatrix(unassignedElementsIdx, unassignedElementsIdx);
			
			// Compute the average correlation of each unassigned element to all the other unassigned elements
			// Computation is done for each row
			var subCorrMatRows = subCorrMat.toRowArray(function(i, j, val) {
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
					// Note: At this stage, the LC element cannot have been assigned to the Hc cluster above if LC <> HC, since
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
 * @file Misc. optimisation functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.lpsolvePrimalDualHybridGradient_ = lpsolvePrimalDualHybridGradient_;
/* End Wrapper private methods - Unit tests usage only */
 
 
 
/**
* @function lpsolvePrimalDualHybridGradient_
*
* @summary Returns an optimal solution to a linear program, using a primal-dual hybrid gradient algorithm.
*
* @description This function computes an optimal solution to a linear program using a 
* preconditioned primal-dual hybrid gradient algorithm, c.f. the first reference.
*
* The linear program to solve is assumed to be provided in the following format:
*
* min <c/x>
*
* s.t. Ae*x = be (equality constraints)
*      Ai*x <= bi (inequality constraints)
*      lb <= x <= ub (bound constraints)
*
* with:
* - c an n by 1 matrix
* - Ae an optional me by n matrix
* - be an optional me by 1 matrix
* - Ai an optional mi by n matrix
* - bi an optional mi by 1 matrix
* - lb an optional n by 1 matrix, which can contain negative infinity values (-Infinity) corresponding to unbounded variables on the negaitve axis
* - ub an optional n by 1 matrix, which can contain positive infinity values (Infinity) corresponding to unbounded variables on the positive axis
*
* and with:
* - lb assumed to be an n by 1 matrix made of zeroes if not provided
* - ub assumed to be an n by 1 matrix made of positive infinity values if not provided
* 
* To be noted that the algorithm used internally requires the linear problem to be feasible and bounded
* (i.e., admit a finite optimal solution) to converge.
* 
* @see <a href="http://ieeexplore.ieee.org/document/6126441/">T. Pock and A. Chambolle, "Diagonal preconditioning for first order primal-dual algorithms in convex optimization" 2011 International Conference on Computer Vision, Barcelona, 2011, pp. 1762-1769.</a>
* @see <a href="https://arxiv.org/abs/1305.0546">Tom Goldstein, Min Li, Xiaoming Yuan, Ernie Esser, Richard Baraniuk, "Adaptive Primal-Dual Hybrid Gradient Methods forSaddle-Point Problems", 05/2013, eprint arXiv:1305.0546</a>
*
* @param {Matrix_} Ae an optional me by n matrix; must be null if not provided.
* @param {Matrix_} be an optional me by 1 matrix; must be null if not provided.
* @param {Matrix_} Ai an optional mi by n matrix; must be null if not provided.
* @param {Matrix_} bi an optional mi by 1 matrix; must be null if not provided.
* @param {Matrix_} c an n by n matrix.
* @param {Matrix_} lb an optional n by 1 matrix; must be null if not provided.
* @param {Matrix_} ub an optional n by 1 matrix; must be null if not provided.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-08.
* @param {number} opt.maxIter maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100000.
* @return {Array<Object>} an array arr containing two elements: 
* - arr[0] an n by 1 matrix containing the optimal solution x^* to the linear program
* - arr[1] the optimal value of the function x -> <c/x>, i.e. <c/x^*>
*
* @example
* lpsolvePrimalDualHybridGradient_(X);
* // X
*/
 function lpsolvePrimalDualHybridGradient_(Ae, be, Ai, bi, c, lb, ub, opt) {
    // ------
    
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-08;
	var maxIterations = opt.maxIter || 100000;
	

	// ------

	// Misc. checks
	var eqContraints = false;
	var ineqContraints = false;
	var boundContraints = false;
	
	if (Ae !== null && be !== null) {
		eqContraints = true;
	}
	else if (Ae !== null && be === null) {
		throw new Error('equality constraints vector is missing');
	}
	else if (Ae === null && be !== null) {
		throw new Error('equality constraints matrix is missing');
	}

	if (Ai !== null && bi !== null) {
		ineqContraints = true;
	}
	else if (Ai !== null && bi === null) {
		throw new Error('inequality constraints vector is missing');
	}
	else if (Ai === null && bi !== null) {
		throw new Error('inequality constraints matrix is missing');
	}
	
	if (lb !== null && ub !== null) {
	    boundContraints = true;
	}
	else if (lb !== null && ub === null ) {
    	throw new Error('upper bounds constraints vector is missing');
	}
	else if (lb === null && ub !== null ) {
    	throw new Error('lower bounds constraints vector is missing');
	}
	
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
	
	if (boundContraints) {
		if (lb.nbRows !== null) {
			if (!(lb instanceof Matrix_)) {
				throw new Error('sixth input must be a matrix');
			}
			if (lb.nbRows !== c.nbRows) {
				throw new Error('sixth input number of rows and fifth input number of rows do not match: ' + lb.nbRows + '-' + c.nbRows);
			}
		}
        if (ub.nbRows !== null) {
			if (!(ub instanceof Matrix_)) {
				throw new Error('seventh input must be a matrix');
			}
			if (ub.nbRows !== c.nbRows) {
				throw new Error('seventh input number of rows and fifth input number of rows do not match: ' + ub.nbRows + '-' + c.nbRows);
			}
		}
	}

	
	// ------
	
	// Initializations
	// Constraints
	var me = 0;
	var ye_k = null;
	var ye_kp = null;
	var res_ye_kp_ye_k = null;
	if (eqContraints) {
	    me = Ae.nbRows; // the number of equality constaints
    	ye_k = Matrix_.zeros(me, 1); // the dual iterate ye_k
    	ye_kp = Matrix_.zeros(me, 1); // the dual iterate ye_k+1
    	res_ye_kp_ye_k = Matrix_.zeros(me, 1); // the residual ye_kp - ye_k
	}
	var mi = 0;
	var yi_k = null;
	var yi_kp = null;
	var res_yi_kp_yi_k = null;
	if (ineqContraints) {
	    mi = Ai.nbRows; // the number of inequality (<=) constaints
    	yi_k = Matrix_.zeros(mi, 1); // the dual iterate yi_k
    	yi_kp = Matrix_.zeros(mi, 1); // the dual iterate yi_k+1
    	res_yi_kp_yi_k = Matrix_.zeros(mi, 1); // the residual yi_kp - yi_k
	}
	var m = me + mi; // the total number of constraints
	
    // Variables
    var n = c.nbRows;
	var x_k = Matrix_.ones(n, 1); // the primal iterate x_k
	var x_kp = Matrix_.ones(n, 1); // the primal iterate x_k+1
	var z_k = Matrix_.ones(n, 1); // the relaxed iterate z_k = 2*x_k+1 - x_k
	var res_x_kp_x_k = Matrix_.zeros(n, 1); // the residual x_kp - x_k

    // Misc.
	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var ttmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var tmp_vec_me = null;
	if (eqContraints) {
		tmp_vec_me = Matrix_.zeros(me, 1); // a temporary placeholder vector of dimension me
	}
	var tmp_vec_mi = null;
	if (ineqContraints) {
		tmp_vec_mi = Matrix_.zeros(mi, 1); // a temporary placeholder vector of dimension mi
	}
	
	// Computation of the diagonal matrices T and S = [Se Si]^t with alpha = 1
	// and mu = 0.9995, nu = 0.9995 (so that mu * nu < 1), c.f. formula 10 and
	// remark 3 of the first reference.
	var mu = 0.9995;
	var nu = 0.9995;
	var T = Matrix_.fill(n, 1, 
				function(i,j) { 
						var columnNorm = 0;
						if (eqContraints) {
						    var aeColNorm = Ae.vectorNorm('one', 'column', i); 
							columnNorm += aeColNorm;
						}
						if (ineqContraints) {
						    var aiColNorm = Ai.vectorNorm('one', 'column', i);
							columnNorm += aiColNorm;
						}
						if (columnNorm == 0) { // in case there is no coefficient in the column, a value of 1 is harmless
							columnNorm = 1;
						}
						return mu * 1/columnNorm;
				});
	
	var Se = null;
	if (eqContraints) {
	    Se = Matrix_.fill(me, 1, 
	                      function(i,j) { 
						    var aeRowNorm = Ae.vectorNorm('one', 'row', i);
							if (aeRowNorm == 0) { // in case there is no coefficient in the row, a value of 1 is harmless
								aeRowNorm = 1;
							}
							return nu * 1/aeRowNorm;
				          });
	}
	var Si = null;
	if (ineqContraints) {
	    Si = Matrix_.fill(mi, 1, 
	                      function(i,j) { 
						    var aiRowNorm = Ai.vectorNorm('one', 'row', i);
							if (aiRowNorm == 0) { // in case there is no coefficient in the row, a value of 1 is harmless
								aiRowNorm = 1;
							}
			        	    return nu * 1/aiRowNorm;
				          });
	}
	
	// ------
	
	// Main loop of the algorithm, c.f. formula 18 of the first reference.
	//
	// The convergence is guaranteed by the theorem 1 of the first reference
	// in case the primal linear program admit a finite optimal solution.
	//
	//
	// To be noted that the formulation used below is slightly different from
	// the formulation in formula 17 of the first reference (equality constraints only).
	//
	// To arrive at this form, the formulation below mixes formula 17 of the first reference,
	// formula 26 of the second reference, and additional bound constraints on x.
	var iter = 0;
	while (true) {
		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}

		// Update the number of iterations
		++iter;

		// Primal update: 
		// - x_k+1 = proj(x_k - T*(A^t*y_k + c))_[0, +infinity[, with A = [Ae Ai]^t and y_k = [ye yi]^t (no bound constraints)
		// - x_k+1 = proj(x_k - T*(A^t*y_k + c))_[lb, ub], with A = [Ae Ai]^t and y_k = [ye yi]^t (bound constraints)
		if (eqContraints && ineqContraints) {
			x_kp = Matrix_.xpy(Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), Matrix_.txy(Ai, yi_k, ttmp_vec_n), tmp_vec_n), c, x_kp);
		}
		else if (eqContraints) {
			x_kp = Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), c, x_kp);
		}
		else if (ineqContraints) {
			x_kp = Matrix_.xpy( Matrix_.txy(Ai, yi_k, tmp_vec_n), c, x_kp);
		}
		x_kp = Matrix_.xmy(x_k, Matrix_.elementwiseProduct(x_kp, T, tmp_vec_n), x_kp);
		if (boundContraints) {
			// Projection on the interval [lb, ub]
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < lb.data[i-1]) {
					x_kp.data[i-1] = lb.data[i-1];
				}
				else if (x_kp.data[i-1] > ub.data[i-1]) {
					x_kp.data[i-1] = ub.data[i-1];
				}
			}
		}
		else {			
			// Projection on the non-negative orthant
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < 0) {
					x_kp.data[i-1] = 0;
				}
			}
		}

		// Relaxed iterate update:
		// - z_k = 2*x_k+1 - x_k
		z_k = Matrix_.axpby(2, x_kp, -1, x_k, z_k);
		
		// Dual update:
		// - ye_k+1 = ye_k + Se*(Ae*z_k - be) (equality constraints)
		// - yi_k+1 = proj(yi_k + Si*(Ai*z_k - bi))_[0, +infinity[ (inequality constraints)
		if (eqContraints) {
			ye_kp = Matrix_.xpy(ye_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ae, z_k, tmp_vec_me), be, tmp_vec_me), Se, tmp_vec_me), ye_kp);
		}
		if (ineqContraints) {
			yi_kp = Matrix_.xpy(yi_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ai, z_k, tmp_vec_mi), bi, tmp_vec_mi), Si, tmp_vec_mi), yi_kp);
			
			// Projection on the non-negative orthant
			for (var i = 1; i <= mi; ++i) {
				if (yi_kp.data[i-1] < 0) {
					yi_kp.data[i-1] = 0;
				}
			}
		}
	
		// Convergence conditions for (x_k, y_k = [ye yi]^t) to be a saddle point of the min-max problem:
		// - Convergence of the primal iterates (relative) x_k: ||x_k+1 - x_k||_inf <= eps * ||x_k+1||_inf
		// - Convergence of the dual iterates (relative) y_k: ||y_k+1 - y_k||_inf <= eps * ||y_k+1||_inf
		res_x_kp_x_k = Matrix_.xmy(x_kp, x_k, res_x_kp_x_k);
		var res_x_kp_x_k_inf_norm = res_x_kp_x_k.vectorNorm('infinity');
		var x_kp_inf_norm = x_kp.vectorNorm('infinity');
		
		var res_ye_kp_ye_k_inf_norm = 0;
		var ye_kp_inf_norm = 0;
		if (eqContraints) {
		    res_ye_kp_ye_k = Matrix_.xmy(ye_kp, ye_k, res_ye_kp_ye_k);
		    res_ye_kp_ye_k_inf_norm = res_ye_kp_ye_k.vectorNorm('infinity');
											 
			ye_kp_inf_norm = ye_kp.vectorNorm('infinity');
		}
		
		var res_yi_kp_yi_k_inf_norm = 0;
		var yi_kp_inf_norm = 0;
		if (ineqContraints) {
		    res_yi_kp_yi_k = Matrix_.xmy(yi_kp, yi_k, res_yi_kp_yi_k);
		    res_yi_kp_yi_k_inf_norm = res_yi_kp_yi_k.vectorNorm('infinity');
			
			yi_kp_inf_norm = yi_kp.vectorNorm('infinity');
		}
		
		if (res_x_kp_x_k_inf_norm <= eps * x_kp_inf_norm  && 
		    res_ye_kp_ye_k_inf_norm <= eps * ye_kp_inf_norm &&
			res_yi_kp_yi_k_inf_norm <= eps * yi_kp_inf_norm) {
			break;
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
	}
	
	// Compute the objective function value.
	var fctVal = Matrix_.vectorDotProduct(c, x_kp);

	// Return the computed primal iterate and the associated objective
	// function value.
	return [x_kp, fctVal];		
 }
 
 

/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.simplexRationalRounding_ = simplexRationalRounding_;
self.simplexRandomSampler_ = simplexRandomSampler_;
self.simplexDeterministicRationalSampler_ = simplexDeterministicRationalSampler_;
self.simplexRationalGirdSearch_ = simplexRationalGirdSearch_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function simplexDeterministicRationalSampler_
*
* @summary Returns a function to generate all the points on a rational grid of the unit simplex of R^n.
*
* @description This function constructs a function to generate all the points on the k-th rational grid
* of the unit simplex of R^n, 1/k * I_n(k), c.f. the first reference.
* 
* The algorithm used is based on the enumeration of all the k-compositions of the integer n, c.f. the second reference.
*
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic
*  optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
* @see <a href="https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to generate points, 
* a natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing all 
* the points on the k-th rational grid of the unit simplex of R^n.
*
* @example
* var mySampler = new simplexDeterministicRationalSampler_(3, 10);
* mySampler.sample(); mySampler.sample(); ...; mySampler.sample();
* // [1, 0, 0]; [0.9, 0.1, 0]; ...; null
*/
function simplexDeterministicRationalSampler_(n, k) {
	// Initializations
	this.n = n;
	this.k = k;
	this.compositionIterator = new compositionsIterator_(k, n);
	this.compositionIteratorStatus = true;
	
	/**
	* @function sample
	*
	* @summary Returns a point on the k-th rational grid of the unit simplex of R^n.
	*
	* @description This function generates a point on the k-th rational grid of the unit simplex of R^n.
	*
	* Each call to this function results in the generation of a new point on the k-th rational grid
	* of the unit simplex of R^n, until exhaustion of all such points.
	*
	* @memberof simplexDeterministicRationalSampler_
	* @return {Array.<number>|null} an array of n real numbers corresponding to the coordinates of the generated point in R^n,
	* or null in case all such points have been generated.
	*
	*/
	this.sample = function() {
		// Return null in case there is no more samples to draw
		if (!this.compositionIteratorStatus) {
			return null;
		}
		
		// Generate a new k-composition of n
		var nextComposition = this.compositionIterator.next();
		
		// Compute the current rational grid point by normalizing the generated k-composition
		var x = nextComposition[1];
		for (var i = 0; i < x.length; ++i) {
			x[i] = x[i] / k;
		}

		// Update the internal iterator status
		this.compositionIteratorStatus = nextComposition[0];	
		
		// Return the point being sampled
		return x;
	}
}
 
 
/**
* @function simplexRandomSampler_
*
* @summary Returns a function to compute random points on the unit simplex of R^n.
*
* @description This function constructs a function to compute random points uniformly distributed on
* the unit simplex of R^n, using the algorithm 2 of the reference.
* 
* @see <a href="https://doi.org/10.1016/0377-2217(82)90161-8">R.Y. Rubinstein, Generating random vectors uniformly distributed inside and on 
* the surface of different regions, In European Journal of Operational Research, Volume 10, Issue 2, 1982, Pages 205-209</a>
*
* @param {number} n the dimension of the unit simplex of R^n, natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing random  
* points on the unit simplex of R^n.
*
* @example
* var mySampler = new simplexRandomSampler_(3);
* mySampler.sample();
* // [0.25, 0, 0.75]
*/
function simplexRandomSampler_(n) {
	// Initializations
	this.n = n;
	this.x = new Array(n); // the coordinates of a point being sampled
	
	/**
	* @function sample
	*
	* @summary Returns a random point on the unit simplex of R^n.
	*
	* @description This function computes a point choosen uniformly at random on the unit simplex of R^n,
	* using the algorithm 2 of the reference, which is linear with n (i.e., which has a time complexity of O(n)).
	*
	* @memberof simplexRandomSampler_
	* @return {Array.<number>} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample = function() {
		// Computation of n random variables from EXP(1), which will form the basis
		// of the coordinates of the point being sampled	
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from EXP(1) using the inverse method, with no need for the minus sign
			// as the negative sign would cancel out at the subsequent normalization step.
			var u = 1 - Math.random(); // belongs to ]0,1], so that the logarithm below is properly defined
			var e = Math.log(u);

			// Set the i-th coordinate of the point being sampled.
			this.x[i] = e;
			
			// Compute the running sum of the exponential variables, for the subsequant normalization step.
			sum += e;
		}

		// Normalization of the computed coordinates of the point being sampled, so that
		// they all belong to [0,1] and sum to 1.
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/sum;
		}
		
		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice();
	}
}
 
 

/**
* @function simplexRationalRounding_
*
* @summary Compute a rational point on the unit simplex that is the closest to a point on the unit simplex.
*
* @description Given a point x = (x_1,...,x_n) on the standard simplex of R^n, this function computes a proximal point xr = (xr_1,...,xr_n) on 
* the r-th rational grid of the unit simplex of R^n, 1/r * I_n(r), with I_n(r) the set of N^n containing the points m = (m_1,...,m_n) 
* statisfying sum m_i = r with r a strictly positive natural integer, so that the computed proximal point xr is one of the closest points to x 
* on this grid with respect to any norm in a large class, c.f. the first reference.
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-013-0126-2">M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258</a>
* @see <a href="https://arxiv.org/abs/1501.00014">Rama Cont, Massoud Heidari, Optimal rounding under integer constraints</a>
* 
* @param {Array.<number>} x a point belonging to the standard simplex of R^n, array of n real numbers.
* @param {number} r the indice of the rational grid of the unit simplex of R^n, 1/r * I_n(r), on which to compute the closest point to x, natural integer greater than or equal to 1.
* @return {Array.<number>} the computed closest point to x on the r-th rational grid of the unit simplex of R^n, array of n real numbers.
*
* @example
* simplexRationalRounding_([0.5759, 0.0671, 0.3570], 20);
* // [0.6, 0.05, 0.35]
*/
function simplexRationalRounding_(x, r) {
	// TODO: Checks, if enabled

	// Compute the integer and fractional parts of the coordinates of the input point multiplied by r, 
	// as described in paragraph 2 of the first reference.
	// In parallel, compute k as also defined in paragraph 2 of the first reference. 
	var k = r;
	var xPartsWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		//
		var rx = r * x[i];
		var integerPart = Math.floor(rx);
		var fractionalPart = rx - integerPart;
		
		//
		k -= integerPart;
		
		//
		xPartsWithIndexes[i] = [integerPart, fractionalPart, i];
	}

	// Re-order the coordinates according to decreasing values of the fractional parts, c.f. theorem 1 of the first reference.
	// In case the fractional parts are equal, re-order the coordinates according to decreasing values of the integer parts,
	// c.f. paragraph 3 of the second reference.
	xPartsWithIndexes.sort(function(a, b) {
		if (b[1] < a[1]) {
			return -1;
		}
		else if (b[1] > a[1]) {
			return 1;
		}
		else { // coordinates have equal fractional parts
			return b[0] - a[0];
		}
	}); 

	// Rounding rule: round the k largest fractional parts up to one and all other fractional parts down to zero,
	// as described in paragraph 2 of the first reference.
	var xr = new Array(x.length);
	for (var i = 0; i < k; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = (xPartsWithIndexes[i][0] + 1) / r;
	}
	for (var i = k; i < xPartsWithIndexes.length; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = xPartsWithIndexes[i][0] / r;
	}

	// Return the computed point
	return xr;
}

/**
* @function simplexRationalGirdSearch_
*
* @summary Compute the point(s) minimizing a real-valued arbitrary function of several real variables
* defined on the unit simplex, using an exhaustive search algorithm on a grid made of rational points belonging to the unit simplex.
*
* @description This function returns the list of points x = (x_1,...,x_n) belonging to the unit simplex of R^n which 
* minimize an arbitrary real-valued function fct of n real variables defined on the unit simplex of R^n, 
* using an exhaustive search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the reference.

* To be noted that per lemma 1 of the reference, the number of points on such a rational grid is equal to
* factorial(n + k - 1) / (factorial(k - 1) * factorial(n)), i.e., binomial(n+k-1, n-1), 
* so that this method can be of limited use, even for small n.
*
* For instance, n=5 and k=100 already result in 4598126 points on which to evaluate fct.
*
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
*
* @param {function} fct a real-valued function of n real variables defined on the unit simplex of R^n, 
* which must take as first input argument an array of n real numbers corresponding to a point on the unit simplex of R^n 
* and which must return as output a real number.
* @param {number} n the number of variables of the function fct, natural integer superior or equal to 1.
* @param {number} k the indice of the rational grid of the unit simplex of R^n on which to minimize the function fct, a natural integer superior or equal to 1.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to a point of R^n 
* minimizing the function fct on the k-th rational grid of the unit simplex of R^n.
*
*/
function simplexRationalGirdSearch_(fct, n, k) {
	// Initialize the current minimum value and the current list of associated grid points
	var minValue = Number.POSITIVE_INFINITY;
	var minValueGridPoints = [];

	// Proceed to an exhaustive grid search on the set 1/k * I_n(k), c.f. the reference.
	var sampler = new simplexDeterministicRationalSampler_(n, k);
	var weights = sampler.sample();
	while (weights !== null) {  
		// Evaluate the function fct at the current grid point
		var fctValue = fct(weights);
	  
		// If the function value at the current grid point is lower than the current minimum value, this value
		// becomes the new minimum value and the current grid point becomes the new (unique for now) associated grid point.
		if (fctValue < minValue) {
			minValue = fctValue;
			minValueGridPoints = [weights];
		}
		// In case of equality of the function value at the current grid point with the current minimum value, 
		// the current grid point is added to the list of grid points associated to the current minimum value.
		else if (fctValue == minValue) {
			minValueGridPoints.push(weights);
		}
		// Otherwise, nothing needs to be done
		
		// Generate a new grid point
		weights = sampler.sample();
	}
	
	// Return the list of grid points associated to the minimum value of fct
	return minValueGridPoints;
}

/**
 * @file Functions related to cluster risk parity portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


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
* @param {object} opt the optional parameters for the algorithms used by the function.
* @param {number} opt.eps the tolerance parameter for the convergence of the ERC algorithms, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the ERC algorithms, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.clusteringMode the method to use for the clusters computation, a string either equals to:
* - 'ftca': usage of the list of clusters automatically constructed by the FTCA algorithm from David Varadi
* - 'manual': usage of a list of clusters provided in input in the opt.clusters option
* ; defaults to 'ftca'.
* @param {number} opt.ftcaThreshold the correlation threshold to use in the FTCA algorithm in case opt.clusteringMode is equal to 'ftca', a real number beloning to [0,1]; defaults to 0.5.
* @param {Array.<Array.<number>>} opt.clusters the list of clusters to use in the algorithm in case opt.clusteringMode is equal to 'manual', an array of m arrays of strictly positive integers representing the indexes of the assets in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
* @return {Array.<number>} the weights corresponding to the cluster risk parity portfolio, array of n real numbers.
*
* @example
* clusterRiskParityWeights([[0.1,0], [0,0.2]], {clusteringMode: 'manual', clusters: [[1], [2]]});
*  
*/
self.clusterRiskParityWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	
	// Clustering options
	var clusteringMode = opt.clusteringMode || 'ftca';
	
	// Convert sigma to matrix format and convert it to acovariance matrix
	var sigma = new Matrix_(sigma).toCovarianceMatrix(sigma);
	var nbAssets = sigma.nbRows;

	
	// ------
	// Clustering logic
	var clusters = [];
	
	// In case clusters are provided, check that they form a partition of the set [1..nbAssets]
	if (clusteringMode === 'manual') {
		// Decode the manual clustering options
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
		// Extract the correlation matrix from the covariance matrix, as a double array
		var corrMat = sigma.getCorrelationMatrix().toRowArray();
		
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
		
		var clusterSigma = sigma.submatrix(assetsIndexes, assetsIndexes);
		
		// Compute ERC weights for these assets
		var assetsWeights = self.equalRiskContributionWeights(clusterSigma, opt);
		
		// - Populate the change of base matrix for the current cluster with the computed ERC weights
		for (var j = 0; j < assetsWeights.length; ++j) {
			assetsToClustersWeights.setValueAt(i+1, assetsIndexes[j], assetsWeights[j]);
		}
	}
	
	// 2(tmp) - Compute the transpose of the matrix mapping the initial assets to the clusters space
	var assetsToClustersWeightsT = assetsToClustersWeights.transpose()
	
	// 3a - Compute the the covariance matrix associated to the weighted clusters space, using formula Var(A*X) = A * Var(X) * A'.
	var clustersSigma = Matrix_.xy(assetsToClustersWeights, Matrix_.xy(sigma, assetsToClustersWeightsT));
	
	// 3b - Compute ERC weights in the clusters space
	var clustersWeights = self.equalRiskContributionWeights(clustersSigma, opt);
	
	// 3c - Compute original assets weights, using formula A' * Y
	var weights = Matrix_.xy(assetsToClustersWeightsT, new Matrix_(clustersWeights));
	
	// Return them (already normalized)
	return weights.toArray();
}


/**
 * @file Functions related to equal risk bounding portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskBoundingWeights
*
* @summary Compute the weights of the equal risk bounding portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk bounding, defined as the portfolio for which the contribution of each asset 
* included in the portfolio to the risk of the portfolio is equal, c.f. the reference.
*
* A property of the equal risk bounding portfolio is that it is either equal to the equal risk contribution portfolio
* (if all the assets are included in the portfolio), or to a portfolio with a smaller variance than the 
* equal risk contribution portfolio and for which the risk contribution of each asset included in the 
* portfolio is strictly smaller than in the equal risk contribution portfolio (if some assets are not included
* in the portfolio).
*
* Hence, the equal risk bounding portfolio is an equal risk contribution portfolio restricted
* to a (possibly strict) subset of the n assets.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* The algorithm used to solve the associated optimisation problem is an exhaustive computation of all the ERC portfolios
* defined on all the 2^n subsets of the n assets.
* This approach is expected to produce a solution within a reasonable amount of time for small n (e.g. n <= 20),
* but due to the combinatorial nature of the problem, the computation for greater values of n will not be tractable.
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-016-0477-6">Cesarone, F. & Tardella F., Equal Risk Bounding is better than Risk Parity for portfolio selection, J Glob Optim (2017) 68: 439</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the underlying ERC algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the underlying ERC algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the equal risk bounding portfolio, array of n real numbers.
*
* @example
* equalRiskBoundingWeights([[1,-9/10, 3/5], [-9/10, 1,-1/5],[ 3/5, -1/5, 4]]);
* // [0.5, 0.5, 0]
*/
self.equalRiskBoundingWeights = function (sigma, opt) {	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Determine the size of the universe
	var nbAssets = sigma.nbRows;
	
	// Initialize the current minimum value of the risk contribution and the current list of associated assets/assets weights
	var minRCValue = Number.MAX_VALUE;
	var minRCAssetsIndexes = [];
	var minRCAssetsWeights = [];

	// Proceed to an exhaustive enumeration of all the subsets of the set {1,...,nbAssets},
	// in order to find the x^ERB as detailled in section 3.3.1, formula 22, of the reference.
	//
	// The empty set is skipped.
	var nextSubsetIterator = new subsetsIterator_(nbAssets);
	var nextSubset = nextSubsetIterator.next();
	do {
		// Generate a new subset	
		var nextSubset = nextSubsetIterator.next();
		
		// Extract the associated assets indexes
		var assetsIndexes = nextSubset[1];
		
		// Extract the covariance matrix of the associated assets
		var subsetSigma = sigma.submatrix(assetsIndexes, assetsIndexes);
		
		// Compute ERC weights for these assets
		var assetsWeights = self.equalRiskContributionWeights(subsetSigma, opt);
		
		// Compute the risk contribution of the first asset
		//
		// Note: all risk contributions are equal, by definition of the ERC portfolio
		assetsWeightsMatrix = new Matrix_(assetsWeights);
		var rcValue =  assetsWeightsMatrix.getValueAt(1,1) * Matrix_.xy(subsetSigma, assetsWeightsMatrix).getValueAt(1,1) // = x_1 * (Sigma*x)_1
		
		// If the risk contribution of the current subset is lower than the current minimum risk contribution, it
		// becomes the new minimum risk contribution and the current subset becomes the new list of associated assets.
		if (rcValue < minRCValue) {
			minRCValue = rcValue;
			minRCAssetsIndexes = assetsIndexes;
			minRCAssetsWeights = assetsWeights;
		}
		// Otherwise, nothing needs to be done
	}
	while (nextSubset[0]);
	
	// Compute the original assets weights, following formula 22 of the reference
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < minRCAssetsIndexes.length; ++i) {
		weights.setValueAt(minRCAssetsIndexes[i], 1, minRCAssetsWeights[i]);
	}
	
	// Return the computed weights (already normalized)
	return weights.toArray();
}


/**
 * @file Functions related to equal risk budget portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskBudgetWeights
*
* @summary Compute the weights of the equal risk budget portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk budgets, defined as w_i = 1/sigma_i / (sum(1/sigma_j), j=1..n), i=1..n, with sigma_i the standard deviation 
* of the asset i.
*
* This portfolio is unique.
*
* This portfolio maximizes the Sharpe ratio if the assets mean returns are proportional to their volatilities and all pair-wise correlations are equal.
* 
* @see <a href="https://ssrn.com/abstract=1949003">Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description (September 13, 2011). The Journal of Portfolio Management, vol. 38, no. 3, Spring 2012.</a>
* 
* @param {Matrix_|<Array.<number>} sigma the variance vector (sigma_i),i=1..n of the n assets in the considered universe, an n by 1 matrix (i.e., vector) or an array of n real numbers statisfying sigma[i-1] = sigma_i.
* @param {object} opt optional parameters for the algorithm, unused.
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
	var weights = new Matrix_(sigma).elemMap(function(i,j,val) { return 1/Math.sqrt(val); })
	weights = weights.normalize(weights);
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to equal risk contributions portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function equalRiskContributionWeights
*
* @summary Compute the weights of the equal risk contribution portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with equal risk contributions.
*
* This portfolio has the property that the contribution of each asset to the risk of the portfolio is equal,
* and is a special case of the more generic risk budgeting portfolio, with all risk budgets
* equal, c.f. the first reference.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
*
* To be noted that the algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/jpm.2010.36.4.060">Maillard, S., Roncalli, T., Teiletche, J.: The properties of equally weighted risk contribution portfolios. J. Portf. Manag. 36, 60–70 (2010)</a>
* @see <a href="https://arxiv.org/abs/1311.4057">Théophile Griveau-Billion, Jean-Charles Richard, Thierry Roncalli; A Fast Algorithm for Computing High-dimensional Risk Parity Portfolios. eprint arXiv:1311.4057</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the equal risk contribution portfolio, array of n real numbers.
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

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


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
	var weights = Matrix_.ones(nbAssets, 1);
	weights = weights.normalize(weights);

	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function globalMinimumVarianceWeights
*
* @summary Compute the weights of the global minimum variance portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only
* global minimum variance portfolio of n assets.
*
* This portfolio is Markowitz-efficient (i.e., it lies on the Markowitz efficient frontier) and is the portfolio
* with the lowest variance among all the feasible portfolios.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* To be noted that the algorithm used internally is a cyclical coordinate descent, c.f. the reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square Matrix or array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {Array.<number>} the weights corresponding to the global minimum variance portfolio, array of n real numbers.
*
* @example
* globalMinimumVarianceWeights([[0.0400, 0.0100], [0.0100, 0.0100]], {eps: 1e-10, maxIter: 10000});
* // XX
*/
self.globalMinimumVarianceWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-8;
	var maxIterations = opt.maxIter || 10000;
	
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1);
	x = x.normalize(x);
	
	// Preparational computations
	var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
	var sigma_x_2 = Matrix_.xy(sigma, x); // SIGMA*x
	var var_x = Matrix_.vectorDotProduct(sigma_x, x); // sigma(x), the portfolio variance
	var obj = 0.5 * var_x - x.sum();
	
	// Main loop until convergence, guaranteed as per hypotheses on sigma
	var iter = 0;
	var converged = false;
	while (!converged) {
        // By default, let's assume the new iteration will lead to the convergence of the algorithm
    	converged = true;
    	
    	for (var i = 1; i <= nbAssets; ++i) {		
			// Define the coefficients of the second order polynomial ax_i^2 + b_ix + c, c.f. the reference
    	    var a = sigma.getValueAt(i,i); // sigma_i^2, always > 0
    	    var b = sigma_x.getValueAt(i,1) - x.getValueAt(i,1) * sigma.getValueAt(i,i) - 1; // (SIGMA*x)_i - x_i*sigma_i^2 - 1, might be any sign
    	    var c = 0;
    	    
    	    // Note: what follows is not detailled in the reference, as lambda_erc is supposed there to be non zero.
    	    //
    	    // The equation ax_i^2 + bx_i + c = 0 has two roots: 0 and another one, possibly strictly negative, 0 or strictly positive:
			// Case #1 - If the other root is 0 or strictly negative, x_i^* is then equal to 0
			// Case #2 - If the other root is strictly positive, x_i^* is equal to the value minimizing the objective function in 1D

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
    	    sigma_x = Matrix_.xy(sigma, x, sigma_x);
    	    var_x = Matrix_.vectorDotProduct(sigma_x, x);
			
			  // Case #1 objective function computation
			var obj_star = 0.5 * var_x - x.sum();
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				sigma_x_2 = Matrix_.xy(sigma, x, sigma_x_2);
				var var_x_2 = Matrix_.vectorDotProduct(sigma_x_2, x);
				
				// Case #2 objective function computation
				var obj_star_2 = 0.5 * var_x_2 - x.sum();
				
				// Selection between the two cases
				// Note: Portfolio sparsity is priviledged in case of equality
				if (obj_star_2 < obj_star) {
					// Case #2 weights update is already done
					
					// Update the products for convergence condition evaluation + next loop evaluation
					sigma_x = Matrix_.copy(sigma_x_2, sigma_x);
					var_x = var_x_2;
					obj_star = obj_star_2;
				}
				else {
					// Case #1 weights re-update, as erased by case #2 test
					x.setValueAt(i, 1, xi_star_1);
					
					// No update on the products, as case #1 is the default case
				}
			}
        }
		
		// Update the convergence condition: |obj - obj*| <= eps
		if (Math.abs(obj_star - obj) > eps) {
			converged = false;
			obj = obj_star;
		}

		// Update the number of iterations
		++iter;
		
		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}
	
	// Return the computed weights, after normalization
	x = x.normalize(x);
	return x.toArray();
}


/**
 * @file Functions related to grid search portfolios.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function gridSearchWeights
*
* @summary Compute the weights of a portfolio minimizing an arbitrary objective function.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only portfolio
* of n assets minimizing an arbitrary real-valued objective function fct of n real variables defined on the unit simplex of R^n, 
* which is the n-1 dimensional set of R^n containing the points x = (x_1,...,x_n) statisfying sum x_i = 1 and x_i >= 0, i = 1..n, 
* c.f. the first reference.
*
* Since such a portfolio might not be unique, all the weights corresponding to the same minimum value of the function fct 
* are provided in output.
*
* The minimization of the function fct is achieved through one of the following optimisation methods:
* - Deterministic search on a grid of rational points belonging to the unit simplex of R^n
*
* To be noted that finding the minimum value(s) of an arbitrary objective function on the unit simplex
* is an NP-hard problem, c.f. the second reference, so that all the optimisation algorithms above are expected to
* be non-polynomial in n.
* 
* @see <a href="https://en.wikipedia.org/wiki/Simplex">Simplex</a>
* @see <a href="http://www.sciencedirect.com/science/article/pii/S0377221707004262">De Klerk, E., Den Hertog, D., Elabwabi, G.: On the complexity of optimization over the standard simplex. Eur. J Oper. Res. 191, 773–785 (2008)</a>
* @see <a href="https://ideas.repec.org/p/cor/louvco/2003071.html">Nesterov, Yurii. Random walk in a simplex and quadratic optimization over convex polytopes. CORE Discussion Papers ; 2003/71 (2003)</a>
*
* @param {number} nbAssets the number of assets in the considered universe, natural integer superior or equal to 1.
* @param {function} fct a real-valued objective function of n real variables to minimize on the unit simplex of R^n,
* which must take as first input argument an array of n real numbers corresponding to the weights w1,...,wn of the n assets 
* in the considered universe and which must return as output a real number.
* @param {Object} opt the optional parameters for the algorithms used by the function.
* @param {string} opt.optimisationMethod the optimisation method to use in order to minimize the function fct, a string either equals to:
* - 'deterministic': usage of a deterministic grid search algorithm on the k-th rational grid of the unit simplex of R^n, 1/k * I_n(k), c.f. the third reference, where k is defined through the parameter opt.rationalGrid.k
* -
* ; defaults to 'deterministic'.
* @param {number} opt.rationalGrid.k the indice k of the k-th rational grid of the unit simplex of R^n to use in case opt.optimisationMethod is equal to 'deterministic', a natural integer greater than or equal to 1; defaults to n.
* @return {Array.<Array.<number>>} an array of possibly several arrays of n real numbers, each array of n real numbers corresponding to
* the weights of a portfolio minimizing the function fct.
*
* @example
* gridSearchWeights(3, function(arr) { return arr[0];}, {optimisationMethod: 'deterministic', rationalGrid: {k: 2}});
* // [[0,1,0],[0,0.5,0.5],[0,0,1]]
*/
self.gridSearchWeights = function (nbAssets, fct, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	opt.optimisationMethod = opt.optimisationMethod || 'deterministic';
	
	// Select the proper optimisation method
	if (opt.optimisationMethod === 'deterministic') {
		// Decode options for the deterministic grid search method
		opt.rationalGrid = opt.rationalGrid || { k: nbAssets };
		var k = opt.rationalGrid.k;
		
		// Call the rational grid search method
		return simplexRationalGirdSearch_(fct, nbAssets, k);
	}
	else {
	    throw new Error('unsupported optimisation method');
	}
}



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
* @param {boolean} opt.contraints.partialInvestment parameter set to true in case the full investment constraint of the portfolio must be replaced
* by a partial investment constraint; defaults to false.
* @return {Array.<number>} the weights corresponding to a minimax portfolio, array of real numbers of length n.
*
* @example
* minimaxWeights([[0.01, -0.02, 0.01], [-0.05, 0.03, 0.01]]);
* // ~[ 0.727, 0.273 ]
*/
self.minimaxWeights = function (assetsReturns, opt) {
	// TODO: Checks, if enabled

	// Decode options
	if (opt === undefined) {
		opt = { contraints: {} };
	}
	var partialInvestmentContraint = false;
	if (opt.contraints.partialInvestment !== undefined) {
		partialInvestmentContraint = opt.contraints.partialInvestment;
	}
	
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
	
		// Build the equality constraints (c.f. formula 1d of the first reference for the inequality equivalent):
		// - Full investment (optional)
	var Ae = null;
	var be = null;
	if (partialInvestmentContraint === false) {
		Ae = Matrix_.fill(1, nbAssets + 1, function(i,j) { return j <= nbAssets ? 1 : 0; }); // Ae = [1,...,1,0]
		be = Matrix_.ones(1,1); // be = [1]
	}
	
		// Build the inequality constraints (c.f. formula 1b of the first reference):
		// - Portfolio return greater than or equal to the minimum portfolio return, for each period
		// - Partial investment (optional)
	var Ai = Matrix_.fill(nbPeriods + (partialInvestmentContraint ? 1 : 0), nbAssets + 1, 
						  function(i,j) { 
								if (i <= nbPeriods) { return j <= nbAssets ? -assetsReturns[j-1][i-1] : 1; }
								else if (i == nbPeriods + 1) { return j <= nbAssets ? 1 : 0; }
						  }); // Ai = [[-ret11,...,-retN1,1]...[-ret1T,...,-retNT,1] (optional: , [1,...,1,0])]
	var bi = Matrix_.fill(nbPeriods + (partialInvestmentContraint ? 1 : 0), 1, 
						  function(i,j) { 
							  if (i <= nbPeriods) { return 0; }
							  else if (i == nbPeriods + 1) { return 1; }
						  }); // bi = [0, ..., 0 (optional: , 1)]
	
		// Build the bound constraints (c.f. formula 1e of the first reference + the definition of M_p):
		// - No short sales
		// - Absence of leverage
		// - "Unbounded" minimum portfolio return
	var lb = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 0 : -Infinity; }); // lb = [0,...,0, -Infinity]
	var ub = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 1 : Infinity; });  // ub = [1,...,1, Infinity]
	
		// Solve the constructed linear program, which is:
		// - Bounded: the portfolio weights belong to the unit simplex
		//            the minimum portfolio return is bounded below by the minimum portfolio return over all the periods
		// - Feasible: the portfolio with the minimum return over all the periods is a solution to the linear program
		//
		// Note: given the assumptions above, the convergence of the primal-dual hybrid gradient algorithm is guaranteed.
	var lpSolution = lpsolvePrimalDualHybridGradient_(Ae, be, Ai, bi, c, lb, ub, {maxIter: -1});
	
	
	// ----
	
	// Extract the computed portfolio weights.
	var weights = Matrix_.fill(nbAssets, 1, function(i,j) {  return lpSolution[0].getValueAt(i, 1); });

	// Return the computed weights.
	return weights.toArray();
}

/**
 * @file Functions related to minimum correlation (heuristic) portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


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
* This portfolio is meant to approximate the most diversified portfolio (MDP).
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
	var invStddevs = variances.elemMap(function(i,j,val) { return 1/Math.sqrt(val); });
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
	var rowsElements = adjustedRho.toRowArray(function(i, j, val) {
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
	weights = Matrix_.xy(adjustedRho, weights).normalize();
	
	// Step 7: Scale portfolio weights by assets standard deviations
	weights = Matrix_.vectorHadamardProduct(weights, invStddevs).normalize();
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to most diversified portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


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
* The algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
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
	var stddevs = sigma.diagonal().elemMap(function(i,j,val) { return Math.sqrt(val); });

	// TODO: Checks, if enabled
	// Check that diagonal entries of sigma are strictly positive
	// Check that sigma is symmetric and positive definite
	// Check that sigma and rb are rows compatible


	// ------
	var nbAssets = sigma.nbRows;
	
	// Initial point for the algorithm is an equal weight vector
	var x = Matrix_.ones(nbAssets, 1);
	x = x.normalize(x);
	
	// Preparational computations
	var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
	var sigma_x_2 = Matrix_.xy(sigma, x); // SIGMA*x
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
			// Case #2 - If the other root is strictly positive, x_i^* is equal to the value maximizing the diversification ratio in 1D

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
    	    sigma_x = Matrix_.xy(sigma, x, sigma_x)
    	    s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x));
			
			  // Case #1 diversification ratio computation
			var dr_star = Matrix_.vectorDotProduct(stddevs, x) / s_x; // Weighted average volatility of the portfolio divided by the volatility of the portfolio
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				sigma_x_2 = Matrix_.xy(sigma, x, sigma_x_2);
				var s_x_2 = Math.sqrt(Matrix_.vectorDotProduct(sigma_x_2, x));
				
				// Case #2 diversification ratio computation
				var dr_star_2 = Matrix_.vectorDotProduct(stddevs, x) / s_x_2;
				
				// Selection between the two cases
				// Note: Portfolio sparsity is priviledged in case of equality
				if (dr_star_2 > dr_star) {
					// Case #2 weights update is already done
					
					// Update the products for convergence condition evaluation + next loop evaluation
					sigma_x = Matrix_.copy(sigma_x_2, sigma_x);
					s_x = s_x_2;
					dr_star = dr_star_2;
				}
				else {
					// Case #1 weights re-update, as erased by case #2 test
					x.setValueAt(i, 1, xi_star_1);
					
					// No update on the products, as case #1 is the default case
				}
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
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}
	
	// Return the computed weights, after normalization
	x = x.normalize(x);
	return x.toArray();
}


/**
 * @file Functions related to proportional minimum variance (heuristic) portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


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
* This portfolio is meant to approximate the global minimum variance (GMV) portfolio.
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
	var rowsElements = sigma.toRowArray();
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	var elementsMean = mean_(rowsAverages);
	var elementsStddev = sampleStddev_(rowsAverages);

	// Step 2: Gaussian convertion, and proportional average covar weigth
	var weights = new Matrix_(rowsAverages, 1).elemMap(function(i, j, val) { 
		return 1 - normcdf_((val - elementsMean)/elementsStddev);
	});
	weights = weights.normalize(weights);
	
	// Step 3: Scale portfolio weights by assets variances
	var invVariancesWeights = sigma.diagonal().elemMap(function(i,j,val) { return 1/val; }).normalize();
	weights = Matrix_.vectorHadamardProduct(weights, invVariancesWeights).normalize();
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to random weights portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function randomWeights
*
* @summary Compute the weights of a randomly generated portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to a fully invested and long-only 
* random portfolio of n assets.
*
* Optionally, the following constraints can be added:
* - Minimum number of assets to include in the portfolio
* - Maximum number of assets to include in the portfolio
*
* This portfolio is not unique.
*
* Random portfolios have several applications in asset allocation, c.f. the first reference, as
* well as in trading strategies evaluation, c.f. the second reference.
*
* To be noted that the algorithms used internally allow the random portfolios to be generated uniformly
* among all the feasible portfolios.
*
* @see <a href="https://arxiv.org/abs/1008.3718">William T. Shaw, Monte Carlo Portfolio Optimization for General Investor Risk-Return Objectives and Arbitrary Return Distributions: a Solution for Long-only Portfolios</a>
* @see <a href="https://ssrn.com/abstract=881735">Burns, Patrick, Random Portfolios for Evaluating Trading Strategies (January 13, 2006)</a>
*
* @param {number} nbAssets the number of assets in the universe, natural integer superior or equal to 1.
* @param {object} opt optional parameters for the algorithm.
* @param {number} opt.contraints.minAssets the minimum number of assets to include in the portfolio, an integer i satisfying 1 <= i <= nbAssets; defaults to 1.
* @param {number} opt.contraints.maxAssets the maximum number of assets to include in the portfolio, an integer j satisfying i <= j <= nbAssets; defaults to nbAssets.
* @return {Array.<number>} the weights corresponding to a random portfolio, array of real numbers of length nbAssets.
*
* @example
* randomWeights(5);
* // ~[0, 0 0.33, 0.33, 0.33]
*/
self.randomWeights = function (nbAssets, opt) {
	// TODO: Checks, if enabled

	// Decode options
	if (opt === undefined) {
		opt = { contraints: {} };
	}
	var nbMinAssets = opt.contraints.minAssets || 1;
	var nbMaxAssets = opt.contraints.maxAssets || nbAssets;
	
	// 1 - Generate the number of assets to include in the portfolio (uniform generation)
	var nbSelectedAssets = Math.floor(Math.random() * (nbMaxAssets - nbMinAssets +1)) + nbMinAssets;
	
	// 2 - Generate the indices of the assets to include in the portfolio (uniform generation)
	var selectedAssetsIdx = new randomKSubsetIterator_(nbAssets, nbSelectedAssets).next();
	
	// 3 - Generate the weights of the assets to include in the portfolio (uniform generation)
	//
	// Extra caution needs to be taken in case one of the weights is zero,
	// as exactly nbSelectedAssets must be included in the portfolio.
	var selectedAssetsWeights;
	var simplexSampler = new simplexRandomSampler_(nbSelectedAssets);
	while (true) {
		// Generate a sample of assets weights
		selectedAssetsWeights = simplexSampler.sample();
		
		// Reject the sample if there is one asset with a zero weight
		var rejectSample = false;
		for (var i = 0; i < nbSelectedAssets; ++i) {
			if (selectedAssetsWeights[i] == 0) {
				rejectSample = true;
				break;
			}
		}
		if (!rejectSample) {
			break;
		}
	}
	
	// Compute the final weights vector:
	// - The weights associated to assets not included in the portfolio at step 2 are set to zero
	// - The weights associated to assets included in the portfolio at step 2 are set to their values generated at step 3
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < nbSelectedAssets; ++i) {
		// Extract included assets information
		var assetIdx = selectedAssetsIdx[i];
		var assetWeight = selectedAssetsWeights[i];
		
		// Update the weights vector
		weights.setValueAt(assetIdx, 1, assetWeight);
	}

	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to risk budgeting portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function riskBudgetingWeights
*
* @summary Compute the weights of the risk budgeting portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only portfolio
* of n assets with risk budgeting contraints.
*
* This portfolio has the property that the total contribution 
* of each asset to the risk of the portfolio is equal to a pre-determined budget weight.
*
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* To be noted that the algorithm used internally is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
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
	var x = Matrix_.ones(nbAssets, 1);
	x = x.normalize(x);
	
	// Preparational computations
	var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
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
    	    sigma_x = Matrix_.xy(sigma, x, sigma_x)
    	    s_x = Math.sqrt(Matrix_.vectorDotProduct(sigma_x, x));
    	    
    	    // Update the convergence condition: |RC_i* - b_i| <= eps, i = 1..nbAssets
    	    var rci_star = x.getValueAt(i,1) * sigma_x.getValueAt(i,1) / s_x;
    	    if (Math.abs(rci_star - rb.getValueAt(i,1)) > eps) {
    	        converged = false;
    	    }
        }
		
		// Update the number of iterations
		++iter;
		
		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}
	
	// Return the computed weights, after normalization
	x = x.normalize(x);
	return x.toArray();
}


/**
* @file Functions related to (rational) rounding of floating-point portfolio weights.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function roundedWeights
*
* @summary Compute the closest rational approximation of the weights of a portfolio.
*
* @description Given n (floating-point) weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, 
* this function returns n (rational) weights wr_1,...,wr_n associated to a fully invested and long-only portfolio of n assets
* satisfying:
* - k * wr_i is a natural integer, i=1..n
* - wr_1,...,wr_n are the closest weights to w_1,...,w_n, in the sense defined in the reference.
*
* To be noted that typical values of k are 10 (rounding to 10%), 20 (rounding to 5%) and 100 (rounding to 1%).
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-013-0126-2">.M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258.</a>
* 
* @param {Array.<number>} originalWeights the weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, array of n real numbers.
* @param {number} k the value to which the rounded weights will be a multiple of the inverse, natural integer greater than or equal to 1.
* @return {Array.<number>} the rounded weights wr_1,...,wr_n, array of n real numbers.
*
* @example
* roundedWeights([0.5759, 0.0671, 0.3570], 10);
* // [0.6, 0.1, 0.3]
* roundedWeights([0.5759, 0.0671, 0.3570], 20);
* // [0.6, 0.05, 0.35]
* roundedWeights([0.5759, 0.0671, 0.3570], 100);
* // [0.57, 0.07, 0.36]
*/
self.roundedWeights = function (originalWeights, k) {
	// ------
	
	// Call to the simplex rational rounding method
	var roundedWeights = simplexRationalRounding_(originalWeights, k);

	// Return the computed weights
	return roundedWeights;	
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
