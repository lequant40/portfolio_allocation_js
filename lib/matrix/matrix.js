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
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j] / sum; 
			}
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
		var r = Matrix_.qrDecomposition(this, {qLess: true});
		
		// By property of the Givens QR decomposition, the determinant of the matrix
		// is then the product of the diagonal elements of the R matrix.
		var det = 1.0;
		for (var i = 0; i < r.nbRows; ++i) {
		        det *= r.data[i * r.nbColumns + i];
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
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					absSum += Math.abs(this.data[i * this.nbColumns + j]);
				}
			}
			return absSum;
		}
		else if (p == 'two') {
			// Compute the l^2 vector norm using an accurate algorithm by S. J. Hammarling
			// C.f. problem 27.5 of the second reference
			var t = 0;
			var s = 1;
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    for (var j = idxStartColumn; j < idxEndColumn; ++j) {

				    var val = this.data[i * this.nbColumns + j];
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
				}
			}
			return t * Math.sqrt(s);
		}
		else if (p == 'infinity') {
			// Compute the largest absolute value of the elements of the matrix
		    var maxAbsVal = 0;
			for (var i = idxStartRow; i < idxEndRow; ++i) {
    		    for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					 maxAbsVal = Math.max(maxAbsVal, Math.abs(this.data[i * this.nbColumns + j]));
				}
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
	var tol = eps || 0;
	for (var i = 0; i < nbRows; ++i) {
		for (var j = 0; j < nbColumns; ++j) {
			if (Math.abs(a.data[i * nbColumns + j] - b.data[i * nbColumns + j]) > tol) {
				return false;
			}
		}
	}

	//
	return true;
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
* @param {Matrix_} out an optional n by m matrix.
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = a * X.data[i * X.nbColumns + j] + b*Y.data[i * Y.nbColumns + j];
		}
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = A.data[i * A.nbColumns + j];
		}
	}
	
	// Return the computed matrix
	return obj;
};



/**
* @function product
*
* @summary Returns the product of a matrix with another matrix.
*
* @description This function computes the product A*B of a matrix A with another matrix B,
* where A is a n by m matrix and B is a m by p matrix.
*
* @param {Matrix_} A a n by m matrix.
* @param {Matrix_} B a m by p matrix.
* @param {Matrix_} out an optional n by p matrix.
* @return {Matrix_} the matrix product A*B, either stored in the matrix out or in a new matrix, a n by p matrix.
*
* @example
* product(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.product = function(a, b, out) {
	// Delegate the computations to the internal function.
	return Matrix_.axy(1, a, b);
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
	if (X.nbRows === Y.nbRows && X.nbColumns === Y.nbColumns) { // full matrix elementwise product
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = X.data[i * X.nbColumns + j] * Y.data[i * Y.nbColumns + j];
			}
		}
	}
	else if (Y.nbRows === X.nbRows && Y.nbColumns === 1) { // row matrix elementwise product
		for (var i = 0; i < obj.nbRows; ++i) {
			var y_i = Y.data[i];
			
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = X.data[i * X.nbColumns + j] * y_i;
			}
		}
	}
	else if (Y.nbRows === 1 && Y.nbColumns === X.nbColumns) { // column matrix elementwise product
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = X.data[i * X.nbColumns + j] * Y.data[j];
			}
		}
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

	// Computation of a*X*Y product in IKJ format, cache aware
	for (var i = 0; i < X.nbRows; ++i) {
		for (var j = 0; j < Y.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
	
		for (var k = 0; k < X.nbColumns; ++k) {
			var a_ik = a * X.data[i * X.nbColumns + k];
			
			for (var j = 0; j < Y.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] += a_ik * Y.data[k * Y.nbColumns + j];
			}
		}
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
* atxy(Matrix_(1, [[1,4], [2,5], [3,6]]), Matrix_([[1,1], [2,2], [3,3]]));
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
	var obj = allocateMatrix_(x.nbRows, x.nbRows, out);
	
	// Result matrix computation
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
		obj.data[i * obj.nbColumns + i] = x.data[i * x.nbColumns];
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < i; ++j) {
			obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
		}
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = fct(i+1 ,j+1);
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = fct(i+1 ,j+1);
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 1;
		}
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
	for (var i = 0; i < obj.nbRows; ++i) {
		for (var j = 0; j < obj.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
		obj.data[i * obj.nbColumns + i] = 1;
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
	for (var i = 0; i < x.nbRows; ++i) {
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
* @description This function computes a QR decomposition of a matrix A, using Givens rotations 
* as described in the algorithm 5.2.4 (and above discussion therein) of the first reference.
* 
* To be noted that the 'givens' internal function used in this function is not the same
* as the one described in the first reference, but is the continuous one described
* in the algorithm 4 of the second reference.
* 
* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
* @see <a href="http://www.netlib.org/lapack/lawnspdf/lawn150.pdf">Anderson, Edward (4 December 2000). Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem. LAPACK Working Note. University of Tennessee at Knoxville and Oak Ridge National Laboratory.</a>
*
* @param {Matrix_} A an m by n matrix, with m >= n.
* @param {object} opt the optional parameters for the algorithm.
* @param {boolean} opt.qLess a boolean parameter enabling to discard the computation of the Q matrix; defaults to false.
* @return {<Array.<Matrix_>|Matrix_} either an R matrix if opt.qLess is set to true or an array of two matrices [Q, R] otherwise, 
* with Q and R satisfying the following properties:
* - Q is an m by m matrix 
* - R is an m by n matrix
* - A = Q*R
* - Q is an orthogonal matrix, with a determinant equals to 1 (i.e., a rotation)
* - R is an upper triangular matrix, with its bottom (m−n) rows consisting entirely of zeroes
*
* @example
* qrDecomposition(Matrix_([[1],[2],[3]]), Matrix_([[1],[2],[3]]));
* // XX
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
    * 
    * @example
    * givens(X, y);
    * // XX
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
	var rr = new Matrix_(A); // represents R
	
	// Create the matrix that will hold Q
	var qq = null;
	if (!isQLessDecomposition) {
		qq = Matrix_.identity(m); // represents Q
	}
	
	// Core of the algorithm
	for (var j = 1; j <= n; ++j) {
		for (var i = m; i >= j+1; --i) {
			// This loop iteration will set R(i,j) to 0 using an appropriate Givens rotation
			
			// Builds G(i-1, i, theta)	
			var a = rr.data[(i-2) * rr.nbColumns + (j-1)]; // R(i-1, j)
			var b = rr.data[(i-1) * rr.nbColumns + (j-1)]; // R(i, j)
			var givensArr = givens(a, b);
			var c = givensArr[0];
			var s = givensArr[1];
			var r = givensArr[2];
			
			// Update R (left multiply R with the transpose of the Givens rotation, c.f. 5.1.9 of the first reference)		
				// Loop below unrolled for k = j, since R(i,j) must be made zero (this can help avoiding "fake" non-zeroes)
			rr.data[(i-2) * rr.nbColumns + (j-1)] = r; // R(i-1, j) = r
			rr.data[(i-1) * rr.nbColumns + (j-1)] = 0; // R(i, j) = 0

				// Loop resumed at k = j+1
			for (var k = j+1; k <= n; ++k) {
				var t1 = rr.data[(i-2) * rr.nbColumns + (k-1)]; // t1 = R(i-1, k)
				var t2 = rr.data[(i-1) * rr.nbColumns + (k-1)]; // t2 = R(i, k)
				rr.data[(i-2) * rr.nbColumns + (k-1)] = c*t1 - s*t2; // R(i-1, k) = ...
				rr.data[(i-1) * rr.nbColumns + (k-1)] = s*t1 + c*t2; // R(i, k) = ...
			}
			
			// Update Q (right multiply Q with the Givens rotation, c.f. 5.1.9 of the first reference)
			if (!isQLessDecomposition) {
				for (var k = 1; k <= m; ++k) {
					var t1 = qq.data[(k-1) * qq.nbColumns + (i-2)] // t1 = Q(k,i-1)
					var t2 = qq.data[(k-1) * qq.nbColumns + (i-1)] // t2 = Q(k,i)
					qq.data[(k-1) * qq.nbColumns + (i-2)] = c*t1 - s*t2; // Q(k,i-1) = ...
					qq.data[(k-1) * qq.nbColumns + (i-1)] = s*t1 + c*t2; // Q(k,i) = ...
				}
			}
		}
	}
	
	// Return either the computed [Q, R] pair or the matrix R
	if (!isQLessDecomposition) {
		return [qq, rr];
	}
	else {
		return rr;
	}
}




//James Demmel and Kresimir Veselic, Jacobi’s Method is More Accurate than QR, SIAM Journal on Matrix Analysis and Applications, 1992, Vol. 13, No. 4 : pp. 1204-1245
//https://doi.org/10.1137/0613074
// Algorithm 4.1
// Assume m >= n; otherwise, trasnpose A and retranspose before the end
// int gsl_linalg_SV_decomp_jacobi
/*
A general rectangular M-by-N matrix A has a singular value decomposition (SVD) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,

A = U S V^T

m >= n, golub 2.4.3

The singulmar values \sigma_i = S_{ii} are all non-negative and are generally chosen to form a non-increasing sequence \sigma_1 >= \sigma_2 >= ... >= \sigma_N >= 0.

The singular value decomposition of a matrix has many practical uses. The condition number of the matrix is given by the ratio of the largest singular value to the smallest singular value. The presence of a zero singular value indicates that the matrix is singular. The number of non-zero singular values indicates the rank of the matrix. In practice singular value decomposition of a rank-deficient matrix will not produce exact zeroes for singular values, due to finite numerical precision. Small singular values should be edited by choosing a suitable tolerance.

For a rank-deficient matrix, the null space of A is given by the columns of V corresponding to the zero singular values. Similarly, the range of A is given by columns of U corresponding to the non-zero singular values.

Note that the routines here compute the “thin” version of the SVD with U as M-by-N orthogonal matrix. This allows in-place computation and is the most commonly-used form in practice. Mathematically, the “full” SVD is defined with U as an M-by-M orthogonal matrix and S as an M-by-N diagonal matrix (with additional rows of zeros). 
*/


//Z = null(A) is an orthonormal basis for the null space of A obtained from the singular value decomposition. That is, A*Z has negligible elements, size(Z,2) is the nullity of A, and Z'*Z = I.
/*
There are a number of ways to compute the rank of a matrix. MATLAB® software uses the method based on the singular value decomposition, or SVD. The SVD algorithm is the most time consuming, but also the most reliable.

The rank algorithm is

s = svd(A);
tol = Math.max(m,n) * eps

r = sum(s > tol);
2.5.2 for null subspace projection

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
		throw new Error('matrix has more columns than rows: ' + '(' + A.nbRows + ') v.s. ' + '(' + A.nbColumns + ')');
	}
	
	
	// ------
	
	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	
	// Create a copy of A so that it is not overwritten
	var uu = new Matrix_(A); // represents U 
	
	// Create the matrix that will hold V
	var vv = Matrix_.identity(n); // represents V
	
	
	// ------
	
	// Core of the algorithm, guaranteed to converge per theorem 4.9 of the reference
	/*var iter = 0;
	while (true) {
		// Update the number of iterations (number of sweeps)
		++iter;

		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
			
		// For all pairs (i, j) with i < j, compute one sweep: n*(n-1)/2 Jacobi rotations
		var converged = true; // convergence is assumed to hold for this sweep
		for (var j = 2; j <= n; ++j) {
			for (var i = 1; i <= j-1; ++i) {
				// Compute the (i, j) 2x2 submatrix [[aa, cc], [cc, bb]] of U^T*U
				var aa = 0;
				var bb = 0;
				var cc = 0;
			    for (var k = 1; k <= m; ++k) {
					aa += uu.data[(k-1) * uu.nbColumns + (i-1)] * uu.data[(k-1) * uu.nbColumns + (i-1)]; // aa = sum U(k,i)^2, k=1..m
					bb += uu.data[(k-1) * uu.nbColumns + (j-1)] * uu.data[(k-1) * uu.nbColumns + (j-1)] ; // bb = sum U(k,j)^2, k=1..m
					cc += uu.data[(k-1) * uu.nbColumns + (i-1)] * uu.data[(k-1) * uu.nbColumns + (j-1)]; // cc = sum U(k,i)*U(k,j), k=1..m
				}

				// Compute the Jacobi rotation which diagonalizes the 2x2 submatrix above
				if (cc == 0) { // in this case, no update is performed as zeta below would be infinite, then, t would be 0 and then cs = 1, sn = 0
					continue;
				}
				var zeta = (bb - aa)/(2 * cc);
				var t = ((zeta >= 0) ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta)); // first part emulates sign(zeta)
				var cs = 1 / Math.sqrt(1 + t*t);
				var sn = cs * t;
				
				// Update columns i and j of G (right multiply G with the Jacobi rotation)
				for (var k = 1; k <= m; ++k) {
					var t1 = uu.data[(k-1) * uu.nbColumns + (i-1)] // t1 = U(k,i)
					var t2 = uu.data[(k-1) * uu.nbColumns + (j-1)] // t2 = U(k,j)
					uu.data[(k-1) * uu.nbColumns + (i-1)] = cs*t1 - sn*t2 // U(k,i) = ...
					uu.data[(k-1) * uu.nbColumns + (j-1)] = sn*t1 + cs*t2 // U(k,j) = ...
				}
				
				// Update the matrix V of right singular vectors (right multiply V with the Jacobi rotation)
				for (var k = 1; k <= n; ++k) {
					var t1 = vv.data[(k-1) * vv.nbColumns + (i-1)] // t1 = V(k,i)
					var t2 = vv.data[(k-1) * vv.nbColumns + (j-1)] // t2 = V(k,j)
					vv.data[(k-1) * vv.nbColumns + (i-1)] = cs*t1 - sn*t2 // V(k,i) = ...
					vv.data[(k-1) * vv.nbColumns + (j-1)] = sn*t1 + cs*t2 // V(k,j) = ...
				}
				
				// Update the convergence criteria
				if (Math.abs(cc) > eps * Math.sqrt(aa * bb)) {
					converged = false;
				}
			}
		}
		
		// In case the covnergence criteria is still true at the end of the sweep, 
		// the algorithm can be stopped
		if (converged == true) {
			break;
		}
	}
	*/
	var iter = 0;
	var u_frob_norm = uu.matrixNorm('frobenius'); 
	var u_columns_two_norm_sq = new Array(n);
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
	var sigmas = new Array(n);
	for (var j = 1; j <= n; ++j) {
		// Compute sigma_j
		var singVal_j = uu.vectorNorm('two', 'column', j);
		sigmas[j-1] = [singVal_j, j];
	}
	sigmas.sort(function(a, b) { return b[0] - a[0]; });
	
	// Compute the thin U,S and V matrices
	var uuu = Matrix_.zeros(m, n); 
	var sss = Matrix_.zeros(n, n);
	var vvv = Matrix_.zeros(n, n);
	for (var j = 1; j <= n; ++j) {
		// Extract singluar values information
		var sigma_j = sigmas[j-1][0];
		var sigma_j_col_idx = sigmas[j-1][1];

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
	
	// Additional work is needed in order to achieve the full singular value decomposition of A with:
    // - A = U*S*V^t
	// - U is an orthogonal square m by m matrix
	// - S a rectangular m by n matrix with the n singular values of A on its diagonal 
	//   verifying S(i,i) = sigma_i and sigma_1 >= ... >= sigma_n >= 0
	// - V is an orthogonal square n by n matrix
	
	// Full U computation is made through a QR decomposition of U and completing the missing
	// columns of U by those of Q
	//
	// Full S computation is made by extending the thin S with zeroes
	var qr = Matrix_.qrDecomposition(uuu);
	var q = qr[0];
	
	var uuuu = Matrix_.zeros(m, m);
	var ssss = Matrix_.zeros(m, n);
	for (var j = 1; j <= n; ++j) {
		// Extract singluar values information
		var sigma_j = sigmas[j-1][0];
		
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


//* @see G.H. Golub and C.F. Van Loan, Matrix Computations, 4th Edition, Johns Hopkins Studies in the Mathematical Sciences
Matrix_.nullSpace = function(A, opt) {
    // ------
   
    // Decode options
    if (opt === undefined) {
        opt = {};
    }
    var eps = opt.eps || 1e-16;

   
    // ------
   
    // Misc. checks
    if (!(A instanceof Matrix_)) {
        throw new Error('first input must be a matrix');
    }   


    // ------
   
    // Initializations
    var m = A.nbRows;
    var n = A.nbColumns;
    var a_ns;
   
    // ------
   
    // Depending on the shape of the matrix A, a different algorithm is used:
    // - m >= n: an orthogonal basis of Ker(A) is directly read from a thin SVD
    // decomposition of the matrix A, c.f. corollary 2.4.6 of the first reference
    //
    // - m < n: an orthogonal basis of Ker(A) is obtained by computing the orthogonal
    // complement of an orthogonal basis of Range(A^t), this orthogonal basis being directly
    // read from a thin SVD decomposition of A^t, c.f. corollary 2.4.6 of the first reference
    //
    // The underlying properties which justifies this approach is that
    // R^m = Range(A^t) _|_ Ker(A) for any m by n real matrix A
    if (m >= n) {
        // Compute a thin SVD decomposition of A: A = U*S*V^t
		var svd = Matrix_.svdDecomposition(A, {maxIter: -1});
        var u = svd[0];
        var s = svd[1];
        var v = svd[2];
       
        // Determine the numerical rank of A
        var r;
        for (r = n; r >= 1; --r) {
            // Comparison of the r-th greatest singular value of A with the tolerance parameter to determine
            // the first "zero" singular value of A
            if (s.data[(r-1)*s.nbColumns + (r-1)] > m * eps) {
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
            a_ns = new Matrix_.zeros(n, n-p);
            for (var j = p+1; j <= n; ++j) {
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
       
        // Determine the numerical rank of A^t
        var r;
        for (r = m; r >= 1; --r) {
            // Comparison of the r-th greatest singular value of A^t with the tolerance parameter to determine
            // the first "zero" singular value of A^t
            if (s.data[(r-1)*s.nbColumns + (r-1)] > n * eps) {
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
* @function linsolveKaczmarz
*
* @summary TODO
*
* @description TODO
* 
* @see <a href="https://doi.org/10.1016/j.matcom.2004.01.021">Constantin Popa and Rafal Zdunek. 2004. Kaczmarz extended algorithm for tomographic image reconstruction from limited-data. Math. Comput. Simul. 65, 6 (May 2004), 579-598.</a>
*
* @param {Matrix_} A an m by n matrix, with m >= n.
* @param {Matrix_} b an m by 1 matrix (e.g., a vector).
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.epsAbs the tolerance parameter for the convergence of the algorithm on the absolute error between two iterates, a strictly positive real number; defaults to 1e-14.
* @param {number} opt.epsRel the tolerance parameter for the convergence of the algorithm on the relative error between two iterates, a strictly positive real number; defaults to 1e-12.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 10000.
* @return {<Matrix_} an n by 1 matrix x^* (e.g., a vector) satisfying the following properties:
* - If the linear system of equations Ax = b is square with A invertible, then x^* is the unique solution of this system
* - If the linear system of equations Ax = b is overdetermined and consistent (i.e., b belongs to Range(A)), then x^* is the minimum euclidian norm solution of the least square problem min ||Ax - b||_2
* - If the linear system of equations Ax = b is overdetermined and inconsistent (i.e., b does not belong to Range(A)), then x^* is the minimum euclidian norm solution of the least square problem min ||Ax - b||_2
*
* @example
* linsolveKaczmarz(Matrix_([[1],[2],[3]]), Matrix_([[1],[2],[3]]));
* // XX
*/

// Kaczmarz Extended Algorithm for Tomographic Image Reconstruction from Limited-Data, algorithm KE, formulas 58-60 
Matrix_.linsolveKaczmarz = function(A, b, opt) {
	// Constants
	machineEps = Number.EPSILON || Math.pow(2, -52);

	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var epsAbs = opt.epsAbs || 1e-14;
	var epsRel = opt.epsRel || 1e-12;
	var maxIterations = opt.maxIter || 10000;
	
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
	if (A.nbRows < A.nbColumns) {
		throw new Error('matrix has more columns than rows: ' + '(' + A.nbRows + ') v.s. ' + '(' + A.nbColumns + ')');
	}	

	// ------
	
	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	var x_km = new Matrix_.zeros(n, 1); // the previous solution vector
	var x_k = new Matrix_.zeros(n, 1); // the current solution vector
	var y_k = Matrix_.copy(b); // the current projection of b on the hyperplanes generated by the columns of A	
	var b_k = Matrix_.copy(b); // the current difference betwen b abd y_k
	
	// ------
	
	// Preliminary computation of the squares of the 2-norms of the rows of A
	var a_rows_two_norm_sq = new Array(m);
	for (var i = 1; i <= m; ++i) {
		var a_i_two_norm = A.vectorNorm('two', 'row', i);
		if (a_i_two_norm <= machineEps) {
			a_rows_two_norm_sq[i-1] = 0; // for all practical purposes, this row of A is null
		}
		else {
			a_rows_two_norm_sq[i-1] = a_i_two_norm * a_i_two_norm;
		}
	}
	
	// Preliminary computation of the squares of the 2-norms of the columns of A
	var a_columns_two_norm_sq = new Array(n);
	for (var j = 1; j <= n; ++j) {
		var alpha_j_two_norm = A.vectorNorm('two', 'column', j);
		if (alpha_j_two_norm <= machineEps) {
			a_columns_two_norm_sq[j-1] = 0; // for all practical purposes, this column of A is null
		}
		else {
			a_columns_two_norm_sq[j-1] = alpha_j_two_norm * alpha_j_two_norm;
		}
	}
	
	// ------
	// Main loop until convergence, guaranteed as per formula 64 of the reference.	
	var iter = 0;
	var delta_x = new Matrix_.zeros(n, 1); // the delta vector (current iteration solution vector minus previous iteration solution vector)
	while (true) {
		// Update the number of iterations
		++iter;

		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		
		// Cycle through the n columns of A and orthogonally project the current iterate y_k onto
		// the hyperplanes generated by the columns A(:,j), j=1..n, c.f. formula 58 of the reference.
		for (var j = 1; j <= n; ++j) {
			// Limit case: skip the projection in case the current column of A is null
			if (a_columns_two_norm_sq[j-1] == 0) {
			    continue;
			}
			
			// Compute <A(:,j)/y_k>
			var a_j_y_k = 0;
			for (var i = 1; i <= m; ++i) {
				a_j_y_k += A.data[(i-1) * A.nbColumns + (j-1)] * y_k.data[(i-1) * y_k.nbColumns];
			}
			
			// Update y_k: y_k = y_k - <A(:,j)/y_k>/||A(:,j)||_2^2 * A(:,j)
			for (var i = 1; i <= m; ++i) {
				y_k.data[(i-1) * y_k.nbColumns] -= a_j_y_k / a_columns_two_norm_sq[j-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
		}
		
		// Update b, c.f. formula 59 of the reference.
		b_k = Matrix_.axpby(1, b, -1, y_k, b_k);
		
		// Cycle through the m rows of A and orthogonally project the current iterate x_k onto
		// the solution hyperplane of <A(i,:)/x_k> = b_i, i=1..m, c.f. formula 60 of the reference.
		for (var i = 1; i <= m; ++i) {
			// Limit case: skip the projection in case the current row of A is null
			if (a_rows_two_norm_sq[i-1] == 0) {
			    continue;
			}
			
			// Compute r_k = <A(i,:)/x_k> - b_k(i)
			var a_i_x_k = 0;
			for (var j = 1; j <= n; ++j) {
				a_i_x_k += A.data[(i-1) * A.nbColumns + (j-1)] * x_k.data[(j-1) * x_k.nbColumns];
			}
			var r_k = a_i_x_k - b_k.data[(i-1) * b.nbColumns];
			
			// Update x_k: x_k = x_k - r_k/||A(i,:)||_2^2 * A(i,:)
			for (var j = 1; j <= n; ++j) {
				x_k.data[(j-1) * x_k.nbColumns] -= r_k / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
		}
		
		// Convergence conditions (not in the references, but usual condition for convergent sequences): 
		// - Absolute error: ||x_k - x_km||_inf <= epsAbs
		// - Relative error: ||x_k - x_km||_inf <= epsRel * ( 1 + ||x_k||_inf)
		var delta_x_inf_norm = Matrix_.axpby(1, x_k, -1, x_km, delta_x).vectorNorm('infinity');
		if (delta_x_inf_norm <= epsAbs && delta_x_inf_norm <= epsRel * (1 + x_k.vectorNorm('infinity'))) {
			break;
		}
		
		// Prepare the next iteration
		// Update the previous iteration solution vector to the current iteration solution vector: x_km = x_k
		x_km = Matrix_.copy(x_k, x_km);
	}
	
	// ------
	
	// Return the computed solution
	return x_k;
}

/**
* @function linsolveRandomizedExtendedKaczmarz
*
* @summary TODO
*
* @description TODO, with convergence depending on the square condition number of the matrix A.
* 
* @see <a href="https://arxiv.org/abs/1205.5770">Zouzias, Anastasios and Freris, Nikolaos, Randomized Extended Kaczmarz for Solving Least-Squares, 05/2012, eprint arXiv:1205.5770</a>
*
* @param {Matrix_} A a m by n matrix.
* @param {Matrix_} b a m by 1 matrix.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-12.
aram {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer or -1 to force an infinite number of iterations; defaults to 100000.
* @return {<Matrix_} an n by 1 matrix x^* satisfying the following properties:
* - If the linear system of equations Ax = b is square with A invertible, then x^* is the unique solution of this system
* - If the linear system of equations Ax = b is not square (i.e., overdetermined or underdetermined), then x^* is the minimum euclidian norm solution of the least square problem min ||Ax - b||_2
*
* @example
* linsolveRandomizedExtendedKaczmarz(X);
* // XX
*/
Matrix_.linsolveRandomizedExtendedKaczmarz = function(A, b, opt) {
	// ------
	
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-12;
	var maxIterations = opt.maxIter || 100000;
	
	
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
	var x_k = new Matrix_.zeros(n, 1); // the current solution
	var z_k = Matrix_.copy(b); // the current "adjusted" b
	
	// Preliminary computation of the Frobenius norm of A
	var a_frob_norm = A.matrixNorm('frobenius');
	var a_frob_norm_sq = a_frob_norm * a_frob_norm;
	
	// Limit case: null matrix
	if (a_frob_norm === 0) {
		return x_k;
	}
	
	// Preliminary computation of the squares of the 2-norms of the rows of A,
	// as well as the probabilities q_i with their associated sampler.
	var a_rows_two_norm_sq = new Array(m);
	var q = new Array(m);
	for (var i = 1; i <= m; ++i) {
		var a_i_two_norm = A.vectorNorm('two', 'row', i);
		a_rows_two_norm_sq[i-1] = a_i_two_norm * a_i_two_norm;
		
		q[i-1] = a_rows_two_norm_sq[i-1]/a_frob_norm_sq;
	}
	var qSampler = new aliasMethodSampler_(q);
	
	// Preliminary computation of the squares of the 2-norms of the columns of A,
	// as well as the probabilities p_j with their associated sampler.
	var a_columns_two_norm_sq = new Array(n);
	var p = new Array(n);
	for (var j = 1; j <= n; ++j) {
		var alpha_j_two_norm = A.vectorNorm('two', 'column', j);
		a_columns_two_norm_sq[j-1] = alpha_j_two_norm * alpha_j_two_norm;
		
		p[j-1] = a_columns_two_norm_sq[j-1]/a_frob_norm_sq;
	}
	var pSampler = new aliasMethodSampler_(p);
	
	
	// ------
	
	// Main loop until convergence, guaranteed as per theorem 8 of the reference.	
	var iter = 0;
	var x_res = new Matrix_.zeros(m, 1); // the x residuals vector
	var b_res = new Matrix_.zeros(m, 1); // the b residuals vector
	var a_x_k = new Matrix_.zeros(m, 1); // A*x_k
	var ta_z_k = new Matrix_.zeros(n, 1); // A^t*z_k
	while (true) {
		// Update the number of iterations
		++iter;

		// Check the number of iterations
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		
		
		// Pick a column index j with probability p_j
		var j = pSampler.sample() + 1;
		
		// Orthogonally project the current iterate y_k onto the hyperplane generated 
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
			a_x_k =  Matrix_.axy(1, A, x_k, a_x_k);
			b_res = Matrix_.axpby(1, b, -1, z_k, b_res);
			x_res = Matrix_.axpby(1, a_x_k, -1, b_res, x_res);
			var x_res_two_norm = x_res.vectorNorm('two');
			
			// Compute ||A^tz_k||_2
			ta_z_k = Matrix_.atxy(1, A, z_k, ta_z_k);
			var ta_z_k_two_norm = ta_z_k.vectorNorm('two');
			
			if (x_res_two_norm <= eps * a_frob_norm * x_k_two_norm && 
				ta_z_k_two_norm <= eps * a_frob_norm_sq * x_k_two_norm) {
				break;
			}
		}
	}
	
	
	// ------
	
	// Return the computed solution
	return x_k;
}
