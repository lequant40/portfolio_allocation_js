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
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
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
		
		// Return the computed sum
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);

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
	* @return {Matrix_} a column matrix (i.e., vector) containing the diagonal elements of the original matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).diagonal();
	* // Matrix_([1,5])
	*/
    diagonal: function () {
    	// Checks
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
    	
    	// Result vector instantiation
    	var obj = Object.create(Matrix_.prototype);
    	obj.nbRows = this.nbRows;
    	obj.nbColumns = 1;
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
    	
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
	* @return {Matrix_} a column matrix (i.e., vector) containing the elements of the i-th row of the original matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).row(1);
	* // Matrix_([1,2])
	*/
    row: function (i) {
		// Bounds check
		if (i < 1 || i > this.nbRows) {
			throw new Error(
			'index out of bounds when getting matrix row, (' + i + 
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
    	
    	// Result vector instantiation
    	var obj = Object.create(Matrix_.prototype);
    	obj.nbRows = this.nbColumns;
    	obj.nbColumns = 1;
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
    	
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
	* @param {Array.<number>} rindexes the row indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..number of matrix rows
    * @param {Array.<number>} cindexes the column indexes of the original matrix elements to keep, array of strictly increasing natural integers belonging to 1..number of matrix columns
	* @return {Matrix_} a matrix whose elements correspond to the elements of the original matrix whose row/column indexes belong to the input lists of row/column indexes to keep.
	*
	* @example
	* Matrix_([[1,2,3], [4,5,6]]).submatrix([1], [2, 3];
	* // Matrix_([[2,3]])
	*/
    submatrix : function(rindexes, cindexes) {
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
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
    	
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
		
		// Compute the (Givens) QR decomposition of the matrix
		var qr = Matrix_.qrDecomposition(this);
		
		// By property of the Givens QR decomposition, the determinant of the matrix
		// is then the product of the diagonal elements of the R matrix.
		var r = qr[1];
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
	*
	* If the parameter subset is provided and is equal to 'row', the parameter idx correspond to the row index of the matrix
	* for which to compute the vector norm.
	*
	* @see <a href="https://en.wikipedia.org/wiki/Norm_(mathematics)">Norm (mathematics)</a>
	* @see Nicholas J. Higham. 2002. Accuracy and Stability of Numerical Algorithms (2nd ed.). Soc. for Industrial and Applied Math., Philadelphia, PA, USA. 
	*
	* @memberof Matrix_
	* @param {string} p the vector norm to compute, a string either equals to 'one', to 'two' or to 'infinity'.
	* @param {string} subset the subset of the matrix for which to compute the vector norm, an optional string either equals to 'matrix' or to 'row'; defaults to 'matrix'.
	* @param {number} idx if subset is equal to 'row', the row index of the matrix for which to compute the vector norm, a natural integer belonging to 1..nbRows.
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
		
		// Supported subsets are 'matrix' and 'row'
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
			if (Math.abs(a.data[i * nbRows + j] - b.data[i * nbRows + j]) > tol) {
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
* The algorithm implemented uses an IKJ form, cache aware.
* 
* @param {Matrix_} A a n by m matrix.
* @param {Matrix_} B a m by p matrix.
* @param {Matrix_} out an optional n by m matrix.
* @return {Matrix_} the matrix product A*B, either stored in the matrix out or in a new matrix, a m by p matrix.
*
* @example
* product(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,1], [2,2], [3,3]]));
* // Matrix_([[14,14], [32,32]])
*/
Matrix_.product = function(a, b, out) {
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
	
	// Result matrix allocation
	var obj = allocateMatrix_(a.nbRows, b.nbColumns, out);

	// Computation of A*B product in IKJ format, cache aware
	for (var i = 0; i < a.nbRows; ++i) {
		for (var j = 0; j < b.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}
	
		for (var k = 0; k < a.nbColumns; ++k) {
			var a_ik = a.data[i * a.nbColumns + k];
			
			for (var j = 0; j < b.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] += a_ik * b.data[k * b.nbColumns + j];
			}
		}
	}
	
	// Return the computed matrix
    return obj;
};

/**
* @function rowColumnProduct
*
* @summary Returns the product of a matrix row and a matrix column.
*
* @description This function computes the product A(i,:) * B(:,j) of the i-th row
* of the matrix A with the j-th column of the matrix B, where A is a n by m matrix
* and B is a m by p matrix.
*
* @param {Matrix_} A a n by m matrix.
* @param {number} i the row index of the matrix A, a natural integer belonging to 1..n.
* @param {Matrix_} B a m by p matrix.
* @param {number} j the column index of the matrix B, a natural integer belonging to 1..m.
* @return {number} the product A(i,:) * B(:,j), a real number.
*
* @example
* rowColumnProduct(Matrix_([[1,2,3], [4,5,6]]), 1, Matrix_([[1,1], [2,2], [3,3]]), 1);
* // 14
*/
Matrix_.rowColumnProduct = function(a, i, b, j) {
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
	
	// Bounds check for a
	if (i < 1 || i > a.nbRows) {
		throw new Error(
		'index out of bounds when getting matrix row, (' + i + 
		') in size (' + a.nbRows + ',' + a.nbColumns + ')');
	}
	
	// Bounds check for b
	if (j < 1 || j > b.nbRows) {
		throw new Error(
		'index out of bounds when getting matrix column, (' + j + 
		') in size (' + b.nbRows + ',' + b.nbColumns + ')');
	}
	
	// Computation of <A(i,:))/B(:,j)>
	var rowColProd = 0;
	for (var k = 0; k < a.nbColumns; ++k) {
		rowColProd += a.data[(i-1) * a.nbColumns + k] * b.data[k * b.nbColumns + (j-1)]; 
	}
	
	// Return the computed value
	return rowColProd;
};


/**
* @function diagonal
*
* @summary Returns a diagonal matrix constructed from an array of numbers.
*
* @description This function computes a diagonal square matrix (a_ij),i=1..n,j=1..n from an array arr of n real numbers, 
* with coefficients a_ij satisfying a_ij = 0, i <> j and a_ij = arr[i-1], i==j.
*
* @param {Array<number>} x an array of real numbers of size n.
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

	// Result matrix computation
	return Matrix_.fillSymetric(arr.length, function(i,j) { if (i == j) { return arr[i-1]; } else { return 0; } });
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
* @return {Matrix_} a n by n symetric matrix with its elements computed by the fonction fct.
*
* @example
* fillSymetric(2, function(i,j) { return 0; });
* // == Matrix_([[0,0], [0,0]])
*/
Matrix_.fillSymetric = function(n, fct) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows and columns out of bounds: ' + n);
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(n, n);
	
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
* @return {Matrix_} a n by m matrix with its elements computed by the fonction fct.
*
* @example
* fill(2, 2, function(i,j) { return 0; });
* // == Matrix_([[0,0], [0,0]])
*/
Matrix_.fill = function(n, m, fct) {
	// Checks
	if (n < 1) {
		throw new Error('input number of rows out of bounds: ' + n);
	}
	if (m < 1) {
		throw new Error('input number of columns out of bounds: ' + m);
	}
	
	// Result matrix allocation
	var obj = allocateMatrix_(n, m);
	
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
* @return {Matrix_} the constructed matrix.
*
* @example
* zeros(3, 2);
* // Matrix_([[0,0], [0,0], [0,0]])
*/
Matrix_.zeros = function(n, m) {
	return Matrix_.fill(n, m, function(i,j) { return 0; });	

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
	return Matrix_.fill(n, m, function(i,j) { return 1; });	
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
	return Matrix_.fillSymetric(n, function(i,j) { if (i == j) { return 1; } else { return 0; } });
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
	
	// Result matrix allocation
	var obj = allocateMatrix_(x.nbRows, 1);
	
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
* @return {<Array.<Matrix_>} an array of two matrices - Q at array index 0 and R at array index 1 - satisfying the following properties:
* - Q is an m by m matrix 
* - R is an m by n matrix
* - A = QR
* - Q is an orthogonal matrix, with a determinant equals to 1 (i.e., a rotation)
* - R is an upper triangular matrix
*
* @example
* qrDecomposition(Matrix_([[1],[2],[3]]), Matrix_([[1],[2],[3]]));
* // XX
*/
Matrix_.qrDecomposition = function(a) {
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
		var c = NaN;
		var s = NaN;
		var r = NaN;
		
		if (b == 0) {
			c = (a >= 0) ? 1 : -1; // Emulates sign(a)
			s = 0;
			r = Math.abs(a);
		}
		else if (a == 0) {
			c = 0;
			s = ((b >= 0) ? 1 : -1); // Emulates sign(b)
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


	// Initializations
	var m = a.nbRows;
	var n = a.nbColumns;

	// Checks
	if (m < n) {
		throw new Error('matrix has more columns than rows: ' + '(' + m + ') v.s. ' + '(' + n + ')');
	}
	
	// Create a copy of A so that it is not overwritten
	var rr = new Matrix_(a); // Represents R
	
	// Create the matrix that will hold Q
	var qq = Matrix_.identity(m); // Represents Q
	
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
			for (var k = 1; k <= m; ++k) {
				var t1 = qq.data[(k-1) * qq.nbColumns + (i-2)] // t1 = Q(k,i-1)
				var t2 = qq.data[(k-1) * qq.nbColumns + (i-1)] // t2 = Q(k,i)
				qq.data[(k-1) * qq.nbColumns + (i-2)] = c*t1 - s*t2; // Q(k,i-1) = ...
				qq.data[(k-1) * qq.nbColumns + (i-1)] = s*t1 + c*t2; // Q(k,i) = ...
			}
			
		}
	}
	
	// Return the computed [Q, R] pair;
	return [qq, rr];
}


/**
* @function linsolveKaczmarz
*
* @summary TODO
*
* @description TODO
* 
* @see <a href="https://www.sciencedirect.com/science/article/pii/002437959090207S">Hanke, M. and Niethammer, W. [1990], On the acceleration of Kaczmarz’smethod for inconsistent linear systems, Linear Algebra and its Applications 130, 83–98.</a>
*
* @param {Matrix_} A an m by n matrix, with m >= n.
* @param {Matrix_} b an m by 1 column matrix (e.g., a vector).
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-16.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 10000.
* @return {<Matrix_} an n by 1 column matrix x^* (e.g., a vector) satisfying the following properties:
* - If the linear system of equations Ax = b is square with A invertible, then x^* is the unique solution of this system
* - If the linear system of equations Ax = b is overdetermined and consistent (i.e., b belongs to Range(A)), then x^* is the minimum euclidian norm solution of the least square problem min ||Ax - b||
* - If the linear system of equations Ax = b is overdetermined and inconsistent (i.e., b does not belong to Range(A)), then x^* is a solution of the least square problem min ||Ax - b||
*
* @example
* linsolveKaczmarz(Matrix_([[1],[2],[3]]), Matrix_([[1],[2],[3]]));
* // XX
*/
Matrix_.linsolveKaczmarz = function(A, b, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-16;
	var maxIterations = opt.maxIter || 10000;
	
	// ------
	
	// Ensure a,b are matrices
	if (!(A instanceof Matrix_)) {
		throw new Error('a must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('b must be a matrix');
	}
	
	// Checks
	if (A.nbColumns !== b.nbRows) {
		throw new Error('Matrices sizes do not match: ' + '(' + A.nbRows + ',' + A.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns !== 1) {
		throw new Error('b is not a vector: ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	
	// Initializations
	var m = A.nbRows;
	var n = A.nbColumns;
	var x_km = new Matrix_.zeros(n, 1); // the previous iteration solution vector; x in the algorithm 2.1 of the reference
	var x_k = new Matrix_.zeros(n, 1); // the current iteration solution vector; u in the algorithm 2.1 of the reference
	
	// Checks
	if (m < n) {
		throw new Error('matrix has more columns than rows: ' + '(' + m + ') v.s. ' + '(' + n + ')');
	}
	
	// Preliminary computation of the squares of the 2-norms of the rows of A
	var a_rows_two_norm_sq = new Array(m);
	for (var i = 1; i <= m; ++i) {
		var a_i_two_norm = A.vectorNorm('two', 'row', i);
		a_rows_two_norm_sq[i-1] = a_i_two_norm * a_i_two_norm;
	}
	
	// Main loop until convergence, guaranteed as per theorem 3.1 of the reference
	var iter = 0;
	var converged = false;
	var delta_x = new Matrix_.zeros(n, 1); // the delta vector (current iteration solution vector minus previous iteration solution vector)
	while (!converged) {
		// Cycle through the m rows of A and orthogonally project the current iterate x_k onto
		// the solution hyperplane of <A(i,:)/x_k> = b_i, i=1..m.
		for (var i = 1; i <= m; ++i) {
			// Limit case: skip the projection in case the current row of A is null
			if (a_rows_two_norm_sq[i-1] == 0) {
			    continue;
			}
			
			// Compute r, c.f. algorithm 2.1 of the reference: r = b(i) - <A(i,:)/x_k>
			var r = b.data[(i-1) * b.nbColumns] - Matrix_.rowColumnProduct(A, i, x_k, 1);
			
			// Update x_k, c.f. algorithm 2.1 of the reference: x_k = x_k + r/||A(i,:)||_2^2 * A(i,:)
			for (var j = 1; j <= n; ++j) {
				x_k.data[(j-1) * x_k.nbColumns] += r / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
		}
		
		// Check the convergence conditions: 
		// - Absoute error: ||x_k - x_km||_abs <= eps 
		// - Relative error ||x_k - x_km||_abs / ||x_k||_abs <= eps
		delta_x = Matrix_.axpby(1, x_k, -1, x_km, delta_x);
		
		var delta_x_inf_norm = delta_x.vectorNorm('infinity');
		var x_k_inf_norm = x_k.vectorNorm('infinity');
		if (delta_x_inf_norm <= eps && delta_x_inf_norm <= eps * x_k_inf_norm) {
			converged = true;
		}
		
		// Prepare the next iteration
		// Update the previous iteration solution vector to the current iteration solution vector: x_km = x_k
		x_km = Matrix_.copy(x_k, x_km);
		
		// Update the number of iterations
		++iter;
		
		// Check the number of iterations
		if (iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
	}
	
	// Return the computed solution
	return x_k;
}
