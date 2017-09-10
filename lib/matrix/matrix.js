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
		that.data = FLOAT64_ARRAY(that.nbRows * that.nbColumns);
		
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
		that.data = FLOAT64_ARRAY(that.nbRows * that.nbColumns);
		
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
		that.data = FLOAT64_ARRAY(that.nbRows * that.nbColumns);
		
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
    	
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
    	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
    	
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
		obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
		
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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

