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
	
	/**
	* @function getDiagonal
	*
	* @summary Returns the diagonal elements of a square matrix.
	*
	* @description This function returns, as a vector of n rows, the diagonal elements (a_ii), i=1..n from the original square matrix (a_ij),i=1..n,j=1..n.
	*
	* @memberof Matrix_
	* @return {Vector_} a vector containing the diagonal coefficients of the original matrix.
	*
	* @example
	* Matrix_([[1,2], [4,5]]).getDiagonal();
	* // Vector_([1,5])
	*/
    getDiagonal: function () {
    	// Checks
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
    	
    	// Result vector instantiation
    	var obj = Object.create(Vector_.prototype);
    	obj.nbRows = this.nbRows;
    	obj.nbColumns = 1;
    	obj.data = new Array(obj.nbRows);
    	
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
		var obj = Object.create(this);
		obj.nbRows = this.nbRows;
		obj.nbColumns = this.nbColumns;
		obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
		
		// Computation of the fonction fct applied to the coefficients of A
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = fct(i+1, j+1, this.data[i * obj.nbColumns + j]);
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
* @return {Vector_} the matrix product a*b, an m by p matrix.
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
* @function areEqual
*
* @summary Determines if two vectors are identical.
*
* @description This function determines if two vectors x and y are identical, i.e. if x_i == y_i, i=1..n where
* n is their common number of rows.
*
* In case an optional tolerance parameter eps is provided, the strict component by component equality condition x_i == y_i
* is relaxed into |x_i - y_i| <= eps.
* 
* @param {Vector_} x, a vector.
* @param {Vector_} y, a vector.
* @param {number} eps, an optional real number; by default, eps is equal to zero.
* @return {boolean} true if x and y are identical, false otherwise.
*
* @example
* areEqual(Matrix_([[1,2,3], [4,5,6]]), Matrix_([[1,2,3], [4,5,6]]));
* // true
*/
Vector_.areEqual = function (x, y, eps) {
	// Ensure x,y are vectors
	if (!(x instanceof Vector_)) {
		throw new Error('x must be a vector');
	}
	if (!(y instanceof Vector_)) {
		throw new Error('y must be a vector');
	}
	
    // Real computation is delegated to the equivalent Matrix method
    return Matrix_.areEqual(x, y, eps);
}


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
* @function vectorProduct
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
* vectorProduct(Matrix_([[1,2,3], [4,5,6]]), Vector_([1,2,3]));
* // Vector_([14,32])
*/
Matrix_.vectorProduct = function(a, x) {
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
* @function hadamardProduct
*
* @summary Returns the Hadamard product of two vectors.
*
* @description This function computes the Hadamard product x*y of two vectors x and y of the same size.
* 
* @param {Vector_} x a vector.
* @param {Vector_} y a vector of same size as x.
* @return {Vector_} the Hadamard product x*y.
*
* @example
* hadamardProduct(Vector_([1,2,3]), Vector_([1,2,3]));
* // Vector_([1,4,9])
*/
Vector_.hadamardProduct = function(x, y) {
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
	
	// Result vector instantiation
	var obj = Object.create(Vector_.prototype);
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
* @function dotProduct
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
* dotProduct(Vector_([1,2,3]), Vector_([1,2,3]));
* // 14
*/
Vector_.dotProduct = function(x, y) {
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
* @function ones
*
* @summary Returns a vector made of ones.
*
* @description This function builds a vector (v_i),i=1..n of size n and satisfying v_i = 1, i=1..n.
* 
* @param {number} n the row length of the vector to construct.
* @return {Vector_} the constructed vector.
*
* @example
* ones(3);
* // Vector_([1,1,1])
*/
Vector_.ones = function(n) {
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
	
	// ------
	
	// The output weights are defined as the normalized inverses of the assets standard deviations.
	var weights = new Vector_(sigma).elemPower(-1/2).normalize();
	
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
* This portfolio is a special case of the more generic risk budgeting portfolio, with all risk budgets
* equal.
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
	var rb = Vector_.ones(nbAssets).normalize();

	// Compute the associated risk budgeting weights
	return self.riskBudgetingWeights(sigma, rb.toArray(), opt);
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
	var weights = Vector_.ones(nbAssets).normalize();

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
* This portfolio maximizes the Sharpe ratio if the Sharpe ratio for each stock is the same.
*
* The algorithm used is a cyclical coordinate descent, c.f. the second reference, whose convergence is guaranteed
* if the covariance matrix of the assets is definite positive.
*
* @see <a href="http://www.iijournals.com/doi/abs/10.3905/JPM.2008.35.1.40">Toward Maximum Diversification by Y. Choueifaty, Y. Coignard, The Journal of Portfolio Management, Fall 2008, Vol. 35, No. 1: pp. 40-51</a>
* @see <a href="https://ssrn.com/abstract=2595051">Richard, Jean-Charles and Roncalli, Thierry, Smart Beta: Managing Diversification of Minimum Variance Portfolios (March 2015)</a>
* 
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm.
* @param {number} opt.eps the tolerance parameter for the convergence of the algorithm, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the algorithm, a strictly positive natural integer; defaults to 1000.
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
	var maxIterations = opt.maxIter || 1000;
	
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
	var x = Vector_.ones(nbAssets).normalize();
	
	// Preparational computations
	var sigma_x = Matrix_.vectorProduct(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Vector_.dotProduct(sigma_x, x)); // sigma(x)
	var dr = Vector_.dotProduct(stddevs, x) / s_x; // diversification ratio DR
	
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
    	    sigma_x = Matrix_.vectorProduct(sigma, x)
    	    s_x = Math.sqrt(Vector_.dotProduct(sigma_x, x));
			
			  // Case #1 diversification ratio computation
			var dr_star = Vector_.dotProduct(stddevs, x) / s_x; // Weighted average volatility of the portfolio divided by the volatility of the portfolio
			
			// Case #2 test
			if (r1 > 0) {
				var xi_star_2 = r1;

				// Case #2 weights update
				x.setValueAt(i, 1, xi_star_2);
				
				// Case #2 updated SIGMA*x and x'*SIGMA*x products
				var sigma_x_2 = Matrix_.vectorProduct(sigma, x)
				var s_x_2 = Math.sqrt(Vector_.dotProduct(sigma_x_2, x));
				
				// Case #2 diversification ratio computation
				var dr_star_2 = Vector_.dotProduct(stddevs, x) / s_x_2;
				
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
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithm, unused.
* @return {Array.<number>} the weights corresponding to theminimum correlation (heuristic) portfolio, array of n real numbers.
*
* @example
* minCorrWeights([[0.1, 0], [0, 0.2]]);
* // ~[0.59, 0.41]
*/
self.minCorrWeights = function (sigma, opt) {
	// Convert sigma to matrix format
	var sigma = new Matrix_(sigma);
	
	// Extract the correlation matrix and the misc. variances-related vectors
	// Note: this could be optimized to O(n) instead of O(n^2), if ever necessary
	var variances = sigma.getDiagonal();
	var invStddevs = variances.elemPower(-1/2);
	var matInvStddevs = Matrix_.diagonal(invStddevs.toArray());
	var rho = Matrix_.product(matInvStddevs, Matrix_.product(sigma, matInvStddevs));

	// TODO: Checks, if enabled	
	
	// ------
	var nbAssets = rho.nbRows;
	
	// Specific case to be filtered out: number of assets is 2 (or 1, but this case is less useful)
	if (nbAssets <= 2) {
		return self.equalRiskBudgetWeights(variances.toArray(), opt);
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
	var weights = new Vector_(ranks).normalize();
	
	// Step 6: Combine rank portfolio weight estimates with Adjusted Correlation Matrix
	weights = Matrix_.vectorProduct(adjustedRho, weights).normalize();
	
	// Step 7: Scale portfolio weights by assets standard deviations
	weights = Vector_.hadamardProduct(weights, invStddevs).normalize();
	
	// Return the computed weights
	return weights.toArray();
}


/**
 * @file Functions related to proportional minimum variance(heuristic) portfolio.
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
* @param {Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, array of n array of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
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
	
	// Specific case to be filtered out: number of assets is 2 (or 1, but this case is less useful)
	// TODO ?

	// Step 1: Average pairwise covariance, and associated mean/standard deviation
	var rowsElements = sigma.toDoubleArray();
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	var elementsMean = mean_(rowsAverages);
	var elementsStddev = sampleStddev_(rowsAverages);

	// Step 2: Gaussian convertion, and proportional average covar weigth
	var weights = new Vector_(rowsAverages).elemMap(function(i, j, val) { 
		return 1 - normcdf_((val - elementsMean)/elementsStddev);
	}).normalize();
	
	// Step 3: Scale portfolio weights by assets variances
	var invVariancesWeights = sigma.getDiagonal().elemPower(-1).normalize();
	weights = Vector_.hadamardProduct(weights, invVariancesWeights).normalize();
	
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
	var x = Vector_.ones(nbAssets).normalize();
	
	// Preparational computations
	var sigma_x = Matrix_.vectorProduct(sigma, x); // SIGMA*x
	var s_x = Math.sqrt(Vector_.dotProduct(sigma_x, x)); // sigma(x)
	
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
			    throw new Error('discriminant not positive during iteration ' + iter + ', covariance matrix might not be psd');
			}
    	    var q = -(b_p + sign_b_p * Math.sqrt(disc));
    	    var r1 = q/a;
    	    var r2 = c/q;
    	    var xi_star = r1 > 0 ? r1 : r2;
    	    
    	    // Update the weights
    	    x.setValueAt(i,1,xi_star);
    	    
    	    // Compute the updated SIGMA*x and x'*SIGMA*x products for convergence condition evaluation + next loop evaluation
    	    sigma_x = Matrix_.vectorProduct(sigma, x)
    	    s_x = Math.sqrt(Vector_.dotProduct(sigma_x, x));
    	    
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
