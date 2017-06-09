/**
 * @file Functions related to vector object.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.Vector_ = function(arr) { return new Vector_(arr); }
self.vectorDotProduct_ = function(x, y) { return vectorDotProduct_(x, y); }
self.matrixVectorProduct_ = function(a, x) { return matrixVectorProduct_(a, x); }
self.vectorOnes_ = function(n) { return vectorOnes_(n); }
self.vectorSum_ = function(x) { return vectorSum_(x); }
/* End Wrapper private methods - Unit tests usage only */


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
	this.data = FLOAT64_ARRAY(this.nbRows * this.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);

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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
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
	obj.data = FLOAT64_ARRAY(obj.nbRows * obj.nbColumns);
	
	// Vector filling
	for (var i = 0; i < obj.nbRows; ++i) {
		obj.data[i] = 1;
	}
	
	// Return
	return obj;
}


