// ------------------------------------------------------------
QUnit.module('Vector internal module', {
  before: function() {
    // Generate a random vector
	 // Generate the dimension
	var max = 10;
	var min = 2;
	var rVect = Math.floor(Math.random()*(max-min+1) + min);
	
	// Generate the values 
	var simpleArray = new Array(rVect);
    for (var i = 0; i < rVect; ++i) {
	  simpleArray[i] = Math.random();
    }
	this.vectValues = simpleArray;	
	
    // Generate a random non square matrix
	 // Generate the dimensions
	var max = 10;
	var min = 2;
	var c = rVect;
	var r = Math.floor(Math.random()*(max-min+1) + min);
	while (r == c) {
	  r = Math.floor(Math.random()*(max-min+1) + min);
	}
	
	// Generate the values 
	var doubleArray = new Array(r);
    for (var i = 0; i < r; ++i) {
	  doubleArray[i] = new Array(c);
	  for (var j = 0; j < c; ++j) {
	    doubleArray[i][j] = Math.random();
	  }
    }
	this.nsMatValues = doubleArray;	
  }
});


QUnit.test('Vector basic manipulations', function(assert) {    
  // Test vector creation using random data
  var vect = PortfolioAllocation.Vector_(this.vectValues);

  // Check the dimensions	
  assert.equal(vect.nbRows, this.vectValues.length, 'Number of rows');
  assert.equal(vect.nbColumns, 1, 'Number of columns');
  
  // Ensure the vector is not square
  assert.equal(vect.isSquare(), false, 'Vector not square');
  
  // Check the values with getValueAt
  for (var i = 0; i < vect.nbRows; ++i) {
	  assert.equal(vect.getValueAt(i+1,1), this.vectValues[i], 'Vector values');
  }
});


QUnit.test('Vector to array', function(assert) {    
  // Test vector to array using static data
  var vec = PortfolioAllocation.Vector_([1,2,3]);
  assert.deepEqual(vec.toArray(), [1,2,3], 'Vector to array');
  
  // Test vector to array using the random vector, using identity vector(vector to array) == vector
  var nsVec = PortfolioAllocation.Vector_(this.vectValues);
  var nsVecArray = nsVec.toArray();
  assert.deepEqual(nsVecArray, this.vectValues, 'Array to vector to array');
  assert.equal(PortfolioAllocation.matrixIdentical_(nsVec, PortfolioAllocation.Vector_(nsVecArray)), true, 'Vector to array to vector');
});


QUnit.test('Vector dot product', function(assert) {    
  // Test vector dot product using static data
  var vec1 = PortfolioAllocation.Vector_([1,2,3]);
  var vec2 = PortfolioAllocation.Vector_([4,5,6]);
  assert.equal(PortfolioAllocation.vectorDotProduct_(vec1, vec2), 32, 'Vector dot product');
  
  // Test vector dot product using formula n(n+1)(2n+1)/6 == dot([1,2,3,...,n],[1,2,3...,n])
  var n_max = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var intArray = new Array(n_max);
  for (var n = 1; n <= n_max; ++n) {
    intArray[n-1] = n;
  }
  var vec = PortfolioAllocation.Vector_(intArray);
  assert.equal(PortfolioAllocation.vectorDotProduct_(vec, vec), n_max*(n_max+1)*(2*n_max+1)/6, 'Vector dot product integer');
});


QUnit.test('Matrix-vector product', function(assert) {    
  // Test matrix-vector product using static data
  var mat = PortfolioAllocation.Matrix_([[1,2,3], [4,5,6]]);
  var vec = PortfolioAllocation.Vector_([1,2,3]);
  var prodVec = PortfolioAllocation.matrixVectorProduct_(mat, vec);

  var expectedResVec = PortfolioAllocation.Vector_([14,32]);
  assert.equal(PortfolioAllocation.matrixIdentical_(prodVec, expectedResVec), true, 'Matrix-vector product');
});


QUnit.test('Vector sum', function(assert) {    
  // Test vector sum using formula n == sum([1,1,1,...,1])
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.vectorOnes_(n); 
  assert.equal(vec.sum(), n, 'Vector sum');
});


QUnit.test('Vector normalisation', function(assert) {    
  // Static data
  var vec = PortfolioAllocation.Vector_([1,2,3]);
  var resvec = PortfolioAllocation.Vector_([1/6,2/6,3/6]);
  assert.equal(PortfolioAllocation.matrixIdentical_(vec.normalize(), resvec), true, 'Vector static normalized');

  // Vector of ones
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.vectorOnes_(n).normalize(); 
  for (var i = 1; i <= n; ++i) {
    assert.equal(vec.getValueAt(i,1), 1/n, 'Vector of ones normalized');
  }
});


QUnit.test('Specific vectors creation', function(assert) {    
  // Vector of ones
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.vectorOnes_(n); 
  
  for (var i = 1; i <= n; ++i) {
    assert.equal(vec.getValueAt(i,1), 1, 'Vector of ones');
  }
});


