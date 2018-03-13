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
  var vect = new PortfolioAllocation.Matrix(this.vectValues);
	
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
  var vec = new PortfolioAllocation.Matrix([1,2,3]);
  assert.deepEqual(vec.toArray(), [1,2,3], 'Vector to array');
  
  // Test vector to array using the random vector, using identity vector(vector to array) == vector
  var nsVec = new PortfolioAllocation.Matrix(this.vectValues);
  var nsVecArray = nsVec.toArray();
  assert.deepEqual(nsVecArray, this.vectValues, 'Array to vector to array');
  assert.equal(PortfolioAllocation.Matrix.areEqual(nsVec, new PortfolioAllocation.Matrix(nsVecArray)), true, 'Vector to array to vector');
});


QUnit.test('Vector dot product', function(assert) {    
  // Test vector dot product using static data
  var vec1 = new PortfolioAllocation.Matrix([1,2,3]);
  var vec2 = new PortfolioAllocation.Matrix([4,5,6]);
  assert.equal(PortfolioAllocation.Matrix.vectorDotProduct(vec1, vec2), 32, 'Vector dot product');
  
  // Test vector dot product using formula n(n+1)(2n+1)/6 == dot([1,2,3,...,n],[1,2,3...,n])
  var n_max = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var intArray = new Array(n_max);
  for (var n = 1; n <= n_max; ++n) {
    intArray[n-1] = n;
  }
  var vec = new PortfolioAllocation.Matrix(intArray);
  assert.equal(PortfolioAllocation.Matrix.vectorDotProduct(vec, vec), n_max*(n_max+1)*(2*n_max+1)/6, 'Vector dot product integer');
});


QUnit.test('Vector Hadamard product', function(assert) {    
  // Test vector Hadamard product using static data
  var vec1 = new PortfolioAllocation.Matrix([1,2,3]);
  var vec2 = new PortfolioAllocation.Matrix([1,2,3]);
  
  var hadamardVec = PortfolioAllocation.Matrix.vectorHadamardProduct(vec1, vec2);
  var expectedVec = new PortfolioAllocation.Matrix([1,4,9]);
  
  assert.equal(PortfolioAllocation.Matrix.areEqual(hadamardVec, expectedVec), true, 'Vector Hadamard product');
});



QUnit.test('Matrix-vector product', function(assert) {    
  // Test matrix-vector product using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var vec = new PortfolioAllocation.Matrix([1,2,3]);
  var prodVec = PortfolioAllocation.Matrix.xy(mat, vec);

  var expectedResVec = new PortfolioAllocation.Matrix([14,32]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(prodVec, expectedResVec), true, 'Matrix-vector product');
});


QUnit.test('Vector sum', function(assert) {    
  // Test vector sum using formula n == sum([1,1,1,...,1])
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.Matrix.ones(n, 1); 
  assert.equal(vec.sum(), n, 'Vector sum');
});


QUnit.test('Vector normalisation', function(assert) {    
  // Static data
  var vec = new PortfolioAllocation.Matrix([1,2,3]);
  var resvec = new PortfolioAllocation.Matrix([1/6,2/6,3/6]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(vec.normalize(), resvec), true, 'Vector static normalized');

  // Vector of ones
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.Matrix.ones(n, 1).normalize(); 
  for (var i = 1; i <= n; ++i) {
    assert.equal(vec.getValueAt(i,1), 1/n, 'Vector of ones normalized');
  }
});


QUnit.test('Specific vectors creation', function(assert) {    
  // Vector of ones
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.Matrix.ones(n, 1); 
  
  for (var i = 1; i <= n; ++i) {
    assert.equal(vec.getValueAt(i,1), 1, 'Vector of ones');
  }
  
  // Vector of zeros
  var n = Math.floor(Math.random()*(100-2+1) + 2); // min = 2, max = 100
  var vec = PortfolioAllocation.Matrix.zeros(n, 1); 
  
  for (var i = 1; i <= n; ++i) {
    assert.equal(vec.getValueAt(i,1), 0, 'Vector of zeros');
  }

});


