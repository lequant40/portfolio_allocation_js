// ------------------------------------------------------------
QUnit.module('Matrix internal module', {
  before: function() {
    // Generate a random non square matrix
	 // Generate the dimensions
	var max = 10;
	var min = 2;
	var r = Math.floor(Math.random()*(max-min+1) + min);
	var c = Math.floor(Math.random()*(max-min+1) + min);
	while (c == r) {
	  c = Math.floor(Math.random()*(max-min+1) + min);
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
	
	// Generate a random square matrix
	var doubleArray2 = new Array(r);
    for (var i = 0; i < r; ++i) {
	  doubleArray2[i] = new Array(r);
	  for (var j = 0; j < r; ++j) {
	    doubleArray2[i][j] = Math.random();
	  }
    }
	this.sMatValues = doubleArray2;	
  }
});


QUnit.test('Matrix basic manipulations', function(assert) {    
  // Test non square matrix using random data
  {
      var nsMat = new PortfolioAllocation.Matrix(this.nsMatValues);
    
      // Check the dimensions	
      assert.equal(nsMat.nbRows, this.nsMatValues.length, 'Number of rows');
      assert.equal(nsMat.nbColumns, this.nsMatValues[0].length, 'Number of columns');
      
      // Ensure the matrix is not square
      assert.equal(nsMat.isSquare(), false, 'Matrix not square');
    
      // Ensure the matrix is equal to itself
      assert.equal(PortfolioAllocation.Matrix.areEqual(nsMat,nsMat), true, 'Matrix equal to itself');
      
      // Check the values with getValueAt
      for (var i = 0; i < nsMat.nbRows; ++i) {
    	for (var j = 0; j < nsMat.nbColumns; ++j) {
    	  assert.equal(nsMat.getValueAt(i+1,j+1), this.nsMatValues[i][j], 'Matrix values');
    	}
      }
  }
  
  // Test equality using the non square matrix
  {
      // Create a new matrix equal to the previous one and replace 
      // the values of the matrix with new values with setValueAt
      var nsMat2 = new PortfolioAllocation.Matrix(this.nsMatValues);
      for (var i = 1; i <= nsMat.nbRows; ++i) {
    	for (var j = 1; j <= nsMat.nbColumns; ++j) {
    	  var newVal = nsMat2.getValueAt(i,j) + Math.random();
    	  nsMat2.setValueAt(i,j, newVal);
    	  assert.equal(nsMat2.getValueAt(i,j), newVal, 'Matrix values #2');
    	}
      }
    
      // Ensure the old matrix is not equal to the new matrix with a 0 tolerance...
      assert.equal(PortfolioAllocation.Matrix.areEqual(nsMat, nsMat2), false, 'Old matrix not strictly equal to new matrix');
    
      // ... but is with a 1 tolerance, as Math.random max value is 1, so that the new matrix coefficients must be 
      // equal to the old matrix coefficients + at most 1 !
      assert.equal(PortfolioAllocation.Matrix.areEqual(nsMat, nsMat2, 1), true, 'Old matrix equal to new matrix within 1');
  }
  
  // Test square matrix using random data
  {
      var sMat = new PortfolioAllocation.Matrix(this.sMatValues);
    
      // Ensure the matrix is square
      assert.equal(sMat.isSquare(), true, 'Matrix is square');
      
      // Extract the matrix diagonal elements
      var diag = sMat.diagonal();
      for (var i = 1; i <= sMat.nbRows; ++i) {
          assert.equal(diag.getValueAt(i, 1), sMat.getValueAt(i, i), 'Matrix diagonal elements');
      }
  }

});

QUnit.test('Matrix-matrix product', function(assert) {    
  // Test matrix-matrix product using static data
  var mat1 = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var mat2 = new PortfolioAllocation.Matrix([[1,2], [2,4], [3,6]]);
  var prodMat = PortfolioAllocation.Matrix.product(mat1, mat2);

  var expectedResMat = new PortfolioAllocation.Matrix([[14,28], [32,64]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(prodMat, expectedResMat), true, 'Matrix-matrix product');
  
  // TODO: use random data
});

QUnit.test('Matrix row - column product', function(assert) {    
  // Test matrix row - column product using static data
  var mat1 = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var mat2 = new PortfolioAllocation.Matrix([[1,2], [2,4], [3,6]]);
  var expectedResMat = new PortfolioAllocation.Matrix([[14,28], [32,64]]);
  
  for (var i = 1; i <= mat1.nbRows; ++i) {
	  for (var j = 1; j <= mat2.nbColumns; ++j) {
			assert.equal(PortfolioAllocation.Matrix.rowColumnProduct(mat1, i, mat2, j), expectedResMat.getValueAt(i,j), 'Matrix row - column product #(' + i + "," + j + ")");
	  }
  }
  
  // TODO: use random data
});


QUnit.test('Matrix element wise function', function(assert) {    
  // Test matrix element wise function using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var powerMat = mat.elemMap(function(i,j,val) { return Math.pow(val, 2); });
 
  var expectedResMat = new PortfolioAllocation.Matrix([[1,4,9], [16,25,36]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(powerMat, expectedResMat), true, 'Matrix element wise function');
});


QUnit.test('Matrix to string', function(assert) {    
  // Test matrix to string using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,10]]);
  assert.equal(mat.toString(), '[  1  2  3 ]\n[  4  5 10 ]\n', 'Matrix to string');
});


QUnit.test('Matrix to row array', function(assert) {    
  // Test using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,10]]);
  assert.deepEqual(mat.toRowArray(), [[1,2,3], [4,5,10]], 'Matrix to double array');
  assert.deepEqual(mat.toRowArray(function(i,j,val) { return i==j; }), [[1], [5]], 'Matrix to doule array with function');
  
  // Test using the random matrix, using identity matrix(matrix to array) == matrix
  var nsMat = new PortfolioAllocation.Matrix(this.nsMatValues);
  var nsMatArray = nsMat.toRowArray();
  assert.deepEqual(nsMatArray, this.nsMatValues, 'Array to matrix to array');
  assert.equal(PortfolioAllocation.Matrix.areEqual(nsMat, new PortfolioAllocation.Matrix(nsMatArray)), true, 'Matrix to double array to matrix');
});


QUnit.test('Matrix to array', function(assert) {    
  // Test using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,10]]);
  assert.deepEqual(mat.toArray(), [1,2,3,4,5,10], 'Matrix to array');
  assert.deepEqual(mat.toArray(function(i,j,val) { return i==j; }), [1, 5], 'Matrix to array with function');  
});


QUnit.test('Diagonal matrix creation', function(assert) {    
  // Test using static data
  var mat = PortfolioAllocation.Matrix.diagonal([1,2,3]);
  var expectedMat = new PortfolioAllocation.Matrix([[1,0,0], [0,2,0], [0,0,3]]);
  assert.deepEqual(mat.toArray(), expectedMat.toArray(), 'Diagonal matrix creation');
});


QUnit.test('Symetric matrix creation', function(assert) {    
  // Test using static data
  var mat = PortfolioAllocation.Matrix.fillSymetric(2, function(i,j) { return i+j; });
  var expectedMat = new PortfolioAllocation.Matrix([[2,3], [3,4]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(mat, expectedMat), true, 'Symetric matrix creation');
});


QUnit.test('Submatrix extraction', function(assert) {    
  // Test using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6], [7,8,9]]);
  var subMat = mat.submatrix([1,3], [1, 3]);
  var expectedMat = new PortfolioAllocation.Matrix([[1,3], [7,9]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(subMat, expectedMat), true, 'Submatrix extraction');
});

QUnit.test('Matrix row extraction', function(assert) {    
  // Test using static data
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var expected_row_1 =  new PortfolioAllocation.Matrix([1,2,3]);
  var expected_row_2 =  new PortfolioAllocation.Matrix([4,5,6]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(mat.row(1), expected_row_1), true, 'Matrix row extraction #1');
  assert.equal(PortfolioAllocation.Matrix.areEqual(mat.row(2), expected_row_2), true, 'Matrix row extraction #2');
});

QUnit.test('Zeros matrix creation', function(assert) {    
  // Test using static data
  var mat = PortfolioAllocation.Matrix.zeros(3, 2);
  var expectedMat = new PortfolioAllocation.Matrix([[0,0], [0,0], [0,0]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(mat, expectedMat), true, 'Zeros matrix creation');
});

QUnit.test('Identity matrix creation', function(assert) {    
  // Test using static data
  var mat = PortfolioAllocation.Matrix.identity(3);
  var expectedMat = new PortfolioAllocation.Matrix([[1,0,0], [0,1,0], [0,0,1]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(mat, expectedMat), true, 'Identity matrix creation');
});

QUnit.test('Transpose matrix', function(assert) {    
  // Test using static data  
  var mat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var transpMat = mat.transpose();
  var expectedMat = new PortfolioAllocation.Matrix([[1,4], [2,5], [3,6]]); 
  assert.equal(PortfolioAllocation.Matrix.areEqual(transpMat, expectedMat), true, 'Transpose matrix');
});


QUnit.test('Givens QR decomposition', function(assert) {    
  // Test using static data  
  {
	  var mat = new PortfolioAllocation.Matrix([[1,2,0], [1,1,1], [2,1,0]]);
	  
	  // Computation of a QR decomposition
	  var qr = PortfolioAllocation.Matrix.qrDecomposition(mat);
	  var q = qr[0];
	  var r = qr[1];
	  var qqp = PortfolioAllocation.Matrix.product(q, q.transpose());  
	  var qtimesr = PortfolioAllocation.Matrix.product(q, r);
	  
	  // Expected matrices were verified with Matlab
	  var expectedQMat = new PortfolioAllocation.Matrix([[1,4], [2,5], [3,6]]); 
	  var expectedRMat = new PortfolioAllocation.Matrix([[1,4], [2,5], [3,6]]);
	  
	  assert.equal(PortfolioAllocation.Matrix.areEqual(expectedQMat, expectedRMat), true, 'Givens QR decomposition');
  }
  
  // TODO: Test using random data: check Q,R dimensions, check Q*R = A, check R upper triangular, check Q orthogonal: Q*q^t = Identity (m)
  
  // TODO: Test error case
});


QUnit.test('Determinant computation', function(assert) {    
  // Test using static data  
  var mat = new PortfolioAllocation.Matrix([[-2,2,-3], [-1,1,3], [2,0,-1]]);
  var expectedValue = 18;

  assert.equal(Math.abs(mat.determinant() - expectedValue) <= 1e-16, true, 'Determinant computation');
  
  // TODO: Test random data
});


QUnit.test('Matrix copy', function(assert) {    
  // Test using static data
  {
    // Original matrix to be copied, and expected copy
    var originalMat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
    var expectedCopyMat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
    
    // First copy: standard one
    var copyMat = PortfolioAllocation.Matrix.copy(originalMat);
    
    // Second copy: using an already existing matrix
    var outputMat = PortfolioAllocation.Matrix.zeros(originalMat.nbRows, originalMat.nbColumns);
    var outputCopyMat = PortfolioAllocation.Matrix.copy(originalMat, outputMat);
    
    // Alteration of the original matrix to check that the copy is a real copy and not a reference
    originalMat.setValueAt(1, 1, 0);
    var originalUpdatedMat = new PortfolioAllocation.Matrix([[0,2,3], [4,5,6]]);
    assert.equal(PortfolioAllocation.Matrix.areEqual(originalMat, originalUpdatedMat), true, 'Matrix copy #1/0');
    
    // Tests of the method behaviour
    assert.equal(PortfolioAllocation.Matrix.areEqual(copyMat, expectedCopyMat), true, 'Matrix copy #1/3');
    assert.equal(PortfolioAllocation.Matrix.areEqual(outputMat, outputCopyMat), true, 'Matrix copy #1/3');
    assert.equal(PortfolioAllocation.Matrix.areEqual(outputCopyMat, expectedCopyMat), true, 'Matrix copy #1/2');
  }
  
  // TODO: Test random data
});


QUnit.test('Matrix AXPBY', function(assert) {    
  // Test using static data
  {
    // 
    var mat1 = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
	var a = -1;
    var mat2 = new PortfolioAllocation.Matrix([[7,8,9], [10,11,12]]);
	var b = 1;
	var expectedMat = new PortfolioAllocation.Matrix([[6,6,6], [6,6,6]]);
	
	// First computation: standard one
	var computedMat = PortfolioAllocation.Matrix.axpby(a, mat1, b, mat2);
	
	// Second computation: using an already existing matrix
    var outputMat = PortfolioAllocation.Matrix.zeros(mat1.nbRows, mat1.nbColumns);
	var outputComputedMat = PortfolioAllocation.Matrix.axpby(a, mat1, b, mat2, outputMat);
	
	// Tests of the method behaviour
	assert.equal(PortfolioAllocation.Matrix.areEqual(computedMat, expectedMat), true, 'Matrix axpby #1/1');
	assert.equal(PortfolioAllocation.Matrix.areEqual(outputMat, outputComputedMat), true, 'Matrix axpby #1/2');
	assert.equal(PortfolioAllocation.Matrix.areEqual(outputComputedMat, expectedMat), true, 'Matrix axpby #1/3');
  }
  
  // TODO: Test random data
});


QUnit.test('Matrix norms computation', function(assert) {    
  // Unsupported norm
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,2,8]]);
      assert.throws(function() { mat.matrixNorm(2) },
		                         new Error('unsupported matrix norm: 2'),
		                         "Matrix norm computation - unsupported norm");
  }
		
  // Test using static data  
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,2,8]]);
      var expectedNorm1 = 19;
      var expectedNormInf = 15;
      var expectedNormFrobenius = Math.sqrt(3*3 + 5*5 + 7*7 + 2*2 + 6*6 + 4*4 + 2*2 + 8*8);

      assert.equal(mat.matrixNorm('one'), expectedNorm1, 'Matrix 1-norm computation #1');
      assert.equal(mat.matrixNorm('infinity'), expectedNormInf, 'Matrix infinity-norm computation #1');
      assert.equal(mat.matrixNorm('frobenius'), expectedNormFrobenius, 'Matrix Frobenius-norm computation #1');
  }

  // Test using static data  
  {
      var mat = new PortfolioAllocation.Matrix([[0,0], [0,0]]);
      var expectedNorm1 = 0;
      var expectedNormInf = 0;
      var expectedNormFrobenius = 0;

      assert.equal(mat.matrixNorm('one'), expectedNorm1, 'Matrix 1-norm computation #2');
      assert.equal(mat.matrixNorm('infinity'), expectedNormInf, 'Matrix infinity-norm computation #2');
      assert.equal(mat.matrixNorm('frobenius'), expectedNormFrobenius, 'Matrix Frobenius-norm computation #2');
  }
  
  // TODO: Test random data
});

QUnit.test('Matrix vector norms computation', function(assert) {    
  // Unsupported norm
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,2,8]]);
      assert.throws(function() { mat.vectorNorm(2) },
		                         new Error('unsupported vector norm: 2'),
		                         "Matrix vector norm computation - unsupported norm");
  }
		
  // Test using static data  
  {
      var mat = new PortfolioAllocation.Matrix([[-1,2], [3,-4]]);
      var expectedNorm1 = 10;
      var expectedNormInf = 4;
      var expectedNorm2 = Math.sqrt(30);

      assert.equal(mat.vectorNorm('one'), expectedNorm1, 'Matrix vector 1-norm computation');
      assert.equal(mat.vectorNorm('infinity'), expectedNormInf, 'Matrix vector infinity-norm computation');
      assert.equal(mat.vectorNorm('two'), expectedNorm2, 'Matrix vector 2-norm computation');
  }
  
  // TODO: Test random data
});


QUnit.test('Matrix row norms computation', function(assert) {    
  // Unsupported norm
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,2,8]]);
      assert.throws(function() { mat.vectorNorm(2, 'row', 1) },
		                         new Error('unsupported vector norm: 2'),
		                         "Matrix row norm computation - unsupported norm");
  }
		
  // Test using static data  
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,0,0]]);
      var expectedRowsNorm1 = [15, 12, 0];
      var expectedRowsNormInf = [7, 6, 0];
      var expectedRowsNorm2 = [Math.sqrt(3*3 + 5*5 + 7*7), Math.sqrt(2*2 + 6*6 +4*4), 0];

	  for(var i = 0; i < mat.nbRows; ++i) {
		assert.equal(mat.vectorNorm('one', 'row', i+1), expectedRowsNorm1[i], 'Matrix row 1-norm computation');
		assert.equal(mat.vectorNorm('infinity', 'row', i+1), expectedRowsNormInf[i], 'Matrix row infinity-norm computation');
		assert.equal(Math.abs(mat.vectorNorm('two', 'row', i+1) - expectedRowsNorm2[i]) <= 1e-15, true, 'Matrix row 2-norm computation');
	  }
  }
  
  // TODO: Test random data
});


QUnit.test('Matrix linear system solving via Kaczmarz algorithms', function(assert) {    
  // Limit case: null matrix; in this case, the least square solution is null as well
  {
	  var A = new PortfolioAllocation.Matrix([[0,0,0], [0,0,0], [0, 0, 0]]);
	  var b = new PortfolioAllocation.Matrix([1,2,3]);
      var expectedX = new PortfolioAllocation.Matrix([0, 0, 0]);

	  var x = PortfolioAllocation.Matrix.linsolveKaczmarz(A, b);
	  assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-14), true, 'Linear system solve - Kaczmarz algorithm #1');
  }
  
  // Test using static data
  {
      // Source: https://en.wikipedia.org/wiki/System_of_linear_equations
	  var A = new PortfolioAllocation.Matrix([[3,2,-1], [2,-2,4], [-1, 0.5, -1]]);
	  var b = new PortfolioAllocation.Matrix([1,-2,0]);
      var expectedX = new PortfolioAllocation.Matrix([1, -2, -2]);

	  var x = PortfolioAllocation.Matrix.linsolveKaczmarz(A, b);
	  assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-14), true, 'Linear system solve - Kaczmarz algorithm #2');
  }
  
  // TODO: Test random data
});
 