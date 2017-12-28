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

  var axyMat = PortfolioAllocation.Matrix.axy(2, mat1, mat2);
  var expectedResMat = new PortfolioAllocation.Matrix([[28,56], [64, 128]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(axyMat, expectedResMat), true, 'Matrix-matrix AXY');

  var mat3 = mat1.transpose();
  var atxyMat = PortfolioAllocation.Matrix.atxy(2, mat3, mat2);
  assert.equal(PortfolioAllocation.Matrix.areEqual(atxyMat, expectedResMat), true, 'Matrix-matrix ATXY');

  var mat4 = mat2.transpose();
  var axtyMat = PortfolioAllocation.Matrix.axty(2, mat1, mat4);
  assert.equal(PortfolioAllocation.Matrix.areEqual(axtyMat, expectedResMat), true, 'Matrix-matrix AXTY');
  
  // TODO: use random data
});


QUnit.test('Matrix-matrix elementwise product', function(assert) {    
  // Test matrix-matrix elementwise product using static data
  var mat1 = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);
  var mat2 = new PortfolioAllocation.Matrix([[1,2,3]]);
  var mat3 = new PortfolioAllocation.Matrix([[1],[4]]);
  var mat4 = new PortfolioAllocation.Matrix([[1,2,3], [4,5,6]]);

  var prodMat = PortfolioAllocation.Matrix.elementwiseProduct(mat1, mat2);
  var expectedResMat = new PortfolioAllocation.Matrix([[1,4,9], [4,10,18]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(prodMat, expectedResMat), true, 'Matrix-matrix elementwise product - Row matrix');

  var prodMat = PortfolioAllocation.Matrix.elementwiseProduct(mat1, mat3);
  var expectedResMat = new PortfolioAllocation.Matrix([[1,2,3], [16,20,24]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(prodMat, expectedResMat), true, 'Matrix-matrix elementwise product - Column matrix');

  var prodMat = PortfolioAllocation.Matrix.elementwiseProduct(mat1, mat4);
  var expectedResMat = new PortfolioAllocation.Matrix([[1,4,9], [16,25,36]]);
  assert.equal(PortfolioAllocation.Matrix.areEqual(prodMat, expectedResMat), true, 'Matrix-matrix elementwise product - Full matrix');

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


QUnit.test('Marix negativity and positivity functions', function(assert) {    
  // Test using static data
  {
	  var nothingMat = new PortfolioAllocation.Matrix([[-1,2,3], [4,5,10]]);
	  var posMat = new PortfolioAllocation.Matrix([[1,2,3], [4,5,10]]);
	  var negMat = new PortfolioAllocation.Matrix([[-1,-2,-3], [-4,-5,-10]]);
	  var nonNegMat = new PortfolioAllocation.Matrix([[0,2,3], [4,5,10]]);
	  var nonPosMat = new PortfolioAllocation.Matrix([[0,-2,-3], [-4,-5,-10]]);
	  
	  assert.equal(posMat.isNonNegative(), true, 'Marix negativity and positivity - Test #1/1');
	  assert.equal(posMat.isPositive(), true, 'Marix negativity and positivity - Test #1/2');
	  assert.equal(posMat.isNonPositive(), false, 'Marix negativity and positivity - Test #1/3');
	  assert.equal(posMat.isNegative(), false, 'Marix negativity and positivity - Test #1/4');
	  
	  assert.equal(nonNegMat.isNonNegative(), true, 'Marix negativity and positivity - Test #2/1');
	  assert.equal(nonNegMat.isPositive(), false, 'Marix negativity and positivity - Test #2/2');
	  assert.equal(nonNegMat.isNonPositive(), false, 'Marix negativity and positivity - Test #2/3');
	  assert.equal(nonNegMat.isNegative(), false, 'Marix negativity and positivity - Test #2/4');  
	  
	  assert.equal(nothingMat.isNonNegative(), false, 'Marix negativity and positivity - Test #3/1');
	  assert.equal(nothingMat.isPositive(), false, 'Marix negativity and positivity - Test #3/2');
	  assert.equal(nothingMat.isNonPositive(), false, 'Marix negativity and positivity - Test #3/3');
	  assert.equal(nothingMat.isNegative(), false, 'Marix negativity and positivity - Test #3/4');
	  
	  assert.equal(negMat.isNonNegative(), false, 'Marix negativity and positivity - Test #4/1');
	  assert.equal(negMat.isPositive(), false, 'Marix negativity and positivity - Test #4/2');
	  assert.equal(negMat.isNonPositive(), true, 'Marix negativity and positivity - Test #4/3');
	  assert.equal(negMat.isNegative(), true, 'Marix negativity and positivity - Test #4/4');

	  assert.equal(nonPosMat.isNonNegative(), false, 'Marix negativity and positivity - Test #5/1');
	  assert.equal(nonPosMat.isPositive(), false, 'Marix negativity and positivity - Test #5/2');
	  assert.equal(nonPosMat.isNonPositive(), true, 'Marix negativity and positivity - Test #5/3');
	  assert.equal(nonPosMat.isNegative(), false, 'Marix negativity and positivity - Test #5/4');
  }

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
  var mat = PortfolioAllocation.Matrix.diagonal(new PortfolioAllocation.Matrix([1,2,3]));
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

QUnit.test('Matrix minimum and maximum values computation', function(assert) {    
  // Test using static data  
  {
      var mat = new PortfolioAllocation.Matrix([[-3,5,7], [2,6,4], [0,2,8]]);
      var expectedMinVal = -3;
	  var expectedMaxVal = 8;

      assert.equal(mat.min(), expectedMinVal, 'Matrix minimum computation #1');
	  assert.equal(mat.max(), expectedMaxVal, 'Matrix maximum computation #1');
  }
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


QUnit.test('Matrix row/columns norms computation', function(assert) {    
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
		assert.equal(Math.abs(mat.vectorNorm('two', 'row', i+1) - expectedRowsNorm2[i]) <= 1e-14, true, 'Matrix row 2-norm computation');
	  }
	  
	  
      var expectedColumnsNorm1 = [5, 11, 11];
      var expectedColumnsNormInf = [3, 6, 7];
      var expectedColumnsNorm2 = [Math.sqrt(3*3 + 2*2), Math.sqrt(5*5 + 6*6), Math.sqrt(7*7 + 4*4)];
	  
	  for(var j = 0; j < mat.nbColumns; ++j) {
		assert.equal(mat.vectorNorm('one', 'column', j+1), expectedColumnsNorm1[j], 'Matrix column 1-norm computation');
		assert.equal(mat.vectorNorm('infinity', 'column', j+1), expectedColumnsNormInf[j], 'Matrix column infinity-norm computation');
		assert.equal(Math.abs(mat.vectorNorm('two', 'column', j+1) - expectedColumnsNorm2[j]) <= 1e-14, true, 'Matrix column 2-norm computation');
	  }
  }
  
  // TODO: Test random data
});


QUnit.test('Linear system solver - Kaczmarz algorithm', function(assert) {    
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
	  assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-10), true, 'Linear system solve - Kaczmarz algorithm #2');
  }
  
  // Reference: Tanabe, K.: An algorithm for the constrained maximization in nonlinear programming, Research Memorandum No. 31, The Institute of Statistical Mathematics, 1969.
  // Test using static data
  {
	// Problem 1
	var A = new PortfolioAllocation.Matrix([[-3.2, 2.9, 1.6, 0.1], [0.0, -1.1, 2.3, 1.0], [5.1, 4.8, 0.2, 4.9], [2.0, 1.1, 1.9, -2.9]]);
	var b = new PortfolioAllocation.Matrix([1.4, 2.2, 15.0, 2.1]);
	var expectedX = new PortfolioAllocation.Matrix([1, 1, 1, 1]);
	
	var x = PortfolioAllocation.Matrix.linsolveKaczmarz(A, b);
	assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-14), true, 'Linear system solve - Kaczmarz algorithm #3');
	
	
	// Problem 4
	var A = new PortfolioAllocation.Matrix([[1, 3, 2, -1], [1, 2, -1, -2], [1, -1, 	2, 3], [2, 1, 1, 1], [5, 5, 4, 1], [4, -1, 5, 7]]);
	var b = new PortfolioAllocation.Matrix([5, 0, 5, 5, 15, 15]);
	// To be noted that the true solution in the reference is [1, 1, 1, 1], with A*x - b == 0 and ||x||_2 == 2
	// The expected solution below verifies ||A*x - b||_inf ~= 1e-13, but ||x||_2 ~= 1.96, which is less than 2, hence why it is preferred to the true solution.
	var expectedX = new PortfolioAllocation.Matrix([1.1538461538462281, 0.7692307692307582, 1.153846153846079, 0.7692307692307796]);

	var x = PortfolioAllocation.Matrix.linsolveKaczmarz(A, b);
	assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-14), true, 'Linear system solve - Kaczmarz algorithm #4');
	
	
	// Problem 6
	var A = PortfolioAllocation.Matrix.fill(84, 84, function(i,j) { if (i == j) { return 6; } else if (j == i+1) { return 1; } else if (i == j+1) { return 8; } else { return 0;} });
	var b = PortfolioAllocation.Matrix.fill(84, 1, function(i,j) { if (i == 1) { return 7; } else if (i == 84) { return 14; } else { return 15;} });
	// To be noted that the true solution in the reference is [1, ..., 1], with A*x - b == 0 and ||x||_2 ~> 9.16
	// The expected solution below verifies ||A*x - b||_inf ~= 1e-13, but ||x||_2 ~> 9.15, which is less than 9.16, hence why it is preferred to the true solution.
	// To also be noted that the solution in the reference is almost the same as the expected solution below.
	var expectedX = new PortfolioAllocation.Matrix([1, 0.9999999999999999, 1.0000000000000002, 0.9999999999999996, 1.0000000000000004, 0.9999999999999996, 1.0000000000000004, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999997, 1.0000000000000002, 0.9999999999999999,                  1,                  1, 0.9999999999999999, 1.0000000000000002, 0.9999999999999998, 1.0000000000000002, 0.9999999999999992,  1.000000000000001, 0.9999999999999987,  1.000000000000002, 0.9999999999999974, 1.0000000000000024, 0.9999999999999983, 0.9999999999999991,  1.000000000000007, 0.9999999999999785, 1.0000000000000522, 0.9999999999998842, 1.0000000000002456, 0.9999999999994922,  1.000000000001035, 0.9999999999979075, 1.0000000000042115,  0.999999999991547, 1.0000000000169387, 0.9999999999660868, 1.0000000000678653, 0.9999999998642277, 1.0000000002715892,  0.999999999456776, 1.0000000010864951, 0.9999999978269619, 1.0000000043461248,  0.999999991307702, 1.0000000173846435,  0.999999965230667, 1.0000000695387083,  0.999999860922549, 1.0000002781549087, 0.9999994436902785, 1.0000011126189468, 0.9999977747641977,  1.000004450463144, 0.9999910991076416, 1.0000178016489207, 0.9999643972454112,  1.000071203336108, 0.9998576020201143, 1.0002847611904067, 0.9994306166966817, 1.0011382102966262, 0.9977258046468135, 1.0045394897460913, 0.9909566243489603, 1.0179443359374987, 0.9646809895833339, 1.0683593749999998, 0.8723958333333336, 1.2187499999999998, 0.7083333333333335]);
	
	var x = PortfolioAllocation.Matrix.linsolveKaczmarz(A, b);
	assert.equal(PortfolioAllocation.Matrix.areEqual(x, expectedX, 1e-14), true, 'Linear system solve - Kaczmarz algorithm #5');
  }
  // TODO: Test random data
});
 