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
  }
});


QUnit.test('Matrix basic manipulations', function(assert) {    
  // Test matrix creation using random data
  var nsMat = PortfolioAllocation.Matrix_(this.nsMatValues);

  // Check the dimensions	
  assert.equal(nsMat.nbRows, this.nsMatValues.length, 'Number of rows');
  assert.equal(nsMat.nbColumns, this.nsMatValues[0].length, 'Number of columns');
  
  // Ensure the matrix is not square
  assert.equal(nsMat.isSquare(), false, 'Matrix not square');
  
  // Ensure the matrix is equal to itself
  assert.equal(PortfolioAllocation.matrixIdentical_(nsMat,nsMat), true, 'Matrix equal to itself');
  
  // Check the values with getValueAt
  for (var i = 0; i < nsMat.nbRows; ++i) {
	for (var j = 0; j < nsMat.nbColumns; ++j) {
	  assert.equal(nsMat.getValueAt(i+1,j+1), this.nsMatValues[i][j], 'Matrix values');
	}
  }
  
  // Create a new matrix equal to the previous one and replace 
  // the values of the matrix with new values with setValueAt
  var nsMat2 = PortfolioAllocation.Matrix_(this.nsMatValues);
  for (var i = 1; i <= nsMat.nbRows; ++i) {
	for (var j = 1; j <= nsMat.nbColumns; ++j) {
	  var newVal = nsMat2.getValueAt(i,j) + Math.random();
	  nsMat2.setValueAt(i,j, newVal);
	  assert.equal(nsMat2.getValueAt(i,j), newVal, 'Matrix values #2');
	}
  }

  // Ensure the old matrix is not equal to the new matrix with a 0 tolerance...
  assert.equal(PortfolioAllocation.matrixIdentical_(nsMat, nsMat2), false, 'Old matrix not strictly equal to new matrix');

  // ... but is with a 1 tolerance, as Math.random max value is 1, so that the new matrix coefficients must be 
  // equal to the old matrix coefficients + at most 1 !
  assert.equal(PortfolioAllocation.matrixIdentical_(nsMat, nsMat2, 1), true, 'Old matrix equal to new matrix within 1');
});


QUnit.test('Matrix-matrix product', function(assert) {    
  // Test matrix-matrix product using static data
  var mat1 = PortfolioAllocation.Matrix_([[1,2,3], [4,5,6]]);
  var mat2 = PortfolioAllocation.Matrix_([[1,2], [2,4], [3,6]]);
  var prodMat = PortfolioAllocation.matrixMatrixProduct_(mat1, mat2);

  var expectedResMat = PortfolioAllocation.Matrix_([[14,28], [32,64]]);
  assert.equal(PortfolioAllocation.matrixIdentical_(prodMat, expectedResMat), true, 'Matrix-matrix product');
});


QUnit.test('Matrix element wise power', function(assert) {    
  // Test matrix element wise power using static data
  var mat = PortfolioAllocation.Matrix_([[1,2,3], [4,5,6]]);
  var powerMat = mat.elementWisePower(2);
 
  var expectedResMat = PortfolioAllocation.Matrix_([[1,4,9], [16,25,36]]);
  assert.equal(PortfolioAllocation.matrixIdentical_(powerMat, expectedResMat), true, 'Matrix element wise power');
});


QUnit.test('Matrix to string', function(assert) {    
  // Test matrix to string using static data
  var mat = PortfolioAllocation.Matrix_([[1,2,3], [4,5,10]]);
  assert.equal(mat.toString(), '[  1  2  3 ]\n[  4  5 10 ]\n', 'Matrix to string');
});


QUnit.test('Matrix to array', function(assert) {    
  // Test matrix to array using static data
  var mat = PortfolioAllocation.Matrix_([[1,2,3], [4,5,10]]);
  assert.deepEqual(mat.toArray(), [[1,2,3], [4,5,10]], 'Matrix to array');
  
  // Test matrix to array using the random matrix, using identity matrix(matrix to array) == matrix
  var nsMat = PortfolioAllocation.Matrix_(this.nsMatValues);
  var nsMatArray = nsMat.toArray();
  assert.deepEqual(nsMatArray, this.nsMatValues, 'Array to matrix to array');
  assert.equal(PortfolioAllocation.matrixIdentical_(nsMat, PortfolioAllocation.Matrix_(nsMatArray)), true, 'Matrix to array to matrix');
});

