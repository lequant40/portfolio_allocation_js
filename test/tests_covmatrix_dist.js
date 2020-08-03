// ------------------------------------------------------------
QUnit.module('Covariance matrix module', {
  before: function() {
    // 
  }
});


QUnit.test('Covariance matrix creation', function(assert) {    
  // Static data
  {
	  var cov = PortfolioAllocation.covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]]);
	  assert.deepEqual(cov.toRowArray(), [[0.0003555555555555556,  -0.0005333333333333334], 
	                                      [-0.0005333333333333334, 0.0010666666666666667]], 'Covariance matrix creation, population covariance');
	  
	  var scov = PortfolioAllocation.covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {method: "sample-covariance"});
	  assert.deepEqual(scov.toRowArray(), [[0.0005333333333333335, -0.0008],  
	                                       [-0.0008, 0.0016]], 'Covariance matrix creation, sample covariance');	  
	  
	  // The computations below were verified using the Matlab script cov1para.m from Ledoit and Wolf
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, 3], [2, 4, 5]], {method: "ledoit-wolf-shrinked-covariance", shrinkageTarget: "constant-variance-null-correlation"});
	  assert.deepEqual(lwcov.toRowArray(), [[0.8193967163039326, 0.6563573883161513],  
	                                        [0.6563573883161513, 1.4028255059182893]], 'Covariance matrix creation, Ledoit-Wolf constant covariance model');
	  
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2], [2, 4], [3, 5]], {method: "ledoit-wolf-shrinked-covariance" , shrinkageTarget: "constant-variance-null-correlation"});
	  assert.deepEqual(lwcov.toRowArray(), [[0.25, 0.5, 0.5], 
	                                        [0.5, 1, 1], 
											[0.5, 1, 1]], 'Covariance matrix creation, Ledoit-Wolf constant covariance model #2');
	  
	  // The computations below were verified using the Matlab script covCor.m from Ledoit and Wolf
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {method: "ledoit-wolf-shrinked-covariance", shrinkageTarget: "constant-correlation"});
	  assert.deepEqual(lwcov.toRowArray(), [[2.5000000000000004, -0.7147951569608585, -3.8283673413445896],
	                                        [-0.7147951569608585, 2.75, -2.6768525072041105 ], 
											[ -3.8283673413445896, -2.6768525072041105, 28.75]], 'Covariance matrix creation, Ledoit-Wolf constant correlation model');
											
	  // The computation below were NOT verified with code from De Nard, pending answer
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {method: "ledoit-wolf-shrinked-covariance", shrinkageTarget: "constant-variance-correlation"});
	  assert.deepEqual(lwcov.toRowArray(), [[6.525679576617755, -1.3292338224681264, -5.4112039932643725],
	                                        [-1.3292338224681264, 6.661745248977629, -2.0095621842675007], 
											[-5.4112039932643725, -2.0095621842675007, 20.812575174404618]], 'Covariance matrix creation, Ledoit-Wolf constant variance covariance model');
											
	  // The computation below were not verified with any code from Schafer
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {method: "ledoit-wolf-shrinked-covariance", shrinkageTarget: "null-correlation"});
	  assert.deepEqual(lwcov.toRowArray(), [[2.5000000000000004, 0, -4.32972972972973],
	                                        [0, 2.75, -0.7216216216216218], 
											[-4.32972972972973, -0.7216216216216218, 28.75]], 'Covariance matrix creation, Ledoit-Wolf null correlation model');
  }

  
});

QUnit.test('Random correlation matrix creation', function(assert) {    
	// Static data
	{
		var n = 3;
		
		// Generate the matrix
		var corr = PortfolioAllocation.randomCorrelationMatrix(n);
	
		// Ensure the correlation matrix is a correlation matrix
		// Check that the matrix is square
		assert.equal(corr.nbRows == n && corr.nbColumns == n, true, 'Random correlation matrix generation, square');
		
		// Check that the matrix is numerically symmetric
		assert.equal(corr.isSymmetric(1e-14), true, 'Random correlation matrix generation, symmetric');
		
		// Check that the matrix has true unit diagonal
		var unit_diagonal = true;
		for (var i = 1; i <= n; ++i) {
			if (corr.getValue(i,i) != 1) {
				unit_diagonal = false;
			}
		}
		assert.equal(unit_diagonal, true, 'Random correlation matrix generation, unit diagonal');

		// Check that the matrix is positive semidefinite
		// Check that the determinant is >= 0
		assert.equal(corr.determinant() >= 0, true, 'Random correlation matrix generation, semi-positive definite');
		
		// Check that all the eigenvalues are >= 0
		// TODO
  }
});

QUnit.test('Random variances vector creation', function(assert) {    
	// Static data
	{
		// Generate the vector
		var variances = PortfolioAllocation.randomVariances(3);
	
		// Ensure the variances are all greater than 0
		assert.equal(variances.isNonNegative(), true, 'Random variances vector creation - Non negative variances');
  }
});

QUnit.test('Perturbed variances vector computation', function(assert) {    
	// Static data
	{
		// Generate the vector
		var variances = PortfolioAllocation.randomVariances(3);
		var variancesVector = variances.toArray();
			
		// The vectors must be of the same size
		var pVariancesVec = PortfolioAllocation.perturbedVariances(variances);
		assert.equal(variances.nbColumns, pVariancesVec.nbColumns, 'Mean vector perturbations, same vector size');
		
		// With k equal to zero, the two vectors must be identical		
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { k: 0});
		var pVariancesVecArray = pVariancesVec.toArray();
		assert.deepEqual(variancesVector, pVariancesVecArray, 'Mean vector perturbations, no perturbation');
		
		// With a small k, the vectors must not be too different (except bad luck)
		// try 1e-3, 1e-2 with relative change !!!
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { k: 0.001});
		pVariancesVecArray = pVariancesVec.toArray();
		var diffOK = true;
		for (var i = 0; i < variancesVector.length; ++i) {
			if (Math.abs((pVariancesVecArray[i] - variancesVector[i])/variancesVector[i]) > 1e-2) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, small perturbation');
		
		// With a big k, the resulting vector must be well above 100
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { k: 50000});
		pVariancesVecArray = pVariancesVec.toArray();
		var diffOK = true;
		for (var i = 0; i < variancesVector.length; ++i) {
			if (Math.abs((pVariancesVecArray[i] - variancesVector[i])/variancesVector[i]) < 100) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, big perturbation');
	
		// Ensure the variances are all greater than 0
		assert.equal(variances.isNonNegative(), true, 'Random variances vector creation - Non negative variances');
  }
});

QUnit.test('Covariance matrix sgv', function(assert) {    
  // Test using random data
	{
		var mat = PortfolioAllocation.randomCorrelationMatrix(4).toCovarianceMatrix();
		assert.equal(mat.standardizedGeneralizedVariance(), Math.pow(mat.determinant(), 1/4), 'Sgv computation');
	}
});