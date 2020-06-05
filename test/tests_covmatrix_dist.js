// ------------------------------------------------------------
QUnit.module('Covariance matrix module', {
  before: function() {
    // 
  }
});


QUnit.test('Covariance matrix creation', function(assert) {    
  // Static data
  {
	  var cov = PortfolioAllocation.covarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
	  assert.deepEqual(cov.toRowArray(), [[0.0003555555555555556,  -0.0005333333333333334], [-0.0005333333333333334, 0.0010666666666666667]], 'Covariance matrix creation #1');
	  
	  var scov = PortfolioAllocation.sampleCovarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
	  assert.deepEqual(scov.toRowArray(), [[0.0005333333333333335, -0.0008], [-0.0008, 0.0016]], 'Covariance matrix creation #2');	  
  }
});

QUnit.test('Random covariance matrix creation', function(assert) {    
	// Static data
	{
		// Generate the matrix
		var cov = PortfolioAllocation.randomCovarianceMatrix(3);
	
		// Extract the correlation matrix and the variances vector
		var corr = cov.getCorrelationMatrix();
		var variances = cov.getVariances();
		
		// Ensure the correlation matrix is a correlation matrix
		// TODO
		
		// Ensure the variances are all greater than 0
		assert.equal(variances.isNonNegative(), true, 'Random covariance matrix creation, #1 - Non negative variances');
  }
});


QUnit.test('Covariance matrix sgv', function(assert) {    
  // Test using random data
	{
		var mat = PortfolioAllocation.randomCovarianceMatrix(4);
		assert.equal(mat.standardizedGeneralizedVariance(), Math.pow(mat.determinant(), 1/4), 'Sgv computation');
	}
});