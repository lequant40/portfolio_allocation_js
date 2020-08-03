// ------------------------------------------------------------
QUnit.module('Covariance matrix module', {
  before: function() {
    // 
  }
});


QUnit.test('Covariance matrix native creation', function(assert) {    
  // Static data
  {
	  var mat =  new PortfolioAllocation.Matrix([[1, 0], [0,1]]);
	  console.log(mat)
	  var matCov = mat.toCovarianceMatrix();
	  assert.equal(PortfolioAllocation.Matrix.areEqual(mat,matCov), true, 'Covariances matrices must be equal');
	  assert.equal(typeof matCov.getCorrelationMatrix == 'function' &&
				   typeof matCov.getVariances == 'function' &&
				   typeof matCov.standardizedGeneralizedVariance == 'function', true, 'Covariances matrices methods added');
  }
  
});


QUnit.test('Covariance matrix functions', function(assert) {    
  // Test using static data
  {
	  var mat = new PortfolioAllocation.Matrix([[1, 1, 8.1], [1, 16, 18], [8.1, 18, 81]]).toCovarianceMatrix();

	  // Test correlation matrix extraction
	  var expectedCorrMat = new PortfolioAllocation.Matrix([[1, 0.25, 0.9], [0.25, 1, 0.5], [0.9, 0.5, 1]]);
	  assert.equal(PortfolioAllocation.Matrix.areEqual(mat.getCorrelationMatrix(), expectedCorrMat, 1e-14), true, 'Covariance matrix functions - Correlation matrix extraction');
	  
	  // Test variances vector extraction
	  var expectedVarianceMat = new PortfolioAllocation.Matrix([1, 16, 81]);
	  assert.equal(PortfolioAllocation.Matrix.areEqual(mat.getVariances(), expectedVarianceMat, 1e-14), true, 'Covariance matrix functions - Variances vector extraction');
	  
	  // Test standard deviations vector extraction
	  var expectedStdDevMat = new PortfolioAllocation.Matrix([1, 4, 9]);
	  assert.equal(PortfolioAllocation.Matrix.areEqual(mat.getStandardDeviations(), expectedStdDevMat, 1e-14), true, 'Covariance matrix functions - Standard deviations vector extraction');
  }
});