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
	  var matCov = mat.toCovarianceMatrix();
	  assert.equal(PortfolioAllocation.Matrix.areEqual(mat,matCov), true, 'Covariances matrices must be equal');
	  assert.equal(typeof matCov.getCorrelationMatrix == 'function' &&
				   typeof matCov.getVariancesVector == 'function' &&
				   typeof matCov.standardizedGeneralizedVariance == 'function', true, 'Covariances matrices methods added');
  }
  
});


QUnit.test('Covariance matrix sgv', function(assert) {    
  // Move to DIST
  // Test using static data  
  var mat = new PortfolioAllocation.Matrix([[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]]).toCovarianceMatrix();
  assert.equal(mat.standardizedGeneralizedVariance(), Math.pow(mat.determinant(), 1/4), true, 'Sgv computation');
  
    // TODO: Test using formula and random data
});