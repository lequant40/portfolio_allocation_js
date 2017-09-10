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
				   typeof matCov.getVariancesVector == 'function', true, 'Covariances matrices methods added');
  }
  
});
