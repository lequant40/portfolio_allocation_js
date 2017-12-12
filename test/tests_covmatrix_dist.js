// ------------------------------------------------------------
QUnit.module('Covariance matrix module', {
  before: function() {
    // 
  }
});


QUnit.test('Covariance matrix native creation', function(assert) {    
  // Static data
  {
	  var cov = PortfolioAllocation.covarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
	  assert.deepEqual(cov.toRowArray(), [[0.0003555555555555556,  -0.0005333333333333334], [-0.0005333333333333334, 0.0010666666666666667]], 'Covariance matrix #1');
	  
	  var scov = PortfolioAllocation.sampleCovarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
	  assert.deepEqual(scov.toRowArray(), [[0.0005333333333333335, -0.0008], [-0.0008, 0.0016]], 'Covariance matrix #2');	  
  }
  
});

