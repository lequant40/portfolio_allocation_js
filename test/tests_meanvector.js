// ------------------------------------------------------------
QUnit.module('Returns computation module', {
  before: function() {
    // 
  }
});


QUnit.test('Mean vector computation', function(assert) {    
  // Static data
  {
	  // Test the sample mean vector computation
	  var meanVec = PortfolioAllocation.meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {method: "sample-mean"}).toArray();
	  var meanVecDefault = PortfolioAllocation.meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]]).toArray();
	  var expectedMeanVec = [0.023333333333333334, -0.010000000000000002];
	  
	  assert.equal(meanVec.length, expectedMeanVec.length, 'Mean vector computation, default parameters #1');
	  assert.equal(meanVecDefault.length, expectedMeanVec.length, 'Mean vector computation, default parameters #2');
	  var resultOK = true;
	  for (var i = 0; i < expectedMeanVec.length; ++i) {
          if (expectedMeanVec[i] != meanVec[i] || expectedMeanVec[i] != meanVecDefault[i]) {
			  resultOK = false;
			  break;
		  }
	  }
	  assert.equal(resultOK, true, 'Mean vector computation, default parameters #3');

	  // Test the DeMiguel shrinkage mean vector computation
	  var meanVec = PortfolioAllocation.meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {method: "demiguel-shrinked-mean"}).toArray();
	  var expectedMeanVec = [0.015659472422062354, -0.0023261390887290207];

	  assert.equal(meanVec.length, expectedMeanVec.length, 'Mean vector computation, DeMiguel shrinkage #1');
	  var resultOK = true;
	  for (var i = 0; i < expectedMeanVec.length; ++i) {
          if (expectedMeanVec[i] != meanVec[i]) {
			  resultOK = false;
			  break;
		  }
	  }
	  assert.equal(resultOK, true, 'Mean vector computation, DeMiguel shrinkage #2');	  
  } 
});


QUnit.test('Mean vector random computation', function(assert) {    
  // Static data
  {
	  //
	  var n = 3;
	  var meanVec = PortfolioAllocation.randomMeanVector(n).toArray();
	  assert.equal(meanVec.length, n, 'Mean vector random computation, vector size');
  } 
});


QUnit.test('Mean vector perturbations', function(assert) {    
	// Static data
	{  
		var meanVec = PortfolioAllocation.meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]])
		var meanVecArray = meanVec.toArray();
		
		// The vectors must be of the same size
		var pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec);
		assert.equal(meanVec.nbColumns, pMeanVec.nbColumns, 'Mean vector perturbations, same vector size');
		
		// With k equal to zero, the two vectors must be identical		
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { k: 0});
		var pMeanVecArray = pMeanVec.toArray();
		assert.deepEqual(meanVecArray, pMeanVecArray, 'Mean vector perturbations, no perturbation');
		
		// With a small k, the vectors must not be too different (except bad luck)
		// try 1e-3, 1e-2 with relative change !!!
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { k: 0.001});
		pMeanVecArray = pMeanVec.toArray();
		var diffOK = true;
		for (var i = 0; i < meanVecArray.length; ++i) {
			if (Math.abs((pMeanVecArray[i] - meanVecArray[i])/meanVecArray[i]) > 1e-2) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, small perturbation');
		
		// With a big k, the resulting vector must be well above 100
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { k: 50000});
		pMeanVecArray = pMeanVec.toArray();
		var diffOK = true;
		for (var i = 0; i < meanVecArray.length; ++i) {
			if (Math.abs((pMeanVecArray[i] - meanVecArray[i])/meanVecArray[i]) < 100) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, big perturbation');
  } 
});
