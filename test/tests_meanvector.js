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
	  var meanVec = PortfolioAllocation.meanVector([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {method: "linear-shrinkage"}).toArray();
	  var expectedMeanVec = [0.015659472422062354, -0.0023261390887290207];

	  assert.equal(meanVec.length, expectedMeanVec.length, 'Mean vector computation, DeMiguel linear shrinkage #1');
	  var resultOK = true;
	  for (var i = 0; i < expectedMeanVec.length; ++i) {
          if (expectedMeanVec[i] != meanVec[i]) {
			  resultOK = false;
			  break;
		  }
	  }
	  assert.equal(resultOK, true, 'Mean vector computation, DeMiguel linear shrinkage #2');	  
  } 
});


QUnit.test('Mean vector random computation', function(assert) {    
  // Static data
  {
	  // Test the size of the generated vector
	  var n = 3;
	  var meanVec = PortfolioAllocation.randomMeanVector(n).toArray();
	  assert.equal(meanVec.length, n, 'Mean vector random computation, vector size');
	  
	  // Test that the standard deviation has an effect on the generated vector
	  // Test with 0 standard deviation
	  var meanVec = PortfolioAllocation.randomMeanVector(n, 0);
	  assert.equal(meanVec.vectorNorm('infinity') == 0, true, 'Mean vector random computation, 0 standard deviation');
	  
	  // Test with default standard deviation (= 1)
	  var meanVec = PortfolioAllocation.randomMeanVector(n);
	  assert.equal(meanVec.vectorNorm('infinity') <= 5 && meanVec.vectorNorm('infinity') > 0, true, 'Mean vector random computation, default standard deviation');

	  // Test with a high standard deviation
	  var meanVec = PortfolioAllocation.randomMeanVector(n, 200);
	  assert.equal(meanVec.vectorNorm('infinity') > 0  && meanVec.vectorNorm('infinity') >= 50, true, 'Mean vector random computation, high standard deviation');
  } 
});


QUnit.test('Mean vector perturbations', function(assert) {    
	// Static data
	{  
		var meanVec = PortfolioAllocation.randomMeanVector(3);
		var meanVecArray = meanVec.toArray();
		
		// The vectors must be of the same size
		var pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec);
		assert.equal(meanVec.nbColumns, pMeanVec.nbColumns, 'Mean vector perturbations, same vector size');
		
		// With sigma equal to zero, the two vectors must be identical		
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { sigma: 0});
		var pMeanVecArray = pMeanVec.toArray();
		assert.deepEqual(meanVecArray, pMeanVecArray, 'Mean vector perturbations, no perturbation');
		
		// With a small sigma, the vectors must not be too different (except bad luck)
		// try 1e-3, 1e-2 with relative change !!!
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { sigma: 0.001});
		pMeanVecArray = pMeanVec.toArray();
		var diffOK = true;
		for (var i = 0; i < meanVecArray.length; ++i) {
			if (Math.abs((pMeanVecArray[i] - meanVecArray[i])/meanVecArray[i]) > 1e-2) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, small perturbation');
		
		// With a big sigma, the resulting vector must be well above 100
		pMeanVec = PortfolioAllocation.perturbedMeanVector(meanVec, { sigma: 50000});
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
