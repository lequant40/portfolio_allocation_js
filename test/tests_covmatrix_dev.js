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


QUnit.test('Indefinite correlation matrix repair - internal tests', function(assert) {     
	// Static data testing the spectral method
	// Reference: The most general methodology to create a valid correlation matrix for risk management and option pricing purposes, Rebonato, Jackel
	{
		// Define the broken correlation matrix, validated thanks to the provided eigenvalues
		var brokenCorrMat = [[1.000, -0.022, -0.081, -0.845, 0.779, 0.066, -0.051, 0.003, -0.086, 0.041, 0.048, -0.067],
		                     [-0.022, 1.000, -0.140, -0.077, 0.114, -0.096, -0.033, -0.083, 0.028, 0.045, 0.060, -1.485],
							 [-0.081, -0.140, 1.000, 0.095, -0.014, -0.046, -0.010, -0.021, -0.049, 0.002, -0.130, 0.007],
							 [-0.845, -0.077, 0.095, 1.000, -0.067, 0.063, 0.056, -0.059, 0.067, 0.040, -0.068, 0.027],
							 [0.779, 0.114, -0.014, -0.067, 1.000, -0.074, 0.042, 0.004, -0.046, 0.024, 0.043, -0.110],
							 [0.066, -0.096, -0.046, 0.063, -0.074, 1.000, 0.141, 0.027, -0.023, 0.056, -0.043, 0.077],
							 [-0.051, -0.033, -0.010, 0.056, 0.042, 0.141, 1.000, 0.318, -0.027, -0.045, 0.026, 0.020],
							 [0.003, -0.083, -0.021, -0.059, 0.004, 0.027, 0.318, 1.000, -0.065, -0.058, 0.010, 0.003],
							 [-0.086, 0.028, -0.049, 0.067, -0.046, -0.023, -0.027, -0.065, 1.000, -0.002, -0.096, -0.018],
							 [0.041, 0.045, 0.002, 0.040, 0.024, 0.056, -0.045, -0.058, -0.002, 1.000, -0.170, -1.072],
							 [0.048, 0.060, -0.130, -0.068, 0.043, -0.043, 0.026, 0.010, -0.096, -0.170, 1.000, -0.035],
							 [-0.067, -1.485, 0.007, 0.027, -0.110, 0.077, 0.020, 0.003, -0.018, -1.072, -0.035, 1.000]];
		
		// Fix it
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "spectral"});	
			
		// Compare the result with the expected result, but since this one is truncated to 4 digits in the reference,
		// the two are necessarily different, so that the norm of the difference between the two matrices is checked.
		var expectedB = new PortfolioAllocation.Matrix([[-0.199, -0.963, 0.038, 0.160, 0.056, -0.011, -0.019, 0.035, -0.032, -0.002, 0.000, 0.000],
                                                        [-0.854, 0.121, -0.028, -0.326, -0.036, -0.098, -0.027, 0.266, 0.253, 0.048, 0.000, 0.000],
						                                [0.090, 0.131, 0.111, 0.452, -0.646, -0.047, -0.001, 0.522, -0.252, 0.055, 0.000, 0.000],
						                                [0.146, 0.755, -0.054, 0.056, -0.127, -0.084, -0.553, -0.217, 0.070, 0.144, 0.000, 0.000],
						                                [-0.223, -0.663, -0.013, 0.143, -0.191, -0.256, -0.581, -0.177, 0.083, 0.115, 0.000, 0.000],
						                                [0.087, 0.010, -0.349, 0.317, 0.547, 0.422, -0.284, 0.404, 0.038, 0.218, 0.000, 0.000],
						                                [0.068, 0.016, -0.795, 0.033, -0.033, -0.210, -0.122, 0.092, -0.021, -0.541, 0.000, 0.000],
						                                [0.071, -0.069, -0.733, -0.007, -0.155, -0.217, 0.346, -0.155, -0.044, 0.486, 0.000, 0.000],
						                                [-0.011, 0.157, 0.201, 0.031, 0.494, -0.731, -0.022, 0.170, -0.349, 0.062, 0.000, 0.000],
						                                [-0.643, 0.170, -0.025, 0.502, 0.090, 0.194, 0.070, -0.359, -0.344, -0.080, 0.000, 0.000],
						                                [-0.012, -0.149, -0.112, -0.690, -0.092, 0.277, -0.242, 0.038, -0.584, 0.050, 0.000, 0.000],
						                                [0.983, -0.163, 0.066, -0.004, 0.030, -0.014, -0.019, -0.036, 0.027, -0.016, 0.000, 0.000]]);
		var expectedFixedCorrMat = PortfolioAllocation.Matrix.axty(1, expectedB, expectedB); // C.f. the formula in the reference.

		assert.equal(PortfolioAllocation.Matrix.xmy(expectedFixedCorrMat, fixedCorrMat).matrixNorm('infinity') <= 6e-3, true, 'Indefinite correlation matrix repair - Rebonato method, internal test #1');
  }
  
  // Test with random data that the linear shrinkage is correctly taking into account the lower bound on eigenvalues
  {
	 // The min eiganvalue of this matrix is ~ -0.0073
	 var brokenCorrMat = [[1, 0.9, 0.7],
						 [0.9, 1, 0.3],
						 [0.7, 0.3, 1]];
							 
	 // Generate a random lower bound on the smallest eigenvalue of the repaired correlation matrix
	 var minEigVal = Math.random(); // belongs to ]0,1[
	 
	 // Fix the broken correlation matrix
	 var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "linear-shrinkage", minEigenvalue: minEigVal})
	
	 // Validate that the smallest eigenvalue of the fixed correlation matrix is greater than 
	 // the lower bound
	 var jacobi = PortfolioAllocation.Matrix.eig(fixedCorrMat, {sortedEigenvalues: true});
	 var minEigenvalue = jacobi[1].data[2];
	 assert.equal(Math.abs(minEigenvalue - minEigVal) <= 1e-4, true, 'Indefinite correlation matrix repair - Shrinkage method, internal test #1');
  }
  
  // Test with random data that the nearest correlation method is correctly taking into account the lower bound on eigenvalues
  {
	 // The min eiganvalue of this matrix is ~ -0.0073
	 var brokenCorrMat = [[1, 0.9, 0.7],
						 [0.9, 1, 0.3],
						 [0.7, 0.3, 1]];
							 
	 // Generate a random lower bound on the smallest eigenvalue of the repaired correlation matrix
	 var minEigVal = Math.random(); // belongs to ]0,1[
	 
	 // Fix the broken correlation matrix
	 var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {minEigenvalue: minEigVal, epsNcm: 1e-6})
	
	 // Validate that the smallest eigenvalue of the fixed correlation matrix is greater than 
	 // the lower bound
	 var jacobi = PortfolioAllocation.Matrix.eig(fixedCorrMat, {sortedEigenvalues: true});
	 var minEigenvalue = jacobi[1].data[2];
	 assert.equal(Math.abs(minEigenvalue - minEigVal) <= 1e-5, true, 'Indefinite correlation matrix repair - Nearest correlation method, internal test #1');
  }
});




QUnit.test('Nearest correlation matrix - internal tests', function(assert) {     
  // Test with random data that the nearest correlation matrix is correctly taking into account the lower bound on eigenvalues
  {
	 // The min eiganvalue of this matrix is ~ -0.0073
	 var brokenCorrMat = [[1, 0.9, 0.7],
						 [0.9, 1, 0.3],
						 [0.7, 0.3, 1]];
							 
	 // Generate a random lower bound on the smallest eigenvalue of the nearest correlation matrix
	 var minEigVal = Math.random(); // belongs to ]0,1[
	 
	 // Compute the nearest correlation matrix to the broken correlation matrix
	 var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(brokenCorrMat, {minEigenvalue: minEigVal, eps: 1e-8, maxIter: 100000})
	
	 // Validate that the smallest eigenvalue of the nearest correlation matrix is greater than 
	 // the lower bound
	 var jacobi = PortfolioAllocation.Matrix.eig(nearestCorrMat, {sortedEigenvalues: true});
	 var minEigenvalue = jacobi[1].data[2];
	 assert.equal(Math.abs(minEigenvalue - minEigVal) <= 1e-6, true, 'Nearest correlation matrix - Internal test #1');
  }
  
  // Test matrices whose nearest correlation matrix were initially NOT a correlation matrix (strictly negative eigenvalue) !
  //
  // The bug was due to an incorrect convergence condition.
  {
	  //
	  var mat = [[1,                    0.618421696258621  ,  0.4874298723751651,  -0.49079033902258856  , 0.11944609519358924 ,  -0.5578129586552424  , -0.4771156006785486 ,  -0.5891536415285011 ,  -0.3033661105868079 ,    0.346969569262485],
				[0.618421696258621,                     1  ,  0.5625061910498927 ,  -0.3066061449896621 ,  0.22036065346794634 ,  -0.2995526904104299 ,  -0.2665281748643263 ,  -0.4384779718501961  ,-0.04157807243431871 ,  0.36441921116789094],
				[0.4874298723751651 ,   0.5625061910498927  ,                   1,  -0.24898733591796773 ,  0.48207046028556727 ,  -0.3253014452697968 , -0.39490635735791035, -0.018501626749401468 , 0.034978192975079626,  -0.13605230462054563],
				[-0.49079033902258856,   -0.3066061449896621 , -0.24898733591796773 ,                    1,    0.3523200216615302 ,   0.5874187797438452,   0.35430037283032667,   0.28094692738663546  ,  0.5704634723157168,   -0.5160250884331637],
				[0.11944609519358924 ,  0.22036065346794634 ,  0.48207046028556727 ,   0.3523200216615302 ,                    1, -0.011807320173530558 ,  -0.1845706823452529 ,  0.48127181356733134 ,   0.5240171002774565 , -0.04905664181384413],
				[-0.5578129586552424 ,  -0.2995526904104299,   -0.3253014452697968 ,   0.5874187797438452, -0.011807320173530558,                     1 ,   0.9395942716104612 ,   0.3019431882153892 ,   0.5712314628087278 ,  -0.5702138190028966],
				[-0.4771156006785486,   -0.2665281748643263,  -0.39490635735791035 ,  0.35430037283032667,   -0.1845706823452529 ,   0.9395942716104612 ,                    1,    0.3594948044435324 ,  0.40258020819009377,  -0.23665358042668136],
				[-0.5891536415285011,   -0.4384779718501961, -0.018501626749401468 ,  0.28094692738663546,   0.48127181356733134 ,   0.3019431882153892,    0.3594948044435324 ,                    1,     0.307071908776885,  -0.27512543932175526],
				[-0.3033661105868079,  -0.04157807243431871,  0.034978192975079626 ,   0.5704634723157168,    0.5240171002774565 ,   0.5712314628087278,   0.40258020819009377,     0.307071908776885,                     1,  -0.12697231187851038],
				[0.346969569262485,   0.36441921116789094 , -0.13605230462054563  , -0.5160250884331637 , -0.04905664181384413  , -0.5702138190028966,  -0.23665358042668136,  -0.27512543932175526 , -0.12697231187851038  ,                   1]]
	  var mat = PortfolioAllocation.Matrix(mat);
	  
	  // Ensure the matrix is not a correlation matrix
	  assert.equal(mat.isCorrelationMatrix(), false, 'Nearest correlation matrix - Internal test #2/1');
		
	  // Compute the nearest correlation matrix, which must be a correlation matrix
	  var nearest_mat = PortfolioAllocation.nearestCorrelationMatrix(mat);
	  
	  // Ensure the nearest correlation matrix is a correlation matrix
	  assert.equal(nearest_mat.isCorrelationMatrix(1e-5), true, 'Nearest correlation matrix - Internal test #2/2');
  }
  

	// Test with static data a small matrix, which needs more than 800 iterations to converge
	// Reference: Anderson acceleration of the alternating projections method for computing the nearest correlation matrix, Nicholas J. Higham, Natasa Strabic
	// n = 6
	{
		var mat = [[0.010712, 0.000654, 0.002391, 0.010059, -0.008321, 0.001738],
					[0.000654, 0.000004, 0.002917, 0.000650, 0.002263, 0.002913],
					[0.002391, 0.002917, 0.013225, -0.000525, 0.010834, 0.010309],
					[0.010059, 0.000650, -0.000525, 0.009409, -0.010584, -0.001175],
					[-0.008321, 0.002263, 0.010834, -0.010584, 0.019155, 0.008571],
					[0.001738, 0.002913, 0.010309, -0.001175, 0.008571, 0.007396]];
		var mat = PortfolioAllocation.Matrix(mat).toCovarianceMatrix().getCorrelationMatrix();

		// Ensure the matrix is not a correlation matrix
		assert.equal(mat.isCorrelationMatrix(), false, 'Nearest correlation matrix - Internal test #3/1');
		  
		// Compute the nearest correlation matrix at 1e-14
		assert.throws(function() { PortfolioAllocation.nearestCorrelationMatrix(mat, {eps:1e-14}); },
						 new Error('maximum number of iterations reached: 100'),
						 "Nearest correlation matrix - Internal test #3/2");
		var nearest_mat = PortfolioAllocation.nearestCorrelationMatrix(mat, {eps:1e-14,maxIter: 900});
		  
		// Ensure the nearest correlation matrix is a correlation matrix at at least 1e-12
		assert.equal(nearest_mat.isCorrelationMatrix(1e-12), true, 'Nearest correlation matrix - Internal test #3/2');
		
		
		// Compute the nearest correlation matrix at full precision
		var nearest_mat = PortfolioAllocation.nearestCorrelationMatrix(mat, {eps:0,maxIter: 900});
		  
		// Ensure the nearest correlation matrix is a strict correlation matrix
		assert.equal(nearest_mat.isCorrelationMatrix(), true, 'Nearest correlation matrix - Internal test #3/3');
	}
});			