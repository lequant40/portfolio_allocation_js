// ------------------------------------------------------------
QUnit.module('Covariance matrix module', {
  before: function() {
    // 
  }
});


QUnit.test('Covariance matrix creation', function(assert) {    
  // Error cases
  {
	assert.throws(function() { 
		var cov = PortfolioAllocation.covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03]]);
	},
	 new Error('inconsistent input matrix dimensions'),
	 "Covariance matrix creation, inconsistent input matrix dimensions");
  }

  // Static data
  {
	  var cov = PortfolioAllocation.covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]]);
	  assert.equal(cov.isCovarianceMatrix(0, "exception"), true, 'Covariance matrix checks');
	  assert.deepEqual(cov.toRowArray(), [[0.0003555555555555556,  -0.0005333333333333334], 
	                                      [-0.0005333333333333334, 0.0010666666666666667]], 'Covariance matrix creation');
	  	  										
	  // The computations below were verified using the Matlab script cov1para.m from Ledoit and Wolf
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, 3], [2, 4, 5]], {regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-variance-null-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[0.8193967163039326, 0.6563573883161513],  
	                                        [0.6563573883161513, 1.4028255059182893]], 'Covariance matrix creation, Ledoit-Wolf constant covariance model');
	  	  
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2], [2, 4], [3, 5]], {regularizationMethod: "linear-shrinkage" , shrinkageTarget: "constant-variance-null-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[0.25, 0.5, 0.5], 
	                                        [0.5, 1, 1], 
											[0.5, 1, 1]], 'Covariance matrix creation, Ledoit-Wolf constant covariance model #2');
	  
	  // The computations below were verified using the Matlab script covCor.m from Ledoit and Wolf
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[2.5000000000000004, -0.7147951569608585, -3.8283673413445896],
	                                        [-0.7147951569608585, 2.75, -2.6768525072041105 ], 
											[ -3.8283673413445896, -2.6768525072041105, 28.75]], 'Covariance matrix creation, Ledoit-Wolf constant correlation model');
											
	  // The computation below were NOT verified with code from De Nard, pending answer
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-variance-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[6.525679576617755, -1.3292338224681264, -5.4112039932643725],
	                                        [-1.3292338224681264, 6.661745248977629, -2.0095621842675007], 
											[-5.4112039932643725, -2.0095621842675007, 20.812575174404618]], 'Covariance matrix creation, Ledoit-Wolf constant variance covariance model');
											
	  // The computation below were not verified with any code from Schafer
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {regularizationMethod: "linear-shrinkage", shrinkageTarget: "null-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[2.5000000000000004, 0, -4.32972972972973],
	                                        [0, 2.75, -0.7216216216216218], 
											[-4.32972972972973, -0.7216216216216218, 28.75]], 'Covariance matrix creation, Ledoit-Wolf null correlation model');
  }
  
  // Static data with non centered data
  {
	  var scov = PortfolioAllocation.covarianceMatrix([[0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]], {assumeZeroMean: false});
	  assert.equal(scov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(scov.toRowArray(), [[0.0005333333333333335, -0.0008],  
	                                       [-0.0008, 0.0016]], 'Covariance matrix creation, no zero mean');

	  // The four linear-shrinkage regularization methods below have NOT been validated against any code
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, 3], [2, 4, 5]], {assumeZeroMean: false, regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-variance-null-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[1.3419243986254297, 0.7306701030927834],  
	                                        [0.7306701030927834	, 1.9914089347079036]], 'Covariance matrix creation, Ledoit-Wolf constant covariance model, no zero mean');
											
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {assumeZeroMean: false, regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[3.333333333333333, -0.9955990054728224, -4.8859840643561725 ],
											[ -0.9955990054728224,   3.666666666666666,  -3.654051336659059 ],
											[ -4.8859840643561725,  -3.654051336659059,   38.33333333333333 ]], 'Covariance matrix creation, Ledoit-Wolf constant correlation model, no zero mean');
											
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {assumeZeroMean: false, regularizationMethod: "linear-shrinkage", shrinkageTarget: "constant-variance-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[   9.302812153826505, -1.9710543275213308,  -6.902628913895051 ],
											[ -1.9710543275213308,   9.467197973372297, -2.7929834252502843 ],
											[  -6.902628913895051, -2.7929834252502843 , 26.563323206134523 ]], 'Covariance matrix creation, Ledoit-Wolf constant variance covariance model, no zero mean');
											
	  var lwcov = PortfolioAllocation.covarianceMatrix([[1, 2, -1, -2], [2, 4, 6, 2], [3, 5, 9, 17]], {assumeZeroMean: false, regularizationMethod: "linear-shrinkage", shrinkageTarget: "null-correlation"});
	  assert.equal(lwcov.isCovarianceMatrix(), true, 'Covariance matrix checks');
	  assert.deepEqual(lwcov.toRowArray(), [[3.333333333333333, 0,  -5.495195195195195 ],
											[0 , 3.6666666666666665, -0.9158658658658657 ],
											[-5.495195195195195, -0.9158658658658657, 38.33333333333333 ]], 'Covariance matrix creation, Ledoit-Wolf null correlation model, no zero mean');
  }

  
});

QUnit.test('Random correlation matrix creation', function(assert) {    
	// Static data
	{
		var n = 3;
		
		// Generate the matrix
		var corr = PortfolioAllocation.randomCorrelationMatrix(n);
	
		// Ensure the correlation matrix is a correlation (and thus also covariance) matrix
		assert.equal(corr.isCorrelationMatrix(), true, 'Random correlation matrix generation, all inclusive test');	
  }
});

QUnit.test('Indefinite correlation matrix repair', function(assert) {    
	// Test error cases
	{
		// Non square
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2, 3], [4, 5, 6]]); },
						 new Error('input matrix must be symmetric'),
						 "Indefinite correlation matrix repair - Non square matrix");

		// Non symmetric 
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [4, 5]]); },
						 new Error('input matrix must be symmetric'),
						 "Indefinite correlation matrix repair - Non symmetric matrix");
		
		// Non unit diagonal
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 2]]); },
						 new Error('input matrix must be unit diagonal'),
						 "Indefinite correlation matrix repair - Non unit diagonal matrix");
		
		// Eigenvalues threshold parameter not belonging to [0,1] interval
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 1]], {minEigenvalue: 1.01}); },
						 new Error('lower bound on the smallest eigenvalue(s) of the correlation matrix must belong to interval [0,1]'),
						 "Indefinite correlation matrix repair - eigenvalues threshold parameter not belonging to [0,1]");				 
						 
		// Unsupported repair method
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 1]], {method: "unknown"}); },
						 new Error('unsupported correlation matrix repair method'),
						 "Indefinite correlation matrix repair - unsupported correlation matrix repair method");
		
		// Shrinking target matrix not a correlation matrix
		var targetMatrix = [[1, 3], [2, 1]];
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 1]], {method: "linear-shrinkage", shrinkageTarget: targetMatrix}); },
						 new Error('shrinkage target matrix must be a correlation matrix'),
						 "Indefinite correlation matrix repair - shrinkage target matrix not symmetric");

		var targetMatrix = [[1, 3], [3, 2]];
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 1]], {method: "linear-shrinkage", shrinkageTarget: targetMatrix}); },
						 new Error('shrinkage target matrix must be a correlation matrix'),
						 "Indefinite correlation matrix repair - shrinkage target matrix not unit diagonal");
						 
		// Shrinking target matrix not satisfying the lower bound on the smallest eigenvalue(s)
		var targetMatrix = [[1, 1], [1, 1]];
		assert.throws(function() { PortfolioAllocation.repairCorrelationMatrix([[1, 2], [2, 1]], {method: "linear-shrinkage", shrinkageTarget: targetMatrix, minEigenvalue: 0.05}); },
						 new Error('smallest eigenvalue of the shrinkage target matrix strictly lower than the desired lower bound on the smallest eigenvalue(s) of the correlation matrix'),
						 "Indefinite correlation matrix repair - shrinkage target matrix not satisfying the lower bound on the smallest eigenvalue(s)");				 

	}
	
	// Static data, with spectral repair
	// Reference: The most general methodology to create a valid correlation matrix for risk management and option pricing purposes, Rebonato, Jackel
	{
		// Define the broken correlation matrix
		var brokenCorrMat = [[1, 0.9, 0.7],
		                     [0.9, 1, 0.3],
							 [0.7, 0.3, 1]];
		
		// Fix it
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "spectral"}).toRowArray();
		
		// Compare the result with the expected result
		var expectedFixedCorrMat = [[1, 0.89402, 0.69632],
		                            [0.89402, 1, 0.30100],
							        [0.69632, 0.30100, 1]];
		
		var matrixEqual = true;
		for (var i = 0; i < 3; ++i) {
			for (var j = 0; j < 3; ++j) {
				if ( Math.abs( fixedCorrMat[i][j] - expectedFixedCorrMat[i][j] ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Spectral method, test #1');
  }
 
 	// Static data, with spectral repair
	// Reference: Fixing a broken correlation matrix, Phil Joubert and Stephen Langdell
	{
		// Define the broken correlation matrix
		var brokenCorrMat = [[1, 0.95, 0],
		                     [0.95, 1, 0.95],
							 [0, 0.95, 1]];
		
		// Fix it
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "spectral"}).toRowArray();
		
		// Compare the result with the expected result
		var expectedFixedCorrMat = [[1, 0.73454, 0.07908],
		                            [0.73454, 1, 0.73454],
							        [0.07908, 0.73454, 1]];
		
		var matrixEqual = true;
		for (var i = 0; i < 3; ++i) {
			for (var j = 0; j < 3; ++j) {
				if ( Math.abs( fixedCorrMat[i][j] - expectedFixedCorrMat[i][j] ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Spectral method, test #2');
  }
  
	// Static data, with shrinkage repair
	// Reference: Restoring Definiteness via Shrinking, with an Application to Correlation Matrices with a Fixed Block, Higham, Strabic, Sego,
	{
		// 
		var n = 5;
		
		// Define the broken correlation matrix
		var brokenCorrMat = [[1.000, 0.900, 0.450, 0.300, 0.225],
                             [0.900, 1.000, 0.900, 0.450, 0.300],
                             [0.450, 0.900, 1.000, 0.900, 0.450],
                             [0.300, 0.450, 0.900, 1.000, 0.900],
                             [0.225, 0.300, 0.450, 0.900, 1.000]];
		
		// Fix the broken correlation matrix by shrinking it towards a target correlation matrix
		var targetCorrMat = [[1.00, 0.90, 0.00, 0.00, 0.00],
                             [0.90, 1.00, 0.00, 0.00, 0.00],
                             [0.00, 0.00, 1.00, 0.00, 0.45],
                             [0.00, 0.00, 0.00, 1.00, 0.45],
                             [0.00, 0.00, 0.45, 0.45, 1.00]];
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "linear-shrinkage", shrinkageTarget: targetCorrMat}).toRowArray();
		
		// Compare the result with the expected result
		var expectedFixedCorrMat = [[1.000, 0.900, 0.343, 0.228, 0.171],
                                    [0.900, 1.000, 0.685, 0.343, 0.228],
                                    [0.343, 0.685, 1.000, 0.685, 0.450],
                                    [0.228, 0.343, 0.685, 1.000, 0.793],
                                    [0.171, 0.228, 0.450, 0.793, 1.000]];
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( fixedCorrMat[i][j] - expectedFixedCorrMat[i][j] ) > 1e-3 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Shrinkage method, test #1 towards a shrinking target matrix');
		
		
		// Fix the broken correlation matrix by shrinking it towards the identity matrix
		var fixedCorrMat2 = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "linear-shrinkage"}).toRowArray();

		// Compare the result with the expected result, which has been calculated using the algorithm,
		// once the matrix above was correct.
		var expectedFixedCorrMat2 = [[1,  0.7656817972961426,  0.3828408986480713, 0.25522726576538085, 0.19142044932403565],
                                     [  0.7656817972961426,                   1,  0.7656817972961426,  0.3828408986480713, 0.25522726576538085 ],
                                     [  0.3828408986480713,  0.7656817972961426,                   1,  0.7656817972961426,  0.3828408986480713 ],
                                     [ 0.25522726576538085,  0.3828408986480713,  0.7656817972961426,                   1,  0.7656817972961426 ],
                                     [ 0.19142044932403565, 0.25522726576538085,  0.3828408986480713,  0.7656817972961426,                   1 ]];
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( fixedCorrMat2[i][j] - expectedFixedCorrMat2[i][j] ) > 1e-6 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Shrinkage method, test #1 towards the identity matrix');
  }
  
  // Static data, with nearest correlation matrix repair
  // Reference: Anderson acceleration of the alternating projections method for computing the nearest correlation matrix, Nicholas J. Higham, Natasa Strabic
  // Same tests as for the Nearest Correlation Matrix tests
  {
		var brokenCorrMat = [[1,     0.18,  -0.13, -0.26, 0.19,  -0.25, -0.12],
							 [0.18,  1,     0.22,  -0.14, 0.31,  0.16,  0.09],
							 [-0.13, 0.22,  1,     0.06,  -0.08, 0.04,  0.04],
							 [-0.26, -0.14, 0.06,  1,     0.85,  0.85,  0.85],
							 [0.19,  0.31,  -0.08, 0.85,  1,     0.85,  0.85],
							 [-0.25, 0.16,  0.04,  0.85,  0.85,  1,     0.85],
							 [-0.12, 0.09,  0.04,  0.85,  0.85,  0.85,  1]];
		var n = brokenCorrMat.length;
		
		// When no minimum eigenvalue is provided nor method, 0 must be the default with nearest correlation matrix method
		var expectedMat =  [[1.000000, 0.183844, -0.131789, -0.251417, 0.178399, -0.247928, -0.119061],
							 [ 0.183844, 1.000000, 0.218241, -0.131564, 0.298599, 0.162037, 0.090923],
							 [-0.131789, 0.218241, 1.000000, 0.056073, -0.074692, 0.039052, 0.039570],
							 [-0.251417, -0.131564, 0.056073, 1.000000, 0.824539, 0.854548, 0.852062],
							 [ 0.178399, 0.298599, -0.074692, 0.824539, 1.000000, 0.843853, 0.847214],
							 [-0.247928, 0.162037, 0.039052, 0.854548, 0.843853, 1.000000, 0.850498],
							 [-0.119061, 0.090923, 0.039570, 0.852062, 0.847214, 0.850498, 1.000000]];
		var fixedCorrMatDefault = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat);
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "nearest-correlation-matrix", minEigenvalue: 0});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - fixedCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ||
				     Math.abs( fixedCorrMat.getValueAt(i+1,j+1) - fixedCorrMatDefault.getValueAt(i+1,j+1) ) > 0 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Nearest correlation matrix method, test #1');
		
		// Test with a minimum eigenvalue of 0.1
		var expectedMat =  [[1.000000, 0.192428, -0.135963, -0.232656, 0.143938, -0.237614, -0.114311],
							[0.192428, 1.000000, 0.213837, -0.106498, 0.269406, 0.163030, 0.093218],
							[-0.135963, 0.213837, 1.000000, 0.045597, -0.060853, 0.037398, 0.038215],
							[-0.232656, -0.106498, 0.045597, 1.000000, 0.759191, 0.851434, 0.855982],
							[0.143938, 0.269406, -0.060853, 0.759191, 1.000000, 0.815831, 0.833260],
							[-0.237614, 0.163030, 0.037398, 0.851434, 0.815831, 1.000000, 0.857818],
							[-0.114311, 0.093218, 0.038215, 0.855982, 0.833260, 0.857818, 1.000000]];
		var fixedCorrMat = PortfolioAllocation.repairCorrelationMatrix(brokenCorrMat, {method: "nearest-correlation-matrix", minEigenvalue: 0.1});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - fixedCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Indefinite correlation matrix repair - Nearest correlation matrix method, test #2');
  }

});


QUnit.test('Perturbed correlation matrix computation', function(assert) {    
	// Static data
	// Test the additive noise method
	{
		// Generate an input correlation matrix
		var n = 10;
		var corrMat = PortfolioAllocation.randomCorrelationMatrix(n);

		// Perturb it
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "additive-noise"});
		
		// Test that the perturbed matrix has the same size as the original matrix
		assert.equal(perturbedCorrMat.nbColumns == perturbedCorrMat.nbColumns && perturbedCorrMat.nbRows == perturbedCorrMat.nbRows, true, 'Perturbed correlation matrix, additive noise method, same matrix size');
		
		// Test that the perturbed matrix is a correlation matrix
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, additive noise method, all inclusive');
		
		// Test that with default max noise level (0.01), the matrix is slightly perturbed
		// Can fail if the matrix had to be repaired
		var matrixEqual = true;
		var matrixNearlyEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( corrMat.getValueAt(i,j) != perturbedCorrMat.getValueAt(i,j) )  {
					matrixEqual = false;
				}

				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 0.01 ) {
					matrixNearlyEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, false, 'Perturbed correlation matrix generation, additive noise method, default standard deviation, matrices different');
		assert.equal(matrixNearlyEqual, true, 'Perturbed correlation matrix generation, additive noise method, default standard deviation, matrices not too different');
		
		
		// Test that with max noise level set to 0, the matrix is not perturbed
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "additive-noise", noiseLevelMax: 0});
		var matrixEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 1e-12 ) {
					matrixEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, true, 'Perturbed correlation matrix generation, additive noise method, 0 max noise level');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, additive noise method, 0 max noise level, all inclusive');
		

		// Test that with higher max noise level, the matrix is more perturbed
		// Can fail if the matrix had to be repaired
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "additive-noise", noiseLevelMax: 0.1});
		var matrixDifferent = false;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 0.05 ) {
					matrixDifferent = true;
					break;
				}
			}
		}	
		assert.equal(matrixDifferent, true, 'Perturbed correlation matrix generation, additive noise method, high max noise level, matrices different');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(1e-5), true, 'Perturbed correlation matrix generation, additive noise method, high max noise level, all inclusive');
  }
  
	
	// Static data
	// Test the spectral method
	{
		// Generate an input correlation matrix
		var n = 10;
		var corrMat = PortfolioAllocation.randomCorrelationMatrix(n);

		// Perturb it
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "spectral-noise"});
		
		// Test that the perturbed matrix has the same size as the original matrix
		assert.equal(perturbedCorrMat.nbColumns == perturbedCorrMat.nbColumns && perturbedCorrMat.nbRows == perturbedCorrMat.nbRows, true, 'Perturbed correlation matrix, spectral method, same matrix size');
		
		// Test that the perturbed matrix is a correlation matrix
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, spectral method, all inclusive');
		
		// Test that with default standard deviations, the matrix is slightly perturbed
		// Can fail if no luck
		var matrixEqual = true;
		var matrixNearlyEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( corrMat.getValueAt(i,j) != perturbedCorrMat.getValueAt(i,j) )  {
					matrixEqual = false;
				}

				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 0.2 ) {
					matrixNearlyEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, false, 'Perturbed correlation matrix generation, spectral method, default standard deviation, matrices different');
		assert.equal(matrixNearlyEqual, true, 'Perturbed correlation matrix generation, spectral method, default standard deviation, matrices not too different');
		
		
		// Test that with standard deviation set to 0, the matrix is not perturbed
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "spectral-noise", sigma: 0});
		var matrixEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 1e-12 ) {
					matrixEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, true, 'Perturbed correlation matrix generation, spectral method, 0 standard deviation');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, spectral method, 0 standard deviation, all inclusive');
		

		// Test that with higher standard deviation, the matrix is more perturbed
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "spectral-noise", sigma: 2});
		var matrixDifferent = false;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 0.1 ) {
					matrixDifferent = true;
					break;
				}
			}
		}	
		assert.equal(matrixDifferent, true, 'Perturbed correlation matrix generation, spectral method, high standard deviation, matrices different');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, spectral method, high standard deviation, all inclusive');
  }
  
	// Static data
	// Test the multiplicative noise method
	{
		// Generate an input correlation matrix
		var n = 30;
		var corrMat = PortfolioAllocation.randomCorrelationMatrix(n);
		
		// Perturb it
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "multiplicative-noise"});
		
		// Test that the perturbed matrix has the same size as the original matrix
		assert.equal(perturbedCorrMat.nbColumns == perturbedCorrMat.nbColumns && perturbedCorrMat.nbRows == perturbedCorrMat.nbRows, true, 'Perturbed correlation matrix, multiplicative noise method, same matrix size');
		
		// Test that the perturbed matrix is a correlation matrix
		assert.equal(perturbedCorrMat.isCorrelationMatrix(1e-5), true, 'Perturbed correlation matrix generation, multiplicative noise method, all inclusive');
		
		// Test that with default standard deviations, the matrix is slightly perturbed
		// Can fail if no luck
		var matrixEqual = true;
		var matrixNearlyEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( corrMat.getValueAt(i,j) != perturbedCorrMat.getValueAt(i,j) )  {
					matrixEqual = false;
				}

				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 0.2 ) {
					matrixNearlyEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, false, 'Perturbed correlation matrix generation, multiplicative noise method, default standard deviation, matrices different');
		assert.equal(matrixNearlyEqual, true, 'Perturbed correlation matrix generation, multiplicative noise method, default standard deviation, matrices not too different');
		
		
		// Test that with standard deviation set to 0, the matrix is not perturbed
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "multiplicative-noise", sigma: 0});
		var matrixEqual = true;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 1e-12 ) {
					matrixEqual = false;
					break;
				}
			}
		}	
		assert.equal(matrixEqual, true, 'Perturbed correlation matrix generation, multiplicative noise method, 0 standard deviation');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(), true, 'Perturbed correlation matrix generation, multiplicative noise method, 0 standard deviation, all inclusive');

		// Test that with higher standard deviation, the matrix is more perturbed
		var perturbedCorrMat = PortfolioAllocation.perturbedCorrelationMatrix(corrMat, {method: "multiplicative-noise", sigma: 0.1});
		var matrixDifferent = false;
		for (var i = 1; i <= n; ++i) {
			for (var j = 1; j <= n; ++j) {
				if ( Math.abs( corrMat.getValueAt(i,j) - perturbedCorrMat.getValueAt(i,j) ) > 1e-2 ) {
					matrixDifferent = true;
					break;
				}
			}
		}
		assert.equal(matrixDifferent, true, 'Perturbed correlation matrix generation, multiplicative noise method, high standard deviation, matrices different');
		assert.equal(perturbedCorrMat.isCorrelationMatrix(1e-5), true, 'Perturbed correlation matrix generation, multiplicative noise method, high standard deviation, all inclusive');
  }
});



QUnit.test('Random variances vector creation', function(assert) {    
	// Static data
	{
		var n = 3;
		
		// Test the size of the generated vector
		// and ensure the variances are all greater than 0
		var variances = PortfolioAllocation.randomVariances(n);
		assert.equal(variances.nbRows, n, 'Random variances vector creation, vector size');
		assert.equal(variances.isNonNegative(), true, 'Random variances vector creation - Non negative variances'); 
	  
		// Test that the standard deviation has an effect on the generated vector
		// Test with 0 standard deviation
		var variances = PortfolioAllocation.randomVariances(n, 0);
		assert.equal(variances.vectorNorm('infinity') == 0, true, 'Random variances vector creation, 0 standard deviation');

		// Test with default standard deviation (= 1)
		var variances = PortfolioAllocation.randomVariances(n);
		assert.equal(variances.vectorNorm('infinity') <= 5 && variances.vectorNorm('infinity') > 0, true, 'Random variances vector creation, default standard deviation');

		// Test with a high standard deviation
		var variances = PortfolioAllocation.randomVariances(n, 500);
		assert.equal(variances.vectorNorm('infinity') > 0  && variances.vectorNorm('infinity') >= 50, true, 'Random variances vector creation, high standard deviation');
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
		
		// With sigma equal to zero, the two vectors must be identical		
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { sigma: 0});
		var pVariancesVecArray = pVariancesVec.toArray();
		assert.deepEqual(variancesVector, pVariancesVecArray, 'Mean vector perturbations, no perturbation');
		
		// With a small sigma, the vectors must not be too different (except bad luck)
		// try 1e-3, 1e-2 with relative change !!!
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { sigma: 0.001});
		pVariancesVecArray = pVariancesVec.toArray();
		var diffOK = true;
		for (var i = 0; i < variancesVector.length; ++i) {
			if (Math.abs((pVariancesVecArray[i] - variancesVector[i])/variancesVector[i]) > 1e-2) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, small perturbation');
		
		// With a big sigma, the resulting vector must be well above 10
		pVariancesVec = PortfolioAllocation.perturbedVariances(variances, { sigma: 100000});
		pVariancesVecArray = pVariancesVec.toArray();
		var diffOK = true;
		for (var i = 0; i < variancesVector.length; ++i) {
			if (Math.abs((pVariancesVecArray[i] - variancesVector[i])/variancesVector[i]) < 10) {
				diffOK = false;
				break;
			}
		}
		assert.equal(diffOK, true, 'Mean vector perturbations, big perturbation');
	
		// Ensure the variances are all greater than 0
		assert.equal(variances.isNonNegative(), true, 'Random variances vector creation - Non negative variances');
  }

});


QUnit.test('Nearest correlation matrix', function(assert) {    
  // Test using static data
  // Reference: Computing the Nearest Correlation Matrix|A Problem from Finance, Nicholas J. Higham
	{
		var mat = [[1, 1, 0],[1, 1, 1],[0, 1, 1]];
		var n = mat.length;
		var expectedMat = [[1.0000, 0.7607, 0.1573], [0.7607, 1.0000, 0.7607], [0.1573, 0.7607, 1.0000]];
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat);
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix - Test #1');
	}
	
  // Test using static data
  // Reference: Computing the Nearest Correlation Matrix|A Problem from Finance, Nicholas J. Higham
	{
		var mat = [[2, -1, 0, 0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]];
		var n = mat.length;
		var expectedMat = [[1.0000, -0.8084, 0.1916, 0.1068], [-0.8084, 1.0000, -0.6562, 0.1916], [0.1916, -0.6562, 1.0000, -0.8084], [0.1068 ,0.1916, -0.8084, 1.000]];
		
		// Even with default precision, the generated matrix will be a numerically exact correlation matrix,
		// thanks to the congruent transformation implemented as a polishing step.
		//
		// Check this
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat);
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix - Test #2/1');
		assert.equal(nearestCorrMat.isCorrelationMatrix(), true, 'Nearest correlation matrix - Test #2, default precision');		
	}
	
  // Test using static data
  // Reference: Computing the nearest correlation matrix By Rick Wicklin on The DO Loop November 28, 2012
	{
	var mat = [[1.0,  0.99, 0.35], [0.99, 1.0,  0.80], [0.35, 0.80, 1.0]];
		var n = mat.length;
		var expectedMat = [[1.0000, 0.9044, 0.3896], [0.9044, 1.0000, 0.7453], [0.3896, 0.7453, 1.0000]];
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat, {eps:1e-8});
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix - Test #3');
	}
	
  // Test using static data
  // Reference: Computing Nearest Covariance and Correlation Matrices, Craig Lucas Thesis
  {
		var mat = [[1.0000, -0.3250, 0.1881, 0.5760, 0.0064, -0.6111, -0.0724, -0.1589],
					[-0.3250, 1.0000, 0.2048, 0.2436, 0.4058, 0.2730, 0.2869, 0.4241],
					[0.1881, 0.2048, 1.0000, -0.1325, 0.7658, 0.2765, -0.6172, 0.9006],
					[0.5760, 0.2436, -0.1325, 1.0000, 0.3041, 0.0126, 0.6452, -0.3210],
					[0.0064, 0.4058, 0.7658, 0.3041, 1.0000, 0.6652, -0.3293, 0.9939],
					[-0.6111, 0.2730, 0.2765, 0.0126, 0.6652, 1.0000, 0.0492, 0.5964],
					[-0.0724, 0.2869, -0.6172, 0.6452, -0.3293, 0.0492, 1.0000, -0.3983],
					[-0.1589, 0.4241, 0.9006, -0.3210, 0.9939, 0.5964, -0.3983, 1.0000]];
		var n = mat.length;
		var expectedMat = [[1.0000, -0.3112, 0.1889, 0.5396, 0.0268, -0.5925, -0.0621, -0.1921],
							[-0.3112, 1.0000, 0.2050, 0.2265, 0.4148, 0.2822, 0.2915, 0.4088],
							[0.1889, 0.2050, 1.0000, -0.1468, 0.7880, 0.2727, -0.6085, 0.8802],
							[0.5396, 0.2265, -0.1468, 1.0000, 0.2137, 0.0015, 0.6069, -0.2208],
							[0.0268, 0.4148, 0.7880, 0.2137, 1.0000, 0.6580, -0.2812, 0.8762],
							[-0.5925, 0.2822, 0.2727, 0.0015, 0.6580, 1.0000, 0.0479, 0.5932],
							[-0.0621, 0.2915, -0.6085, 0.6069, -0.2812, 0.0479, 1.0000, -0.4470],
							[-0.1921, 0.4088, 0.8802, -0.2208, 0.8762, 0.5932, -0.4470, 1.0000]];
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat);

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix - Test #4');
	}
	
  // Test using static data
  // Reference: Anderson acceleration of the alternating projections method for computing the nearest correlation matrix, Nicholas J. Higham, Natasa Strabic
  // Calculations below validated using nearcorr_new.m Matlab program from the authors with maximum accuracy
  {
		var mat = [[1,     0.18,  -0.13, -0.26, 0.19,  -0.25, -0.12],
					 [0.18,  1,     0.22,  -0.14, 0.31,  0.16,  0.09],
					 [-0.13, 0.22,  1,     0.06,  -0.08, 0.04,  0.04],
					 [-0.26, -0.14, 0.06,  1,     0.85,  0.85,  0.85],
					 [0.19,  0.31,  -0.08, 0.85,  1,     0.85,  0.85],
					 [-0.25, 0.16,  0.04,  0.85,  0.85,  1,     0.85],
					 [-0.12, 0.09,  0.04,  0.85,  0.85,  0.85,  1]];
		var n = mat.length;
		
		// When no minimum eigenvalue is provided, 1e-8 must be the default
		var expectedMat =  [[1.000000, 0.183844, -0.131789, -0.251417, 0.178399, -0.247928, -0.119061],
							 [ 0.183844, 1.000000, 0.218241, -0.131564, 0.298599, 0.162037, 0.090923],
							 [-0.131789, 0.218241, 1.000000, 0.056073, -0.074692, 0.039052, 0.039570],
							 [-0.251417, -0.131564, 0.056073, 1.000000, 0.824539, 0.854548, 0.852062],
							 [ 0.178399, 0.298599, -0.074692, 0.824539, 1.000000, 0.843853, 0.847214],
							 [-0.247928, 0.162037, 0.039052, 0.854548, 0.843853, 1.000000, 0.850498],
							 [-0.119061, 0.090923, 0.039570, 0.852062, 0.847214, 0.850498, 1.000000]];
		var nearestCorrMatDefault = PortfolioAllocation.nearestCorrelationMatrix(mat);
		var nearestCorrMatDelta = PortfolioAllocation.nearestCorrelationMatrix(mat, {minEigenvalue: 1e-8});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMatDefault.getValueAt(i+1,j+1) ) > 1e-4 || 
				     Math.abs( nearestCorrMatDefault.getValueAt(i+1,j+1) - nearestCorrMatDelta.getValueAt(i+1,j+1) ) > 0 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix, positive definite - Delta = 0');
		
		// Test with a minimum eigenvalue of 1e-8, the matrix must not change
		var expectedMat =  [[1.000000, 0.183844, -0.131789, -0.251417, 0.178399, -0.247928, -0.119061],
							 [ 0.183844, 1.000000, 0.218241, -0.131564, 0.298599, 0.162037, 0.090923],
							 [-0.131789, 0.218241, 1.000000, 0.056073, -0.074692, 0.039052, 0.039570],
							 [-0.251417, -0.131564, 0.056073, 1.000000, 0.824539, 0.854548, 0.852062],
							 [ 0.178399, 0.298599, -0.074692, 0.824539, 1.000000, 0.843853, 0.847214],
							 [-0.247928, 0.162037, 0.039052, 0.854548, 0.843853, 1.000000, 0.850498],
							 [-0.119061, 0.090923, 0.039570, 0.852062, 0.847214, 0.850498, 1.000000]];
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat, {minEigenvalue: 1e-8});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix, positive definite - Delta = 1e-8');
		
		// Test with a minimum eigenvalue of 0.1
		var expectedMat =  [[1.000000, 0.192428, -0.135963, -0.232656, 0.143938, -0.237614, -0.114311],
							[0.192428, 1.000000, 0.213837, -0.106498, 0.269406, 0.163030, 0.093218],
							[-0.135963, 0.213837, 1.000000, 0.045597, -0.060853, 0.037398, 0.038215],
							[-0.232656, -0.106498, 0.045597, 1.000000, 0.759191, 0.851434, 0.855982],
							[0.143938, 0.269406, -0.060853, 0.759191, 1.000000, 0.815831, 0.833260],
							[-0.237614, 0.163030, 0.037398, 0.851434, 0.815831, 1.000000, 0.857818],
							[-0.114311, 0.093218, 0.038215, 0.855982, 0.833260, 0.857818, 1.000000]];
		var nearestCorrMat = PortfolioAllocation.nearestCorrelationMatrix(mat, {minEigenvalue: 0.1});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedMat[i][j] - nearestCorrMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Nearest correlation matrix, positive definite - Delta = 0.1');
	}
	
	// Test using static data on a matrix which takes ~20000 iterations to converge to full precision
	{
		var mat = [[                      1,   0.03362837328256285,    0.1889627792864824,    0.2057437574498187,-0.0044618267644031704,    0.0975573943197419,  -0.12560468451285348,   0.17445948346484055,   0.08091522232237142,  -0.38645611295971993,  -0.07744236763031459,   0.10544055441000916,   -0.3304416390282105,    0.5786386214395042,   0.08154866142160715,    0.3111696478005631,  -0.06500097193573268,   -0.8410632620400488, -0.012748426648486667,   0.31833088826251565,    0.5495462420516737,  -0.12136514741700843,   -0.4819462208558094,   -0.3168422503925507,    0.0694765260708388,     -0.84383330246617,  0.033925645085928255,   -0.2369898943772112,  0.015032758726772899,   0.08710227399541885],
					[    0.03362837328256285,                     1,  -0.06964962319293763, -0.013595236591766162,  -0.40299134364956146,  -0.04110639493842262,   0.18151727257372285,    0.3976721834684879,   0.17510685313091154,  -0.05908194180910702, -0.028255303569252162,   -0.2562631940880807,    0.5208356036252447, -0.049393850575032705,  0.003120878037483506, -0.042496508355601044,   0.03308056704461644,   -0.0634236992392466,   0.05823509287341512, -0.004725844949174946,   -0.5636194966526595,  -0.04073761682678527,   0.21944576481354322, -0.010350813662266659,    0.3634716254189384,  0.013683665383489514,    0.0666749578383111,   0.11904196598439208,   -0.2654284061714957,    0.4730379093934179],
					[     0.1889627792864824,  -0.06964962319293763,                     1,   -0.3526726700003423,   0.04978419830658411,   0.06392330839310646,  -0.12529924087027172,    0.0442683475771076,   -0.1954839835539136,    0.1817869724525281,   0.39739472154698086,   -0.1303473566860232,    0.6848085954190859,   0.09760724971562995,   -0.7120705331979182,  -0.08545324322060366,   0.19719796172177154,   -0.8294571303359755,     0.697269253708484,  0.025404700858405686,   0.06630194899010945,   0.31397639049163356,   0.36753878221352404,    0.1525595838014024,   0.26600778408203346, -0.016148718125606948,   0.20064969655798662,   0.30039482785073884,   -0.2777363372965186,   0.35617772156815514],
					[     0.2057437574498187, -0.013595236591766162,   -0.3526726700003423,                     1,0.00006152655390854771,   0.06268871139170942,  0.053466682769447375,   0.05215592426598065,   -0.3709258916036329, -0.050762684592276244,-0.0028927105990979468,  -0.13265214387080135,  -0.10803595064969912,  -0.05668674305399061,  -0.02654984266654007,   0.01665971801211744,  -0.05167651977381764,  0.034349543467026844, -0.046200784502070415, -0.011258789399801884,  -0.19874772326733112,  0.006536490324764677,  -0.01728558434073849,    0.5940689330822759,  -0.04677991500894784,   0.02428865859473822,    -0.172671693296398,   -0.2392659445949045,  -0.23754577245707648,    -0.275762559323563],
					[ -0.0044618267644031704,  -0.40299134364956146,   0.04978419830658411,0.00006152655390854771,                     1,  -0.36123602402080895,  -0.03268485976827452,    0.2943780065184315,   0.19619025214544394,   -0.2552799035670831,   -0.7081591732927086,   0.15859330685027603, -0.025570439744845994,   0.04299217738141946,   0.12743815296251407,   0.15099402367169068,   -0.4976329138583712,  -0.13034754752618524,    -0.384489092274836,-0.0018288529508970648,  -0.16532027189344853, -0.043864205765321074,   -0.2487918490903405, 0.0012027977885858332,  0.022494166963966222,   -0.2869172452093622,  -0.32816731378686254,   -0.0384365047539285,  0.028185866289586504,  -0.28342265951749956],
					[     0.0975573943197419,  -0.04110639493842262,   0.06392330839310646,   0.06268871139170942,  -0.36123602402080895,                     1,  -0.19450985849624064,  -0.29761632819999506,    0.5469116889875193,   0.36538541010299114,  -0.11737588161349842,   -0.4814314609019866,   0.19172018239447977,  -0.20769699227242583,     0.610353313985166,    0.1529056249512108,  -0.05624930193791267,   0.22201827499702917,   0.30450367624673386,  -0.10000484624014182, -0.022914355240576976,   0.18768640522719002, -0.023209166871707706,   0.10152110743360489,  0.009283593076093243,    0.0752398165886433,   0.09862759498249037,  -0.21921735455072244,   -0.2136808201707631,   0.07462720947319375],
					[   -0.12560468451285348,   0.18151727257372285,  -0.12529924087027172,  0.053466682769447375,  -0.03268485976827452,  -0.19450985849624064,                     1,  0.000342737205531862, 0.0013339493197323676, -0.009583962538139106,   0.05873106369639609,  -0.01638202084940745,    0.1990075947604494,   -0.4063169912987429,  0.015105051647033239,  0.047552813957069145,   0.23055411775553347,    0.0329377888086339,   -0.3959834400280367,  -0.04521783406957059,   0.09154183073843764,  -0.08866383842786077,    0.3679536689748993,    0.5115212527586548,   0.03333047250170801,   0.17740308072588865,  -0.04652168452863328,  -0.19779051550087912,      0.48416358823977,   0.12380359058118116],
					[    0.17445948346484055,    0.3976721834684879,    0.0442683475771076,   0.05215592426598065,    0.2943780065184315,  -0.29761632819999506,  0.000342737205531862,                     1, -0.007201659094295932,   -0.3932700929329244,   -0.6164675517139614,  0.018471399305595316,   0.03475627652619073,   0.13703971574926513,  -0.29792942264432887,   0.21418470856733743, 0.0029489682188216892,    0.1066046844009827,  -0.14371420260483023,  0.011949701609853437,  -0.17815350178184353,   0.37559730564678273,  -0.20216879264319224,   -0.1312944090861684,   0.17340981313031542,  0.007872315510893325,   0.04847497922596598,   0.07512117575433298,    0.1909367151944158,    0.0323811242415796],
					[    0.08091522232237142,   0.17510685313091154,   -0.1954839835539136,   -0.3709258916036329,   0.19619025214544394,    0.5469116889875193, 0.0013339493197323676, -0.007201659094295932,                     1,  0.000771188145206407,  -0.19833816161521087,  -0.01959478024967369,   0.14237848206818346,  -0.03243933484038869,  -0.03194339346831293,   0.12978345879504077,   0.01626595934518267,   -0.3111113831084434, -0.009907988622202468,  -0.05668875181855866,   0.35922551090016835,  -0.12073772388670044,  -0.34298538556866776,  -0.23433729443231405,    0.2687949995026267,   -0.6875192705977949,  0.022488278013989352,    0.2667287225373032,  -0.14950035330945455,    0.1301184202252945],
					[   -0.38645611295971993,  -0.05908194180910702,    0.1817869724525281, -0.050762684592276244,   -0.2552799035670831,   0.36538541010299114, -0.009583962538139106,   -0.3932700929329244,  0.000771188145206407,                     1,  -0.07329766854273546,   0.19923620526264568,   0.05525038269469989, -0.016697733533609524,  -0.18181052283018823,  -0.01178410024202195,   0.26582414355924927,  -0.24216585304794294,    0.8154049779258034,  -0.15552210894479584,0.00036606808670074716, 0.0001038187111231005,  -0.16361610845226018,    0.3699741358824692,  -0.09071585462254621,   0.30121030510421154,  0.018755549124796957,  -0.31029522474521315,   -0.3657930466497671,  -0.20207413829413484],
					[   -0.07744236763031459, -0.028255303569252162,   0.39739472154698086,-0.0028927105990979468,   -0.7081591732927086,  -0.11737588161349842,   0.05873106369639609,   -0.6164675517139614,  -0.19833816161521087,  -0.07329766854273546,                     1,  -0.21753839731245536,  -0.03711194928741433, -0.023349301059714327,  -0.08888929221022504,   -0.8372417918532487,   0.20339086663796976, -0.006179144319545693,    0.8473815776217516,   0.03175732389893964,  0.019402162126520305, -0.007535561323084914,   0.37491502067553284,  -0.06784171500570993, -0.006151934521117095,    0.2915956745947752,   0.23934311092924576,   -0.3517456301017262,  -0.15479667043964349,   0.05214490025342785],
					[    0.10544055441000916,   -0.2562631940880807,   -0.1303473566860232,  -0.13265214387080135,   0.15859330685027603,   -0.4814314609019866,  -0.01638202084940745,  0.018471399305595316,  -0.01959478024967369,   0.19923620526264568,  -0.21753839731245536,                     1,  0.025676171463345755,  -0.24353830707512056,   0.20457220100192736,   0.09181965401786074,  -0.11774442785177347,   0.14854229919978537,  -0.14419619842294013, 0.0032016216459880282,  -0.04786288893371655,  0.017641162849298453, -0.054121304225692037,  -0.13672385232046186, -0.010958019690700368,   0.06013944824125146,  0.021605142818977185,  0.007674282582095672,    0.7351440677796907,  0.006559576479502226],
					[    -0.3304416390282105,    0.5208356036252447,    0.6848085954190859,  -0.10803595064969912, -0.025570439744845994,   0.19172018239447977,    0.1990075947604494,   0.03475627652619073,   0.14237848206818346,   0.05525038269469989,  -0.03711194928741433,  0.025676171463345755,                     1,  -0.13928832129079055,  -0.15288436707677142,  -0.21558045657443575,   0.09162937523595518, -0.006360647568300344,   0.44275596793948885,   -0.0227190733948297,  0.022868859829516825,     0.301594067038583,  -0.21923048817615526,  -0.09290014037834599,     0.334083962333145,    0.9370922102744562,  -0.06666823330426179,   -0.5954602539509863,  -0.21880328909120453,   0.06500213262810271],
					[     0.5786386214395042, -0.049393850575032705,   0.09760724971562995,  -0.05668674305399061,   0.04299217738141946,  -0.20769699227242583,   -0.4063169912987429,   0.13703971574926513,  -0.03243933484038869, -0.016697733533609524, -0.023349301059714327,  -0.24353830707512056,  -0.13928832129079055,                     1,   0.13720609792912936,  -0.01835996342898162,  0.053517476625399414,  -0.07952769922799299,  -0.30325742710241377,   0.46796596439404897,    0.6264842079934234,  -0.01686693234527459,    -0.899255527324118,  -0.45153750781238866,   0.05782049578670538,   -0.3620234652916817,   0.16079996159378404,    0.2808445006715864,    0.0106509924685626,     0.337778017544219],
					[    0.08154866142160715,  0.003120878037483506,   -0.7120705331979182,  -0.02654984266654007,   0.12743815296251407,     0.610353313985166,  0.015105051647033239,  -0.29792942264432887,  -0.03194339346831293,  -0.18181052283018823,  -0.08888929221022504,   0.20457220100192736,  -0.15288436707677142,   0.13720609792912936,                     1,   0.38041000024669963,    0.0578308533736533,    0.3774346318581759, -0.047687515977000286,   -0.0867880983842151,  -0.27330012865028797,   -0.7047670599702823,   -0.1976250775955519,    0.7372350741465112,  -0.11077163584879184,  -0.04186470867001903,0.00040586576968277253, -0.020357386478963527,  0.036245329574606584,    0.0511825390367075],
					[     0.3111696478005631, -0.042496508355601044,  -0.08545324322060366,   0.01665971801211744,   0.15099402367169068,    0.1529056249512108,  0.047552813957069145,   0.21418470856733743,   0.12978345879504077,  -0.01178410024202195,   -0.8372417918532487,   0.09181965401786074,  -0.21558045657443575,  -0.01835996342898162,   0.38041000024669963,                     1,   -0.3797432205210555,   0.19275216598421915, -0.046575519557213015,  -0.10903490586431608,   0.25464449480970086,  -0.01940919422141176,   -0.1200206989070653,   0.11508587856606452, -0.030686125998387284,  0.020846349165428796,  0.050992843774932135,  -0.13538976954715845,    0.2044675189510433,   -0.4376182613065415],
					[   -0.06500097193573268,   0.03308056704461644,   0.19719796172177154,  -0.05167651977381764,   -0.4976329138583712,  -0.05624930193791267,   0.23055411775553347, 0.0029489682188216892,   0.01626595934518267,   0.26582414355924927,   0.20339086663796976,  -0.11774442785177347,   0.09162937523595518,  0.053517476625399414,    0.0578308533736533,   -0.3797432205210555,                     1,  -0.10896910152698858,  0.007786975343206232, -0.014597311019823296,  0.001665451710272586,   0.04788204634775381,    0.3641163005561989,   -0.2106273326307495,   0.37208001647324745,    0.1989916087085751,  -0.13653601413685823,  -0.04157874776157013,   0.14030243108762938, -0.002452335191418184],
					[    -0.8410632620400488,   -0.0634236992392466,   -0.8294571303359755,  0.034349543467026844,  -0.13034754752618524,   0.22201827499702917,    0.0329377888086339,    0.1066046844009827,   -0.3111113831084434,  -0.24216585304794294, -0.006179144319545693,   0.14854229919978537, -0.006360647568300344,  -0.07952769922799299,    0.3774346318581759,   0.19275216598421915,  -0.10896910152698858,                     1,  -0.31605109653305374,  -0.11794916086870087,  -0.06983744897569738,  -0.06311513430980038,-0.0016481403638468287,    0.1344942213117672,  0.017187184385583115,  -0.03709822059681859,   0.14170606365288696,   -0.0411706193108413,    0.6055357951557646,  -0.13616091604900568],
					[  -0.012748426648486667,   0.05823509287341512,     0.697269253708484, -0.046200784502070415,    -0.384489092274836,   0.30450367624673386,   -0.3959834400280367,  -0.14371420260483023, -0.009907988622202468,    0.8154049779258034,    0.8473815776217516,  -0.14419619842294013,   0.44275596793948885,  -0.30325742710241377, -0.047687515977000286, -0.046575519557213015,  0.007786975343206232,  -0.31605109653305374,                     1,  -0.03309757203740167,   0.06637886538682056,    0.4453931297759861,   0.09464198682174094, -0.051582407442756516,  0.043680115818458765,    0.2501976320896041,   0.21393108339367176,  -0.21893243582913247,   -0.5621474670890907,  -0.15351885400011267],
					[    0.31833088826251565, -0.004725844949174946,  0.025404700858405686, -0.011258789399801884,-0.0018288529508970648,  -0.10000484624014182,  -0.04521783406957059,  0.011949701609853437,  -0.05668875181855866,  -0.15552210894479584,   0.03175732389893964, 0.0032016216459880282,   -0.0227190733948297,   0.46796596439404897,   -0.0867880983842151,  -0.10903490586431608, -0.014597311019823296,  -0.11794916086870087,  -0.03309757203740167,                     1,-0.0006262992420206521,  -0.48521657148899827,   0.03130326715821451,    -0.077210590017826,  -0.09327658302394438,  -0.02408218064317969,  -0.06625520391842574,-0.0060813421850797736,   0.06189745043264944,   0.06933578750874124],
					[     0.5495462420516737,   -0.5636194966526595,   0.06630194899010945,  -0.19874772326733112,  -0.16532027189344853, -0.022914355240576976,   0.09154183073843764,  -0.17815350178184353,   0.35922551090016835,0.00036606808670074716,  0.019402162126520305,  -0.04786288893371655,  0.022868859829516825,    0.6264842079934234,  -0.27330012865028797,   0.25464449480970086,  0.001665451710272586,  -0.06983744897569738,   0.06637886538682056,-0.0006262992420206521,                     1,  -0.21901765708235735,   -0.4432611641882161,  -0.10338618253708498,   -0.1660293790719887,   -0.5998520750577832,  -0.16769199903981316,  -0.06517554206559666,   0.04945026233886099,   0.06942588225083968],
					[   -0.12136514741700843,  -0.04073761682678527,   0.31397639049163356,  0.006536490324764677, -0.043864205765321074,   0.18768640522719002,  -0.08866383842786077,   0.37559730564678273,  -0.12073772388670044, 0.0001038187111231005, -0.007535561323084914,  0.017641162849298453,     0.301594067038583,  -0.01686693234527459,   -0.7047670599702823,  -0.01940919422141176,   0.04788204634775381,  -0.06311513430980038,    0.4453931297759861,  -0.48521657148899827,  -0.21901765708235735,                     1,   0.43518098878498573,    0.1363536299121965,   0.36201373100107886,    0.5730709805002415,   0.17358316153769748,    0.1073037962740902,  -0.04661921312641353,    0.2660736899939656],
					[    -0.4819462208558094,   0.21944576481354322,   0.36753878221352404,  -0.01728558434073849,   -0.2487918490903405, -0.023209166871707706,    0.3679536689748993,  -0.20216879264319224,  -0.34298538556866776,  -0.16361610845226018,   0.37491502067553284, -0.054121304225692037,  -0.21923048817615526,    -0.899255527324118,   -0.1976250775955519,   -0.1200206989070653,    0.3641163005561989,-0.0016481403638468287,   0.09464198682174094,   0.03130326715821451,   -0.4432611641882161,   0.43518098878498573,                     1,    0.2371052596554239,   0.18883416248534676,    0.8735710431504162,   0.40752509091722366,  -0.06760117991476561,   -0.0861808614511945,   0.23250158080841152],
					[    -0.3168422503925507, -0.010350813662266659,    0.1525595838014024,    0.5940689330822759, 0.0012027977885858332,   0.10152110743360489,    0.5115212527586548,   -0.1312944090861684,  -0.23433729443231405,    0.3699741358824692,  -0.06784171500570993,  -0.13672385232046186,  -0.09290014037834599,  -0.45153750781238866,    0.7372350741465112,   0.11508587856606452,   -0.2106273326307495,    0.1344942213117672, -0.051582407442756516,    -0.077210590017826,  -0.10338618253708498,    0.1363536299121965,    0.2371052596554239,                     1,    0.3484413193935787,   0.09739987200907463,  -0.21979905692213114,    0.3607693551943086,   0.07439040078283819,   -0.3132710378563707],
					[     0.0694765260708388,    0.3634716254189384,   0.26600778408203346,  -0.04677991500894784,  0.022494166963966222,  0.009283593076093243,   0.03333047250170801,   0.17340981313031542,    0.2687949995026267,  -0.09071585462254621, -0.006151934521117095, -0.010958019690700368,     0.334083962333145,   0.05782049578670538,  -0.11077163584879184, -0.030686125998387284,   0.37208001647324745,  0.017187184385583115,  0.043680115818458765,  -0.09327658302394438,   -0.1660293790719887,   0.36201373100107886,   0.18883416248534676,    0.3484413193935787,                     1,   0.16340469107068892,  -0.04516429283760329,   0.03591342436268478, -0.057567882682352314,   0.12930985038319273],
					[      -0.84383330246617,  0.013683665383489514, -0.016148718125606948,   0.02428865859473822,   -0.2869172452093622,    0.0752398165886433,   0.17740308072588865,  0.007872315510893325,   -0.6875192705977949,   0.30121030510421154,    0.2915956745947752,   0.06013944824125146,    0.9370922102744562,   -0.3620234652916817,  -0.04186470867001903,  0.020846349165428796,    0.1989916087085751,  -0.03709822059681859,    0.2501976320896041,  -0.02408218064317969,   -0.5998520750577832,    0.5730709805002415,    0.8735710431504162,   0.09739987200907463,   0.16340469107068892,                     1,   0.17028166009224072,   -0.9205685482449166,   0.11821883860981837,   0.03316783689365472],
					[   0.033925645085928255,    0.0666749578383111,   0.20064969655798662,    -0.172671693296398,  -0.32816731378686254,   0.09862759498249037,  -0.04652168452863328,   0.04847497922596598,  0.022488278013989352,  0.018755549124796957,   0.23934311092924576,  0.021605142818977185,  -0.06666823330426179,   0.16079996159378404,0.00040586576968277253,  0.050992843774932135,  -0.13653601413685823,   0.14170606365288696,   0.21393108339367176,  -0.06625520391842574,  -0.16769199903981316,   0.17358316153769748,   0.40752509091722366,  -0.21979905692213114,  -0.04516429283760329,   0.17028166009224072,                     1,  -0.18792754584037022,   -0.3268740011378931,  0.045110047369108744],
					[    -0.2369898943772112,   0.11904196598439208,   0.30039482785073884,   -0.2392659445949045,   -0.0384365047539285,  -0.21921735455072244,  -0.19779051550087912,   0.07512117575433298,    0.2667287225373032,  -0.31029522474521315,   -0.3517456301017262,  0.007674282582095672,   -0.5954602539509863,    0.2808445006715864, -0.020357386478963527,  -0.13538976954715845,  -0.04157874776157013,   -0.0411706193108413,  -0.21893243582913247,-0.0060813421850797736,  -0.06517554206559666,    0.1073037962740902,  -0.06760117991476561,    0.3607693551943086,   0.03591342436268478,   -0.9205685482449166,  -0.18792754584037022,                     1,    0.2537646064466628,  -0.12446104070927157],
					[   0.015032758726772899,   -0.2654284061714957,   -0.2777363372965186,  -0.23754577245707648,  0.028185866289586504,   -0.2136808201707631,      0.48416358823977,    0.1909367151944158,  -0.14950035330945455,   -0.3657930466497671,  -0.15479667043964349,    0.7351440677796907,  -0.21880328909120453,    0.0106509924685626,  0.036245329574606584,    0.2044675189510433,   0.14030243108762938,    0.6055357951557646,   -0.5621474670890907,   0.06189745043264944,   0.04945026233886099,  -0.04661921312641353,   -0.0861808614511945,   0.07439040078283819, -0.057567882682352314,   0.11821883860981837,   -0.3268740011378931,    0.2537646064466628,                     1,  -0.13226489541600775],
					[    0.08710227399541885,    0.4730379093934179,   0.35617772156815514,    -0.275762559323563,  -0.28342265951749956,   0.07462720947319375,   0.12380359058118116,    0.0323811242415796,    0.1301184202252945,  -0.20207413829413484,   0.05214490025342785,  0.006559576479502226,   0.06500213262810271,     0.337778017544219,    0.0511825390367075,   -0.4376182613065415, -0.002452335191418184,  -0.13616091604900568,  -0.15351885400011267,   0.06933578750874124,   0.06942588225083968,    0.2660736899939656,   0.23250158080841152,   -0.3132710378563707,   0.12930985038319273,   0.03316783689365472,  0.045110047369108744,  -0.12446104070927157,  -0.13226489541600775,                     1]]
		  
		  var nearest_mat = PortfolioAllocation.nearestCorrelationMatrix(mat, {eps: 0});
		  assert.equal(nearest_mat.isCorrelationMatrix(), true, 'Nearest correlation matrix - Test high number of iterations to full precision');	
	}

});

QUnit.test('Covariance matrix from correlation matrix', function(assert) {    
	// Test error cases
	{  
		assert.throws(function() { 
			var stddev = [1, 2];
			var corrMat = [[1, 0.1, 3], [0.1, 1, 2]];
			var covMat = PortfolioAllocation.covarianceMatrixFromCorrelationMatrix(corrMat, stddev, {diagonalVectorType: "standard-deviations"});
		},
		 new Error('input correlation matrix dimensions not compatible with input diagonal vector dimension'),
		 "Covariance matrix creation from correlation matrix, correlation matrix not square");
		 
		assert.throws(function() { 
			var stddev = [1, 2, 2];
			var corrMat = [[1, 0.1], [0.1, 1]];
			var covMat = PortfolioAllocation.covarianceMatrixFromCorrelationMatrix(corrMat, stddev, {diagonalVectorType: "standard-deviations"});
		},
		 new Error('input correlation matrix dimensions not compatible with input diagonal vector dimension'),
		 "Covariance matrix creation from correlation matrix, incompatible dimensions matrix - vector");
	}
  
  // Test using static data
  // Reference: Table 9 from the paper Constrained Risk Budgeting Portfolios https://arxiv.org/abs/1902.05710
	{
		var stddev = [0.05,0.05,0.07,0.1,0.15,0.15,0.15,0.18];
		var corrMat = [[100,  80,  60, -20, -10, -20, -20, -20],
						[ 80, 100,  40, -20, -20, -10, -20, -20],
						[ 60,  40, 100,  50,  30,  20,  20,  30],
						[-20, -20,  50, 100,  60,  60,  50,  60],
						[-10, -20,  30,  60, 100,  90,  70,  70],
						[-20, -10,  20,  60,  90, 100,  60,  70],
						[-20, -20,  20,  50,  70,  60, 100,  70],
						[-20, -20,  30,  60,  70,  70,  70, 100]];
		var n = corrMat.length;
		
		// Normalize the correlation matrix
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				corrMat[i][j] = corrMat[i][j] / 100;
			}
		}
		
		var expectedCovMat = [[ 0.0025 ,  0.002  ,  0.0021 , -0.001  , -0.00075, -0.0015 , -0.0015 , -0.0018 ],
							  [ 0.002  ,  0.0025 ,  0.0014 , -0.001  , -0.0015 , -0.00075, -0.0015 , -0.0018 ],
							  [ 0.0021 ,  0.0014 ,  0.0049 ,  0.0035 ,  0.00315,  0.0021 , 0.0021 ,  0.00378],
							  [-0.001  , -0.001  ,  0.0035 ,  0.01   ,  0.009  ,  0.009  , 0.0075 ,  0.0108 ],
							  [-0.00075, -0.0015 ,  0.00315,  0.009  ,  0.0225 ,  0.02025, 0.01575,  0.0189 ],
							  [-0.0015 , -0.00075,  0.0021 ,  0.009  ,  0.02025,  0.0225 , 0.0135 ,  0.0189 ],
							  [-0.0015 , -0.0015 ,  0.0021 ,  0.0075 ,  0.01575,  0.0135 , 0.0225 ,  0.0189 ],
							  [-0.0018 , -0.0018 ,  0.00378,  0.0108 ,  0.0189 ,  0.0189 , 0.0189 ,  0.0324 ]];
		var covMat = PortfolioAllocation.covarianceMatrixFromCorrelationMatrix(corrMat, stddev, {diagonalVectorType: "standard-deviations"});

		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedCovMat[i][j] - covMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Covariance matrix computation from correlation matrix, test #1');
	}
	

	// Test using static data
	// Reference: Chopra, Vijay K; Ziemba, William T, The effect of errors in means, variances, and covariances on optimal portfolio choice
	{
		var std = [8.8308, 8.4585, 10.040, 8.6215, 5.9886, 6.8767, 5.8162, 5.6385, 8.0047, 8.2125];
		var corr = [[1,      0.3660, 0.3457, 0.1606, 0.2279, 0.5133, 0.5203, 0.2176, 0.3267, 0.5101],
				   [0.3660, 1,      0.5379, 0.2165, 0.4986, 0.5823, 0.5569, 0.4760, 0.6517, 0.5853],
				   [0.3457, 0.5379, 1,      0.2218, 0.4283, 0.4051, 0.4492, 0.3867, 0.4883, 0.6569],
				   [0.1606, 0.2165, 0.2218, 1,      0.0569, 0.3609, 0.2325, 0.2289, 0.1726, 0.3814],
				   [0.2279, 0.4986, 0.4283, 0.0569, 1,      0.3619, 0.4811, 0.5952, 0.4378, 0.4368],
				   [0.5133, 0.5823, 0.4051, 0.3609, 0.3619, 1,      0.6167, 0.4996, 0.5811, 0.5644],
				   [0.5203, 0.5569, 0.4492, 0.2325, 0.4811, 0.6167, 1,      0.6037, 0.5671, 0.6032],
				   [0.2176, 0.4760, 0.3867, 0.2289, 0.5952, 0.4996, 0.6037, 1,      0.5012, 0.4772],
				   [0.3267, 0.6517, 0.4883, 0.1726, 0.4378, 0.5811, 0.5671, 0.5012, 1,      0.6039],
				   [0.5101, 0.5853, 0.6569, 0.3814, 0.4368, 0.5644, 0.6032, 0.4772, 0.6039, 1]];
		var n = corr.length;
		
		// Compute the variances
		var variances = new Array(n);
		for (var i = 0; i < n; ++i) {
			variances[i] = std[i] * std[i];
		}
		
		var expectedCovMat = [[77.98302864,      27.3384877788 , 30.650187902399995 ,12.227239597319999  ,  12.052292971752, 31.171047119387996  ,  26.723491968888 ,10.834840558079998, 23.093738485091997  ,    36.9939542445],
							   [27.338487778800005,  71.54622225000001,       45.680264586, 15.788253352875001  ,   25.25637014766 , 33.87039193498501  ,   27.39743489613    ,   22.701988071  ,  44.125143900915, 40.658116910625004],
							   [30.6501879024,       45.680264586, 100.80159999999998    ,   19.198976948   ,   25.7517704952, 27.968941746799995  ,    26.2308758816, 21.891295817999996, 39.243297900399995    ,    54.16370415],
							   [12.227239597319999 ,   15.788253352875, 19.198976947999995 , 74.33026224999999 ,     2.93778767781, 21.396847580144996, 11.658565629749999, 11.127361821974999, 11.911561133229998 ,    27.00467182125],
							   [12.052292971752001 ,    25.25637014766,      25.7517704952 ,     2.93778767781  ,      35.86332996 ,   14.903695453878 ,   16.757143738452, 20.097952398719997 ,20.986795142675998, 21.482425692000003],
							   [31.171047119388  ,  33.870391934985, 27.968941746799995, 21.396847580144996, 14.903695453877999 , 47.28900288999999  ,  24.665695108418 ,19.371626765819997, 31.987184396738993, 31.874432854499997],
							   [26.723491968888  ,   27.39743489613  ,    26.2308758816, 11.658565629749999  ,  16.757143738452  ,  24.665695108418, 33.828182440000006 ,    19.79812640169 ,26.402438484994004   ,    28.812175236],
							   [10.834840558079998 ,22.701988071000002, 21.891295817999996 ,11.127361821974999, 20.097952398719997, 19.371626765819997 ,    19.79812640169, 31.792682249999995, 22.621411876139998 ,     22.0973096925],
							   [23.093738485091997 ,   44.125143900915 ,39.243297900399995 ,11.911561133229998 ,   20.986795142676, 31.987184396738993, 26.402438484994004, 22.621411876139998  ,      64.07522209, 39.699539785125005],
							   [36.9939542445 ,40.658116910625004  ,      54.16370415  ,   27.00467182125    ,   21.482425692   ,   31.8744328545   ,    28.812175236  ,    22.0973096925  ,  39.699539785125 , 67.44515625000001]];
		var covMat = PortfolioAllocation.covarianceMatrixFromCorrelationMatrix(corr, variances);
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( expectedCovMat[i][j] - covMat.getValueAt(i+1,j+1) ) > 1e-4 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Covariance matrix computation from correlation matrix, test #2');		   
	}
});

QUnit.test('Covariance matrix from covariance matrix', function(assert) {    
	// Test error cases
	{  
		assert.throws(function() { 
			var covMat = [[1, 0.1], [0.1, 1, 2]];
			var covMat = PortfolioAllocation.covarianceMatrixFromCovarianceMatrix(covMat);
		},
		 new Error('unsupported input type'),
		 "Covariance matrix creation from covariance matrix, covariance matrix dimensions not consistent");
	}
  
  // Test using static data
	{
		var n = 2;
		var covMat = [[1, 0.1], [0.1, 1]];
		var covMatMatrix = PortfolioAllocation.covarianceMatrixFromCovarianceMatrix(covMat);
		
		var matrixEqual = true;
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < n; ++j) {
				if ( Math.abs( covMat[i][j] - covMatMatrix.getValueAt(i+1,j+1) ) > 0 ) {
					matrixEqual = false;
					break;
				}
			}
		}
		assert.equal(matrixEqual, true, 'Covariance matrix computation from covariance matrix, test #1');
	}
});

QUnit.test('Covariance matrix sgv', function(assert) {    
  // Test using random data
	{
		var mat = PortfolioAllocation.randomCorrelationMatrix(4).toCovarianceMatrix();
		assert.equal(mat.standardizedGeneralizedVariance(), Math.pow(mat.determinant(), 1/4), 'Sgv computation');
	}
});