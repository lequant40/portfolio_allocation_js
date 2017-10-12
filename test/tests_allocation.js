// ------------------------------------------------------------
QUnit.module('Assets allocation module', {
  before: function() {
    // Generate a random 2*2 covariance matrix 
	this.randomCovarianceMatrix = function() {
		// Generate a random 2*2 covariance matrix
		var var1 = 1 - Math.random(); // belongs to ]0,1]
		var var2 = 1 - Math.random(); // belongs to ]0,1]
		var cov = 2*Math.random() - 1; // belongs to ]-1,1[
		while (var1*var2 - cov*cov <= 0) { // The covariance matrix must be DP, so that the covariance cannot be any value: determinant must be strictly positive
			cov = 2*Math.random() - 1;
		}
		
		return [[var1, cov], [cov, var2]];
	};
  }
});


QUnit.test('Equal weights portfolio', function(assert) {    
  // Random data, n assets case
  {
    // Generate a random number of assets
	var nbAssets = Math.floor(Math.random()*(50-1+1) + 1); // max 50 min 1
    
	// Compute EW weights
	var weights = PortfolioAllocation.equalWeights(nbAssets);
  
    // Check the number of output weights
	assert.equal(weights.length, nbAssets, 'EW - Number of weights');
	
	// Compare EW weights to exact closed-form formula
	var expectedWeights = [];
	for (var i = 0; i < nbAssets; ++i) {
		expectedWeights[i] = 1/nbAssets;
	}
    for (var i = 0; i < weights.length; ++i) {
	  assert.equal(weights[i], expectedWeights[i], 'EW - Values ' + i);
	}
  }  
});


QUnit.test('Equal risk budget portfolio', function(assert) {     
  // Reference: Carvalho, Raul Leote de and Xiao, Lu and Moulin, Pierre, Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description
  
  // Random data, n assets case
  {
	// Generate a random number of assets
	var nbAssets = Math.floor(Math.random()*(50-1+1) + 1); // max 50 min 1
	
	// Generate n random variances
	var sigma = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		sigma[i] = 1 - Math.random(); // belongs to ]0,1];
	}
	
	// Compute ERB weights
	var weights = PortfolioAllocation.equalRiskBudgetWeights(sigma);
	
	// Check the number of output weights
	assert.equal(weights.length, nbAssets, 'ERB - Number of weights');
	
	// Compare ERB weights to exact closed-form formula, c.f. formula 3 of the reference
	var denom = 0;
	for (var i = 0; i < nbAssets; ++i) { 
		denom += 1/Math.sqrt(sigma[i]); 
	}
	var expectedWeights = [];
	for (var i = 0; i < nbAssets; ++i) {
		expectedWeights[i] = 1/Math.sqrt(sigma[i]) / denom;
	}
	for (var i = 0; i < nbAssets; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'ERB - Values ' + i);  
	}
  }
});


QUnit.test('Equal risk contributions portfolio', function(assert) {    
  // Reference: Maillard, S., Roncalli, T., Teiletche, J.: The properties of equally weighted risk contribution portfolios. J. Portf. Manag. 36, 60–70 (2010)
  
  // Example 1 p. 65 of reference with constant correlation equals to 0
  {
	  var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]]); 
	  var expectedWeights =  [48, 24, 16, 12];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC - Values #1 ' + i);
	  }
  }
  
  
  // Example 2 p. 65 of reference (note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed)
  {
	  var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]]); 
	  var expectedWeights =  [38.4, 19.2, 24.3, 18.2];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		  assert.equal(Math.round(1000 * weights[i])/10, expectedWeights[i], 'ERC - Values #2 ' + i);
	  }
  }
  
  // Random data, 2 assets case, c.f. p 62 of reference
  {
	  // Generate a random 2*2 covariance matrix
	  var covMat = this.randomCovarianceMatrix();
	  
	  // Compute ERC weights
	  var weights = PortfolioAllocation.equalRiskContributionWeights(covMat);
	  
	  // Compare ERC weights to exact closed-form formula, c.f. p 62 of reference
	  var var1 = covMat[0][0];
	  var var2 = covMat[1][1];
	  var std1 = Math.sqrt(var1);
	  var std2 = Math.sqrt(var2);
	  var expectedWeights =  [(1/std1) / (1/std1 + 1/std2), (1/std2) / (1/std1 + 1/std2)];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'ERC - Values #3 ' + i);  
	  }
  }
});


QUnit.test('Risk budgeting portfolio', function(assert) {	
	// Reference: Bruder, Benjamin and Roncalli, Thierry, Managing Risk Exposures Using the Risk Budgeting Approach
    
    // Random data, 2 assets case, c.f. p 6 of reference
	{
	    // Generate a random 2*2 covariance matrix
	    var covMat = this.randomCovarianceMatrix();
	  
		// Generate random risk budgets
		var b = 1 - Math.random(); // belongs to ]0,1]

		// Compute RB weights
		var weights = PortfolioAllocation.riskBudgetingWeights(covMat, [b, 1 - b]);
		
		// Compare RB weights to exact closed-form formula, c.f. p 6 of reference
	    var var1 = covMat[0][0];
	    var var2 = covMat[1][1];
		var cov = covMat[0][1];
		var std1 = Math.sqrt(var1);
		var std2 = Math.sqrt(var2);
		var corr = cov/(std1*std2);
		var w_star = (b - 0.5)*corr*std1*std2 - b*var2 + std1*std2*Math.sqrt( (b - 0.5)*(b - 0.5)*corr*corr + b*(1 - b) );
		w_star /= (1 - b)*var1 - b*var2 + 2*(b - 0.5)*corr*std1*std2;	
		var expectedWeights = [w_star, 1 - w_star];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'RB - Values #1 ' + i);  
		}
	}
	
	
    // Random data, n assets case, 0 correlation, c.f. p 7 of reference
    // The associated matrix is always definite positive as long as all variances are strictly positive
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(20-2+1) + 2); // max 20 min 2
		
		// Generate a random n*n covariance matrix with 0 correlation
		var sigma = new Array(nbAssets);
		for (var i = 0; i < nbAssets; ++i) {
			sigma[i] = new Array(nbAssets);
			for (var j = 0; j < nbAssets; ++j) { 
				sigma[i][j] = 0; // 0 correlation
			}
			sigma[i][i] = 1 - Math.random(); // belongs to ]0,1];
		}
		
		// Generate n random risk budgets
		var rb = new Array(nbAssets);
		var sum_rb = 0;
		for (var i = 0; i < nbAssets; ++i) { 
			rb[i] = 1 - Math.random(); // belongs to ]0,1]; so that rb[i] > 0
			sum_rb += rb[i]; 
		}
		for (var i = 0; i < nbAssets; ++i) {
			rb[i] /= sum_rb; // normalization, so that sum(rb[i]) == 1
		}

		// Compute RB weights
		var weights = PortfolioAllocation.riskBudgetingWeights(sigma, rb, {eps: 1e-10, maxIter: 10000});
		
		// Compare RB weights to exact closed-form formula, c.f. p 7 of reference
		var denom = 0;
		for (var i = 0; i < nbAssets; ++i) { 
			denom += Math.sqrt(rb[i]) * 1/Math.sqrt(sigma[i][i]); 
		}
		for (var i = 0; i < nbAssets; ++i) {
			assert.equal(Math.abs(weights[i] - Math.sqrt(rb[i]) * 1/Math.sqrt(sigma[i][i]) / denom) <= 1e-8, true, 'RB - Values #2 ' + i);  
		}
	}
});


QUnit.test('Most diversified portfolio', function(assert) {    
	// Reference: Choueifaty, Yves and Froidure, Tristan and Reynier, Julien, Properties of the Most Diversified Portfolio (July 6, 2011). Journal of Investment Strategies, Vol.2(2), Spring 2013, pp.49-70.

	// Simple 2 assets case, p. 9 of the reference
	// Note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed
	{
		// Covariance matrix
		var sigma = [[0.0400, 0.0100], [0.0100, 0.0100]];
		
		// Compute MDP weights
		var weights = PortfolioAllocation.mostDiversifiedWeights(sigma);
		
		// Compare MDP weights to values provided p. 9 of the reference
		var expectedWeights = [0.33, 0.67];
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP - Values #1 0');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP - Values #1 1');
	}
	
	// Duplicate 3 assets is skipped, as not a real example (asset A duplication is a polico, and thus would never be selected)
	
	// Leverage invariance 2 assets case, p. 10 of the reference
	// Note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed
	{
		// Covariance matrix
		var sigma = [[0.0025, 0.0025], [0.0025, 0.0100]];
		
		// Compute MDP weights
		var weights = PortfolioAllocation.mostDiversifiedWeights(sigma);
		
		// Compare MDP weights to values provided p. 11 of the reference
		var expectedWeights = [0.67, 0.33];
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP - Values #2 0');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP - Values #2 1');
	}
	
	// Polico invariance 3 assets case, p. 10 of the reference
	// Note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed
	{
		// Covariance matrix
		var sigma = [[0.0400, 0.0100, 0.0225], [0.0100, 0.0100, 0.0075], [0.0225, 0.0075, 0.0131]];
		
		// Compute MDP weights
		var weights = PortfolioAllocation.mostDiversifiedWeights(sigma);
		
		// Compare MDP weights to values provided p. 11 of the reference
		var expectedWeights = [0.33, 0.67, 0];
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP - Values #3 0');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP - Values #3 1');
		assert.equal(Math.abs(weights[2] - expectedWeights[2]) <= 1e-2, true, 'MDP - Values #3 2');
	}
	
	
	// Reference: Toward Maximum Diversification by Y. Choueifaty, Y. Coignard, The Journal of Portfolio Management, Fall 2008, Vol. 35, No. 1: pp. 40-51
	
	// Example 2 p.41 of the reference
	// Note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed
	{
		// Covariance matrix
		var sigma = [[0.0100, 0.0090, 0.0010], [0.0090, 0.0100, 0.0010], [0.0010, 0.0010, 0.0100]];
		
		// Compute MDP weights
		var weights = PortfolioAllocation.mostDiversifiedWeights(sigma);

		// Compare MDP weights to values provided p. 11 of the reference
		var expectedWeights = [25.7, 25.7, 48.6];
		assert.equal(Math.round(weights[0]*1000)/10, expectedWeights[0], true, 'MDP - Values #4 0');
		assert.equal(Math.round(weights[1]*1000)/10, expectedWeights[1], true, 'MDP - Values #4 1');
		assert.equal(Math.round(weights[2]*1000)/10,  expectedWeights[2], true, 'MDP - Values #4 2');
	}
	
	
	// Static data: correlation matrix built with randcorr Matlab gallery function
	{
		var corr = [[1.0000, 0.1966, -0.2071, -0.2084, 0.2021], [0.1966, 1.0000, 0.1993, 0.2008, -0.1963], [-0.2071, 0.1993, 1.0000, -0.2037, 0.1960], [-0.2084, 0.2008, -0.2037, 1.0000, 0.1978], [0.2021, -0.1963, 0.1960, 0.1978, 1.0000]];
		var nbAssets = corr.length;
		
		{	
			// In case of a correlation matrix, maximizing the diversification ratio
			// is equivalent to minimizing the variance, c.f. formula 2 of the reference
			//
			// Variance minimization with the correlation matrix corr has been implemented
			// in Matlab, with CVX library, and output weights copy/pasted here.

			// Compute MDP weights
			var weights = PortfolioAllocation.mostDiversifiedWeights(corr, {eps: 1e-10});
			
			// Compare MDP weights to values obtained through Matlab CVX library usage
			var expectedWeights = [0.334077696834404, 0.000000171869935, 0.332780568181053, 0.333141372759862, 0.000000190356559];
			for (var i = 0; i < nbAssets; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-6, true, 'MDP - Values #5 ' + i);
			}
		}
		
		{
			// In case of a covariance matrix, maximizing the diversification ratio
			// is equivalent to minimizing the variance using the correlation matrix,
			// plus rescaling the obtained weights by the assets volatilities and normalizing
			// them if no cash, c.f. discussion above formula 4 of the reference.
			
			// Generate random standard deviations
			var stddevs = new Array(nbAssets);
			for (var i = 0; i < nbAssets; ++i) { 
				stddevs[i] = 1 - Math.random(); // belongs to ]0,1]
			}
			
			// Generate the associated covariance matrix from the correlation matrix above
			var sigma = corr.slice();
			for (var i = 0; i < nbAssets; ++i) { 
				for (var j = 0; j < nbAssets; ++j) { 
					sigma[i][j] *= stddevs[i]*stddevs[j];
				}
			}
			
			// Compute MDP weights
			var weights = PortfolioAllocation.mostDiversifiedWeights(sigma);
			
			// Compute expected weights by rescaling + normalizing the weights obtained from the Matlab CVX library usage
			// on the correlation matrix above
			var expectedWeights = [0.334077696834404, 0.000000171869935, 0.332780568181053, 0.333141372759862, 0.000000190356559]; // correlation weights
			var sumWeights = 0.0;
			for (var i = 0; i < nbAssets; ++i) { 
				expectedWeights[i] /= stddevs[i]; // correlation rescaled weights
				sumWeights += expectedWeights[i];
			}
			for (var i = 0; i < nbAssets; ++i) { 
				expectedWeights[i] /= sumWeights; // correlation rescaled + normalized weights = expected weights
			}
			
			// Compare MDP weights to expected weights
			for (var i = 0; i < nbAssets; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-4, true, 'MDP - Values #6 ' + i);
			}
		}
	}
});



QUnit.test('Minimum correlation heuristic portfolio', function(assert) {    
  // Reference: David Varadi, Michael Kapler, Henry Bee, Corey Rittenhouse, Sep. 2012, The Minimum Correlation Algorithm: A Practical Diversification Tool
  
  // Example of the appendix, with MinCorr 2 example being inversed with MinCorr example...
  {
	  // Data
	  var sigma = [[0.019600000000000, 0.022680000000000, 0.026180000000000], [0.022680000000000, 0.032400000000000, 0.027720000000000], [0.026180000000000, 0.027720000000000, 0.048400000000000]];
	  var nbAssets = sigma.length;

	  // Compute MinCorr weights
	  var weights = PortfolioAllocation.minCorrWeights(sigma);
	  
	  // Compare MinCorr weights to expected weights
	  var expectedWeights = [0.21, 0.31, 0.48];
	  for (var i = 0; i < nbAssets; ++i) { 
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-2, true, 'MinCorr - Values #1 ' + i);
	  }
  }

  // Limit case: 2 assets, must be equal to ERB
  {
	    // Generate a random 2*2 covariance matrix
	    var covMat = this.randomCovarianceMatrix();
		
		// Compute MinCorr weights
		var weights = PortfolioAllocation.minCorrWeights(covMat);
		
		// Compute ERB weights
	    var var1 = covMat[0][0];
	    var var2 = covMat[1][1];
		var expectedWeights = PortfolioAllocation.equalRiskBudgetWeights([var1, var2]);
		
		// Compare the weights
		for (var i = 0; i < 2; ++i) { 
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'MinCorr - 2 assets ' + i);
		}
  }  
});


QUnit.test('Proportional minimum variance heuristic portfolio', function(assert) {    
  // Reference: Minimum Variance Algorithm (MVA) Excel Sheet
  
  // Example of the Excel sheet
  {
	  // Data
	  var cov = [[0.000090, 0.000044, 0.000028, 0.000034],[0.000044, 0.000084, 0.000068, 0.000039],[0.000028, 0.000068, 0.000101, 0.000036],[0.000034, 0.000039, 0.000036, 0.000039]];
	  var nbAssets = cov.length;

	  // Compute MVA weights
	  var weights = PortfolioAllocation.minVarWeights(cov);
	  
	  // Compare MVA weights to expected weights
	  var expectedWeights = [0.18, 0.07, 0.07, 0.68];
	  for (var i = 0; i < nbAssets; ++i) { 
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-2, true, 'MinVar - Values #1 ' + i);
	  }
  }
});



QUnit.test('Cluster risk parity portfolio', function(assert) {    
  // Error cases
  {
	// Dummy covariance matrix
	var sigma = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

	// Unsupported clustering method
    assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'none' }) },
		new Error('unsupported clustering method'),
		"CRP - Manual clustering checks unsupported clustering method");

    // In case manual clusters are provided for n assets, they must form a partition of the set [1..n]
	// Test no clusters
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual' }) }, 
		new Error('missing asset index: 1'),
		"CRP - Manual clustering checks no clusters");
	
	// Test the presence of an empty cluster
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual', clusters: [[1, 2, 3, 4], []] }) }, 
		new Error('empty cluster at index: 1'),
		"CRP - Manual clustering checks empty cluster");
	
	// Test asset index out of bounds #1
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual', clusters: [[0, 1, 2, 3]] }) }, 
		new Error('asset index out of bounds: 0'),
		"CRP - Manual clustering checks out of bounds #1");
		
	// Test asset index out of bounds #2
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual', clusters: [[1, 2, 3, 5]] }) }, 
		new Error('asset index out of bounds: 5'),
		"CRP - Manual clustering checks out of bounds #2")

	// Test missing asset index
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual', clusters: [[1, 2, 3]] }) }, 
		new Error('missing asset index: 4'),
		"CRP - Manual clustering checks missing asset index")
	
	// Test duplicate asset index
	assert.throws(function() { 
		PortfolioAllocation.clusterRiskParityWeights(sigma, { clusteringMode: 'manual', clusters: [[1, 1, 2, 3]] }) }, 
		new Error('duplicate asset index: 1'),
		"CRP - Manual clustering checks duplicate asset index")
  } 

  
  // Limit cases
  {
	// Example taken from the ERC
	var sigma = [[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]];
	var ercWeights = [38.4, 19.2, 24.3, 18.2];

	// Limit case #1: in case all assets belong to a unique cluster, the CRP weights must be equal to the ERC weights
	{
		var opt = { clusteringMode: 'manual', clusters: [[1, 2, 3, 4]] };
		var weights = PortfolioAllocation.clusterRiskParityWeights(sigma, opt); 
		var expectedWeights = ercWeights;
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(1000 * weights[i])/10, expectedWeights[i], 'CRP - Values #1 ' + i);
		}
	}
  
	  // Limit case #2: in case each assets belong to one cluster, the CRP weights must be equal to the ERC weights
	  {
		  var opt = { clusteringMode: 'manual', clusters: [[1], [2], [3], [4]] };
		  var weights = PortfolioAllocation.clusterRiskParityWeights(sigma, opt); 
		  var expectedWeights = ercWeights;
		  for (var i = 0; i < expectedWeights.length; ++i) {
			  assert.equal(Math.round(1000 * weights[i])/10, expectedWeights[i], 'CRP - Values #2 ' + i);
		  }
	  }
  }
  
  
  // In case of 2 clusters of 2 assets each, the final weights must verify w_1f = p * w_1erc, w_2f = p * w_2erc,
  // w_3f = (1-p) * w_3erc, w_4f = (1-p) * w_4erc with w_ierc the weight for asset i obtained though the ERC applied
  // intra clusters and p a real number belonging to ]0,1[ dependant on the clusters variances/correlation.
  {
      // Generate two random 2*2 covariance matrices, one for each of the two clusters
	  var sigma_c1 = this.randomCovarianceMatrix();
	  var sigma_c2 = this.randomCovarianceMatrix();

	  // Create a 4x4 bloc matrix using the two 2x2 matrices above
	  var sigma = [[sigma_c1[0][0], sigma_c1[0][1], 0, 0], [sigma_c1[1][0], sigma_c1[1][1], 0, 0], [0, 0, sigma_c2[0][0], sigma_c2[0][1]], [0, 0, sigma_c2[1][0], sigma_c2[1][1]]];
		
      // Compute CRP weights for the whole portfolio
	  var opt = { clusteringMode: 'manual', clusters: [[1,2], [3,4]] };
	  var weights = PortfolioAllocation.clusterRiskParityWeights(sigma, opt); 
	  
  	  // Compute ERC weights for each cluster
	  var weights_c1 = PortfolioAllocation.equalRiskContributionWeights(sigma_c1); 
	  var weights_c2 = PortfolioAllocation.equalRiskContributionWeights(sigma_c2); 

	  // Ensure the relationships described above are satisfied
	  var p = weights[0]/weights_c1[0];
	  assert.equal(Math.abs(p - weights[1]/weights_c1[1]) <= 1e-8, true, 'CRP - Values #3 0');
	  assert.equal(Math.abs(1-p - weights[2]/weights_c2[0]) <= 1e-8, true, 'CRP - Values #3 1');
	  assert.equal(Math.abs(1-p - weights[3]/weights_c2[1]) <= 1e-8, true, 'CRP - Values #3 2');
  }
  

  
  // CRP using FTCA, using static data
  {
	// Example taken from the ERC
	var sigma = [[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]];
	
	// Call CRP with no options, which should defaults to FTCA
	var expectedWeights = [ 0.3262379197337176, 0.16311896054014643, 0.29179607083498427, 0.21884704889115175 ];
	var weights = PortfolioAllocation.clusterRiskParityWeights(sigma);
	for (var i = 0; i < sigma[0].length; ++i) { 
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'CRP - Values #4 ' + i);
	}
	  
	// Call CRP with FTCA mode, which should output the same result as above
	var opt = { clusteringMode: 'ftca' };
	weights = PortfolioAllocation.clusterRiskParityWeights(sigma, opt); 		
	for (var i = 0; i < sigma[0].length; ++i) { 
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'CRP - Values #5 ' + i);
	}

  }
});


QUnit.test('Global minimum variance portfolio', function(assert) {    
  // Reference: Portfolio Optimization versus Risk-Budgeting Allocation, Thierry Roncalli, WG RISK ESSEC, January 18, 2012
  // Examples using static data
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0396, 0.0398], [0.0398, 0.0400]]); 
	  var expectedWeights =  [1, 0];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #1 ' + i);
	  }
  }
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0396, 0.0414], [0.0396, 0.0484, 0.0455], [0.0414, 0.0455, 0.0529]]); 
	  var expectedWeights =  [0.9565927119697761, 0.043407288030223846, 0];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #2 ' + i);
	  }
  }
  // TODO: This one does not pass; pending the Markowitz algorithm implementation to decide if the reference is false,
  // as the volatility obtained by the computed portfolio is 0.1987534683949592, which is lower than the volatility of 0.2017 in the reference !
  /*
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0374, 0.0391], [0.0374, 0.0484, 0.0430], [0.0391, 0.0430, 0.0529]]); 
	  console.log(weights)
	  var expectedWeights =  [0.7009, 0.2378, 0.0613];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #3 ' + i);
	  }
  }
  */
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0418, 0.0437], [0.0418, 0.0484, 0.0481], [0.0437, 0.0481, 0.0529]]); 
	  var expectedWeights =  [1, 0, 0];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #4 ' + i);
	  }
  }
  
  // Reference: Understanding the Impact of Weights Constraints in Portfolio Theory, Thierry Roncalli
  // Example using static data
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0225,0.0030,0.0150,0.0225], [0.0030,0.0400,0.0350,0.0240], [0.0150,0.0350,0.0625,0.0600], [0.0225,0.0240,0.0600,0.0900]]); 
	  var expectedWeights =  [0.6548673025499347, 0.3451326974500652, 0, 0];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #5 ' + i);
	  }
  }
  	
  // Reference: Private communication with TrendXplorer
  // Example using static data
  {
	  var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.03401428,0.0333167,-0.00614739,0.00926415,-0.0064081],[0.0333167,0.06323421,-0.00855552,0.02245369,-0.00480642],[-0.00614739,-0.00855552,0.01444902,-0.00432445,0.00690744],[0.00926415,0.02245369,-0.00432445,0.02622712,0.0016983],[-0.0064081,-0.00480642,0.00690744,0.00169834,0.0116492]]); 
	  var expectedWeights =  [0.22372329306378916, 0, 0.30887622610697146, 0.13382010508052292, 0.33358037574871646];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #6 ' + i);
	  }
  }
  
});



QUnit.test('Grid search portfolio', function(assert) {    
	// Test unsupported optimisation method
	{
	    
        assert.throws(function() { 
    		PortfolioAllocation.gridSearchWeights(2,  function (arr) { return 1; }, {optimisationMethod: 'none'}) },
    		new Error('unsupported optimisation method'),
    		"Grid search portfolio - Unsupported optimisation method");
	}
		
	// Test rational grid search using static data
	{
		// Objective function: portfolio variance, for three assets
		// Covariance matrix taken from GMV test case #2
		function portfolio_variance_three_assets(arr) { 
			var covMat = [[0.0400, 0.0396, 0.0414], [0.0396, 0.0484, 0.0455], [0.0414, 0.0455, 0.0529]];
			
			return arr[0]*arr[0]*covMat[0][0] + arr[1]*arr[1]*covMat[1][1] + arr[2]*arr[2]*covMat[2][2] +
			       2*arr[0]*arr[1]*covMat[0][1] + 2*arr[0]*arr[2]*covMat[0][2] + 2*arr[1]*arr[2]*covMat[1][2];
		}
		
		// From GMV test case #2, expected (exact) weights are [0.9565927119697761, 0.043407288030223846, 0];
		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets), 
	                	PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'rationalGridSearch', rationalGridSearch: {k: 3}}), 
		                'Grid search portfolio - Default values');

		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'rationalGridSearch', rationalGridSearch: {k: 3}}), 
	                	[[1,0,0]], 
		                'Grid search portfolio - Values #1');

		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'rationalGridSearch', rationalGridSearch: {k: 10}}), 
	                	[[1,0,0]], 
		                'Grid search portfolio - Values #2');
		                
		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'rationalGridSearch', rationalGridSearch: {k: 100}}), 
	                	[[0.96,0.04,0]], 
		                'Grid search portfolio - Values #3');
		                
		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'rationalGridSearch', rationalGridSearch: {k: 1000}}), 
	                	[[0.957,0.043,0]], 
		                'Grid search portfolio - Values #4');
	}

});
