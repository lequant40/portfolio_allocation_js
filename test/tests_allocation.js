// ------------------------------------------------------------
QUnit.module('Assets allocation module', {
  before: function() {
    // 
  }
});


QUnit.test('Equal weights portfolio', function(assert) {    
  // Limit case, one asset
  {
	  var weights = PortfolioAllocation.equalWeights(1);
	  assert.equal(weights.length, 1, 'Equal weights with one asset');
	  assert.equal(weights[0], 1, 'Equal weights with one asset #2');
  }

  // General case, generate random values for testing
  {
    var nbAssets = Math.floor(Math.random()*(50-2+1) + 2); // max 50 min 2
    var weights = PortfolioAllocation.equalWeights(nbAssets);
  
    assert.equal(weights.length, nbAssets, 'Equal weights');
    for (var j = 0; j < weights.length; ++j) {
	  assert.equal(weights[j], 1/nbAssets, 'Equal weights #2');
	}
  }  
});


QUnit.test('Equal risk budget portfolio', function(assert) {     
  // Random data, n assets case
  {
	var nbAssets = Math.floor(Math.random()*(50-3+1) + 3); // max 50 min 3
	
	// Generate n variances
	var sigma = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		sigma[i] = 1 - Math.random(); // belongs to ]0,1];
	}
	
	// Compute ERB weights
	var weights = PortfolioAllocation.equalRiskBudgetWeights(sigma);
	
	// Compare ERB weights to exact closed-form formula, c.f. formula 3 of the reference
	var denom = 0;
	for (var i = 0; i < nbAssets; ++i) { 
		denom += 1/Math.sqrt(sigma[i]); 
	}
	for (var i = 0; i < nbAssets; ++i) {
		assert.equal(Math.abs(weights[i] - 1/Math.sqrt(sigma[i]) / denom) <= 1e-8, true, 'ERB #1/' + i);  
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
		assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC #1');
	  }
  }
  
  
  // Example 2 p. 65 of reference (note: covariance matrix is provided here v.s. correlation matrix in the reference => matrix operations were needed)
  {
	  var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]]); 
	  var expectedWeights =  [38.4, 19.2, 24.3, 18.2];
	  for (var i = 0; i < expectedWeights.length; ++i) {
		  assert.equal(Math.round(1000 * weights[i])/10, expectedWeights[i], 'ERC #2');
	  }
  }
  
  // Random data, 2 assets case, c.f. p 62 of reference
  {
	  // Generate a 2*2 covariance matrix
	  var var1 = 1 - Math.random(); // belongs to ]0,1]
	  var var2 = 1 - Math.random(); // belongs to ]0,1]
	  var cov = 2*Math.random() - 1; // belongs to ]-1,1[
	  while (var1*var2 - cov*cov <= 0) { // The covariance matrix must be DP, so that the covariance cannot be any value: determinant must be strictly positive
	    cov = 2*Math.random() - 1;
	  }
	  
	  // Compute ERC weights
	  var weights = PortfolioAllocation.equalRiskContributionWeights([[var1, cov], [cov, var2]]);
	  
	  // Compare ERC weights to exact closed-form formula, c.f. p 62 of reference
	  var std1 = Math.sqrt(var1);
	  var std2 = Math.sqrt(var2);
	  var expectedWeights =  [(1/std1) / (1/std1 + 1/std2), (1/std2) / (1/std1 + 1/std2)];
	  assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-8, true, 'ERC #3/1 ' + i);  
	  assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-8, true, 'ERC #3/2 ' + i);
  }
});


QUnit.test('Risk budgeting portfolio', function(assert) {	
	// Reference: Bruder, Benjamin and Roncalli, Thierry, Managing Risk Exposures Using the Risk Budgeting Approach
    
    // Random data, 2 assets case, c.f. p 6 of reference
	{
		// Generate a 2*2 covariance matrix
		var var1 = 1 - Math.random(); // belongs to ]0,1]
		var var2 = 1 - Math.random(); // belongs to ]0,1]
		var cov = 2*Math.random() - 1; // belongs to ]-1,1[
		while (var1*var2 - cov*cov <= 0) { // The covariance matrix must be DP, so that the correlation cannot be any value: determinant must be strictly positive
		  cov = 2*Math.random() - 1;
		}
	  
		// Generate risk budgets
		var b = 1 - Math.random(); // belongs to ]0,1]

		// Compute RB weights
		var weights = PortfolioAllocation.riskBudgetingWeights([[var1, cov], [cov, var2]], [b, 1 - b]);
		
		// Compare RB weights to exact closed-form formula, c.f. p 6 of reference
		var std1 = Math.sqrt(var1);
		var std2 = Math.sqrt(var2);
		var corr = cov/(std1*std2);
		var w_star = (b - 0.5)*corr*std1*std2 - b*var2 + std1*std2*Math.sqrt( (b - 0.5)*(b - 0.5)*corr*corr + b*(1 - b) );
		w_star /= (1 - b)*var1 - b*var2 + 2*(b - 0.5)*corr*std1*std2;	
		var expectedWeights = [w_star, 1 - w_star];
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-8, true, 'RB #1/1');  
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-8, true, 'RB #1/2');
	}
	
	
    // Random data, n assets case, 0 correlation, c.f. p 7 of reference
    // The associated matrix is always definite positive as long as all variances are strictly positive
	{
		var nbAssets = Math.floor(Math.random()*(50-3+1) + 3); // max 50 min 3
		
		// Generate a n*n covariance matrix with 0 correlation
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
			assert.equal(Math.abs(weights[i] - Math.sqrt(rb[i]) * 1/Math.sqrt(sigma[i][i]) / denom) <= 1e-8, true, 'RB #2/' + i);  
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
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP #1/1');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP #1/2');
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
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP #2/1');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP #2/2');
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
		assert.equal(Math.abs(weights[0] - expectedWeights[0]) <= 1e-2, true, 'MDP #3/1');
		assert.equal(Math.abs(weights[1] - expectedWeights[1]) <= 1e-2, true, 'MDP #3/2');
		assert.equal(Math.abs(weights[2] - expectedWeights[2]) <= 1e-2, true, 'MDP #3/3');
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
		assert.equal(Math.round(weights[0]*1000)/10, expectedWeights[0], true, 'MDP #4/1');
		assert.equal(Math.round(weights[1]*1000)/10, expectedWeights[1], true, 'MDP #4/2');
		assert.equal(Math.round(weights[2]*1000)/10,  expectedWeights[2], true, 'MDP #4/3');
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
				assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-6, true, 'MDP #5/' + i);
			}
		}
		
		{
			// In case of a covariance matrix, maximizing the diversification ratio
			// is equivalent to minimizing the variance, using the correlation matrix,
			// plus rescaling the obtained weights by the assets volatilities + normalizing
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
				assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-4, true, 'MDP #6/' + i);
			}
		}
	}
});



QUnit.test('Minimum correlation heuristic portfolio', function(assert) {    
  // Reference: David Varadi, Michael Kapler, Henry Bee, Corey Rittenhouse, Sep. 2012, The Minimum Correlation Algorithm: A Practical Diversification Tool
  
  // Example of the appendix, with MinCorr 2 example being inversed with MinCorr example...
  {
	  // Data
	  //var corr = [[1.0000, 0.90, 0.85], [0.90, 1.0000, 0.70], [0.85, 0.70, 1.0000]];
	  //var vars = [0.14*0.14, 0.18*0.18, 0.22*0.22];
	  var sigma = [[0.019600000000000, 0.022680000000000, 0.026180000000000], [0.022680000000000, 0.032400000000000, 0.027720000000000], [0.026180000000000, 0.027720000000000, 0.048400000000000]];
	  var nbAssets = sigma.length;

	  // Compute MinCorr weights
	  var weights = PortfolioAllocation.minCorrWeights(sigma);
	  
	  // Compare MinCorr weights to expected weights
	  var expectedWeights = [0.21, 0.31, 0.48];
	  for (var i = 0; i < nbAssets; ++i) { 
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-2, true, 'MinCorr weights ' + i);
	  }
  }

  // Limit case: 2 assets, must be equal to ERB
  {
		// Generate a 2*2 covariance matrix
		var var1 = 1 - Math.random(); // belongs to ]0,1]
		var var2 = 1 - Math.random(); // belongs to ]0,1]
		var cov = 2*Math.random() - 1; // belongs to ]-1,1[
		while (var1*var2 - cov*cov < 0) { // The covariance matrix must be DSP, so that the correlation cannot be any value: determinant must be positive
		  cov = 2*Math.random() - 1;
		}
		
		// Compute MinCorr weights
		var weights = PortfolioAllocation.minCorrWeights([[var1, cov], [cov, var2]]);
		
		// Compute ERB weights
		var expectedWeights = PortfolioAllocation.equalRiskBudgetWeights([var1, var2]);
		
		// Compare the weights
		for (var i = 0; i < 2; ++i) { 
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'MinCorr weights 2 assets ' + i);
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
		assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-2, true, 'MinVar weights ' + i);
	  }
  }
});