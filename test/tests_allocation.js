// ------------------------------------------------------------
QUnit.module('Assets allocation module', {
  before: function() {
    // 
  }
});


QUnit.test('Equal weights portfolio', function(assert) {    
  // Limit case, one asset
  var weights = PortfolioAllocation.equalWeights(1);
  assert.equal(weights.length, 1, 'Equal weights with one asset');
  assert.equal(weights[0], 1, 'Equal weights with one asset #2');

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
