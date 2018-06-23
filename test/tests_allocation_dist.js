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
	
	// Reference: Xi Bai, Katya Scheinberg, Reha Tutuncu, Least-squares approach to risk parity in portfolio selection
	//
	// To be noted that this test also test the option of providing the portfolio volatility in output.
	{
		var covMar = [[94.868,33.750,12.325,-1.178,8.778],
					[33.750,445.642,98.955,-7.901,84.954],
					[12.325,98.955,117.265,0.503,45.184],
					[-1.178,-7.901,0.503,5.460,1.057],
					[8.778,84.954,45.184,1.057,34.126]];
		
		var weights = PortfolioAllocation.equalRiskContributionWeights(covMar, {outputPortfolioVolatility: true}); 
		var expectedWeights =  [0.125, 0.047, 0.083, 0.613, 0.132];
		var expectedVolatility = 3.0406150258182514;
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[0][i] - expectedWeights[i]) <= 1e-3, true, 'ERC - Values #4 ' + i);
		}
		assert.equal(Math.abs(weights[1] - expectedVolatility) <= 1e-3, true, 'ERC - Values #4, volatility');
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
		assert.equal(Math.abs(Math.round(weights[0]*1000)/10 - expectedWeights[0]) <= 1e-2, true, 'MDP - Values #4 0');
		assert.equal(Math.abs(Math.round(weights[1]*1000)/10 - expectedWeights[1]) <= 1e-2, true, 'MDP - Values #4 1');
		assert.equal(Math.abs(Math.round(weights[2]*1000)/10- expectedWeights[2]) <= 1e-2, true, 'MDP - Values #4 2');
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
	{
		var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0396, 0.0398], [0.0398, 0.0400]]); 
		var expectedWeights =  [1, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #1 ' + i);
		}
	}
	{
		var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0396, 0.0414], [0.0396, 0.0484, 0.0455], [0.0414, 0.0455, 0.0529]]); 
		var expectedWeights =  [0.9565217391304353, 0.04347826086956469, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #2 ' + i);
		}
	}
	// Note: The volatility obtained by the computed portfolio is 0.019751470588235294, which is lower than the volatility of 0.2017 in the reference,
	// associated with the commented weights !
	{
		var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0374, 0.0391], [0.0374, 0.0484, 0.0430], [0.0391, 0.0430, 0.0529]]); 
		//var expectedWeights =  [0.7009, 0.2378, 0.0613];//
		var expectedWeights =  [0.8088235294117648, 0.19117647058823511, 0 ];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #3 ' + i);
		}
	}

	{
		var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0400, 0.0418, 0.0437], [0.0418, 0.0484, 0.0481], [0.0437, 0.0481, 0.0529]]); 
		var expectedWeights =  [1, 0, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #4 ' + i);
		}
	}

	// Reference: Understanding the Impact of Weights Constraints in Portfolio Theory, Thierry Roncalli
	{
		var weights = PortfolioAllocation.globalMinimumVarianceWeights([[0.0225,0.0030,0.0150,0.0225], [0.0030,0.0400,0.0350,0.0240], [0.0150,0.0350,0.0625,0.0600], [0.0225,0.0240,0.0600,0.0900]]); 
		var expectedWeights =  [0.6548672566371683, 0.34513274336283173, 0, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #5 ' + i);
		}
	}
	
	// Reference: Private communication with TrendXplorer
	{
		var covMat = [[0.03401428,0.0333167,-0.00614739,0.00926415,-0.0064081],
					  [0.0333167,0.06323421,-0.00855552,0.02245369,-0.00480642],
					  [-0.00614739,-0.00855552,0.01444902,-0.00432445,0.00690744],
					  [0.00926415,0.02245369,-0.00432445,0.02622712,0.0016983],
					  [-0.0064081,-0.00480642,0.00690744,0.00169834,0.0116492]];
		
		var weights = PortfolioAllocation.globalMinimumVarianceWeights(covMat); 
		var expectedWeights =  [0.22146166011081053, 0, 0.3167829128512612, 0.13743827218475788, 0.32431715485317053];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #6 ' + i);
		}
	}
	
	// Reference: Xi Bai, Katya Scheinberg, Reha Tutuncu, Least-squares approach to risk parity in portfolio selection
	{
		var covMat = [[94.868,33.750,12.325,-1.178,8.778],
					  [33.750,445.642,98.955,-7.901,84.954],
					  [12.325,98.955,117.265,0.503,45.184],
					  [-1.178,-7.901,0.503,5.460,1.057],
					  [8.778,84.954,45.184,1.057,34.126]];
			
		var weights = PortfolioAllocation.globalMinimumVarianceWeights(covMat); 
		var expectedWeights =  [0.050, 0.006, 0, 0.862, 0.082];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-3, true, 'GMV - Values #7 ' + i);
		}
	}
	
	// Reference: Xi Bai, Katya Scheinberg, Reha Tutuncu, Least-squares approach to risk parity in portfolio selection
	// Constraints on assets weights
	{
		var covMat = [[94.868,33.750,12.325,-1.178,8.778],
					  [33.750,445.642,98.955,-7.901,84.954],
					  [12.325,98.955,117.265,0.503,45.184],
					  [-1.178,-7.901,0.503,5.460,1.057],
					  [8.778,84.954,45.184,1.057,34.126]];
			
		var minWeights = [0.05,0.05,0.05,0.05,0.05];
		var maxWeights = [0.35,0.35,0.35,0.35,0.35];
		var weights = PortfolioAllocation.globalMinimumVarianceWeights(covMat, { constraints: {minWeights: minWeights, maxWeights: maxWeights} }); 
		var expectedWeights =  [0.200,0.050,0.050,0.350,0.350];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-3, true, 'GMV - Values #8, constrained ' + i);
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
		PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'deterministic', rationalGrid: {k: 3}}), 
		'Grid search portfolio - Default values');

		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'deterministic', rationalGrid: {k: 3}}), 
		[[1,0,0]], 
		'Grid search portfolio - Values #1');

		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'deterministic', rationalGrid: {k: 10}}), 
		[[1,0,0]], 
		'Grid search portfolio - Values #2');
		
		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'deterministic', rationalGrid: {k: 100}}), 
		[[0.96,0.04,0]], 
		'Grid search portfolio - Values #3');
		
		assert.deepEqual(PortfolioAllocation.gridSearchWeights(3, portfolio_variance_three_assets, {optimisationMethod: 'deterministic', rationalGrid: {k: 1000}}), 
		[[0.957,0.043,0]], 
		'Grid search portfolio - Values #4');
	}
});



QUnit.test('Equal risk bounding portfolio', function(assert) {    
	// Reference: Cesarone, F. & Tardella F., Equal Risk Bounding is better than Risk Parity for portfolio selection, J Glob Optim (2017) 68: 439

	// Example of Table 1
	{
		// Data
		var sigma = [[1,-9/10, 3/5], [-9/10, 1,-1/5],[ 3/5, -1/5, 4]];
		var nbAssets = sigma.length;

		// Compute ERB weights
		var weights = PortfolioAllocation.equalRiskBoundingWeights(sigma);
		
		// Compare ERB weights to expected weights
		var expectedWeights = [0.5, 0.5, 0];
		for (var i = 0; i < nbAssets; ++i) { 
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'ERB - Values #1 ' + i);
		}
	}

	// Example of Table 2
	{
		// Data
		var sigma = [[1/100, 21/5000, -147/10000, 4/625, -31/2000], [21/5000, 1/25, -57/5000, -21/1250, 21/500], [-147/10000, -57/5000, 9/100, 36/625, 3/125], [4/625, -21/1250, 36/625, 4/25, 1/250], [-31/2000, 21/500, 3/125, 1/250,  1/4]];
		var nbAssets = sigma.length;

		// Compute ERB weights
		var weights = PortfolioAllocation.equalRiskBoundingWeights(sigma);
		
		// Compare ERB weights to expected weights
		var expectedWeights = [0.583, 0.157, 0.186, 0, 0.074];
		for (var i = 0; i < nbAssets; ++i) { 
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-3, true, 'ERB - Values #2 ' + i);
		}
	}

/*
	// Static example with 20 assets, to show the tractability of the ERB algorithm with this number of assets
	// This takes around 30 seconds on a 2012 computer
	{
	// Data
	var sigma = [[1.0000,0.3010,-0.2531,0.0497,-0.1286,0.0689,-0.0366,-0.0950,0.0502,-0.0342,0.0074,0.0107,-0.0971,0.2335,0.0807,-0.0024,-0.1245,0.0835,0.1783,-0.0989],
				[0.3010,1.0000,0.1298,0.0349,-0.0470,-0.1580,0.0884,-0.0551,0.0455,-0.0042,-0.1111,-0.0481,0.2739,-0.1916,-0.0774,0.0735,0.0239,0.0936,0.0030,-0.2237],
				[-0.2531,0.1298,1.0000,-0.0888,0.0696,-0.0881,0.1338,-0.1019,-0.2108,0.0422,-0.1493,-0.2719,0.1637,-0.2037,-0.2762,-0.0429,-0.1202,0.1024,0.1826,-0.2397],
				[0.0497,0.0349,-0.0888,1.0000,-0.0711,0.0303,-0.2361,-0.2027,0.2024,0.3856,-0.0313,0.1930,0.1567,-0.2357,0.0109,-0.1681,-0.0752,0.2626,0.0212,0.0834],
				[-0.1286,-0.0470,0.0696,-0.0711,1.0000,-0.2109,0.0067,0.0439,-0.2620,-0.2619,0.0822,0.1585,-0.1585,0.0025,-0.3311,-0.1409,-0.0738,-0.0535,-0.0684,0.2921],
				[0.0689,-0.1580,-0.0881,0.0303,-0.2109,1.0000,-0.1430,-0.2636,-0.0971,0.1451,-0.1107,-0.0961,-0.1054,-0.0189,0.0006,0.0678,0.1373,0.3409,-0.0739,0.1559],
				[-0.0366,0.0884,0.1338,-0.2361,0.0067,-0.1430,1.0000,-0.1179,0.1548,-0.0829,0.1414,-0.3829,-0.0553,0.0796,0.0092,0.0328,-0.1051,0.0139,-0.0700,0.1420],
				[-0.0950,-0.0551,-0.1019,-0.2027,0.0439,-0.2636,-0.1179,1.0000,-0.1681,0.1354,-0.0915,-0.1071,-0.1485,-0.0956,0.0881,-0.0646,0.1750,-0.2343,0.0476,-0.0378],
				[0.0502,0.0455,-0.2108,0.2024,-0.2620,-0.0971,0.1548,-0.1681,1.0000,0.0967,0.1938,-0.0046,0.3219,0.1068,0.0268,0.0781,0.0529,0.0346,0.1081,0.1386],
				[-0.0342,-0.0042,0.0422,0.3856,-0.2619,0.1451,-0.0829,0.1354,0.0967,1.0000,0.0773,-0.1698,-0.1185,-0.1201,0.1164,0.0458,-0.0040,-0.0958,-0.2451,0.0366],
				[0.0074,-0.1111,-0.1493,-0.0313,0.0822,-0.1107,0.1414,-0.0915,0.1938,0.0773,1.0000,-0.1579,-0.2839,0.0309,-0.1498,0.3240,-0.0849,-0.1001,0.0279,0.1015],
				[0.0107,-0.0481,-0.2719,0.1930,0.1585,-0.0961,-0.3829,-0.1071,-0.0046,-0.1698,-0.1579,1.0000,0.2488,0.1018,-0.0556,-0.0657,0.0473,-0.1634,-0.1715,0.0021],
				[-0.0971,0.2739,0.1637,0.1567,-0.1585,-0.1054,-0.0553,-0.1485,0.3219,-0.1185,-0.2839,0.2488,1.0000,0.0940,-0.1194,0.1852,-0.1272,-0.1893,0.0399,-0.0084],
				[0.2335,-0.1916,-0.2037,-0.2357,0.0025,-0.0189,0.0796,-0.0956,0.1068,-0.1201,0.0309,0.1018,0.0940,1.0000,0.1020,0.0104,-0.1762,-0.2871,0.1063,0.0213],
				[0.0807,-0.0774,-0.2762,0.0109,-0.3311,0.0006,0.0092,0.0881,0.0268,0.1164,-0.1498,-0.0556,-0.1194,0.1020,1.0000,0.0257,-0.1306,0.1000,-0.0951,-0.2423],
				[-0.0024,0.0735,-0.0429,-0.1681,-0.1409,0.0678,0.0328,-0.0646,0.0781,0.0458,0.3240,-0.0657,0.1852,0.0104,0.0257,1.0000,-0.2438,-0.2582,-0.2274,0.0310],
				[-0.1245,0.0239,-0.1202,-0.0752,-0.0738,0.1373,-0.1051,0.1750,0.0529,-0.0040,-0.0849,0.0473,-0.1272,-0.1762,-0.1306,-0.2438,1.0000,-0.0231,0.1164,0.0304],
				[0.0835,0.0936,0.1024,0.2626,-0.0535,0.3409,0.0139,-0.2343,0.0346,-0.0958,-0.1001,-0.1634,-0.1893,-0.2871,0.1000,-0.2582,-0.0231,1.0000,0.0632,0.0336],
				[0.1783,0.0030,0.1826,0.0212,-0.0684,-0.0739,-0.0700,0.0476,0.1081,-0.2451,0.0279,-0.1715,0.0399,0.1063,-0.0951,-0.2274,0.1164,0.0632,1.0000,-0.3556],
				[-0.0989,-0.2237,-0.2397,0.0834,0.2921,0.1559,0.1420,-0.0378,0.1386,0.0366,0.1015,0.0021,-0.0084,0.0213,-0.2423,0.0310,0.0304,0.0336,-0.3556,1.0000]]
	var nbAssets = sigma.length;
				
	// Compute ERB weights
	var weights = PortfolioAllocation.equalRiskBoundingWeights(sigma);				
	}
*/
});


QUnit.test('Random weights portfolio', function(assert) {    
	// Test with random data, no constraints
	{
		// Setup static parameters of the random test
		var nbTests = 50;
		var nbAssetsMin = 1;
		var nbAssetsMax = 50;

		// Aim of these tests is to check that for a portfolio of n assets:
		// - The number of weights computed is n
		// - The weights belong to the interval [0, 1]
		// - The weights sum to 1	  
		for (var i = 0; i < nbTests; ++i) {
			// Generate a random number of assets
			var nbAssets = Math.floor(Math.random()*(nbAssetsMax - nbAssetsMin + 1) + nbAssetsMin);
			
			// Generate a random portfolio for this number of assets
			var randomWeights = PortfolioAllocation.randomWeights(nbAssets);
			
			// Check that the number of weights corresponds to the number of assets
			var nbWeights = randomWeights.length;
			assert.equal(nbWeights, nbAssets, "Random weights portfolio, number of weights - Test " + i);
			
			// Check that the weights belong to the unit interval
			var weightsBelongToUnitInterval = true;
			for (var k = 0; k < randomWeights.length; ++k) {
				if (randomWeights[k] > 1 || randomWeights[k] < 0) {
					weightsBelongToUnitInterval = false;
					break;
				}
			}
			assert.equal(weightsBelongToUnitInterval, true, "Random weights portfolio, weights in unit interval - Test " + i);

			// Check that the sum of the weights is 1, near to machine precision
			var weightsSum = 0;
			for (var k = 0; k < randomWeights.length; ++k) {
				weightsSum += randomWeights[k];
			}
			assert.equal(Math.abs(weightsSum - 1) <= 1e-15, true, "Random weights portfolio, weights sum to one - Test " + i);
		}
	}

	// Test with random data, cardinality constraints
	{
		// Setup static parameters of the random test
		var nbTests = 50;
		var nbAssetsMin = 1;
		var nbAssetsMax = 50;

		// Aim of these tests is to check that for a portfolio of n assets constrained to hold betwen i and j assets:
		// - The number of non-zero weights is between i and j
		for (var i = 0; i < nbTests; ++i) {
			// Generate a random number of assets
			var nbAssets = Math.floor(Math.random()*(nbAssetsMax - nbAssetsMin + 1) + nbAssetsMin);
			
			// Generate random cardinality constraints
			var maxAssets = Math.floor(Math.random()*(nbAssets - 1 + 1) + 1);
			var minAssets = Math.floor(Math.random()*(maxAssets - 1 + 1) + 1);
			
			// Generate a random portfolio for this number of assets and these cardinality constraints
			var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { minAssets: minAssets, maxAssets: maxAssets } });
			
			// Check that the number of non-zero weights is between minAssets and maxAssets
			var nbNonZeroAssets = 0;
			for (var k = 0; k < randomWeights.length; ++k) {
				if (randomWeights[k] != 0) {
					nbNonZeroAssets = nbNonZeroAssets + 1;
				}
			}
			assert.equal(maxAssets >= nbNonZeroAssets, true, "Random weights portfolio with cardinality constraints, number of non-zero weights - Test " + i + "/1");
			assert.equal(minAssets <= nbNonZeroAssets, true, "Random weights portfolio with cardinality constraints, number of non-zero weights - Test " + i + "/2");
		}
	}
});


QUnit.test('Minimax portfolio', function(assert) {    
	// Static data
	// Reference: Christos Papahristodoulou, Optimal portfolios using Linear Programming models (preliminary version on Internet)
	// (b) Maximin formulation example
	{
		// Define the assets returns
		var assetsReturns = [[0.054,  0.045,  -0.03, -0.018,  0.043,  0.047,  0.055,  0.036, -0.039, -0.043,  0.046, 0.052], // Asset 1 returns
		[0.032,  0.055, -0.036,  0.052,  0.047,  0.034,  0.063,  0.048,  0.025,   0.04,  0.036, -0.017], // ...
		[0.064,  0.056,  0.048,  0.007,  0.053,  0.036,  0.017,  0.047, -0.059,  0.047,   0.04, 0.032],
		[0.038,  0.062, -0.037,   0.05,  0.065, -0.043,  0.062,  0.034,  0.035,  0.056,  0.057, 0.025],
		[0.049,  0.067, -0.039,  0.051,  0.049,  0.037,  0.055,  0.025,  0.052,  0.02,   0.045,  0.04]];
		
		// Compute the associated minimax portfolio weights
		var minimaxWeights = PortfolioAllocation.minimaxWeights(assetsReturns);
		
		// Compare minimax weights to expected weights
		var expectedWeights = [0, 0, 0.459596, 0, 0.540404];
		for (var i = 0; i < assetsReturns.length; ++i) { 
			assert.equal(Math.abs(minimaxWeights[i] - expectedWeights[i]) <= 1e-3, true, 'Minimax - Values ' + i);
		}
		
		// Compute the minimum portfolio return
		var minimaxWeights = PortfolioAllocation.minimaxWeights(assetsReturns, {outputMinimumPortfolioReturn: true});
		var expectedMinPortfolioReturn = 0.000985;
		assert.equal(Math.abs(minimaxWeights[1] - expectedMinPortfolioReturn) <= 1e-3, true, 'Minimax - Minimum Portfolio Return');
	} 
	
	// Static data, no full investment
	// Reference: Christos Papahristodoulou, Optimal portfolios using Linear Programming models (preliminary version on Internet)
	// (b) Maximin formulation example, with March return of security C changed from 0.048 to -0.03
	{
		// Define the assets returns
		var assetsReturns = [[0.054,  0.045,  -0.03, -0.018,  0.043,  0.047,  0.055,  0.036, -0.039, -0.043,  0.046, 0.052], // Asset 1 returns
		[0.032,  0.055, -0.036,  0.052,  0.047,  0.034,  0.063,  0.048,  0.025,   0.04,  0.036, -0.017], // ...
		[0.064,  0.056,  -0.03,  0.007,  0.053,  0.036,  0.017,  0.047, -0.059,  0.047,   0.04, 0.032], 
		[0.038,  0.062, -0.037,   0.05,  0.065, -0.043,  0.062,  0.034,  0.035,  0.056,  0.057, 0.025],
		[0.049,  0.067, -0.039,  0.051,  0.049,  0.037,  0.055,  0.025,  0.052,  0.02,   0.045,  0.04]];
		
		// Compute the associated minimax portfolio weights
		var minimaxWeights = PortfolioAllocation.minimaxWeights(assetsReturns, {constraints: {fullInvestment: false}});
		
		// Compare minimax weights to expected weights
		var expectedWeights = [0, 0, 0, 0, 0];
		for (var i = 0; i < assetsReturns.length; ++i) { 
			assert.equal(Math.abs(minimaxWeights[i] - expectedWeights[i]) <= 0, true, 'Minimax #2 - Values ' + i);
		}
		
		// Compute the minimum portfolio return
		var minimaxWeights = PortfolioAllocation.minimaxWeights(assetsReturns, {outputMinimumPortfolioReturn: true});
		var expectedMinPortfolioReturn =  -0.030794545037876233;
		assert.equal(Math.abs(minimaxWeights[1] - expectedMinPortfolioReturn) <= 1e-3, true, 'Minimax #2 - Minimum Portfolio Return');
	} 
});


QUnit.test('Mean variance portfolio - corner portfolios computation', function(assert) {    
	// Test using static data
	// Test limit cases for determining the critical line algorithm starting portfolio:
	{
		// Lower bounds binding (sum lb_i == 1)
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios([0.1, 0.2], [[1,0],[0,1]], { constraints: {minWeights: [0.4, 0.6]} });
		assert.deepEqual(cornerPortfolios, [[[0.4, 0.6],  0.16, 0.7211102550927979]], 'Mean variance portfolio - Corner portfolios, lower bounds binding');
		
		// Upper bounds binding (sum ub_i == 1)
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios([0.1, 0.2], [[1,0],[0,1]], { constraints: {maxWeights: [0.6, 0.4]} });
		assert.deepEqual(cornerPortfolios, [[[0.6, 0.4],  0.14, 0.7211102550927979]], 'Mean variance portfolio - Corner portfolios, upper bounds binding');
	}

	// Test using static data
	// Reference: https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm
	// Note: this test also allows checking that numerically equal corner portfolios are filtered out by the algorithm
	{
		var expectedCornerPortfolios = [[[0.2, 0.30000000000000004, 0.5], 7.8500000000000005, 8.777323054325846], 
										[[0.2, 0.5, 0.30000000000000004], 6.950000000000001, 6.92166165021088], 
										[[0.22180737780348653, 0.5, 0.27819262219651353], 6.775540977572108, 6.6431091200292105],
										[[0.451915610952186, 0.348084389047814, 0.2], 5.618295361667349, 4.81952189552386],
										[[0.5, 0.2999999999999999, 0.2], 5.449999999999999, 4.560824486866382]].reverse();
		
		var covMat = [[1, 2.96, 2.31],
					[2.96, 54.76, 39.886],
					[2.31,	39.886,	237.16]];
		var returns = [2.8000, 6.3000, 10.8000];
		
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat, { constraints: {minWeights: [0.2, 0.2, 0.2], maxWeights: [0.5, 0.5, 0.5]} });
		assert.deepEqual(cornerPortfolios, expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #1');
	}
	
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter VIII "The computing procedure"
	{
		var expectedCornerPortfolios = [[[0, 1, 0], 0.146, 0.29223278392404917], 
										[[0, 0.22496808316614988, 0.7750319168338501], 0.1320494254969907, 0.15908574925905583], 
										[[0.8414051841746248, 0, 0.15859481582537516], 0.07246725784447476, 0.12220061216349064], 
										[[0.9931034482758623, 0, 0.006896551724137813], 0.0624551724137931, 0.12082760588883482]].reverse();
		
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat);
		assert.deepEqual(cornerPortfolios, expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #2');
	}
	
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter II "Illustrative portfolio analyses"
	{
		var expectedCornerPortfolios = [[[0, 0, 0, 0, 100, 0, 0, 0, 0, 0], 0.1981111111, 0.35761725182099363], 
										[[0, 0, 0, 0, 96, 0, 0, 4, 0, 0], 0.19776400166660185, 0.35451400571130365], 
										[[0, 0, 0, 8, 92, 0, 0, 0, 0, 0], 0.19602068236429498, 0.3403126771859089], 
										[[0, 0, 0, 42, 58, 0, 0, 0, 0, 0], 0.18785503614831273, 0.29187554881080124], // In Markowitz book, this is actually [0, 0, 0, 41, 59, 0, 0, 0, 0, 0], but this seems a wrong rounding 
										[[0, 0, 0, 25, 32, 0, 43, 0, 0, 0], 0.16161551062944268, 0.2069944757322166],
										[[0, 0, 16, 3, 12, 0, 69, 0, 0, 0], 0.14035962162442164, 0.16472622029650955],
										[[0, 0, 0, 0, 0, 0, 0, 0, 0, 100], 0, 0]].reverse();
	
		var covMat = [[0.05338816358,0.02149069753,0.02865533642,0.04896485802,0.01624895062,0.03223945062,0.02425553395,0.03999812963,0.0361509784, 0],
					[0.02149069753,0.01468446914,0.01878391358,0.02441658642,0.008041938272,0.01002193827,0.01448993827,0.02536259259,0.02083593827, 0],
					[0.02865533642,0.01878391358,0.08550016358,0.06260714198,0.04439938272,0.01328671605,0.01043991049,0.06864603704,0.0420215216, 0],
					[0.04896485802,0.02441658642,0.06260714198,0.09546446914,0.05153806173,0.02902461728,0.02077028395,0.09002012963,0.03664589506, 0],
					[0.01624895062,0.008041938272,0.04439938272,0.05153806173,0.1278900988,0.0128384321,0.02091715432,0.1015344074,0.04497232099, 0],
					[0.03223945062,0.01002193827,0.01328671605,0.02902461728,0.0128384321,0.04125832099,0.01127854321,0.02960762963,0.02165332099, 0],
					[0.02425553395,0.01448993827,0.01043991049,0.02077028395,0.02091715432,0.01127854321,0.02883379321,0.02913762963,0.01739445988, 0],
					[0.03999812963,0.02536259259,0.06864603704,0.09002012963,0.1015344074,0.02960762963,0.02913762963,0.1467278889,0.05284057407, 0],
					[0.0361509784,0.02083593827,0.0420215216,0.03664589506,0.04497232099,0.02165332099,0.01739445988,0.05284057407,0.07926979321, 0],
					[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
		var returns = [0.06594444444,0.06155555556,0.1460555556,0.1734444444,0.1981111111,0.05511111111,0.1276111111,0.1903333333,0.1156111111, 0];

		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat);
		
		for (var i = 0; i < expectedCornerPortfolios.length; ++i) {
			var expectedWeights = expectedCornerPortfolios[i][0];
			var expectedReturn = expectedCornerPortfolios[i][1];
			var expectedVolatility = expectedCornerPortfolios[i][2];
			
			var computedWeights = cornerPortfolios[i][0];
			var roundedComputedWeights = PortfolioAllocation.roundedWeights(computedWeights, 100);
			var computedReturn = cornerPortfolios[i][1];
			var computedVolatility = cornerPortfolios[i][2];
			
			assert.equal(Math.abs(computedReturn - expectedReturn) <= 1e-14, true, 'Mean variance portfolio - Corner portfolios #3, #' + i + ", return");
			assert.equal(Math.abs(computedVolatility - expectedVolatility) <= 1e-14, true, 'Mean variance portfolio - Corner portfolios #3, #' + i + ", volatility");
			
			for (var j = 0; j < expectedWeights.length; ++j) {
				assert.equal(expectedWeights[j]/100, roundedComputedWeights[j], 'Mean variance portfolio - Corner portfolios #3, #' + i + ",  " + j + "/" + expectedWeights.length);
			}
		}
	}
	
	// Test using static data
	// Reference: A Simple Spreadsheet-Based Exposition of the Markowitz Critical Line Method for Portfolio Selection, Clarence C. Kwan
	{
		var expectedCornerPortfolios = [[[0, 0, 1], 0.12, 0.1], 
										[[0, 0.6485013623978204, 0.3514986376021796], 0.0940599455040872, 0.05237167799176803], 
										[[0.9754098360655736, 0, 0.024590163934426246], 0.051721311475409835, 0.019905041975576163], 
										[[0.9799999999999999, 0, 0.02000000000000001], 0.051399999999999994, 0.019899748742132396]].reverse();
		
		var covMat = [[0.0004, 0.0004, 0.0002],
					 [0.0004, 0.0025,0.001],
					  [0.0002, 0.001, 0.01]];
		var returns = [0.05, 0.08, 0.12];
		
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat);
		assert.deepEqual(cornerPortfolios, expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #4');
	}
	
	// Test using static data (upper bounds)
	// Reference: A Simple Spreadsheet-Based Exposition of the Markowitz Critical Line Method for Portfolio Selection, Clarence C. Kwan
	{
		var expectedCornerPortfolios = [[[0, 0.30000000000000004, 0.7], 0.108, 0.07446475676452585], 
										[[0, 0.6485013623978203, 0.3514986376021798], 0.0940599455040872, 0.05237167799176804], 
										[[0.7, 0.18310626702997274, 0.11689373297002724], 0.06367574931880109, 0.02438316869432441], 
										[[0.7, 0.2438095238095238, 0.05619047619047619], 0.061247619047619044, 0.02357642082775965]].reverse();
		
		var covMat = [[0.0004, 0.0004, 0.0002],
					 [0.0004, 0.0025,0.001],
					  [0.0002, 0.001, 0.01]];
		var returns = [0.05, 0.08, 0.12];
		
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat, { constraints: {maxWeights: [0.7, 0.7, 0.7]} });
		assert.deepEqual(cornerPortfolios, expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #5');
	}
	
	// Test using static data
	// Reference: An Open-Source Implementation of the Critical-Line Algorithm for Portfolio Optimization, David H. Bailey and Marcos Lopez de Prado
	{
		var expectedCornerPortfolios = [[[0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 1.19, 0.9520003676469878], 
										[[0.6493694070931811, 0.3506305929068189, 0, 0, 0, 0, 0, 0, 0, 0], 1.1802594588936024, 0.5456568742682376], 
										[[0.4339841341086239, 0.23124750065448754, 0, 0.334768365236889, 0, 0, 0, 0, 0, 0], 1.1600564524217891, 0.4172556503733193],
										[[0.12688785385570883, 0.07234334721032556, 0, 0.28125374926334057, 0, 0, 0, 0, 0, 0.5195150496706249], 1.1112622642799619, 0.2667196142113227],
										[[0.12320100405906734, 0.07044407130753655, 0, 0.2789935668090118, 0, 0, 0, 0.006435564362887149, 0, 0.5209257934614971], 1.1083602383747904, 0.26501699530147116],
										[[0.0869215492990579, 0.050451042268558385, 0, 0.22359401742288823, 0, 0.17383161507156486, 0, 0.03017301555135618, 0, 0.4350287603865743], 1.0224838894431951, 0.2296800735960956],
										[[0.0846709411996219, 0.049253858741118525, 0, 0.21963390336360733, 0, 0.18003923464176064, 0, 0.03102980185535347, 0.006485702415438152, 0.42888655778310003], 1.0153059205224644, 0.22798274842085042],
										[[0.07378925302280315, 0.043828660769718863, 0, 0.19897560805881487, 0.026158159857441972, 0.19815187227970524, 0, 0.03341958639919798, 0.027902966026643668, 0.3977738935856743], 0.9727204340349849, 0.21955488008267854],
										[[0.06834400480527462, 0.041387026820649334, 0.015215259551836627, 0.18813443107045838, 0.03416248599274816, 0.20231943214747125, 0, 0.0339293235595669, 0.03363264959172938, 0.38287538646026537], 0.9499368157550272, 0.21602457168689138],
										[[0.03696858147921504, 0.02690083780081047, 0.0949424305647986, 0.1257759521946726, 0.0767460810325476, 0.21935567131616898, 0.029987096882220312, 0.035963284621386274, 0.06134983772972688, 0.29201022637845325], 0.8032153598765337, 0.20523761981386526]].reverse();
										
		var covMat = [[0.40755159,0.03175842,0.05183923,0.05663904,0.0330226,0.00827775,0.02165938,0.01332419,0.0343476,0.02249903],
					[0.03175842,0.9063047,0.03136385,0.02687256,0.01917172,0.00934384,0.02495043,0.00761036,0.02874874,0.01336866],
					[0.05183923,0.03136385,0.19490901,0.04408485,0.03006772,0.01322738,0.03525971,0.0115493,0.0427563,0.02057303],
					[0.05663904,0.02687256,0.04408485,0.19528471,0.02777345,0.00526665,0.01375808,0.00780878,0.02914176,0.01640377],
					[0.0330226,0.01917172,0.03006772,0.02777345,0.34059105,0.00777055,0.02067844,0.00736409,0.02542657,0.01284075],
					[0.00827775,0.00934384,0.01322738,0.00526665,0.00777055,0.15983874,0.02105575,0.00518686,0.01723737,0.00723779],
					[0.02165938,0.02495043,0.03525971,0.01375808,0.02067844,0.02105575,0.68056711,0.01377882,0.04627027,0.01926088],
					[0.01332419,0.00761036,0.0115493,0.00780878,0.00736409,0.00518686,0.01377882,0.95526918,0.0106553,0.00760955],
					[0.0343476,0.02874874,0.0427563,0.02914176,0.02542657,0.01723737,0.04627027,0.0106553,0.31681584,0.01854318],
					[0.02249903,0.01336866,0.02057303,0.01640377,0.01284075,0.00723779,0.01926088,0.00760955,0.01854318,0.11079287]];
		var returns = [1.175,1.19,0.396,1.12,0.346,0.679,0.089,0.73,0.481,1.08];

		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat);
		assert.deepEqual(cornerPortfolios, expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #6');
	}
});



QUnit.test('Mean variance portfolio - target return weights portfolio', function(assert) {    
	// Test using random data
	// Test unreachable cases
	{
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		var maxReturn = returns[1];
		var unreachableMaxReturn = maxReturn + (1 - Math.random()); // > maxReturn
		assert.throws(function() { 
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: unreachableMaxReturn}}) },
			new Error('target return not reachable'),
			"Mean variance portfolio - Target return weights portfolio, unreachable target return #1");
			
		var minReturn = returns[0];
		var unreachableMinReturn = minReturn - (1 - Math.random()); // < minReturn
		assert.throws(function() { 
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: unreachableMinReturn}}) },
			new Error('target return not reachable'),
			"Mean variance portfolio - Target return weights portfolio, unreachable target return #2");
	}
	
	// Test using static data
	{
		var covMat = [[1, 0.3], [0.3, 1]];
		var returns = [0.1, 0.2];
		
		var targetReturn = returns[1]; // the second asset
		var expectedWeights = [0, 1]; 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: targetReturn}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - target return weights portfolio #0');
		
		var targetReturn = 0.15; // exact mix of the two assets
		var expectedWeights = [0.5, 0.5];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: targetReturn}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - target return weights portfolio #0 bis');
	}
		
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter II "Illustrative portfolio analyses"
	{	
		var covMat = [[0.05338816358,0.02149069753,0.02865533642,0.04896485802,0.01624895062,0.03223945062,0.02425553395,0.03999812963,0.0361509784, 0],
					[0.02149069753,0.01468446914,0.01878391358,0.02441658642,0.008041938272,0.01002193827,0.01448993827,0.02536259259,0.02083593827, 0],
					[0.02865533642,0.01878391358,0.08550016358,0.06260714198,0.04439938272,0.01328671605,0.01043991049,0.06864603704,0.0420215216, 0],
					[0.04896485802,0.02441658642,0.06260714198,0.09546446914,0.05153806173,0.02902461728,0.02077028395,0.09002012963,0.03664589506, 0],
					[0.01624895062,0.008041938272,0.04439938272,0.05153806173,0.1278900988,0.0128384321,0.02091715432,0.1015344074,0.04497232099, 0],
					[0.03223945062,0.01002193827,0.01328671605,0.02902461728,0.0128384321,0.04125832099,0.01127854321,0.02960762963,0.02165332099, 0],
					[0.02425553395,0.01448993827,0.01043991049,0.02077028395,0.02091715432,0.01127854321,0.02883379321,0.02913762963,0.01739445988, 0],
					[0.03999812963,0.02536259259,0.06864603704,0.09002012963,0.1015344074,0.02960762963,0.02913762963,0.1467278889,0.05284057407, 0],
					[0.0361509784,0.02083593827,0.0420215216,0.03664589506,0.04497232099,0.02165332099,0.01739445988,0.05284057407,0.07926979321, 0],
					[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
		var returns = [0.06594444444,0.06155555556,0.1460555556,0.1734444444,0.1981111111,0.05511111111,0.1276111111,0.1903333333,0.1156111111, 0];
		
		var expectedRoundedWeights = [0, 0, 0.08, 0.015, 0.06, 0, 0.345, 0, 0, 0.5];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: 0.07}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #1');
		
		var expectedRoundedWeights = [0, 0, 0, 0.335, 0.455, 0, 0.21, 0, 0, 0]; // In Markowitz book, this is actually [0, 0, 0, 0.33, 0.455, 0, 0.215, 0, 0, 0], but this seems a wrong rounding 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: 0.175}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #2');
		
		var expectedRoundedWeights = [0, 0, 0, 0.3, 0.39, 0, 0.31, 0, 0, 0]; // In Markowitz book, this is actually [0, 0, 0, 0.29, 0.39, 0, 0.32, 0, 0, 0], but this seems a wrong rounding 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: 0.169}});
		weights = PortfolioAllocation.roundedWeights(weights, 100);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #3');
		
		var expectedRoundedWeights = [0, 0, 0.09, 0.13, 0.21, 0, 0.57, 0, 0, 0];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetReturn: 0.150}});
		weights = PortfolioAllocation.roundedWeights(weights, 100);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #4');
	}
});	


QUnit.test('Mean variance portfolio - target volatility weights portfolio', function(assert) {    
	function generateRandomValue(minVal, maxVal) {	
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test using static data
	// Test unreachable cases
	{
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		
		var nbSamples = 10;
		for (var k = 0; k < nbSamples; ++k) {
			var maxVolatility = Math.sqrt(covMat[1][1]);
			var unreachableMaxVolatility = generateRandomValue(maxVolatility - 1e-6, 2 * maxVolatility); // > maxVolatility
			assert.throws(function() { 
			PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetVolatility: unreachableMaxVolatility}}) },
				new Error('target volatility not reachable'),
				"Mean variance portfolio - Target volatility weights portfolio, unreachable target volatility #1 - " + k + "/" + nbSamples);
		
			var minVolatility = 0.12082760588883482; // computed thanks to the efficient frontier
			var unreachableMinVolatility = generateRandomValue(0, minVolatility - 1e-6); // < minVolatility
			assert.throws(function() { 
			PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetVolatility: unreachableMinVolatility}}) },
				new Error('target volatility not reachable'),
				"Mean variance portfolio - Target volatility weights portfolio, unreachable target volatility #2 - " + k + "/" + nbSamples);
		}
	}
	
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter II "Illustrative portfolio analyses"
	{	
		var covMat = [[0.05338816358,0.02149069753,0.02865533642,0.04896485802,0.01624895062,0.03223945062,0.02425553395,0.03999812963,0.0361509784, 0],
					[0.02149069753,0.01468446914,0.01878391358,0.02441658642,0.008041938272,0.01002193827,0.01448993827,0.02536259259,0.02083593827, 0],
					[0.02865533642,0.01878391358,0.08550016358,0.06260714198,0.04439938272,0.01328671605,0.01043991049,0.06864603704,0.0420215216, 0],
					[0.04896485802,0.02441658642,0.06260714198,0.09546446914,0.05153806173,0.02902461728,0.02077028395,0.09002012963,0.03664589506, 0],
					[0.01624895062,0.008041938272,0.04439938272,0.05153806173,0.1278900988,0.0128384321,0.02091715432,0.1015344074,0.04497232099, 0],
					[0.03223945062,0.01002193827,0.01328671605,0.02902461728,0.0128384321,0.04125832099,0.01127854321,0.02960762963,0.02165332099, 0],
					[0.02425553395,0.01448993827,0.01043991049,0.02077028395,0.02091715432,0.01127854321,0.02883379321,0.02913762963,0.01739445988, 0],
					[0.03999812963,0.02536259259,0.06864603704,0.09002012963,0.1015344074,0.02960762963,0.02913762963,0.1467278889,0.05284057407, 0],
					[0.0361509784,0.02083593827,0.0420215216,0.03664589506,0.04497232099,0.02165332099,0.01739445988,0.05284057407,0.07926979321, 0],
					[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
		var returns = [0.06594444444,0.06155555556,0.1460555556,0.1734444444,0.1981111111,0.05511111111,0.1276111111,0.1903333333,0.1156111111, 0];
		
		var expectedRoundedWeights = [0, 0, 0.08, 0.015, 0.055, 0, 0.335, 0, 0, 0.515]; // In Markowitz book, this is actually [0, 0, 0.08, 0.015, 0.06, 0, 0.345, 0, 0, 0.5], but this seems a wrong rounding 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { constraints: {targetVolatility: 0.08}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target volatility weights portfolio #1');
	}
});	


QUnit.test('Mean variance portfolio - efficient portfolios computation', function(assert) {    
	// Test using static data
	// Test the limit case of only one corner portfolio
	{
		// Lower bounds binding (sum lb_i == 1)
		assert.throws(function() { 
			PortfolioAllocation.meanVarianceEfficientFrontierPortfolios([0.1, 0.2], [[1,0],[0,1]], { constraints: {minWeights: [0.4, 0.6]} }) },
			new Error('efficient frontier made of only one corner portfolio: only one efficient portfolio can be computed'),
			"Mean variance portfolio - Efficient portfolios, lower bounds binding KO");
			
		var efficientPortfolios = PortfolioAllocation.meanVarianceEfficientFrontierPortfolios([0.1, 0.2], [[1,0],[0,1]], { nbPortfolios: 1, constraints: {minWeights: [0.4, 0.6]} });
		assert.deepEqual(efficientPortfolios, [[[0.4, 0.6],  0.16, 0.7211102550927979]], 'Mean variance portfolio - Efficient portfolios, lower bounds binding OK');
	}

	// Test using static data
	// Test the limit case of only one efficient portfolio to be computed
	// Reference: Portfolio Selection, H. Markowitz example, chapter VIII "The computing procedure"
	{
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		assert.throws(function() { 
			PortfolioAllocation.meanVarianceEfficientFrontierPortfolios(returns, covMat, { nbPortfolios: 1 }) },
			new Error('efficient frontier made of several corner portfolios: at least two efficient portfolios must be computed'),
			"Mean variance portfolio - Efficient portfolios, one efficient portfolio KO");
	}
	
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter VIII "The computing procedure"
	{		
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Test that by default, 100 efficient portfolios are computed
		var efficientPortfolios = PortfolioAllocation.meanVarianceEfficientFrontierPortfolios(returns, covMat);
		assert.equal(efficientPortfolios.length, 100, 'Mean variance portfolio - Efficient portfolios, default number of efficient portfolios computed');
		
		// Test a small number of efficient portfolios
		var expectedEfficientPortfolios = [[[0.9931034482758623, 0, 0.006896551724137813], 0.0624551724137931, 0.12082760588883482], 
											 [[0, 0.39839572192513373, 0.6016042780748663], 0.1351711229946524, 0.1702926860745384], 
											 [[0, 0.5989304812834225, 0.4010695187165775], 0.13878074866310158, 0.20069797992103716], 
											 [[0, 0.7994652406417113, 0.20053475935828874], 0.1423903743315508, 0.24306339262469084], 
											 [[0, 1, 0], 0.146, 0.29223278392404917]];
										
		var efficientPortfolios = PortfolioAllocation.meanVarianceEfficientFrontierPortfolios(returns, covMat, { nbPortfolios: 5 });
		assert.deepEqual(efficientPortfolios, expectedEfficientPortfolios, 'Mean variance portfolio - Efficient portfolios #1');
	}
});


QUnit.test('Mean variance portfolio - maximum Sharpe ratio computation', function(assert) {    
	// Test using static data
	// Reference: An Open-Source Implementation of the Critical-Line Algorithm for Portfolio Optimization, David H. Bailey and Marcos Lopez de Prado
	{
		// Note: this portfolio has a Sharpe ratio of ~4.4535(3...), with a null risk free rate, which
		// is exactly the same value as in the reference.
		var expectedWeights = [0.08397318948217423, 0.04890598613711377, 0, 0.21830925954049477, 0.0016773041709357821,
							   0.18120064671441183, 0, 0.031183038765169507, 0.007859012532850177, 0.42689156265685];
										
		var covMat = [[0.40755159,0.03175842,0.05183923,0.05663904,0.0330226,0.00827775,0.02165938,0.01332419,0.0343476,0.02249903],
					[0.03175842,0.9063047,0.03136385,0.02687256,0.01917172,0.00934384,0.02495043,0.00761036,0.02874874,0.01336866],
					[0.05183923,0.03136385,0.19490901,0.04408485,0.03006772,0.01322738,0.03525971,0.0115493,0.0427563,0.02057303],
					[0.05663904,0.02687256,0.04408485,0.19528471,0.02777345,0.00526665,0.01375808,0.00780878,0.02914176,0.01640377],
					[0.0330226,0.01917172,0.03006772,0.02777345,0.34059105,0.00777055,0.02067844,0.00736409,0.02542657,0.01284075],
					[0.00827775,0.00934384,0.01322738,0.00526665,0.00777055,0.15983874,0.02105575,0.00518686,0.01723737,0.00723779],
					[0.02165938,0.02495043,0.03525971,0.01375808,0.02067844,0.02105575,0.68056711,0.01377882,0.04627027,0.01926088],
					[0.01332419,0.00761036,0.0115493,0.00780878,0.00736409,0.00518686,0.01377882,0.95526918,0.0106553,0.00760955],
					[0.0343476,0.02874874,0.0427563,0.02914176,0.02542657,0.01723737,0.04627027,0.0106553,0.31681584,0.01854318],
					[0.02249903,0.01336866,0.02057303,0.01640377,0.01284075,0.00723779,0.01926088,0.00760955,0.01854318,0.11079287]];
		var returns = [1.175,1.19,0.396,1.12,0.346,0.679,0.089,0.73,0.481,1.08];
		var rf = 0;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf);
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #1');
	}
	
	// Test using static data
	// Reference: Optimal Portfolio Selection Without Short Sales Under the Full-Information Covariance Structure: A Pedagogic Consideration, Clarence C. Y. Kwan and Yufei Yuan
	{
		var expectedRoundedWeights = [0.3318, 0, 0.2443, 0, 0.4239, 0];

		// Note: the covariance matrix below is the result of Diag(STD)*CORR*Diag(STD), because
		// the covariance matrix as displayed in the reference has some digits truncated.
		var covMat = [[0.010, 0.010, -0.0018, 0.0024, 0.0016, 0.0048],
		              [0.010, 0.0625, 0.0135, 0.009, 0.002, 0.008],
					  [-0.0018, 0.0135, 0.0324, 0.00432, -0.00288, 0.00864],
					  [0.0024, 0.009, 0.00432, 0.0144, 0.0096, 0.00192],
					  [0.0016, 0.002, -0.00288, 0.0096, 0.0064, 0.00256],
					  [0.0048, 0.008, 0.00864, 0.00192, 0.00256, 0.0256]];
		var returns = [0.15, 0.18, 0.20, 0.11, 0.13, 0.12];
		var rf = 0.08;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf);
		maxSharpeRatioPortfolioWeights = PortfolioAllocation.roundedWeights(maxSharpeRatioPortfolioWeights, 10000);
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedRoundedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #2');
	}
	
	// Test using static data
	// Reference: Maximizing the Sharpe Ratio and Information Ratio in the Barra Optimizer, Leonid Kopman, Scott Liu
	{
		var expectedRoundedWeights = [0.2727, 0.7273];
		
		var covMat = [[0.05, 0.02],
		              [0.02, 0.07]];
		var returns = [0.05, 0.1];
		var rf = 0;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf);
		maxSharpeRatioPortfolioWeights = PortfolioAllocation.roundedWeights(maxSharpeRatioPortfolioWeights, 10000);
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedRoundedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #3');
	}

	// Test using static data
	// Reference: Optimal mean-variance portfolio selection using Cauchy–Schwarz maximization, Hsin-Hung Chen, Hsien-Tang Tsai, Dennis K. J. Lin
	{
		// Note: the weights in the reference do not sum to 1 (they sum to 0.999)
		//
		// Note 2: the weights in the reference do not reach the maximum Sharpe ratio, 
		// because they reach ~0.743369... compared to ~0.743391... for the computed weights
		// below, plus there is a mismatch between the computed Sharpe ratio and the 
		// Sharpe ratio provided in the reference.
		//
		// Most probably, this is due to truncated numerical data displayed in the reference
		// v.s. non truncated data used by the researchers in their algorithm.
		
		//var referenceRoundedWeights = [0, 0.1355, 0, 0.5153, 0.2292, 0.119];
		var expectedRoundedWeights = [0, 0.1355, 0, 0.5153, 0.2290, 0.1202];
		
		/* The covariance matrix below has been obtained through Diag(STD)*CORR*Diag(STD)
		var stdDev = [0.2929, 0.3252, 0.3278, 0.2555, 0.3943, 0.5698];
		var corrMat = [[1.0000, 0.1484, 0.3076, 0.3685, 0.4482, 0.3802],
		              [0.1484, 1.0000, 0.0382, -0.1451, 0.2714, 0.1025],
					  [0.3076, 0.0382, 1.0000, 0.2990, 0.3634, 0.4388],
					  [0.3685, -0.1451, 0.2990, 1.0000, 0.2756, 0.2796],
					  [0.4482, 0.2714, 0.3634, 0.2756, 1.0000, 0.7542],
					  [0.3802, 0.1025, 0.4388, 0.2796, 0.7542, 1.0000]];
		*/
		var covMat = [[0.08579041, 0.014135260272, 0.029533481911999993, 0.027577047575, 0.05176282865399999, 0.06345325848399999 ],
						[0.014135260271999999, 0.10575504, 0.004072141391999999, -0.01205615586, 0.034800634104, 0.018993143399999995 ],
						[0.029533481911999997, 0.004072141392, 0.10745284, 0.0250421171, 0.046970009635999996, 0.081959257072 ],
						[0.027577047575000003, -0.01205615586, 0.0250421171, 0.06528025, 0.02776494994, 0.04070525844 ],
						[0.051762828654, 0.03480063410399999, 0.046970009635999996, 0.027764949940000002, 0.15547249, 0.16944772798799998 ],
						[0.06345325848399999, 0.018993143399999995, 0.081959257072, 0.04070525844, 0.16944772798799998, 0.32467204 ]];
		var returns = [0.0912, 0.1087, 0.1283, 0.1882, 0.2896, 0.3739];
		var rf = 0.05;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf);
		maxSharpeRatioPortfolioWeights = PortfolioAllocation.roundedWeights(maxSharpeRatioPortfolioWeights, 10000);
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedRoundedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #4');
	}
	
	// Test using static data
	// Test the case of a positive semi-definite covariance matrix with a 0 volatility point, as well as a 0 return point, on the efficient frontier
	// The result has been validated with a grid search
	// Reference: Portfolio Selection, H. Markowitz example, chapter II "Illustrative portfolio analyses"
	{		
		var expectedWeights = [0, 0, 0.16171761768982862, 0.0316824768157338, 0.1179234805759474, 0, 0.6886764249184905, 0, 0, 0];
		
		var covMat = [[0.05338816358,0.02149069753,0.02865533642,0.04896485802,0.01624895062,0.03223945062,0.02425553395,0.03999812963,0.0361509784, 0],
					[0.02149069753,0.01468446914,0.01878391358,0.02441658642,0.008041938272,0.01002193827,0.01448993827,0.02536259259,0.02083593827, 0],
					[0.02865533642,0.01878391358,0.08550016358,0.06260714198,0.04439938272,0.01328671605,0.01043991049,0.06864603704,0.0420215216, 0],
					[0.04896485802,0.02441658642,0.06260714198,0.09546446914,0.05153806173,0.02902461728,0.02077028395,0.09002012963,0.03664589506, 0],
					[0.01624895062,0.008041938272,0.04439938272,0.05153806173,0.1278900988,0.0128384321,0.02091715432,0.1015344074,0.04497232099, 0],
					[0.03223945062,0.01002193827,0.01328671605,0.02902461728,0.0128384321,0.04125832099,0.01127854321,0.02960762963,0.02165332099, 0],
					[0.02425553395,0.01448993827,0.01043991049,0.02077028395,0.02091715432,0.01127854321,0.02883379321,0.02913762963,0.01739445988, 0],
					[0.03999812963,0.02536259259,0.06864603704,0.09002012963,0.1015344074,0.02960762963,0.02913762963,0.1467278889,0.05284057407, 0],
					[0.0361509784,0.02083593827,0.0420215216,0.03664589506,0.04497232099,0.02165332099,0.01739445988,0.05284057407,0.07926979321, 0],
					[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
		var returns = [0.06594444444,0.06155555556,0.1460555556,0.1734444444,0.1981111111,0.05511111111,0.1276111111,0.1903333333,0.1156111111, 0];
		var rf = 0;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf);
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #5');
	}
	
	// Test using static data
	// Test the case of only one portfolio on the efficient frontier
	{
		var expectedWeights = [0.4, 0.6];
		
		var covMat = [[1,0],[0,1]];
		var returns = [0.1, 0.2];
		var rf = 0;
		
		var maxSharpeRatioPortfolioWeights = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf, { constraints: {minWeights: [0.4, 0.6]} });
		assert.deepEqual(maxSharpeRatioPortfolioWeights, expectedWeights, 'Mean variance portfolio - Maximum Sharpe Ratio portfolio #6');
	}
	
	
	// Test using random data
	// Test the case of a too high risk free rate
	// Reference: Maximizing the Sharpe Ratio and Information Ratio in the Barra Optimizer, Leonid Kopman, Scott Liu
	{
		var covMat = [[0.05, 0.02],
		              [0.02, 0.07]];
		var returns = [0.05, 0.1];
		var rf = 0.1 + (1 - Math.random()); // rf will be > 0.1, which is the maximum return attainable
		
		assert.throws(function() { 
			PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf) },
			new Error('no corner portfolio with a strictly positive excess return'),
			"Mean variance portfolio - Maximum Sharpe Ratio portfolio #7");
	}
});


QUnit.test('Rounded weights portfolio', function(assert) {    
	// Example with static data
	{
		var testValues = [[0.7373, 0.2627, 0], [0.5759, 0.0671, 0.3570], [0.22, 0.66, 0.11, 0.01], [0.22, 0.66, 0.11, 0.01], [0.5, 0.5]];
		var testGridIndices = [10, 10, 1, 5, 1];
		var expectedValues = [[0.70, 0.30, 0], [0.60, 0.10, 0.30], [0, 1, 0, 0], [0.2, 0.6, 0.2, 0], [1, 0]];
		
		for (var i = 0; i < testValues.length; ++i) {
			assert.deepEqual(PortfolioAllocation.roundedWeights(testValues[i], testGridIndices[i]), expectedValues[i], 'Rounded weights #' + i);
		} 
	}  

});