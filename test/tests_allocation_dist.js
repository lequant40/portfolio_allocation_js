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
		var expectedWeights = [0.125, 0.047, 0.083, 0.613, 0.132];
		var expectedVolatility = 3.0406150258182514;
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[0][i] - expectedWeights[i]) <= 1e-3, true, 'ERC - Values #4 ' + i);
		}
		assert.equal(Math.abs(weights[1] - expectedVolatility) <= 1e-3, true, 'ERC - Values #4, volatility');
	}
	
	// Reference: https://github.com/lequant40/portfolio_allocation_js/issues/3
	//
	// This is a test of the ERC method wit a semi-positive definite covariance matrix
	{
		var cov = PortfolioAllocation.sampleCovarianceMatrix([-0.001697998787, 0.001427530069, -0.0005054947277, 0.0001213800916, -0.0001618204804],
                                                   [-0.002961208738, 0.001640864801, -0.0001441519189, 0.0004640260466, 0.000206138145],
                                                   [-0.001196934304, 0.0002291811055, -0.0004775812854, 0.001433428435, 0.0009680408618],
                                                   [0.0007749231464, -0.0002552713534, 0.0003744935825, 0.00247583719, 0.002435774483],
                                                   [-0.01361067487, -0.009173715734, -0.01167442776, -0.002090384233, -0.01495151011],
                                                   [-0.02573117139, 0.002675293971, 0.003048918071, -0.01685867172, -0.008688999321],
                                                   [-0.007799506409, 0.001570415597, -0.008377497996, -0.002298475244, -0.02687699473],
                                                   [-0.009589945124, -0.02074629016, 0.001121818449, -0.003655601888, 0.01279545706],
                                                   [-0.0001174153557, -0.0252073843, 0.008193945345, -0.006319268471, 0.002470726766],
                                                   [0.01103327496, 0.003395115191, -0.003901529538, 0.002079722704, -0.005188516084]);
		
		var weights = PortfolioAllocation.equalRiskContributionWeights(cov, {outputPortfolioVolatility: true}); 
		
		var expectedWeights = [0.27653353547900505, 0.08992490010322292, 0.07660067915254527, 0.10281081303356884, 0.00627339817004538, 0.07548939567866735, 0.004543021673884871, 0.0858506918300374, 0.02204534860588049, 0.25992821627314233];
		var expectedVolatility = 0;
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[0][i] - expectedWeights[i]) <= 1e-8, true, 'ERC - Values #5 ' + i);
		}
		assert.equal(Math.abs(weights[1] - expectedVolatility) <= 1e-8, true, 'ERC - Values #5, volatility');
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
				expectedWeights[i] /= sumWeights; // correlation rescaled + normalized weights = final expected weights
			}
			
			// Compare MDP weights to final expected weights
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
		var weights = PortfolioAllocation.minimumCorrelationWeights(sigma);
		
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
		var weights = PortfolioAllocation.minimumCorrelationWeights(covMat);
		
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
		var weights = PortfolioAllocation.proportionalMinimumVarianceWeights(cov);
		
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
		var expectedWeights =  [0.22372361188383866, 0, 0.308875157539891, 0.1338200074030125, 0.333581223173258];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #6 ' + i);
		}
	}
	
	// Reference: Private communication
	{
		var covMat = [[0.00008167944062,	0.00002381649417,	0.0000270928852,	0.00004867187199,	0.00003320921106],
					[0.00002381649417,	0.00002359207427,	0.000005565882395,	0.0000188629976,	0.00001333414531],
					[0.0000270928852,	0.000005565882395,	0.00001048947312,	0.0000168927258,	0.00001058070153],
					[0.00004867187199,	0.0000188629976,	0.0000168927258,	0.00003645410697,	0.00002279528797],
					[0.00003320921106,	0.00001333414531,	0.00001058070153,	0.00002279528797,	0.00001607412275]];
					
		// True expected values
		var weights = PortfolioAllocation.globalMinimumVarianceWeights(covMat); 
		var expectedWeights =  [0, 0.2145375758374286, 0.7854624241625715, 0, 0]
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - Values #6 ' + i);
		}
		
		// Wrong values, obtained by setting eps with a too big value
		var weights = PortfolioAllocation.globalMinimumVarianceWeights(covMat, {eps: 1e-4}); 
		var expectedWeights =  [0.2, 0.2, 0.2, 0.2, 0.2]
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'GMV - (Wrong) Values #7 ' + i);
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



QUnit.test('Numerical optimization portfolio', function(assert) {    
	// Test unsupported optimisation method
	{
		assert.throws(function() { 
			PortfolioAllocation.numericalOptimizationWeights(3,  function (arr) { return 1; }, {optimizationMethod: 'none'}) },
			new Error('unsupported optimisation method'),
			"Grid search portfolio - Unsupported optimisation method");
	}
	
	// Test grid search using static data
	{
		// Objective function: portfolio variance, for three assets
		// Covariance matrix taken from GMV test case #2
		// From GMV test case #2, expected (exact) weights are [0.9565927119697761, 0.043407288030223846, 0];
		function portfolio_variance_three_assets(arr) { 
			var covMat = [[0.0400, 0.0396, 0.0414], [0.0396, 0.0484, 0.0455], [0.0414, 0.0455, 0.0529]];
			
			return arr[0]*arr[0]*covMat[0][0] + arr[1]*arr[1]*covMat[1][1] + arr[2]*arr[2]*covMat[2][2] +
			2*arr[0]*arr[1]*covMat[0][1] + 2*arr[0]*arr[2]*covMat[0][2] + 2*arr[1]*arr[2]*covMat[1][2];
		}
		
		// Test default values for the grid search
		assert.deepEqual(PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets), 
		PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets, {optimizationMethod: 'grid-search', optimizationMethodParams: {k: 3}}), 
		'Numerical optimization portfolio - Default values');

		// Test finer and finer rational grids
		var expectedWeights = [[1,0,0]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets, {optimizationMethodParams: {k: 3}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio - Values #1' + i);
		}
		
		var expectedWeights = [[1,0,0]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets, {optimizationMethodParams: {k: 10}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio - Values #2' + i);
		}
		
		var expectedWeights = [[0.96,0.04,0]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets, {optimizationMethodParams: {k: 100}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio - Values #3' + i);
		}
		
		var expectedWeights = [[0.957,0.043,0]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_variance_three_assets, {optimizationMethodParams: {k: 1000}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio - Values #4' + i);
		}
	}

	// Test grid search with bounds contraints using static data
	{
		// Objective function: portfolio return, for three assets
		function portfolio_return_three_assets(arr) { 
			var ret = [-0.1, 0.3, .1];
			
			return -(arr[0]*ret[0] + arr[1]*ret[1] + arr[2]*ret[2]);
		}
		
		var expectedWeights = [[0.5,0.5,0]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_return_three_assets, {optimizationMethodParams: {k: 10}, constraints: {minWeights: [0.5, 0, 0]}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio, bounds contraints - Values #5' + i);
		}
		
		var expectedWeights = [[0,0.5,0.5]];
		var weights = PortfolioAllocation.numericalOptimizationWeights(3, portfolio_return_three_assets, {optimizationMethodParams: {k: 10}, constraints: {maxWeights: [1, 0.5, 1]}});
		for (var i = 0; i < weights.length; ++i) {
			var weightsOk = true;
			for (var j = 0; j < 3; ++j) {
			   if (expectedWeights[i][j] != weights[i][j]) {
				 weightsOk = false;
				 break;
			   }
			}	  
			assert.equal(weightsOk, true, 'Numerical optimization portfolio, bounds contraints - Values #6' + i);
		}
	}
	
	// TODO: Use random data with return as function + random contraints v.s. solving the same problem via linear programming
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




QUnit.test('Mean variance portfolio - error cases', function(assert) {    
	// Test using static data
	// Test the unsupported cases:
	// - Missing mean-variance optimization method
	// - Unsupported mean-variance optimization method
	{
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		assert.throws(function() { 
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat) },
			new Error('missing mean-variance optimization method'),
			"Mean variance portfolio - Missing mean-variance optimization method");

		assert.throws(function() { 
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'test', constraints: {return: 0.10}}) },
			new Error('unsupported mean-variance optimization method'),
			"Mean variance portfolio - Unsupported mean-variance optimization method");
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
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: unreachableMaxReturn}}) },
			new Error('return not reachable'),
			"Mean variance portfolio - Target return weights portfolio, unreachable target return #1");
			
		var minReturn = returns[0];
		var unreachableMinReturn = minReturn - (1 - Math.random()); // < minReturn
		assert.throws(function() { 
		PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: unreachableMinReturn}}) },
			new Error('return not reachable'),
			"Mean variance portfolio - Target return weights portfolio, unreachable target return #2");
	}
	
	// Test using static data
	{
		var covMat = [[1, 0.3], [0.3, 1]];
		var returns = [0.1, 0.2];
		
		var targetReturn = returns[1]; // the second asset
		var expectedWeights = [0, 1]; 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: targetReturn}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - target return weights portfolio #0');
		
		var targetReturn = 0.15; // exact mix of the two assets
		var expectedWeights = [0.5, 0.5];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: targetReturn}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - target return weights portfolio #0 bis');
	}
		
	// Test using static data (positive semi-definite covariance matrix)
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
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: 0.07}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #1');
		
		var expectedRoundedWeights = [0, 0, 0, 0.335, 0.455, 0, 0.21, 0, 0, 0]; // In Markowitz book, this is actually [0, 0, 0, 0.33, 0.455, 0, 0.215, 0, 0, 0], but this seems a wrong rounding 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: 0.175}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #2');
		
		var expectedRoundedWeights = [0, 0, 0, 0.3, 0.39, 0, 0.31, 0, 0, 0]; // In Markowitz book, this is actually [0, 0, 0, 0.29, 0.39, 0, 0.32, 0, 0, 0], but this seems a wrong rounding 
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: 0.169}});
		weights = PortfolioAllocation.roundedWeights(weights, 100);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target return weights portfolio #3');
		
		var expectedRoundedWeights = [0, 0, 0.09, 0.13, 0.21, 0, 0.57, 0, 0, 0];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetReturn', constraints: {return: 0.150}});
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
		var minVolatility = 0.12082760588883482; // computed thanks to the efficient frontier
		var maxVolatility = Math.sqrt(covMat[1][1]);
		for (var k = 0; k < nbSamples; ++k) {
			var unreachableMaxVolatility = generateRandomValue(maxVolatility - 1e-6, 2 * maxVolatility); // > maxVolatility
			assert.throws(function() { 
			PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetVolatility', constraints: {volatility: unreachableMaxVolatility}}) },
				new Error('volatility not reachable'),
				"Mean variance portfolio - Target volatility weights portfolio, unreachable target volatility #1 - " + k + "/" + nbSamples);
		
			var unreachableMinVolatility = generateRandomValue(0, minVolatility - 1e-6); // < minVolatility
			assert.throws(function() { 
			PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetVolatility', constraints: {volatility: unreachableMinVolatility}}) },
				new Error('volatility not reachable'),
				"Mean variance portfolio - Target volatility weights portfolio, unreachable target volatility #2 - " + k + "/" + nbSamples);
		}
	}
	
	// Test using static data (positive semi-definite covariance matrix)
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
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'targetVolatility', constraints: {volatility: 0.08}});
		weights = PortfolioAllocation.roundedWeights(weights, 200);
		assert.deepEqual(weights, expectedRoundedWeights, 'Mean variance portfolio - target volatility weights portfolio #1');
	}
	
});	


QUnit.test('Mean variance portfolio - maximum target volatility weights portfolio', function(assert) {    
	function generateRandomValue(minVal, maxVal) {	
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test using static data
	{
		var covMat = [[0.0146, 0.0187, 0.0145],
					 [0.0187, 0.0854, 0.0104],
					  [0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		var nbSamples = 100;
		var minAttainableVolatility = 0.12082760588883482; // computed thanks to the efficient frontier
		var maxAttainableVolatility = Math.sqrt(covMat[1][1]);
		for (var k = 0; k < nbSamples; ++k) {	
			// In case the maximum target volatility is lower than the mimimum attainable volatility,
			// the portfolio cannot be fully invested and must correspond to the efficient portfolio with cash
			// and with a target volatility constaint 
			var maxTargetVolatility = generateRandomValue(0, minAttainableVolatility - 1e-6); // < minAttainableVolatility
			var weights_mtv = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat,
			                                                                      { optimizationMethod: 'maximumTargetVolatility', constraints: {maxVolatility: maxTargetVolatility}});
			
			var investment = weights_mtv[0] + weights_mtv[1] + weights_mtv[2];
			assert.equal(investment < 1, true, 'Mean variance portfolio - maximum target volatility weights portfolio #1/1');
			
			var weights_tv = PortfolioAllocation.meanVarianceOptimizationWeights([returns[0], returns[1], returns[2], 0], 
			                                                                     [[covMat[0][0], covMat[0][1], covMat[0][2], 0],
                                                             					 [covMat[1][0], covMat[1][1], covMat[1][2], 0],
					                                                             [covMat[2][0], covMat[2][1], covMat[2][2], 0],
					                                                             [0, 0, 0, 0]], 
			                                                                     { optimizationMethod: 'targetVolatility', constraints: {volatility: maxTargetVolatility}});
			for (var i = 0; i < weights_mtv.length; ++i) {
				assert.equal(Math.abs(weights_mtv[i] - weights_tv[i]) <= 1e-8, true, 'Mean variance portfolio - maximum target volatility weights portfolio #1/2 ' + i);
			}
			assert.equal(Math.abs(1-investment - weights_tv[i]) <= 1e-8, true, 'Mean variance portfolio - maximum target volatility weights portfolio #1/2 ' + i);	
			
			// In case the maximum target volatility is greater than the maximum attainable volatility,
			// the portfolio must be the efficient portfolio with the highest return
			var maxTargetVolatility = generateRandomValue(maxAttainableVolatility - 1e-6, 2 * maxAttainableVolatility); // > maxAttainableVolatility
			var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, 
			                                                                  { optimizationMethod: 'maximumTargetVolatility', constraints: {maxVolatility: maxTargetVolatility}});
			assert.deepEqual(weights, [0, 1, 0], 'Mean variance portfolio - maximum target volatility weights portfolio #2');
			
			// Otherwise, the portfolio must correspond to the efficient portfolio with a target volatility constaint
			var maxTargetVolatility = generateRandomValue(minAttainableVolatility + 1e-6, maxAttainableVolatility - 1e-6); // > minAttainableVolatility && < maxAttainableVolatility
			var weights_mtv = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, 
			                                                                      { optimizationMethod: 'maximumTargetVolatility', constraints: {maxVolatility: maxTargetVolatility}});
			var weights_tv = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, 
			                                                                     { optimizationMethod: 'targetVolatility', constraints: {volatility: maxTargetVolatility}});
			assert.deepEqual(weights_mtv, weights_tv, 'Mean variance portfolio - maximum target volatility weights portfolio #3');
		}
	}
	
	// Test using static data
	// Test that a target maximum volatility of zero produces a zero-weights portfolio
	{
		// Problem data
		var covMat =[[0.0146, 0.0187, 0.0145],
					[0.0187, 0.0854, 0.0104],
					[0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Compute the desired maximum volatility
		var maxVolatility = 0;
		
		// Test that a 0 maximum volatility target is properly supported
		var expectedWeights = [0, 0, 0];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'maximumTargetVolatility', constraints: {maxVolatility: maxVolatility}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - maximum target volatility weights portfolio #4');
	}
	
	
	// Test using static data
	// Test the algorithm when all returns are negative; the associated portfolio must be null.
	{
		// Problem data
		var covMat = [[0.002533106361, 0.002161135405, 0.002605003686],
                     [0.002161135405, 0.002252008852, 0.002439575054],
                     [0.002605003686, 0.002439575054, 0.003315580611]];
		var returns = [-0.0000244756209, -0.003612424187, -0.004045817996];
		
		// Compute the desired maximum volatility
		var maxVolatility = 0.10/Math.sqrt(12);
		
		// Test that the algorithm is behaving properly
		var expectedWeights = [0, 0, 0];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'maximumTargetVolatility', constraints: {maxVolatility: maxVolatility}});
		assert.deepEqual(weights, expectedWeights, 'Mean variance portfolio - Test negative returns');
	}
});	


QUnit.test('Mean variance portfolio - minimum variance weights portfolio', function(assert) {    
	// Test that the GMV portfolio as computed by MVO algorithm is efficient, 
	// using a semi-definite positive covariance matrix
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02], 
		                                                                  [[0.0400, 0.0400], [0.0400, 0.0400]], 
																		  { optimizationMethod: 'minimumVariance'});
		var expectedWeights =  [0, 1];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance, efficient portfolio #1 ' + i);
		}

		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.02, 0.01], 
		                                                                  [[0.0400, 0.0400], [0.0400, 0.0400]], 
																		  { optimizationMethod: 'minimumVariance'});
		var expectedWeights =  [1, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance, efficient portfolio #2 ' + i);
		}
	}

	// Re-use the tests from the GMV portfolio
	// As the solution of the GMV portfolio is unique when the covariance matrix is positive definite, 
	// fake returns data are used below.
	
	// Reference: Portfolio Optimization versus Risk-Budgeting Allocation, Thierry Roncalli, WG RISK ESSEC, January 18, 2012
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02], 
		                                                                  [[0.0396, 0.0398], [0.0398, 0.0400]], 
																		  { optimizationMethod: 'minimumVariance'});
		var expectedWeights =  [1, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #1 ' + i);
		}
	}
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03], 
		                                                                  [[0.0400, 0.0396, 0.0414], [0.0396, 0.0484, 0.0455], [0.0414, 0.0455, 0.0529]], 
																		  { optimizationMethod: 'minimumVariance'});
		var expectedWeights =  [0.9565217391304353, 0.04347826086956469, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #2 ' + i);
		}
	}
	// Note: The volatility obtained by the computed portfolio is 0.019751470588235294, which is lower than the volatility of 0.2017 in the reference,
	// associated with the commented weights !
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03], 
		                                                                 [[0.0400, 0.0374, 0.0391], [0.0374, 0.0484, 0.0430], [0.0391, 0.0430, 0.0529]],
																		 { optimizationMethod: 'minimumVariance'}); 
		//var expectedWeights =  [0.7009, 0.2378, 0.0613];//
		var expectedWeights =  [0.8088235294117648, 0.19117647058823511, 0 ];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #3 ' + i);
		}
	}
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03], 
		                                                                  [[0.0400, 0.0418, 0.0437], [0.0418, 0.0484, 0.0481], [0.0437, 0.0481, 0.0529]],
																		  { optimizationMethod: 'minimumVariance'}); 
		var expectedWeights =  [1, 0, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #4 ' + i);
		}
	}
	
	// Reference: Understanding the Impact of Weights Constraints in Portfolio Theory, Thierry Roncalli
	{
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03, 0.04], 
		                                                                  [[0.0225,0.0030,0.0150,0.0225], [0.0030,0.0400,0.0350,0.0240], [0.0150,0.0350,0.0625,0.0600], [0.0225,0.0240,0.0600,0.0900]],
																	      { optimizationMethod: 'minimumVariance'}); 
		var expectedWeights =  [0.6548672566371683, 0.34513274336283173, 0, 0];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #5 ' + i);
		}
	}
	
	// Reference: Private communication with TrendXplorer
	{
		var covMat = [[0.03401428,0.0333167,-0.00614739,0.00926415,-0.0064081],
					  [0.0333167,0.06323421,-0.00855552,0.02245369,-0.00480642],
					  [-0.00614739,-0.00855552,0.01444902,-0.00432445,0.00690744],
					  [0.00926415,0.02245369,-0.00432445,0.02622712,0.0016983],
					  [-0.0064081,-0.00480642,0.00690744,0.00169834,0.0116492]];
		
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03, 0.04, 0.05], 
		                                                                  covMat,
																		  { optimizationMethod: 'minimumVariance'}); 
		//var expectedWeights = [0.22372361188383866, 0, 0.308875157539891, 0.1338200074030125, 0.333581223173258]; // these weights were obtained though the GMV method above, but are associated to a higher volatility: 0.004813123264651144 v.s. 0.004813123264643914 for the weights computed through MVO
		var expectedWeights = [0.2237239682833401, 0, 0.30887475456633273, 0.13381893271287823, 0.33358234443744894];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-8, true, 'Mean variance portfolio - minimum variance #6 ' + i);
		}
	}
	
	// Reference: Xi Bai, Katya Scheinberg, Reha Tutuncu, Least-squares approach to risk parity in portfolio selection
	{
		var covMat = [[94.868,33.750,12.325,-1.178,8.778],
					  [33.750,445.642,98.955,-7.901,84.954],
					  [12.325,98.955,117.265,0.503,45.184],
					  [-1.178,-7.901,0.503,5.460,1.057],
					  [8.778,84.954,45.184,1.057,34.126]];
			
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03, 0.04, 0.05], 
		                                                                  covMat,
																		  { optimizationMethod: 'minimumVariance'}); 
		var expectedWeights =  [0.050, 0.006, 0, 0.862, 0.082];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-3, true, 'Mean variance portfolio - minimum variance #7 ' + i);
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
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights([0.01, 0.02, 0.03, 0.04, 0.05], 
		                                                                  covMat, 
																		  { optimizationMethod: 'minimumVariance', constraints: {minWeights: minWeights, maxWeights: maxWeights} }); 
		var expectedWeights =  [0.200,0.050,0.050,0.350,0.350];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[i] - expectedWeights[i]) <= 1e-3, true, 'Mean variance portfolio - minimum variance #8, constrained ' + i);
		}
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


QUnit.test('Random subspace optimization method - error and limit cases', function(assert) {    
	// Test using static data
	// Test the unsupported cases:
	// - Unsupported subsets generation method
	// - Unsupported subsets aggregation method
	// - Wrong size of subsets of assets to generate
	// - No feasible portfolio generated
	// - Number of infeasible portfolios generated greater than the admissible ratio
	// - Wrong number of weights returned by the portfolio optimization method
	{
		// Problem data
		var nbAssets = 3;

		// Test an unsupported subsets generation method
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, null, { subsetsGenerationMethod: 'unknown' })},
			new Error('unsupported subsets of assets generation method'),
			"Random subspace optimization method - Unsupported subsets generation method");

		// Test an unsupported subsets aggregation method
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, null, { subsetsAggregationMethod: 'unknown' }) },
			new Error('unsupported aggregation method'),
			"Random subspace optimization method - Unsupported aggregation method");

		// Test wrong size of subsets of assets to generate: < 2
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, null, { sizeSubsets: 1 }) },
			new Error('the size of the subsets of assets must be between 2 and 3'),
			"Random subspace optimization method - Size of subsets of assets to generate < 2");

		// Test wrong size of subsets of assets to generate: > n
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, null, { sizeSubsets: 4 }) },
			new Error('the size of the subsets of assets must be between 2 and 3'),
			"Random subspace optimization method - Size of subsets of assets to generate > n");
		
		// Test no feasible portfolio generated
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, 
			                                                      function(idx, opt) { throw new Error('infeasible portfolio optimization problem'); }) },
			new Error('no feasible portfolio generated'),
			"Random subspace optimization method - No feasible portfolio generated");	
		
		// Test number of feasible portfolios generated less than the admissible ratio
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, 
			                                                      function(idx, opt) { if (idx[0] == 1) { throw new Error('infeasible portfolio optimization problem'); } else { return [0, 1]; } }, 
																  { subsetsGenerationMethod: 'deterministic',
                                                                    maxInfeasibleSubsetsRatio: 2/3}) },
			new Error('too many infeasible portfolios generated'),
			"Random subspace optimization method - Too many infeasible portfolios generated");
		
		// Test wrong number of weights returned by the portfolio optimization method
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, 
			                                                      function(idx, opt) { return [0, 0, 0, 1]; }) },
			new Error('internal error: the portfolio optimization method did not return the expected number of weights'),
			"Random subspace optimization method - Wrong number of weights returned by the portfolio optimization method");	
	}
	
	// Test using random data
	// Test the limit cases:
	// - 1 asset
	// - Number of assets equal to the size of the subsets to generate
	{
		// Test 1 asset
		{
			var expectedWeights = [1];
			var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(1, null);
			assert.deepEqual(weights, expectedWeights, 'Random subspace optimization method - 1 asset');	
		}
		
		// Test number of assets equal to the size of the subsets to generate
		//
		// For performances reasons, it must be tested that only one call to the portfolio optimization method takes place
		{
			// Generate a random number of assets
			var nbAssets = Math.floor(Math.random()*(50-2+1) + 2); // max 50 min 2
			
			// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
			function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
				// Count the number of function calls
				subsetEqualWeightsOptimization.called++;
				
				return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
			}
			subsetEqualWeightsOptimization.called = 0;
			
			// Compute the RSO equally weighted portfolio
			var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {sizeSubsets: nbAssets});
			
			// Compare the weights as well as the number of calls of the portfolio optimization method
			var expectedWeights = PortfolioAllocation.equalWeights(nbAssets);
			assert.deepEqual(weights, expectedWeights, 'Random subspace optimization method -  Number of assets equal to the size of the subsets to generate, #1');	
			
			var expectedFunctionCalls = 1;
			assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Number of assets equal to the size of the subsets to generate, #2');	
		}
	}
});	

QUnit.test('Random subspace optimization method', function(assert) {    
	// Tests the properties of the number of subsets to generate in case of random subsets generation, using random data 
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(25-2+1) + 2); // max 25 min 2
			
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// Count the number of function calls
			subsetEqualWeightsOptimization.called++;
			
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}
		
		// Compute the RSO equally weighted portfolio
		// The default number of subsets in case of random subsets generation method to generate must be equal to 128
		var expectedFunctionCalls = 128;
		subsetEqualWeightsOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization);
		assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Default number of random subsets to generate');
		
		// Compute the RSO equally weighted portfolio
		// The number of subsets to generate in case of random subsets generation method must be equal to the specified one
		var nbSubsetsToGenerate = Math.floor(Math.random()*(1000-1+1) + 1); // max 1000 min 1
		var expectedFunctionCalls = nbSubsetsToGenerate;
		subsetEqualWeightsOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {nbRandomSubsets: nbSubsetsToGenerate});
		assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Specified number of random subsets to generate');					
	}

	// Test the properties of the size of the generated subsets in case of random subsets generation, using random data  
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(50-2+1) + 2); // max 50 min 2
			
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// In case the size of the generated subsets is not correct, throw
			assert.equal(subsetAssetsIdx.length, expectedSizeSubsets, 'Random subspace optimization method - Size of random subsets to generate');
			
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}
		
		// Compute the RSO equally weighted portfolio
		// The default size of the subsets to generate must be FLOOR(SQRT(nbAssets))
		var expectedSizeSubsets = Math.max(2, Math.floor(Math.sqrt(nbAssets)));
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization);	
		
		// Compute the RSO equally weighted portfolio
		// The size of the generated subsets must be equal to the specified one
		var sizeSubsetsToGenerate = Math.floor(Math.random()*(nbAssets-2+1) + 2); // max nbAssets min 2
		var expectedSizeSubsets = sizeSubsetsToGenerate;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {sizeSubsets: sizeSubsetsToGenerate});			
	}
	
	// Tests the properties of the number of the subsets to generate in case of deterministic subsets generation, using random data 
	{
		function binomial(n, k) { if (k === 0) return 1; return (n * binomial(n-1, k-1)) / k; }
		
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(25-2+1) + 2); // max 25 min 2
			
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// Count the number of function calls
			subsetEqualWeightsOptimization.called++;
			
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}

		// Compute the RSO equally weighted portfolio
		// The number of subsets to generate in case of deterministic subsets generation method must be equal to Binomial(nbAssets, sizeSubsets)
		var sizeSubsetsToGenerate = Math.floor(Math.random()*(nbAssets-2+1) + 2); // max nbAssets min 2		
		var expectedFunctionCalls = binomial(nbAssets, sizeSubsetsToGenerate);
		subsetEqualWeightsOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {subsetsGenerationMethod: 'deterministic', sizeSubsets: sizeSubsetsToGenerate});
		assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Number of deterministic subsets to generate');
	}
	
	// Test the behaviour of the deterministic subsets generation method, using random data 
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(25-2+1) + 2); // max 25 min 2
			
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}

		// Compute the RSO equally weighted portfolio
		// In the specific case of the deterministic subsets generation method, an RSO equally weighted portfolio
		// must be equal to an equally weighted portfolio over all the original assets
		var expectedWeights = PortfolioAllocation.equalWeights(nbAssets);
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {subsetsGenerationMethod: 'deterministic'});
		
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-8) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Random subspace optimization method - Deterministic subsets generation method');
	}

	// Test the behaviour of the automated subsetting in case minWeights/maxWeights constraints are provided, using random data 
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(50-2+1) + 2); // max 50 min 2
				
		// Generate random lower and upper bounds on these assets 
		// (feasibility of the associated portfolio optimization problem is not relevant).
		var lowerBounds = new Array(nbAssets);
		var upperBounds = new Array(nbAssets);
		for (var i = 0; i < nbAssets; ++i) {
			lowerBounds[i] = Math.random();
			upperBounds[i] = Math.random();
		}
		
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// In case the automated subsetting of minWeights/maxWeights is not correct, throw
			var boundsOK = true;

			for (var i = 0; i < subsetAssetsIdx.length; ++i) {
				if (subsetPortfolioOptimizationMethodParams.constraints.minWeights[i] != lowerBounds[subsetAssetsIdx[i]-1]) {
					boundsOK = false;
					break;
				}
				if (subsetPortfolioOptimizationMethodParams.constraints.maxWeights[i] != upperBounds[subsetAssetsIdx[i]-1]) {
					boundsOK = false;
					break;
				}
			}
			assert.equal(boundsOK, true, 'Random subspace optimization method -  Automated subsetting of minWeights/maxWeights constraints');
			
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}
		
		// Compute the RSO equally weighted portfolio
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, { constraints: { minWeights: lowerBounds, maxWeights: upperBounds}});	
	}
	
	// Test the subsets aggregation method properties, using static data
	{
		var nbAssets = 3;
		
		// Define the portfolio optimization method to test the subsets aggregation method
		function subsetOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// Use the number of function calls to compute the subset portfolio weights
			subsetOptimization.called++;

			if (subsetOptimization.called === 1) {
				return [1, 0]; // subsetAssetsIdx is equal to [1, 2] -> [1, 0, 0] for original portfolio weights
			}
			else if (subsetOptimization.called === 2) {
				return [0, 1]; // subsetAssetsIdx is equal to [1, 3] -> [0, 0, 1] for original portfolio weights
			}
			else {
				return [0, 0]; // subsetAssetsIdx is equal to [2, 3] -> [0, 0, 0] for original portfolio weights
			}
		}
		
		// Compute the RSO portfolio
		// Test the 'average' subsets aggregation method
		var expectedWeights = [0.3333333333333333, 0, 0.3333333333333333];
		subsetOptimization.called = 0;	
		var averageWeights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetOptimization, {subsetsGenerationMethod: 'deterministic', subsetsAggregationMethod: 'average'});
		assert.deepEqual(averageWeights, expectedWeights, 'Random subspace optimization method - Average subsets aggregation method');

		// Compute the RSO portfolio
		// The default subsets aggregation method must be 'average'		
		subsetOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetOptimization, {subsetsGenerationMethod: 'deterministic'});
		assert.deepEqual(weights, averageWeights, 'Random subspace optimization method - Default subsets aggregation method');

		// Compute the RSO portfolio
		// Test the 'median' subsets aggregation method
		var expectedWeights = [0.21134201925334137, 0, 0.21134201925334137];
		subsetOptimization.called = 0;	
		var medianWeights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetOptimization, {subsetsGenerationMethod: 'deterministic', subsetsAggregationMethod: 'median'});
		assert.deepEqual(medianWeights, expectedWeights, 'Random subspace optimization method - Median subsets aggregation method');		
	}
});	



QUnit.test('Random subspace mean variance portfolio - error cases', function(assert) {    	
	// Test using static data
	// Test the unsupported cases:
	// - Missing mean-variance optimization method
	// - Target volatility not provided
	// - Unsupported mean-variance optimization method
	{
		// Problem data
		var covMat =[[0.0146, 0.0187, 0.0145],
					[0.0187, 0.0854, 0.0104],
					[0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Compute the target maximum volatility
		var maxVolatility = covMat[0][0];
		
		// Test missing mean-variance optimization method
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat)},
			new Error('missing mean-variance optimization method'),
			"Random subspace mean variance portfolio - Missing mean-variance optimization method");

		// Test no desired maximum volatility provided
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, { subsetPortfolioOptimizationMethodParams: { optimizationMethod: 'maximumTargetVolatility' } }) },
			new Error('missing mean-variance optimization method maximum volatility constraint'),
			"Random subspace mean variance portfolio - No desired maximum volatility");

		// Test an unsupported mean-variance optimization method
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, { subsetPortfolioOptimizationMethodParams: { optimizationMethod: 'unknown' } }) },
			new Error('unsupported mean-variance optimization method'),
			"Random subspace mean variance portfolio - Unsupported mean-variance optimization method");
	}
});	


QUnit.test('Random subspace mean variance portfolio', function(assert) {    
	// Test using static data
	{
		var covMat = [[1, 0.3, -0.2], 
		              [0.3, 1, 0.6], 
					  [-0.2, 0.6, 1]];
		var returns = [0.1, 0.2, 0.15];
		
		var expectedWeights = [0.09395200108235617, 0.1853986157617026, 0.11749197001156553];
		var weights = PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, 
		                                                                                { subsetsGenerationMethod: 'deterministic', 
																			              subsetPortfolioOptimizationMethodParams: {
																							  optimizationMethod: 'maximumTargetVolatility',
																							  constraints: {maxVolatility: Math.sqrt(0.10)}
																						  }
																						});
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-8) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Random subspace mean variance portfolio, #1');
	}
	
	// Test using static data volatility target, which can lead to infeasible portfolios
	// TODO

	// Test using static data return target, which can lead to infeasible portfolios
	// TODO

	// Test using static data
	// Reference: http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html, 100 Portfolios Formed on Size and Book-to-Market (10 x 10)
	// Returns data from 1963/07 to 1973/06 used to computed the expected returns/covariance matrix
	// 4 portfolios have been removed due to missing data
	//
	// This example is used to stress-test the algorithm, but it is not possible to compare the results to some expected weights 
	{
		var covMat = [[0.008626335034, 0.006715399917, 0.007042419992, 0.007048459, 0.0066366781, 0.006594382192, 0.006188482528, 0.006144514577, 0.005581524855, 0.005791371247, 0.006839503165, 0.006026119572, 0.00556042548, 0.004834846375, 0.004511975974, 0.00502810301, 0.004273827084, 0.004391065668, 0.004333287312, 0.005263026397, 0.006093797699, 0.005328759167, 0.004885570623, 0.005006898428, 0.004375058515, 0.004536227401, 0.00382383205, 0.004188885336, 0.004227284984, 0.004480470259, 0.005652719041, 0.004684925202, 0.004824278, 0.004694296476, 0.0038311798, 0.004131720831, 0.003918816918, 0.00363012853, 0.004303746538, 0.004208510197, 0.00500884843, 0.004497613403, 0.003925394458, 0.003826366745, 0.003479102582, 0.003217506668, 0.003092539967, 0.003387487433, 0.003395715859, 0.004000762357, 0.004921900893, 0.004177801985, 0.003654297707, 0.003598637829, 0.00347358152, 0.003243638023, 0.003030619633, 0.003440298944, 0.003931812779, 0.004004083997, 0.00379606579, 0.00338809299, 0.003208075825, 0.002935523324, 0.002579952426, 0.002747787374, 0.002993945136, 0.003157948431, 0.003715079122, 0.003337110603, 0.002848833667, 0.002949900901, 0.002982481082, 0.002356035502, 0.002656567017, 0.002544831954, 0.002931137112, 0.003179676958, 0.00410448382, 0.003094311644, 0.002412252584, 0.002356759221, 0.001942123379, 0.002101900141, 0.002560515946, 0.002299870077, 0.002083781247, 0.002701932892, 0.002757407576, 0.002009863366, 0.001980655134, 0.001876568152, 0.001624975344, 0.001797528616, 0.001651044294, 0.001615803773],
						[0.006715399917, 0.006683590368, 0.006014896374, 0.005908565288, 0.005722001587, 0.005670592412, 0.005404122902, 0.005228753424, 0.004781393518, 0.004867670897, 0.006057294555, 0.005460206075, 0.004834687069, 0.00427960842, 0.004094676601, 0.004430417541, 0.003713099085, 0.003918644507, 0.003805239886, 0.004724829646, 0.005547769669, 0.004641512869, 0.004275172741, 0.004535712695, 0.004002068481, 0.003998855101, 0.003461115736, 0.003868338392, 0.003743206232, 0.00397440804, 0.005132938989, 0.004315534403, 0.004334464152, 0.004243317135, 0.003586062754, 0.003673558324, 0.003599684906, 0.003328464186, 0.003895523432, 0.003933823625, 0.00454823586, 0.004173927655, 0.00350670524, 0.003414971713, 0.003313012908, 0.003017905708, 0.002835420141, 0.00304936103, 0.003091043629, 0.00371479407, 0.004518149563, 0.003944750725, 0.003314415963, 0.00343389213, 0.003213881416, 0.003034389219, 0.002929474526, 0.00334420473, 0.003452367593, 0.003555839981, 0.003597330441, 0.003228763006, 0.002824046695, 0.002719042358, 0.002595143005, 0.00264760704, 0.002929607457, 0.003064194774, 0.003378001969, 0.003206598288, 0.002600477361, 0.002858297419, 0.002724947082, 0.002440528653, 0.002583794193, 0.002526042736, 0.002966731622, 0.003089654503, 0.003539376119, 0.002988376321, 0.002327715382, 0.0023059398, 0.0020451527, 0.001968497738, 0.002535640801, 0.002346819807, 0.002105827515, 0.00248426969, 0.002719471677, 0.002091166113, 0.00204061465, 0.001945355621, 0.00177718764, 0.001833740965, 0.001833126743, 0.001802919713],
						[0.007042419992, 0.006014896374, 0.006439832216, 0.006144856481, 0.005759584222, 0.005845477915, 0.005529484399, 0.005390914697, 0.004908853124, 0.005103830919, 0.006073314853, 0.005499701853, 0.004969641429, 0.004398641598, 0.004134976123, 0.004525198382, 0.003883988522, 0.003985939094, 0.003937507735, 0.004796383717, 0.00554862731, 0.00472970036, 0.004383809539, 0.004507738633, 0.004099612817, 0.004120392418, 0.003465223346, 0.003916806442, 0.003893286734, 0.00427363146, 0.005206899053, 0.004323028253, 0.004343625848, 0.00432711529, 0.003575259214, 0.003722534095, 0.003569795438, 0.003341274835, 0.00395494306, 0.004053391677, 0.0045785797, 0.004162307351, 0.003616132532, 0.003500002623, 0.00324854751, 0.003014847819, 0.002964497291, 0.003052088947, 0.003176130571, 0.003789134736, 0.004571681602, 0.003864248166, 0.003371414581, 0.003366216157, 0.003231814936, 0.003016512653, 0.002922320365, 0.003274787097, 0.003685314298, 0.003723157506, 0.003694396668, 0.003228902136, 0.002950115286, 0.002751010898, 0.002432745554, 0.002705575637, 0.002880348639, 0.002975663554, 0.003502498329, 0.00322837047, 0.002724240891, 0.002751144754, 0.002720678945, 0.00232456853, 0.00253747154, 0.002549888275, 0.002757467982, 0.003098168353, 0.003999211129, 0.002958428789, 0.002268515803, 0.002255255898, 0.00193368615, 0.00201404631, 0.00253937852, 0.002258501074, 0.002197152045, 0.002615141548, 0.002802102426, 0.001993131787, 0.001941433241, 0.001881035982, 0.001629664591, 0.001757523692, 0.001675342356, 0.001705850964],
						[0.007048459, 0.005908565288, 0.006144856481, 0.006721824469, 0.005852814111, 0.005856465086, 0.005599641385, 0.005466069036, 0.005048488162, 0.005242939938, 0.006111949473, 0.005465037035, 0.005129873692, 0.004416730131, 0.004205399686, 0.004458656717, 0.003886593217, 0.003969078367, 0.003980490598, 0.004776577303, 0.005403948677, 0.004734314166, 0.004410779587, 0.004359942384, 0.004043738715, 0.00403985557, 0.003443041304, 0.003961335731, 0.003888049826, 0.004038458787, 0.004970339251, 0.004186782856, 0.004351604312, 0.004236631295, 0.003515810312, 0.003803261423, 0.003629734042, 0.003345952371, 0.003877361097, 0.003927757425, 0.004493401784, 0.004028057711, 0.003580285672, 0.00341839726, 0.003306099986, 0.003046845036, 0.002896100679, 0.00306826418, 0.003125160756, 0.003719892798, 0.004420189562, 0.003735997127, 0.003210599149, 0.003358889658, 0.003149259674, 0.002974762136, 0.002953982871, 0.003297600332, 0.003627668097, 0.003716680489, 0.003460502947, 0.003201833917, 0.002927183154, 0.002807933105, 0.002393663894, 0.002612167028, 0.002891811709, 0.002789078706, 0.00345071268, 0.003015226454, 0.002634119186, 0.002712130997, 0.002692079815, 0.00220023004, 0.002432468817, 0.002479192684, 0.002699152581, 0.002967319487, 0.003811419029, 0.002813450776, 0.002223566, 0.002173940325, 0.001910269587, 0.001976848994, 0.002485448666, 0.002206599569, 0.002148219717, 0.002382785761, 0.002707287202, 0.001804422913, 0.001881753917, 0.001755480548, 0.001556326854, 0.001713310797, 0.001644961803, 0.001681966033],
						[0.0066366781, 0.005722001587, 0.005759584222, 0.005852814111, 0.005804198968, 0.005518711262, 0.005212667932, 0.00520579585, 0.004695260663, 0.004884578153, 0.00563581538, 0.0051420724, 0.004752966669, 0.004155890204, 0.003897907344, 0.004250814324, 0.00363196973, 0.003733774514, 0.003756053426, 0.004551709944, 0.005102927184, 0.004374288211, 0.004156728923, 0.004290280092, 0.003757206817, 0.00387302364, 0.003394904912, 0.003736503159, 0.00370489262, 0.003989563848, 0.004710173158, 0.003965921161, 0.004097876261, 0.003968104092, 0.003341734566, 0.003518722423, 0.003444342914, 0.00311780955, 0.003726171349, 0.003807930121, 0.004239428375, 0.003945152337, 0.003383915216, 0.003340113351, 0.003251743842, 0.002902764019, 0.002782632546, 0.002923215582, 0.003007653059, 0.003567505347, 0.004212690159, 0.003554662776, 0.003178857495, 0.003158288587, 0.003031807835, 0.002869155025, 0.002756593414, 0.003070416077, 0.003440003388, 0.003495357866, 0.003285121985, 0.003009518249, 0.002794968285, 0.002665138781, 0.002399572175, 0.002482747036, 0.002766082028, 0.002807242675, 0.003305060105, 0.002907239952, 0.002521817319, 0.002658128393, 0.002581899568, 0.002227657241, 0.002470283393, 0.002403485387, 0.002673820064, 0.002835663113, 0.003627688117, 0.002699749688, 0.002162000371, 0.002126479925, 0.001833840701, 0.001935762069, 0.002350195792, 0.0021989561, 0.002055965995, 0.002359758516, 0.002504526251, 0.001829998069, 0.001840005844, 0.001776080995, 0.001555161991, 0.001703873611, 0.001733502399, 0.001681532412],
						[0.006594382192, 0.005670592412, 0.005845477915, 0.005856465086, 0.005518711262, 0.005885907481, 0.005249993788, 0.005186831761, 0.004725555913, 0.00494187006, 0.005657857003, 0.005163692584, 0.004706632845, 0.004173076723, 0.003902254463, 0.004284994782, 0.00371978044, 0.003800660873, 0.003827532969, 0.004633361599, 0.005034250859, 0.004381039134, 0.004140177808, 0.004271141624, 0.003815087583, 0.003888776041, 0.003412560729, 0.003715946693, 0.003694375629, 0.003998670233, 0.00466125509, 0.003956994273, 0.004064539484, 0.003977163851, 0.003358623301, 0.003518226553, 0.003373116109, 0.003176878082, 0.003775985543, 0.003852669217, 0.004223583686, 0.003867575054, 0.003372837448, 0.003276297438, 0.003134634966, 0.002941816614, 0.002837975836, 0.002966669732, 0.00306358095, 0.003607251209, 0.004047749529, 0.003518512219, 0.003220634491, 0.003225088129, 0.003074806453, 0.002869951923, 0.002794455413, 0.00314857206, 0.003445206259, 0.003484276388, 0.003242065918, 0.00296896478, 0.002867250372, 0.002673288533, 0.002378738807, 0.002476341968, 0.00278931335, 0.002744695175, 0.003419811664, 0.002865320479, 0.00253830332, 0.002672007559, 0.002546935189, 0.002227989655, 0.00248874816, 0.002472291663, 0.002701562387, 0.002937912347, 0.003746374369, 0.002725652958, 0.002111424879, 0.002183330207, 0.001826697256, 0.001950773694, 0.00231638808, 0.002196664573, 0.002116575774, 0.002469559865, 0.002616001643, 0.001804342367, 0.00183879967, 0.001747813931, 0.001571649058, 0.001717569968, 0.001702643639, 0.001731275443],
						[0.006188482528, 0.005404122902, 0.005529484399, 0.005599641385, 0.005212667932, 0.005249993788, 0.00529751027, 0.00488858012, 0.004479292431, 0.004661042738, 0.005244524452, 0.004844508527, 0.004519702335, 0.003901269768, 0.003814645021, 0.004000603532, 0.003560416181, 0.003645410684, 0.003622276933, 0.004387740333, 0.004834011926, 0.004163941856, 0.003798179843, 0.003995557898, 0.00373019354, 0.003728689545, 0.00322625399, 0.003643686193, 0.003498595594, 0.00389183741, 0.00449380729, 0.003763685986, 0.003833958135, 0.003889181773, 0.003237279705, 0.003384040493, 0.003371544744, 0.003103820812, 0.003638987516, 0.003674315677, 0.003997987643, 0.00374307051, 0.003281052417, 0.003160486446, 0.003002804477, 0.002773694236, 0.002668950943, 0.002792739793, 0.0028610154, 0.003470426784, 0.003914053963, 0.003418586327, 0.003024243577, 0.003035802487, 0.002865340173, 0.002746671583, 0.002713610852, 0.002952582409, 0.003384751458, 0.003398783072, 0.003164531892, 0.002893749342, 0.002576146487, 0.002534087762, 0.002298139155, 0.002369198989, 0.002650016671, 0.002596581078, 0.003333011815, 0.002746955289, 0.0024088781, 0.002490168387, 0.002476504628, 0.002163291487, 0.002345133043, 0.002378529183, 0.002524251707, 0.002792157251, 0.003547327024, 0.002582676068, 0.001961667199, 0.00198535973, 0.001802418959, 0.00179555948, 0.002302544639, 0.002102420911, 0.001963682724, 0.002249388048, 0.002506105964, 0.001728073759, 0.001709216314, 0.001684785524, 0.0015206781, 0.00160284881, 0.001599105457, 0.001639530125],
						[0.006144514577, 0.005228753424, 0.005390914697, 0.005466069036, 0.00520579585, 0.005186831761, 0.00488858012, 0.00512539599, 0.004447562299, 0.004656774587, 0.005183246677, 0.00477351631, 0.004412606969, 0.003868620735, 0.003722629721, 0.003997654346, 0.00351483726, 0.003593304073, 0.00361034607, 0.004428704724, 0.004702079826, 0.004092405047, 0.003907918293, 0.004086110471, 0.003608937357, 0.003750278469, 0.003359302132, 0.003640954706, 0.003612722983, 0.003832686285, 0.004280689398, 0.003780650207, 0.003853131458, 0.003750471713, 0.003221510183, 0.003375329421, 0.00337204551, 0.003078967411, 0.003631111805, 0.003705470881, 0.004051691901, 0.003758263597, 0.003219548501, 0.00322410537, 0.003053520724, 0.00282870831, 0.002747014377, 0.002873165543, 0.002847918425, 0.003448010781, 0.003931554582, 0.003313178608, 0.003028926938, 0.003032511577, 0.002911464346, 0.002745567374, 0.002658174324, 0.003029978077, 0.003456979638, 0.003418816348, 0.003089516401, 0.002823881199, 0.002694148822, 0.002581815013, 0.002274152315, 0.002452890052, 0.002639344627, 0.002777070666, 0.003271405217, 0.002731006516, 0.002417733311, 0.002545041644, 0.002445071567, 0.002144393999, 0.002389854799, 0.002344599062, 0.002619809282, 0.002852917885, 0.003539833913, 0.002638401879, 0.002021217408, 0.002049355105, 0.001773477959, 0.001877138149, 0.002263133226, 0.00219391231, 0.002079094009, 0.002384093743, 0.002466890693, 0.001729753224, 0.001775353392, 0.001711371032, 0.001494355257, 0.001666919889, 0.001640984047, 0.001597138892],
						[0.005581524855, 0.004781393518, 0.004908853124, 0.005048488162, 0.004695260663, 0.004725555913, 0.004479292431, 0.004447562299, 0.004275795316, 0.004278451337, 0.004720826626, 0.004336785326, 0.004055297277, 0.003554139547, 0.003381220835, 0.003618877409, 0.003181156489, 0.003217439503, 0.003259972434, 0.004035653591, 0.004306432777, 0.003736776838, 0.003600043099, 0.003660395382, 0.003208342404, 0.003419811102, 0.002923043046, 0.003287871386, 0.003183279271, 0.003487995237, 0.003933100761, 0.003328081687, 0.003459011173, 0.003390186576, 0.002954394167, 0.002954282232, 0.003005882414, 0.002767666317, 0.003242079448, 0.003357423284, 0.003561854758, 0.003334818137, 0.002942575524, 0.002935013377, 0.002794287977, 0.002568939824, 0.002487098668, 0.002538668543, 0.002648547527, 0.003152359262, 0.003521043502, 0.003050485455, 0.00263296343, 0.002732898344, 0.002639358885, 0.00245577162, 0.00240993133, 0.002745503964, 0.003105609008, 0.003080811196, 0.002765467253, 0.002570024624, 0.002403539014, 0.002365928298, 0.002060583994, 0.002128802882, 0.002378697448, 0.002375537224, 0.002979944128, 0.002466248967, 0.002183305605, 0.00233221606, 0.002247460908, 0.00188411888, 0.002077564129, 0.002123604214, 0.002368484015, 0.002484697379, 0.003263079647, 0.002315163065, 0.001846701912, 0.00181648464, 0.001622591534, 0.001671580843, 0.002067975693, 0.001944108977, 0.001812298457, 0.002112557084, 0.002276282331, 0.001537089641, 0.001600344549, 0.001539839272, 0.001380351461, 0.001445855309, 0.001467154418, 0.001461220643],
						[0.005791371247, 0.004867670897, 0.005103830919, 0.005242939938, 0.004884578153, 0.00494187006, 0.004661042738, 0.004656774587, 0.004278451337, 0.004640905014, 0.00482602218, 0.004458280353, 0.004183449769, 0.003645954911, 0.003425163222, 0.003676113678, 0.003263922036, 0.003322619244, 0.00337732856, 0.004206406044, 0.004412656627, 0.003819775592, 0.003648447356, 0.003715759146, 0.003384254134, 0.003439829175, 0.003038721767, 0.003383124035, 0.003306777089, 0.003614084622, 0.003957043704, 0.003394430058, 0.003520707926, 0.00340538183, 0.002993569578, 0.003161183487, 0.002984117306, 0.002825546864, 0.003317706912, 0.003465333124, 0.003620403799, 0.00336627656, 0.002999951399, 0.002888635953, 0.002803079025, 0.002559037881, 0.002522217349, 0.002613303501, 0.002654684242, 0.003261252232, 0.003583525581, 0.003009493202, 0.002728739572, 0.00276611704, 0.00257947995, 0.002479049424, 0.002497650293, 0.00278677727, 0.003201191769, 0.003171901734, 0.002852639794, 0.002556324433, 0.002395752497, 0.002348819196, 0.002024066215, 0.002145725437, 0.002443764267, 0.002445315495, 0.002970631164, 0.002471499072, 0.00217780669, 0.002264720672, 0.002206492379, 0.001912792792, 0.002093213831, 0.00213031204, 0.00237844092, 0.002584543946, 0.003406758634, 0.002316135476, 0.001841431471, 0.001807023938, 0.001595657068, 0.001697851274, 0.002067148721, 0.001983828019, 0.001884783965, 0.002081992134, 0.00232940124, 0.001503426019, 0.001564091681, 0.001499557562, 0.001293417718, 0.001433147338, 0.001459815065, 0.001456114868],
						[0.006839503165, 0.006057294555, 0.006073314853, 0.006111949473, 0.00563581538, 0.005657857003, 0.005244524452, 0.005183246677, 0.004720826626, 0.00482602218, 0.00715380894, 0.00568062383, 0.005189433557, 0.004630976216, 0.004338145108, 0.004523508387, 0.003894140167, 0.003935601645, 0.003681198312, 0.004681323607, 0.005995902923, 0.005062114042, 0.004514844807, 0.004573083939, 0.004188643902, 0.003980118185, 0.003462811048, 0.003809565319, 0.003914667109, 0.004022603493, 0.005713294222, 0.004837360272, 0.00468190493, 0.004517144876, 0.003730300128, 0.003878401571, 0.003657651758, 0.003235475906, 0.003778763194, 0.00390359752, 0.005197194723, 0.00435597295, 0.003811533041, 0.003647806326, 0.003519599805, 0.003250822993, 0.002851594993, 0.003133866747, 0.003116543132, 0.003691230013, 0.004951209144, 0.004333963728, 0.003606982422, 0.003729182959, 0.003338004342, 0.003044495193, 0.002987117227, 0.003204243619, 0.003555873242, 0.003764693424, 0.004334541133, 0.003553878345, 0.003231203484, 0.002943099406, 0.002696793804, 0.002640047186, 0.003154011213, 0.0030178547, 0.003233047781, 0.003705899965, 0.003057049188, 0.003098922292, 0.002903858749, 0.00246176927, 0.002599881832, 0.002592969402, 0.002729180908, 0.003147871729, 0.003568750149, 0.003247957188, 0.002593234347, 0.002508700781, 0.002074163928, 0.00205201753, 0.002574940276, 0.002281943433, 0.002136746694, 0.002578710274, 0.002542355418, 0.002360605143, 0.002220554732, 0.002113770986, 0.001821547472, 0.001864275149, 0.001779315548, 0.00185400446],
						[0.006026119572, 0.005460206075, 0.005499701853, 0.005465037035, 0.0051420724, 0.005163692584, 0.004844508527, 0.00477351631, 0.004336785326, 0.004458280353, 0.00568062383, 0.005716769838, 0.0045563704, 0.004052643741, 0.003908646935, 0.004102625285, 0.003502287588, 0.003641064362, 0.003561828604, 0.004279809736, 0.005264020749, 0.004513935274, 0.004092159872, 0.004119086638, 0.003783692589, 0.003874680981, 0.003291156127, 0.003635289158, 0.003635098529, 0.003767603955, 0.005032640564, 0.00427678359, 0.004236813741, 0.004152959872, 0.003422740202, 0.003607572409, 0.003390717311, 0.003054099004, 0.003663946352, 0.003692743513, 0.004504288936, 0.004016010617, 0.003455978013, 0.003492681159, 0.003262636438, 0.003004556681, 0.002798426204, 0.002855994378, 0.002949423992, 0.003455392502, 0.004442682175, 0.003868226071, 0.003236383866, 0.003332541527, 0.003093094158, 0.002851603943, 0.002816757742, 0.003150937277, 0.003371518839, 0.003421802213, 0.003810412801, 0.003295890373, 0.002979283605, 0.002720910432, 0.002541390204, 0.002527176121, 0.002777588485, 0.00282129975, 0.003289166972, 0.003287562718, 0.00280330556, 0.002900054774, 0.00263994758, 0.002390474875, 0.002500464799, 0.002486099168, 0.002786945071, 0.002980381819, 0.003553381874, 0.002975489791, 0.002443998787, 0.002360782909, 0.002067739944, 0.001998907565, 0.002534972683, 0.002236704862, 0.00217467103, 0.002359074018, 0.002538706295, 0.002178968325, 0.002083663409, 0.001987904628, 0.001794990531, 0.001780425059, 0.001803667429, 0.001750839599],
						[0.00556042548, 0.004834687069, 0.004969641429, 0.005129873692, 0.004752966669, 0.004706632845, 0.004519702335, 0.004412606969, 0.004055297277, 0.004183449769, 0.005189433557, 0.0045563704, 0.005036636106, 0.0038638914, 0.003697488713, 0.003872703062, 0.003331785478, 0.003494649964, 0.00343602649, 0.004112031878, 0.004721710351, 0.004106976616, 0.003784774051, 0.003711100379, 0.00356760794, 0.003569918691, 0.003247660017, 0.003570598902, 0.00352122161, 0.003795115997, 0.004639492099, 0.003769156389, 0.00392091161, 0.003720303187, 0.003305639622, 0.003329627441, 0.00330156276, 0.002915103266, 0.00347520935, 0.003564945967, 0.004115633012, 0.003837504568, 0.003139800334, 0.003134726059, 0.003276090681, 0.002948850289, 0.002638145651, 0.002842945493, 0.003000883638, 0.003455188678, 0.00405233259, 0.003689967188, 0.002976667966, 0.00317508292, 0.003002616169, 0.002716058382, 0.002865345579, 0.002911171213, 0.0033207877, 0.003383236234, 0.003420739267, 0.003050888097, 0.002846153878, 0.002658557789, 0.00248201665, 0.002494794761, 0.002796370631, 0.002825190452, 0.003221488961, 0.002953910851, 0.002647777308, 0.00267266428, 0.002675630159, 0.002214654882, 0.002358356459, 0.002414975934, 0.002506934412, 0.002874428286, 0.003548025327, 0.002676243201, 0.00227159661, 0.002203276806, 0.002040575238, 0.001982869769, 0.002484016339, 0.002223537119, 0.002091956136, 0.002527216656, 0.00245712161, 0.00192903563, 0.002014819832, 0.00197991251, 0.001669107432, 0.001871868301, 0.001839191813, 0.001779332353],
						[0.004834846375, 0.00427960842, 0.004398641598, 0.004416730131, 0.004155890204, 0.004173076723, 0.003901269768, 0.003868620735, 0.003554139547, 0.003645954911, 0.004630976216, 0.004052643741, 0.0038638914, 0.003873833194, 0.003250860823, 0.003363765299, 0.002882109091, 0.003108371262, 0.002917879079, 0.003545322717, 0.004265575331, 0.00365789709, 0.003297331834, 0.003420793192, 0.003095560522, 0.003088608654, 0.002741008061, 0.00297496778, 0.002989693088, 0.003099461678, 0.003987789659, 0.003450262327, 0.003397950778, 0.00336390694, 0.002851727225, 0.002891509442, 0.002818125281, 0.002423383822, 0.002985284814, 0.003095657796, 0.003702221795, 0.003339200871, 0.002930283919, 0.002769878876, 0.002819545433, 0.002524422025, 0.002300083628, 0.002393712985, 0.002549847288, 0.002984198221, 0.003592822144, 0.00315543466, 0.002677606164, 0.002848600657, 0.002637211277, 0.002388155994, 0.002352838588, 0.002637722498, 0.00285052315, 0.002949664632, 0.003039306391, 0.002737754808, 0.002598694665, 0.002320043343, 0.002149391811, 0.002141692279, 0.002455331697, 0.002372667954, 0.002748106064, 0.002744602577, 0.002357861268, 0.002469088971, 0.002266383371, 0.001915814607, 0.002137951166, 0.002032070219, 0.002162047775, 0.002401755863, 0.003083914722, 0.002498508202, 0.002029350921, 0.002015332046, 0.001718203143, 0.001682071098, 0.002112971385, 0.001897921416, 0.001773218193, 0.002041876751, 0.002086498421, 0.001765418608, 0.001761942727, 0.001704274891, 0.00152144715, 0.001538600425, 0.001536824264, 0.00156334151],
						[0.004511975974, 0.004094676601, 0.004134976123, 0.004205399686, 0.003897907344, 0.003902254463, 0.003814645021, 0.003722629721, 0.003381220835, 0.003425163222, 0.004338145108, 0.003908646935, 0.003697488713, 0.003250860823, 0.00382247693, 0.003213676416, 0.002929488868, 0.002993233232, 0.002912453456, 0.003561962143, 0.004050883743, 0.003683104325, 0.003225346346, 0.003335041919, 0.003178531074, 0.003193442983, 0.002712119593, 0.003110599618, 0.003033533277, 0.003223937921, 0.003787378884, 0.003347239651, 0.003280235516, 0.003345649319, 0.002756067678, 0.002767110209, 0.002914685529, 0.002613713399, 0.002927980931, 0.003180273915, 0.003624700079, 0.003335812529, 0.002856828165, 0.002905681077, 0.002759159271, 0.002564853416, 0.002352480996, 0.002384606717, 0.002531822685, 0.002957736965, 0.003428895014, 0.003150561594, 0.00263577638, 0.002772655784, 0.00257720391, 0.002387000641, 0.002357687016, 0.002638124, 0.002886483586, 0.002862711554, 0.002969202312, 0.002639004767, 0.002423556742, 0.002302659107, 0.002211500634, 0.002164797252, 0.002527220919, 0.002339467906, 0.002896492621, 0.002636127797, 0.002371647666, 0.002451257932, 0.002351451581, 0.001999883551, 0.002172617411, 0.002162700965, 0.002268925825, 0.002563441418, 0.003174658432, 0.002401444852, 0.001947953835, 0.001925131204, 0.001891678493, 0.001711701281, 0.002146911221, 0.002034977418, 0.001895268085, 0.002260311334, 0.002183371878, 0.001703799545, 0.001711400505, 0.001713412, 0.001534514796, 0.00171764512, 0.001629135281, 0.00166261154],
						[0.00502810301, 0.004430417541, 0.004525198382, 0.004458656717, 0.004250814324, 0.004284994782, 0.004000603532, 0.003997654346, 0.003618877409, 0.003676113678, 0.004523508387, 0.004102625285, 0.003872703062, 0.003363765299, 0.003213676416, 0.003883552949, 0.002967310677, 0.003138069179, 0.002950611747, 0.0036818911, 0.004142008143, 0.003493402694, 0.003368136036, 0.003444248522, 0.003129152572, 0.003263735276, 0.002896374512, 0.003062725198, 0.003016611448, 0.003231478037, 0.003999979567, 0.003306005668, 0.00343324518, 0.003199677237, 0.002782050118, 0.002870923467, 0.002859287353, 0.00266264615, 0.003115850747, 0.003073912489, 0.003631119945, 0.003343044549, 0.00289164188, 0.002883360707, 0.002722695521, 0.002568207708, 0.002351255139, 0.002444673168, 0.002568367833, 0.003064479536, 0.003569656968, 0.003122008957, 0.00274180904, 0.002790751101, 0.00260473893, 0.002401464726, 0.002339342052, 0.002678514626, 0.002819510233, 0.002970506193, 0.002821127856, 0.002624749646, 0.002501434852, 0.002245179492, 0.002192825626, 0.002164558022, 0.002394731624, 0.002495798699, 0.002804875205, 0.002463355151, 0.002245615832, 0.002447919046, 0.002278303381, 0.001985647385, 0.002228671992, 0.002114276384, 0.002285185976, 0.002486016887, 0.003060032788, 0.002355436976, 0.001903773391, 0.001912987547, 0.001729006647, 0.001644036295, 0.002084742402, 0.00191297072, 0.001819593425, 0.002199419717, 0.002243951131, 0.001705130056, 0.001712051313, 0.001694201176, 0.001474786573, 0.001594540806, 0.001590322296, 0.001530771977],
						[0.004273827084, 0.003713099085, 0.003883988522, 0.003886593217, 0.00363196973, 0.00371978044, 0.003560416181, 0.00351483726, 0.003181156489, 0.003263922036, 0.003894140167, 0.003502287588, 0.003331785478, 0.002882109091, 0.002929488868, 0.002967310677, 0.003044465794, 0.002733450929, 0.002635010927, 0.003274962376, 0.003681445239, 0.003179269931, 0.002838947865, 0.003071383135, 0.002865501922, 0.002941029134, 0.002605314467, 0.002783878304, 0.002705720095, 0.00292749477, 0.003485988704, 0.00309340092, 0.003076533575, 0.002884361566, 0.002539860536, 0.002578019495, 0.002610227175, 0.002350060347, 0.002642787593, 0.002800009414, 0.003194270936, 0.002960725557, 0.002526066962, 0.002600968519, 0.002373548708, 0.00233798433, 0.002110325803, 0.002179840156, 0.002177178217, 0.002647144259, 0.002943346366, 0.002783759962, 0.00247288644, 0.002423331136, 0.002272073889, 0.002157631265, 0.002150761998, 0.002429601622, 0.00258608262, 0.002629567419, 0.002489536208, 0.002321161174, 0.00216214058, 0.002093913888, 0.001967676888, 0.001976377986, 0.00212531014, 0.002052453022, 0.002483985727, 0.00221991072, 0.002045275133, 0.002103865601, 0.00200061205, 0.001752657322, 0.001914621865, 0.001904786973, 0.00196953169, 0.002308013374, 0.002725449108, 0.002072097471, 0.00164384809, 0.001742680996, 0.001600301608, 0.001561717078, 0.001859878712, 0.001712852285, 0.001679136737, 0.001911216259, 0.001958052651, 0.001415231144, 0.00147673357, 0.001456591754, 0.001216366287, 0.001374620282, 0.001367862053, 0.001422734186],
						[0.004391065668, 0.003918644507, 0.003985939094, 0.003969078367, 0.003733774514, 0.003800660873, 0.003645410684, 0.003593304073, 0.003217439503, 0.003322619244, 0.003935601645, 0.003641064362, 0.003494649964, 0.003108371262, 0.002993233232, 0.003138069179, 0.002733450929, 0.00321693965, 0.0027661938, 0.003352791662, 0.003773810179, 0.003292673967, 0.002995705642, 0.003130947967, 0.002948265548, 0.002959804007, 0.002697129301, 0.002903473022, 0.002869045667, 0.003031698175, 0.003634588359, 0.003086208451, 0.003203120949, 0.003092637923, 0.002672967873, 0.002738575719, 0.002763725286, 0.002485791671, 0.00294973082, 0.002929667055, 0.003251648963, 0.003135708876, 0.002704731628, 0.002625437784, 0.002637326282, 0.002359066446, 0.002255987615, 0.002265448903, 0.002420139337, 0.002924262235, 0.003350293021, 0.002951902598, 0.002461926978, 0.002655913987, 0.002527600965, 0.002327256475, 0.002300116212, 0.002474336877, 0.002739514386, 0.002810675895, 0.002664617975, 0.002577909445, 0.002402551657, 0.002241335738, 0.002078977645, 0.002098506638, 0.002261782071, 0.002340604819, 0.002765602047, 0.002390797271, 0.002148585308, 0.002367405345, 0.002224109391, 0.001910957354, 0.002059282986, 0.002030608105, 0.002143121664, 0.002335886857, 0.002960454181, 0.002227701645, 0.001910514222, 0.00190527573, 0.001734247778, 0.001685271276, 0.002121516629, 0.001922968469, 0.001818487049, 0.002148142537, 0.0021263257, 0.001657059318, 0.001654215713, 0.001605786675, 0.001515116099, 0.001592331939, 0.001566729393, 0.001568749327],
						[0.004333287312, 0.003805239886, 0.003937507735, 0.003980490598, 0.003756053426, 0.003827532969, 0.003622276933, 0.00361034607, 0.003259972434, 0.00337732856, 0.003681198312, 0.003561828604, 0.00343602649, 0.002917879079, 0.002912453456, 0.002950611747, 0.002635010927, 0.0027661938, 0.003139940126, 0.003427624549, 0.003578906624, 0.00307687328, 0.002969666323, 0.00302098847, 0.002743564708, 0.002933038244, 0.002606332896, 0.002895006413, 0.002805153446, 0.003141477769, 0.003278536375, 0.002795899555, 0.002969012724, 0.002875247054, 0.002543491766, 0.00256979023, 0.002599392839, 0.002394934838, 0.002872689712, 0.002960508812, 0.002908189163, 0.002926999021, 0.002450992361, 0.002432798607, 0.00239761976, 0.002189184483, 0.002227756109, 0.002331667699, 0.002421487646, 0.002875366694, 0.00298844369, 0.002606981108, 0.002317270052, 0.002405332647, 0.002385423734, 0.002234018339, 0.002249965789, 0.002484491958, 0.002751166494, 0.002687481463, 0.002383835896, 0.002218656927, 0.002133343147, 0.002062347094, 0.001937510594, 0.001986508049, 0.002164657603, 0.002328611956, 0.002693200054, 0.002124462034, 0.001949411141, 0.002055621468, 0.002065006753, 0.001803744266, 0.002000441696, 0.001965647569, 0.002239790169, 0.002399122799, 0.002997783046, 0.002008603922, 0.00166251154, 0.001645494544, 0.001594663869, 0.00156716494, 0.001981218007, 0.001904040506, 0.001840492434, 0.002101416757, 0.002207618824, 0.001389582413, 0.001484180669, 0.001422802845, 0.001374981132, 0.001391746185, 0.00151998509, 0.001501016057],
						[0.005263026397, 0.004724829646, 0.004796383717, 0.004776577303, 0.004551709944, 0.004633361599, 0.004387740333, 0.004428704724, 0.004035653591, 0.004206406044, 0.004681323607, 0.004279809736, 0.004112031878, 0.003545322717, 0.003561962143, 0.0036818911, 0.003274962376, 0.003352791662, 0.003427624549, 0.004713592054, 0.004412996599, 0.003691369524, 0.003495881485, 0.003740263973, 0.003380891517, 0.003446890802, 0.003165427335, 0.003558889199, 0.003415026759, 0.003849844066, 0.004096787328, 0.003443607935, 0.003573020979, 0.003413438583, 0.003179537078, 0.003132294152, 0.003151896462, 0.002894221711, 0.003341524005, 0.003697354707, 0.003778983674, 0.003473167568, 0.003193445134, 0.003074032778, 0.003001363495, 0.002813582463, 0.002728806664, 0.002756412871, 0.002729660122, 0.003574148845, 0.003742135969, 0.003261816697, 0.002881400693, 0.00295569575, 0.002771345251, 0.00263992019, 0.002789369528, 0.003027644353, 0.003339619025, 0.003391684513, 0.003009037919, 0.0027456759, 0.002538829406, 0.00246199063, 0.002241391942, 0.002453759042, 0.002658168093, 0.002796176512, 0.003167364769, 0.002702996108, 0.002418751103, 0.002507110088, 0.002455918705, 0.002200525629, 0.002357724883, 0.002368100476, 0.002607530716, 0.002848028378, 0.003700040109, 0.002599035893, 0.002003765691, 0.002031277, 0.001906474095, 0.001893518972, 0.002302300689, 0.002221521965, 0.002142081755, 0.002591435225, 0.002621234068, 0.001826813062, 0.001795320974, 0.001810262683, 0.001647121456, 0.001765786251, 0.001753677596, 0.00178581935],
						[0.006093797699, 0.005547769669, 0.00554862731, 0.005403948677, 0.005102927184, 0.005034250859, 0.004834011926, 0.004702079826, 0.004306432777, 0.004412656627, 0.005995902923, 0.005264020749, 0.004721710351, 0.004265575331, 0.004050883743, 0.004142008143, 0.003681445239, 0.003773810179, 0.003578906624, 0.004412996599, 0.006400274475, 0.004829108861, 0.004184353511, 0.004320974139, 0.003946196099, 0.003911279088, 0.003315031289, 0.003730880704, 0.003706376531, 0.003906636903, 0.005361479462, 0.004625380525, 0.004510725187, 0.004359761436, 0.003687701075, 0.003642718224, 0.003533366195, 0.003147419208, 0.003727981077, 0.003698115635, 0.0048289164, 0.004370629445, 0.003732576263, 0.003622061423, 0.003367869088, 0.003209881211, 0.002915577328, 0.003030936925, 0.003163595894, 0.003710364719, 0.004846839608, 0.004181164132, 0.003579984805, 0.003684863634, 0.003303553928, 0.002993021396, 0.003079061485, 0.003216678089, 0.003516415998, 0.003563033468, 0.004092782743, 0.003510117601, 0.003131965539, 0.002865249967, 0.002725611012, 0.002733212066, 0.002940756862, 0.002965537976, 0.003303248719, 0.003620976929, 0.003034990502, 0.003093187533, 0.002901151042, 0.002514304441, 0.002581041583, 0.002572983051, 0.002693267864, 0.003012297037, 0.003726597402, 0.00317777977, 0.002611472402, 0.00252861915, 0.002266167678, 0.002117982045, 0.002668235806, 0.002302778545, 0.00226003792, 0.002571792485, 0.002718526671, 0.002277131577, 0.002209301972, 0.002131248267, 0.001875955335, 0.00190288067, 0.001968385308, 0.00192768926],
						[0.005328759167, 0.004641512869, 0.00472970036, 0.004734314166, 0.004374288211, 0.004381039134, 0.004163941856, 0.004092405047, 0.003736776838, 0.003819775592, 0.005062114042, 0.004513935274, 0.004106976616, 0.00365789709, 0.003683104325, 0.003493402694, 0.003179269931, 0.003292673967, 0.00307687328, 0.003691369524, 0.004829108861, 0.004846353964, 0.003666264902, 0.003736210123, 0.003435652026, 0.003327434877, 0.002786963435, 0.003327749918, 0.003165502457, 0.003293778884, 0.00458887544, 0.003868415051, 0.003810806758, 0.003812135294, 0.003031235293, 0.00315819996, 0.003127765307, 0.002679971527, 0.003145197477, 0.003229609154, 0.004040609008, 0.003625382924, 0.003093686213, 0.003046358195, 0.002982998964, 0.002713429037, 0.00249017059, 0.002583465192, 0.002667415222, 0.002995624971, 0.004059070863, 0.003637699749, 0.002998394955, 0.003056517742, 0.002810718864, 0.002538625765, 0.0024774441, 0.00270785449, 0.003038493889, 0.003132510879, 0.003491058531, 0.00292253278, 0.002812505344, 0.002501828039, 0.002316819949, 0.002253495925, 0.002516935413, 0.002451923454, 0.002803321185, 0.003086465809, 0.002570900195, 0.002649884076, 0.002471890206, 0.002036089535, 0.002110086306, 0.002118461054, 0.002273008844, 0.002606396206, 0.003256030568, 0.002670710051, 0.002162432991, 0.002133878537, 0.00189527577, 0.001799982984, 0.00224771988, 0.001895963612, 0.001801342803, 0.002218445191, 0.00220851091, 0.001925310451, 0.001872333325, 0.001733466947, 0.001460972446, 0.001628559069, 0.001487323111, 0.001518940644],
						[0.004885570623, 0.004275172741, 0.004383809539, 0.004410779587, 0.004156728923, 0.004140177808, 0.003798179843, 0.003907918293, 0.003600043099, 0.003648447356, 0.004514844807, 0.004092159872, 0.003784774051, 0.003297331834, 0.003225346346, 0.003368136036, 0.002838947865, 0.002995705642, 0.002969666323, 0.003495881485, 0.004184353511, 0.003666264902, 0.004179783337, 0.003471414173, 0.003023854828, 0.003356610161, 0.002698519832, 0.003038886897, 0.003003883224, 0.003341485813, 0.003794583205, 0.003293611224, 0.003492744597, 0.003326876146, 0.002819414503, 0.002903537259, 0.002864338184, 0.002573591133, 0.002974747996, 0.003097156618, 0.003585739896, 0.003359651379, 0.002702284528, 0.00280620731, 0.002800027585, 0.002587252992, 0.00237875066, 0.002318001535, 0.002577506406, 0.003046000888, 0.003738788513, 0.003156040757, 0.002682901711, 0.002738199901, 0.002771310142, 0.002361977066, 0.00228325972, 0.00259919058, 0.00286215559, 0.002856169711, 0.002929151204, 0.002730446441, 0.002602943627, 0.002393241764, 0.002125213045, 0.002038210391, 0.002278688507, 0.002451203521, 0.002800159684, 0.002582754777, 0.002317755746, 0.002533026189, 0.002351341128, 0.001900481007, 0.002128911057, 0.002028752638, 0.002231438776, 0.002492451223, 0.003146094171, 0.002458733781, 0.001977984212, 0.002032270797, 0.001658662239, 0.001736552634, 0.002087295719, 0.001880784241, 0.001813132445, 0.00215663356, 0.002089571752, 0.00168525719, 0.001848983201, 0.001633653508, 0.00153148607, 0.001521509988, 0.001577008256, 0.001443522175],
						[0.005006898428, 0.004535712695, 0.004507738633, 0.004359942384, 0.004290280092, 0.004271141624, 0.003995557898, 0.004086110471, 0.003660395382, 0.003715759146, 0.004573083939, 0.004119086638, 0.003711100379, 0.003420793192, 0.003335041919, 0.003444248522, 0.003071383135, 0.003130947967, 0.00302098847, 0.003740263973, 0.004320974139, 0.003736210123, 0.003471414173, 0.004219814752, 0.003074347462, 0.00335558002, 0.002853020663, 0.003163075181, 0.003162460188, 0.003383226237, 0.003946854786, 0.003358245728, 0.003492581255, 0.00347055143, 0.002864913462, 0.002947059973, 0.002977872901, 0.002697692301, 0.00306046926, 0.003209536799, 0.003709842156, 0.003330488186, 0.002923999009, 0.003021807234, 0.002858419046, 0.002557295966, 0.002518513222, 0.00260775945, 0.002569833295, 0.003149431155, 0.003752861444, 0.003185902301, 0.002801078981, 0.002773726189, 0.002783920954, 0.002526022159, 0.002509305232, 0.00270304282, 0.002986993227, 0.003066809247, 0.003032034749, 0.002679082237, 0.002584795124, 0.00242750741, 0.002219893885, 0.002286300568, 0.002416626155, 0.002570982108, 0.002968969687, 0.002657204096, 0.002334412356, 0.002524391039, 0.00233510927, 0.002032861451, 0.002332175229, 0.002109761204, 0.002414349499, 0.002598142746, 0.003200476161, 0.002519373609, 0.00194844727, 0.002044466161, 0.001777307171, 0.001766203399, 0.002226324109, 0.002049160193, 0.00191435273, 0.002332399788, 0.002324302701, 0.001714808166, 0.001772016561, 0.001672395879, 0.001495796545, 0.001619156643, 0.001669654434, 0.001534071611],
						[0.004375058515, 0.004002068481, 0.004099612817, 0.004043738715, 0.003757206817, 0.003815087583, 0.00373019354, 0.003608937357, 0.003208342404, 0.003384254134, 0.004188643902, 0.003783692589, 0.00356760794, 0.003095560522, 0.003178531074, 0.003129152572, 0.002865501922, 0.002948265548, 0.002743564708, 0.003380891517, 0.003946196099, 0.003435652026, 0.003023854828, 0.003074347462, 0.003587875816, 0.002988766884, 0.002676701988, 0.002935990214, 0.002847540267, 0.003027385698, 0.003833575022, 0.003323279761, 0.003323942424, 0.003132949639, 0.002732966861, 0.002828162298, 0.002733380803, 0.002420330457, 0.002918191526, 0.002943742574, 0.003517684979, 0.003152410878, 0.002764708831, 0.002673637997, 0.002574437589, 0.002454666916, 0.002211056038, 0.002238484936, 0.002225958147, 0.002883740895, 0.00337803957, 0.003022747245, 0.002552899392, 0.002698310427, 0.002457512633, 0.0023063138, 0.002281525885, 0.002496782109, 0.002708649286, 0.002886268458, 0.002825941231, 0.002597498699, 0.002309029674, 0.002188957955, 0.002129636242, 0.002074964307, 0.002305354962, 0.002290907059, 0.002604838112, 0.002458205216, 0.002175005627, 0.002335369335, 0.002151563889, 0.001885691664, 0.002024016797, 0.001981590044, 0.002175490999, 0.002433038059, 0.002925505835, 0.002296204885, 0.001843445266, 0.001861205512, 0.001682774242, 0.001597012896, 0.002031397406, 0.001888393315, 0.001786999549, 0.002043348134, 0.002036454966, 0.001631303246, 0.001626562968, 0.001603387422, 0.00138335841, 0.001554231646, 0.00142310159, 0.001518029966],
						[0.004536227401, 0.003998855101, 0.004120392418, 0.00403985557, 0.00387302364, 0.003888776041, 0.003728689545, 0.003750278469, 0.003419811102, 0.003439829175, 0.003980118185, 0.003874680981, 0.003569918691, 0.003088608654, 0.003193442983, 0.003263735276, 0.002941029134, 0.002959804007, 0.002933038244, 0.003446890802, 0.003911279088, 0.003327434877, 0.003356610161, 0.00335558002, 0.002988766884, 0.003799319481, 0.002803502387, 0.003093195063, 0.003051601905, 0.003287937879, 0.003713693103, 0.003127171505, 0.003227872339, 0.003230790424, 0.002656555739, 0.002746432141, 0.002948685667, 0.002492862527, 0.003131702677, 0.003070495267, 0.003482263014, 0.003260850024, 0.002760310073, 0.00289823961, 0.002656191464, 0.002533025109, 0.002453741056, 0.002358155389, 0.002498554716, 0.002957772661, 0.003434297459, 0.003012799755, 0.002586533688, 0.002659368438, 0.002638276807, 0.002365156932, 0.002266873072, 0.002727073654, 0.00298885925, 0.002962162727, 0.002724667161, 0.002603957511, 0.00239513157, 0.002366851537, 0.002204672391, 0.002124099997, 0.002254771995, 0.002485193819, 0.002957245289, 0.002441650783, 0.002273750322, 0.002472290253, 0.002295598751, 0.001996875499, 0.002265404134, 0.002150671015, 0.002292409218, 0.002491206316, 0.003125602712, 0.00234600001, 0.001961811834, 0.001976809919, 0.001834543289, 0.001731717203, 0.002172197335, 0.0020091459, 0.001939627161, 0.002327275906, 0.002234689403, 0.001603869682, 0.001739008661, 0.001730536194, 0.00160289107, 0.00161558132, 0.001635188623, 0.001538851026],
						[0.00382383205, 0.003461115736, 0.003465223346, 0.003443041304, 0.003394904912, 0.003412560729, 0.00322625399, 0.003359302132, 0.002923043046, 0.003038721767, 0.003462811048, 0.003291156127, 0.003247660017, 0.002741008061, 0.002712119593, 0.002896374512, 0.002605314467, 0.002697129301, 0.002606332896, 0.003165427335, 0.003315031289, 0.002786963435, 0.002698519832, 0.002853020663, 0.002676701988, 0.002803502387, 0.00306700713, 0.002735058614, 0.002695636045, 0.002991305127, 0.00319868205, 0.002879138605, 0.002890740144, 0.00277773913, 0.002525763465, 0.002601870775, 0.002549273479, 0.002315693113, 0.002736115921, 0.00279314337, 0.002947256641, 0.002949088362, 0.002438530769, 0.002494613245, 0.002427808488, 0.0023698177, 0.002180214972, 0.002209502928, 0.002311303426, 0.002783458454, 0.002905864464, 0.002661201218, 0.002421328055, 0.00247621901, 0.002356458479, 0.002146446492, 0.002187583781, 0.002492489763, 0.002553266905, 0.002613346258, 0.002381112847, 0.002299599325, 0.002138305566, 0.002185859516, 0.002131839155, 0.002043831445, 0.002271061995, 0.00229176301, 0.00263039294, 0.00211106795, 0.002025041307, 0.002161650698, 0.002002517374, 0.001885994037, 0.002042060296, 0.001988821433, 0.002202094963, 0.002329232732, 0.002684325948, 0.002051566, 0.001725660176, 0.001823756082, 0.001687102396, 0.001621796456, 0.001933316665, 0.001898223352, 0.001877129076, 0.002004063264, 0.002049152968, 0.001490542702, 0.001573055978, 0.001556844155, 0.001397123182, 0.001458963187, 0.001493556021, 0.001504941644],
						[0.004188885336, 0.003868338392, 0.003916806442, 0.003961335731, 0.003736503159, 0.003715946693, 0.003643686193, 0.003640954706, 0.003287871386, 0.003383124035, 0.003809565319, 0.003635289158, 0.003570598902, 0.00297496778, 0.003110599618, 0.003062725198, 0.002783878304, 0.002903473022, 0.002895006413, 0.003558889199, 0.003730880704, 0.003327749918, 0.003038886897, 0.003163075181, 0.002935990214, 0.003093195063, 0.002735058614, 0.003367045221, 0.00296287914, 0.003316453212, 0.003559494386, 0.003022089541, 0.00315814943, 0.003016359076, 0.002667549281, 0.002677872445, 0.002914688297, 0.002509079189, 0.003004795826, 0.003058055976, 0.003248632596, 0.003158675801, 0.00260980252, 0.002709653171, 0.002696348858, 0.002456309315, 0.002423594261, 0.002423369531, 0.002442380546, 0.003036829971, 0.003284227606, 0.0029268648, 0.002456698341, 0.002705616909, 0.002486510885, 0.002321681675, 0.002466251737, 0.002640528653, 0.002961271301, 0.002871738576, 0.002652598675, 0.002495345023, 0.002318905424, 0.0022650826, 0.002093950813, 0.002171420609, 0.002351905388, 0.002521800233, 0.002877079814, 0.002384031672, 0.002151030603, 0.002317802626, 0.00226400016, 0.002005291199, 0.002147263607, 0.00215961835, 0.002237437268, 0.002500962852, 0.0032285676, 0.00226710001, 0.001868700425, 0.001834793041, 0.001821025008, 0.001726301118, 0.002169956628, 0.002032204841, 0.001974984732, 0.002326966181, 0.002348328682, 0.001651025562, 0.001691110869, 0.001726021726, 0.001591942106, 0.001656064673, 0.001621731298, 0.001688645496],
						[0.004227284984, 0.003743206232, 0.003893286734, 0.003888049826, 0.00370489262, 0.003694375629, 0.003498595594, 0.003612722983, 0.003183279271, 0.003306777089, 0.003914667109, 0.003635098529, 0.00352122161, 0.002989693088, 0.003033533277, 0.003016611448, 0.002705720095, 0.002869045667, 0.002805153446, 0.003415026759, 0.003706376531, 0.003165502457, 0.003003883224, 0.003162460188, 0.002847540267, 0.003051601905, 0.002695636045, 0.00296287914, 0.003230447816, 0.003234096768, 0.003466110855, 0.003037932947, 0.003089193132, 0.00300929479, 0.00257306987, 0.00266034213, 0.00273167138, 0.002420467153, 0.002915252933, 0.002969772584, 0.003252968055, 0.003136741139, 0.002526157358, 0.002702039418, 0.002609354892, 0.002335447056, 0.002329560571, 0.002357951196, 0.00246716127, 0.002996461714, 0.003250542802, 0.002808323347, 0.002478598141, 0.002593966541, 0.002468548448, 0.00228161994, 0.002341879803, 0.002538539563, 0.002828931675, 0.002913652699, 0.002624850473, 0.002433160548, 0.002332568494, 0.002191589414, 0.00211624525, 0.00214557857, 0.002334825057, 0.002464649154, 0.002830893306, 0.002365428219, 0.002215903769, 0.002292992308, 0.002186510036, 0.001928013396, 0.002127059326, 0.002026009422, 0.002247039578, 0.002468833821, 0.003034514466, 0.002220333684, 0.001885137482, 0.001871014565, 0.001755383966, 0.001768774387, 0.002103557194, 0.001956164834, 0.001946371307, 0.002220914779, 0.002292293272, 0.00153782559, 0.001667038649, 0.001664785053, 0.001486070658, 0.001637206048, 0.001584095463, 0.001656696165],
						[0.004480470259, 0.00397440804, 0.00427363146, 0.004038458787, 0.003989563848, 0.003998670233, 0.00389183741, 0.003832686285, 0.003487995237, 0.003614084622, 0.004022603493, 0.003767603955, 0.003795115997, 0.003099461678, 0.003223937921, 0.003231478037, 0.00292749477, 0.003031698175, 0.003141477769, 0.003849844066, 0.003906636903, 0.003293778884, 0.003341485813, 0.003383226237, 0.003027385698, 0.003287937879, 0.002991305127, 0.003316453212, 0.003234096768, 0.004257539481, 0.00360892735, 0.003089312147, 0.00314914342, 0.003293199865, 0.002775515786, 0.002848306622, 0.002993796676, 0.002688962135, 0.003248114478, 0.003503084115, 0.003305595067, 0.003361283934, 0.002744912459, 0.002780883803, 0.002735180879, 0.002597800854, 0.002580293557, 0.002562534438, 0.002716587516, 0.003329739, 0.003374641407, 0.002900609376, 0.002693202523, 0.002619729093, 0.002685051147, 0.00240451075, 0.002498561864, 0.002693658519, 0.003142076927, 0.003097697868, 0.00259576905, 0.002535927788, 0.002366890196, 0.002410639485, 0.002165431495, 0.002328296499, 0.00245450379, 0.002711650031, 0.003238324681, 0.002348861861, 0.002259480343, 0.002309906809, 0.002355579841, 0.002038839365, 0.002239321229, 0.002284708105, 0.002384424834, 0.002667066797, 0.003507441082, 0.002274031605, 0.001883094939, 0.001964218654, 0.00181361643, 0.001834717014, 0.002292201244, 0.002216786869, 0.002177387041, 0.002613953357, 0.002567909878, 0.001524680113, 0.001779855608, 0.001690825931, 0.001564463116, 0.001679989757, 0.001671031498, 0.001728721221],
						[0.005652719041, 0.005132938989, 0.005206899053, 0.004970339251, 0.004710173158, 0.00466125509, 0.00449380729, 0.004280689398, 0.003933100761, 0.003957043704, 0.005713294222, 0.005032640564, 0.004639492099, 0.003987789659, 0.003787378884, 0.003999979567, 0.003485988704, 0.003634588359, 0.003278536375, 0.004096787328, 0.005361479462, 0.00458887544, 0.003794583205, 0.003946854786, 0.003833575022, 0.003713693103, 0.00319868205, 0.003559494386, 0.003466110855, 0.00360892735, 0.006389629496, 0.004543286095, 0.00435339964, 0.00421183499, 0.003643843907, 0.003552145708, 0.003377381468, 0.002940465183, 0.003674445074, 0.003546939204, 0.004832200731, 0.004167746276, 0.003702399965, 0.003501187289, 0.003368893116, 0.003111353431, 0.002770718558, 0.002793733171, 0.002945230267, 0.003532627186, 0.004823624868, 0.004206178584, 0.003434773421, 0.003554204283, 0.00317215526, 0.002848030125, 0.002739330765, 0.003111718767, 0.003134480965, 0.003584368242, 0.004337202597, 0.003567835962, 0.003097652087, 0.002903839306, 0.002862427116, 0.002742707281, 0.00303943435, 0.002918479447, 0.002895614018, 0.003558485885, 0.003049905475, 0.003085204278, 0.0028736148, 0.002450573985, 0.002543276964, 0.002552567112, 0.002644966347, 0.002905540983, 0.003571788009, 0.003231873214, 0.002752491792, 0.00253094309, 0.002247082258, 0.002097977091, 0.002585262202, 0.002213443299, 0.002178675095, 0.002510933318, 0.002546967054, 0.002496554737, 0.002295747956, 0.002140649214, 0.00184150236, 0.001894217899, 0.001779130572, 0.00175253825],
						[0.004684925202, 0.004315534403, 0.004323028253, 0.004186782856, 0.003965921161, 0.003956994273, 0.003763685986, 0.003780650207, 0.003328081687, 0.003394430058, 0.004837360272, 0.00427678359, 0.003769156389, 0.003450262327, 0.003347239651, 0.003306005668, 0.00309340092, 0.003086208451, 0.002795899555, 0.003443607935, 0.004625380525, 0.003868415051, 0.003293611224, 0.003358245728, 0.003323279761, 0.003127171505, 0.002879138605, 0.003022089541, 0.003037932947, 0.003089312147, 0.004543286095, 0.004619143534, 0.003587010759, 0.003561213594, 0.003084867702, 0.003114246027, 0.002929868051, 0.002477946386, 0.003101819533, 0.003145192983, 0.004121189692, 0.003653686039, 0.003180433981, 0.003098678454, 0.002956622019, 0.002822098896, 0.002300704258, 0.002436448376, 0.002488366376, 0.002998644227, 0.003872081124, 0.003577494119, 0.003077587403, 0.003050500899, 0.002748571156, 0.002444597683, 0.002435739537, 0.002681119544, 0.002759195476, 0.003072875462, 0.003471787293, 0.00303231373, 0.002760249885, 0.002524300375, 0.002394184833, 0.002314679862, 0.002563165524, 0.002435194682, 0.002696053478, 0.003113548988, 0.002701328888, 0.00269269277, 0.002425354578, 0.002117866581, 0.002233633219, 0.002230699024, 0.002236322138, 0.002683642961, 0.002926002114, 0.002786528801, 0.002314710289, 0.002238214903, 0.001959224997, 0.001948823312, 0.002238119388, 0.001904738534, 0.001925369836, 0.002235240477, 0.002217031467, 0.002028930962, 0.002013238493, 0.001896815714, 0.001624917218, 0.001691962259, 0.001543014166, 0.001629675397],
						[0.004824278, 0.004334464152, 0.004343625848, 0.004351604312, 0.004097876261, 0.004064539484, 0.003833958135, 0.003853131458, 0.003459011173, 0.003520707926, 0.00468190493, 0.004236813741, 0.00392091161, 0.003397950778, 0.003280235516, 0.00343324518, 0.003076533575, 0.003203120949, 0.002969012724, 0.003573020979, 0.004510725187, 0.003810806758, 0.003492744597, 0.003492581255, 0.003323942424, 0.003227872339, 0.002890740144, 0.00315814943, 0.003089193132, 0.00314914342, 0.00435339964, 0.003587010759, 0.00433071444, 0.003369536117, 0.003050768156, 0.002964649715, 0.003030567194, 0.002635540325, 0.00313070882, 0.003219439668, 0.003833987856, 0.003641186078, 0.002890780361, 0.00295686402, 0.003012434425, 0.002712519187, 0.002492592931, 0.002536415273, 0.002606521851, 0.003156837627, 0.003945092908, 0.003427932124, 0.002807858163, 0.003140722747, 0.002845330975, 0.002578218747, 0.002494993998, 0.002687851375, 0.002980382082, 0.002966540382, 0.00324307985, 0.003040830314, 0.002747754518, 0.002462205185, 0.00236962802, 0.002326459033, 0.002444427051, 0.002492011945, 0.002907310129, 0.002862370156, 0.002476277853, 0.002678459174, 0.002505757031, 0.00214030932, 0.002277384155, 0.002153952799, 0.00231883584, 0.002832604928, 0.003215172404, 0.00257810375, 0.002168545284, 0.002214754088, 0.001950404018, 0.00188735876, 0.002251596033, 0.002023847758, 0.001959518691, 0.002326900403, 0.002387760302, 0.001960869109, 0.00197660144, 0.001872576902, 0.00167473677, 0.001639575669, 0.001682717651, 0.001722636547],
						[0.004694296476, 0.004243317135, 0.00432711529, 0.004236631295, 0.003968104092, 0.003977163851, 0.003889181773, 0.003750471713, 0.003390186576, 0.00340538183, 0.004517144876, 0.004152959872, 0.003720303187, 0.00336390694, 0.003345649319, 0.003199677237, 0.002884361566, 0.003092637923, 0.002875247054, 0.003413438583, 0.004359761436, 0.003812135294, 0.003326876146, 0.00347055143, 0.003132949639, 0.003230790424, 0.00277773913, 0.003016359076, 0.00300929479, 0.003293199865, 0.00421183499, 0.003561213594, 0.003369536117, 0.00399116046, 0.002839814268, 0.002955843653, 0.002934780701, 0.002577714692, 0.003126073595, 0.00308664757, 0.003817000885, 0.003431157463, 0.002981412052, 0.002907614761, 0.002769647323, 0.002650533834, 0.002431642779, 0.002420464412, 0.002619850514, 0.002904114795, 0.003772310698, 0.003295540396, 0.002766997167, 0.002774585934, 0.002778439634, 0.002491512264, 0.002358714498, 0.00263069605, 0.002811431806, 0.002873848244, 0.003285817904, 0.002866599622, 0.002514428298, 0.002443524259, 0.002274656068, 0.00221991396, 0.002444334685, 0.002513984598, 0.002864591507, 0.002849514422, 0.002484707529, 0.002510756872, 0.002307803728, 0.002045824989, 0.002175205296, 0.002223016351, 0.002255465655, 0.002509502702, 0.003086664146, 0.002559313469, 0.002125851234, 0.002042135601, 0.001865060743, 0.001717740157, 0.002234850987, 0.001997505152, 0.001918653967, 0.002221537306, 0.002176170363, 0.001858518206, 0.001821790762, 0.001711593746, 0.001589475474, 0.001632520931, 0.001637769712, 0.001580531026],
						[0.0038311798, 0.003586062754, 0.003575259214, 0.003515810312, 0.003341734566, 0.003358623301, 0.003237279705, 0.003221510183, 0.002954394167, 0.002993569578, 0.003730300128, 0.003422740202, 0.003305639622, 0.002851727225, 0.002756067678, 0.002782050118, 0.002539860536, 0.002672967873, 0.002543491766, 0.003179537078, 0.003687701075, 0.003031235293, 0.002819414503, 0.002864913462, 0.002732966861, 0.002656555739, 0.002525763465, 0.002667549281, 0.00257306987, 0.002775515786, 0.003643843907, 0.003084867702, 0.003050768156, 0.002839814268, 0.003113261711, 0.002605304699, 0.002523949341, 0.002256129179, 0.002623905379, 0.002786969656, 0.003148736152, 0.002983505359, 0.002586300099, 0.002513571863, 0.002510273881, 0.002352049805, 0.002093244918, 0.002053063802, 0.002270429913, 0.002740250691, 0.00312637474, 0.002822365865, 0.002439855205, 0.002520047139, 0.002422558995, 0.002079851963, 0.002235407413, 0.002373286904, 0.002521135822, 0.002633358523, 0.002729913204, 0.00243657735, 0.002218249685, 0.002148550824, 0.00206209881, 0.002026451581, 0.0021809998, 0.002184221162, 0.002537095063, 0.002429949287, 0.002107394936, 0.002216443786, 0.002131532637, 0.001853938575, 0.001872324134, 0.001909036257, 0.002066530273, 0.00233571516, 0.002790613686, 0.002256992943, 0.001784834378, 0.001836535011, 0.001683397023, 0.001668607442, 0.001948746689, 0.001817722133, 0.001742905187, 0.001849463061, 0.001900591909, 0.001683628808, 0.001661515879, 0.00162246157, 0.001411743762, 0.001469458031, 0.001513598958, 0.001448566505],
						[0.004131720831, 0.003673558324, 0.003722534095, 0.003803261423, 0.003518722423, 0.003518226553, 0.003384040493, 0.003375329421, 0.002954282232, 0.003161183487, 0.003878401571, 0.003607572409, 0.003329627441, 0.002891509442, 0.002767110209, 0.002870923467, 0.002578019495, 0.002738575719, 0.00256979023, 0.003132294152, 0.003642718224, 0.00315819996, 0.002903537259, 0.002947059973, 0.002828162298, 0.002746432141, 0.002601870775, 0.002677872445, 0.00266034213, 0.002848306622, 0.003552145708, 0.003114246027, 0.002964649715, 0.002955843653, 0.002605304699, 0.003246127894, 0.002597081723, 0.002289770099, 0.002681355027, 0.002807109181, 0.003251010168, 0.002906715171, 0.002690006606, 0.002534592074, 0.002472799134, 0.002380507514, 0.002121420707, 0.002153921171, 0.002199153858, 0.00274069719, 0.003199500798, 0.002891564344, 0.002373193014, 0.002513681059, 0.002319327304, 0.002084102308, 0.00224190194, 0.002394257036, 0.00256134582, 0.00252863471, 0.002803873853, 0.002419776366, 0.002180843066, 0.002103332855, 0.002042384937, 0.001968265032, 0.002240327534, 0.002214681019, 0.002462364389, 0.002424836448, 0.002049107198, 0.002083366443, 0.00197987481, 0.001875720046, 0.001887109246, 0.001949640749, 0.002115343694, 0.002307184219, 0.002697538545, 0.002179623843, 0.001839841472, 0.001812945473, 0.001639023288, 0.001532550038, 0.001956169542, 0.00177166698, 0.001680075186, 0.001911598533, 0.001958579231, 0.001715632971, 0.001600549369, 0.001576340268, 0.00143020767, 0.001363441344, 0.001368196308, 0.001363541314],
						[0.003918816918, 0.003599684906, 0.003569795438, 0.003629734042, 0.003444342914, 0.003373116109, 0.003371544744, 0.00337204551, 0.003005882414, 0.002984117306, 0.003657651758, 0.003390717311, 0.00330156276, 0.002818125281, 0.002914685529, 0.002859287353, 0.002610227175, 0.002763725286, 0.002599392839, 0.003151896462, 0.003533366195, 0.003127765307, 0.002864338184, 0.002977872901, 0.002733380803, 0.002948685667, 0.002549273479, 0.002914688297, 0.00273167138, 0.002993796676, 0.003377381468, 0.002929868051, 0.003030567194, 0.002934780701, 0.002523949341, 0.002597081723, 0.003241193577, 0.002403651445, 0.00292431403, 0.002867665129, 0.003037304181, 0.003101957178, 0.002501104067, 0.002572921285, 0.002675836699, 0.002472353739, 0.002332287125, 0.002256236355, 0.002316121493, 0.002808901081, 0.003053367321, 0.002762400551, 0.002330442759, 0.002681173546, 0.002527505865, 0.002247129268, 0.002221462671, 0.002480365542, 0.002739803121, 0.002757168903, 0.00244546032, 0.002450922864, 0.002309054272, 0.002160272368, 0.002095941077, 0.002156213221, 0.002307032847, 0.002425733908, 0.00276051117, 0.002288477239, 0.002059309865, 0.002395086863, 0.002338798784, 0.001965452817, 0.002133015049, 0.002127609696, 0.002180562054, 0.002459078114, 0.003023566438, 0.002211358705, 0.001910315612, 0.001850515018, 0.001865905876, 0.001676086598, 0.002156986478, 0.001938585802, 0.001892057444, 0.002211490898, 0.002181747593, 0.00161786628, 0.001671386618, 0.001651221385, 0.001634855501, 0.001659412043, 0.001620250311, 0.001693878156],
						[0.00363012853, 0.003328464186, 0.003341274835, 0.003345952371, 0.00311780955, 0.003176878082, 0.003103820812, 0.003078967411, 0.002767666317, 0.002825546864, 0.003235475906, 0.003054099004, 0.002915103266, 0.002423383822, 0.002613713399, 0.00266264615, 0.002350060347, 0.002485791671, 0.002394934838, 0.002894221711, 0.003147419208, 0.002679971527, 0.002573591133, 0.002697692301, 0.002420330457, 0.002492862527, 0.002315693113, 0.002509079189, 0.002420467153, 0.002688962135, 0.002940465183, 0.002477946386, 0.002635540325, 0.002577714692, 0.002256129179, 0.002289770099, 0.002403651445, 0.002585858455, 0.00253544555, 0.002491761963, 0.002702254236, 0.002568549821, 0.002301293365, 0.002340550062, 0.002103647344, 0.002006869365, 0.001983430368, 0.002017001486, 0.002145451789, 0.002531211976, 0.002653443198, 0.002399420782, 0.002150288726, 0.002176605807, 0.002159345857, 0.001969783828, 0.001983946335, 0.002136173418, 0.002299743107, 0.002245431262, 0.002152025938, 0.00202899222, 0.001868377558, 0.001880645303, 0.001805562653, 0.001860821108, 0.001986037911, 0.002040183949, 0.002325492097, 0.001988182248, 0.001757660917, 0.001900188593, 0.001913647432, 0.001699751327, 0.001726585286, 0.001847907701, 0.001886751418, 0.002194595047, 0.002489518257, 0.001822208334, 0.001510692335, 0.001534000403, 0.001554887457, 0.001469200418, 0.001768324436, 0.001737591391, 0.001611170221, 0.001913410436, 0.001930658194, 0.001297908725, 0.001312630625, 0.001291381265, 0.001203881102, 0.001290540791, 0.001357172477, 0.001286333176],
						[0.004303746538, 0.003895523432, 0.00395494306, 0.003877361097, 0.003726171349, 0.003775985543, 0.003638987516, 0.003631111805, 0.003242079448, 0.003317706912, 0.003778763194, 0.003663946352, 0.00347520935, 0.002985284814, 0.002927980931, 0.003115850747, 0.002642787593, 0.00294973082, 0.002872689712, 0.003341524005, 0.003727981077, 0.003145197477, 0.002974747996, 0.00306046926, 0.002918191526, 0.003131702677, 0.002736115921, 0.003004795826, 0.002915252933, 0.003248114478, 0.003674445074, 0.003101819533, 0.00313070882, 0.003126073595, 0.002623905379, 0.002681355027, 0.00292431403, 0.00253544555, 0.003536369603, 0.003104430334, 0.0032369198, 0.003232893531, 0.002606564187, 0.002637783682, 0.002598238767, 0.002373926728, 0.002323770748, 0.002325965117, 0.002503271213, 0.002921510826, 0.003188651794, 0.002843995701, 0.002491374206, 0.002673806675, 0.002571918771, 0.002379817153, 0.002362586563, 0.002660356607, 0.002844804395, 0.002827361159, 0.002658962851, 0.002517867529, 0.0022987776, 0.002230930061, 0.002179070837, 0.002236853026, 0.002385015081, 0.002537127192, 0.002868935179, 0.002362173787, 0.002143727368, 0.002323640337, 0.002270978414, 0.001988655607, 0.002219449236, 0.002179725325, 0.002221286989, 0.002421063393, 0.003063830896, 0.002320725035, 0.001971615854, 0.001932172136, 0.001822227792, 0.001709856006, 0.00215968941, 0.002086821575, 0.001997052557, 0.002346739526, 0.002276102439, 0.001649338559, 0.001713164575, 0.001690114187, 0.00157166069, 0.001657481338, 0.001697401999, 0.001620050265],
						[0.004208510197, 0.003933823625, 0.004053391677, 0.003927757425, 0.003807930121, 0.003852669217, 0.003674315677, 0.003705470881, 0.003357423284, 0.003465333124, 0.00390359752, 0.003692743513, 0.003564945967, 0.003095657796, 0.003180273915, 0.003073912489, 0.002800009414, 0.002929667055, 0.002960508812, 0.003697354707, 0.003698115635, 0.003229609154, 0.003097156618, 0.003209536799, 0.002943742574, 0.003070495267, 0.00279314337, 0.003058055976, 0.002969772584, 0.003503084115, 0.003546939204, 0.003145192983, 0.003219439668, 0.00308664757, 0.002786969656, 0.002807109181, 0.002867665129, 0.002491761963, 0.003104430334, 0.003862731071, 0.003180756389, 0.003252099621, 0.00273071594, 0.002610622681, 0.002728077169, 0.002511774866, 0.002401525174, 0.002438591233, 0.002548323656, 0.003072985728, 0.003316450463, 0.002800158523, 0.002509262718, 0.002616284985, 0.002480641143, 0.002301242311, 0.002420236356, 0.002662738394, 0.003062524055, 0.003058764787, 0.002719144755, 0.002530666059, 0.002268505579, 0.00224192054, 0.002096645124, 0.002253077345, 0.00241413202, 0.002519864315, 0.002964938217, 0.002489035842, 0.002126379973, 0.002259785774, 0.002215060521, 0.001987783614, 0.002186204378, 0.002162518963, 0.002429174419, 0.002663709273, 0.003367061096, 0.002292200105, 0.001857039652, 0.001928067481, 0.001756471334, 0.001761635758, 0.002084969638, 0.002156477303, 0.002065400695, 0.002442677105, 0.00243369268, 0.001676106267, 0.001761323688, 0.001715023211, 0.001546563915, 0.001547746646, 0.0016355327, 0.001723242545],
						[0.00500884843, 0.00454823586, 0.0045785797, 0.004493401784, 0.004239428375, 0.004223583686, 0.003997987643, 0.004051691901, 0.003561854758, 0.003620403799, 0.005197194723, 0.004504288936, 0.004115633012, 0.003702221795, 0.003624700079, 0.003631119945, 0.003194270936, 0.003251648963, 0.002908189163, 0.003778983674, 0.0048289164, 0.004040609008, 0.003585739896, 0.003709842156, 0.003517684979, 0.003482263014, 0.002947256641, 0.003248632596, 0.003252968055, 0.003305595067, 0.004832200731, 0.004121189692, 0.003833987856, 0.003817000885, 0.003148736152, 0.003251010168, 0.003037304181, 0.002702254236, 0.0032369198, 0.003180756389, 0.004867089328, 0.003781505789, 0.003375505154, 0.00328495968, 0.003029895698, 0.00294335548, 0.002422654082, 0.002583255508, 0.002618984466, 0.003097829045, 0.004262926982, 0.003772977814, 0.003133892725, 0.003135532692, 0.002847280697, 0.002614391589, 0.002467966524, 0.00279316539, 0.002985339315, 0.003151541273, 0.003871364107, 0.003256375384, 0.002828057206, 0.002597729241, 0.002459102341, 0.002355242093, 0.002675446452, 0.002647022174, 0.00282064572, 0.003334871745, 0.002860552143, 0.002801038174, 0.002492115849, 0.002183844354, 0.002363437762, 0.002239527586, 0.002268189282, 0.00267658867, 0.00312213084, 0.00295891475, 0.00243252121, 0.002317128886, 0.001976554273, 0.001845401642, 0.00234040907, 0.002079558915, 0.001874360924, 0.002396075867, 0.002140531648, 0.002277918698, 0.002093178431, 0.002042077396, 0.001717972348, 0.001779658493, 0.001629050738, 0.001584058446],
						[0.004497613403, 0.004173927655, 0.004162307351, 0.004028057711, 0.003945152337, 0.003867575054, 0.00374307051, 0.003758263597, 0.003334818137, 0.00336627656, 0.00435597295, 0.004016010617, 0.003837504568, 0.003339200871, 0.003335812529, 0.003343044549, 0.002960725557, 0.003135708876, 0.002926999021, 0.003473167568, 0.004370629445, 0.003625382924, 0.003359651379, 0.003330488186, 0.003152410878, 0.003260850024, 0.002949088362, 0.003158675801, 0.003136741139, 0.003361283934, 0.004167746276, 0.003653686039, 0.003641186078, 0.003431157463, 0.002983505359, 0.002906715171, 0.003101957178, 0.002568549821, 0.003232893531, 0.003252099621, 0.003781505789, 0.004323497568, 0.00299665323, 0.003046330957, 0.00307135378, 0.002827228213, 0.002427048034, 0.002554118485, 0.002790492181, 0.003161271791, 0.003790943132, 0.003395870217, 0.002897485451, 0.003176837063, 0.002927014815, 0.002613668375, 0.002481251664, 0.002773992033, 0.002958353609, 0.002879259816, 0.003210875327, 0.003026247196, 0.002743137921, 0.002572044672, 0.002498455262, 0.002438483166, 0.002639653815, 0.002611863186, 0.00311687834, 0.002884994166, 0.00258528776, 0.002802613883, 0.002626465312, 0.002368972886, 0.002417322856, 0.002343205445, 0.002340431675, 0.002674308167, 0.003246847673, 0.002728429499, 0.002414224911, 0.002314975858, 0.002088830654, 0.002000047114, 0.002459099058, 0.00217047674, 0.00203353407, 0.002300070477, 0.002407624321, 0.002036565539, 0.002083728244, 0.002002175109, 0.001875702655, 0.001848791536, 0.001887412237, 0.001908320189],
						[0.003925394458, 0.00350670524, 0.003616132532, 0.003580285672, 0.003383915216, 0.003372837448, 0.003281052417, 0.003219548501, 0.002942575524, 0.002999951399, 0.003811533041, 0.003455978013, 0.003139800334, 0.002930283919, 0.002856828165, 0.00289164188, 0.002526066962, 0.002704731628, 0.002450992361, 0.003193445134, 0.003732576263, 0.003093686213, 0.002702284528, 0.002923999009, 0.002764708831, 0.002760310073, 0.002438530769, 0.00260980252, 0.002526157358, 0.002744912459, 0.003702399965, 0.003180433981, 0.002890780361, 0.002981412052, 0.002586300099, 0.002690006606, 0.002501104067, 0.002301293365, 0.002606564187, 0.00273071594, 0.003375505154, 0.00299665323, 0.003236821869, 0.002687041388, 0.002555955639, 0.002407418934, 0.002158106374, 0.002110355545, 0.00219950667, 0.00265617112, 0.003335814834, 0.002947628681, 0.002472716885, 0.002548346256, 0.002319075618, 0.00211918581, 0.002132957306, 0.002253955311, 0.002451245485, 0.002538588326, 0.002866070977, 0.002630237528, 0.002252036354, 0.002199742604, 0.002015571044, 0.001971782777, 0.002144784117, 0.002030323869, 0.002408010983, 0.002520117346, 0.002221175975, 0.002228063953, 0.002056039595, 0.001878906302, 0.00193734829, 0.001977574331, 0.001972020328, 0.002174230852, 0.002690242371, 0.002311024277, 0.001885365619, 0.00180390305, 0.001602129373, 0.001579517001, 0.001933646242, 0.001659412079, 0.001606100416, 0.001862895769, 0.001994029146, 0.001724103763, 0.001573455973, 0.001526185446, 0.001367473141, 0.001388784454, 0.001310831026, 0.001294276419],
						[0.003826366745, 0.003414971713, 0.003500002623, 0.00341839726, 0.003340113351, 0.003276297438, 0.003160486446, 0.00322410537, 0.002935013377, 0.002888635953, 0.003647806326, 0.003492681159, 0.003134726059, 0.002769878876, 0.002905681077, 0.002883360707, 0.002600968519, 0.002625437784, 0.002432798607, 0.003074032778, 0.003622061423, 0.003046358195, 0.00280620731, 0.003021807234, 0.002673637997, 0.00289823961, 0.002494613245, 0.002709653171, 0.002702039418, 0.002780883803, 0.003501187289, 0.003098678454, 0.00295686402, 0.002907614761, 0.002513571863, 0.002534592074, 0.002572921285, 0.002340550062, 0.002637783682, 0.002610622681, 0.00328495968, 0.003046330957, 0.002687041388, 0.003172722247, 0.002703964364, 0.002380647836, 0.002275645923, 0.002110789127, 0.002225447566, 0.002748700439, 0.003147833761, 0.002897016079, 0.002525782909, 0.002544187618, 0.002423581281, 0.002150737582, 0.002264520511, 0.002413970608, 0.002546536325, 0.002455990775, 0.00277449007, 0.002477180466, 0.002369015847, 0.002247473907, 0.00215330362, 0.002032713052, 0.002175111753, 0.002230589291, 0.002510170656, 0.002475605517, 0.002287571624, 0.002359078671, 0.002171076485, 0.00200555453, 0.00202633518, 0.00214252625, 0.001987810504, 0.002343696514, 0.002688894305, 0.00223512425, 0.001907502635, 0.001900859149, 0.001811115973, 0.001725504202, 0.002020156092, 0.001792794597, 0.001736940424, 0.002063233178, 0.002008129068, 0.001679055023, 0.001674292826, 0.001687232765, 0.001571458355, 0.001584612812, 0.001583883056, 0.001552290241],
						[0.003479102582, 0.003313012908, 0.00324854751, 0.003306099986, 0.003251743842, 0.003134634966, 0.003002804477, 0.003053520724, 0.002794287977, 0.002803079025, 0.003519599805, 0.003262636438, 0.003276090681, 0.002819545433, 0.002759159271, 0.002722695521, 0.002373548708, 0.002637326282, 0.00239761976, 0.003001363495, 0.003367869088, 0.002982998964, 0.002800027585, 0.002858419046, 0.002574437589, 0.002656191464, 0.002427808488, 0.002696348858, 0.002609354892, 0.002735180879, 0.003368893116, 0.002956622019, 0.003012434425, 0.002769647323, 0.002510273881, 0.002472799134, 0.002675836699, 0.002103647344, 0.002598238767, 0.002728077169, 0.003029895698, 0.00307135378, 0.002555955639, 0.002703964364, 0.003326836059, 0.002469461262, 0.002218134497, 0.002089632921, 0.002262581838, 0.002719719359, 0.00317050846, 0.002873997188, 0.002347548213, 0.002657807666, 0.002457468729, 0.002198946237, 0.002201418296, 0.002347567932, 0.002524656713, 0.00247925764, 0.002722687972, 0.002562334361, 0.002485296261, 0.002264635978, 0.002122171843, 0.002073310735, 0.002203683292, 0.00230487603, 0.00248902096, 0.002451727614, 0.002269753022, 0.002469453739, 0.002207256383, 0.002010964229, 0.002105815658, 0.00210722721, 0.002005792894, 0.002344255713, 0.002866565122, 0.002265267561, 0.002038896069, 0.001958905783, 0.001872698697, 0.001695447636, 0.00206957304, 0.00182965049, 0.001730261696, 0.002064534697, 0.001966684595, 0.001747999998, 0.001856325133, 0.001757823844, 0.001656088474, 0.001698686823, 0.001684434505, 0.00171185129],
						[0.003217506668, 0.003017905708, 0.003014847819, 0.003046845036, 0.002902764019, 0.002941816614, 0.002773694236, 0.00282870831, 0.002568939824, 0.002559037881, 0.003250822993, 0.003004556681, 0.002948850289, 0.002524422025, 0.002564853416, 0.002568207708, 0.00233798433, 0.002359066446, 0.002189184483, 0.002813582463, 0.003209881211, 0.002713429037, 0.002587252992, 0.002557295966, 0.002454666916, 0.002533025109, 0.0023698177, 0.002456309315, 0.002335447056, 0.002597800854, 0.003111353431, 0.002822098896, 0.002712519187, 0.002650533834, 0.002352049805, 0.002380507514, 0.002472353739, 0.002006869365, 0.002373926728, 0.002511774866, 0.00294335548, 0.002827228213, 0.002407418934, 0.002380647836, 0.002469461262, 0.00273628785, 0.002050882213, 0.002019935754, 0.002077264058, 0.002482654811, 0.002798265508, 0.002574299859, 0.002235745093, 0.002374428454, 0.002322647411, 0.002029250779, 0.00195440405, 0.002153424168, 0.002335183798, 0.002326288792, 0.002413520365, 0.002282572339, 0.00210209921, 0.002039529681, 0.002037678444, 0.001915372075, 0.002023075094, 0.002057182916, 0.00232234439, 0.002193428133, 0.002100723821, 0.002191852187, 0.002013411673, 0.001817086302, 0.001956677058, 0.001922936359, 0.001926255715, 0.002134125577, 0.002562093626, 0.002117743289, 0.001801102602, 0.001774330601, 0.001647109971, 0.001577969709, 0.001900135167, 0.001711745192, 0.001660298011, 0.001848218608, 0.001831724253, 0.001566536254, 0.001659315375, 0.001530303126, 0.001421089374, 0.001542052019, 0.001504801496, 0.00145597054],
						[0.003092539967, 0.002835420141, 0.002964497291, 0.002896100679, 0.002782632546, 0.002837975836, 0.002668950943, 0.002747014377, 0.002487098668, 0.002522217349, 0.002851594993, 0.002798426204, 0.002638145651, 0.002300083628, 0.002352480996, 0.002351255139, 0.002110325803, 0.002255987615, 0.002227756109, 0.002728806664, 0.002915577328, 0.00249017059, 0.00237875066, 0.002518513222, 0.002211056038, 0.002453741056, 0.002180214972, 0.002423594261, 0.002329560571, 0.002580293557, 0.002770718558, 0.002300704258, 0.002492592931, 0.002431642779, 0.002093244918, 0.002121420707, 0.002332287125, 0.001983430368, 0.002323770748, 0.002401525174, 0.002422654082, 0.002427048034, 0.002158106374, 0.002275645923, 0.002218134497, 0.002050882213, 0.002292387056, 0.001877801099, 0.001959666183, 0.002394099937, 0.002641668335, 0.002242256865, 0.001981723824, 0.002117546945, 0.002090233351, 0.001864608629, 0.001948433257, 0.002073591892, 0.00233437301, 0.002300474881, 0.002038243146, 0.002012237245, 0.001985023885, 0.001852543988, 0.001776673135, 0.001841625137, 0.00188748672, 0.001934718146, 0.002294160273, 0.001867289702, 0.001752422095, 0.001966253632, 0.001836744414, 0.00163411733, 0.001779050864, 0.001877644921, 0.001820894591, 0.002145394252, 0.002643949387, 0.001795900975, 0.001530803387, 0.001533680315, 0.001455476451, 0.001490387539, 0.001746800454, 0.001655796305, 0.001665854964, 0.001962286642, 0.001977779456, 0.001304139577, 0.001386751283, 0.001319070256, 0.001340951321, 0.00135730855, 0.001383488268, 0.001367615045],
						[0.003387487433, 0.00304936103, 0.003052088947, 0.00306826418, 0.002923215582, 0.002966669732, 0.002792739793, 0.002873165543, 0.002538668543, 0.002613303501, 0.003133866747, 0.002855994378, 0.002842945493, 0.002393712985, 0.002384606717, 0.002444673168, 0.002179840156, 0.002265448903, 0.002331667699, 0.002756412871, 0.003030936925, 0.002583465192, 0.002318001535, 0.00260775945, 0.002238484936, 0.002358155389, 0.002209502928, 0.002423369531, 0.002357951196, 0.002562534438, 0.002793733171, 0.002436448376, 0.002536415273, 0.002420464412, 0.002053063802, 0.002153921171, 0.002256236355, 0.002017001486, 0.002325965117, 0.002438591233, 0.002583255508, 0.002554118485, 0.002110355545, 0.002110789127, 0.002089632921, 0.002019935754, 0.001877801099, 0.002255582718, 0.002068414027, 0.002404775849, 0.002593694069, 0.002293415019, 0.001993140333, 0.002173903779, 0.002054921614, 0.001893534002, 0.001927413423, 0.002076444179, 0.002311205946, 0.002274966194, 0.002120337891, 0.001995792556, 0.001938888603, 0.001855814132, 0.001722785106, 0.001832567826, 0.001856676975, 0.002010200721, 0.002260962916, 0.001973444816, 0.001803630977, 0.001894799551, 0.001817529201, 0.001639164356, 0.00178368941, 0.001739598415, 0.001934769853, 0.002056383569, 0.002448597599, 0.001835817712, 0.001554517899, 0.001528020481, 0.001494005094, 0.001475744626, 0.001733468366, 0.001701345065, 0.001608682258, 0.001935430794, 0.001909424209, 0.001308641757, 0.001388652916, 0.001359905258, 0.001246940377, 0.001404625531, 0.001371462861, 0.001372454729],
						[0.003395715859, 0.003091043629, 0.003176130571, 0.003125160756, 0.003007653059, 0.00306358095, 0.0028610154, 0.002847918425, 0.002648547527, 0.002654684242, 0.003116543132, 0.002949423992, 0.003000883638, 0.002549847288, 0.002531822685, 0.002568367833, 0.002177178217, 0.002420139337, 0.002421487646, 0.002729660122, 0.003163595894, 0.002667415222, 0.002577506406, 0.002569833295, 0.002225958147, 0.002498554716, 0.002311303426, 0.002442380546, 0.00246716127, 0.002716587516, 0.002945230267, 0.002488366376, 0.002606521851, 0.002619850514, 0.002270429913, 0.002199153858, 0.002316121493, 0.002145451789, 0.002503271213, 0.002548323656, 0.002618984466, 0.002790492181, 0.00219950667, 0.002225447566, 0.002262581838, 0.002077264058, 0.001959666183, 0.002068414027, 0.002620710447, 0.002603254346, 0.002640121329, 0.002404387829, 0.002117499327, 0.00231447102, 0.002240947952, 0.002012175463, 0.001954114917, 0.002216503475, 0.002286007278, 0.002303013543, 0.002179407599, 0.002098648833, 0.002033042901, 0.001928023816, 0.001856704911, 0.001840407715, 0.002022820615, 0.002108084727, 0.002461211042, 0.002013526607, 0.001899826248, 0.00206525987, 0.00192052398, 0.001766285481, 0.001905770945, 0.00183325678, 0.001966713259, 0.002198692526, 0.002631403393, 0.001867927429, 0.001659350044, 0.001623427959, 0.001643625112, 0.001512696411, 0.001903089234, 0.001771103512, 0.001686312902, 0.001977935665, 0.001958590849, 0.001336237643, 0.001430388361, 0.001420022489, 0.001365864643, 0.00148274205, 0.001519772785, 0.001519431052],
						[0.004000762357, 0.00371479407, 0.003789134736, 0.003719892798, 0.003567505347, 0.003607251209, 0.003470426784, 0.003448010781, 0.003152359262, 0.003261252232, 0.003691230013, 0.003455392502, 0.003455188678, 0.002984198221, 0.002957736965, 0.003064479536, 0.002647144259, 0.002924262235, 0.002875366694, 0.003574148845, 0.003710364719, 0.002995624971, 0.003046000888, 0.003149431155, 0.002883740895, 0.002957772661, 0.002783458454, 0.003036829971, 0.002996461714, 0.003329739, 0.003532627186, 0.002998644227, 0.003156837627, 0.002904114795, 0.002740250691, 0.00274069719, 0.002808901081, 0.002531211976, 0.002921510826, 0.003072985728, 0.003097829045, 0.003161271791, 0.00265617112, 0.002748700439, 0.002719719359, 0.002482654811, 0.002394099937, 0.002404775849, 0.002603254346, 0.003507874446, 0.003261877243, 0.002892954422, 0.002485275888, 0.002805093691, 0.002547038788, 0.002374874245, 0.00248914718, 0.002706362721, 0.00290711498, 0.002931542597, 0.002570636166, 0.002504877642, 0.00234094259, 0.002328788354, 0.002196628596, 0.002275081176, 0.002412950967, 0.002613051957, 0.00293017241, 0.002370038322, 0.002156137142, 0.002405012624, 0.002390016422, 0.002112982315, 0.002240209112, 0.002241507991, 0.002418662726, 0.002577853902, 0.003224991981, 0.002298418854, 0.001922046187, 0.001935690547, 0.001934885145, 0.001863184141, 0.002231743966, 0.002133414729, 0.002105693351, 0.002335114657, 0.002455373533, 0.00165968462, 0.001753961274, 0.001731624546, 0.0017056357, 0.001694330807, 0.001752820866, 0.001814234354],
						[0.004921900893, 0.004518149563, 0.004571681602, 0.004420189562, 0.004212690159, 0.004047749529, 0.003914053963, 0.003931554582, 0.003521043502, 0.003583525581, 0.004951209144, 0.004442682175, 0.00405233259, 0.003592822144, 0.003428895014, 0.003569656968, 0.002943346366, 0.003350293021, 0.00298844369, 0.003742135969, 0.004846839608, 0.004059070863, 0.003738788513, 0.003752861444, 0.00337803957, 0.003434297459, 0.002905864464, 0.003284227606, 0.003250542802, 0.003374641407, 0.004823624868, 0.003872081124, 0.003945092908, 0.003772310698, 0.00312637474, 0.003199500798, 0.003053367321, 0.002653443198, 0.003188651794, 0.003316450463, 0.004262926982, 0.003790943132, 0.003335814834, 0.003147833761, 0.00317050846, 0.002798265508, 0.002641668335, 0.002593694069, 0.002640121329, 0.003261877243, 0.005037470377, 0.003785148817, 0.003021184311, 0.003144722717, 0.002826997808, 0.002607740669, 0.002574420361, 0.002691890067, 0.003095029334, 0.003187072618, 0.003814743727, 0.003287373498, 0.002838411383, 0.002605977434, 0.002375888774, 0.002473344017, 0.002512097295, 0.002746438852, 0.002892256126, 0.003257937869, 0.002711491509, 0.002887254079, 0.002549283435, 0.002204302719, 0.00230289056, 0.002175761723, 0.002439105183, 0.002668846883, 0.003424183921, 0.002966632996, 0.002436055038, 0.002355746166, 0.002025684896, 0.001889482664, 0.002335763509, 0.002037849418, 0.002005744056, 0.002458236182, 0.002416220994, 0.002355888284, 0.002182515245, 0.002002751956, 0.001770018254, 0.001689529892, 0.001693807964, 0.001633407093],
						[0.004177801985, 0.003944750725, 0.003864248166, 0.003735997127, 0.003554662776, 0.003518512219, 0.003418586327, 0.003313178608, 0.003050485455, 0.003009493202, 0.004333963728, 0.003868226071, 0.003689967188, 0.00315543466, 0.003150561594, 0.003122008957, 0.002783759962, 0.002951902598, 0.002606981108, 0.003261816697, 0.004181164132, 0.003637699749, 0.003156040757, 0.003185902301, 0.003022747245, 0.003012799755, 0.002661201218, 0.0029268648, 0.002808323347, 0.002900609376, 0.004206178584, 0.003577494119, 0.003427932124, 0.003295540396, 0.002822365865, 0.002891564344, 0.002762400551, 0.002399420782, 0.002843995701, 0.002800158523, 0.003772977814, 0.003395870217, 0.002947628681, 0.002897016079, 0.002873997188, 0.002574299859, 0.002242256865, 0.002293415019, 0.002404387829, 0.002892954422, 0.003785148817, 0.003795038486, 0.00274344516, 0.002846294208, 0.002643327318, 0.002373651119, 0.002450064418, 0.002525335517, 0.002638444729, 0.002732705094, 0.003402279329, 0.002927043135, 0.002580653227, 0.002380309535, 0.002335352466, 0.002231337129, 0.0024369971, 0.002456676516, 0.002660254894, 0.002915525926, 0.002531438179, 0.002613211472, 0.002426121036, 0.002052740998, 0.002158481331, 0.0021140707, 0.002201615384, 0.00245250601, 0.002904037504, 0.002652139789, 0.002216415693, 0.002137932195, 0.001988579044, 0.001773262991, 0.002244430763, 0.001871709159, 0.001750546507, 0.002199654503, 0.002098263899, 0.002135683637, 0.002002634763, 0.001942345413, 0.001702335031, 0.001709835884, 0.001653003911, 0.001587142632],
						[0.003654297707, 0.003314415963, 0.003371414581, 0.003210599149, 0.003178857495, 0.003220634491, 0.003024243577, 0.003028926938, 0.00263296343, 0.002728739572, 0.003606982422, 0.003236383866, 0.002976667966, 0.002677606164, 0.00263577638, 0.00274180904, 0.00247288644, 0.002461926978, 0.002317270052, 0.002881400693, 0.003579984805, 0.002998394955, 0.002682901711, 0.002801078981, 0.002552899392, 0.002586533688, 0.002421328055, 0.002456698341, 0.002478598141, 0.002693202523, 0.003434773421, 0.003077587403, 0.002807858163, 0.002766997167, 0.002439855205, 0.002373193014, 0.002330442759, 0.002150288726, 0.002491374206, 0.002509262718, 0.003133892725, 0.002897485451, 0.002472716885, 0.002525782909, 0.002347548213, 0.002235745093, 0.001981723824, 0.001993140333, 0.002117499327, 0.002485275888, 0.003021184311, 0.00274344516, 0.002907841859, 0.002360326866, 0.002201432646, 0.00197900765, 0.002028213155, 0.002182706382, 0.002323935946, 0.002457778017, 0.002650656099, 0.00229871417, 0.002193479256, 0.002106481104, 0.001988621986, 0.001907748692, 0.002067338217, 0.002054173723, 0.002371889688, 0.002318416621, 0.002158656476, 0.002257302626, 0.001995026144, 0.001761774637, 0.001884896651, 0.001870519581, 0.001795517347, 0.002097052813, 0.002458630052, 0.00220955019, 0.001801943668, 0.001859650234, 0.00158182962, 0.00155812396, 0.001802007614, 0.00162487014, 0.001567712893, 0.001873815258, 0.001786576186, 0.001529337994, 0.001611606875, 0.001470394496, 0.001221941677, 0.001405332314, 0.001389306134, 0.001394237652],
						[0.003598637829, 0.00343389213, 0.003366216157, 0.003358889658, 0.003158288587, 0.003225088129, 0.003035802487, 0.003032511577, 0.002732898344, 0.00276611704, 0.003729182959, 0.003332541527, 0.00317508292, 0.002848600657, 0.002772655784, 0.002790751101, 0.002423331136, 0.002655913987, 0.002405332647, 0.00295569575, 0.003684863634, 0.003056517742, 0.002738199901, 0.002773726189, 0.002698310427, 0.002659368438, 0.00247621901, 0.002705616909, 0.002593966541, 0.002619729093, 0.003554204283, 0.003050500899, 0.003140722747, 0.002774585934, 0.002520047139, 0.002513681059, 0.002681173546, 0.002176605807, 0.002673806675, 0.002616284985, 0.003135532692, 0.003176837063, 0.002548346256, 0.002544187618, 0.002657807666, 0.002374428454, 0.002117546945, 0.002173903779, 0.00231447102, 0.002805093691, 0.003144722717, 0.002846294208, 0.002360326866, 0.003149523875, 0.002415837636, 0.002200086391, 0.002162371878, 0.002324971658, 0.002461091217, 0.002459801104, 0.00268770369, 0.00263931478, 0.002375925292, 0.002159224955, 0.002075736788, 0.002131126522, 0.002320906714, 0.002170903329, 0.002468893541, 0.002472894677, 0.002205323436, 0.002454136129, 0.002299884068, 0.00205816416, 0.002082048619, 0.002101663354, 0.00213781531, 0.002254171392, 0.002674545428, 0.002328507026, 0.002047874052, 0.001979971247, 0.001877199216, 0.001696418522, 0.002103916017, 0.001831609016, 0.001815691681, 0.002105651673, 0.002104577206, 0.001764384174, 0.001754955104, 0.001726516932, 0.00164620921, 0.00164119908, 0.001632717802, 0.00170226126],
						[0.00347358152, 0.003213881416, 0.003231814936, 0.003149259674, 0.003031807835, 0.003074806453, 0.002865340173, 0.002911464346, 0.002639358885, 0.00257947995, 0.003338004342, 0.003093094158, 0.003002616169, 0.002637211277, 0.00257720391, 0.00260473893, 0.002272073889, 0.002527600965, 0.002385423734, 0.002771345251, 0.003303553928, 0.002810718864, 0.002771310142, 0.002783920954, 0.002457512633, 0.002638276807, 0.002356458479, 0.002486510885, 0.002468548448, 0.002685051147, 0.00317215526, 0.002748571156, 0.002845330975, 0.002778439634, 0.002422558995, 0.002319327304, 0.002527505865, 0.002159345857, 0.002571918771, 0.002480641143, 0.002847280697, 0.002927014815, 0.002319075618, 0.002423581281, 0.002457468729, 0.002322647411, 0.002090233351, 0.002054921614, 0.002240947952, 0.002547038788, 0.002826997808, 0.002643327318, 0.002201432646, 0.002415837636, 0.002767315247, 0.002043651733, 0.002034993437, 0.002204876043, 0.00234345839, 0.002419818446, 0.002411671324, 0.002288762897, 0.002288017767, 0.002072733425, 0.002008139257, 0.001982204687, 0.002128501004, 0.002193425371, 0.002539332658, 0.002208228764, 0.002025280393, 0.002246620993, 0.002133769892, 0.001897569735, 0.002027384031, 0.001938074558, 0.001924333968, 0.002256134662, 0.002770495738, 0.00210142722, 0.001832840001, 0.001837111486, 0.001635391534, 0.001627968138, 0.002071366733, 0.00171055106, 0.001684877891, 0.001950284629, 0.002041913208, 0.0014921953, 0.001613936507, 0.001541106991, 0.001471749075, 0.001556784332, 0.001587092448, 0.001462726316],
						[0.003243638023, 0.003034389219, 0.003016512653, 0.002974762136, 0.002869155025, 0.002869951923, 0.002746671583, 0.002745567374, 0.00245577162, 0.002479049424, 0.003044495193, 0.002851603943, 0.002716058382, 0.002388155994, 0.002387000641, 0.002401464726, 0.002157631265, 0.002327256475, 0.002234018339, 0.00263992019, 0.002993021396, 0.002538625765, 0.002361977066, 0.002526022159, 0.0023063138, 0.002365156932, 0.002146446492, 0.002321681675, 0.00228161994, 0.00240451075, 0.002848030125, 0.002444597683, 0.002578218747, 0.002491512264, 0.002079851963, 0.002084102308, 0.002247129268, 0.001969783828, 0.002379817153, 0.002301242311, 0.002614391589, 0.002613668375, 0.00211918581, 0.002150737582, 0.002198946237, 0.002029250779, 0.001864608629, 0.001893534002, 0.002012175463, 0.002374874245, 0.002607740669, 0.002373651119, 0.00197900765, 0.002200086391, 0.002043651733, 0.002258678707, 0.001848273189, 0.002082300372, 0.002138891732, 0.002093305463, 0.002110123955, 0.002090690778, 0.001923976682, 0.001807177388, 0.001846959803, 0.001777795269, 0.001862867801, 0.002014353121, 0.002167454362, 0.00190608691, 0.001789749635, 0.001991950908, 0.001904169139, 0.001637374176, 0.00181615008, 0.001761109708, 0.001886317215, 0.001973870058, 0.002307852321, 0.001854531318, 0.001563540612, 0.001555407534, 0.001528674383, 0.001389151699, 0.001813067579, 0.00165444263, 0.001623472287, 0.00179482723, 0.001723097603, 0.001368238167, 0.001434811152, 0.001338967163, 0.001364436889, 0.001395907181, 0.001442945758, 0.001410140658],
						[0.003030619633, 0.002929474526, 0.002922320365, 0.002953982871, 0.002756593414, 0.002794455413, 0.002713610852, 0.002658174324, 0.00240993133, 0.002497650293, 0.002987117227, 0.002816757742, 0.002865345579, 0.002352838588, 0.002357687016, 0.002339342052, 0.002150761998, 0.002300116212, 0.002249965789, 0.002789369528, 0.003079061485, 0.0024774441, 0.00228325972, 0.002509305232, 0.002281525885, 0.002266873072, 0.002187583781, 0.002466251737, 0.002341879803, 0.002498561864, 0.002739330765, 0.002435739537, 0.002494993998, 0.002358714498, 0.002235407413, 0.00224190194, 0.002221462671, 0.001983946335, 0.002362586563, 0.002420236356, 0.002467966524, 0.002481251664, 0.002132957306, 0.002264520511, 0.002201418296, 0.00195440405, 0.001948433257, 0.001927413423, 0.001954114917, 0.00248914718, 0.002574420361, 0.002450064418, 0.002028213155, 0.002162371878, 0.002034993437, 0.001848273189, 0.002538786963, 0.00222124262, 0.002253509776, 0.002253188285, 0.002299911894, 0.0020225735, 0.001944254201, 0.001863789412, 0.001819545462, 0.001845145489, 0.001922706237, 0.002048916468, 0.00233515705, 0.002032495395, 0.001825373988, 0.001888469975, 0.001866935206, 0.001704091112, 0.001790801096, 0.001876568938, 0.001794484491, 0.002072880281, 0.00253783427, 0.001899653245, 0.001507957861, 0.001604302022, 0.001605942258, 0.001543079407, 0.001825533626, 0.001679398022, 0.001668945698, 0.001919123145, 0.001938805242, 0.001453504753, 0.001421291285, 0.001527167595, 0.001413216884, 0.001448010909, 0.001478431473, 0.001485237838],
						[0.003440298944, 0.00334420473, 0.003274787097, 0.003297600332, 0.003070416077, 0.00314857206, 0.002952582409, 0.003029978077, 0.002745503964, 0.00278677727, 0.003204243619, 0.003150937277, 0.002911171213, 0.002637722498, 0.002638124, 0.002678514626, 0.002429601622, 0.002474336877, 0.002484491958, 0.003027644353, 0.003216678089, 0.00270785449, 0.00259919058, 0.00270304282, 0.002496782109, 0.002727073654, 0.002492489763, 0.002640528653, 0.002538539563, 0.002693658519, 0.003111718767, 0.002681119544, 0.002687851375, 0.00263069605, 0.002373286904, 0.002394257036, 0.002480365542, 0.002136173418, 0.002660356607, 0.002662738394, 0.00279316539, 0.002773992033, 0.002253955311, 0.002413970608, 0.002347567932, 0.002153424168, 0.002073591892, 0.002076444179, 0.002216503475, 0.002706362721, 0.002691890067, 0.002525335517, 0.002182706382, 0.002324971658, 0.002204876043, 0.002082300372, 0.00222124262, 0.002888158837, 0.00246768224, 0.002545452725, 0.002354012316, 0.002171205185, 0.002126558897, 0.002073529808, 0.002061390344, 0.002037461461, 0.00216609843, 0.00225437478, 0.002504191509, 0.002126676968, 0.001914553843, 0.002159707946, 0.002024030435, 0.001879695348, 0.00199321726, 0.002005570746, 0.002129976208, 0.002232949935, 0.002727763283, 0.002024908208, 0.001716670092, 0.001777444928, 0.001822485902, 0.001672337281, 0.00202609524, 0.001939476904, 0.001907641035, 0.002103609887, 0.002081524935, 0.001497014391, 0.001545671162, 0.001640709803, 0.001567264754, 0.001563318738, 0.001658344814, 0.001638997416],
						[0.003931812779, 0.003452367593, 0.003685314298, 0.003627668097, 0.003440003388, 0.003445206259, 0.003384751458, 0.003456979638, 0.003105609008, 0.003201191769, 0.003555873242, 0.003371518839, 0.0033207877, 0.00285052315, 0.002886483586, 0.002819510233, 0.00258608262, 0.002739514386, 0.002751166494, 0.003339619025, 0.003516415998, 0.003038493889, 0.00286215559, 0.002986993227, 0.002708649286, 0.00298885925, 0.002553266905, 0.002961271301, 0.002828931675, 0.003142076927, 0.003134480965, 0.002759195476, 0.002980382082, 0.002811431806, 0.002521135822, 0.00256134582, 0.002739803121, 0.002299743107, 0.002844804395, 0.003062524055, 0.002985339315, 0.002958353609, 0.002451245485, 0.002546536325, 0.002524656713, 0.002335183798, 0.00233437301, 0.002311205946, 0.002286007278, 0.00290711498, 0.003095029334, 0.002638444729, 0.002323935946, 0.002461091217, 0.00234345839, 0.002138891732, 0.002253509776, 0.00246768224, 0.003350583833, 0.002772584344, 0.002545890912, 0.00228588627, 0.002134542835, 0.002105589144, 0.001972762346, 0.002065197164, 0.002134455168, 0.002349684212, 0.002885162927, 0.002390285755, 0.001977518945, 0.002179275575, 0.002159351368, 0.001854111527, 0.002045975845, 0.002050252042, 0.002150942185, 0.002416861729, 0.00308709855, 0.002117514202, 0.001769438899, 0.001741036688, 0.001707515433, 0.001644695146, 0.00200538994, 0.001967585467, 0.001912712957, 0.002272632476, 0.002188558276, 0.001539832924, 0.001644256099, 0.001611616989, 0.001492172521, 0.001583504098, 0.001624889665, 0.001665134061],
						[0.004004083997, 0.003555839981, 0.003723157506, 0.003716680489, 0.003495357866, 0.003484276388, 0.003398783072, 0.003418816348, 0.003080811196, 0.003171901734, 0.003764693424, 0.003421802213, 0.003383236234, 0.002949664632, 0.002862711554, 0.002970506193, 0.002629567419, 0.002810675895, 0.002687481463, 0.003391684513, 0.003563033468, 0.003132510879, 0.002856169711, 0.003066809247, 0.002886268458, 0.002962162727, 0.002613346258, 0.002871738576, 0.002913652699, 0.003097697868, 0.003584368242, 0.003072875462, 0.002966540382, 0.002873848244, 0.002633358523, 0.00252863471, 0.002757168903, 0.002245431262, 0.002827361159, 0.003058764787, 0.003151541273, 0.002879259816, 0.002538588326, 0.002455990775, 0.00247925764, 0.002326288792, 0.002300474881, 0.002274966194, 0.002303013543, 0.002931542597, 0.003187072618, 0.002732705094, 0.002457778017, 0.002459801104, 0.002419818446, 0.002093305463, 0.002253188285, 0.002545452725, 0.002772584344, 0.004079755781, 0.002537002164, 0.002413232926, 0.002321551262, 0.002189272422, 0.002113071136, 0.002099741406, 0.002304073073, 0.002393506739, 0.002693485986, 0.002321912481, 0.002075753427, 0.002323698689, 0.00217477644, 0.00173595849, 0.002037136712, 0.002029942012, 0.002329128467, 0.002497822871, 0.002997682812, 0.002208828426, 0.0018639162, 0.001807356548, 0.001639376038, 0.001659224599, 0.002125126885, 0.001869063588, 0.00195522756, 0.002256425015, 0.002297224771, 0.001523107793, 0.001675158753, 0.001569474671, 0.001440842522, 0.001625692915, 0.001511484483, 0.00162087977],
						[0.00379606579, 0.003597330441, 0.003694396668, 0.003460502947, 0.003285121985, 0.003242065918, 0.003164531892, 0.003089516401, 0.002765467253, 0.002852639794, 0.004334541133, 0.003810412801, 0.003420739267, 0.003039306391, 0.002969202312, 0.002821127856, 0.002489536208, 0.002664617975, 0.002383835896, 0.003009037919, 0.004092782743, 0.003491058531, 0.002929151204, 0.003032034749, 0.002825941231, 0.002724667161, 0.002381112847, 0.002652598675, 0.002624850473, 0.00259576905, 0.004337202597, 0.003471787293, 0.00324307985, 0.003285817904, 0.002729913204, 0.002803873853, 0.00244546032, 0.002152025938, 0.002658962851, 0.002719144755, 0.003871364107, 0.003210875327, 0.002866070977, 0.00277449007, 0.002722687972, 0.002413520365, 0.002038243146, 0.002120337891, 0.002179407599, 0.002570636166, 0.003814743727, 0.003402279329, 0.002650656099, 0.00268770369, 0.002411671324, 0.002110123955, 0.002299911894, 0.002354012316, 0.002545890912, 0.002537002164, 0.004115982628, 0.00285063903, 0.002468721806, 0.002265342652, 0.002177720248, 0.002024486222, 0.002346885696, 0.002171199726, 0.002353119122, 0.003107641537, 0.002544776951, 0.002392017713, 0.0021142444, 0.002042435891, 0.001953488988, 0.00201758249, 0.002020718018, 0.002274620913, 0.00277005577, 0.002698335303, 0.002290962724, 0.002012465449, 0.00185434367, 0.001648297176, 0.002036498587, 0.001800994314, 0.001637370346, 0.001989093522, 0.001834252413, 0.00224975273, 0.001895068104, 0.001898527721, 0.001539200166, 0.001481003002, 0.001471806275, 0.001417212846],
						[0.00338809299, 0.003228763006, 0.003228902136, 0.003201833917, 0.003009518249, 0.00296896478, 0.002893749342, 0.002823881199, 0.002570024624, 0.002556324433, 0.003553878345, 0.003295890373, 0.003050888097, 0.002737754808, 0.002639004767, 0.002624749646, 0.002321161174, 0.002577909445, 0.002218656927, 0.0027456759, 0.003510117601, 0.00292253278, 0.002730446441, 0.002679082237, 0.002597498699, 0.002603957511, 0.002299599325, 0.002495345023, 0.002433160548, 0.002535927788, 0.003567835962, 0.00303231373, 0.003040830314, 0.002866599622, 0.00243657735, 0.002419776366, 0.002450922864, 0.00202899222, 0.002517867529, 0.002530666059, 0.003256375384, 0.003026247196, 0.002630237528, 0.002477180466, 0.002562334361, 0.002282572339, 0.002012237245, 0.001995792556, 0.002098648833, 0.002504877642, 0.003287373498, 0.002927043135, 0.00229871417, 0.00263931478, 0.002288762897, 0.002090690778, 0.0020225735, 0.002171205185, 0.00228588627, 0.002413232926, 0.00285063903, 0.002928884545, 0.002343570116, 0.002190788495, 0.002009566027, 0.001967943614, 0.002072171205, 0.002006980689, 0.002390876113, 0.002517607928, 0.002260702104, 0.002359948953, 0.002099432758, 0.001857530485, 0.001919568355, 0.001893251837, 0.001914514415, 0.002093109033, 0.002478238425, 0.002328326893, 0.002022264844, 0.001956531438, 0.00171365323, 0.001608642565, 0.001983211094, 0.001738711241, 0.001653203216, 0.001914985238, 0.001870729067, 0.00181912707, 0.001775680814, 0.001709574839, 0.001544124909, 0.001523394415, 0.001502743129, 0.001442843411],
						[0.003208075825, 0.002824046695, 0.002950115286, 0.002927183154, 0.002794968285, 0.002867250372, 0.002576146487, 0.002694148822, 0.002403539014, 0.002395752497, 0.003231203484, 0.002979283605, 0.002846153878, 0.002598694665, 0.002423556742, 0.002501434852, 0.00216214058, 0.002402551657, 0.002133343147, 0.002538829406, 0.003131965539, 0.002812505344, 0.002602943627, 0.002584795124, 0.002309029674, 0.00239513157, 0.002138305566, 0.002318905424, 0.002332568494, 0.002366890196, 0.003097652087, 0.002760249885, 0.002747754518, 0.002514428298, 0.002218249685, 0.002180843066, 0.002309054272, 0.001868377558, 0.0022987776, 0.002268505579, 0.002828057206, 0.002743137921, 0.002252036354, 0.002369015847, 0.002485296261, 0.00210209921, 0.001985023885, 0.001938888603, 0.002033042901, 0.00234094259, 0.002838411383, 0.002580653227, 0.002193479256, 0.002375925292, 0.002288017767, 0.001923976682, 0.001944254201, 0.002126558897, 0.002134542835, 0.002321551262, 0.002468721806, 0.002343570116, 0.002667728836, 0.002059014354, 0.001875570377, 0.001815770394, 0.001926379175, 0.001910367164, 0.002259947412, 0.002210713674, 0.002105304897, 0.00227225667, 0.002006248791, 0.001678812607, 0.00184701393, 0.00179998138, 0.001704258851, 0.002055336166, 0.002590616688, 0.002073480187, 0.001874480495, 0.001817398768, 0.00158261204, 0.00160578189, 0.001834642637, 0.001570663998, 0.001592933608, 0.002057573856, 0.001761267425, 0.001508379831, 0.001597609687, 0.001539238963, 0.001439716463, 0.001489612646, 0.001449216874, 0.001435943811],
						[0.002935523324, 0.002719042358, 0.002751010898, 0.002807933105, 0.002665138781, 0.002673288533, 0.002534087762, 0.002581815013, 0.002365928298, 0.002348819196, 0.002943099406, 0.002720910432, 0.002658557789, 0.002320043343, 0.002302659107, 0.002245179492, 0.002093913888, 0.002241335738, 0.002062347094, 0.00246199063, 0.002865249967, 0.002501828039, 0.002393241764, 0.00242750741, 0.002188957955, 0.002366851537, 0.002185859516, 0.0022650826, 0.002191589414, 0.002410639485, 0.002903839306, 0.002524300375, 0.002462205185, 0.002443524259, 0.002148550824, 0.002103332855, 0.002160272368, 0.001880645303, 0.002230930061, 0.00224192054, 0.002597729241, 0.002572044672, 0.002199742604, 0.002247473907, 0.002264635978, 0.002039529681, 0.001852543988, 0.001855814132, 0.001928023816, 0.002328788354, 0.002605977434, 0.002380309535, 0.002106481104, 0.002159224955, 0.002072733425, 0.001807177388, 0.001863789412, 0.002073529808, 0.002105589144, 0.002189272422, 0.002265342652, 0.002190788495, 0.002059014354, 0.002239974825, 0.001860978705, 0.001763013832, 0.001851172919, 0.001813111348, 0.002145930135, 0.002033279861, 0.001984916688, 0.002110346607, 0.001896770696, 0.001675683985, 0.001751742541, 0.001808930336, 0.001737045766, 0.001859046805, 0.002208657034, 0.001916434512, 0.001698026684, 0.001711879501, 0.001552914235, 0.001560586611, 0.001794451241, 0.001666964874, 0.001556821933, 0.001761721129, 0.001674310159, 0.001338224906, 0.001469942223, 0.001372897973, 0.001312547653, 0.001355054415, 0.001392112944, 0.001311553054],
						[0.002579952426, 0.002595143005, 0.002432745554, 0.002393663894, 0.002399572175, 0.002378738807, 0.002298139155, 0.002274152315, 0.002060583994, 0.002024066215, 0.002696793804, 0.002541390204, 0.00248201665, 0.002149391811, 0.002211500634, 0.002192825626, 0.001967676888, 0.002078977645, 0.001937510594, 0.002241391942, 0.002725611012, 0.002316819949, 0.002125213045, 0.002219893885, 0.002129636242, 0.002204672391, 0.002131839155, 0.002093950813, 0.00211624525, 0.002165431495, 0.002862427116, 0.002394184833, 0.00236962802, 0.002274656068, 0.00206209881, 0.002042384937, 0.002095941077, 0.001805562653, 0.002179070837, 0.002096645124, 0.002459102341, 0.002498455262, 0.002015571044, 0.00215330362, 0.002122171843, 0.002037678444, 0.001776673135, 0.001722785106, 0.001856704911, 0.002196628596, 0.002375888774, 0.002335352466, 0.001988621986, 0.002075736788, 0.002008139257, 0.001846959803, 0.001819545462, 0.002061390344, 0.001972762346, 0.002113071136, 0.002177720248, 0.002009566027, 0.001875570377, 0.001860978705, 0.00217574472, 0.001695211453, 0.001905841881, 0.001938044055, 0.001955343717, 0.001947203413, 0.001835683003, 0.00199579587, 0.001878570542, 0.001637405318, 0.001789838079, 0.001715748504, 0.001829044759, 0.001879697817, 0.00210540579, 0.001800355523, 0.0016475462, 0.001614763026, 0.001592853435, 0.001466974956, 0.001730962735, 0.001611239802, 0.001517674255, 0.001656282469, 0.00164398893, 0.001406303153, 0.001505474531, 0.001409622401, 0.001346738712, 0.001380300963, 0.001465920283, 0.001417331467],
						[0.002747787374, 0.00264760704, 0.002705575637, 0.002612167028, 0.002482747036, 0.002476341968, 0.002369198989, 0.002452890052, 0.002128802882, 0.002145725437, 0.002640047186, 0.002527176121, 0.002494794761, 0.002141692279, 0.002164797252, 0.002164558022, 0.001976377986, 0.002098506638, 0.001986508049, 0.002453759042, 0.002733212066, 0.002253495925, 0.002038210391, 0.002286300568, 0.002074964307, 0.002124099997, 0.002043831445, 0.002171420609, 0.00214557857, 0.002328296499, 0.002742707281, 0.002314679862, 0.002326459033, 0.00221991396, 0.002026451581, 0.001968265032, 0.002156213221, 0.001860821108, 0.002236853026, 0.002253077345, 0.002355242093, 0.002438483166, 0.001971782777, 0.002032713052, 0.002073310735, 0.001915372075, 0.001841625137, 0.001832567826, 0.001840407715, 0.002275081176, 0.002473344017, 0.002231337129, 0.001907748692, 0.002131126522, 0.001982204687, 0.001777795269, 0.001845145489, 0.002037461461, 0.002065197164, 0.002099741406, 0.002024486222, 0.001967943614, 0.001815770394, 0.001763013832, 0.001695211453, 0.002226843776, 0.001831680806, 0.001941214445, 0.002092000333, 0.001879618498, 0.001749958111, 0.001914182171, 0.001800106007, 0.001640827494, 0.001703409595, 0.00178709596, 0.001800410688, 0.001956771025, 0.002342896959, 0.001864830712, 0.001542859377, 0.00157997703, 0.001570629667, 0.001497027997, 0.001760054456, 0.001638850111, 0.001596252795, 0.001828888462, 0.001783426135, 0.001414236235, 0.001439126627, 0.001421393745, 0.001257311173, 0.001417501875, 0.001434609122, 0.001413057928],
						[0.002993945136, 0.002929607457, 0.002880348639, 0.002891811709, 0.002766082028, 0.00278931335, 0.002650016671, 0.002639344627, 0.002378697448, 0.002443764267, 0.003154011213, 0.002777588485, 0.002796370631, 0.002455331697, 0.002527220919, 0.002394731624, 0.00212531014, 0.002261782071, 0.002164657603, 0.002658168093, 0.002940756862, 0.002516935413, 0.002278688507, 0.002416626155, 0.002305354962, 0.002254771995, 0.002271061995, 0.002351905388, 0.002334825057, 0.00245450379, 0.00303943435, 0.002563165524, 0.002444427051, 0.002444334685, 0.0021809998, 0.002240327534, 0.002307032847, 0.001986037911, 0.002385015081, 0.00241413202, 0.002675446452, 0.002639653815, 0.002144784117, 0.002175111753, 0.002203683292, 0.002023075094, 0.00188748672, 0.001856676975, 0.002022820615, 0.002412950967, 0.002512097295, 0.0024369971, 0.002067338217, 0.002320906714, 0.002128501004, 0.001862867801, 0.001922706237, 0.00216609843, 0.002134455168, 0.002304073073, 0.002346885696, 0.002072171205, 0.001926379175, 0.001851172919, 0.001905841881, 0.001831680806, 0.002362971986, 0.002013099909, 0.002211800006, 0.002037217978, 0.001839230725, 0.001984006171, 0.001901105125, 0.001751261424, 0.001834798411, 0.001852889709, 0.001891929825, 0.002045778129, 0.002521603449, 0.001912894773, 0.001712256401, 0.001621937332, 0.001611890165, 0.001470146558, 0.001863226423, 0.001699113779, 0.001664173048, 0.001824116362, 0.001925407783, 0.001432696674, 0.001447393032, 0.001488967379, 0.001353134112, 0.001441163181, 0.001423381553, 0.001498832196],
						[0.003157948431, 0.003064194774, 0.002975663554, 0.002789078706, 0.002807242675, 0.002744695175, 0.002596581078, 0.002777070666, 0.002375537224, 0.002445315495, 0.0030178547, 0.00282129975, 0.002825190452, 0.002372667954, 0.002339467906, 0.002495798699, 0.002052453022, 0.002340604819, 0.002328611956, 0.002796176512, 0.002965537976, 0.002451923454, 0.002451203521, 0.002570982108, 0.002290907059, 0.002485193819, 0.00229176301, 0.002521800233, 0.002464649154, 0.002711650031, 0.002918479447, 0.002435194682, 0.002492011945, 0.002513984598, 0.002184221162, 0.002214681019, 0.002425733908, 0.002040183949, 0.002537127192, 0.002519864315, 0.002647022174, 0.002611863186, 0.002030323869, 0.002230589291, 0.00230487603, 0.002057182916, 0.001934718146, 0.002010200721, 0.002108084727, 0.002613051957, 0.002746438852, 0.002456676516, 0.002054173723, 0.002170903329, 0.002193425371, 0.002014353121, 0.002048916468, 0.00225437478, 0.002349684212, 0.002393506739, 0.002171199726, 0.002006980689, 0.001910367164, 0.001813111348, 0.001938044055, 0.001941214445, 0.002013099909, 0.002913229606, 0.002342612304, 0.002031034186, 0.001806357272, 0.001986177557, 0.001948085404, 0.001738131697, 0.001961482042, 0.001842175714, 0.002085970703, 0.002257239072, 0.002693510647, 0.001937579107, 0.001624053584, 0.001595164909, 0.001618440704, 0.001438808775, 0.001889528393, 0.001777879072, 0.001792062365, 0.00216390251, 0.001970141906, 0.001485800171, 0.001609694131, 0.001519896088, 0.001452871784, 0.001512179815, 0.001595801621, 0.00157436436],
						[0.003715079122, 0.003378001969, 0.003502498329, 0.00345071268, 0.003305060105, 0.003419811664, 0.003333011815, 0.003271405217, 0.002979944128, 0.002970631164, 0.003233047781, 0.003289166972, 0.003221488961, 0.002748106064, 0.002896492621, 0.002804875205, 0.002483985727, 0.002765602047, 0.002693200054, 0.003167364769, 0.003303248719, 0.002803321185, 0.002800159684, 0.002968969687, 0.002604838112, 0.002957245289, 0.00263039294, 0.002877079814, 0.002830893306, 0.003238324681, 0.002895614018, 0.002696053478, 0.002907310129, 0.002864591507, 0.002537095063, 0.002462364389, 0.00276051117, 0.002325492097, 0.002868935179, 0.002964938217, 0.00282064572, 0.00311687834, 0.002408010983, 0.002510170656, 0.00248902096, 0.00232234439, 0.002294160273, 0.002260962916, 0.002461211042, 0.00293017241, 0.002892256126, 0.002660254894, 0.002371889688, 0.002468893541, 0.002539332658, 0.002167454362, 0.00233515705, 0.002504191509, 0.002885162927, 0.002693485986, 0.002353119122, 0.002390876113, 0.002259947412, 0.002145930135, 0.001955343717, 0.002092000333, 0.002211800006, 0.002342612304, 0.003658817031, 0.00212495448, 0.001982733737, 0.002274283814, 0.002126321301, 0.001906339992, 0.002098460066, 0.002175184895, 0.002123685539, 0.002434289261, 0.003287733483, 0.002143422865, 0.001787322454, 0.001853956393, 0.001761340961, 0.001761927511, 0.002223767524, 0.00202875252, 0.002015049614, 0.002359241942, 0.002349290251, 0.001502287893, 0.001690429149, 0.001701696033, 0.001637735267, 0.001694576331, 0.001727829921, 0.001714672068],
						[0.003337110603, 0.003206598288, 0.00322837047, 0.003015226454, 0.002907239952, 0.002865320479, 0.002746955289, 0.002731006516, 0.002466248967, 0.002471499072, 0.003705899965, 0.003287562718, 0.002953910851, 0.002744602577, 0.002636127797, 0.002463355151, 0.00221991072, 0.002390797271, 0.002124462034, 0.002702996108, 0.003620976929, 0.003086465809, 0.002582754777, 0.002657204096, 0.002458205216, 0.002441650783, 0.00211106795, 0.002384031672, 0.002365428219, 0.002348861861, 0.003558485885, 0.003113548988, 0.002862370156, 0.002849514422, 0.002429949287, 0.002424836448, 0.002288477239, 0.001988182248, 0.002362173787, 0.002489035842, 0.003334871745, 0.002884994166, 0.002520117346, 0.002475605517, 0.002451727614, 0.002193428133, 0.001867289702, 0.001973444816, 0.002013526607, 0.002370038322, 0.003257937869, 0.002915525926, 0.002318416621, 0.002472894677, 0.002208228764, 0.00190608691, 0.002032495395, 0.002126676968, 0.002390285755, 0.002321912481, 0.003107641537, 0.002517607928, 0.002210713674, 0.002033279861, 0.001947203413, 0.001879618498, 0.002037217978, 0.002031034186, 0.00212495448, 0.002986835437, 0.002274080231, 0.002230011678, 0.002045158491, 0.001799474408, 0.001782427378, 0.00181265829, 0.001825797102, 0.002054788709, 0.002457843811, 0.002349393317, 0.00199023081, 0.001857663227, 0.001724921252, 0.001551595006, 0.001880673805, 0.001614511409, 0.001502270361, 0.001976858326, 0.001714422288, 0.00193244399, 0.00173933189, 0.001693475392, 0.0014430677, 0.001442626283, 0.001374532828, 0.001388198613],
						[0.002848833667, 0.002600477361, 0.002724240891, 0.002634119186, 0.002521817319, 0.00253830332, 0.0024088781, 0.002417733311, 0.002183305605, 0.00217780669, 0.003057049188, 0.00280330556, 0.002647777308, 0.002357861268, 0.002371647666, 0.002245615832, 0.002045275133, 0.002148585308, 0.001949411141, 0.002418751103, 0.003034990502, 0.002570900195, 0.002317755746, 0.002334412356, 0.002175005627, 0.002273750322, 0.002025041307, 0.002151030603, 0.002215903769, 0.002259480343, 0.003049905475, 0.002701328888, 0.002476277853, 0.002484707529, 0.002107394936, 0.002049107198, 0.002059309865, 0.001757660917, 0.002143727368, 0.002126379973, 0.002860552143, 0.00258528776, 0.002221175975, 0.002287571624, 0.002269753022, 0.002100723821, 0.001752422095, 0.001803630977, 0.001899826248, 0.002156137142, 0.002711491509, 0.002531438179, 0.002158656476, 0.002205323436, 0.002025280393, 0.001789749635, 0.001825373988, 0.001914553843, 0.001977518945, 0.002075753427, 0.002544776951, 0.002260702104, 0.002105304897, 0.001984916688, 0.001835683003, 0.001749958111, 0.001839230725, 0.001806357272, 0.001982733737, 0.002274080231, 0.00229457775, 0.002117357162, 0.001845375086, 0.001633187581, 0.001749012712, 0.001681082922, 0.001634979053, 0.00184047917, 0.002192893711, 0.002051811528, 0.001802743262, 0.001717172202, 0.001572679132, 0.001451319933, 0.001753953548, 0.001512238654, 0.001449994849, 0.001779310054, 0.001611715347, 0.001544583527, 0.00155642601, 0.001497240486, 0.001291858116, 0.001423933613, 0.001357856058, 0.001340842036],
						[0.002949900901, 0.002858297419, 0.002751144754, 0.002712130997, 0.002658128393, 0.002672007559, 0.002490168387, 0.002545041644, 0.00233221606, 0.002264720672, 0.003098922292, 0.002900054774, 0.00267266428, 0.002469088971, 0.002451257932, 0.002447919046, 0.002103865601, 0.002367405345, 0.002055621468, 0.002507110088, 0.003093187533, 0.002649884076, 0.002533026189, 0.002524391039, 0.002335369335, 0.002472290253, 0.002161650698, 0.002317802626, 0.002292992308, 0.002309906809, 0.003085204278, 0.00269269277, 0.002678459174, 0.002510756872, 0.002216443786, 0.002083366443, 0.002395086863, 0.001900188593, 0.002323640337, 0.002259785774, 0.002801038174, 0.002802613883, 0.002228063953, 0.002359078671, 0.002469453739, 0.002191852187, 0.001966253632, 0.001894799551, 0.00206525987, 0.002405012624, 0.002887254079, 0.002613211472, 0.002257302626, 0.002454136129, 0.002246620993, 0.001991950908, 0.001888469975, 0.002159707946, 0.002179275575, 0.002323698689, 0.002392017713, 0.002359948953, 0.00227225667, 0.002110346607, 0.00199579587, 0.001914182171, 0.001984006171, 0.001986177557, 0.002274283814, 0.002230011678, 0.002117357162, 0.002621401631, 0.002126901078, 0.001825069475, 0.002014683763, 0.001882770839, 0.001953755803, 0.002056603229, 0.00240569299, 0.002109080628, 0.001936116943, 0.001858294563, 0.001781608032, 0.001627779459, 0.001970454534, 0.001773160458, 0.001636927628, 0.001953381381, 0.00179774555, 0.001620992162, 0.001691735074, 0.001597226743, 0.0015912288, 0.001611998182, 0.001645564968, 0.00160054969],
						[0.002982481082, 0.002724947082, 0.002720678945, 0.002692079815, 0.002581899568, 0.002546935189, 0.002476504628, 0.002445071567, 0.002247460908, 0.002206492379, 0.002903858749, 0.00263994758, 0.002675630159, 0.002266383371, 0.002351451581, 0.002278303381, 0.00200061205, 0.002224109391, 0.002065006753, 0.002455918705, 0.002901151042, 0.002471890206, 0.002351341128, 0.00233510927, 0.002151563889, 0.002295598751, 0.002002517374, 0.00226400016, 0.002186510036, 0.002355579841, 0.0028736148, 0.002425354578, 0.002505757031, 0.002307803728, 0.002131532637, 0.00197987481, 0.002338798784, 0.001913647432, 0.002270978414, 0.002215060521, 0.002492115849, 0.002626465312, 0.002056039595, 0.002171076485, 0.002207256383, 0.002013411673, 0.001836744414, 0.001817529201, 0.00192052398, 0.002390016422, 0.002549283435, 0.002426121036, 0.001995026144, 0.002299884068, 0.002133769892, 0.001904169139, 0.001866935206, 0.002024030435, 0.002159351368, 0.00217477644, 0.0021142444, 0.002099432758, 0.002006248791, 0.001896770696, 0.001878570542, 0.001800106007, 0.001901105125, 0.001948085404, 0.002126321301, 0.002045158491, 0.001845375086, 0.002126901078, 0.002276511839, 0.001678201675, 0.001833440004, 0.001792319511, 0.001811427404, 0.001902093582, 0.002396847716, 0.001904328088, 0.001718775095, 0.001651573889, 0.001642356721, 0.001562520092, 0.00190059034, 0.001684161069, 0.001583373631, 0.001914626162, 0.001781562354, 0.001410880098, 0.001494793329, 0.001431770042, 0.001430295881, 0.001462799528, 0.001531357426, 0.001467325825],
						[0.002356035502, 0.002440528653, 0.00232456853, 0.00220023004, 0.002227657241, 0.002227989655, 0.002163291487, 0.002144393999, 0.00188411888, 0.001912792792, 0.00246176927, 0.002390474875, 0.002214654882, 0.001915814607, 0.001999883551, 0.001985647385, 0.001752657322, 0.001910957354, 0.001803744266, 0.002200525629, 0.002514304441, 0.002036089535, 0.001900481007, 0.002032861451, 0.001885691664, 0.001996875499, 0.001885994037, 0.002005291199, 0.001928013396, 0.002038839365, 0.002450573985, 0.002117866581, 0.00214030932, 0.002045824989, 0.001853938575, 0.001875720046, 0.001965452817, 0.001699751327, 0.001988655607, 0.001987783614, 0.002183844354, 0.002368972886, 0.001878906302, 0.00200555453, 0.002010964229, 0.001817086302, 0.00163411733, 0.001639164356, 0.001766285481, 0.002112982315, 0.002204302719, 0.002052740998, 0.001761774637, 0.00205816416, 0.001897569735, 0.001637374176, 0.001704091112, 0.001879695348, 0.001854111527, 0.00173595849, 0.002042435891, 0.001857530485, 0.001678812607, 0.001675683985, 0.001637405318, 0.001640827494, 0.001751261424, 0.001738131697, 0.001906339992, 0.001799474408, 0.001633187581, 0.001825069475, 0.001678201675, 0.001901617019, 0.001689389323, 0.001681978974, 0.001659577292, 0.001815411915, 0.002077210542, 0.001680707469, 0.001479750176, 0.001444464685, 0.001531887011, 0.001363783562, 0.001697949351, 0.001514333937, 0.001467214757, 0.001568259676, 0.001741996881, 0.001342113321, 0.001277767174, 0.001379396749, 0.001353480725, 0.001308016213, 0.001373687721, 0.001352817337],
						[0.002656567017, 0.002583794193, 0.00253747154, 0.002432468817, 0.002470283393, 0.00248874816, 0.002345133043, 0.002389854799, 0.002077564129, 0.002093213831, 0.002599881832, 0.002500464799, 0.002358356459, 0.002137951166, 0.002172617411, 0.002228671992, 0.001914621865, 0.002059282986, 0.002000441696, 0.002357724883, 0.002581041583, 0.002110086306, 0.002128911057, 0.002332175229, 0.002024016797, 0.002265404134, 0.002042060296, 0.002147263607, 0.002127059326, 0.002239321229, 0.002543276964, 0.002233633219, 0.002277384155, 0.002175205296, 0.001872324134, 0.001887109246, 0.002133015049, 0.001726585286, 0.002219449236, 0.002186204378, 0.002363437762, 0.002417322856, 0.00193734829, 0.00202633518, 0.002105815658, 0.001956677058, 0.001779050864, 0.00178368941, 0.001905770945, 0.002240209112, 0.00230289056, 0.002158481331, 0.001884896651, 0.002082048619, 0.002027384031, 0.00181615008, 0.001790801096, 0.00199321726, 0.002045975845, 0.002037136712, 0.001953488988, 0.001919568355, 0.00184701393, 0.001751742541, 0.001789838079, 0.001703409595, 0.001834798411, 0.001961482042, 0.002098460066, 0.001782427378, 0.001749012712, 0.002014683763, 0.001833440004, 0.001689389323, 0.00208699488, 0.001712450073, 0.001813296546, 0.001918209543, 0.002270302748, 0.001759219865, 0.001542083706, 0.001524916872, 0.001519306473, 0.001411680124, 0.00173158236, 0.001670207893, 0.001543731544, 0.001747236502, 0.001804869635, 0.001302159289, 0.001408221757, 0.001359346906, 0.001375020272, 0.001435210999, 0.001468558918, 0.001440521518],
						[0.002544831954, 0.002526042736, 0.002549888275, 0.002479192684, 0.002403485387, 0.002472291663, 0.002378529183, 0.002344599062, 0.002123604214, 0.00213031204, 0.002592969402, 0.002486099168, 0.002414975934, 0.002032070219, 0.002162700965, 0.002114276384, 0.001904786973, 0.002030608105, 0.001965647569, 0.002368100476, 0.002572983051, 0.002118461054, 0.002028752638, 0.002109761204, 0.001981590044, 0.002150671015, 0.001988821433, 0.00215961835, 0.002026009422, 0.002284708105, 0.002552567112, 0.002230699024, 0.002153952799, 0.002223016351, 0.001909036257, 0.001949640749, 0.002127609696, 0.001847907701, 0.002179725325, 0.002162518963, 0.002239527586, 0.002343205445, 0.001977574331, 0.00214252625, 0.00210722721, 0.001922936359, 0.001877644921, 0.001739598415, 0.00183325678, 0.002241507991, 0.002175761723, 0.0021140707, 0.001870519581, 0.002101663354, 0.001938074558, 0.001761109708, 0.001876568938, 0.002005570746, 0.002050252042, 0.002029942012, 0.00201758249, 0.001893251837, 0.00179998138, 0.001808930336, 0.001715748504, 0.00178709596, 0.001852889709, 0.001842175714, 0.002175184895, 0.00181265829, 0.001681082922, 0.001882770839, 0.001792319511, 0.001681978974, 0.001712450073, 0.002175490603, 0.001743625395, 0.00199189188, 0.00226151043, 0.001711165998, 0.001540525879, 0.001519947769, 0.001537543115, 0.00147923634, 0.001773098228, 0.00170197913, 0.001674688296, 0.001831954788, 0.001804142466, 0.001318801308, 0.001365465744, 0.001388263461, 0.00141299199, 0.001423903284, 0.001504110078, 0.001507816599],
						[0.002931137112, 0.002966731622, 0.002757467982, 0.002699152581, 0.002673820064, 0.002701562387, 0.002524251707, 0.002619809282, 0.002368484015, 0.00237844092, 0.002729180908, 0.002786945071, 0.002506934412, 0.002162047775, 0.002268925825, 0.002285185976, 0.00196953169, 0.002143121664, 0.002239790169, 0.002607530716, 0.002693267864, 0.002273008844, 0.002231438776, 0.002414349499, 0.002175490999, 0.002292409218, 0.002202094963, 0.002237437268, 0.002247039578, 0.002384424834, 0.002644966347, 0.002236322138, 0.00231883584, 0.002255465655, 0.002066530273, 0.002115343694, 0.002180562054, 0.001886751418, 0.002221286989, 0.002429174419, 0.002268189282, 0.002340431675, 0.001972020328, 0.001987810504, 0.002005792894, 0.001926255715, 0.001820894591, 0.001934769853, 0.001966713259, 0.002418662726, 0.002439105183, 0.002201615384, 0.001795517347, 0.00213781531, 0.001924333968, 0.001886317215, 0.001794484491, 0.002129976208, 0.002150942185, 0.002329128467, 0.002020718018, 0.001914514415, 0.001704258851, 0.001737045766, 0.001829044759, 0.001800410688, 0.001891929825, 0.002085970703, 0.002123685539, 0.001825797102, 0.001634979053, 0.001953755803, 0.001811427404, 0.001659577292, 0.001813296546, 0.001743625395, 0.002560726752, 0.002083073614, 0.002325012524, 0.001761420977, 0.001542287283, 0.00150474773, 0.001612232676, 0.001417474375, 0.001773750525, 0.001785032259, 0.001751504465, 0.001892191381, 0.001904908508, 0.001385364582, 0.001454037065, 0.001371315787, 0.001356353647, 0.001369128949, 0.001459390899, 0.001517580729],
						[0.003179676958, 0.003089654503, 0.003098168353, 0.002967319487, 0.002835663113, 0.002937912347, 0.002792157251, 0.002852917885, 0.002484697379, 0.002584543946, 0.003147871729, 0.002980381819, 0.002874428286, 0.002401755863, 0.002563441418, 0.002486016887, 0.002308013374, 0.002335886857, 0.002399122799, 0.002848028378, 0.003012297037, 0.002606396206, 0.002492451223, 0.002598142746, 0.002433038059, 0.002491206316, 0.002329232732, 0.002500962852, 0.002468833821, 0.002667066797, 0.002905540983, 0.002683642961, 0.002832604928, 0.002509502702, 0.00233571516, 0.002307184219, 0.002459078114, 0.002194595047, 0.002421063393, 0.002663709273, 0.00267658867, 0.002674308167, 0.002174230852, 0.002343696514, 0.002344255713, 0.002134125577, 0.002145394252, 0.002056383569, 0.002198692526, 0.002577853902, 0.002668846883, 0.00245250601, 0.002097052813, 0.002254171392, 0.002256134662, 0.001973870058, 0.002072880281, 0.002232949935, 0.002416861729, 0.002497822871, 0.002274620913, 0.002093109033, 0.002055336166, 0.001859046805, 0.001879697817, 0.001956771025, 0.002045778129, 0.002257239072, 0.002434289261, 0.002054788709, 0.00184047917, 0.002056603229, 0.001902093582, 0.001815411915, 0.001918209543, 0.00199189188, 0.002083073614, 0.002980570262, 0.002776905698, 0.001883144273, 0.001638469279, 0.001646822523, 0.001668124132, 0.001633452578, 0.001904942811, 0.00181302171, 0.001790193869, 0.002071404173, 0.002108539649, 0.001503615307, 0.001489054779, 0.001519445517, 0.001495595492, 0.001462829871, 0.001541539547, 0.001597668889],
						[0.00410448382, 0.003539376119, 0.003999211129, 0.003811419029, 0.003627688117, 0.003746374369, 0.003547327024, 0.003539833913, 0.003263079647, 0.003406758634, 0.003568750149, 0.003553381874, 0.003548025327, 0.003083914722, 0.003174658432, 0.003060032788, 0.002725449108, 0.002960454181, 0.002997783046, 0.003700040109, 0.003726597402, 0.003256030568, 0.003146094171, 0.003200476161, 0.002925505835, 0.003125602712, 0.002684325948, 0.0032285676, 0.003034514466, 0.003507441082, 0.003571788009, 0.002926002114, 0.003215172404, 0.003086664146, 0.002790613686, 0.002697538545, 0.003023566438, 0.002489518257, 0.003063830896, 0.003367061096, 0.00312213084, 0.003246847673, 0.002690242371, 0.002688894305, 0.002866565122, 0.002562093626, 0.002643949387, 0.002448597599, 0.002631403393, 0.003224991981, 0.003424183921, 0.002904037504, 0.002458630052, 0.002674545428, 0.002770495738, 0.002307852321, 0.00253783427, 0.002727763283, 0.00308709855, 0.002997682812, 0.00277005577, 0.002478238425, 0.002590616688, 0.002208657034, 0.00210540579, 0.002342896959, 0.002521603449, 0.002693510647, 0.003287733483, 0.002457843811, 0.002192893711, 0.00240569299, 0.002396847716, 0.002077210542, 0.002270302748, 0.00226151043, 0.002325012524, 0.002776905698, 0.004495284329, 0.002342090592, 0.00196423921, 0.001968095145, 0.001963940144, 0.001932027363, 0.002236250352, 0.002092354005, 0.002225072321, 0.002734847721, 0.002739584572, 0.001727048959, 0.001749568877, 0.001809098989, 0.001689155895, 0.001757563505, 0.001734222474, 0.001758184786],
						[0.003094311644, 0.002988376321, 0.002958428789, 0.002813450776, 0.002699749688, 0.002725652958, 0.002582676068, 0.002638401879, 0.002315163065, 0.002316135476, 0.003247957188, 0.002975489791, 0.002676243201, 0.002498508202, 0.002401444852, 0.002355436976, 0.002072097471, 0.002227701645, 0.002008603922, 0.002599035893, 0.00317777977, 0.002670710051, 0.002458733781, 0.002519373609, 0.002296204885, 0.00234600001, 0.002051566, 0.00226710001, 0.002220333684, 0.002274031605, 0.003231873214, 0.002786528801, 0.00257810375, 0.002559313469, 0.002256992943, 0.002179623843, 0.002211358705, 0.001822208334, 0.002320725035, 0.002292200105, 0.00295891475, 0.002728429499, 0.002311024277, 0.00223512425, 0.002265267561, 0.002117743289, 0.001795900975, 0.001835817712, 0.001867927429, 0.002298418854, 0.002966632996, 0.002652139789, 0.00220955019, 0.002328507026, 0.00210142722, 0.001854531318, 0.001899653245, 0.002024908208, 0.002117514202, 0.002208828426, 0.002698335303, 0.002328326893, 0.002073480187, 0.001916434512, 0.001800355523, 0.001864830712, 0.001912894773, 0.001937579107, 0.002143422865, 0.002349393317, 0.002051811528, 0.002109080628, 0.001904328088, 0.001680707469, 0.001759219865, 0.001711165998, 0.001761420977, 0.001883144273, 0.002342090592, 0.00238108428, 0.00184448747, 0.001772725629, 0.001556963915, 0.001482927742, 0.001793079798, 0.001564612875, 0.001518886766, 0.001774290456, 0.001573060605, 0.001739458821, 0.001649662703, 0.001562606651, 0.001332975416, 0.00145065916, 0.001354462311, 0.00131376431],
						[0.002412252584, 0.002327715382, 0.002268515803, 0.002223566, 0.002162000371, 0.002111424879, 0.001961667199, 0.002021217408, 0.001846701912, 0.001841431471, 0.002593234347, 0.002443998787, 0.00227159661, 0.002029350921, 0.001947953835, 0.001903773391, 0.00164384809, 0.001910514222, 0.00166251154, 0.002003765691, 0.002611472402, 0.002162432991, 0.001977984212, 0.00194844727, 0.001843445266, 0.001961811834, 0.001725660176, 0.001868700425, 0.001885137482, 0.001883094939, 0.002752491792, 0.002314710289, 0.002168545284, 0.002125851234, 0.001784834378, 0.001839841472, 0.001910315612, 0.001510692335, 0.001971615854, 0.001857039652, 0.00243252121, 0.002414224911, 0.001885365619, 0.001907502635, 0.002038896069, 0.001801102602, 0.001530803387, 0.001554517899, 0.001659350044, 0.001922046187, 0.002436055038, 0.002216415693, 0.001801943668, 0.002047874052, 0.001832840001, 0.001563540612, 0.001507957861, 0.001716670092, 0.001769438899, 0.0018639162, 0.002290962724, 0.002022264844, 0.001874480495, 0.001698026684, 0.0016475462, 0.001542859377, 0.001712256401, 0.001624053584, 0.001787322454, 0.00199023081, 0.001802743262, 0.001936116943, 0.001718775095, 0.001479750176, 0.001542083706, 0.001540525879, 0.001542287283, 0.001638469279, 0.00196423921, 0.00184448747, 0.00192243273, 0.001573558801, 0.001467854587, 0.001348935717, 0.001656951629, 0.001443941498, 0.001344356872, 0.001572495244, 0.001432723256, 0.001472780098, 0.001464959819, 0.00137748601, 0.001316196758, 0.001285898014, 0.001326462412, 0.001307264791],
						[0.002356759221, 0.0023059398, 0.002255255898, 0.002173940325, 0.002126479925, 0.002183330207, 0.00198535973, 0.002049355105, 0.00181648464, 0.001807023938, 0.002508700781, 0.002360782909, 0.002203276806, 0.002015332046, 0.001925131204, 0.001912987547, 0.001742680996, 0.00190527573, 0.001645494544, 0.002031277, 0.00252861915, 0.002133878537, 0.002032270797, 0.002044466161, 0.001861205512, 0.001976809919, 0.001823756082, 0.001834793041, 0.001871014565, 0.001964218654, 0.00253094309, 0.002238214903, 0.002214754088, 0.002042135601, 0.001836535011, 0.001812945473, 0.001850515018, 0.001534000403, 0.001932172136, 0.001928067481, 0.002317128886, 0.002314975858, 0.00180390305, 0.001900859149, 0.001958905783, 0.001774330601, 0.001533680315, 0.001528020481, 0.001623427959, 0.001935690547, 0.002355746166, 0.002137932195, 0.001859650234, 0.001979971247, 0.001837111486, 0.001555407534, 0.001604302022, 0.001777444928, 0.001741036688, 0.001807356548, 0.002012465449, 0.001956531438, 0.001817398768, 0.001711879501, 0.001614763026, 0.00157997703, 0.001621937332, 0.001595164909, 0.001853956393, 0.001857663227, 0.001717172202, 0.001858294563, 0.001651573889, 0.001444464685, 0.001524916872, 0.001519947769, 0.00150474773, 0.001646822523, 0.001968095145, 0.001772725629, 0.001573558801, 0.001742285217, 0.001398927468, 0.00134751789, 0.001524758113, 0.001405406504, 0.001348093794, 0.001599095327, 0.001447074662, 0.001360521947, 0.001412518689, 0.001338489337, 0.001231215769, 0.001236032472, 0.001274523524, 0.001226998722],
						[0.001942123379, 0.0020451527, 0.00193368615, 0.001910269587, 0.001833840701, 0.001826697256, 0.001802418959, 0.001773477959, 0.001622591534, 0.001595657068, 0.002074163928, 0.002067739944, 0.002040575238, 0.001718203143, 0.001891678493, 0.001729006647, 0.001600301608, 0.001734247778, 0.001594663869, 0.001906474095, 0.002266167678, 0.00189527577, 0.001658662239, 0.001777307171, 0.001682774242, 0.001834543289, 0.001687102396, 0.001821025008, 0.001755383966, 0.00181361643, 0.002247082258, 0.001959224997, 0.001950404018, 0.001865060743, 0.001683397023, 0.001639023288, 0.001865905876, 0.001554887457, 0.001822227792, 0.001756471334, 0.001976554273, 0.002088830654, 0.001602129373, 0.001811115973, 0.001872698697, 0.001647109971, 0.001455476451, 0.001494005094, 0.001643625112, 0.001934885145, 0.002025684896, 0.001988579044, 0.00158182962, 0.001877199216, 0.001635391534, 0.001528674383, 0.001605942258, 0.001822485902, 0.001707515433, 0.001639376038, 0.00185434367, 0.00171365323, 0.00158261204, 0.001552914235, 0.001592853435, 0.001570629667, 0.001611890165, 0.001618440704, 0.001761340961, 0.001724921252, 0.001572679132, 0.001781608032, 0.001642356721, 0.001531887011, 0.001519306473, 0.001537543115, 0.001612232676, 0.001668124132, 0.001963940144, 0.001556963915, 0.001467854587, 0.001398927468, 0.001698701726, 0.001327940391, 0.001616634301, 0.001494835972, 0.001437515905, 0.001637855863, 0.001541159024, 0.001283218879, 0.001300200301, 0.001356399067, 0.001310637426, 0.001323083041, 0.001367549979, 0.001368370588],
						[0.002101900141, 0.001968497738, 0.00201404631, 0.001976848994, 0.001935762069, 0.001950773694, 0.00179555948, 0.001877138149, 0.001671580843, 0.001697851274, 0.00205201753, 0.001998907565, 0.001982869769, 0.001682071098, 0.001711701281, 0.001644036295, 0.001561717078, 0.001685271276, 0.00156716494, 0.001893518972, 0.002117982045, 0.001799982984, 0.001736552634, 0.001766203399, 0.001597012896, 0.001731717203, 0.001621796456, 0.001726301118, 0.001768774387, 0.001834717014, 0.002097977091, 0.001948823312, 0.00188735876, 0.001717740157, 0.001668607442, 0.001532550038, 0.001676086598, 0.001469200418, 0.001709856006, 0.001761635758, 0.001845401642, 0.002000047114, 0.001579517001, 0.001725504202, 0.001695447636, 0.001577969709, 0.001490387539, 0.001475744626, 0.001512696411, 0.001863184141, 0.001889482664, 0.001773262991, 0.00155812396, 0.001696418522, 0.001627968138, 0.001389151699, 0.001543079407, 0.001672337281, 0.001644695146, 0.001659224599, 0.001648297176, 0.001608642565, 0.00160578189, 0.001560586611, 0.001466974956, 0.001497027997, 0.001470146558, 0.001438808775, 0.001761927511, 0.001551595006, 0.001451319933, 0.001627779459, 0.001562520092, 0.001363783562, 0.001411680124, 0.00147923634, 0.001417474375, 0.001633452578, 0.001932027363, 0.001482927742, 0.001348935717, 0.00134751789, 0.001327940391, 0.001499097178, 0.001464277672, 0.001374150364, 0.001351491444, 0.001519239066, 0.001488590535, 0.001073216745, 0.001183648946, 0.001176539898, 0.001128847482, 0.00120403918, 0.001211319905, 0.001150033874],
						[0.002560515946, 0.002535640801, 0.00253937852, 0.002485448666, 0.002350195792, 0.00231638808, 0.002302544639, 0.002263133226, 0.002067975693, 0.002067148721, 0.002574940276, 0.002534972683, 0.002484016339, 0.002112971385, 0.002146911221, 0.002084742402, 0.001859878712, 0.002121516629, 0.001981218007, 0.002302300689, 0.002668235806, 0.00224771988, 0.002087295719, 0.002226324109, 0.002031397406, 0.002172197335, 0.001933316665, 0.002169956628, 0.002103557194, 0.002292201244, 0.002585262202, 0.002238119388, 0.002251596033, 0.002234850987, 0.001948746689, 0.001956169542, 0.002156986478, 0.001768324436, 0.00215968941, 0.002084969638, 0.00234040907, 0.002459099058, 0.001933646242, 0.002020156092, 0.00206957304, 0.001900135167, 0.001746800454, 0.001733468366, 0.001903089234, 0.002231743966, 0.002335763509, 0.002244430763, 0.001802007614, 0.002103916017, 0.002071366733, 0.001813067579, 0.001825533626, 0.00202609524, 0.00200538994, 0.002125126885, 0.002036498587, 0.001983211094, 0.001834642637, 0.001794451241, 0.001730962735, 0.001760054456, 0.001863226423, 0.001889528393, 0.002223767524, 0.001880673805, 0.001753953548, 0.001970454534, 0.00190059034, 0.001697949351, 0.00173158236, 0.001773098228, 0.001773750525, 0.001904942811, 0.002236250352, 0.001793079798, 0.001656951629, 0.001524758113, 0.001616634301, 0.001464277672, 0.002158109585, 0.00167480383, 0.001628784164, 0.001809829016, 0.001810650856, 0.001328240001, 0.001424361999, 0.001421208669, 0.00141370992, 0.001458385528, 0.001509848185, 0.001468393513],
						[0.002299870077, 0.002346819807, 0.002258501074, 0.002206599569, 0.0021989561, 0.002196664573, 0.002102420911, 0.00219391231, 0.001944108977, 0.001983828019, 0.002281943433, 0.002236704862, 0.002223537119, 0.001897921416, 0.002034977418, 0.00191297072, 0.001712852285, 0.001922968469, 0.001904040506, 0.002221521965, 0.002302778545, 0.001895963612, 0.001880784241, 0.002049160193, 0.001888393315, 0.0020091459, 0.001898223352, 0.002032204841, 0.001956164834, 0.002216786869, 0.002213443299, 0.001904738534, 0.002023847758, 0.001997505152, 0.001817722133, 0.00177166698, 0.001938585802, 0.001737591391, 0.002086821575, 0.002156477303, 0.002079558915, 0.00217047674, 0.001659412079, 0.001792794597, 0.00182965049, 0.001711745192, 0.001655796305, 0.001701345065, 0.001771103512, 0.002133414729, 0.002037849418, 0.001871709159, 0.00162487014, 0.001831609016, 0.00171055106, 0.00165444263, 0.001679398022, 0.001939476904, 0.001967585467, 0.001869063588, 0.001800994314, 0.001738711241, 0.001570663998, 0.001666964874, 0.001611239802, 0.001638850111, 0.001699113779, 0.001777879072, 0.00202875252, 0.001614511409, 0.001512238654, 0.001773160458, 0.001684161069, 0.001514333937, 0.001670207893, 0.00170197913, 0.001785032259, 0.00181302171, 0.002092354005, 0.001564612875, 0.001443941498, 0.001405406504, 0.001494835972, 0.001374150364, 0.00167480383, 0.001943278547, 0.001557263926, 0.001799886207, 0.001633126472, 0.001208674032, 0.001284975258, 0.001287052074, 0.001306314015, 0.001310153149, 0.001438865939, 0.001397350548],
						[0.002083781247, 0.002105827515, 0.002197152045, 0.002148219717, 0.002055965995, 0.002116575774, 0.001963682724, 0.002079094009, 0.001812298457, 0.001884783965, 0.002136746694, 0.00217467103, 0.002091956136, 0.001773218193, 0.001895268085, 0.001819593425, 0.001679136737, 0.001818487049, 0.001840492434, 0.002142081755, 0.00226003792, 0.001801342803, 0.001813132445, 0.00191435273, 0.001786999549, 0.001939627161, 0.001877129076, 0.001974984732, 0.001946371307, 0.002177387041, 0.002178675095, 0.001925369836, 0.001959518691, 0.001918653967, 0.001742905187, 0.001680075186, 0.001892057444, 0.001611170221, 0.001997052557, 0.002065400695, 0.001874360924, 0.00203353407, 0.001606100416, 0.001736940424, 0.001730261696, 0.001660298011, 0.001665854964, 0.001608682258, 0.001686312902, 0.002105693351, 0.002005744056, 0.001750546507, 0.001567712893, 0.001815691681, 0.001684877891, 0.001623472287, 0.001668945698, 0.001907641035, 0.001912712957, 0.00195522756, 0.001637370346, 0.001653203216, 0.001592933608, 0.001556821933, 0.001517674255, 0.001596252795, 0.001664173048, 0.001792062365, 0.002015049614, 0.001502270361, 0.001449994849, 0.001636927628, 0.001583373631, 0.001467214757, 0.001543731544, 0.001674688296, 0.001751504465, 0.001790193869, 0.002225072321, 0.001518886766, 0.001344356872, 0.001348093794, 0.001437515905, 0.001351491444, 0.001628784164, 0.001557263926, 0.001869148927, 0.001821622391, 0.001754296788, 0.001100220997, 0.001245637172, 0.001224408789, 0.001235302959, 0.001315514696, 0.001423332905, 0.001436971018],
						[0.002701932892, 0.00248426969, 0.002615141548, 0.002382785761, 0.002359758516, 0.002469559865, 0.002249388048, 0.002384093743, 0.002112557084, 0.002081992134, 0.002578710274, 0.002359074018, 0.002527216656, 0.002041876751, 0.002260311334, 0.002199419717, 0.001911216259, 0.002148142537, 0.002101416757, 0.002591435225, 0.002571792485, 0.002218445191, 0.00215663356, 0.002332399788, 0.002043348134, 0.002327275906, 0.002004063264, 0.002326966181, 0.002220914779, 0.002613953357, 0.002510933318, 0.002235240477, 0.002326900403, 0.002221537306, 0.001849463061, 0.001911598533, 0.002211490898, 0.001913410436, 0.002346739526, 0.002442677105, 0.002396075867, 0.002300070477, 0.001862895769, 0.002063233178, 0.002064534697, 0.001848218608, 0.001962286642, 0.001935430794, 0.001977935665, 0.002335114657, 0.002458236182, 0.002199654503, 0.001873815258, 0.002105651673, 0.001950284629, 0.00179482723, 0.001919123145, 0.002103609887, 0.002272632476, 0.002256425015, 0.001989093522, 0.001914985238, 0.002057573856, 0.001761721129, 0.001656282469, 0.001828888462, 0.001824116362, 0.00216390251, 0.002359241942, 0.001976858326, 0.001779310054, 0.001953381381, 0.001914626162, 0.001568259676, 0.001747236502, 0.001831954788, 0.001892191381, 0.002071404173, 0.002734847721, 0.001774290456, 0.001572495244, 0.001599095327, 0.001637855863, 0.001519239066, 0.001809829016, 0.001799886207, 0.001821622391, 0.002979244709, 0.002074992858, 0.001353083366, 0.001470737906, 0.001491283712, 0.001517181308, 0.001570266791, 0.001544666878, 0.001551373396],
						[0.002757407576, 0.002719471677, 0.002802102426, 0.002707287202, 0.002504526251, 0.002616001643, 0.002506105964, 0.002466890693, 0.002276282331, 0.00232940124, 0.002542355418, 0.002538706295, 0.00245712161, 0.002086498421, 0.002183371878, 0.002243951131, 0.001958052651, 0.0021263257, 0.002207618824, 0.002621234068, 0.002718526671, 0.00220851091, 0.002089571752, 0.002324302701, 0.002036454966, 0.002234689403, 0.002049152968, 0.002348328682, 0.002292293272, 0.002567909878, 0.002546967054, 0.002217031467, 0.002387760302, 0.002176170363, 0.001900591909, 0.001958579231, 0.002181747593, 0.001930658194, 0.002276102439, 0.00243369268, 0.002140531648, 0.002407624321, 0.001994029146, 0.002008129068, 0.001966684595, 0.001831724253, 0.001977779456, 0.001909424209, 0.001958590849, 0.002455373533, 0.002416220994, 0.002098263899, 0.001786576186, 0.002104577206, 0.002041913208, 0.001723097603, 0.001938805242, 0.002081524935, 0.002188558276, 0.002297224771, 0.001834252413, 0.001870729067, 0.001761267425, 0.001674310159, 0.00164398893, 0.001783426135, 0.001925407783, 0.001970141906, 0.002349290251, 0.001714422288, 0.001611715347, 0.00179774555, 0.001781562354, 0.001741996881, 0.001804869635, 0.001804142466, 0.001904908508, 0.002108539649, 0.002739584572, 0.001573060605, 0.001432723256, 0.001447074662, 0.001541159024, 0.001488590535, 0.001810650856, 0.001633126472, 0.001754296788, 0.002074992858, 0.002635131973, 0.001167564371, 0.00130376408, 0.001329507691, 0.001380001238, 0.001321609962, 0.001336970271, 0.001462465311],
						[0.002009863366, 0.002091166113, 0.001993131787, 0.001804422913, 0.001829998069, 0.001804342367, 0.001728073759, 0.001729753224, 0.001537089641, 0.001503426019, 0.002360605143, 0.002178968325, 0.00192903563, 0.001765418608, 0.001703799545, 0.001705130056, 0.001415231144, 0.001657059318, 0.001389582413, 0.001826813062, 0.002277131577, 0.001925310451, 0.00168525719, 0.001714808166, 0.001631303246, 0.001603869682, 0.001490542702, 0.001651025562, 0.00153782559, 0.001524680113, 0.002496554737, 0.002028930962, 0.001960869109, 0.001858518206, 0.001683628808, 0.001715632971, 0.00161786628, 0.001297908725, 0.001649338559, 0.001676106267, 0.002277918698, 0.002036565539, 0.001724103763, 0.001679055023, 0.001747999998, 0.001566536254, 0.001304139577, 0.001308641757, 0.001336237643, 0.00165968462, 0.002355888284, 0.002135683637, 0.001529337994, 0.001764384174, 0.0014921953, 0.001368238167, 0.001453504753, 0.001497014391, 0.001539832924, 0.001523107793, 0.00224975273, 0.00181912707, 0.001508379831, 0.001338224906, 0.001406303153, 0.001414236235, 0.001432696674, 0.001485800171, 0.001502287893, 0.00193244399, 0.001544583527, 0.001620992162, 0.001410880098, 0.001342113321, 0.001302159289, 0.001318801308, 0.001385364582, 0.001503615307, 0.001727048959, 0.001739458821, 0.001472780098, 0.001360521947, 0.001283218879, 0.001073216745, 0.001328240001, 0.001208674032, 0.001100220997, 0.001353083366, 0.001167564371, 0.001732174103, 0.001355419126, 0.001356884386, 0.001206983927, 0.001050330198, 0.001120576959, 0.001069502127],
						[0.001980655134, 0.00204061465, 0.001941433241, 0.001881753917, 0.001840005844, 0.00183879967, 0.001709216314, 0.001775353392, 0.001600344549, 0.001564091681, 0.002220554732, 0.002083663409, 0.002014819832, 0.001761942727, 0.001711400505, 0.001712051313, 0.00147673357, 0.001654215713, 0.001484180669, 0.001795320974, 0.002209301972, 0.001872333325, 0.001848983201, 0.001772016561, 0.001626562968, 0.001739008661, 0.001573055978, 0.001691110869, 0.001667038649, 0.001779855608, 0.002295747956, 0.002013238493, 0.00197660144, 0.001821790762, 0.001661515879, 0.001600549369, 0.001671386618, 0.001312630625, 0.001713164575, 0.001761323688, 0.002093178431, 0.002083728244, 0.001573455973, 0.001674292826, 0.001856325133, 0.001659315375, 0.001386751283, 0.001388652916, 0.001430388361, 0.001753961274, 0.002182515245, 0.002002634763, 0.001611606875, 0.001754955104, 0.001613936507, 0.001434811152, 0.001421291285, 0.001545671162, 0.001644256099, 0.001675158753, 0.001895068104, 0.001775680814, 0.001597609687, 0.001469942223, 0.001505474531, 0.001439126627, 0.001447393032, 0.001609694131, 0.001690429149, 0.00173933189, 0.00155642601, 0.001691735074, 0.001494793329, 0.001277767174, 0.001408221757, 0.001365465744, 0.001454037065, 0.001489054779, 0.001749568877, 0.001649662703, 0.001464959819, 0.001412518689, 0.001300200301, 0.001183648946, 0.001424361999, 0.001284975258, 0.001245637172, 0.001470737906, 0.00130376408, 0.001355419126, 0.001475355593, 0.001258830108, 0.001160337151, 0.001174410199, 0.001188517155, 0.001154916643],
						[0.001876568152, 0.001945355621, 0.001881035982, 0.001755480548, 0.001776080995, 0.001747813931, 0.001684785524, 0.001711371032, 0.001539839272, 0.001499557562, 0.002113770986, 0.001987904628, 0.00197991251, 0.001704274891, 0.001713412, 0.001694201176, 0.001456591754, 0.001605786675, 0.001422802845, 0.001810262683, 0.002131248267, 0.001733466947, 0.001633653508, 0.001672395879, 0.001603387422, 0.001730536194, 0.001556844155, 0.001726021726, 0.001664785053, 0.001690825931, 0.002140649214, 0.001896815714, 0.001872576902, 0.001711593746, 0.00162246157, 0.001576340268, 0.001651221385, 0.001291381265, 0.001690114187, 0.001715023211, 0.002042077396, 0.002002175109, 0.001526185446, 0.001687232765, 0.001757823844, 0.001530303126, 0.001319070256, 0.001359905258, 0.001420022489, 0.001731624546, 0.002002751956, 0.001942345413, 0.001470394496, 0.001726516932, 0.001541106991, 0.001338967163, 0.001527167595, 0.001640709803, 0.001611616989, 0.001569474671, 0.001898527721, 0.001709574839, 0.001539238963, 0.001372897973, 0.001409622401, 0.001421393745, 0.001488967379, 0.001519896088, 0.001701696033, 0.001693475392, 0.001497240486, 0.001597226743, 0.001431770042, 0.001379396749, 0.001359346906, 0.001388263461, 0.001371315787, 0.001519445517, 0.001809098989, 0.001562606651, 0.00137748601, 0.001338489337, 0.001356399067, 0.001176539898, 0.001421208669, 0.001287052074, 0.001224408789, 0.001491283712, 0.001329507691, 0.001356884386, 0.001258830108, 0.001477072021, 0.001255149956, 0.001200739301, 0.001226329989, 0.001201041421],
						[0.001624975344, 0.00177718764, 0.001629664591, 0.001556326854, 0.001555161991, 0.001571649058, 0.0015206781, 0.001494355257, 0.001380351461, 0.001293417718, 0.001821547472, 0.001794990531, 0.001669107432, 0.00152144715, 0.001534514796, 0.001474786573, 0.001216366287, 0.001515116099, 0.001374981132, 0.001647121456, 0.001875955335, 0.001460972446, 0.00153148607, 0.001495796545, 0.00138335841, 0.00160289107, 0.001397123182, 0.001591942106, 0.001486070658, 0.001564463116, 0.00184150236, 0.001624917218, 0.00167473677, 0.001589475474, 0.001411743762, 0.00143020767, 0.001634855501, 0.001203881102, 0.00157166069, 0.001546563915, 0.001717972348, 0.001875702655, 0.001367473141, 0.001571458355, 0.001656088474, 0.001421089374, 0.001340951321, 0.001246940377, 0.001365864643, 0.0017056357, 0.001770018254, 0.001702335031, 0.001221941677, 0.00164620921, 0.001471749075, 0.001364436889, 0.001413216884, 0.001567264754, 0.001492172521, 0.001440842522, 0.001539200166, 0.001544124909, 0.001439716463, 0.001312547653, 0.001346738712, 0.001257311173, 0.001353134112, 0.001452871784, 0.001637735267, 0.0014430677, 0.001291858116, 0.0015912288, 0.001430295881, 0.001353480725, 0.001375020272, 0.00141299199, 0.001356353647, 0.001495595492, 0.001689155895, 0.001332975416, 0.001316196758, 0.001231215769, 0.001310637426, 0.001128847482, 0.00141370992, 0.001306314015, 0.001235302959, 0.001517181308, 0.001380001238, 0.001206983927, 0.001160337151, 0.001255149956, 0.001555627002, 0.001175399778, 0.001255432815, 0.001297431844],
						[0.001797528616, 0.001833740965, 0.001757523692, 0.001713310797, 0.001703873611, 0.001717569968, 0.00160284881, 0.001666919889, 0.001445855309, 0.001433147338, 0.001864275149, 0.001780425059, 0.001871868301, 0.001538600425, 0.00171764512, 0.001594540806, 0.001374620282, 0.001592331939, 0.001391746185, 0.001765786251, 0.00190288067, 0.001628559069, 0.001521509988, 0.001619156643, 0.001554231646, 0.00161558132, 0.001458963187, 0.001656064673, 0.001637206048, 0.001679989757, 0.001894217899, 0.001691962259, 0.001639575669, 0.001632520931, 0.001469458031, 0.001363441344, 0.001659412043, 0.001290540791, 0.001657481338, 0.001547746646, 0.001779658493, 0.001848791536, 0.001388784454, 0.001584612812, 0.001698686823, 0.001542052019, 0.00135730855, 0.001404625531, 0.00148274205, 0.001694330807, 0.001689529892, 0.001709835884, 0.001405332314, 0.00164119908, 0.001556784332, 0.001395907181, 0.001448010909, 0.001563318738, 0.001583504098, 0.001625692915, 0.001481003002, 0.001523394415, 0.001489612646, 0.001355054415, 0.001380300963, 0.001417501875, 0.001441163181, 0.001512179815, 0.001694576331, 0.001442626283, 0.001423933613, 0.001611998182, 0.001462799528, 0.001308016213, 0.001435210999, 0.001423903284, 0.001369128949, 0.001462829871, 0.001757563505, 0.00145065916, 0.001285898014, 0.001236032472, 0.001323083041, 0.00120403918, 0.001458385528, 0.001310153149, 0.001315514696, 0.001570266791, 0.001321609962, 0.001050330198, 0.001174410199, 0.001200739301, 0.001175399778, 0.00148087213, 0.001291575366, 0.001280782044],
						[0.001651044294, 0.001833126743, 0.001675342356, 0.001644961803, 0.001733502399, 0.001702643639, 0.001599105457, 0.001640984047, 0.001467154418, 0.001459815065, 0.001779315548, 0.001803667429, 0.001839191813, 0.001536824264, 0.001629135281, 0.001590322296, 0.001367862053, 0.001566729393, 0.00151998509, 0.001753677596, 0.001968385308, 0.001487323111, 0.001577008256, 0.001669654434, 0.00142310159, 0.001635188623, 0.001493556021, 0.001621731298, 0.001584095463, 0.001671031498, 0.001779130572, 0.001543014166, 0.001682717651, 0.001637769712, 0.001513598958, 0.001368196308, 0.001620250311, 0.001357172477, 0.001697401999, 0.0016355327, 0.001629050738, 0.001887412237, 0.001310831026, 0.001583883056, 0.001684434505, 0.001504801496, 0.001383488268, 0.001371462861, 0.001519772785, 0.001752820866, 0.001693807964, 0.001653003911, 0.001389306134, 0.001632717802, 0.001587092448, 0.001442945758, 0.001478431473, 0.001658344814, 0.001624889665, 0.001511484483, 0.001471806275, 0.001502743129, 0.001449216874, 0.001392112944, 0.001465920283, 0.001434609122, 0.001423381553, 0.001595801621, 0.001727829921, 0.001374532828, 0.001357856058, 0.001645564968, 0.001531357426, 0.001373687721, 0.001468558918, 0.001504110078, 0.001459390899, 0.001541539547, 0.001734222474, 0.001354462311, 0.001326462412, 0.001274523524, 0.001367549979, 0.001211319905, 0.001509848185, 0.001438865939, 0.001423332905, 0.001544666878, 0.001336970271, 0.001120576959, 0.001188517155, 0.001226329989, 0.001255432815, 0.001291575366, 0.001652398889, 0.001400097182],
						[0.001615803773, 0.001802919713, 0.001705850964, 0.001681966033, 0.001681532412, 0.001731275443, 0.001639530125, 0.001597138892, 0.001461220643, 0.001456114868, 0.00185400446, 0.001750839599, 0.001779332353, 0.00156334151, 0.00166261154, 0.001530771977, 0.001422734186, 0.001568749327, 0.001501016057, 0.00178581935, 0.00192768926, 0.001518940644, 0.001443522175, 0.001534071611, 0.001518029966, 0.001538851026, 0.001504941644, 0.001688645496, 0.001656696165, 0.001728721221, 0.00175253825, 0.001629675397, 0.001722636547, 0.001580531026, 0.001448566505, 0.001363541314, 0.001693878156, 0.001286333176, 0.001620050265, 0.001723242545, 0.001584058446, 0.001908320189, 0.001294276419, 0.001552290241, 0.00171185129, 0.00145597054, 0.001367615045, 0.001372454729, 0.001519431052, 0.001814234354, 0.001633407093, 0.001587142632, 0.001394237652, 0.00170226126, 0.001462726316, 0.001410140658, 0.001485237838, 0.001638997416, 0.001665134061, 0.00162087977, 0.001417212846, 0.001442843411, 0.001435943811, 0.001311553054, 0.001417331467, 0.001413057928, 0.001498832196, 0.00157436436, 0.001714672068, 0.001388198613, 0.001340842036, 0.00160054969, 0.001467325825, 0.001352817337, 0.001440521518, 0.001507816599, 0.001517580729, 0.001597668889, 0.001758184786, 0.00131376431, 0.001307264791, 0.001226998722, 0.001368370588, 0.001150033874, 0.001468393513, 0.001397350548, 0.001436971018, 0.001551373396, 0.001462465311, 0.001069502127, 0.001154916643, 0.001201041421, 0.001297431844, 0.001280782044, 0.001400097182, 0.001685882495]];

		var returns = [0.005873325, 0.007285383333, 0.0088992, 0.01417630833, 0.009235308333, 0.01318765, 0.01406648333, 0.01360380833, 0.01448148333, 0.01661586667, 0.004312425, 0.008158858333, 0.006681016667, 0.0005192, 0.006546558333, 0.005799741667, 0.009362716667, 0.005585858333, 0.01105469167, 0.01514394167, 0.0008603833333, 0.005390983333, 0.004285958333, 0.002965683333, 0.00882215, 0.01205829167, 0.007865941667, 0.009536483333, 0.01536361667, 0.01206999167, 0.001890775, 0.001219958333, 0.005129775, 0.00944175, 0.008676658333, 0.01137755, 0.008224466667, 0.01020840833, 0.01157376667, 0.007151091667, 0.005009416667, 0.005049083333, 0.006737966667, 0.01528215, 0.004935233333, 0.01058240833, 0.008808275, 0.0079628, 0.01132546667, 0.009103333333, 0.002207241667, 0.001660125, 0.00703745, 0.003528141667, 0.0058985, 0.00719915, 0.01374685833, 0.01278683333, 0.0099985, 0.008261508333, 0.0080636, 0.002264716667, 0.004111808333, 0.0021296, 0.0056057, 0.005727558333, 0.009533016667, 0.01281698333, 0.010599625, 0.006604125, 0.0015463, 0.001751933333, 0.003734616667, 0.008640916667, 0.007365616667, 0.0103935, 0.01234144167, 0.009268191667, 0.0082898, 0.004892658333, 0.004707458333, 0.005013725, 0.004047533333, 0.0034394, 0.006425441667, 0.007648291667, 0.005843025, 0.00818545, 0.0080509, 0.007949991667, 0.0058026, 0.003725691667, 0.001433233333, 0.002909541667, 0.006006575, 0.007056141667];
		
		// Compute the maximum desired volatility
		var maxVolatility = 0.10;
		
		// Compute the RSO-MVO weights
		var weights = PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, { nbRandomSubsets: 16384*1, 
		                                                                                                   subsetPortfolioOptimizationMethodParams: { 
																										       optimizationMethod: 'maximumTargetVolatility', 
																											   constraints: {
																												   maxVolatility: maxVolatility 
																											    }
																										     }
																										  });
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