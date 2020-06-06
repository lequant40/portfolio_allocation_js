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
		// Test default parameters
		var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]]); 
		var expectedWeights = [48, 24, 16, 12];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC - Values #1 ' + i);
		}
		
		// Test a number of cycles equal to 0, which must output
		// the initial internal weights, corresponding to an equally weighted portfolio
		var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]], {nbCycles: 0}); 
		var expectedWeights = [25, 25, 25, 25];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC - Values #1, no cycle ' + i);
		}
		
		// Test a number of cycles equal to 1
		//
		// Note: weights taken from running the algorithm once
		var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]], {nbCycles: 1}); 
		var expectedWeights = [0.38407763703774006, 0.24824686014490746, 0.1979048506362147, 0.16977065218113782];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(100 * weights[i]), Math.round(100 * expectedWeights[i]), 'ERC - Values #1, one cycle ' + i);
		}
		
		// Test a number of cycles equal to 10, which must converge to the
		// same solution as the default parameters
		var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]], {nbCycles: 10}); 
		var expectedWeights = [48, 24, 16, 12];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC - Values #1, ten cycles ' + i);
		}
		
		// Test a number of cycles equal to 100, which must converge to the
		// same solution as the default parameters
		var weights = PortfolioAllocation.equalRiskContributionWeights([[0.01, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.09, 0], [0, 0, 0, 0.16]], {nbCycles: 100}); 
		var expectedWeights = [48, 24, 16, 12];
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.round(100 * weights[i]), expectedWeights[i], 'ERC - Values #1, hundred cycles ' + i);
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
		var covMat = [[94.868,33.750,12.325,-1.178,8.778],
					[33.750,445.642,98.955,-7.901,84.954],
					[12.325,98.955,117.265,0.503,45.184],
					[-1.178,-7.901,0.503,5.460,1.057],
					[8.778,84.954,45.184,1.057,34.126]];
		
		var weights = PortfolioAllocation.equalRiskContributionWeights(covMat, {outputPortfolioVolatility: true}); 
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
	//
	// The ERC method converges in this case to a MV portfolio, and risk contributions
	// cannot be equal (there is no risk anymore).
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
		
		// Check the expected weights and the expected volatility
		var expectedWeights = [0.27653353547900505, 0.08992490010322292, 0.07660067915254527, 0.10281081303356884, 0.00627339817004538, 0.07548939567866735, 0.004543021673884871, 0.0858506918300374, 0.02204534860588049, 0.25992821627314233];
		var expectedVolatility = 0;
		for (var i = 0; i < expectedWeights.length; ++i) {
			assert.equal(Math.abs(weights[0][i] - expectedWeights[i]) <= 1e-4, true, 'ERC, semi-positive definite covariance matrix - Weight #' + i);
		}
		assert.equal(Math.abs(weights[1] - expectedVolatility) <= 1e-8, true, 'ERC, semi-positive definite covariance matrix - Volatility');
	}
	
	// Reference: Table 9 from the paper Constrained Risk Budgeting Portfolios https://arxiv.org/abs/1902.05710
	//
	// Test different coordinates sampling methods
	{
		/* The covariance matrix below has been obtained through Diag(STD)*CORR*Diag(STD)
		var stdDev = [0.05,0.05,0.07,0.1,0.15,0.15,0.15,0.18];
		var corrMat = [[100,  80,  60, -20, -10, -20, -20, -20],
						[ 80, 100,  40, -20, -20, -10, -20, -20],
						[ 60,  40, 100,  50,  30,  20,  20,  30],
						[-20, -20,  50, 100,  60,  60,  50,  60],
						[-10, -20,  30,  60, 100,  90,  70,  70],
						[-20, -10,  20,  60,  90, 100,  60,  70],
						[-20, -20,  20,  50,  70,  60, 100,  70],
						[-20, -20,  30,  60,  70,  70,  70, 100]]/100;
		*/
		var covMat = [[ 0.0025 ,  0.002  ,  0.0021 , -0.001  , -0.00075, -0.0015 , -0.0015 , -0.0018 ],
                      [ 0.002  ,  0.0025 ,  0.0014 , -0.001  , -0.0015 , -0.00075, -0.0015 , -0.0018 ],
                      [ 0.0021 ,  0.0014 ,  0.0049 ,  0.0035 ,  0.00315,  0.0021 , 0.0021 ,  0.00378],
                      [-0.001  , -0.001  ,  0.0035 ,  0.01   ,  0.009  ,  0.009  , 0.0075 ,  0.0108 ],
                      [-0.00075, -0.0015 ,  0.00315,  0.009  ,  0.0225 ,  0.02025, 0.01575,  0.0189 ],
                      [-0.0015 , -0.00075,  0.0021 ,  0.009  ,  0.02025,  0.0225 , 0.0135 ,  0.0189 ],
                      [-0.0015 , -0.0015 ,  0.0021 ,  0.0075 ,  0.01575,  0.0135 , 0.0225 ,  0.0189 ],
                      [-0.0018 , -0.0018 ,  0.00378,  0.0108 ,  0.0189 ,  0.0189 , 0.0189 ,  0.0324 ]];
					  
       var expectedWeights = [26.8306, 28.6769, 11.4095,  9.7985,  5.6135,  5.9029,  6.656,   5.1121];
	   
	   
	   // Check that the default coordinates sampler is 'cyclic', and that it converges
	   var defaultSamplingWeights = PortfolioAllocation.equalRiskContributionWeights(covMat);
	   var cyclicWeights = PortfolioAllocation.equalRiskContributionWeights(covMat, {coordinatesSampler: 'cyclic'});
	   assert.deepEqual(defaultSamplingWeights, cyclicWeights, 'ERC - Default coordinates sampler');
	   
	   for (var i = 0; i < expectedWeights.length; ++i) {
	      assert.equal(Math.abs(defaultSamplingWeights[i] - expectedWeights[i]/100) <= 1e-5, true, 'ERC - Default coordinates sampler weights ' + i);
	   }
	   
	   
	   // Check that the shuffled cyclic coordinates sampler converges
	   var shuffledCyclicWeights = PortfolioAllocation.equalRiskContributionWeights(covMat, {coordinatesSampler: 'shuffledCyclic'});
	   for (var i = 0; i < expectedWeights.length; ++i) {
	      assert.equal(Math.abs(defaultSamplingWeights[i] - expectedWeights[i]/100) <= 1e-5, true, 'ERC - Shuffled cyclic coordinates sampler weights ' + i);
	   }
	   
	   
	   // Check that the randomized coordinates sampler converges
	   //
	   // Note: There is no theoretical guarantees that the algorithm converges to a solution of the problem,
	   // but in practice with these data, it seems to converge.
	   var randomizedWeights = PortfolioAllocation.equalRiskContributionWeights(covMat, {coordinatesSampler: 'randomized'});
	   for (var i = 0; i < expectedWeights.length; ++i) {
	      assert.equal(Math.abs(randomizedWeights[i] - expectedWeights[i]/100) <= 1e-5, true, 'ERC - Randomized coordinates sampler weights ' + i);
	   }
	   
	   
	   // Check that the ACF coordinates samples converges
	   var acfWeights = PortfolioAllocation.equalRiskContributionWeights(covMat, {coordinatesSampler: 'acf'});
	   for (var i = 0; i < expectedWeights.length; ++i) {
	      assert.equal(Math.abs(acfWeights[i] - expectedWeights[i]/100) <= 1e-5, true, 'ERC - ACF coordinates sampler weights ' + i);
	   }
	}
	
	// Reference: Table 6 from the paper Constrained Risk Budgeting Portfolios https://arxiv.org/abs/1902.05710
	//
	// Test weights constraints on ERC portfolios
	{
		var covMat = [[0.0225, 0.003, 0.015, 0.0225, 0.0075],
		              [0.003, 0.04, 0.035, 0.024, 0.008],
					  [0.015, 0.035, 0.0625, 0.06, 0.00125],
					  [0.0225, 0.024, 0.06, 0.09, 0.003],
					  [0.0075, 0.008, 0.00125, 0.003, 0.01]];
		
		// Compute the min/max weights constraints
		var initWeights = [25.00, 25.00, 10.00, 10.00, 30.00];
		var minWeights = new Array(5);
		var maxWeights = new Array(5);
		for (var i = 0; i < minWeights.length; ++i) {
			minWeights[i] = initWeights[i]/100 - 0.05;
			maxWeights[i] = initWeights[i]/100 + 0.05;
		}
		
		// Compute the min/max weights-constrained ERC portfolio
		var weights = PortfolioAllocation.equalRiskContributionWeights(covMat, {constraints: { minWeights: minWeights, maxWeights: maxWeights } } );

		//
		var expectedWeights = [22.89, 20.00, 11.69, 10.42, 35.00];
	    for (var i = 0; i < expectedWeights.length; ++i) {
	       assert.equal(Math.abs(weights[i] - expectedWeights[i]/100) <= 1e-4, true, 'ERC - Bounds constraints weights ' + i);
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
	
	// Reference: Table 1 from the paper Constrained Risk Budgeting Portfolios https://arxiv.org/abs/1902.05710
	//
	// Test weights constraints on ERC and RB portfolios
	{
		var covMat = [[0.01,0.0075,0.01,0.015],
		              [0.0075,0.0225,0.015,0.0225],
					  [ 0.01,0.015,0.04,0.045],
					  [0.015,0.0225,0.045,0.09]];
					  
		// Compute ERC weights with no bounds constraints
		{
			var expectedWeights = [41.01, 27.34, 18.99, 12.66];
			
			var rb = [0.25, 0.25, 0.25, 0.25];
			var weights = PortfolioAllocation.riskBudgetingWeights(covMat, rb);

			for (var i = 0; i < expectedWeights.length; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]/100) <= 1e-3, true, 'RB - Test #3 ERC, ' + i);
			}
		}

		// Compute ERC weights with bounds constraints x_i <= 30%, i=1..4
		{
			var expectedWeights = [30, 30, 24.00, 16.00]; // the weights in the reference paper are [30, 30, 24.57, 15.43]
			
			var rb = [0.25, 0.25, 0.25, 0.25];
			var weights = PortfolioAllocation.riskBudgetingWeights(covMat, rb, {constraints: { maxWeights: [0.3, 0.3, 0.3, 0.3] } } );
			
			for (var i = 0; i < expectedWeights.length; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]/100) <= 1e-3, true, 'RB - Test #3 ERC with bounds constraints, ' + i);
			}
		}
		
		// Compute RB weights with no bounds constraints
		{
			var expectedWeights = [45.05, 30.04, 14.67, 10.24];
			
			var rb = [0.3, 0.3, 0.195, 0.205];
			var weights = PortfolioAllocation.riskBudgetingWeights(covMat, rb);
			
			for (var i = 0; i < expectedWeights.length; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]/100) <= 1e-3, true, 'RB - Test #3 RB, ' + i);
			}
		}
		
		// Compute RB weights with bounds constraints x_i <= 30%, i=1..4
		{
			var expectedWeights = [30, 30, 23.56, 16.43]; // the weights in the reference are [30, 30, 24.43, 15.57]
			
			var rb = [0.3, 0.3, 0.195, 0.205];
			var weights = PortfolioAllocation.riskBudgetingWeights(covMat, rb, {constraints: { maxWeights: [0.3, 0.3, 0.3, 0.3] } } );

			for (var i = 0; i < expectedWeights.length; ++i) {
				assert.equal(Math.abs(weights[i] - expectedWeights[i]/100) <= 1e-3, true, 'RB - Test #3 RB with bounds constraints, ' + i);
			}
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
			var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { minNbAssets: minAssets, maxNbAssets: maxAssets } });
			
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

	// Test with random data, exposure constraints
	{
		// Setup static parameters of the random test
		var nbTests = 10;
		var nbAssetsMin = 1;
		var nbAssetsMax = 50;

		// Aim of these tests is to check that for a portfolio of n assets with a constrained exposure,
		// the exposure is properly constrainted.
		for (var i = 0; i < nbTests; ++i) {
			// Generate a random number of assets
			var nbAssets = Math.floor(Math.random()*(nbAssetsMax - nbAssetsMin + 1) + nbAssetsMin);
			
			// Generate a random exposure constraint
			var maxExposure = Math.random();
			var minExposure = Math.min(Math.random(), maxExposure);

			// Generate a random portfolio for this number of assets/exposure
			var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { minExposure: minExposure, maxExposure: maxExposure } });

			// Check that the sum of the weights is belongs to the exposure interval, near to machine precision
			var weightsSum = 0;
			for (var k = 0; k < randomWeights.length; ++k) {
				weightsSum += randomWeights[k];
			}
			assert.equal(weightsSum <= maxExposure + 1e-14, true, "Random weights portfolio with exposure constraints, weights lower than max exposure - Test " + i);
			assert.equal(weightsSum >= minExposure - 1e-14, true, "Random weights portfolio with exposure constraints, weights greater than min exposure - Test " + i);
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
	
	// Test using static data in order to fix bug https://github.com/lequant40/portfolio_allocation_js/issues/6
	// The covariance matrix here has a determinant ~1e-22, and condition number ~5e6
	{
		// Problem data
		var covMat = [[0.04902254557, 0.04255604021, 0.04501327517, 0.04330779376, 0.03019691712, 0.02548665503, -0.01073017105, -0.0004006905689],
					[0.04682328959, 0.04455486658, 0.03786639547, 0.03747189194, 0.02769367774, 0.02256710184, -0.007460602423, -0.000360821725],
					[0.04501327517, 0.03441543433, 0.05846027012, 0.04801847343, 0.02887413717, 0.02797183226, -0.01440997349, -0.0003895354954],
					[0.04330779376, 0.03405688396, 0.04801847343, 0.04558680387, 0.03111517718, 0.02477230838, -0.01272882784, -0.0003624611328],
					[0.03019691712, 0.02516980916, 0.02887413717, 0.03111517718, 0.02614411029, 0.01475643353, -0.008794983792, -0.0002160623154],
					[0.02548665503, 0.02051044473, 0.02797183226, 0.02477230838, 0.01475643353, 0.01618991115, -0.006014483461, -0.0002507995642],
					[-0.01073017105, -0.006780679005, -0.01440997349, -0.01272882784, -0.008794983792, -0.006014483461, 0.005138124692, 0.00007878547574],
					[-0.0004006905689, -0.0003279381686, -0.0003895354954, -0.0003624611328, -0.0002160623154, -0.0002507995642, 0.00007878547574, 0.000007405024165]];
		var returns = [0.01807438883, 0.03238795043, 0.007555801824, 0.007427269126, 0.009034317809, 0.006707731718, 0.007769863126, 0.0007622417915];
		
		var minWeights = [0, 0, 0, 0, 0, 0, 0, 0];
		var maxWeights = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 1];
		
		var expectedCornerPortfolios = [[[0.0008860751090370526, 0, 0, 0.003345695950238685, 0.00007663523943489258, 0.012841518880073127, 0.01192011965466721, 0.97092995516609], 0.0009603955475062189, 0.0018269888034165023],
										[[0.0010092608916674394, 0, 0, 0.004008702066904519, 0.00033764763422769743, 0.012496184295674449, 0.01453557759610034, 0.9676126275149418], 0.0009853811817269642, 0.0018392482173261258], 
										[[0.0017655557593660093, 0, 0.008931031081884725, 0, 0.008263675566746306, 0.0063916733651148624, 0.04309168778365437, 0.9315563764426312], 0.0012618104609474637, 0.002243038039717897], 
										[[0.0012756123767303507, 0, 0.0032588378760465275, 0, 0.004072943929600115, 0.012021095630643686, 0.02160446204652544, 0.957767048138976], 0.0010630233852118463, 0.00190635190598824], 
										[[0, 0.004345917988265768, 0.006849418265594814, 0, 0.004007054238260629, 0.005038387387097269, 0.0284123197082556, 0.9513469024125258], 0.0012084215782098653, 0.0021408491446326377],
										[[0, 0.007621949314523461, 0.01148109322047808, 0, 0.00653096875456001, 0, 0.04935684701395931, 0.9250091416964791], 0.0014811875999548919, 0.0027647693325342387],
										[[0, 0.061938538867079784, 0.10166770263898256, 0, 0.11119983640121779, 0, 0.7251939220927198, 0], 0.009413515516400071, 0.033011575897091465],
										[[0, 0.25, 0.0405985104054761, 0, 0.009674272143084611, 0, 0.6997272174514392, 0], 0.013927927060719282, 0.053425941689894445],
										[[0, 0.25, 0.03158049051256315, 0, 0.023733980661979903, 0, 0.6946855288254571, 0], 0.013947635334077356, 0.05356625537194339],
										[[0.046211968331688, 0.25, 0.008523783192019019, 0, -2.7755575615628914e-17, 0, 0.6952642484762929, 0], 0.014398752755378262, 0.05712938291215708],
										[[0.05663429563220421, 0.25, 0, 0, -2.7755575615628914e-17, 0, 0.6933657043677959, 0], 0.014507974507069983, 0.05801956227886409],
										[[0.25, 0.25, 0, 0, -2.7755575615628914e-17, 0, 0.5, 0], 0.016500516378, 0.09086597054934813],
										[[0.25, 0.25, 0, 0, 0.25, 0, 0.25, 0], 0.01681663004875, 0.1309113982241424]];
				
		// Test that the algorithm is behaving properly
		var cornerPortfolios = PortfolioAllocation.meanVarianceCornerPortfolios(returns, covMat, { constraints: { minWeights: minWeights, maxWeights: maxWeights }});
		
		var cornerPortfoliosOK = true;
		for (var i = 0; i < expectedCornerPortfolios.length; ++i) {
			var cornerPortfolio = cornerPortfolios[i];
			var expectedCornerPortfolio = expectedCornerPortfolios[i];

			//Weights
			var weights = cornerPortfolio[0];
			var expectedWeights = expectedCornerPortfolio[0];
			for (var j = 0; j < expectedWeights.length; ++j) {
				if (Math.abs(weights[j] - expectedWeights[j]) > 1e-6) {
					cornerPortfoliosOK = false;
					break;
				}
			}
			
			// Return
			if (Math.abs(cornerPortfolio[1] - expectedCornerPortfolio[1]) > 1e-6) {
				cornerPortfoliosOK = false;
				break;
			}
			
			// Volatilities
			if (Math.abs(cornerPortfolio[2] - expectedCornerPortfolio[2]) > 1e-6) {
				cornerPortfoliosOK = false;
				break;
			}

		}
		assert.equal(cornerPortfoliosOK, true, 'Mean variance portfolio - Corner portfolios #7');
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
			// In case the maximum target volatility is lower than the minimum attainable volatility,
			// the portfolio cannot be fully invested and must correspond to the efficient portfolio with cash
			// and with a target volatility constraint 
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
			
			// Otherwise, the portfolio must correspond to the efficient portfolio with a target volatility constraint
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
	
	
	// Test using static data for bug fix https://github.com/lequant40/portfolio_allocation_js/issues/5
	// Test that min/max weights constraints are properly managed in case a risk free asset is added internally
	{
		// Problem data
		var covMat = [[0.0004246616877,-0.00005520852069,0.0003954256465,0.00000001152437917,0.0001590470407,0.0002580332644,0.0003335881244,0.0003605784739],
					  [-0.00005520852069,0.00002480135555,-0.00005059666822,-0.00000001082496581,-0.00001673167975,-0.00003073486553,-0.00003900214181,-0.00004548279667],
					  [0.0003954256465,-0.00005059666822,0.0003976604822,-0.000000008144301394,0.0001740534739,0.0002593777442,0.0003420673729,0.0003593307083],
					  [0.00000001152437917,-0.00000001082496581,-0.000000008144301394,0.00000001790999881,-0.00000008021468402,-0.00000006461917657,-0.00000007037516421,0.00000004991265269],
					  [0.0001590470407,-0.00001673167975,0.0001740534739,-0.00000008021468402,0.000110482962,0.0001154225601,0.000157800705,0.0001589926655],
					  [0.0002580332644,-0.00003073486553,0.0002593777442,-0.00000006461917657,0.0001154225601,0.0002185122484,0.0002506289478,0.0002558291246],
					  [0.0003335881244,-0.00003900214181,0.0003420673729,-0.00000007037516421,0.000157800705,0.0002506289478,0.0003278326867,0.0003411302813],
					  [0.0003605784739,-0.00004548279667,0.0003593307083,0.00000004991265269,0.0001589926655,0.0002558291246,0.0003411302813,0.0004078706675]];
		var returns = [0.02836270809, 0.01289860823, 0.003013461519, 0.0009824011995, -0.002791817978, -0.003228292398, -0.01390104713, -0.01415404063];
		
		var minWeights = [0, 0, 0, 0, 0, 0, 0, 0];
		var maxWeights = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 1];
		
		// Compute the desired maximum volatility
		var maxVolatility = 0.05/Math.sqrt(252);
		
		// Test that the algorithm is behaving properly
		var expectedWeights = [0.17660422287563507, 0.25, 0, 0.25, 0, 0, 0, 0];
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, { optimizationMethod: 'maximumTargetVolatility', 
		                                                                                     constraints: {maxVolatility: maxVolatility, minWeights: minWeights, maxWeights: maxWeights}});
		var weightsOK = true;
		for (var i = 0; i < expectedWeights.length; ++i) {
			if (Math.abs(weights[i] - expectedWeights[i]) > 1e-6) {
				weightsOK = false;
				break;
			}
		}
		assert.equal(weightsOK, true, 'Mean variance portfolio - maximum target volatility weights portfolio #5');
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
											 [[0.8793297501999342, 0, 0.12067024980006583], 0.06996423648680435, 0.12160182948150385], 
											 [[0.4207025920873124, 0.11248404158307494, 0.4668133663296126], 0.10225834167073272, 0.13608702325065056], 
											 [[0, 0.4187260623746124,0.5812739376253876], 0.13553706912274302, 0.17262859662810146], 
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
	// Tests the properties of the number of subsets to generate in case of random subsets generation with 2 assets, using random data 
	{
		// Limit case to test
		var nbAssets = 2;
		
		// Define the portfolio optimization method: an equal weigths portfolio optimization method, for simplicity
		function subsetEqualWeightsOptimization(subsetAssetsIdx, subsetPortfolioOptimizationMethodParams) {
			// Count the number of function calls
			subsetEqualWeightsOptimization.called++;
			
			return PortfolioAllocation.equalWeights(subsetAssetsIdx.length);
		}
		
		// Compute the RSO equally weighted portfolio
		// The default number of subsets in case of random subsets generation method to generate must be equal to 128
		var expectedFunctionCalls = 1;
		subsetEqualWeightsOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization);
		assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Default number of random subsets to generate with 2 assets');
		
		// Compute the RSO equally weighted portfolio
		// The number of subsets to generate in case of random subsets generation method must be equal to the specified one
		var nbSubsetsToGenerate = Math.floor(Math.random()*(1000-1+1) + 1); // max 1000 min 1
		var expectedFunctionCalls = 1;
		subsetEqualWeightsOptimization.called = 0;
		var weights = PortfolioAllocation.randomSubspaceOptimizationWeights(nbAssets, subsetEqualWeightsOptimization, {nbRandomSubsets: nbSubsetsToGenerate});
		assert.equal(subsetEqualWeightsOptimization.called, expectedFunctionCalls, 'Random subspace optimization method - Specified number of random subsets to generate with 2 assets');							
	}
	
	// Tests the properties of the number of subsets to generate in case of random subsets generation, using random data 
	{
		// Generate a random number of assets
		var nbAssets = Math.floor(Math.random()*(25-3+1) + 3); // max 25 min 3
			
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
		var expectedWeights = [0.21134794394954964, 0, 0.21134794394954964];
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
	
	// Test using static data infeasible cases
	{
		var covMat = [[1, 0.3, -0.2], 
		              [0.3, 1, 0.6], 
					  [-0.2, 0.6, 1]];
		var returns = [0.1, 0.2, 0.15];
		
		// Test a non-reachable target return
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, 
		                                                                                { subsetsGenerationMethod: 'deterministic', 
																			              subsetPortfolioOptimizationMethodParams: {
																							  optimizationMethod: 'targetReturn',
																							  constraints: {return: 0.3}
																						  }
																						}) },
			new Error('no feasible portfolio generated'),
			"Random subspace mean variance portfolio - No feasible portfolio generated #1");
			
		// Test a non-reachable target volatility
		assert.throws(function() { 
			PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, 
		                                                                                { subsetsGenerationMethod: 'deterministic', 
																			              subsetPortfolioOptimizationMethodParams: {
																							  optimizationMethod: 'targetVolatility',
																							  constraints: {volatility: 0.1}
																						  }
																						}) },
			new Error('no feasible portfolio generated'),
			"Random subspace mean variance portfolio - No feasible portfolio generated #1");
	}

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


QUnit.test('Best constantly rebalanced portfolio', function(assert) {    
	// Test with static data
	{
		// Define the assets price relatives
		var priceRelatives = [[1 - 0.05, 1 - 0.05, 1 - 0.05], // Asset 1 returns, all negatives
		                      [1, 1, 1], // Asset 2 returns, constants
							  [1 + 0.05, 1 + 0.05, 1 + 0.05]]; // Asset 3 returns, all positives
		
		// The BRCP should be the portfolio fully invested in the asset 3, as any other
		// portfolio would lead to a lower terminal wealth.
		var expectedWeights = [0, 0, 1];
		
		// Compute the weights of the BCRP
		var weights = PortfolioAllocation.bestConstantlyRebalancedWeights(priceRelatives);

		//
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
		assert.equal(weightsOK, true, 'Best constantly rebalanced portfolio, #1');
	}  
	

	// Test with random data, comparing the BCRP algorithm with an exhaustive search algorithm
	{
		// Define random assets price relatives for 5 assets and 24 periods
		var nbAssets = 5;
		var nbPeriods = 24;
		var priceRelativeMax = 1.20;
		var priceRelativeMin = 0.80;
		
		var priceRelatives = [];
		for (var i = 0; i < nbAssets; ++i) {
			// Initialization
			priceRelatives[i] = new Array(nbPeriods);
			
			//
			for (var j = 0; j < nbPeriods; ++j) {
				priceRelatives[i][j] = Math.random() * (priceRelativeMax - priceRelativeMin) + priceRelativeMin; // belongs to (priceRelativeMin, priceRelativeMax)
			}
		}
		
		// Compute the weights of the associated BCRP
		var weights = PortfolioAllocation.bestConstantlyRebalancedWeights(priceRelatives);
		
		// Compute the weights of the BRCP using a grid search algorithm:
		//
		// - Objective function to minimize: opposite of the portfolio cumulative wealth
		function portfolio_cumulative_wealth(arr) { 
			// Initialization
			var cumWealthPortfolio = 1;
			
			// Compute the evolution of the portfolio wealth
			for (var j = 0; j < nbPeriods; ++j) {
				// Compute the wealth relative for the current period
				var cumWealthPeriod = 0;
				for (var i = 0; i < nbAssets; ++i) {
					cumWealthPeriod += arr[i] * priceRelatives[i][j];
				}
				
				// Update the portfolio cumulative wealth
				cumWealthPortfolio *= cumWealthPeriod;
			}
			
			// Return it
			return -cumWealthPortfolio;
		}
		
		// - Numerical weights computation
		var numericalWeights = PortfolioAllocation.numericalOptimizationWeights(nbAssets, portfolio_cumulative_wealth, {optimizationMethod: 'grid-search', optimizationMethodParams: {k: 10}})
		
		// Compare the resulting cumulative portfolio wealths, which must be numerically equal
		assert.equal(Math.abs(portfolio_cumulative_wealth(weights) - portfolio_cumulative_wealth(numericalWeights[0])) <= 1e-2, true, 'Best constantly rebalanced portfolio, #2');
	}

});


QUnit.test('Minimum tracking error portfolio', function(assert) {    
	// Test with static data: 1 EEM ETF, 1 WLD ETF, 1 ACWI ETF as benchmark
	{
		// Define the assets returns
		var assetsReturns = [[0.078129363,0.03663478,0.028787879,0.031491869,-0.019533386,-0.042308288,-0.062133001,-0.103910104,-0.027053038,0.082179315,0.004579357,-0.049970429], // Asset 1 returns in 2015 (EEM)
		                     [0.05289784,0.065017266,0.028012939,-0.019089773,0.025592817,-0.03893898,0.026564726,-0.07927768,-0.033253707,0.090633334,0.04075,-0.044859476]]; // Asset 2 returns in 2015 (WLD)
		
		// Define the benchmark returns in 2015 (ACWI)
		var benchmarkReturns = [0.055518502,0.062009871,0.028003081,-0.013945023,0.020509392,-0.03935567,0.017029792,-0.081809115,-0.032750097,0.089738414,0.037148934,-0.045498531];

		// Prior to any numerical computation, since ACWI index in 2015 ~= 90% World index + 10% Emerging index,
		// the minimum tracking error portfolio should have weights close to [0.1, 0.9], which has been
		// validated by the numerical computation below.
		var expectedWeights = [0.10614062381555761, 0.8938593761844424];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns);

		//
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, #1');
	}  
	
	// Test with static data for redundancy: 2 EEM ETFs, 2 WLD ETFs, 1 ACWI ETF as benchmark
	{
		// Define the assets returns
		var assetsReturns = [[0.078129363,0.03663478,0.028787879,0.031491869,-0.019533386,-0.042308288,-0.062133001,-0.103910104,-0.027053038,0.082179315,0.004579357,-0.049970429], // Asset 1 returns in 2015 (EEM)
		                     [0.05289784,0.065017266,0.028012939,-0.019089773,0.025592817,-0.03893898,0.026564726,-0.07927768,-0.033253707,0.090633334,0.04075,-0.044859476], // Asset 2 returns in 2015 (WLD)
							 [0.055323183,0.065000264,0.027985908,-0.019018198,0.025635979,-0.038860104,0.026504942,-0.079552638,-0.033599239,0.090526431,0.041004562,-0.044589974], // Asset 3 returns in 2015 (WLD #2)
							 [0.078085397,0.036893028,0.02912706,0.031861707,-0.019134434,-0.041778299,-0.0617737,-0.103650587,-0.026909091,0.082105488,0.005130735,-0.049278492]]; // Asset 4 returns in 2015 (EEM #2)
		
		// Define the benchmark returns in 2015 (ACWI)
		var benchmarkReturns = [0.055518502,0.062009871,0.028003081,-0.013945023,0.020509392,-0.03935567,0.017029792,-0.081809115,-0.032750097,0.089738414,0.037148934,-0.045498531];

		
		// The default numerical computation allocates evenly between the 2 WLD ETFs and the 2 EEM ETFs.
		var expectedWeightsDefault = [0.05152357319477277, 0.4486638440288671, 0.44666368775168186, 0.053148895024678164];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns);

		//
		var weightsOK = true;
		if (expectedWeightsDefault.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsDefault.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsDefault[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, no cardinality constraints, duplicate assets');
		
		
		// The weights constrained numerical computation allocates between the 2 WLD ETFs and the 2 EEM ETFs according to the weights constraints.
		var expectedWeightsMinMaxWeights = [0.25, 0.5, 0.15, 0.1];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minWeights: [0.25, 0, 0, 0.1], maxWeights: [1, 0.5, 1, 1] }});

		//
		var weightsOK = true;
		if (expectedWeightsMinMaxWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsMinMaxWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsMinMaxWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, weights constraints, duplicate assets #1');
		
		
		// The 2-2 cardinality constrained numerical computation allocates 1 WLD ETF and 1 EEM ETF.
		var expectedWeightsTwoTwo = [0.10614062381555761, 0.8938593761844424, 0, 0];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minNbAssets: 2, maxNbAssets: 2 }});

		//
		var weightsOK = true;
		if (expectedWeightsTwoTwo.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsTwoTwo.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsTwoTwo[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality constraints, duplicate assets #1');

		
		// The 2-2 cardinality constrained numerical computation allocates to the same EEM ETF as above, 
		// but to the other WLD ETF in case maximum weights constraints are additionally imposed on the first WLD ETF.
		var expectedWeightsMaxTwoTwo = [0.10292768040365563, 0, 0.8970723195963444, 0];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { maxWeights: [1,0.25,1,1], minNbAssets: 2, maxNbAssets: 2 }});

		//
		var weightsOK = true;
		if (expectedWeightsMaxTwoTwo.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsMaxTwoTwo.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsMaxTwoTwo[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality and max weights constraints, duplicate assets #1');
		
		
		// The 2-3 cardinality constrained numerical computation allocates to the same EEM ETF as above, 
		// but dispatch to the two WLD ETFs in case maximum weights constraints are additionally imposed on the first WLD ETF.
		var expectedWeightsMaxTwoThree = [0.09775016142607804, 0.2441086935007033, 0.6581411450732186, 0];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { maxWeights: [1,0.25,1,1], minNbAssets: 2, maxNbAssets: 3 }});

		//
		var weightsOK = true;
		if (expectedWeightsMaxTwoThree.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsMaxTwoThree.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsMaxTwoThree[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, relaxed cardinality and max weights constraints, duplicate assets #1');
				
			
		// The 2-2 cardinality constrained with min weights numerical computation allocates 
		// to a different WLD ETF than above in case the min weight constraint would force
		// the ideal WLD ETF to have too much weight (here, 95% v.s. its ideal allocation of ~90%)
		var expectedWeightsMinTwoTwo = [0, 0, 0.8, 0.2];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minWeights: [0.25, 0.95, 0, 0.2], minNbAssets: 2, maxNbAssets: 2 }});

		//
		var weightsOK = true;
		if (expectedWeightsMinTwoTwo.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsMinTwoTwo.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsMinTwoTwo[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality and min weights constraints, duplicate assets #1');
		
		
		// Infeasible 2-2 cardinality constrained computation
		assert.throws(function() { 
			PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minWeights: [0.6, 0.6, 0.6, 0.6], minNbAssets: 2, maxNbAssets: 2 }}); },
			new Error('infeasible problem detected'),
			"Minimum tracking error portfolio, infeasible cardinality constraints, duplicate assets #1");
		
		
		// The 2-2 cardinality constrained with soft inequality constraints numerical computation 
		// must allocate the same as the 2-2, with a small lambda.
		var expectedWeights = expectedWeightsTwoTwo;
		
		// Compute the weights of the minimum tracking error portfolio
		//
		// Soft inequality constraints are w_1 + w_2 <= 0.5 and w_3 + w_4 <= 0.5 
		// (i.e. less than half portfolio allocation in the 2 groups of ETFs, 
		// for instance supposed to be issued by different emitters)
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, 
		                                                              {constraints: { minNbAssets: 2, maxNbAssets: 2 }, 
																	  softConstraints: { Ai:[[1,1,0,0],[0,0,1,1]], bi:[0.5, 0.5], lambdai:1e-8 }} );

		//
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality and inequality constraints, duplicate assets #2');
		
		
		// The 2-2 cardinality constrained with soft inequality constraints numerical computation 
		// must start to balance the ETFs emitters with a bigger lambda.
		var expectedWeights = [0, 0.8909507489391898, 0, 0.10904925106081023];
		
		// Compute the weights of the minimum tracking error portfolio
		//
		// Soft inequality constraints are w_1 + w_2 <= 0.5 and w_3 + w_4 <= 0.5 
		// (i.e. less than half portfolio allocation in the 2 groups of ETFs, 
		// for instance supposed to be issued by different emitters)
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, 
		                                                              {constraints: { minNbAssets: 2, maxNbAssets: 2},
																	   softConstraints: { Ai:[[1,1,0,0],[0,0,1,1]], bi:[0.5, 0.5], lambdai:1e-4 }});

		//
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality and inequality constraints, duplicate assets #3');
		
		
		// The 2-2 cardinality constrained with soft inequality constraints numerical computation 
		// must drop the EEM ETFs completely in favor of balancing the ETFs emitters with an even
		// bigger lambda.
		var expectedWeights = [0, 0.49999999999999983, 0.5000000000000002, 0];
		
		// Compute the weights of the minimum tracking error portfolio
		//
		// Soft inequality constraints are w_1 + w_2 <= 0.5 and w_3 + w_4 <= 0.5 
		// (i.e. less than half portfolio allocation in the 2 groups of ETFs, 
		// for instance supposed to be issued by different emitters)
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, 
		                                                              {constraints: { minNbAssets: 2, maxNbAssets: 2}, 
																	   softConstraints: { Ai:[[1,1,0,0],[0,0,1,1]], bi:[0.5, 0.5], lambdai:1 }});

		//
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality and inequality constraints, duplicate assets #4');
		
		
		// The 1-4 cardinality constrained numerical computation allocates the same as the 2-2.
		var expectedWeightsOneFour = expectedWeightsTwoTwo;
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minNbAssets: 1, maxNbAssets: 4 }});

		//
		var weightsOK = true;
		if (expectedWeightsOneFour.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsOneFour.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsOneFour[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality constraints, duplicate assets #5');
		
		
		// The 4-4 cardinality constrained numerical computation must allocate the same
		// as the non cardinality constrained one, since the non-constrained one is made of 4 assets.
		var expectedWeightsFourFour = expectedWeightsDefault;
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, benchmarkReturns, {constraints: { minNbAssets: 4, maxNbAssets: 4 }});

		//
		var weightsOK = true;
		if (expectedWeightsFourFour.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeightsFourFour.length; ++i) {
				if (Math.abs(weights[i] - expectedWeightsFourFour[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, cardinality constraints, duplicate assets #6');
	} 
	
	
	// INDTRACK1 (Hang Seng Index and Constituents)
	// Reference: OR-library http://people.brunel.ac.uk/~mastjjb/jeb/orlib/indtrackinfo.html
	{
		// Define the index prices
		var indexPrices = [8749.31759356, 8713.53262793, 9001.60515134, 8903.30299873, 8888.69435453, 8727.37800153, 8735.20596277, 9022.01820217, 9461.30705604, 9734.73543429, 9954.80587273, 10125.76073581, 10401.21266867, 10792.5574788, 10713.83410451, 10329.12797352, 10273.26721616, 10851.04530705, 10760.07410226, 10617.20049693, 10877.60002361, 10246.92550535, 10439.28744609, 10385.68099956, 10335.37614212, 9779.43114037, 9991.03460033, 10131.10362998, 9828.95497674, 10098.81550871, 10093.20635736, 9932.17401202, 9934.42832285, 10624.24743709, 11116.30846497, 10988.77127348, 11209.4452282, 11300.94894735, 10434.05105479, 10626.62600127, 9351.12983299, 9324.64611835, 9217.22021955, 9659.82841299, 9784.75628407, 9815.23385648, 10422.72624919, 10498.32553918, 10209.01048224, 10201.57303155, 10398.49684538, 10952.20528678, 11275.08649948, 11542.7992235, 10952.75555162, 10807.09512109, 11256.02248505, 11368.18776174, 11157.36532286, 11887.76203205, 11983.2418585, 12124.48241977, 12067.37912888, 12434.72529143, 12727.02242884, 13131.64459726, 12704.88758155, 12898.81156444, 12781.67615363, 12450.32796246, 12812.17147652, 12550.67142015, 12292.59720621, 11982.15907928, 12405.73875925, 13129.62104265, 13121.59782615, 13393.37540986, 13176.14504807, 13334.3905682, 13478.73746329, 13167.5005648, 13357.32418705, 13625.65817783, 14210.25245265, 15556.48427968, 15477.21064054, 16559.58159748, 15969.98168759, 17220.18342331, 16443.93722691, 16462.53972888, 16497.9164335, 18155.385159, 18760.23273035, 19596.17378757, 21102.46168358, 19528.01194828, 19124.72106578, 20340.93063426, 20196.19322863, 21580.26907684, 20420.20427339, 19216.45554117, 17928.42753473, 17605.26231437, 17583.0209642, 16210.28767206, 16391.16505289, 16028.52276727, 16504.78586886, 16926.96326099, 16253.84734748, 15915.20370944, 15301.39214621, 16214.56553749, 17096.60458863, 16809.93435306, 16391.27155577, 16172.74540901, 16177.71554312, 15764.20038474, 15326.42032158, 14968.90792439, 16183.14718969, 16246.99566259, 16832.44196041, 17044.38267946, 16800.04733626, 16693.33145662, 16683.81719988, 17575.74326782, 18007.89642914, 17694.60047509, 17098.09562887, 16900.65705114, 16481.1067299, 16953.35822324, 16576.44455252, 16649.0085106, 16916.91648989, 16628.38245402, 16734.15755828, 15369.83799317, 14593.68054917, 13825.97233315, 14495.73340614, 14753.98512486, 14539.48833683, 13638.13676456, 12873.2508747, 12918.97610856, 12952.73751959, 13275.4412275, 14223.1393004, 14276.72799645, 14589.02992367, 14529.03330471, 14110.54801221, 15149.51904804, 15060.67790073, 15243.61433712, 15035.01070812, 15375.44714453, 15345.98134942, 14841.22872906, 14790.00084672, 16362.08976832, 15999.07472264, 16559.74135179, 16968.99639522, 16450.46940317, 16532.7073723, 16373.68083109, 16342.04947755, 17097.70511833, 17267.59495267, 16703.60898395, 16777.20246981, 16619.47171357, 15987.73216657, 15790.50659459, 16118.12718518, 16324.17474521, 16673.39766872, 17390.85427873, 16937.27628928, 17122.71554322, 17526.64544295, 17544.18291619, 17564.5249651, 17183.7949414, 17494.51707599, 16706.48456154, 16486.46737455, 16842.98574492, 17506.49864931, 17508.38020008, 17499.5759625, 17630.13073542, 17880.74974818, 18691.07686423, 18709.02259848, 19106.77533153, 19724.1014896, 20358.73436468, 20076.2887431, 20581.66263022, 20218.52333119, 19871.55471851, 19912.11456299, 18740.21019006, 19572.97391154, 19449.55483117, 19773.82058123, 19258.91468689, 19203.3379372, 19051.21633232, 19053.79015177, 18811.47836317, 19200.42685865, 19559.53679895, 19995.43531133, 19874.41254563, 19285.87766447, 19268.65969986, 19562.62538229, 19839.94111544, 19175.27442993, 19250.92697135, 19002.89952852, 19458.02180965, 19710.18511408, 19837.89981036, 20279.28322075, 19807.794998, 19570.95035694, 20180.5905576, 20576.99425425, 20873.4805047, 21132.8505036, 21688.2452405, 22205.93796003, 21989.96788225, 22240.05438064, 22633.91975879, 22879.24912881, 23282.93052184, 23774.86729638, 23257.97334839, 22692.33658512, 23308.88172212, 23792.99053542, 23471.08559906, 23415.54435033, 24595.7736979, 23749.3421076, 23646.81534099, 24248.04181463, 23276.66460276, 23865.25273536, 23783.36977581, 23674.43508629, 22607.95080804, 22169.10571615, 22249.01837253, 21663.73182902, 22217.56452377, 22261.19520111, 22446.82971031, 23220.64409109, 24727.83726153, 24961.38031351, 25439.418463, 26195.8196239, 26013.5577057, 25050.45221705, 26899.71486797, 26975.03015029, 26311.48174495, 27025.61901539, 27638.20579557, 27793.91299721, 29073.90003667, 29550.18088874, 28572.73301307, 27388.54530868, 25090.74580434, 25850.99881918, 25685.75961033, 25532.51972527];
		var indexReturns = new Array(indexPrices.length - 1);
		for (var i = 1; i < indexPrices.length; ++i) {
			indexReturns[i-1] = (indexPrices[i] - indexPrices[i-1])/indexPrices[i-1];
		}
		
		// Define the 31 assets prices
		var assetsPrices = [[9.33675195, 9.86926631, 9.95801871, 9.15924716, 9.24799955, 8.98174237, 9.07049476, 8.73323566, 9.24799955, 9.33675195, 10.93429505, 10.93429505, 10.93429505, 11.37805703, 12.07032571, 11.44905894, 11.37805703, 11.53781134, 11.11179984, 10.93429505, 11.02304745, 10.31302829, 10.0467711, 10.0467711, 10.31302829, 9.51425673, 9.69176152, 9.78051392, 9.60300913, 10.75679026, 10.93429505, 10.49053308, 10.57928547, 11.28930463, 11.37805703, 10.93429505, 11.37805703, 12.1590781, 10.93429505, 11.28930463, 9.33675195, 8.82198805, 8.89298997, 9.51425673, 9.69176152, 10.0467711, 10.66803787, 11.02304745, 10.1355235, 10.40178068, 10.40178068, 10.57928547, 11.20055224, 11.44905894, 11.02304745, 11.02304745, 11.71531613, 11.89282092, 11.37805703, 11.89282092, 11.89282092, 12.51408768, 13.13535445, 12.78034487, 12.60284008, 13.04660205, 12.60284008, 13.13535445, 12.95784966, 12.95784966, 12.95784966, 12.78034487, 12.1590781, 11.44905894, 11.89282092, 12.42533529, 11.89282092, 11.62656373, 12.69159247, 12.69159247, 12.78034487, 12.86909726, 13.75662121, 13.66786882, 14.00512792, 15.24766145, 16.1351854, 18.33624479, 16.56119689, 18.15874, 18.15874, 17.62622563, 17.53747324, 19.38352305, 19.04626395, 21.67333484, 26.27070889, 21.15857095, 20.62605658, 22.91586837, 23.62588753, 25.38318495, 24.495661, 22.73836358, 22.02834442, 21.67333484, 21.15857095, 20.44855179, 19.91603742, 20.09354221, 20.44855179, 21.85083963, 20.98106616, 20.98106616, 20.09354221, 21.85083963, 21.85083963, 21.67333484, 21.67333484, 21.33607574, 20.98106616, 20.62605658, 20.89231376, 20.44855179, 20.71480897, 20.62605658, 20.80356137, 20.80356137, 20.80356137, 21.33607574, 20.1822946, 20.98106616, 21.67333484, 21.06981855, 20.44855179, 20.1822946, 19.56102784, 19.38352305, 18.60250197, 18.42499718, 19.11726586, 18.15874, 17.89248281, 15.86892821, 14.53764229, 13.57911642, 15.69142342, 15.42516624, 15.78017582, 14.62639468, 13.31285924, 13.13535445, 12.69159247, 13.31285924, 14.44888989, 14.53764229, 15.15890905, 15.24766145, 15.42516624, 15.60267103, 15.86892821, 16.046433, 16.64994929, 16.64994929, 16.56119689, 16.56119689, 16.1351854, 17.71497802, 17.27121605, 17.80373042, 18.24749239, 16.91620647, 17.18246366, 17.00495887, 16.82745408, 17.53747324, 17.44872084, 16.91620647, 17.18246366, 17.35996845, 16.49019498, 16.1351854, 16.1351854, 16.64994929, 17.09371126, 17.09371126, 16.91620647, 17.00495887, 17.18246366, 17.00495887, 17.00495887, 16.56119689, 17.35996845, 16.31269019, 16.046433, 15.95768061, 17.44872084, 18.15874, 18.95751155, 19.47227544, 19.73853263, 21.06981855, 20.44855179, 21.85083963, 22.02834442, 21.51358053, 21.06981855, 22.383354, 22.02834442, 21.67333484, 22.02834442, 21.06981855, 22.47210639, 22.2946016, 22.11709681, 22.20584921, 23.27087795, 23.27087795, 23.00462076, 22.82711597, 22.82711597, 23.44838274, 24.0696495, 23.62588753, 23.62588753, 23.35963034, 23.27087795, 23.53713513, 23.27087795, 23.27087795, 23.35963034, 23.27087795, 23.35963034, 23.18212555, 23.98089711, 23.53713513, 23.35963034, 24.0696495, 23.98089711, 24.31815621, 24.22940381, 24.67316579, 25.11692776, 24.4069086, 25.47193734, 24.58441339, 24.76191818, 25.47193734, 27.33573763, 26.53696608, 24.93942297, 25.73819453, 26.0932041, 25.02817537, 25.38318495, 26.27070889, 25.29443255, 25.29443255, 25.29443255, 24.22940381, 25.20568016, 24.22940381, 24.495661, 22.36560352, 21.12306999, 21.74433675, 20.9455652, 21.47807957, 20.50180323, 19.88053646, 20.50180323, 21.03431759, 21.30057478, 21.21182238, 21.30057478, 20.8568128, 21.38932717, 22.18809873, 21.92184154, 20.76806041, 20.9455652, 20.50180323, 19.34802209, 18.46049814, 19.1705173, 18.63800293, 17.57297419, 16.86295503, 17.04045982, 17.04045982, 17.12921222],
							[14.37788798, 14.18263271, 14.59089372, 15.74467486, 16.29493971, 15.6204215, 15.74467486, 16.98720839, 18.97526203, 18.76225628, 18.90426012, 18.97526203, 20.12904317, 22.08159585, 22.25910064, 21.92184154, 22.25910064, 24.31815621, 26.02220219, 25.68494309, 26.35946129, 24.14065142, 24.47791052, 23.96314663, 24.31815621, 23.28862843, 23.11112364, 23.80339232, 23.28862843, 23.44838274, 23.62588753, 23.44838274, 23.96314663, 25.1701792, 25.68494309, 24.65541531, 24.8329201, 24.8329201, 22.77386454, 24.14065142, 20.41305083, 20.71480897, 21.92184154, 22.59635975, 23.28862843, 23.96314663, 24.47791052, 23.96314663, 24.47791052, 24.65541531, 25.47193734, 25.27668207, 26.0932041, 27.31798715, 26.0932041, 25.47193734, 26.0932041, 26.0932041, 25.88019836, 27.93925392, 28.34751494, 29.37704272, 29.18178745, 29.18178745, 29.37704272, 31.43609828, 30.81483151, 31.22309253, 30.81483151, 29.78530373, 30.4065705, 29.78530373, 29.37704272, 28.75577595, 27.74399865, 28.15225967, 27.74399865, 28.15225967, 26.90972614, 26.71447087, 26.71447087, 26.28845937, 26.90972614, 26.71447087, 28.34751494, 29.78530373, 29.37704272, 31.22309253, 30.19356475, 32.25262031, 32.46562606, 34.09867013, 36.15772569, 39.03330328, 40.8793531, 41.50061986, 43.55967542, 45.60098051, 46.84351404, 45.19271949, 43.96793644, 46.026992, 44.78445847, 41.85562944, 35.18144934, 32.87388708, 34.15192156, 32.092866, 34.93294264, 33.63715767, 33.63715767, 34.66668545, 33.90341486, 33.38865097, 32.87388708, 34.15192156, 37.7552688, 38.51853939, 37.24050491, 36.21097713, 36.97424772, 36.21097713, 34.205173, 31.73785642, 34.31167587, 32.55437845, 32.98038995, 33.17564522, 32.7673842, 32.25262031, 31.54260115, 33.28214809, 36.15772569, 35.53645892, 34.205173, 33.90341486, 33.17564522, 34.205173, 33.90341486, 33.79691198, 34.82643976, 34.82643976, 35.12819791, 33.49515384, 31.8443593, 28.45401781, 29.99830948, 32.05736504, 31.73785642, 30.56632481, 28.86227883, 28.70252452, 24.14065142, 23.87439423, 26.19970698, 26.19970698, 25.82694692, 25.63169165, 24.70866674, 26.07545363, 27.05172997, 26.92747662, 26.12870506, 26.67896991, 26.57246704, 25.95120027, 25.40093542, 27.92150344, 27.79725009, 28.41851685, 28.66702356, 28.59602164, 28.84452835, 28.72027499, 28.72027499, 31.73785642, 31.98636313, 31.56035163, 32.28812127, 32.05736504, 30.93908487, 29.94505804, 30.56632481, 30.56632481, 31.06333822, 31.86210977, 30.19356475, 30.88583343, 32.11061648, 32.78513468, 33.65490815, 33.33539953, 33.90341486, 33.03364139, 33.28214809, 33.45965288, 34.6311845, 34.15192156, 33.83241294, 33.95666629, 34.205173, 35.44770653, 35.87371802, 36.97424772, 38.83804801, 40.68409783, 38.83804801, 40.93260453, 40.31133777, 40.68409783, 40.4533416, 38.0925279, 39.63681957, 39.19305759, 39.63681957, 38.74929562, 39.56581765, 38.02152598, 38.16352981, 37.4180097, 38.53628987, 39.49481574, 40.38233969, 42.08638567, 40.82610166, 40.82610166, 41.85562944, 42.08638567, 40.38233969, 40.96810549, 40.96810549, 41.57162178, 43.18691537, 41.19886172, 41.85562944, 41.48286938, 40.75509974, 41.78462753, 42.15738758, 41.71362561, 42.53014764, 43.41767159, 44.81995943, 44.6779556, 44.51820129, 46.4530035, 47.18077314, 49.39958301, 48.95582103, 50.00309929, 49.11557535, 49.55933732, 50.73086893, 50.28710696, 51.33438522, 55.02648485, 52.50591683, 51.61839288, 50.73086893, 49.55933732, 50.00309929, 50.14510313, 48.10379804, 46.23999775, 46.23999775, 46.68375973, 45.97374057, 46.59500733, 46.68375973, 46.77251212, 47.30502649, 50.05635073, 49.96759834, 50.14510313, 50.05635073, 51.56514145, 50.14510313, 53.25143695, 57.33404711, 56.44652317, 56.44652317, 56.80153275, 56.62402796, 56.97903754, 55.38149443, 55.55899922, 53.78395132, 49.08007439, 48.28130283, 48.8138172, 48.28130283],
							[7.36644878, 7.40194974, 7.77470979, 7.63270596, 7.49070213, 7.54395357, 7.66820692, 8.0232165, 8.66223374, 9.40775386, 9.23024907, 9.85151584, 9.76276344, 9.94026823, 9.67401105, 9.67401105, 9.31900147, 9.94026823, 9.85151584, 9.58525865, 9.76276344, 9.31900147, 9.31900147, 9.49650626, 9.49650626, 8.66223374, 9.05274428, 9.05274428, 9.05274428, 9.23024907, 8.82198805, 8.7864871, 8.66223374, 8.7864871, 9.40775386, 8.96399189, 8.73323566, 8.7864871, 7.98771554, 8.34272512, 7.27769638, 7.18894399, 7.27769638, 7.66820692, 7.66820692, 7.81021075, 8.34272512, 8.37822608, 7.93446411, 7.81021075, 7.98771554, 8.66223374, 8.96399189, 9.05274428, 8.25397273, 8.46697847, 8.6977347, 8.6977347, 8.16522033, 8.73323566, 8.87523949, 9.14149668, 8.82198805, 9.14149668, 9.14149668, 9.40775386, 8.82198805, 8.96399189, 8.60898231, 8.34272512, 8.28947369, 7.98771554, 7.72145836, 7.81021075, 8.16522033, 8.60898231, 8.34272512, 8.37822608, 8.66223374, 8.96399189, 9.40775386, 8.96399189, 9.58525865, 10.73903978, 11.53781134, 11.89282092, 11.62656373, 13.75662121, 13.49036403, 15.97543108, 14.82164995, 15.6204215, 15.35416432, 15.79792629, 15.35416432, 15.35416432, 17.92798377, 16.24168827, 15.26541193, 16.59669785, 15.7091739, 15.7091739, 14.82164995, 13.84537361, 13.04660205, 13.13535445, 12.78034487, 12.33658289, 13.31285924, 12.86909726, 12.86909726, 13.22410684, 12.69159247, 12.2478305, 11.62656373, 11.80406852, 12.07032571, 12.07032571, 11.89282092, 11.89282092, 12.07032571, 11.53781134, 11.18280176, 10.68578835, 12.01707427, 11.98157331, 12.78034487, 12.33658289, 12.21232954, 11.92832188, 11.80406852, 12.37208385, 12.90459822, 12.72709343, 12.86909726, 12.86909726, 13.04660205, 13.66786882, 13.31285924, 13.57911642, 13.934126, 13.75662121, 14.20038319, 12.60284008, 11.18280176, 10.73903978, 11.71531613, 11.62656373, 11.36030655, 10.56153499, 9.49650626, 9.76276344, 9.49650626, 10.02902063, 10.82779218, 11.09404936, 10.29527781, 10.29527781, 9.76276344, 9.67401105, 9.76276344, 9.67401105, 9.58525865, 9.58525865, 10.11777302, 9.76276344, 9.85151584, 11.36030655, 11.00529697, 11.36030655, 11.98157331, 11.89282092, 11.80406852, 11.53781134, 11.62656373, 11.80406852, 12.42533529, 11.89282092, 11.89282092, 11.44905894, 9.85151584, 10.11777302, 9.94026823, 9.76276344, 10.11777302, 11.00529697, 9.85151584, 9.31900147, 9.23024907, 9.40775386, 9.05274428, 8.73323566, 8.66223374, 8.46697847, 7.40194974, 7.66820692, 7.57945453, 7.40194974, 7.18894399, 7.13569255, 7.89896315, 7.81021075, 8.7864871, 8.87523949, 9.14149668, 9.23024907, 9.40775386, 10.4727826, 10.11777302, 10.02902063, 9.67401105, 8.6977347, 8.55573087, 7.98771554, 8.55573087, 8.16522033, 8.07646794, 8.21847177, 8.1119689, 7.98771554, 7.89896315, 7.98771554, 8.25397273, 7.98771554, 7.81021075, 7.72145836, 7.45520117, 7.57945453, 7.36644878, 7.63270596, 7.63270596, 7.63270596, 7.77470979, 7.63270596, 7.89896315, 7.33094782, 7.27769638, 7.40194974, 7.57945453, 7.57945453, 7.81021075, 8.7864871, 9.40775386, 9.49650626, 8.96399189, 8.82198805, 8.82198805, 9.23024907, 9.05274428, 8.87523949, 8.82198805, 8.82198805, 8.87523949, 8.73323566, 8.55573087, 8.52022991, 8.25397273, 8.28947369, 8.34272512, 7.72145836, 7.89896315, 8.21847177, 7.93446411, 7.63270596, 7.27769638, 7.27769638, 7.0114392, 7.22444495, 7.13569255, 7.18894399, 7.04694016, 7.18894399, 7.10019159, 7.13569255, 7.66820692, 7.36644878, 7.10019159, 7.63270596, 8.6977347, 7.77470979, 8.7864871, 9.14149668, 8.37822608, 8.60898231, 8.34272512, 8.82198805, 8.46697847, 8.21847177, 7.77470979, 8.28947369, 8.25397273],
							[29.46579511, 27.33573763, 28.75577595, 28.75577595, 27.868252, 26.62571847, 26.98072805, 27.51324242, 29.11078553, 30.88583343, 32.30587175, 31.24084301, 31.4183478, 35.32345318, 34.43592923, 32.66088133, 31.59585259, 32.83838612, 31.77335738, 30.88583343, 31.06333822, 28.40076637, 28.93328074, 28.22326158, 28.04575679, 27.15823284, 27.69074721, 27.868252, 26.80322326, 27.868252, 27.33573763, 26.44821368, 26.44821368, 27.51324242, 31.24084301, 29.46579511, 28.57827116, 29.82080469, 27.868252, 28.75577595, 25.38318495, 25.20568016, 24.85067058, 26.98072805, 26.80322326, 26.80322326, 29.28829032, 29.82080469, 28.93328074, 28.40076637, 28.93328074, 30.88583343, 31.4183478, 31.95086217, 28.93328074, 28.75577595, 30.17581427, 30.53082385, 29.82080469, 31.24084301, 32.83838612, 35.32345318, 36.21097713, 35.50095797, 36.92099628, 39.58356813, 37.09850107, 39.05105376, 37.63101544, 36.03347234, 38.16352981, 37.09850107, 35.85596755, 34.61343402, 37.27600586, 39.05105376, 39.76107292, 40.82610166, 40.29358729, 42.77865435, 44.02118788, 41.00360645, 41.53612082, 42.06863519, 44.02118788, 47.92629325, 48.10379804, 51.65389384, 48.8138172, 53.07393216, 55.02648485, 57.24529472, 55.47024682, 65.67677224, 67.45182014, 67.45182014, 68.33934408, 62.57043841, 61.68291447, 70.11439198, 71.44567791, 74.55201173, 68.78310606, 64.34548631, 59.46410459, 60.35162854, 58.13281867, 53.69519892, 56.80153275, 55.9140088, 57.24529472, 58.57658064, 58.13281867, 56.35777077, 51.12137947, 54.58272287, 59.02034262, 59.46410459, 56.80153275, 56.80153275, 56.35777077, 57.68905669, 55.38149443, 55.02648485, 59.46410459, 61.0616477, 63.01420039, 63.90172434, 65.14425787, 64.25673392, 63.54671476, 68.16183929, 67.80682971, 67.9843345, 66.20928661, 64.78924829, 63.36920997, 62.30418123, 60.35162854, 61.59416207, 62.48168602, 60.52913333, 61.77166686, 55.55899922, 52.363913, 51.65389384, 55.20398964, 56.80153275, 55.55899922, 51.38763666, 45.08621662, 46.86126451, 48.01504565, 49.70134115, 54.31646569, 55.73650401, 57.33404711, 56.09151359, 55.20398964, 59.64160938, 57.86656148, 60.52913333, 59.46410459, 61.41665728, 62.12667644, 59.64160938, 60.35162854, 68.16183929, 67.27431535, 70.29189677, 70.11439198, 64.6117435, 65.85427703, 65.85427703, 66.3867914, 70.46940156, 70.82441114, 69.04936324, 69.22686803, 68.16183929, 67.27431535, 66.03178182, 65.85427703, 66.74180098, 68.69435366, 73.13197341, 72.77696383, 74.37450694, 76.50456442, 75.61704047, 75.97205005, 74.55201173, 76.14955484, 73.66448778, 73.48698299, 74.55201173, 78.27961231, 78.81212668, 79.52214584, 82.0072129, 83.60475601, 89.19615689, 87.68736618, 86.44483265, 88.39738533, 88.75239491, 90.97120479, 89.63991886, 88.57489012, 88.04237576, 87.68736618, 80.40966979, 85.91231828, 86.62233744, 87.50986139, 86.26732786, 84.49227996, 84.49227996, 83.42725122, 82.53972727, 84.66978475, 85.73481349, 88.21988055, 87.68736618, 86.44483265, 85.73481349, 86.44483265, 86.44483265, 84.66978475, 84.66978475, 82.71723206, 83.60475601, 83.60475601, 83.24974643, 85.5573087, 83.07224164, 83.7822608, 86.62233744, 89.63991886, 90.97120479, 90.97120479, 94.07753861, 95.40882453, 94.52130058, 95.85258651, 101.1777302, 106.05911192, 107.39039785, 106.05911192, 104.28406402, 98.95892033, 106.94663587, 108.72168377, 103.84030205, 100.73396823, 109.16544574, 108.72168377, 103.84030205, 106.05911192, 99.84644428, 103.39654008, 104.727826, 103.84030205, 96.29634848, 97.62763441, 103.39654008, 98.07139638, 100.73396823, 101.62149218, 100.73396823, 104.28406402, 111.82801759, 109.16544574, 113.15930352, 114.49058944, 108.72168377, 107.83415982, 114.04682747, 118.92820919, 113.15930352, 118.48444721, 118.48444721, 123.80959091, 137.12245014, 143.33511779, 141.56006989, 135.79116422, 114.49058944, 126.02840078, 126.02840078, 123.80959091],
							[5.55589992, 5.44939705, 5.80440663, 5.50264848, 5.50264848, 5.44939705, 5.32514369, 5.32514369, 5.62690184, 5.71565423, 6.12391525, 6.336921, 6.78068297, 7.45520117, 7.33094782, 7.22444495, 7.27769638, 7.54395357, 7.36644878, 7.27769638, 7.18894399, 6.74518201, 6.69193058, 6.69193058, 6.74518201, 6.336921, 6.2481686, 6.21266764, 6.03516285, 6.12391525, 6.12391525, 5.80440663, 6.07066381, 6.69193058, 6.95818776, 6.95818776, 7.04694016, 6.95818776, 6.51442579, 6.78068297, 5.68015327, 5.94641046, 5.85765806, 5.89315902, 6.03516285, 5.9996619, 6.21266764, 6.44342387, 6.03516285, 6.2481686, 6.336921, 6.74518201, 7.22444495, 7.13569255, 6.9226868, 6.56767722, 6.74518201, 6.51442579, 6.83393441, 7.0114392, 6.88718585, 6.88718585, 6.88718585, 7.18894399, 8.1119689, 8.43147752, 8.25397273, 7.81021075, 7.72145836, 7.54395357, 7.49070213, 7.27769638, 6.88718585, 6.88718585, 7.0114392, 7.0114392, 7.10019159, 7.13569255, 7.18894399, 7.13569255, 7.10019159, 7.04694016, 7.04694016, 7.33094782, 7.36644878, 8.46697847, 9.40775386, 9.76276344, 9.23024907, 11.62656373, 11.44905894, 11.44905894, 11.80406852, 12.86909726, 12.78034487, 13.13535445, 15.6204215, 14.46664037, 13.40161163, 14.64414516, 14.64414516, 15.6204215, 14.82164995, 13.66786882, 12.2478305, 12.2478305, 12.1590781, 10.73903978, 11.80406852, 11.62656373, 11.71531613, 12.1590781, 10.91654457, 10.20652542, 9.49650626, 10.11777302, 11.09404936, 10.73903978, 10.65028739, 10.73903978, 10.73903978, 11.00529697, 10.33077877, 10.29527781, 11.04079793, 11.44905894, 11.92832188, 12.01707427, 11.98157331, 12.10582667, 11.62656373, 12.33658289, 13.17085541, 12.95784966, 13.04660205, 12.01707427, 11.62656373, 11.98157331, 11.71531613, 12.07032571, 12.2478305, 11.80406852, 11.53781134, 9.94026823, 8.82198805, 8.1119689, 9.31900147, 9.23024907, 9.05274428, 8.43147752, 7.84571171, 7.54395357, 7.40194974, 7.54395357, 8.82198805, 9.14149668, 9.40775386, 9.31900147, 8.73323566, 9.40775386, 9.58525865, 10.38403021, 9.85151584, 10.20652542, 10.38403021, 9.76276344, 9.76276344, 11.00529697, 10.91654457, 11.09404936, 11.44905894, 11.36030655, 11.09404936, 11.18280176, 11.27155415, 11.98157331, 11.98157331, 11.44905894, 11.53781134, 11.27155415, 10.73903978, 10.73903978, 10.82779218, 10.38403021, 11.00529697, 11.71531613, 11.09404936, 11.44905894, 11.98157331, 11.89282092, 12.07032571, 11.80406852, 11.18280176, 10.82779218, 10.82779218, 11.36030655, 11.62656373, 11.36030655, 11.36030655, 11.18280176, 11.27155415, 12.33658289, 12.33658289, 12.78034487, 12.78034487, 12.86909726, 12.69159247, 13.13535445, 12.78034487, 12.78034487, 13.04660205, 12.42533529, 12.78034487, 12.60284008, 12.60284008, 12.2478305, 12.33658289, 12.07032571, 11.62656373, 11.36030655, 11.44905894, 11.80406852, 12.1590781, 11.98157331, 11.89282092, 11.98157331, 12.33658289, 12.42533529, 12.95784966, 13.22410684, 12.78034487, 13.04660205, 13.22410684, 13.934126, 14.99915474, 14.73289756, 14.46664037, 14.37788798, 14.91040235, 15.53166911, 15.35416432, 15.6204215, 15.6204215, 15.53166911, 14.99915474, 14.99915474, 14.82164995, 16.68545024, 16.68545024, 16.33044066, 16.06418348, 15.97543108, 16.50794545, 15.88667869, 15.6204215, 16.59669785, 16.33044066, 15.88667869, 15.97543108, 15.6204215, 15.7091739, 15.44291672, 15.53166911, 14.91040235, 13.40161163, 14.20038319, 13.31285924, 13.84537361, 13.75662121, 14.0228784, 15.44291672, 15.6204215, 15.6204215, 15.35416432, 16.06418348, 15.53166911, 14.28913558, 15.7091739, 15.17665953, 14.73289756, 14.99915474, 15.35416432, 15.17665953, 15.97543108, 15.7091739, 16.06418348, 15.88667869, 14.55539277, 15.53166911, 14.99915474, 15.97543108],
							[8.55573087, 8.50247943, 8.73323566, 8.85748901, 9.21249859, 9.15924716, 9.05274428, 8.76873662, 9.56750817, 9.709512, 9.78051392, 9.47875578, 9.99351967, 9.78051392, 9.78051392, 9.40775386, 9.40775386, 9.709512, 9.56750817, 9.47875578, 9.78051392, 9.33675195, 9.85151584, 9.63851009, 9.40775386, 9.33675195, 9.40775386, 9.709512, 9.40775386, 9.63851009, 9.56750817, 9.33675195, 9.19474811, 10.34852925, 10.84554266, 11.20055224, 11.83956948, 11.76856757, 10.84554266, 11.20055224, 10.06452158, 10.27752733, 10.49053308, 10.91654457, 11.20055224, 11.27155415, 11.12955032, 11.55556182, 11.05854841, 10.77454074, 11.34255607, 11.55556182, 12.26558098, 12.12357715, 11.69756565, 11.4845599, 11.76856757, 11.76856757, 12.05257523, 12.62059056, 12.40758481, 12.40758481, 12.1590781, 12.07032571, 12.40758481, 12.40758481, 12.5850896, 12.40758481, 12.31883241, 12.07032571, 12.2478305, 12.07032571, 11.89282092, 11.55556182, 11.73306661, 12.40758481, 12.2478305, 12.2478305, 12.1590781, 12.67384199, 12.5850896, 12.40758481, 12.67384199, 12.67384199, 12.9223487, 14.12938127, 15.15890905, 16.09968444, 15.22991097, 17.37771892, 17.11146174, 16.61444833, 16.61444833, 17.46647132, 17.28896653, 17.62622563, 19.1705173, 18.8332582, 17.62622563, 18.65575341, 19.2592697, 19.43677449, 18.56700102, 17.46647132, 16.77420264, 16.34819114, 16.52569593, 15.49616815, 15.58492055, 16.52569593, 16.33044066, 17.35996845, 17.46647132, 17.0582103, 15.40741576, 16.63219881, 17.1469627, 16.54344641, 16.33044066, 15.92217965, 15.51391863, 15.19441001, 15.10565761, 13.66786882, 14.32463654, 14.64414516, 14.9459033, 14.6973966, 14.6973966, 14.64414516, 14.53764229, 14.99915474, 15.97543108, 15.56717007, 14.89265187, 14.89265187, 14.59089372, 14.83940043, 14.37788798, 14.43113941, 14.59089372, 14.48439085, 14.32463654, 13.66786882, 12.74484391, 11.30705511, 12.69159247, 13.29510876, 12.8335963, 12.17682858, 12.12357715, 11.76856757, 11.66206469, 12.28333146, 13.3483602, 13.2596078, 13.72112025, 14.21813367, 13.80987265, 14.21813367, 14.07612983, 14.18263271, 13.80987265, 13.98737744, 13.80987265, 13.86312409, 13.61461738, 15.15890905, 15.10565761, 15.53166911, 15.7091739, 15.22991097, 15.2831624, 15.35416432, 15.22991097, 15.35416432, 15.65592246, 15.2831624, 15.22991097, 15.0346557, 14.66189564, 14.66189564, 14.73289756, 15.10565761, 15.15890905, 15.65592246, 15.35416432, 15.35416432, 15.58492055, 15.58492055, 15.58492055, 15.35416432, 15.35416432, 15.35416432, 15.2831624, 15.2831624, 15.65592246, 15.53166911, 15.65592246, 15.53166911, 15.35416432, 15.78017582, 16.57894737, 16.95170743, 16.63219881, 17.69722755, 17.87473234, 19.09951539, 18.17649048, 18.79775724, 18.79775724, 18.05223713, 18.62025245, 18.79775724, 19.04626395, 18.85100868, 18.63800293, 18.26524287, 18.26524287, 17.75047898, 17.96348473, 17.75047898, 18.19424096, 18.19424096, 17.89248281, 18.12323904, 18.26524287, 18.19424096, 17.8214809, 18.19424096, 18.19424096, 18.33624479, 18.40724671, 18.19424096, 19.01076299, 18.4959991, 18.93976107, 19.22376874, 19.29477065, 19.22376874, 19.22376874, 19.38352305, 19.82728502, 19.66753071, 20.11129269, 21.06981855, 21.30057478, 22.40110448, 22.70286262, 21.88634059, 21.60233292, 22.11709681, 22.18809873, 21.88634059, 21.88634059, 22.11709681, 21.44257861, 21.60233292, 22.0460949, 21.30057478, 21.60233292, 21.74433675, 22.18809873, 21.81533867, 21.67333484, 21.51358053, 21.15857095, 21.67333484, 22.18809873, 22.01059394, 21.83308915, 22.45435591, 22.18809873, 22.09934633, 23.96314663, 24.31815621, 25.47193734, 27.42449003, 27.51324242, 27.15823284, 26.0932041, 28.66702356, 28.31201398, 29.28829032, 29.55454751, 29.6432999, 29.11078553, 25.56068974, 27.33573763, 27.42449003, 26.98072805],
							[16.41919306, 16.41919306, 16.59669785, 16.59669785, 16.77420264, 16.50794545, 15.79792629, 16.50794545, 17.3954694, 18.28299335, 20.05804125, 21.83308915, 23.96314663, 24.495661, 24.67316579, 23.43063226, 23.07562268, 25.56068974, 24.85067058, 24.85067058, 24.495661, 23.78564184, 22.89811789, 22.18809873, 22.36560352, 20.59055562, 20.76806041, 22.54310831, 22.01059394, 23.78564184, 23.78564184, 23.96314663, 23.78564184, 25.20568016, 26.27070889, 27.33573763, 26.80322326, 26.98072805, 25.38318495, 25.20568016, 22.01059394, 22.01059394, 21.12306999, 22.54310831, 22.01059394, 21.83308915, 23.43063226, 24.495661, 23.07562268, 23.78564184, 24.85067058, 26.44821368, 27.33573763, 27.33573763, 25.56068974, 24.14065142, 25.73819453, 25.73819453, 24.31815621, 25.91569932, 26.0932041, 25.73819453, 26.80322326, 29.28829032, 29.11078553, 29.6432999, 29.11078553, 29.6432999, 29.6432999, 28.75577595, 28.75577595, 29.11078553, 28.57827116, 27.51324242, 28.22326158, 29.99830948, 28.40076637, 29.46579511, 28.75577595, 29.46579511, 30.35331906, 29.11078553, 30.88583343, 32.12836696, 32.83838612, 35.67846276, 36.03347234, 39.76107292, 39.93857771, 42.24613998, 38.51853939, 36.5659867, 38.87354897, 45.26372141, 46.68375973, 49.70134115, 54.58272287, 47.39377888, 45.08621662, 48.99132199, 49.70134115, 55.47024682, 51.65389384, 48.8138172, 44.90871183, 44.19869267, 46.15124536, 44.19869267, 44.90871183, 43.84368309, 43.84368309, 43.84368309, 41.53612082, 40.1160825, 38.3410346, 42.60114956, 44.37619746, 41.53612082, 39.22855855, 37.63101544, 37.80852023, 37.80852023, 36.29972952, 35.05719599, 38.78479658, 38.3410346, 39.31731095, 41.35861603, 40.64859687, 40.47109208, 38.96230137, 41.44736843, 44.10994027, 42.24613998, 40.73734927, 38.3410346, 36.5659867, 36.5659867, 35.05719599, 35.50095797, 35.76721515, 35.05719599, 34.9684436, 31.50710019, 28.13450919, 25.47193734, 27.06948045, 27.24698524, 27.24698524, 25.29443255, 24.58441339, 24.58441339, 24.495661, 24.67316579, 26.89197566, 27.51324242, 27.60199482, 27.69074721, 26.80322326, 29.02203314, 31.24084301, 30.53082385, 29.90955709, 29.90955709, 30.08706188, 29.28829032, 28.75577595, 33.54840528, 32.39462414, 32.39462414, 32.66088133, 33.01589091, 32.21711935, 31.59585259, 31.4183478, 33.1933957, 33.54840528, 33.90341486, 33.45965288, 32.92713851, 29.6432999, 28.75577595, 29.28829032, 28.75577595, 30.26456667, 32.21711935, 31.77335738, 32.92713851, 33.90341486, 35.32345318, 36.29972952, 35.32345318, 35.14594839, 33.72591007, 33.45965288, 33.28214809, 35.50095797, 35.67846276, 34.79093881, 35.94471994, 36.29972952, 38.87354897, 38.07477742, 38.96230137, 41.35861603, 42.77865435, 41.802378, 43.39992111, 41.97988279, 40.55984448, 42.33489237, 39.31731095, 44.46494985, 44.28744506, 44.28744506, 42.15738758, 42.86740674, 42.60114956, 42.24613998, 40.55984448, 41.71362561, 43.04491153, 44.28744506, 43.39992111, 41.71362561, 42.42364477, 42.06863519, 41.8911304, 39.22855855, 39.22855855, 38.16352981, 40.29358729, 41.09235885, 41.71362561, 44.02118788, 42.33489237, 41.00360645, 41.35861603, 41.97988279, 43.04491153, 43.22241632, 43.5774259, 44.19869267, 42.95615914, 45.26372141, 47.39377888, 48.72506481, 49.96759834, 52.45266539, 50.85512229, 51.03262708, 51.21013187, 52.89642737, 52.54141779, 51.74264624, 51.83139863, 51.21013187, 50.05635073, 50.05635073, 48.01504565, 48.54756002, 47.39377888, 46.41750254, 42.77865435, 40.29358729, 41.18111124, 40.47109208, 40.73734927, 39.67232053, 37.89727263, 40.2048349, 41.8911304, 41.8911304, 44.73120704, 45.88498817, 43.48867351, 39.93857771, 42.06863519, 40.55984448, 39.49481574, 41.35861603, 41.71362561, 42.68990195, 44.64245464, 47.48253128, 45.08621662, 42.51239716, 39.93857771, 41.53612082, 41.53612082, 41.35861603],
							[9.95801871, 9.28350051, 8.85748901, 8.76873662, 8.87523949, 8.6977347, 8.43147752, 8.66223374, 9.14149668, 9.14149668, 9.67401105, 9.58525865, 10.56153499, 11.27155415, 10.91654457, 10.73903978, 10.4727826, 11.00529697, 10.65028739, 10.82779218, 10.65028739, 10.02902063, 10.38403021, 10.65028739, 10.38403021, 10.02902063, 10.11777302, 10.11777302, 10.20652542, 10.73903978, 11.09404936, 11.18280176, 11.89282092, 12.51408768, 13.13535445, 13.04660205, 12.86909726, 13.04660205, 11.62656373, 12.1590781, 11.09404936, 10.29527781, 10.20652542, 10.73903978, 11.00529697, 11.18280176, 12.33658289, 11.98157331, 11.36030655, 11.36030655, 11.89282092, 12.60284008, 13.934126, 13.22410684, 12.33658289, 12.78034487, 13.22410684, 13.13535445, 12.86909726, 13.66786882, 13.13535445, 13.84537361, 13.31285924, 13.84537361, 14.0228784, 14.82164995, 14.46664037, 14.37788798, 14.20038319, 13.40161163, 13.40161163, 13.04660205, 12.78034487, 13.04660205, 13.40161163, 14.20038319, 13.75662121, 13.57911642, 14.11163079, 14.82164995, 15.08790714, 14.37788798, 15.08790714, 15.44291672, 15.79792629, 16.33044066, 16.59669785, 17.75047898, 17.4842218, 17.92798377, 17.30671701, 17.66172659, 17.75047898, 18.81550772, 19.88053646, 21.12306999, 22.7206131, 21.12306999, 21.47807957, 22.54310831, 22.36560352, 23.78564184, 22.7206131, 21.83308915, 20.59055562, 19.52552688, 20.05804125, 18.81550772, 19.34802209, 18.81550772, 19.1705173, 18.99301251, 17.92798377, 18.10548856, 17.57297419, 19.52552688, 19.88053646, 18.99301251, 18.28299335, 17.75047898, 17.21796461, 16.50794545, 15.92217965, 15.38966528, 16.27718923, 17.16471318, 17.4842218, 17.53747324, 17.16471318, 16.86295503, 15.92217965, 16.54344641, 18.46049814, 18.10548856, 18.46049814, 17.60847515, 17.21796461, 17.66172659, 16.95170743, 16.77420264, 16.77420264, 16.68545024, 16.06418348, 14.11163079, 13.04660205, 11.80406852, 13.04660205, 12.69159247, 12.42533529, 11.89282092, 11.62656373, 11.62656373, 11.98157331, 11.98157331, 12.60284008, 12.60284008, 12.60284008, 12.86909726, 12.1590781, 12.78034487, 12.78034487, 13.13535445, 11.89282092, 12.33658289, 12.07032571, 11.44905894, 11.36030655, 12.86909726, 12.42533529, 12.42533529, 12.51408768, 12.60284008, 12.33658289, 12.1590781, 12.07032571, 13.04660205, 13.66786882, 12.78034487, 12.86909726, 12.86909726, 12.2478305, 11.98157331, 12.2478305, 12.2478305, 12.51408768, 13.66786882, 12.95784966, 12.95784966, 13.40161163, 13.66786882, 13.75662121, 13.40161163, 13.66786882, 12.78034487, 12.69159247, 12.60284008, 13.49036403, 13.49036403, 13.49036403, 13.84537361, 13.66786882, 14.55539277, 15.08790714, 15.08790714, 16.06418348, 17.3954694, 16.77420264, 16.77420264, 16.50794545, 16.15293587, 16.77420264, 15.6204215, 16.15293587, 16.06418348, 15.97543108, 15.26541193, 15.26541193, 15.35416432, 15.26541193, 14.73289756, 15.44291672, 16.59669785, 17.04045982, 16.68545024, 16.50794545, 16.33044066, 16.59669785, 16.77420264, 16.33044066, 15.79792629, 15.6204215, 16.06418348, 16.06418348, 16.15293587, 16.50794545, 15.97543108, 15.44291672, 16.06418348, 16.15293587, 15.7091739, 16.24168827, 16.77420264, 17.04045982, 16.77420264, 16.86295503, 17.04045982, 18.28299335, 19.70303167, 19.1705173, 19.08176491, 18.90426012, 19.1705173, 19.96928886, 19.52552688, 19.2592697, 19.08176491, 18.01673617, 17.4842218, 17.75047898, 17.57297419, 17.04045982, 16.41919306, 16.41919306, 15.44291672, 14.91040235, 14.55539277, 14.20038319, 13.84537361, 13.66786882, 13.40161163, 14.0228784, 14.73289756, 14.91040235, 15.79792629, 15.6204215, 14.82164995, 13.84537361, 14.82164995, 15.17665953, 14.28913558, 14.73289756, 14.28913558, 14.37788798, 15.79792629, 17.21796461, 15.88667869, 15.26541193, 13.84537361, 13.84537361, 13.66786882, 14.55539277],
							[5.85765806, 6.2481686, 6.51442579, 6.30142004, 6.30142004, 6.39017243, 6.21266764, 6.39017243, 6.39017243, 6.69193058, 6.74518201, 7.04694016, 7.13569255, 7.04694016, 7.22444495, 6.78068297, 6.83393441, 7.36644878, 7.63270596, 7.27769638, 7.54395357, 7.36644878, 7.27769638, 7.45520117, 7.57945453, 7.22444495, 7.10019159, 7.0114392, 6.9226868, 7.13569255, 7.36644878, 7.45520117, 7.10019159, 7.36644878, 7.36644878, 7.33094782, 7.33094782, 7.36644878, 7.40194974, 7.22444495, 7.10019159, 6.74518201, 6.83393441, 7.10019159, 6.88718585, 7.10019159, 7.13569255, 7.27769638, 7.13569255, 7.22444495, 7.33094782, 7.40194974, 7.04694016, 7.54395357, 7.27769638, 7.04694016, 7.04694016, 7.33094782, 7.04694016, 7.10019159, 7.27769638, 7.13569255, 7.10019159, 7.33094782, 7.81021075, 7.98771554, 7.89896315, 7.98771554, 8.28947369, 8.07646794, 8.34272512, 8.37822608, 8.52022991, 8.37822608, 8.46697847, 8.7864871, 8.96399189, 8.87523949, 8.52022991, 9.49650626, 9.49650626, 8.28947369, 8.16522033, 8.60898231, 8.28947369, 7.77470979, 8.34272512, 8.25397273, 8.25397273, 8.25397273, 7.84571171, 7.66820692, 7.84571171, 8.43147752, 7.77470979, 7.98771554, 8.34272512, 8.52022991, 8.25397273, 8.37822608, 8.16522033, 8.37822608, 7.93446411, 7.81021075, 7.27769638, 7.40194974, 7.89896315, 7.72145836, 7.66820692, 7.66820692, 7.72145836, 7.77470979, 8.0232165, 8.28947369, 7.98771554, 8.25397273, 8.21847177, 8.25397273, 8.16522033, 8.0232165, 8.07646794, 7.98771554, 7.98771554, 7.89896315, 8.07646794, 8.07646794, 8.34272512, 8.21847177, 8.3604756, 8.16522033, 8.34272512, 7.98771554, 8.52022991, 8.48472895, 8.20072129, 8.449228, 8.34272512, 8.52022991, 8.52022991, 8.60898231, 8.52022991, 8.34272512, 8.52022991, 8.43147752, 7.93446411, 7.72145836, 8.07646794, 7.98771554, 8.0232165, 8.16522033, 7.77470979, 7.81021075, 7.63270596, 7.72145836, 7.54395357, 7.33094782, 7.72145836, 7.45520117, 7.72145836, 7.93446411, 7.81021075, 7.66820692, 7.81021075, 7.98771554, 7.84571171, 8.16522033, 8.07646794, 8.6977347, 8.82198805, 8.6977347, 8.21847177, 8.34272512, 8.07646794, 8.1119689, 8.25397273, 8.25397273, 8.34272512, 8.1119689, 8.0232165, 7.89896315, 7.54395357, 7.49070213, 7.84571171, 7.63270596, 7.81021075, 8.16522033, 8.1119689, 8.1119689, 8.1119689, 7.98771554, 7.98771554, 7.89896315, 7.98771554, 7.98771554, 8.07646794, 8.0232165, 8.25397273, 8.28947369, 8.37822608, 8.34272512, 8.37822608, 8.55573087, 8.43147752, 8.82198805, 9.40775386, 9.76276344, 9.85151584, 9.58525865, 9.49650626, 9.58525865, 9.58525865, 9.14149668, 9.05274428, 8.96399189, 8.96399189, 8.82198805, 8.87523949, 9.05274428, 9.23024907, 8.96399189, 9.94026823, 10.11777302, 9.94026823, 9.85151584, 9.49650626, 9.94026823, 9.40775386, 9.14149668, 8.87523949, 8.96399189, 9.05274428, 9.67401105, 9.49650626, 9.76276344, 10.38403021, 9.67401105, 9.58525865, 9.67401105, 9.40775386, 9.94026823, 10.65028739, 11.09404936, 11.89282092, 11.09404936, 11.98157331, 13.13535445, 12.1590781, 12.07032571, 12.51408768, 11.71531613, 11.53781134, 11.71531613, 11.36030655, 11.18280176, 11.44905894, 11.44905894, 11.44905894, 11.00529697, 11.27155415, 11.27155415, 13.04660205, 12.42533529, 12.60284008, 12.07032571, 12.1590781, 11.71531613, 11.53781134, 11.71531613, 11.98157331, 11.80406852, 11.80406852, 12.42533529, 12.60284008, 12.51408768, 13.22410684, 12.78034487, 12.86909726, 13.04660205, 13.49036403, 12.78034487, 13.04660205, 13.934126, 13.57911642, 13.934126, 13.934126, 12.78034487, 12.42533529, 11.09404936, 11.44905894, 11.53781134, 11.53781134],
							[1.79279838, 1.81054886, 1.88155077, 2.07680604, 1.97030317, 1.89930125, 1.88155077, 1.84604981, 1.82829934, 2.27206131, 2.20105939, 2.23656035, 2.30756227, 2.69807281, 2.71582328, 2.48506706, 2.37856418, 2.50281754, 2.4495661, 2.39631466, 2.50281754, 2.30756227, 2.23656035, 2.20105939, 2.27206131, 2.09455652, 2.27206131, 2.23656035, 2.13005748, 2.18330891, 2.16555844, 2.112307, 2.09455652, 2.112307, 2.37856418, 2.20105939, 2.27206131, 2.27206131, 2.07680604, 2.16555844, 1.86380029, 1.86380029, 1.88155077, 1.89930125, 1.93480221, 2.04130508, 2.112307, 2.21880987, 2.112307, 2.20105939, 2.23656035, 2.28981179, 2.4495661, 2.94657951, 2.80457568, 2.7868252, 3.14183478, 3.2128367, 3.26608813, 3.53234532, 3.37259101, 3.5500958, 3.85185394, 4.02935873, 3.85185394, 3.81635298, 3.78085202, 3.72760059, 3.63884819, 3.4613434, 3.63884819, 3.63884819, 3.67434915, 3.63884819, 4.08261017, 4.26011496, 4.29561591, 4.38436831, 4.43761975, 4.61512454, 4.4731207, 4.29561591, 4.29561591, 4.5618731, 4.5618731, 5.18313986, 5.32514369, 5.62690184, 5.76890567, 7.04694016, 6.83393441, 7.04694016, 6.69193058, 7.27769638, 7.77470979, 7.63270596, 7.49070213, 7.84571171, 7.04694016, 7.66820692, 6.95818776, 7.04694016, 6.47892483, 6.15941621, 6.21266764, 6.21266764, 5.76890567, 5.85765806, 5.62690184, 5.76890567, 6.07066381, 6.07066381, 6.44342387, 6.15941621, 6.07066381, 7.0114392, 7.36644878, 7.22444495, 7.33094782, 6.95818776, 7.0114392, 7.10019159, 7.36644878, 7.63270596, 7.79246027, 7.65045644, 8.07646794, 8.91074045, 8.67998422, 8.83973853, 8.96399189, 9.40775386, 9.7982644, 10.15327398, 9.99351967, 10.02902063, 9.49650626, 9.76276344, 9.40775386, 9.14149668, 9.58525865, 10.02902063, 10.65028739, 9.94026823, 9.31900147, 8.60898231, 9.14149668, 9.67401105, 10.02902063, 9.58525865, 8.34272512, 8.43147752, 8.25397273, 8.60898231, 9.23024907, 9.05274428, 9.14149668, 9.67401105, 9.05274428, 9.76276344, 9.49650626, 10.02902063, 11.18280176, 11.36030655, 11.00529697, 11.27155415, 11.36030655, 12.86909726, 12.51408768, 12.60284008, 13.49036403, 12.69159247, 12.2478305, 12.51408768, 12.1590781, 12.86909726, 13.22410684, 13.13535445, 14.64414516, 14.20038319, 13.934126, 14.11163079, 14.37788798, 14.91040235, 14.20038319, 14.64414516, 14.28913558, 14.64414516, 16.06418348, 15.97543108, 16.06418348, 15.53166911, 15.79792629, 14.99915474, 14.55539277, 15.08790714, 15.08790714, 14.91040235, 15.26541193, 15.08790714, 15.26541193, 15.44291672, 16.50794545, 17.04045982, 16.95170743, 20.05804125, 18.46049814, 18.90426012, 18.72675533, 20.9455652, 19.79178407, 19.1705173, 18.99301251, 19.52552688, 19.34802209, 18.54925054, 18.28299335, 17.83923138, 18.46049814, 17.92798377, 18.46049814, 18.46049814, 18.90426012, 19.88053646, 19.70303167, 19.96928886, 21.12306999, 21.30057478, 21.38932717, 22.54310831, 20.76806041, 21.74433675, 21.65558436, 21.92184154, 22.27685112, 22.27685112, 20.76806041, 21.12306999, 20.8568128, 20.9455652, 20.50180323, 20.32429844, 20.41305083, 19.61427928, 18.46049814, 18.72675533, 18.10548856, 19.1705173, 19.08176491, 17.92798377, 16.68545024, 17.21796461, 17.3954694, 17.75047898, 16.95170743, 19.34802209, 18.81550772, 18.63800293, 18.99301251, 19.52552688, 20.05804125, 19.2592697, 20.05804125, 19.08176491, 17.21796461, 17.4842218, 17.12921222, 17.83923138, 17.3954694, 16.50794545, 16.68545024, 16.68545024, 16.95170743, 16.77420264, 17.57297419, 17.57297419, 16.59669785, 17.04045982, 17.57297419, 16.68545024, 16.15293587, 16.50794545, 16.24168827, 16.33044066, 17.04045982, 16.50794545, 16.06418348, 13.75662121, 13.40161163, 14.20038319, 14.11163079],
							[14.91040235, 14.99915474, 14.99915474, 14.99915474, 15.08790714, 14.91040235, 14.91040235, 15.17665953, 16.59669785, 16.50794545, 17.57297419, 17.66172659, 17.92798377, 17.30671701, 16.24168827, 15.7091739, 15.08790714, 15.79792629, 15.6204215, 16.06418348, 16.06418348, 15.17665953, 15.44291672, 15.7091739, 16.15293587, 15.44291672, 15.53166911, 15.44291672, 15.08790714, 15.26541193, 14.91040235, 14.64414516, 14.28913558, 15.7091739, 16.24168827, 16.06418348, 17.04045982, 17.3954694, 16.68545024, 17.30671701, 15.44291672, 15.97543108, 15.97543108, 16.68545024, 17.12921222, 16.95170743, 18.10548856, 18.28299335, 17.75047898, 17.30671701, 17.4842218, 18.10548856, 17.92798377, 17.57297419, 17.12921222, 17.12921222, 17.66172659, 17.75047898, 17.3954694, 18.10548856, 17.66172659, 18.28299335, 18.28299335, 19.1705173, 19.88053646, 20.76806041, 19.52552688, 19.1705173, 18.99301251, 18.81550772, 19.34802209, 18.63800293, 18.10548856, 18.63800293, 19.34802209, 20.41305083, 20.05804125, 21.12306999, 21.30057478, 22.54310831, 21.65558436, 22.18809873, 22.18809873, 22.54310831, 24.31815621, 29.6432999, 27.51324242, 29.6432999, 28.04575679, 29.99830948, 27.15823284, 26.44821368, 26.0932041, 28.40076637, 27.69074721, 28.40076637, 28.93328074, 26.98072805, 26.98072805, 28.75577595, 28.04575679, 28.75577595, 27.868252, 26.44821368, 25.20568016, 24.85067058, 24.85067058, 22.7206131, 22.89811789, 22.7206131, 23.96314663, 24.67316579, 26.0932041, 26.27070889, 24.67316579, 26.62571847, 27.15823284, 27.33573763, 26.80322326, 26.62571847, 26.98072805, 26.80322326, 25.47193734, 24.85067058, 26.44821368, 26.71447087, 27.33573763, 28.04575679, 27.69074721, 28.84452835, 29.11078553, 29.55454751, 28.75577595, 28.40076637, 27.51324242, 27.42449003, 26.80322326, 28.66702356, 27.77949961, 27.9570044, 28.31201398, 27.60199482, 27.9570044, 26.71447087, 26.44821368, 24.85067058, 25.73819453, 26.62571847, 26.1819565, 24.85067058, 23.60813705, 24.4069086, 24.495661, 24.58441339, 25.02817537, 24.67316579, 24.58441339, 24.76191818, 24.31815621, 26.0932041, 26.27070889, 26.71447087, 26.35946129, 26.98072805, 27.06948045, 26.89197566, 26.71447087, 28.40076637, 27.42449003, 28.75577595, 29.37704272, 27.77949961, 28.04575679, 27.51324242, 27.15823284, 27.9570044, 27.15823284, 26.00445171, 25.20568016, 24.58441339, 23.60813705, 23.69688944, 24.58441339, 24.76191818, 24.93942297, 25.47193734, 24.76191818, 24.93942297, 24.93942297, 24.85067058, 24.4069086, 23.60813705, 23.96314663, 22.54310831, 22.18809873, 22.7206131, 23.25312747, 24.05189902, 23.87439423, 24.31815621, 24.495661, 25.47193734, 26.1819565, 26.00445171, 25.73819453, 27.77949961, 28.13450919, 29.11078553, 28.13450919, 26.98072805, 27.15823284, 26.71447087, 26.53696608, 27.42449003, 27.06948045, 25.73819453, 25.82694692, 25.64944213, 26.1819565, 26.62571847, 26.71447087, 26.00445171, 25.47193734, 25.11692776, 23.69688944, 23.78564184, 24.67316579, 24.495661, 23.43063226, 22.98687028, 22.6318607, 23.34187986, 22.80936549, 22.80936549, 23.16437507, 22.98687028, 22.6318607, 22.80936549, 23.60813705, 24.31815621, 25.02817537, 24.58441339, 24.22940381, 23.87439423, 24.05189902, 22.89811789, 23.78564184, 23.34187986, 23.78564184, 23.25312747, 22.89811789, 22.7206131, 22.36560352, 21.83308915, 23.60813705, 23.87439423, 22.6318607, 23.78564184, 24.22940381, 23.16437507, 23.96314663, 23.69688944, 24.85067058, 24.85067058, 25.56068974, 23.51938465, 23.34187986, 24.22940381, 23.78564184, 23.07562268, 23.43063226, 25.73819453, 27.60199482, 28.40076637, 30.44207146, 33.72591007, 29.02203314, 32.21711935, 32.83838612, 31.59585259, 32.21711935, 35.76721515, 35.32345318, 35.14594839, 34.43592923, 32.74963372, 31.59585259, 28.93328074, 30.61957625, 30.35331906, 29.7320523],
							[22.54310831, 23.60813705, 24.495661, 24.495661, 23.96314663, 23.78564184, 23.60813705, 24.31815621, 25.02817537, 26.44821368, 26.27070889, 27.868252, 30.17581427, 31.95086217, 31.95086217, 31.06333822, 30.35331906, 31.59585259, 30.70832864, 30.70832864, 30.70832864, 29.11078553, 29.11078553, 27.868252, 27.51324242, 26.27070889, 26.27070889, 26.80322326, 26.80322326, 28.93328074, 28.57827116, 28.40076637, 28.57827116, 31.06333822, 31.24084301, 31.24084301, 31.77335738, 32.30587175, 30.35331906, 30.53082385, 25.73819453, 26.27070889, 26.27070889, 27.15823284, 27.15823284, 27.868252, 30.17581427, 31.24084301, 30.35331906, 30.35331906, 31.06333822, 31.95086217, 32.48337654, 31.95086217, 29.99830948, 30.88583343, 31.77335738, 32.66088133, 31.95086217, 34.79093881, 36.03347234, 34.61343402, 34.61343402, 35.67846276, 36.74349149, 37.45351065, 36.03347234, 37.09850107, 37.27600586, 35.67846276, 36.5659867, 35.50095797, 35.14594839, 33.1933957, 35.32345318, 37.09850107, 36.38848191, 36.92099628, 35.85596755, 35.50095797, 36.74349149, 35.85596755, 37.09850107, 37.27600586, 39.40606334, 44.90871183, 47.21627409, 50.5888651, 50.76636989, 52.89642737, 48.8138172, 48.45880762, 49.70134115, 53.69519892, 58.13281867, 60.79539052, 65.23301026, 59.90786657, 57.24529472, 62.12667644, 61.68291447, 67.45182014, 61.68291447, 55.9140088, 53.69519892, 54.58272287, 55.02648485, 52.18640821, 56.35777077, 52.89642737, 56.80153275, 58.13281867, 55.02648485, 52.89642737, 48.28130283, 52.71892258, 57.24529472, 57.24529472, 53.25143695, 53.07393216, 51.29888426, 50.5888651, 49.79009355, 48.28130283, 55.20398964, 55.02648485, 57.68905669, 55.9140088, 56.09151359, 56.09151359, 55.38149443, 58.39907585, 61.23915249, 59.10909501, 56.62402796, 55.20398964, 53.25143695, 54.1389609, 51.92015103, 51.65389384, 53.78395132, 52.80767497, 53.78395132, 48.63631241, 45.97374057, 40.29358729, 45.08621662, 46.32875015, 46.32875015, 41.26986364, 38.3410346, 39.84982532, 41.18111124, 42.33489237, 45.61873099, 46.32875015, 47.0387693, 44.99746422, 42.42364477, 46.95001691, 44.90871183, 44.81995943, 43.39992111, 44.46494985, 43.39992111, 41.18111124, 41.09235885, 46.41750254, 44.99746422, 44.02118788, 45.61873099, 45.52997859, 45.4412262, 44.64245464, 44.81995943, 46.41750254, 46.23999775, 43.84368309, 43.6661783, 42.42364477, 40.47109208, 40.02733011, 40.64859687, 39.22855855, 41.44736843, 44.64245464, 42.51239716, 42.77865435, 44.81995943, 45.4412262, 46.50625494, 46.15124536, 46.41750254, 43.6661783, 42.42364477, 43.13366393, 45.61873099, 45.52997859, 45.88498817, 46.15124536, 45.70748338, 51.92015103, 50.6776175, 52.71892258, 52.80767497, 55.9140088, 54.1389609, 55.9140088, 56.26901838, 53.96145611, 54.31646569, 50.32260792, 53.78395132, 51.83139863, 52.89642737, 50.76636989, 50.50011271, 50.41136031, 50.32260792, 49.79009355, 50.5888651, 51.92015103, 52.98517976, 50.76636989, 49.70134115, 49.52383636, 49.16882678, 50.14510313, 47.30502649, 47.92629325, 48.28130283, 50.14510313, 50.50011271, 51.47638905, 52.71892258, 51.83139863, 51.83139863, 53.78395132, 55.20398964, 55.38149443, 56.26901838, 57.33404711, 56.97903754, 54.67147527, 57.68905669, 62.12667644, 61.77166686, 66.91930577, 71.00191593, 66.3867914, 66.03178182, 68.51684887, 69.93688719, 66.03178182, 67.45182014, 69.04936324, 64.96675308, 62.12667644, 63.72421955, 58.75408543, 63.90172434, 61.0616477, 59.46410459, 55.73650401, 51.74264624, 52.63017018, 52.80767497, 54.84898006, 53.25143695, 52.18640821, 51.74264624, 54.84898006, 54.31646569, 59.46410459, 61.59416207, 59.81911417, 56.80153275, 61.77166686, 59.64160938, 56.62402796, 57.5115519, 56.26901838, 54.49397048, 58.75408543, 59.64160938, 55.73650401, 52.45266539, 49.79009355, 51.29888426, 48.99132199, 47.83754086],
							[26.44821368, 25.73819453, 26.80322326, 26.0932041, 26.44821368, 26.27070889, 26.0932041, 26.62571847, 29.46579511, 31.06333822, 31.59585259, 31.77335738, 32.83838612, 32.48337654, 32.66088133, 31.24084301, 31.06333822, 33.72591007, 33.54840528, 31.59585259, 31.4183478, 30.17581427, 30.35331906, 28.93328074, 28.93328074, 26.0932041, 27.69074721, 28.57827116, 26.80322326, 29.46579511, 29.82080469, 29.82080469, 29.6432999, 32.12836696, 32.48337654, 32.66088133, 34.08091965, 34.61343402, 31.06333822, 32.48337654, 28.93328074, 27.15823284, 26.98072805, 28.57827116, 29.11078553, 29.82080469, 31.06333822, 30.35331906, 29.11078553, 29.46579511, 29.82080469, 32.48337654, 34.43592923, 35.14594839, 33.1933957, 32.83838612, 34.25842444, 33.54840528, 32.83838612, 35.32345318, 36.38848191, 36.03347234, 36.5659867, 37.98602502, 38.51853939, 40.29358729, 37.45351065, 38.16352981, 37.45351065, 34.9684436, 35.67846276, 34.9684436, 33.72591007, 33.01589091, 33.90341486, 36.5659867, 36.03347234, 36.92099628, 34.9684436, 34.43592923, 34.9684436, 33.72591007, 35.67846276, 36.92099628, 38.16352981, 42.42364477, 42.06863519, 47.57128367, 47.57128367, 50.5888651, 48.99132199, 48.45880762, 48.8138172, 54.58272287, 60.79539052, 65.23301026, 72.77696383, 61.68291447, 59.02034262, 65.23301026, 64.78924829, 68.33934408, 63.45796236, 59.46410459, 55.47024682, 55.02648485, 55.9140088, 50.41136031, 50.41136031, 47.39377888, 48.10379804, 49.16882678, 42.06863519, 42.06863519, 40.1160825, 42.42364477, 45.61873099, 44.55370225, 42.60114956, 42.95615914, 41.71362561, 39.58356813, 37.54226305, 36.38848191, 40.91485406, 40.91485406, 44.28744506, 44.99746422, 45.4412262, 44.81995943, 44.10994027, 47.30502649, 50.14510313, 49.61258876, 47.39377888, 47.30502649, 45.26372141, 45.97374057, 43.75493069, 42.68990195, 43.5774259, 43.48867351, 44.55370225, 39.49481574, 37.18725347, 33.90341486, 36.03347234, 37.18725347, 36.6547391, 33.81466246, 31.95086217, 30.08706188, 29.6432999, 30.97458583, 33.90341486, 34.52468162, 36.6547391, 35.85596755, 34.9684436, 38.16352981, 38.3410346, 37.36475826, 36.6547391, 37.63101544, 36.29972952, 35.67846276, 33.99216725, 42.33489237, 40.64859687, 42.42364477, 43.84368309, 42.51239716, 44.10994027, 44.37619746, 45.70748338, 48.72506481, 49.34633157, 49.08007439, 49.96759834, 50.14510313, 47.39377888, 48.28130283, 50.23385552, 50.23385552, 51.47638905, 55.55899922, 53.25143695, 54.1389609, 55.9140088, 55.9140088, 56.80153275, 52.71892258, 52.98517976, 53.25143695, 51.65389384, 54.67147527, 58.04406627, 56.97903754, 57.33404711, 58.93159022, 59.81911417, 61.77166686, 60.88414291, 61.23915249, 66.03178182, 69.40437282, 64.96675308, 69.7593824, 68.33934408, 66.91930577, 67.09681056, 61.59416207, 65.85427703, 63.90172434, 64.78924829, 61.23915249, 61.41665728, 59.99661896, 59.10909501, 58.93159022, 60.52913333, 62.12667644, 65.14425787, 66.20928661, 62.65919081, 62.65919081, 63.72421955, 65.49926745, 62.30418123, 63.54671476, 62.48168602, 63.19170518, 63.36920997, 65.49926745, 68.51684887, 66.56429619, 66.20928661, 68.33934408, 69.93688719, 71.17942072, 71.00191593, 73.66448778, 76.68206921, 77.56959316, 81.65220332, 81.82970811, 82.53972727, 85.02479433, 92.74625269, 87.3323566, 82.89473685, 86.44483265, 91.85872874, 92.74625269, 89.19615689, 94.07753861, 87.3323566, 84.84728954, 89.19615689, 83.7822608, 86.44483265, 85.20229912, 82.71723206, 74.90702131, 73.13197341, 74.19700215, 73.3094782, 72.42195425, 70.64690635, 73.48698299, 81.11968895, 86.79984223, 86.44483265, 87.86487097, 87.3323566, 85.5573087, 79.87715542, 86.08982307, 82.0072129, 82.18471769, 85.37980391, 85.73481349, 87.86487097, 95.85258651, 104.28406402, 91.41496676, 92.30249071, 85.73481349, 84.84728954, 84.13727038, 80.58717458],
							[21.30057478, 20.9455652, 22.36560352, 21.12306999, 20.9455652, 20.76806041, 20.76806041, 20.59055562, 21.12306999, 21.12306999, 22.54310831, 22.54310831, 24.31815621, 22.89811789, 22.7206131, 21.12306999, 20.9455652, 22.01059394, 20.9455652, 21.12306999, 22.18809873, 21.12306999, 21.83308915, 21.30057478, 20.76806041, 19.1705173, 18.63800293, 19.52552688, 18.81550772, 19.34802209, 19.34802209, 18.81550772, 18.81550772, 19.34802209, 20.41305083, 19.70303167, 19.88053646, 19.70303167, 18.10548856, 17.92798377, 15.7091739, 16.50794545, 16.06418348, 16.41919306, 16.77420264, 16.41919306, 16.68545024, 16.59669785, 16.50794545, 16.50794545, 16.68545024, 17.92798377, 17.75047898, 17.75047898, 17.04045982, 16.06418348, 16.59669785, 16.77420264, 16.59669785, 16.86295503, 17.04045982, 18.81550772, 18.63800293, 18.99301251, 19.52552688, 19.52552688, 19.52552688, 19.1705173, 19.34802209, 17.75047898, 17.57297419, 18.28299335, 17.92798377, 17.30671701, 17.3954694, 17.92798377, 18.10548856, 19.1705173, 18.28299335, 18.28299335, 19.1705173, 18.63800293, 18.81550772, 18.81550772, 19.88053646, 20.76806041, 21.30057478, 22.36560352, 22.01059394, 23.43063226, 21.30057478, 21.12306999, 21.12306999, 22.01059394, 22.89811789, 23.07562268, 26.62571847, 23.60813705, 22.54310831, 23.78564184, 23.60813705, 24.85067058, 24.31815621, 23.78564184, 21.30057478, 21.65558436, 21.47807957, 20.41305083, 19.88053646, 19.34802209, 19.52552688, 20.41305083, 19.70303167, 19.88053646, 19.70303167, 19.88053646, 20.23554604, 19.52552688, 19.52552688, 19.34802209, 19.88053646, 19.52552688, 20.23554604, 21.47807957, 20.9455652, 21.47807957, 22.18809873, 22.54310831, 22.09934633, 22.18809873, 22.27685112, 22.7206131, 23.34187986, 22.98687028, 22.36560352, 21.74433675, 21.12306999, 20.9455652, 20.14679365, 20.05804125, 20.41305083, 19.52552688, 19.70303167, 19.52552688, 19.52552688, 18.46049814, 19.34802209, 19.43677449, 19.96928886, 19.61427928, 18.99301251, 18.72675533, 18.37174575, 18.90426012, 20.9455652, 21.83308915, 22.45435591, 21.92184154, 21.83308915, 21.92184154, 21.38932717, 21.21182238, 20.76806041, 20.8568128, 20.23554604, 19.2592697, 18.10548856, 19.70303167, 19.88053646, 20.05804125, 20.76806041, 20.32429844, 20.32429844, 19.61427928, 20.05804125, 20.50180323, 20.23554604, 19.79178407, 19.52552688, 19.52552688, 20.9455652, 20.41305083, 20.67930802, 20.59055562, 21.38932717, 21.12306999, 20.67930802, 20.8568128, 21.30057478, 21.30057478, 21.47807957, 19.79178407, 20.14679365, 19.2592697, 18.99301251, 19.70303167, 20.8568128, 21.65558436, 21.30057478, 20.9455652, 20.9455652, 22.18809873, 22.54310831, 22.54310831, 23.87439423, 24.93942297, 24.31815621, 24.93942297, 24.76191818, 25.02817537, 25.47193734, 23.16437507, 24.14065142, 24.14065142, 23.96314663, 22.54310831, 22.18809873, 22.89811789, 24.14065142, 24.22940381, 24.4069086, 25.20568016, 25.38318495, 26.00445171, 24.58441339, 24.85067058, 25.20568016, 24.85067058, 23.51938465, 24.31815621, 22.27685112, 22.98687028, 22.89811789, 22.6318607, 23.07562268, 22.80936549, 22.36560352, 22.54310831, 22.36560352, 22.45435591, 22.54310831, 22.6318607, 22.36560352, 21.38932717, 21.21182238, 20.76806041, 21.03431759, 20.8568128, 22.6318607, 21.38932717, 21.03431759, 21.83308915, 21.38932717, 21.38932717, 21.03431759, 21.38932717, 20.8568128, 21.30057478, 21.21182238, 20.9455652, 21.47807957, 21.38932717, 21.21182238, 21.12306999, 20.8568128, 20.59055562, 20.76806041, 21.21182238, 21.12306999, 20.8568128, 20.8568128, 21.12306999, 21.30057478, 20.8568128, 20.67930802, 21.03431759, 22.36560352, 25.73819453, 28.48951877, 24.85067058, 25.56068974, 26.27070889, 26.98072805, 25.73819453, 24.67316579, 24.85067058, 23.16437507, 22.45435591, 22.54310831, 22.54310831, 22.18809873],
							[75.88329765, 75.43953568, 80.3209174, 73.2207258, 71.44567791, 68.33934408, 69.67063001, 73.66448778, 72.77696383, 73.66448778, 77.21458358, 78.5458695, 78.5458695, 84.75853714, 85.64606109, 78.5458695, 80.76467937, 85.20229912, 88.75239491, 86.97734702, 94.96506256, 92.30249071, 95.85258651, 94.07753861, 94.96506256, 89.63991886, 96.74011046, 98.51515836, 94.96506256, 96.74011046, 98.51515836, 98.51515836, 96.74011046, 102.9527781, 110.05296969, 110.94049364, 116.26563734, 116.26563734, 108.2779218, 106.5028739, 97.62763441, 97.62763441, 92.30249071, 95.85258651, 98.51515836, 99.4026823, 106.5028739, 106.5028739, 106.5028739, 108.2779218, 109.16544574, 113.60306549, 115.37811339, 120.70325708, 117.15316129, 116.26563734, 119.81573313, 122.47830498, 121.59078103, 129.57849658, 129.57849658, 126.91592473, 122.47830498, 125.14087683, 127.80344868, 132.24106842, 127.80344868, 129.57849658, 130.46602052, 127.80344868, 132.24106842, 132.24106842, 129.57849658, 126.02840078, 131.35354447, 141.11630791, 145.55392766, 148.21649951, 147.32897556, 146.44145161, 147.32897556, 143.77887976, 142.00383186, 144.66640371, 149.10402346, 155.3166911, 152.65411925, 158.8667869, 153.5416432, 158.8667869, 152.65411925, 153.5416432, 155.3166911, 169.51707429, 177.50478983, 188.15507722, 204.1305083, 189.93012512, 189.93012512, 202.35546041, 200.58041251, 232.53127468, 220.10593939, 209.455652, 197.03031671, 179.27983773, 181.05488563, 164.19193059, 158.8667869, 154.42916715, 158.8667869, 163.30440664, 156.20421505, 151.7665953, 147.32897556, 151.7665953, 160.6418348, 157.091739, 152.65411925, 150.87907135, 154.42916715, 148.21649951, 149.54778543, 145.55392766, 155.3166911, 158.8667869, 165.07945454, 169.96083626, 162.41688269, 156.20421505, 156.20421505, 160.6418348, 165.52321652, 161.08559677, 155.76045308, 153.09788123, 152.21035728, 161.52935874, 160.19807282, 159.75431085, 161.52935874, 160.19807282, 161.97312072, 152.65411925, 148.66026148, 144.66640371, 146.44145161, 150.87907135, 148.21649951, 143.77887976, 138.00997409, 133.57235435, 132.6848304, 134.4598783, 142.89135581, 142.00383186, 145.11016569, 143.33511779, 138.45373607, 150.43530938, 154.87292913, 154.87292913, 156.64797702, 161.52935874, 161.97312072, 159.31054887, 161.08559677, 170.40459824, 168.18578836, 175.72974193, 182.82993352, 181.05488563, 178.39231378, 179.27983773, 176.17350391, 187.26755327, 189.04260117, 181.05488563, 185.49250537, 187.26755327, 180.16736168, 177.50478983, 182.82993352, 184.60498142, 189.04260117, 195.25526881, 190.81764907, 190.81764907, 197.03031671, 197.03031671, 202.35546041, 197.03031671, 202.35546041, 194.36774486, 196.14279276, 197.91784066, 203.24298435, 205.01803225, 205.9055562, 205.01803225, 207.6806041, 210.34317595, 208.56812805, 217.44336754, 222.76851124, 226.31860703, 224.54355913, 227.20613098, 225.43108308, 220.10593939, 220.10593939, 203.24298435, 211.2306999, 205.9055562, 209.455652, 204.1305083, 202.35546041, 202.35546041, 200.58041251, 197.03031671, 198.80536461, 203.24298435, 207.6806041, 209.455652, 205.01803225, 203.24298435, 207.6806041, 219.21841544, 213.89327174, 214.78079569, 219.21841544, 225.43108308, 236.08137047, 239.63146627, 241.40651417, 236.96889442, 233.41879863, 244.06908602, 251.16927761, 254.7193734, 254.7193734, 261.819565, 275.13242424, 277.79499608, 283.12013978, 287.55775952, 284.00766373, 276.01994818, 285.78271163, 280.45756793, 274.24490029, 284.89518768, 291.10785532, 292.88290322, 298.20804691, 324.83376539, 315.07100195, 318.62109774, 326.60881329, 325.72128934, 335.48405278, 335.48405278, 337.25910068, 328.38386118, 313.29595405, 319.50862169, 313.29595405, 326.60881329, 331.04643303, 336.37157673, 350.57195991, 385.18539393, 390.51053762, 395.83568132, 417.1362561, 404.71092081, 408.26101661, 413.5861603, 413.5861603, 427.78654349, 434.88673508, 470.38769305, 465.06254935, 475.71283674, 484.58807623, 482.81302834, 457.96235776, 418.911304, 424.23644769, 408.26101661, 410.03606451],
							[3.23058717, 3.30158909, 3.5500958, 4.11811112, 4.11811112, 4.43761975, 4.38436831, 5.14763891, 4.88138172, 5.00563507, 5.11213795, 5.44939705, 5.80440663, 5.71565423, 5.85765806, 5.71565423, 5.80440663, 6.21266764, 6.74518201, 6.74518201, 6.44342387, 5.76890567, 6.12391525, 6.03516285, 5.89315902, 5.50264848, 5.50264848, 5.62690184, 5.36064465, 5.14763891, 5.18313986, 4.34886735, 5.00563507, 5.2363913, 5.44939705, 4.88138172, 5.2363913, 5.00563507, 4.66837597, 4.66837597, 3.72760059, 3.67434915, 3.85185394, 4.08261017, 4.17136256, 4.02935873, 4.11811112, 4.17136256, 3.90510538, 3.90510538, 4.43761975, 4.91688268, 4.97013412, 4.97013412, 4.66837597, 4.38436831, 4.66837597, 5.00563507, 5.18313986, 5.62690184, 5.50264848, 5.59140088, 5.68015327, 5.71565423, 5.41389609, 5.44939705, 4.82813028, 4.97013412, 4.79262933, 4.52637214, 4.91688268, 4.61512454, 4.61512454, 4.43761975, 4.34886735, 4.79262933, 4.91688268, 4.97013412, 4.79262933, 4.88138172, 5.27189226, 5.14763891, 5.36064465, 5.50264848, 5.76890567, 6.21266764, 6.30142004, 7.63270596, 7.10019159, 8.6977347, 8.6977347, 8.60898231, 9.58525865, 10.20652542, 10.20652542, 9.85151584, 10.29527781, 9.23024907, 9.14149668, 8.87523949, 8.52022991, 8.66223374, 8.7864871, 8.43147752, 7.45520117, 7.22444495, 7.49070213, 6.9226868, 8.0232165, 7.89896315, 7.49070213, 8.21847177, 7.72145836, 7.72145836, 7.22444495, 7.89896315, 8.0232165, 7.93446411, 7.98771554, 8.07646794, 7.98771554, 7.81021075, 8.00546602, 7.65045644, 8.50247943, 8.67998422, 9.05274428, 9.17699763, 9.31900147, 9.17699763, 8.91074045, 9.14149668, 9.62075961, 9.709512, 9.26575003, 9.14149668, 8.43147752, 8.6977347, 8.16522033, 8.37822608, 8.55573087, 8.55573087, 8.6977347, 7.72145836, 7.10019159, 6.74518201, 7.10019159, 7.0114392, 6.78068297, 6.78068297, 5.9996619, 5.89315902, 5.71565423, 5.68015327, 6.44342387, 6.56767722, 6.39017243, 6.2481686, 5.9996619, 6.44342387, 6.60317818, 6.65642962, 6.47892483, 6.51442579, 6.47892483, 6.336921, 6.15941621, 6.83393441, 7.04694016, 7.36644878, 7.72145836, 7.40194974, 7.40194974, 7.36644878, 7.49070213, 7.98771554, 8.1119689, 8.16522033, 8.1119689, 8.16522033, 7.45520117, 7.54395357, 7.84571171, 7.84571171, 7.93446411, 8.37822608, 8.07646794, 8.21847177, 8.34272512, 8.55573087, 8.28947369, 8.07646794, 8.34272512, 8.0232165, 8.1119689, 7.72145836, 8.16522033, 8.25397273, 8.34272512, 8.28947369, 8.25397273, 9.40775386, 10.02902063, 9.85151584, 9.76276344, 9.58525865, 9.40775386, 9.40775386, 9.40775386, 8.96399189, 9.23024907, 8.37822608, 9.05274428, 8.6977347, 9.05274428, 8.73323566, 8.6977347, 8.34272512, 8.37822608, 8.34272512, 8.66223374, 8.66223374, 8.96399189, 8.82198805, 8.66223374, 8.66223374, 8.6977347, 8.52022991, 8.55573087, 8.7864871, 8.7864871, 8.96399189, 9.85151584, 9.94026823, 10.02902063, 9.85151584, 9.76276344, 9.76276344, 9.67401105, 9.67401105, 10.20652542, 9.94026823, 9.94026823, 10.02902063, 9.67401105, 9.49650626, 10.02902063, 10.91654457, 11.44905894, 10.4727826, 10.73903978, 11.09404936, 11.53781134, 13.13535445, 12.69159247, 12.78034487, 12.07032571, 12.33658289, 12.60284008, 11.89282092, 12.51408768, 13.40161163, 12.95784966, 12.51408768, 11.44905894, 11.80406852, 11.00529697, 11.80406852, 12.2478305, 13.13535445, 12.86909726, 12.78034487, 12.95784966, 15.26541193, 18.19424096, 17.30671701, 16.77420264, 19.61427928, 20.67930802, 20.41305083, 19.70303167, 19.61427928, 19.08176491, 19.70303167, 19.70303167, 18.99301251, 18.28299335, 16.59669785, 14.55539277, 16.15293587, 14.99915474],
							[8.12971937, 8.07646794, 8.25397273, 8.16522033, 8.12971937, 8.07646794, 8.25397273, 8.82198805, 8.73323566, 9.21249859, 8.96399189, 9.72726248, 10.40178068, 10.22427589, 9.97576919, 9.72726248, 10.1355235, 9.88701679, 9.97576919, 9.7982644, 9.7982644, 8.6977347, 8.96399189, 8.6977347, 9.21249859, 8.7864871, 9.1237462, 8.96399189, 8.53798039, 9.21249859, 9.4610053, 9.39000338, 9.30125099, 9.63851009, 10.31302829, 10.4727826, 10.82779218, 10.98754649, 9.7982644, 10.1355235, 9.4610053, 8.96399189, 9.21249859, 9.21249859, 9.30125099, 9.30125099, 9.4610053, 9.30125099, 8.96399189, 8.7864871, 9.39000338, 9.21249859, 10.06452158, 9.7982644, 9.4610053, 10.06452158, 9.88701679, 9.88701679, 9.97576919, 9.97576919, 10.40178068, 11.07629889, 10.65028739, 10.73903978, 10.73903978, 11.2360532, 10.98754649, 11.5733123, 12.08807619, 11.9105714, 12.85134678, 12.76259439, 12.60284008, 11.83956948, 12.51408768, 12.42533529, 12.51408768, 12.85134678, 12.60284008, 12.76259439, 12.76259439, 12.60284008, 12.76259439, 12.67384199, 12.85134678, 13.61461738, 14.62639468, 15.54941959, 15.05240618, 16.64994929, 15.38966528, 15.30091288, 15.30091288, 16.56119689, 15.79792629, 16.56119689, 17.07596078, 15.05240618, 14.20038319, 15.30091288, 15.38966528, 15.46066719, 14.11163079, 14.53764229, 14.28913558, 14.78614899, 15.63817198, 14.71514708, 15.21216049, 14.71514708, 14.53764229, 14.71514708, 14.44888989, 13.95187648, 13.52586498, 14.28913558, 14.0228784, 14.28913558, 13.95187648, 13.77437169, 13.6856193, 12.94009918, 13.52586498, 13.13535445, 13.77437169, 12.85134678, 12.85134678, 12.79809535, 12.54958864, 10.98754649, 11.87507044, 11.2360532, 11.62656373, 11.66206469, 11.07629889, 10.98754649, 10.82779218, 11.83956948, 11.41355799, 11.50231038, 11.32480559, 11.07629889, 10.82779218, 10.56153499, 9.05274428, 8.0232165, 9.1237462, 8.87523949, 9.30125099, 8.7864871, 8.07646794, 7.18894399, 7.43745069, 7.73920884, 7.52620309, 7.31319734, 7.40194974, 7.27769638, 6.79843345, 7.10019159, 7.43745069, 7.73920884, 7.89896315, 7.89896315, 8.28947369, 7.95221458, 7.52620309, 8.87523949, 8.6977347, 9.05274428, 9.05274428, 9.05274428, 9.88701679, 9.88701679, 10.40178068, 10.56153499, 11.07629889, 10.82779218, 10.98754649, 11.16505128, 10.56153499, 10.22427589, 10.40178068, 10.65028739, 10.40178068, 10.4727826, 10.4727826, 10.31302829, 10.1355235, 10.56153499, 10.1355235, 10.1355235, 10.22427589, 9.88701679, 9.05274428, 8.6977347, 8.96399189, 9.21249859, 9.21249859, 8.87523949, 9.21249859, 9.7982644, 9.7982644, 9.63851009, 10.06452158, 10.22427589, 9.7982644, 9.88701679, 9.63851009, 9.72726248, 9.7982644, 9.05274428, 9.30125099, 9.1237462, 9.39000338, 9.4610053, 9.4610053, 9.4610053, 9.30125099, 8.73323566, 8.6977347, 9.05274428, 9.14149668, 8.96399189, 8.52022991, 8.60898231, 8.43147752, 8.43147752, 8.25397273, 8.21847177, 8.07646794, 8.34272512, 8.46697847, 8.16522033, 8.1119689, 8.07646794, 7.84571171, 7.89896315, 7.89896315, 9.14149668, 8.52022991, 9.23024907, 9.31900147, 9.14149668, 9.14149668, 9.40775386, 9.23024907, 9.40775386, 9.58525865, 9.14149668, 9.14149668, 9.23024907, 9.23024907, 9.14149668, 9.14149668, 9.58525865, 9.49650626, 9.40775386, 10.02902063, 9.58525865, 9.49650626, 9.67401105, 9.76276344, 9.40775386, 9.14149668, 8.96399189, 9.05274428, 8.82198805, 8.73323566, 8.7864871, 8.87523949, 8.7864871, 8.60898231, 8.87523949, 9.05274428, 8.96399189, 8.46697847, 8.60898231, 8.43147752, 8.34272512, 8.46697847, 8.34272512, 8.0232165, 8.21847177, 8.60898231, 8.34272512, 8.16522033, 7.77470979, 7.89896315, 8.21847177, 8.82198805],
							[49.70134115, 49.34633157, 51.12137947, 49.34633157, 50.5888651, 50.76636989, 51.12137947, 52.18640821, 55.47024682, 56.80153275, 60.35162854, 59.02034262, 61.68291447, 64.78924829, 62.12667644, 59.02034262, 61.23915249, 65.23301026, 65.23301026, 63.01420039, 64.34548631, 57.24529472, 59.90786657, 62.57043841, 61.68291447, 56.35777077, 59.46410459, 62.12667644, 59.90786657, 64.34548631, 64.34548631, 64.78924829, 63.45796236, 67.45182014, 65.23301026, 66.56429619, 65.67677224, 64.34548631, 58.13281867, 57.24529472, 49.52383636, 49.87884594, 48.28130283, 51.65389384, 52.00890342, 52.18640821, 55.9140088, 55.02648485, 54.1389609, 54.1389609, 54.58272287, 58.57658064, 59.90786657, 61.68291447, 55.9140088, 55.9140088, 58.13281867, 60.35162854, 61.23915249, 65.23301026, 66.12053421, 67.89558211, 67.89558211, 68.33934408, 70.11439198, 73.2207258, 70.55815396, 73.66448778, 70.55815396, 68.78310606, 73.66448778, 71.44567791, 67.89558211, 64.34548631, 66.56429619, 70.11439198, 69.22686803, 71.44567791, 65.23301026, 67.00805816, 70.11439198, 68.33934408, 71.00191593, 71.88943988, 77.21458358, 82.53972727, 82.0959653, 91.41496676, 88.30863294, 98.51515836, 94.07753861, 96.74011046, 93.19001466, 109.16544574, 110.05296969, 112.71554154, 123.36582893, 105.61534995, 104.727826, 113.60306549, 111.82801759, 118.04068524, 113.60306549, 102.06525415, 94.96506256, 95.85258651, 94.07753861, 93.19001466, 100.29020625, 94.96506256, 99.4026823, 102.06525415, 100.29020625, 98.51515836, 94.07753861, 102.9527781, 103.84030205, 106.5028739, 103.84030205, 103.84030205, 102.06525415, 101.1777302, 98.07139638, 98.95892033, 104.727826, 106.05911192, 112.71554154, 116.70939931, 104.28406402, 106.5028739, 110.49673167, 112.27177957, 115.37811339, 112.71554154, 110.05296969, 107.39039785, 103.39654008, 102.50901613, 97.18387243, 101.1777302, 101.62149218, 101.62149218, 100.73396823, 92.30249071, 87.68736618, 80.05466021, 83.07224164, 84.84728954, 85.5573087, 77.03707879, 72.24444946, 71.00191593, 73.66448778, 77.39208837, 83.60475601, 87.68736618, 91.41496676, 92.74625269, 90.08368084, 95.85258651, 93.19001466, 93.63377663, 93.63377663, 95.40882453, 96.74011046, 91.85872874, 89.19615689, 100.73396823, 96.74011046, 102.9527781, 108.2779218, 103.84030205, 105.61534995, 102.9527781, 104.727826, 109.60920772, 110.49673167, 106.94663587, 107.39039785, 106.5028739, 100.73396823, 99.84644428, 102.06525415, 104.28406402, 106.5028739, 109.16544574, 107.39039785, 108.72168377, 112.27177957, 112.71554154, 108.72168377, 102.50901613, 104.28406402, 99.4026823, 95.40882453, 98.07139638, 105.17158797, 106.5028739, 105.17158797, 106.94663587, 106.5028739, 111.82801759, 110.94049364, 113.15930352, 118.48444721, 125.5846388, 122.92206696, 127.3596867, 125.5846388, 121.14701906, 120.70325708, 115.37811339, 120.25949511, 120.70325708, 120.25949511, 112.27177957, 111.82801759, 110.94049364, 113.15930352, 113.15930352, 114.49058944, 120.70325708, 122.03454301, 119.81573313, 113.60306549, 114.04682747, 117.59692326, 120.25949511, 118.04068524, 119.81573313, 117.59692326, 119.37197116, 118.92820919, 118.92820919, 121.59078103, 122.03454301, 118.48444721, 121.14701906, 120.70325708, 121.59078103, 120.70325708, 124.69711486, 127.3596867, 123.80959091, 120.25949511, 120.70325708, 122.03454301, 133.57235435, 130.02225855, 126.47216275, 125.5846388, 130.02225855, 130.9097825, 129.57849658, 125.5846388, 127.3596867, 130.9097825, 126.47216275, 126.91592473, 115.82187536, 119.37197116, 118.04068524, 117.59692326, 109.16544574, 103.84030205, 108.2779218, 107.83415982, 106.94663587, 104.28406402, 104.727826, 106.05911192, 111.38425562, 110.94049364, 115.37811339, 115.37811339, 121.59078103, 114.93435141, 122.47830498, 123.80959091, 118.04068524, 122.03454301, 120.70325708, 122.03454301, 130.46602052, 127.80344868, 121.14701906, 115.82187536, 105.17158797, 114.49058944, 111.82801759, 108.72168377],
							[5.85765806, 5.76890567, 6.01741238, 6.01741238, 5.89315902, 5.76890567, 5.68015327, 6.17716669, 7.06469064, 7.18894399, 8.23622225, 8.39597656, 9.28350051, 9.67401105, 9.51425673, 8.7864871, 8.62673279, 9.76276344, 9.35450242, 9.19474811, 9.85151584, 9.0349938, 9.35450242, 9.28350051, 9.0349938, 8.71548518, 8.7864871, 8.71548518, 8.55573087, 8.62673279, 8.87523949, 8.46697847, 8.7864871, 9.0349938, 9.60300913, 9.35450242, 9.51425673, 9.28350051, 8.23622225, 8.39597656, 6.53217627, 6.37242195, 6.78068297, 7.06469064, 7.54395357, 7.45520117, 7.89896315, 8.14746985, 7.57945453, 7.50845261, 7.79246027, 8.30722416, 8.71548518, 8.7864871, 8.14746985, 7.63270596, 7.79246027, 7.73920884, 7.63270596, 8.30722416, 8.30722416, 8.96399189, 8.87523949, 9.28350051, 9.35450242, 10.01127015, 9.67401105, 10.17102446, 10.17102446, 9.67401105, 9.92251775, 9.51425673, 9.19474811, 8.55573087, 8.96399189, 9.67401105, 9.28350051, 9.28350051, 9.0349938, 9.19474811, 9.28350051, 8.96399189, 9.0349938, 9.1237462, 9.35450242, 10.33077877, 10.33077877, 11.28930463, 11.28930463, 11.94607236, 11.37805703, 11.37805703, 11.21830272, 12.90459822, 15.17665953, 15.97543108, 19.52552688, 18.72675533, 18.38949623, 19.52552688, 18.54925054, 19.52552688, 17.91023329, 15.88667869, 13.79212217, 14.3601375, 14.76839851, 12.99335062, 13.95187648, 13.47261355, 13.56136594, 13.95187648, 12.74484391, 12.51408768, 11.69756565, 12.90459822, 13.31285924, 13.15310493, 12.74484391, 12.74484391, 12.42533529, 12.51408768, 11.85731996, 11.00529697, 12.42533529, 12.17682858, 12.74484391, 13.06435253, 12.8335963, 12.86909726, 12.46083625, 13.06435253, 13.50811451, 13.15310493, 13.40161163, 14.44888989, 14.11163079, 14.28913558, 14.28913558, 14.28913558, 14.44888989, 14.67964612, 14.20038319, 12.69159247, 10.91654457, 10.02902063, 10.91654457, 10.73903978, 10.65028739, 9.85151584, 8.82198805, 8.7864871, 9.49650626, 9.94026823, 10.56153499, 9.76276344, 10.56153499, 10.65028739, 10.11777302, 10.38403021, 10.11777302, 10.02902063, 9.76276344, 9.76276344, 9.58525865, 9.49650626, 9.40775386, 10.56153499, 10.02902063, 10.11777302, 10.20652542, 10.02902063, 10.02902063, 10.02902063, 10.02902063, 10.20652542, 10.56153499, 10.20652542, 10.20652542, 9.94026823, 9.67401105, 9.49650626, 9.58525865, 9.58525865, 9.76276344, 9.94026823, 9.67401105, 9.40775386, 10.29527781, 10.4727826, 10.91654457, 10.38403021, 10.73903978, 10.20652542, 10.02902063, 10.20652542, 10.91654457, 10.56153499, 10.29527781, 10.4727826, 10.65028739, 11.98157331, 12.42533529, 13.57911642, 14.73289756, 15.08790714, 14.73289756, 15.17665953, 14.99915474, 14.55539277, 14.55539277, 12.95784966, 13.66786882, 13.66786882, 14.0228784, 13.40161163, 13.66786882, 13.57911642, 13.13535445, 12.86909726, 13.22410684, 13.66786882, 13.934126, 13.84537361, 13.66786882, 13.49036403, 13.57911642, 14.0228784, 13.57911642, 13.49036403, 13.22410684, 13.22410684, 13.66786882, 13.934126, 14.82164995, 14.28913558, 13.84537361, 14.46664037, 14.46664037, 14.91040235, 15.26541193, 15.35416432, 15.53166911, 15.17665953, 15.44291672, 15.88667869, 16.33044066, 17.57297419, 17.4842218, 16.86295503, 17.04045982, 17.04045982, 16.86295503, 16.41919306, 15.6204215, 16.50794545, 16.15293587, 16.15293587, 16.15293587, 15.44291672, 15.88667869, 15.88667869, 15.35416432, 14.55539277, 14.55539277, 14.20038319, 13.22410684, 13.40161163, 12.78034487, 12.60284008, 14.11163079, 14.73289756, 14.46664037, 14.99915474, 16.06418348, 15.88667869, 14.91040235, 15.17665953, 14.91040235, 13.934126, 14.0228784, 13.31285924, 13.22410684, 14.64414516, 14.91040235, 14.11163079, 13.57911642, 13.04660205, 12.78034487, 12.78034487, 12.78034487],
							[16.95170743, 16.41919306, 16.41919306, 16.50794545, 16.50794545, 16.86295503, 16.50794545, 16.50794545, 17.3954694, 17.92798377, 18.99301251, 18.99301251, 19.88053646, 22.18809873, 22.01059394, 20.76806041, 19.88053646, 21.65558436, 20.9455652, 20.41305083, 20.41305083, 18.63800293, 18.46049814, 18.46049814, 18.46049814, 17.66172659, 17.3954694, 17.75047898, 17.4842218, 17.66172659, 18.28299335, 18.46049814, 19.34802209, 20.23554604, 21.47807957, 21.30057478, 20.76806041, 20.59055562, 18.28299335, 18.46049814, 15.53166911, 14.99915474, 15.26541193, 15.97543108, 16.41919306, 15.7091739, 16.77420264, 16.77420264, 16.59669785, 16.59669785, 17.21796461, 18.99301251, 20.05804125, 20.41305083, 18.99301251, 18.10548856, 19.1705173, 18.46049814, 18.10548856, 19.70303167, 19.70303167, 20.9455652, 20.59055562, 22.01059394, 23.07562268, 23.43063226, 22.18809873, 22.54310831, 22.01059394, 20.76806041, 21.30057478, 20.05804125, 20.23554604, 20.23554604, 20.9455652, 21.65558436, 20.59055562, 21.65558436, 20.9455652, 21.12306999, 21.30057478, 20.9455652, 20.9455652, 21.47807957, 22.36560352, 24.495661, 24.14065142, 26.44821368, 26.27070889, 27.868252, 25.38318495, 24.67316579, 24.85067058, 28.04575679, 29.6432999, 30.53082385, 33.54840528, 33.1933957, 32.12836696, 34.25842444, 33.01589091, 33.90341486, 32.12836696, 29.28829032, 27.868252, 27.69074721, 28.04575679, 26.0932041, 26.27070889, 26.27070889, 26.0932041, 26.27070889, 24.67316579, 23.60813705, 21.65558436, 23.78564184, 24.85067058, 25.02817537, 23.60813705, 23.25312747, 22.89811789, 21.83308915, 20.8568128, 21.12306999, 23.60813705, 23.78564184, 24.495661, 25.29443255, 24.85067058, 24.05189902, 23.69688944, 25.20568016, 26.62571847, 25.64944213, 25.20568016, 25.29443255, 24.4069086, 24.76191818, 24.14065142, 24.31815621, 24.4069086, 24.58441339, 23.51938465, 20.9455652, 20.14679365, 19.1705173, 20.05804125, 20.05804125, 19.52552688, 18.46049814, 18.46049814, 17.66172659, 17.21796461, 17.57297419, 19.1705173, 19.52552688, 19.88053646, 19.70303167, 18.99301251, 20.41305083, 20.23554604, 20.9455652, 19.88053646, 19.96928886, 20.41305083, 20.23554604, 20.14679365, 22.09934633, 21.92184154, 22.18809873, 22.36560352, 22.27685112, 22.45435591, 21.65558436, 21.83308915, 23.16437507, 23.51938465, 22.6318607, 22.6318607, 22.09934633, 20.67930802, 21.03431759, 21.03431759, 20.9455652, 21.74433675, 23.34187986, 22.27685112, 22.09934633, 22.89811789, 22.89811789, 22.7206131, 22.54310831, 22.89811789, 20.9455652, 20.9455652, 20.8568128, 22.01059394, 21.12306999, 21.47807957, 21.47807957, 21.83308915, 23.34187986, 23.60813705, 25.56068974, 26.53696608, 27.69074721, 26.62571847, 27.42449003, 26.62571847, 26.00445171, 26.62571847, 25.29443255, 26.1819565, 26.1819565, 26.27070889, 26.0932041, 25.11692776, 25.38318495, 25.29443255, 24.76191818, 25.73819453, 25.82694692, 26.71447087, 26.27070889, 25.82694692, 25.73819453, 25.82694692, 26.00445171, 24.93942297, 25.11692776, 25.38318495, 26.27070889, 26.1819565, 26.0932041, 26.53696608, 25.91569932, 25.47193734, 26.27070889, 26.35946129, 25.82694692, 26.27070889, 26.98072805, 27.15823284, 26.71447087, 27.24698524, 27.60199482, 28.48951877, 29.46579511, 30.61957625, 29.55454751, 28.75577595, 29.82080469, 30.17581427, 29.37704272, 28.04575679, 29.82080469, 28.31201398, 28.22326158, 28.31201398, 27.868252, 27.69074721, 27.69074721, 28.22326158, 26.53696608, 25.02817537, 24.93942297, 25.47193734, 25.29443255, 25.11692776, 25.20568016, 25.64944213, 26.44821368, 26.0932041, 26.62571847, 26.35946129, 26.00445171, 25.20568016, 25.47193734, 25.20568016, 24.85067058, 24.85067058, 24.85067058, 25.11692776, 27.60199482, 27.60199482, 26.0932041, 25.82694692, 24.85067058, 25.02817537, 24.85067058, 26.0932041],
							[38.87354897, 37.27600586, 39.05105376, 41.18111124, 39.40606334, 37.27600586, 37.45351065, 36.92099628, 39.58356813, 41.53612082, 42.42364477, 42.42364477, 42.60114956, 47.39377888, 47.0387693, 45.08621662, 44.19869267, 47.39377888, 45.79623578, 44.55370225, 45.26372141, 43.13366393, 43.48867351, 42.06863519, 42.06863519, 38.69604418, 39.22855855, 39.76107292, 39.05105376, 38.87354897, 37.63101544, 36.74349149, 36.5659867, 38.69604418, 42.60114956, 40.82610166, 40.47109208, 40.64859687, 36.92099628, 37.63101544, 33.01589091, 31.95086217, 31.4183478, 33.1933957, 33.37090049, 33.1933957, 34.25842444, 36.5659867, 34.43592923, 35.85596755, 36.38848191, 38.69604418, 39.76107292, 41.35861603, 39.22855855, 38.87354897, 40.1160825, 40.1160825, 39.22855855, 40.64859687, 43.48867351, 46.32875015, 47.74878846, 47.92629325, 49.16882678, 49.70134115, 48.8138172, 49.16882678, 48.99132199, 46.50625494, 48.8138172, 46.86126451, 44.90871183, 43.48867351, 44.73120704, 49.34633157, 49.52383636, 50.23385552, 49.52383636, 50.05635073, 49.52383636, 47.92629325, 48.28130283, 48.8138172, 50.94387468, 59.02034262, 59.02034262, 64.78924829, 62.12667644, 67.89558211, 62.57043841, 64.34548631, 63.45796236, 72.33320186, 75.43953568, 78.10210752, 83.87101319, 80.76467937, 78.5458695, 83.87101319, 84.31477517, 88.75239491, 84.75853714, 78.98963147, 75.43953568, 77.65834555, 75.43953568, 66.56429619, 72.33320186, 70.11439198, 71.44567791, 71.00191593, 64.78924829, 64.78924829, 62.57043841, 67.00805816, 70.55815396, 67.00805816, 67.00805816, 66.12053421, 64.34548631, 61.68291447, 58.39907585, 54.67147527, 62.12667644, 62.48168602, 67.27431535, 67.9843345, 67.09681056, 65.14425787, 64.78924829, 70.11439198, 72.24444946, 70.82441114, 68.51684887, 66.74180098, 64.6117435, 65.85427703, 63.90172434, 64.6117435, 65.32176266, 65.14425787, 65.32176266, 59.2865998, 55.20398964, 54.1389609, 56.62402796, 57.5115519, 55.9140088, 51.92015103, 46.68375973, 46.15124536, 47.92629325, 50.14510313, 55.20398964, 56.62402796, 59.10909501, 57.68905669, 55.73650401, 61.23915249, 57.15654233, 59.81911417, 60.35162854, 62.30418123, 62.12667644, 57.86656148, 58.39907585, 68.16183929, 65.85427703, 68.16183929, 68.87185845, 65.49926745, 68.69435366, 68.69435366, 67.9843345, 71.88943988, 73.13197341, 70.29189677, 69.22686803, 68.51684887, 66.74180098, 66.03178182, 67.09681056, 67.9843345, 71.00191593, 75.43953568, 72.95446862, 74.72951652, 77.74709795, 77.56959316, 76.859574, 76.50456442, 77.74709795, 74.90702131, 74.37450694, 75.61704047, 78.63462189, 78.27961231, 79.34464105, 80.76467937, 83.60475601, 90.52744281, 90.52744281, 91.85872874, 100.29020625, 99.84644428, 96.29634848, 101.1777302, 98.07139638, 95.40882453, 95.40882453, 91.41496676, 96.74011046, 96.74011046, 101.1777302, 97.62763441, 94.96506256, 94.52130058, 94.96506256, 91.41496676, 94.52130058, 97.62763441, 101.62149218, 99.84644428, 97.18387243, 97.18387243, 98.95892033, 101.62149218, 96.29634848, 96.29634848, 93.19001466, 97.18387243, 97.62763441, 99.4026823, 100.73396823, 96.29634848, 96.29634848, 98.95892033, 100.73396823, 102.06525415, 105.17158797, 107.39039785, 110.05296969, 108.72168377, 110.05296969, 115.82187536, 117.59692326, 121.14701906, 120.70325708, 118.92820919, 116.70939931, 117.15316129, 123.36582893, 120.25949511, 126.47216275, 134.90364027, 130.9097825, 128.24721065, 134.4598783, 127.3596867, 127.3596867, 131.35354447, 132.24106842, 122.47830498, 121.59078103, 121.14701906, 114.49058944, 116.26563734, 119.37197116, 120.25949511, 125.5846388, 132.24106842, 137.12245014, 136.67868817, 140.67254594, 135.34740224, 128.69097263, 138.45373607, 135.79116422, 127.3596867, 137.56621212, 132.24106842, 134.4598783, 154.87292913, 170.84836021, 163.30440664, 155.3166911, 145.55392766, 149.54778543, 149.99154741, 149.99154741],
							[29.28829032, 28.93328074, 28.57827116, 27.868252, 28.93328074, 28.22326158, 28.22326158, 28.75577595, 30.53082385, 29.99830948, 29.99830948, 31.95086217, 32.83838612, 32.66088133, 34.43592923, 33.72591007, 32.48337654, 34.43592923, 34.79093881, 34.08091965, 34.61343402, 32.48337654, 32.48337654, 32.12836696, 31.24084301, 29.82080469, 30.17581427, 29.99830948, 28.93328074, 28.93328074, 29.11078553, 27.33573763, 27.69074721, 29.99830948, 30.88583343, 30.70832864, 31.59585259, 31.06333822, 28.75577595, 29.99830948, 26.27070889, 26.98072805, 26.80322326, 27.69074721, 27.69074721, 27.15823284, 28.22326158, 29.6432999, 28.22326158, 27.69074721, 28.75577595, 29.82080469, 30.88583343, 31.4183478, 29.46579511, 28.57827116, 29.46579511, 29.82080469, 29.99830948, 31.24084301, 31.4183478, 31.4183478, 30.35331906, 31.4183478, 33.01589091, 33.01589091, 31.95086217, 31.95086217, 30.88583343, 31.4183478, 31.59585259, 31.24084301, 31.24084301, 30.35331906, 31.24084301, 33.01589091, 32.48337654, 34.61343402, 34.79093881, 36.5659867, 36.5659867, 35.67846276, 36.92099628, 37.63101544, 38.87354897, 42.06863519, 41.35861603, 44.37619746, 42.06863519, 48.28130283, 46.86126451, 45.97374057, 45.26372141, 47.39377888, 48.10379804, 50.94387468, 57.24529472, 50.05635073, 49.34633157, 50.5888651, 52.54141779, 56.35777077, 51.29888426, 48.99132199, 42.77865435, 42.24613998, 41.53612082, 37.27600586, 39.40606334, 37.80852023, 38.3410346, 41.00360645, 41.71362561, 40.47109208, 39.76107292, 41.35861603, 43.6661783, 42.95615914, 42.77865435, 43.13366393, 43.6661783, 42.60114956, 40.29358729, 39.93857771, 41.8911304, 41.62487321, 41.53612082, 41.71362561, 43.22241632, 43.48867351, 45.26372141, 47.21627409, 47.30502649, 45.52997859, 46.41750254, 44.99746422, 42.24613998, 43.31116872, 41.35861603, 41.44736843, 41.97988279, 41.802378, 41.62487321, 37.18725347, 34.9684436, 33.81466246, 37.80852023, 38.51853939, 37.54226305, 36.47723431, 37.54226305, 36.29972952, 36.5659867, 37.63101544, 40.55984448, 40.55984448, 40.64859687, 40.02733011, 40.47109208, 44.02118788, 44.02118788, 43.93243548, 43.6661783, 44.64245464, 45.3524738, 42.15738758, 41.00360645, 44.99746422, 46.15124536, 48.54756002, 49.52383636, 47.48253128, 48.8138172, 47.92629325, 46.68375973, 47.92629325, 49.61258876, 47.30502649, 49.34633157, 48.19255044, 46.15124536, 45.17496901, 45.88498817, 47.21627409, 47.74878846, 48.37005523, 46.50625494, 45.88498817, 48.10379804, 48.28130283, 47.57128367, 46.23999775, 47.57128367, 45.88498817, 45.3524738, 45.70748338, 46.32875015, 44.37619746, 45.08621662, 43.84368309, 44.99746422, 46.06249296, 45.97374057, 46.68375973, 45.88498817, 47.1275217, 48.01504565, 48.28130283, 47.48253128, 46.50625494, 47.21627409, 44.81995943, 45.17496901, 44.64245464, 45.3524738, 44.28744506, 44.19869267, 44.28744506, 43.48867351, 42.51239716, 43.48867351, 42.95615914, 44.02118788, 43.5774259, 42.42364477, 41.8911304, 41.8911304, 43.93243548, 42.51239716, 41.97988279, 40.55984448, 41.18111124, 40.91485406, 40.55984448, 42.06863519, 41.00360645, 41.26986364, 41.62487321, 42.60114956, 44.46494985, 43.6661783, 44.02118788, 44.37619746, 44.19869267, 44.28744506, 43.5774259, 43.75493069, 43.75493069, 44.02118788, 45.70748338, 45.4412262, 45.3524738, 45.79623578, 44.81995943, 48.01504565, 48.37005523, 47.57128367, 48.9025696, 49.52383636, 47.74878846, 47.74878846, 47.57128367, 47.39377888, 48.37005523, 48.63631241, 48.54756002, 48.10379804, 48.8138172, 48.10379804, 47.74878846, 46.59500733, 48.99132199, 49.08007439, 48.10379804, 49.52383636, 50.5888651, 52.2751606, 58.04406627, 55.38149443, 53.25143695, 53.42894174, 54.84898006, 56.62402796, 55.9140088, 55.73650401, 53.60644653, 50.05635073, 48.10379804, 50.41136031, 50.50011271, 50.50011271],
							[12.07032571, 12.07032571, 12.78034487, 12.60284008, 13.49036403, 13.13535445, 12.86909726, 13.84537361, 14.91040235, 14.46664037, 14.37788798, 14.37788798, 14.20038319, 14.82164995, 15.08790714, 16.33044066, 15.97543108, 16.15293587, 16.24168827, 17.3954694, 17.75047898, 16.68545024, 17.4842218, 16.33044066, 16.41919306, 16.06418348, 15.97543108, 16.24168827, 16.15293587, 15.26541193, 16.41919306, 16.59669785, 16.59669785, 17.04045982, 17.57297419, 17.66172659, 18.10548856, 21.30057478, 19.34802209, 20.41305083, 18.99301251, 19.88053646, 19.52552688, 20.23554604, 20.9455652, 20.76806041, 23.96314663, 21.65558436, 22.54310831, 22.36560352, 23.96314663, 25.20568016, 26.44821368, 26.80322326, 26.0932041, 29.11078553, 28.22326158, 28.75577595, 28.57827116, 29.99830948, 30.53082385, 30.70832864, 30.70832864, 32.48337654, 34.61343402, 36.21097713, 35.50095797, 37.80852023, 37.98602502, 37.63101544, 35.32345318, 35.85596755, 37.27600586, 36.38848191, 34.79093881, 39.22855855, 44.37619746, 44.19869267, 45.97374057, 45.79623578, 43.84368309, 44.90871183, 47.21627409, 45.97374057, 49.34633157, 49.70134115, 48.10379804, 50.05635073, 46.86126451, 47.0387693, 46.86126451, 45.79623578, 51.83139863, 52.00890342, 52.18640821, 53.25143695, 56.80153275, 59.02034262, 55.02648485, 58.13281867, 55.9140088, 59.02034262, 56.35777077, 50.76636989, 44.37619746, 50.94387468, 50.76636989, 50.23385552, 51.65389384, 48.45880762, 52.89642737, 49.52383636, 49.70134115, 51.12137947, 51.47638905, 51.29888426, 52.89642737, 54.1389609, 55.02648485, 55.47024682, 52.89642737, 50.23385552, 52.98517976, 55.02648485, 55.9140088, 54.84898006, 60.88414291, 61.77166686, 61.77166686, 61.23915249, 61.94917165, 64.43423871, 65.32176266, 61.41665728, 61.41665728, 63.72421955, 63.19170518, 63.72421955, 63.36920997, 60.35162854, 59.46410459, 57.68905669, 59.81911417, 56.44652317, 53.78395132, 53.42894174, 59.64160938, 57.33404711, 54.84898006, 48.8138172, 45.4412262, 50.05635073, 51.12137947, 52.71892258, 57.33404711, 49.70134115, 51.12137947, 51.56514145, 47.0387693, 48.72506481, 48.54756002, 46.86126451, 47.1275217, 49.96759834, 53.25143695, 51.12137947, 51.65389384, 52.54141779, 53.42894174, 52.54141779, 54.84898006, 53.16268455, 53.07393216, 47.83754086, 48.28130283, 50.41136031, 50.05635073, 49.70134115, 51.12137947, 51.12137947, 49.16882678, 49.34633157, 51.29888426, 50.94387468, 49.43508397, 53.25143695, 51.83139863, 55.20398964, 53.25143695, 55.38149443, 54.84898006, 54.84898006, 55.38149443, 52.363913, 51.12137947, 50.76636989, 52.363913, 50.76636989, 48.10379804, 49.08007439, 48.9025696, 53.78395132, 53.78395132, 55.73650401, 54.1389609, 53.60644653, 50.85512229, 50.23385552, 50.85512229, 52.45266539, 50.94387468, 48.45880762, 49.96759834, 50.50011271, 49.52383636, 51.56514145, 53.78395132, 54.31646569, 55.02648485, 53.25143695, 52.363913, 52.89642737, 54.31646569, 52.98517976, 52.18640821, 50.32260792, 51.56514145, 52.80767497, 50.94387468, 50.94387468, 50.50011271, 50.94387468, 52.54141779, 51.21013187, 52.00890342, 48.8138172, 49.70134115, 50.05635073, 51.12137947, 51.12137947, 49.16882678, 47.48253128, 48.63631241, 47.21627409, 48.54756002, 50.5888651, 48.63631241, 49.70134115, 52.71892258, 52.09765581, 51.29888426, 51.83139863, 54.1389609, 52.54141779, 55.73650401, 57.68905669, 58.57658064, 60.70663812, 61.23915249, 61.59416207, 60.35162854, 58.57658064, 58.57658064, 57.86656148, 56.44652317, 55.9140088, 55.20398964, 56.26901838, 57.5115519, 56.80153275, 56.62402796, 57.86656148, 57.68905669, 58.04406627, 59.99661896, 57.86656148, 59.2865998, 59.46410459, 61.77166686, 59.81911417, 61.23915249, 61.23915249, 58.75408543, 51.47638905, 50.94387468, 49.96759834, 49.43508397, 48.63631241, 48.9025696, 48.99132199, 48.45880762],
							[17.3954694, 17.57297419, 17.75047898, 18.99301251, 18.63800293, 17.92798377, 18.28299335, 18.63800293, 19.88053646, 20.76806041, 20.23554604, 20.41305083, 20.76806041, 22.89811789, 22.89811789, 22.7206131, 22.36560352, 23.60813705, 23.78564184, 23.60813705, 24.31815621, 22.54310831, 23.60813705, 23.43063226, 22.89811789, 22.18809873, 22.18809873, 22.01059394, 21.83308915, 22.01059394, 22.7206131, 22.18809873, 22.01059394, 23.78564184, 24.14065142, 23.60813705, 24.67316579, 24.85067058, 22.89811789, 24.31815621, 22.18809873, 21.83308915, 21.83308915, 23.96314663, 24.31815621, 23.78564184, 25.56068974, 25.02817537, 24.495661, 24.495661, 25.02817537, 26.44821368, 28.04575679, 29.28829032, 28.22326158, 29.11078553, 30.17581427, 31.06333822, 30.35331906, 33.72591007, 33.90341486, 33.37090049, 32.66088133, 32.83838612, 32.83838612, 32.66088133, 30.70832864, 31.06333822, 30.70832864, 29.99830948, 30.53082385, 27.868252, 28.04575679, 27.69074721, 27.51324242, 28.75577595, 29.82080469, 29.99830948, 28.93328074, 28.75577595, 30.35331906, 29.6432999, 30.53082385, 30.88583343, 32.48337654, 34.43592923, 34.43592923, 37.27600586, 35.14594839, 37.63101544, 36.74349149, 36.74349149, 35.14594839, 39.76107292, 39.76107292, 42.95615914, 45.26372141, 41.53612082, 43.13366393, 44.37619746, 44.37619746, 46.15124536, 43.6661783, 42.06863519, 39.05105376, 36.74349149, 39.05105376, 36.5659867, 39.58356813, 38.51853939, 41.18111124, 42.24613998, 38.51853939, 38.16352981, 36.92099628, 38.51853939, 39.93857771, 39.93857771, 39.22855855, 38.87354897, 39.58356813, 38.51853939, 36.74349149, 35.67846276, 39.76107292, 41.26986364, 42.60114956, 42.24613998, 42.68990195, 42.33489237, 42.33489237, 42.60114956, 43.22241632, 43.48867351, 43.13366393, 42.42364477, 41.53612082, 41.97988279, 40.02733011, 40.1160825, 42.24613998, 40.82610166, 40.29358729, 36.92099628, 31.95086217, 31.50710019, 33.99216725, 33.63715767, 33.1046433, 31.15209062, 26.27070889, 28.84452835, 30.08706188, 31.95086217, 34.08091965, 33.45965288, 33.63715767, 33.28214809, 31.68460498, 34.08091965, 32.48337654, 33.90341486, 34.34717683, 34.52468162, 34.8796912, 33.63715767, 33.72591007, 38.69604418, 34.9684436, 36.12222473, 35.94471994, 36.12222473, 34.61343402, 34.16967204, 34.52468162, 36.03347234, 38.16352981, 36.74349149, 36.47723431, 37.36475826, 36.21097713, 35.94471994, 37.71976784, 39.31731095, 40.64859687, 40.82610166, 39.40606334, 41.44736843, 41.8911304, 42.33489237, 42.15738758, 41.71362561, 43.48867351, 40.1160825, 41.09235885, 42.42364477, 43.39992111, 44.90871183, 44.73120704, 45.26372141, 46.95001691, 48.19255044, 49.25757918, 51.65389384, 52.45266539, 55.73650401, 53.25143695, 57.68905669, 55.55899922, 54.31646569, 54.67147527, 50.41136031, 52.00890342, 53.25143695, 54.49397048, 52.89642737, 53.96145611, 54.31646569, 52.09765581, 50.50011271, 53.16268455, 54.1389609, 56.44652317, 55.9140088, 53.96145611, 54.31646569, 55.55899922, 55.9140088, 55.20398964, 56.80153275, 55.9140088, 57.33404711, 59.99661896, 60.52913333, 61.77166686, 60.35162854, 59.10909501, 59.64160938, 61.41665728, 62.12667644, 63.01420039, 64.96675308, 66.91930577, 66.20928661, 66.56429619, 67.45182014, 68.51684887, 70.11439198, 71.5344303, 70.46940156, 71.00191593, 71.88943988, 77.39208837, 75.61704047, 70.11439198, 72.06694467, 68.16183929, 67.80682971, 70.82441114, 66.3867914, 68.51684887, 70.11439198, 71.35692551, 67.62932492, 67.45182014, 68.16183929, 67.09681056, 67.80682971, 68.69435366, 73.84199257, 75.0845261, 78.81212668, 75.97205005, 79.16713626, 78.63462189, 84.66978475, 80.58717458, 84.13727038, 85.91231828, 84.13727038, 84.66978475, 83.42725122, 84.31477517, 86.62233744, 90.97120479, 85.20229912, 80.05466021, 73.3094782, 76.14955484, 78.4571171, 78.63462189],
							[24.61991435, 24.85067058, 26.39496225, 26.62571847, 25.73819453, 25.06367632, 24.17615237, 24.61991435, 26.62571847, 27.2824862, 28.84452835, 31.06333822, 32.60762989, 35.05719599, 33.72591007, 31.72010594, 33.05139187, 35.50095797, 32.60762989, 32.83838612, 33.72591007, 31.50710019, 31.50710019, 33.28214809, 33.49515384, 29.96280852, 32.1816184, 32.39462414, 29.96280852, 31.06333822, 32.60762989, 30.83258199, 31.50710019, 33.28214809, 34.82643976, 33.72591007, 33.93891582, 35.27020174, 30.38882002, 32.1816184, 28.17001015, 26.83872422, 26.83872422, 27.2824862, 27.06948045, 27.51324242, 29.7320523, 29.7320523, 28.17001015, 27.51324242, 27.74399865, 30.17581427, 30.61957625, 31.72010594, 29.28829032, 27.51324242, 29.7320523, 29.0575341, 28.61377212, 30.17581427, 29.7320523, 31.50710019, 32.1816184, 34.82643976, 35.50095797, 35.71396371, 34.61343402, 35.71396371, 34.61343402, 33.49515384, 33.93891582, 32.83838612, 31.50710019, 29.96280852, 31.95086217, 31.95086217, 32.83838612, 33.05139187, 31.95086217, 32.39462414, 32.83838612, 31.50710019, 31.95086217, 31.72010594, 34.82643976, 36.38848191, 38.16352981, 44.14544123, 44.37619746, 48.37005523, 46.59500733, 46.15124536, 48.8138172, 55.9140088, 52.80767497, 55.47024682, 62.57043841, 55.9140088, 57.24529472, 60.79539052, 62.12667644, 66.12053421, 63.01420039, 59.46410459, 51.92015103, 50.5888651, 52.363913, 46.59500733, 46.15124536, 43.27566776, 43.93243548, 44.81995943, 41.26986364, 40.38233969, 39.93857771, 45.26372141, 46.59500733, 47.0387693, 44.81995943, 41.48286938, 40.82610166, 38.16352981, 37.71976784, 37.89727263, 41.09235885, 39.49481574, 41.8911304, 41.97988279, 42.06863519, 42.06863519, 41.18111124, 42.24613998, 43.75493069, 42.60114956, 41.00360645, 39.76107292, 37.48901161, 39.05105376, 37.48901161, 37.71976784, 37.93277359, 37.04524964, 37.04524964, 32.60762989, 28.17001015, 26.83872422, 28.40076637, 27.74399865, 27.51324242, 25.06367632, 24.85067058, 25.73819453, 28.84452835, 27.9570044, 26.83872422, 28.40076637, 27.9570044, 27.74399865, 25.95120027, 27.15823284, 27.77949961, 28.13450919, 27.51324242, 26.71447087, 26.00445171, 24.67316579, 25.11692776, 29.7320523, 31.24084301, 31.15209062, 30.26456667, 29.28829032, 29.46579511, 29.02203314, 29.28829032, 31.15209062, 33.45965288, 31.95086217, 32.83838612, 31.77335738, 28.93328074, 28.57827116, 28.75577595, 28.40076637, 28.84452835, 33.28214809, 33.01589091, 33.90341486, 37.45351065, 37.36475826, 36.6547391, 36.74349149, 37.63101544, 35.94471994, 34.61343402, 35.14594839, 36.6547391, 35.41220557, 36.83224389, 36.21097713, 35.50095797, 39.40606334, 38.87354897, 39.22855855, 41.97988279, 42.95615914, 38.3410346, 41.71362561, 41.35861603, 38.51853939, 40.29358729, 39.05105376, 42.42364477, 40.55984448, 40.47109208, 38.78479658, 39.76107292, 38.96230137, 38.78479658, 37.54226305, 37.54226305, 39.05105376, 42.06863519, 41.35861603, 38.78479658, 38.78479658, 40.2048349, 39.84982532, 39.49481574, 38.96230137, 38.429787, 39.93857771, 39.67232053, 40.38233969, 41.44736843, 39.76107292, 39.22855855, 39.40606334, 40.91485406, 40.29358729, 42.06863519, 43.13366393, 43.84368309, 43.75493069, 44.99746422, 51.29888426, 51.29888426, 55.02648485, 56.97903754, 55.55899922, 54.84898006, 55.73650401, 58.22157106, 56.62402796, 58.22157106, 63.01420039, 59.64160938, 55.73650401, 58.22157106, 55.73650401, 53.25143695, 51.56514145, 52.45266539, 50.32260792, 46.06249296, 45.4412262, 42.60114956, 41.97988279, 41.71362561, 41.00360645, 43.13366393, 47.48253128, 46.06249296, 47.92629325, 48.63631241, 47.48253128, 42.51239716, 45.88498817, 45.3524738, 43.31116872, 44.02118788, 42.06863519, 42.33489237, 48.45880762, 52.71892258, 47.30502649, 43.31116872, 40.47109208, 42.24613998, 41.71362561, 41.00360645],
							[43.70167926, 44.78445847, 47.78428942, 48.88481912, 49.96759834, 50.251606, 50.51786319, 57.68905669, 58.0263158, 59.05584358, 59.05584358, 59.05584358, 60.75988956, 62.12667644, 62.46393554, 60.75988956, 62.46393554, 68.26834217, 71.67643413, 68.26834217, 71.67643413, 71.00191593, 73.04322101, 75.77679478, 73.7354897, 70.32739773, 69.63512905, 71.00191593, 68.96061085, 69.63512905, 68.96061085, 68.26834217, 67.23881439, 72.36870281, 74.4100079, 73.04322101, 76.46906346, 79.20263722, 73.04322101, 75.10227658, 64.8602502, 65.1975093, 67.23881439, 71.00191593, 72.36870281, 73.04322101, 77.83585034, 79.20263722, 77.83585034, 77.83585034, 79.20263722, 81.2439423, 86.01882115, 91.48596868, 90.1191818, 88.07787671, 94.21954244, 99.4026823, 96.74011046, 103.84030205, 103.84030205, 101.1777302, 101.1777302, 102.06525415, 100.29020625, 104.727826, 103.84030205, 103.84030205, 104.727826, 101.1777302, 103.84030205, 102.06525415, 102.9527781, 100.29020625, 103.84030205, 109.16544574, 107.39039785, 102.06525415, 97.62763441, 95.85258651, 95.85258651, 93.19001466, 94.96506256, 94.96506256, 97.62763441, 102.9527781, 101.1777302, 104.727826, 101.1777302, 105.61534995, 103.84030205, 106.5028739, 106.5028739, 114.49058944, 129.57849658, 128.69097263, 134.01611632, 126.91592473, 126.91592473, 131.35354447, 131.35354447, 142.00383186, 137.56621212, 129.57849658, 119.81573313, 102.06525415, 101.1777302, 93.19001466, 93.19001466, 89.63991886, 91.41496676, 93.19001466, 92.30249071, 91.41496676, 88.30863294, 91.41496676, 100.29020625, 98.51515836, 95.85258651, 93.19001466, 96.74011046, 89.63991886, 91.85872874, 89.63991886, 96.74011046, 95.85258651, 97.18387243, 97.18387243, 96.74011046, 93.63377663, 93.19001466, 100.73396823, 99.4026823, 98.95892033, 96.74011046, 96.29634848, 94.07753861, 94.96506256, 98.51515836, 98.07139638, 101.62149218, 102.50901613, 105.17158797, 101.62149218, 100.29020625, 94.07753861, 98.07139638, 100.29020625, 98.51515836, 93.19001466, 84.84728954, 83.24974643, 75.97205005, 75.97205005, 81.65220332, 81.47469853, 82.36222248, 85.5573087, 83.60475601, 91.41496676, 94.07753861, 96.74011046, 90.52744281, 93.19001466, 94.07753861, 90.52744281, 90.52744281, 100.73396823, 101.1777302, 104.28406402, 107.83415982, 103.84030205, 103.39654008, 103.39654008, 104.727826, 110.05296969, 114.49058944, 111.82801759, 113.60306549, 114.49058944, 114.04682747, 108.2779218, 110.94049364, 111.82801759, 113.60306549, 115.37811339, 111.82801759, 113.15930352, 113.60306549, 113.15930352, 115.37811339, 113.60306549, 116.70939931, 113.60306549, 111.38425562, 114.93435141, 120.70325708, 122.03454301, 122.03454301, 122.03454301, 122.92206696, 124.25335288, 123.36582893, 126.91592473, 130.9097825, 130.9097825, 131.35354447, 133.12859237, 130.9097825, 137.56621212, 135.34740224, 133.12859237, 140.67254594, 140.22878397, 138.89749804, 140.67254594, 142.00383186, 135.34740224, 135.79116422, 134.4598783, 139.34126002, 140.22878397, 143.33511779, 141.11630791, 138.00997409, 138.45373607, 138.45373607, 138.45373607, 134.90364027, 136.23492619, 136.67868817, 139.34126002, 141.11630791, 140.67254594, 144.22264174, 140.67254594, 138.45373607, 143.33511779, 140.67254594, 143.33511779, 148.66026148, 155.76045308, 165.07945454, 161.97312072, 164.19193059, 168.18578836, 168.18578836, 161.52935874, 165.52321652, 160.19807282, 159.75431085, 163.30440664, 164.19193059, 164.19193059, 161.08559677, 169.51707429, 161.97312072, 161.97312072, 165.52321652, 157.091739, 161.97312072, 158.8667869, 149.54778543, 145.99768963, 145.11016569, 142.00383186, 140.67254594, 144.66640371, 147.32897556, 152.21035728, 156.64797702, 166.41074046, 164.19193059, 160.19807282, 165.07945454, 158.8667869, 157.091739, 194.36774486, 196.14279276, 196.14279276, 197.91784066, 197.03031671, 202.35546041, 199.69288856, 182.82993352, 185.49250537, 181.94240957, 165.96697849, 169.51707429, 166.85450244, 160.6418348],
							[44.21644315, 43.73718021, 44.69570608, 43.89693452, 43.89693452, 43.08041249, 43.2401668, 45.17496901, 46.7902626, 50.42911079, 48.8138172, 49.61258876, 51.22788234, 54.0502085, 53.65969797, 51.63614336, 53.25143695, 57.28079568, 54.86673054, 54.0502085, 54.0502085, 48.40555619, 50.02084977, 50.02084977, 49.22207822, 47.60678463, 47.76653894, 50.02084977, 48.08604756, 49.22207822, 49.22207822, 50.42911079, 52.04440438, 54.45846952, 57.68905669, 55.27499155, 56.07376311, 54.86673054, 50.02084977, 50.42911079, 44.69570608, 42.12188663, 42.60114956, 44.69570608, 45.0152147, 44.69570608, 46.31099967, 46.95001691, 44.53595177, 44.05668884, 46.63050829, 50.83737181, 53.25143695, 54.45846952, 49.61258876, 49.22207822, 53.25143695, 52.45266539, 51.63614336, 56.89028514, 56.89028514, 58.48782825, 60.51138285, 62.53493746, 62.92544799, 63.74197003, 61.71841542, 63.33370901, 62.53493746, 61.32790489, 63.74197003, 61.71841542, 59.7126113, 57.28079568, 60.10312184, 63.33370901, 62.12667644, 62.92544799, 61.32790489, 60.51138285, 62.92544799, 61.71841542, 62.53493746, 65.76552463, 68.58785079, 75.84779669, 77.85360082, 85.52180774, 85.52180774, 96.01234082, 90.3676885, 87.86487097, 88.75239491, 98.51515836, 103.84030205, 118.04068524, 126.02840078, 123.36582893, 115.37811339, 117.15316129, 115.37811339, 122.47830498, 112.71554154, 107.39039785, 101.1777302, 104.727826, 104.727826, 94.07753861, 98.51515836, 94.07753861, 98.51515836, 95.85258651, 82.98348925, 82.53972727, 82.0959653, 86.97734702, 94.96506256, 90.52744281, 89.63991886, 88.75239491, 86.08982307, 82.0959653, 76.50456442, 75.79454526, 85.37980391, 84.49227996, 89.19615689, 90.97120479, 90.97120479, 92.74625269, 94.52130058, 100.73396823, 106.05911192, 106.5028739, 102.50901613, 102.06525415, 98.51515836, 101.62149218, 99.84644428, 101.62149218, 102.50901613, 101.62149218, 100.29020625, 88.75239491, 85.5573087, 81.65220332, 83.07224164, 83.42725122, 82.0072129, 73.3094782, 68.51684887, 70.82441114, 74.37450694, 76.859574, 83.95976559, 84.13727038, 89.19615689, 91.41496676, 87.3323566, 96.29634848, 92.74625269, 93.63377663, 92.30249071, 94.07753861, 91.85872874, 87.68736618, 88.21988055, 99.4026823, 94.52130058, 98.07139638, 102.50901613, 100.29020625, 102.50901613, 101.62149218, 101.62149218, 105.17158797, 106.05911192, 102.9527781, 102.9527781, 101.1777302, 94.52130058, 93.63377663, 96.29634848, 102.06525415, 102.9527781, 110.94049364, 109.16544574, 111.38425562, 111.82801759, 111.38425562, 108.2779218, 107.39039785, 110.05296969, 103.84030205, 101.62149218, 106.5028739, 110.94049364, 108.2779218, 106.94663587, 109.16544574, 112.27177957, 118.48444721, 114.93435141, 114.93435141, 124.25335288, 131.35354447, 125.5846388, 129.57849658, 126.91592473, 124.25335288, 123.36582893, 117.15316129, 125.14087683, 122.92206696, 131.79730645, 128.69097263, 128.24721065, 126.02840078, 128.24721065, 128.69097263, 132.24106842, 138.00997409, 140.22878397, 140.67254594, 136.67868817, 138.00997409, 138.89749804, 138.89749804, 132.24106842, 133.57235435, 129.1347346, 131.79730645, 134.01611632, 134.90364027, 138.45373607, 134.01611632, 132.24106842, 137.12245014, 140.22878397, 143.33511779, 146.44145161, 153.98540518, 153.98540518, 152.65411925, 154.42916715, 157.091739, 160.6418348, 167.74202639, 170.40459824, 171.29212219, 161.08559677, 166.41074046, 169.96083626, 167.29826441, 161.97312072, 166.41074046, 159.75431085, 155.3166911, 165.07945454, 157.97926295, 160.6418348, 158.8667869, 157.53550097, 148.21649951, 145.99768963, 145.55392766, 138.45373607, 138.00997409, 138.89749804, 143.77887976, 155.3166911, 162.86064467, 166.41074046, 168.18578836, 169.07331231, 164.63569257, 152.21035728, 169.07331231, 165.52321652, 153.09788123, 164.19193059, 157.97926295, 161.97312072, 179.27983773, 185.49250537, 170.84836021, 166.41074046, 156.20421505, 156.64797702, 158.8667869, 163.74816862],
							[41.1278598, 41.41186747, 43.18691537, 45.84948721, 45.49447763, 45.84948721, 46.59500733, 48.44105714, 50.65986702, 49.55933732, 49.55933732, 50.28710696, 51.03262708, 51.03262708, 52.15090725, 51.03262708, 49.93209738, 51.03262708, 51.03262708, 49.93209738, 50.65986702, 46.96776739, 47.34052745, 48.06829709, 48.06829709, 46.96776739, 47.34052745, 46.59500733, 45.49447763, 46.59500733, 47.34052745, 47.34052745, 49.55933732, 53.25143695, 53.25143695, 53.99695707, 52.15090725, 53.25143695, 48.44105714, 49.93209738, 44.74895752, 44.74895752, 45.49447763, 47.71328751, 48.8138172, 48.8138172, 50.65986702, 50.65986702, 51.40538713, 50.28710696, 51.40538713, 51.40538713, 54.7247267, 54.7247267, 52.15090725, 52.15090725, 54.35196665, 53.25143695, 52.15090725, 53.62419701, 54.35196665, 53.99695707, 53.25143695, 55.84300688, 56.94353658, 58.78958639, 58.06181675, 58.43457681, 58.06181675, 58.43457681, 59.90786657, 58.43457681, 57.68905669, 55.47024682, 58.06181675, 61.75391638, 62.12667644, 63.97272625, 62.87219656, 63.97272625, 64.34548631, 62.12667644, 64.71824637, 66.56429619, 69.88363576, 74.69401556, 73.96624592, 76.91282543, 73.2207258, 80.62267554, 79.87715542, 82.84148541, 82.0959653, 92.30249071, 91.41496676, 93.19001466, 100.29020625, 89.63991886, 88.30863294, 87.42110899, 85.64606109, 86.53358504, 82.53972727, 78.98963147, 72.77696383, 76.32705963, 75.88329765, 70.55815396, 71.44567791, 70.55815396, 72.33320186, 78.10210752, 77.65834555, 71.44567791, 71.00191593, 73.2207258, 75.88329765, 74.55201173, 73.66448778, 73.2207258, 74.55201173, 72.77696383, 68.69435366, 66.56429619, 69.40437282, 70.46940156, 72.06694467, 70.29189677, 68.87185845, 69.04936324, 69.40437282, 72.59945904, 74.37450694, 72.06694467, 70.11439198, 69.7593824, 69.40437282, 71.5344303, 69.7593824, 70.11439198, 69.40437282, 67.9843345, 68.16183929, 63.19170518, 56.62402796, 55.02648485, 58.93159022, 60.88414291, 58.57658064, 53.96145611, 55.20398964, 55.73650401, 57.68905669, 59.81911417, 66.20928661, 67.09681056, 66.56429619, 64.07922913, 63.72421955, 66.56429619, 67.09681056, 66.56429619, 64.96675308, 66.20928661, 64.96675308, 64.78924829, 64.78924829, 69.22686803, 70.11439198, 73.13197341, 75.43953568, 72.95446862, 72.59945904, 71.17942072, 70.64690635, 74.37450694, 71.88943988, 69.93688719, 72.24444946, 70.46940156, 66.74180098, 67.27431535, 66.56429619, 67.62932492, 69.7593824, 72.42195425, 70.82441114, 71.00191593, 72.77696383, 72.77696383, 72.42195425, 72.42195425, 73.48698299, 65.14425787, 62.48168602, 64.25673392, 65.14425787, 63.54671476, 62.8366956, 61.94917165, 63.19170518, 64.25673392, 64.78924829, 67.45182014, 66.20928661, 67.45182014, 66.91930577, 67.09681056, 66.20928661, 64.43423871, 64.43423871, 60.17412375, 62.8366956, 61.94917165, 63.19170518, 61.23915249, 62.30418123, 63.54671476, 63.72421955, 65.14425787, 63.36920997, 63.90172434, 65.32176266, 64.96675308, 63.01420039, 62.65919081, 62.30418123, 61.41665728, 59.64160938, 59.46410459, 58.39907585, 57.5115519, 57.33404711, 56.97903754, 59.46410459, 58.39907585, 57.5115519, 60.17412375, 63.54671476, 63.19170518, 63.01420039, 65.49926745, 66.03178182, 63.54671476, 64.25673392, 61.0616477, 59.2865998, 60.35162854, 58.93159022, 59.46410459, 59.99661896, 59.46410459, 59.2865998, 60.35162854, 59.81911417, 60.35162854, 58.75408543, 64.43423871, 64.78924829, 62.12667644, 61.41665728, 63.54671476, 63.19170518, 62.48168602, 61.59416207, 60.52913333, 59.10909501, 62.65919081, 61.23915249, 61.94917165, 61.23915249, 63.90172434, 65.67677224, 67.80682971, 68.87185845, 79.34464105, 75.0845261, 79.34464105, 77.92460274, 75.43953568, 75.97205005, 78.81212668, 78.98963147, 78.4571171, 78.10210752, 74.55201173, 68.33934408, 63.72421955, 69.58187761, 74.55201173, 73.3094782],
							[2.43181562, 2.34306323, 2.30756227, 3.81635298, 3.67434915, 4.08261017, 4.11811112, 4.4731207, 4.4731207, 4.43761975, 4.61512454, 4.43761975, 4.29561591, 4.34886735, 4.38436831, 4.11811112, 4.91688268, 5.11213795, 5.85765806, 6.03516285, 6.07066381, 5.36064465, 5.62690184, 5.89315902, 6.39017243, 6.336921, 6.15941621, 6.07066381, 6.12391525, 5.89315902, 5.68015327, 5.50264848, 5.76890567, 5.9996619, 6.336921, 6.336921, 6.95818776, 6.65642962, 6.07066381, 6.21266764, 5.36064465, 5.11213795, 5.14763891, 5.36064465, 5.41389609, 5.71565423, 5.85765806, 5.94641046, 5.59140088, 5.55589992, 5.62690184, 5.80440663, 5.80440663, 6.44342387, 6.47892483, 6.07066381, 5.9996619, 6.51442579, 6.21266764, 6.39017243, 6.56767722, 6.56767722, 6.336921, 6.30142004, 6.39017243, 6.39017243, 6.15941621, 6.15941621, 6.21266764, 5.89315902, 5.89315902, 5.62690184, 5.50264848, 5.55589992, 5.44939705, 5.32514369, 5.32514369, 5.27189226, 5.18313986, 5.2363913, 5.55589992, 5.14763891, 5.36064465, 5.18313986, 5.32514369, 5.80440663, 6.56767722, 6.78068297, 6.44342387, 6.95818776, 6.78068297, 6.56767722, 6.74518201, 7.04694016, 8.37822608, 7.98771554, 7.89896315, 8.16522033, 7.33094782, 7.36644878, 7.40194974, 7.36644878, 7.18894399, 7.10019159, 6.30142004, 6.12391525, 6.39017243, 5.80440663, 6.21266764, 5.62690184, 5.62690184, 5.62690184, 5.2363913, 5.18313986, 5.00563507, 5.05888651, 5.62690184, 5.62690184, 5.85765806, 5.50264848, 5.32514369, 5.32514369, 5.27189226, 5.25414178, 5.60915136, 5.30739322, 5.46714753, 5.9109095, 5.83990759, 5.64465232, 5.71565423, 5.76890567, 6.15941621, 5.89315902, 5.80440663, 5.37839513, 5.27189226, 5.11213795, 4.97013412, 5.32514369, 5.2363913, 5.05888651, 5.00563507, 4.5618731, 3.81635298, 3.53234532, 3.72760059, 3.47909388, 3.5500958, 3.24833765, 2.59156993, 2.92882903, 2.82232616, 2.71582328, 3.24833765, 3.1240843, 3.1240843, 2.98208047, 2.89332807, 3.03533191, 2.96432999, 3.4613434, 3.30158909, 3.37259101, 3.44359292, 3.44359292, 3.58559675, 4.4731207, 4.43761975, 4.82813028, 4.88138172, 4.70387693, 4.52637214, 4.5618731, 4.52637214, 4.5618731, 5.00563507, 4.82813028, 4.91688268, 5.41389609, 5.41389609, 5.55589992, 5.50264848, 5.62690184, 6.15941621, 6.30142004, 6.2481686, 5.89315902, 5.89315902, 6.12391525, 6.21266764, 6.07066381, 6.07066381, 5.85765806, 5.62690184, 6.03516285, 6.07066381, 6.336921, 7.04694016, 6.9226868, 7.10019159, 8.1119689, 8.43147752, 7.98771554, 7.66820692, 8.1119689, 7.93446411, 8.1119689, 7.93446411, 7.84571171, 7.81021075, 7.36644878, 8.21847177, 7.89896315, 7.72145836, 7.93446411, 8.07646794, 8.52022991, 8.55573087, 8.60898231, 9.23024907, 9.76276344, 10.02902063, 10.4727826, 9.49650626, 9.31900147, 9.67401105, 10.73903978, 11.09404936, 11.36030655, 11.00529697, 11.27155415, 12.07032571, 12.1590781, 12.78034487, 12.33658289, 11.98157331, 12.51408768, 12.95784966, 13.40161163, 15.17665953, 15.08790714, 14.73289756, 14.99915474, 15.79792629, 16.95170743, 16.41919306, 17.57297419, 19.70303167, 23.34187986, 24.85067058, 24.05189902, 27.24698524, 31.50710019, 29.02203314, 30.35331906, 28.04575679, 29.37704272, 29.46579511, 28.22326158, 28.93328074, 30.88583343, 31.24084301, 30.70832864, 28.48951877, 29.6432999, 28.66702356, 31.86210977, 34.8796912, 41.18111124, 38.78479658, 39.84982532, 38.3410346, 46.50625494, 48.28130283, 47.57128367, 46.50625494, 54.1389609, 67.45182014, 63.72421955, 61.41665728, 63.54671476, 63.72421955, 67.62932492, 67.62932492, 66.56429619, 63.01420039, 57.68905669, 54.84898006, 59.99661896, 56.80153275],
							[29.11078553, 27.868252, 29.28829032, 28.40076637, 28.57827116, 28.93328074, 28.75577595, 29.82080469, 31.06333822, 33.1933957, 32.66088133, 33.1933957, 33.72591007, 34.08091965, 34.9684436, 33.72591007, 33.90341486, 34.9684436, 33.90341486, 33.54840528, 33.72591007, 31.77335738, 30.35331906, 29.82080469, 30.35331906, 28.75577595, 29.11078553, 28.93328074, 27.868252, 28.22326158, 28.04575679, 28.04575679, 27.868252, 30.53082385, 32.30587175, 32.30587175, 32.83838612, 32.66088133, 29.82080469, 29.99830948, 25.02817537, 23.96314663, 23.43063226, 24.85067058, 25.73819453, 25.56068974, 27.33573763, 28.04575679, 27.51324242, 28.04575679, 28.75577595, 31.06333822, 32.48337654, 32.83838612, 31.77335738, 30.17581427, 32.48337654, 32.66088133, 32.12836696, 33.72591007, 34.08091965, 34.61343402, 34.43592923, 38.69604418, 39.40606334, 41.8911304, 39.05105376, 38.16352981, 38.3410346, 38.16352981, 39.05105376, 37.98602502, 36.38848191, 34.79093881, 37.27600586, 39.22855855, 37.45351065, 38.16352981, 36.38848191, 35.85596755, 36.74349149, 36.5659867, 37.63101544, 42.06863519, 42.95615914, 51.47638905, 53.69519892, 56.35777077, 52.363913, 60.35162854, 61.23915249, 60.35162854, 63.01420039, 69.67063001, 72.33320186, 78.10210752, 101.1777302, 89.63991886, 88.30863294, 96.74011046, 94.07753861, 97.62763441, 94.96506256, 84.75853714, 79.87715542, 83.87101319, 84.31477517, 82.53972727, 79.43339345, 74.55201173, 78.10210752, 77.65834555, 67.89558211, 67.45182014, 65.67677224, 71.88943988, 75.88329765, 71.88943988, 71.00191593, 71.00191593, 68.33934408, 65.67677224, 62.30418123, 59.81911417, 69.58187761, 68.33934408, 72.42195425, 71.5344303, 70.46940156, 73.3094782, 73.48698299, 79.34464105, 86.62233744, 87.15485181, 85.02479433, 85.20229912, 83.95976559, 87.68736618, 86.26732786, 87.15485181, 90.52744281, 85.91231828, 88.57489012, 76.50456442, 71.00191593, 64.25673392, 63.72421955, 63.72421955, 65.49926745, 60.70663812, 59.64160938, 58.04406627, 59.64160938, 61.23915249, 66.74180098, 65.67677224, 73.3094782, 75.97205005, 73.66448778, 79.69965063, 76.68206921, 77.21458358, 77.03707879, 79.16713626, 77.56959316, 70.82441114, 69.93688719, 79.87715542, 75.97205005, 78.10210752, 78.63462189, 77.39208837, 77.56959316, 76.14955484, 75.26203089, 78.98963147, 81.82970811, 78.98963147, 77.92460274, 76.14955484, 73.3094782, 72.24444946, 73.84199257, 74.90702131, 76.32705963, 81.29719374, 80.58717458, 82.89473685, 86.44483265, 86.79984223, 85.73481349, 85.37980391, 82.89473685, 79.87715542, 78.10210752, 81.29719374, 84.84728954, 81.82970811, 81.82970811, 82.36222248, 82.71723206, 87.15485181, 90.08368084, 92.74625269, 98.95892033, 104.28406402, 102.06525415, 103.39654008, 102.06525415, 99.84644428, 100.29020625, 93.63377663, 98.51515836, 96.29634848, 98.07139638, 97.62763441, 97.18387243, 94.96506256, 95.85258651, 92.30249071, 96.29634848, 98.51515836, 103.84030205, 101.62149218, 98.07139638, 98.07139638, 102.9527781, 105.17158797, 100.73396823, 101.1777302, 97.62763441, 101.62149218, 103.39654008, 107.39039785, 108.72168377, 107.39039785, 106.05911192, 112.71554154, 114.49058944, 115.37811339, 118.92820919, 122.03454301, 124.69711486, 126.47216275, 124.25335288, 125.5846388, 126.47216275, 135.34740224, 138.00997409, 134.90364027, 127.80344868, 131.35354447, 135.34740224, 132.6848304, 128.24721065, 134.4598783, 127.3596867, 127.3596867, 130.46602052, 122.47830498, 126.91592473, 122.92206696, 122.03454301, 111.82801759, 113.60306549, 114.04682747, 109.60920772, 114.04682747, 114.93435141, 114.93435141, 119.81573313, 126.02840078, 127.80344868, 133.57235435, 134.01611632, 127.80344868, 116.26563734, 128.24721065, 122.03454301, 117.15316129, 123.36582893, 118.48444721, 120.25949511, 133.12859237, 133.12859237, 124.69711486, 120.70325708, 114.04682747, 118.48444721, 115.82187536, 117.15316129],
							[10.91654457, 11.27155415, 11.62656373, 11.36030655, 11.36030655, 11.18280176, 10.73903978, 11.18280176, 11.44905894, 11.98157331, 12.07032571, 13.31285924, 14.64414516, 15.08790714, 14.73289756, 14.55539277, 14.46664037, 15.53166911, 14.99915474, 14.91040235, 14.64414516, 13.57911642, 13.84537361, 14.37788798, 14.11163079, 13.49036403, 13.49036403, 13.40161163, 13.40161163, 13.75662121, 14.11163079, 13.22410684, 13.66786882, 14.55539277, 14.99915474, 15.08790714, 15.35416432, 15.53166911, 14.46664037, 14.55539277, 11.71531613, 11.53781134, 11.44905894, 12.07032571, 11.89282092, 11.98157331, 13.04660205, 13.40161163, 12.42533529, 12.60284008, 13.13535445, 13.75662121, 14.46664037, 14.91040235, 13.934126, 13.04660205, 13.31285924, 13.75662121, 13.13535445, 14.46664037, 14.99915474, 15.08790714, 14.82164995, 15.35416432, 16.06418348, 16.50794545, 15.97543108, 16.68545024, 16.95170743, 16.15293587, 16.86295503, 16.77420264, 16.24168827, 16.24168827, 16.68545024, 17.4842218, 17.75047898, 19.88053646, 20.23554604, 19.88053646, 19.88053646, 19.52552688, 20.41305083, 22.18809873, 23.07562268, 26.0932041, 26.62571847, 29.28829032, 29.6432999, 30.88583343, 28.75577595, 29.28829032, 31.24084301, 33.37090049, 36.21097713, 37.80852023, 38.69604418, 36.5659867, 35.85596755, 38.51853939, 39.76107292, 39.93857771, 36.74349149, 34.25842444, 32.83838612, 33.01589091, 32.83838612, 29.6432999, 31.77335738, 28.93328074, 31.24084301, 33.37090049, 32.66088133, 31.59585259, 30.53082385, 33.1933957, 34.08091965, 33.01589091, 32.12836696, 30.53082385, 29.82080469, 29.6432999, 28.75577595, 27.33573763, 29.55454751, 28.57827116, 30.61957625, 31.4183478, 31.3295954, 30.70832864, 29.90955709, 31.86210977, 32.21711935, 31.50710019, 30.17581427, 30.17581427, 29.11078553, 29.55454751, 29.11078553, 28.75577595, 29.02203314, 28.57827116, 28.31201398, 24.76191818, 23.34187986, 21.47807957, 22.98687028, 23.34187986, 22.80936549, 19.79178407, 18.99301251, 18.90426012, 18.72675533, 20.41305083, 22.7206131, 22.45435591, 23.07562268, 22.27685112, 21.56683196, 22.6318607, 22.45435591, 23.25312747, 22.54310831, 22.89811789, 22.7206131, 21.47807957, 21.47807957, 24.76191818, 23.96314663, 24.14065142, 24.31815621, 23.34187986, 23.07562268, 22.09934633, 22.80936549, 23.34187986, 23.69688944, 22.98687028, 22.80936549, 22.36560352, 21.38932717, 20.76806041, 20.67930802, 20.9455652, 20.9455652, 22.18809873, 20.41305083, 20.50180323, 21.47807957, 22.01059394, 22.6318607, 22.36560352, 22.98687028, 21.83308915, 21.03431759, 22.01059394, 23.16437507, 22.98687028, 23.07562268, 23.16437507, 23.51938465, 26.0932041, 27.06948045, 28.48951877, 29.11078553, 29.7320523, 29.11078553, 29.99830948, 29.55454751, 27.69074721, 29.19953793, 26.71447087, 27.51324242, 27.51324242, 27.60199482, 26.89197566, 27.42449003, 27.33573763, 26.71447087, 24.93942297, 25.64944213, 27.06948045, 27.77949961, 27.15823284, 26.35946129, 26.0932041, 27.15823284, 26.89197566, 26.44821368, 26.00445171, 25.64944213, 26.98072805, 26.00445171, 26.27070889, 27.77949961, 27.69074721, 27.24698524, 28.40076637, 29.82080469, 29.99830948, 29.90955709, 30.44207146, 31.86210977, 30.44207146, 31.50710019, 32.74963372, 35.23470078, 37.09850107, 40.29358729, 38.16352981, 37.27600586, 38.16352981, 39.31731095, 38.25228221, 37.54226305, 41.00360645, 38.96230137, 38.25228221, 38.3410346, 35.32345318, 37.89727263, 36.92099628, 35.23470078, 31.59585259, 30.17581427, 30.08706188, 28.93328074, 29.82080469, 29.02203314, 28.57827116, 28.66702356, 31.68460498, 32.57212893, 34.79093881, 33.63715767, 33.28214809, 31.4183478, 31.24084301, 33.81466246, 32.03961456, 31.77335738, 31.4183478, 31.24084301, 33.90341486, 34.25842444, 31.86210977, 31.24084301, 28.40076637, 29.28829032, 28.75577595, 28.31201398]];
		var assetsReturns = new Array(assetsPrices.length);
		for (var k = 0; k < assetsPrices.length; ++k) {
			assetsReturns[k] = new Array(assetsPrices[k].length - 1);
			for (var i = 1; i < assetsPrices[k].length; ++i) {
				assetsReturns[k][i-1] = (assetsPrices[k][i] - assetsPrices[k][i-1])/assetsPrices[k][i-1];
			}
		}
			
		// Expected weights have been computed with the method below, and value of the objective function verified to be of ~6.8e-4
		var expectedWeights = [0.011715036820412077, 0.010982741983461726, 0.014329219601669463, 0.08138764092018093, 0.005174521042167091, 0.03731032609992767, 0.026151093715953894, 0.004251777653722618, 0.006812177865239405, 0.003484984901039915, 0.09715701892213753, 0.02866896503974561, 0.036382118690755604, 0.025350781767040753, 0.17571563315373878, 0.0017205161735978484, 0.0008476736457532089, 0.025294354986590422, 0, 0.022574142290264214, 0.04012875317955057, 0.029974000042811807, 0.012634234356196011, 0.031698301448724105, 0.01762436687114238, 0.0672892678807561, 0.08493973626173573, 0.03733945213481528, 0, 0.034618365316269016, 0.028442797234600012];
		
		// Compute the weights of the minimum tracking error portfolio
		var weights = PortfolioAllocation.minimumTrackingErrorWeights(assetsReturns, indexReturns);

		//
		var weightsOK = true;
		if (expectedWeights.length != weights.length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedWeights.length; ++i) {
				if (Math.abs(weights[i] - expectedWeights[i]) > 1e-4) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, 'Minimum tracking error portfolio, indtrack1 OR data set (Hang Seng Index)');
	}
});