// ------------------------------------------------------------
QUnit.module('Assets allocation module');
	
QUnit.test('Mean variance portfolio - internal corner portfolios computation', function(assert) {    
	// Used for compatibility purposes with old JavaScript engines
	function cornerPortfoliosToArray(cornerPortfolio) {
		for (var i = 0; i < cornerPortfolio.length; ++i) {
			cornerPortfolio[i][0] = cornerPortfolio[i][0].toArray();
		}
		
		return cornerPortfolio;
	}

	// Test using static data
	// Test infeasible/unsupported cases for determining the critical line algorithm starting portfolio:
	{
		// Lower bounds > Upper bounds
		assert.throws(function() { 
			PortfolioAllocation.computeCornerPortfolios_(new PortfolioAllocation.Matrix([0.1, 0.2]), 
														 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														 { constraints: {minWeights: [0.6, 0.3], maxWeights: [0.2, 1]} }) },
			new Error('infeasible problem detected'),
			"Mean variance portfolio - Corner portfolios, lower bounds greater than upper bounds");
			
		//  Sum lower bounds > 1
		assert.throws(function() { 
			PortfolioAllocation.computeCornerPortfolios_(new PortfolioAllocation.Matrix([0.1, 0.2]), 
														 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														 { constraints: {minWeights: [0.6, 0.5], maxWeights: [0.7, 0.5]} }) },
			new Error('infeasible problem detected'),
			"Mean variance portfolio - Corner portfolios, sum of lower bounds greater than one");
		
		//  Sum upper bounds < 1
		assert.throws(function() { 
			PortfolioAllocation.computeCornerPortfolios_(new PortfolioAllocation.Matrix([0.1, 0.2]), 
														new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														{ constraints: {minWeights: [0.4, 0.4], maxWeights: [0.4, 0.5]} }) },
			new Error('infeasible problem detected'),
			"Mean variance portfolio - Corner portfolios, sum of upper bounds lower than one");
		
		// Identical returns
		assert.throws(function() { 
			PortfolioAllocation.computeCornerPortfolios_(new PortfolioAllocation.Matrix([0.1, 0.1]), 
														 new PortfolioAllocation.Matrix([[1,0],[0,1]])) },
			new Error('unsupported problem detected'),
			"Mean variance portfolio - Corner portfolios, identical returns");
	}
	
	// Test using static data
	// Reference: https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm
	// Note: this test also allows checking that numerically equal corner portfolios are filtered out by the algorithm
	{
		var expectedCornerPortfolios = [[[0.2, 0.30000000000000004, 0.5], 20.898844444444443], 
										[[0.2, 0.5, 0.30000000000000004], 11.1475], 
										[[0.22180737780348653, 0.5, 0.27819262219651353], 10.51088812347172],
										[[0.451915610952186, 0.348084389047814, 0.2], 7.55192170004087],
										[[0.5, 0.2999999999999999, 0.2], 0]];
		
		var covMat = new PortfolioAllocation.Matrix([[1, 2.96, 2.31],
													[2.96, 54.76, 39.886],
													[2.31,	39.886,	237.16]]);
		var returns = new PortfolioAllocation.Matrix([2.8000, 6.3000, 10.8000]);
		
		var cornerPortfolios = PortfolioAllocation.computeCornerPortfolios_(returns, covMat, { constraints: {minWeights: [0.2, 0.2, 0.2], maxWeights: [0.5, 0.5, 0.5]} });
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #1');
	}
		
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter VIII "The computing procedure"
	{
		var expectedCornerPortfolios = [[[0, 1, 0], 4.16666666666667], 
										[[0, 0.22496808316614988, 0.7750319168338501], 0.1408064320019454], 
										[[0.8414051841746248, 0, 0.15859481582537516], 0.03332764893133244], 
										[[0.9931034482758623, 0, 0.006896551724137813], 0]];
				
		var covMat = new PortfolioAllocation.Matrix([[0.0146, 0.0187, 0.0145],
													 [0.0187, 0.0854, 0.0104],
													  [0.0145, 0.0104, 0.0289]]);
		var returns = new PortfolioAllocation.Matrix([0.062, 0.146, 0.128]);
		
		var cornerPortfolios = PortfolioAllocation.computeCornerPortfolios_(returns, covMat);
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #2');
	}
	
	// Test using static data
	// Reference: A Simple Spreadsheet-Based Exposition of the Markowitz Critical Line Method for Portfolio Selection, Clarence C. Kwan
	{
		var expectedCornerPortfolios = [[[0, 0, 1], 0.22500000000000006], 
										[[0, 0.6485013623978204, 0.3514986376021796], 0.05476839237057218], 
										[[0.9754098360655736, 0, 0.024590163934426246], 0.0006557377049180337], 
										[[0.9799999999999999, 0, 0.02000000000000001], 0]];
		
		var covMat = new PortfolioAllocation.Matrix([[0.0004, 0.0004, 0.0002],
													 [0.0004, 0.0025,0.001],
													  [0.0002, 0.001, 0.01]]);
		var returns = new PortfolioAllocation.Matrix([0.05, 0.08, 0.12]);
		
		var cornerPortfolios = PortfolioAllocation.computeCornerPortfolios_(returns, covMat);
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #3');
	}
	
	// Test using static data (upper bounds)
	// Reference: A Simple Spreadsheet-Based Exposition of the Markowitz Critical Line Method for Portfolio Selection, Clarence C. Kwan
	{
		var expectedCornerPortfolios = [[[0, 0.30000000000000004, 0.7], 0.14625], 
										[[0, 0.6485013623978203, 0.3514986376021798], 0.05476839237057221], 
										[[0.7, 0.18310626702997274, 0.11689373297002724], 0.015934604904632152], 
										[[0.7, 0.2438095238095238, 0.05619047619047619], 0]];
		
		var covMat = new PortfolioAllocation.Matrix([[0.0004, 0.0004, 0.0002],
													 [0.0004, 0.0025,0.001],
													  [0.0002, 0.001, 0.01]]);
		var returns = new PortfolioAllocation.Matrix([0.05, 0.08, 0.12]);
		
		var cornerPortfolios = PortfolioAllocation.computeCornerPortfolios_(returns, covMat, { constraints: {maxWeights: [0.7, 0.7, 0.7]} });
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #4');
	}
	
	// Test using static data
	// Reference: An Open-Source Implementation of the Critical-Line Algorithm for Portfolio Optimization, David H. Bailey and Marcos Lopez de Prado
	{
		var expectedCornerPortfolios = [[[0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 58.30308533333371], 
										[[0.6493694070931811, 0.3506305929068189, 0, 0, 0, 0, 0, 0, 0, 0], 4.1742728458857385], 
										[[0.4339841341086239, 0.23124750065448754, 0, 0.334768365236889, 0, 0, 0, 0, 0, 0], 1.9455661414558894], 
										[[0.12688785385570883, 0.07234334721032556, 0, 0.28125374926334057, 0, 0, 0, 0, 0, 0.5195150496706249], 0.16458117494477595],
										[[0.12320100405906734, 0.07044407130753655, 0, 0.2789935668090118, 0, 0, 0, 0.006435564362887149, 0, 0.5209257934614971], 0.1473887508934171],
										[[0.0869215492990579, 0.050451042268558385, 0, 0.22359401742288823, 0, 0.17383161507156486, 0, 0.03017301555135618, 0, 0.4350287603865743], 0.056172204002751545],
										[[0.0846709411996219, 0.049253858741118525, 0, 0.21963390336360733, 0, 0.18003923464176064, 0, 0.03102980185535347, 0.006485702415438152, 0.42888655778310003], 0.05204819067458028],
										[[0.07378925302280315, 0.043828660769718863, 0, 0.19897560805881487, 0.026158159857441972, 0.19815187227970524, 0, 0.03341958639919798, 0.027902966026643668, 0.3977738935856743], 0.03652161374727064],
										[[0.06834400480527462, 0.041387026820649334, 0.015215259551836627, 0.18813443107045838, 0.03416248599274816, 0.20231943214747125, 0, 0.0339293235595669, 0.03363264959172938, 0.38287538646026537], 0.030971168861678777],
										[[0.03696858147921504, 0.02690083780081047, 0.0949424305647986, 0.1257759521946726, 0.0767460810325476, 0.21935567131616898, 0.029987096882220312, 0.035963284621386274, 0.06134983772972688, 0.29201022637845325], 0]];
										
		var covMat = new PortfolioAllocation.Matrix([[0.40755159,0.03175842,0.05183923,0.05663904,0.0330226,0.00827775,0.02165938,0.01332419,0.0343476,0.02249903],
													[0.03175842,0.9063047,0.03136385,0.02687256,0.01917172,0.00934384,0.02495043,0.00761036,0.02874874,0.01336866],
													[0.05183923,0.03136385,0.19490901,0.04408485,0.03006772,0.01322738,0.03525971,0.0115493,0.0427563,0.02057303],
													[0.05663904,0.02687256,0.04408485,0.19528471,0.02777345,0.00526665,0.01375808,0.00780878,0.02914176,0.01640377],
													[0.0330226,0.01917172,0.03006772,0.02777345,0.34059105,0.00777055,0.02067844,0.00736409,0.02542657,0.01284075],
													[0.00827775,0.00934384,0.01322738,0.00526665,0.00777055,0.15983874,0.02105575,0.00518686,0.01723737,0.00723779],
													[0.02165938,0.02495043,0.03525971,0.01375808,0.02067844,0.02105575,0.68056711,0.01377882,0.04627027,0.01926088],
													[0.01332419,0.00761036,0.0115493,0.00780878,0.00736409,0.00518686,0.01377882,0.95526918,0.0106553,0.00760955],
													[0.0343476,0.02874874,0.0427563,0.02914176,0.02542657,0.01723737,0.04627027,0.0106553,0.31681584,0.01854318],
													[0.02249903,0.01336866,0.02057303,0.01640377,0.01284075,0.00723779,0.01926088,0.00760955,0.01854318,0.11079287]]);
		var returns = new PortfolioAllocation.Matrix([1.175,1.19,0.396,1.12,0.346,0.679,0.089,0.73,0.481,1.08]);

		var cornerPortfolios = PortfolioAllocation.computeCornerPortfolios_(returns, covMat);
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #5');
	}
});


QUnit.test('Mean variance portfolio - internal target return weights portfolio', function(assert) {    
	function generateRandomValue(minVal, maxVal) {	
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test using random data
	{
		// Problem data
		var covMat =[[0.0146, 0.0187, 0.0145],
					[0.0187, 0.0854, 0.0104],
					[0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Compute the target return at random
		var minReturn = returns[0];
		var maxReturn = returns[1];
		var targetReturn = generateRandomValue(minReturn, maxReturn);
		
		// Compute the associated portfolio weights
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, {  optimizationMethod: 'targetReturn', constraints: {return: targetReturn}});
		
		// Compare the computed portfolio return with the target return
		var portfolioReturn = PortfolioAllocation.Matrix.vectorDotProduct(new PortfolioAllocation.Matrix(returns), new PortfolioAllocation.Matrix(weights));
		assert.equal(Math.abs(portfolioReturn - targetReturn) <= 1e-8, true, 'Mean variance portfolio - internal target return weights portfolio #1');
	}
});	


QUnit.test('Mean variance portfolio - internal target volatility weights portfolio', function(assert) {    
	function generateRandomValue(minVal, maxVal) {	
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test using random data
	{
		// Problem data
		var covMat = [[0.0146, 0.0187, 0.0145],
					[0.0187, 0.0854, 0.0104],
					[0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Compute the target return at random
		var maxVolatility = Math.sqrt(covMat[1][1]);
		var minVolatility = Math.sqrt(0.014599310344827589); // computed thanks to the efficient frontier
		var targetVolatility = generateRandomValue(minVolatility, maxVolatility);
		
		// Compute the associated portfolio weights
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, {  optimizationMethod: 'targetVolatility', constraints: {volatility: targetVolatility}});
		
		// Compare the computed portfolio volatility with the target volatility
		var portfolioVolatility = Math.sqrt(PortfolioAllocation.Matrix.vectorDotProduct(PortfolioAllocation.Matrix.xy(new PortfolioAllocation.Matrix(covMat), new PortfolioAllocation.Matrix(weights)), 
																			            new PortfolioAllocation.Matrix(weights)));
		assert.equal(Math.abs(portfolioVolatility - targetVolatility) <= 1e-8, true, 'Mean variance portfolio - internal target volatility weights portfolio #1');
	}
});	


QUnit.test('Random subspace mean variance portfolio - internal target max volatility weights portfolio', function(assert) {    
	function generateRandomValue(minVal, maxVal) {	
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test using random data
	{
		// Problem data
		var covMat = [[0.0146, 0.0187, 0.0145],
					[0.0187, 0.0854, 0.0104],
					[0.0145, 0.0104, 0.0289]];
		var returns = [0.062, 0.146, 0.128];
		
		// Compute the target return at random
		var maxVolatility = Math.sqrt(covMat[1][1]);
		var minVolatility = 0; // in case of RSO-MVO, the minimum volatility can be 0; in this case, the final portfolio will be empty
		var targetMaxVolatility = generateRandomValue(minVolatility, maxVolatility);
		
		// Compute the associated portfolio weights
		var weights = PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, { subsetPortfolioOptimizationMethodParams: { 
		                                                                                                      optimizationMethod: 'maximumTargetVolatility', 
																											  constraints: {
																												  maxVolatility: targetMaxVolatility 
																											  }
		                                                                                                    }} );
		
		// Compare the computed portfolio volatility with the target volatility
		var portfolioVolatility = Math.sqrt(PortfolioAllocation.Matrix.vectorDotProduct(PortfolioAllocation.Matrix.xy(new PortfolioAllocation.Matrix(covMat), new PortfolioAllocation.Matrix(weights)), 
																			            new PortfolioAllocation.Matrix(weights)));
		assert.equal(portfolioVolatility <= targetMaxVolatility, true, 'Random subspace mean variance portfolio - internal max target volatility weights portfolio #1');
	}
});	
