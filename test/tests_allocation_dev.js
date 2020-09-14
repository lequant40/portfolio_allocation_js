// ------------------------------------------------------------
QUnit.module('Assets allocation module');
	
QUnit.test('Mean variance portfolio - internal corner portfolios computation', function(assert) {    
	// Used for compatibility purposes with old JavaScript engines
	function cornerPortfoliosToArray(cornerPortfolio) {
		for (var i = 0; i < cornerPortfolio.length; ++i) {
			// WeightscornerPortfolio
			cornerPortfolio[i][0] = cornerPortfolio[i][0].toArray();
		}
		
		return cornerPortfolio;
	}

	// Test using static data
	// Test infeasible/unsupported cases for determining the critical line algorithm starting portfolio:
	{
		// Lower bounds > Upper bounds
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														             { constraints: {minWeights: [0.6, 0.3], maxWeights: [0.2, 1]} }) },
			new Error('infeasible problem detected: lower bound strictly greater than upper bound'),
			"Mean variance portfolio - Corner portfolios, lower bounds greater than upper bounds");
			
		//  Sum lower bounds > 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.2]), 
														 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														 { constraints: {minWeights: [0.6, 0.5], maxWeights: [0.7, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Mean variance portfolio - Corner portfolios, sum of lower bounds greater than one");
		
		//  Sum upper bounds < 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.2]), 
														new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
														{ constraints: {minWeights: [0.4, 0.4], maxWeights: [0.4, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Mean variance portfolio - Corner portfolios, sum of upper bounds lower than one");
		
		// Identical returns, same returns, equal returns
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.1]), 
														 new PortfolioAllocation.Matrix([[1,0],[0,1]])) },
			new Error('unsupported problem: equal returns are not supported by the critical line algorithm'),
			"Mean variance portfolio - Corner portfolios, identical returns");
	}
	
	// Test using static data
	// Test limit cases for determining the critical line algorithm starting portfolio:
	{
		// Lower bounds binding (sum lb_i == 1)
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																						new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																						{ constraints: {minWeights: [0.4, 0.6]} }).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), [[[0.4, 0.6],  0]], 'Mean variance portfolio - Corner portfolios, lower bounds binding');
		
		// Upper bounds binding (sum ub_i == 1)
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																						new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																						{ constraints: {maxWeights: [0.6, 0.4]} }).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), [[[0.6, 0.4],  0]], 'Mean variance portfolio - Corner portfolios, upper bounds binding');
	}
	
	// Test using static data
	// Reference: https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm
	//
	// Note: this test also allows checking that numerically equal corner portfolios are NOT filtered out by the algorithm,
	// c.f. note: "over the range of risk tolerances from 22.30 to 22.94 the optimal composition remains the same. This is
    // not a rounding error. Since the efficient frontier is piecewise quadratic, there is always the possibility that there 
	// is a kink at the point corresponding to a specific level of risk and return. In such a case indifference curves with different 
	// slopes (risk tolerances) can be tangent to the efficient frontier at the same point, giving the same optimal portfolio."
	{
		var expectedCornerPortfolios = [[[0.2, 0.30000000000000004, 0.5], 20.898844444444443], 
		                                [[0.2,0.5,0.3000000000000001],11.470044444444447],
										[[0.2, 0.5, 0.30000000000000004], 11.1475], 
										[[0.22180737780348653, 0.5, 0.27819262219651353], 10.51088812347172],
										[[0.451915610952186, 0.348084389047814, 0.2], 7.55192170004087],
										[[0.5, 0.30000000000000004, 0.2], 6.8672],
										[[0.5, 0.2999999999999999, 0.2], 0]];
		
		var covMat = new PortfolioAllocation.Matrix([[1, 2.96, 2.31],
													[2.96, 54.76, 39.886],
													[2.31,	39.886,	237.16]]);
		var returns = new PortfolioAllocation.Matrix([2.8000, 6.3000, 10.8000]);
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, { constraints: {minWeights: [0.2, 0.2, 0.2], maxWeights: [0.5, 0.5, 0.5]} }).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #1');
	}
		
	// Test using static data
	// Reference: Portfolio Selection, H. Markowitz example, chapter VIII, section "The computing procedure"
	{
		var expectedCornerPortfolios = [[[0, 1, 0], 4.16666666666667], 
										[[0, 0.22496808316614988, 0.7750319168338501], 0.1408064320019454], 
										[[0.8414051841746248, 0, 0.15859481582537516], 0.03332764893133244], 
										[[0.9931034482758623, 0, 0.006896551724137813], 0]];
				
		var covMat = new PortfolioAllocation.Matrix([[0.0146, 0.0187, 0.0145],
													 [0.0187, 0.0854, 0.0104],
													  [0.0145, 0.0104, 0.0289]]);
		var returns = new PortfolioAllocation.Matrix([0.062, 0.146, 0.128]);
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #2');
	}
	
	// Test using static data that partial investment constraint is properly managed
	{
		var expectedCornerPortfolios = [[[0, 1, 0], 4.16666666666667], 
										[[0, 0.23479391919356146, 0.7652060808064386], 0.19184619136655556],
										[[0, 0, 0], 0]];
		
		var covMat = new PortfolioAllocation.Matrix([[0.0146, 0.0187, 0.0145],
													 [0.0187, 0.0854, 0.0104],
													  [0.0145, 0.0104, 0.0289]]);
		var returns = new PortfolioAllocation.Matrix([0.062, 0.146, 0.128]);
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, {constraints: {fullInvestment: false}}).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #2, partial investment constraint');
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
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat).cornerPortfolios;
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
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, { constraints: {maxWeights: [0.7, 0.7, 0.7]} }).cornerPortfolios;
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

		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #5');
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
		
		var expectedCornerPortfolios = [[[0.0008860751090370526, 0, 0, 0.003345695950238685, 0.00007663523943489258, 0.012841518880073127, 0.01192011965466721, 0.97092995516609], 0],
										[[0.0010092608916674394, 0, 0, 0.004008702066904519, 0.00033764763422769743, 0.012496184295674449, 0.01453557759610034, 0.9676126275149418],  0.00046956768316411246], 
										[[0.0017655557593660093, 0, 0.008931031081884725, 0, 0.008263675566746306, 0.0063916733651148624, 0.04309168778365437, 0.9315563764426312], 0.0041508392817256245], 
										[[0.0012756123767303507, 0, 0.0032588378760465275, 0, 0.004072943929600115, 0.012021095630643686, 0.02160446204652544, 0.957767048138976], 0.0014586035429539558], 
										[[0, 0.004345917988265768, 0.006849418265594814, 0, 0.004007054238260629, 0.005038387387097269, 0.0284123197082556, 0.9513469024125258], 0.002851415009160898],
										[[0, 0.007621949314523461, 0.01148109322047808, 0, 0.00653096875456001, 0, 0.04935684701395931, 0.9250091416964791], 0.006507814480398172],
										[[0, 0.061938538867079784, 0.10166770263898256, 0, 0.11119983640121779, 0, 0.7251939220927198, 0],  0.12699435118749036],
										[[0, 0.25, 0.0405985104054761, 0, 0.009674272143084611, 0, 0.6997272174514392, 0], 0.2637718357009199],
										[[0, 0.25, 0.03158049051256315, 0, 0.023733980661979903, 0, 0.6946855288254571, 0], 0.39912101396530497],
										[[0.046211968331688, 0.25, 0.008523783192019019, 0, -2.7755575615628914e-17, 0, 0.6952642484762929, 0], 0.44342871837702696],
										[[0.05663429563220421, 0.25, 0, 0, -2.7755575615628914e-17, 0, 0.6933657043677959, 0], 0.4576339467556026],
										[[0.25, 0.25, 0, 0, -2.7755575615628914e-17, 0, 0.5, 0], 1.8766716337311202],
										[[0.25,0.25,0,0,-2.7755575615628914e-17,0,0.5,0], 8.899361909160644],
										[[0.25, 0.25, 0, 0, 0.25, 0, 0.25, 0], 18.562065370001065]].reverse();
						
		// Test that the algorithm is behaving properly
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																						new PortfolioAllocation.Matrix(covMat), 
																						{ constraints: { minWeights: minWeights, maxWeights: maxWeights }}).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #6');
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
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, {  constraints: {return: targetReturn}});
		
		// Compare the computed portfolio return with the target return
		var portfolioReturn = PortfolioAllocation.Matrix.vectorDotProduct(new PortfolioAllocation.Matrix(returns), new PortfolioAllocation.Matrix(weights));
		assert.equal(Math.abs(portfolioReturn - targetReturn) <= 1e-4, true, 'Internal target return weights portfolio #1');
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
		var weights = PortfolioAllocation.meanVarianceOptimizationWeights(returns, covMat, {  constraints: {volatility: targetVolatility}});
		
		// Compare the computed portfolio volatility with the target volatility
		var portfolioVolatility = Math.sqrt(PortfolioAllocation.Matrix.vectorDotProduct(PortfolioAllocation.Matrix.xy(new PortfolioAllocation.Matrix(covMat), new PortfolioAllocation.Matrix(weights)), 
																			            new PortfolioAllocation.Matrix(weights)));
		assert.equal(Math.abs(portfolioVolatility - targetVolatility) <= 1e-4, true, 'Internal target volatility weights portfolio #1');
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
		var weights = PortfolioAllocation.randomSubspaceMeanVarianceOptimizationWeights(returns, covMat, { subsetsOpt: { 
																											  constraints: {
																												  fullInvestment: false,
																												  maxVolatility: targetMaxVolatility
																											  }
		                                                                                                    }} );
		
		// Compare the computed portfolio volatility with the target volatility
		var portfolioVolatility = Math.sqrt(PortfolioAllocation.Matrix.vectorDotProduct(PortfolioAllocation.Matrix.xy(new PortfolioAllocation.Matrix(covMat), new PortfolioAllocation.Matrix(weights)), 
																			            new PortfolioAllocation.Matrix(weights)));
		assert.equal(portfolioVolatility <= targetMaxVolatility, true, 'Random subspace mean variance portfolio - internal max target volatility weights portfolio #1');
	}
});	



QUnit.test('Random weights portfolio - internal tests requiring access to internal functions', function(assert) {    
	// Test with random data, min and max weights constraints
	{
	  // Setup static parameters of the random test
	  var nbTests = 50;
	  var nbAssetsMin = 1;
	  var nbAssetsMax = 50;
	  
	  // Aim of these tests is to check that the generated weights are compatible with lower, upper and 
	  // lower and upper bounds.
	  for (var i = 0; i < nbTests; ++i) {
		  // Generate a random number of assets
		  var nbAssets = Math.floor(Math.random()*(nbAssetsMax - nbAssetsMin + 1) + nbAssetsMin);

		  
		  // Generate random feasible lower bounds with sum l_i = k with k < 1.
		  var integerBounds = new PortfolioAllocation.randomCompositionsIterator_(100, nbAssets).next();
		  var feasibleLowerBounds = new Array(nbAssets);
		  var feasibleLowerBoundsFactor = Math.random()*(0.9 - 0.1) + 0.1;
		  for (var j = 0; j < nbAssets; ++j) {
			  feasibleLowerBounds[j] = integerBounds[j] * feasibleLowerBoundsFactor/100;
		  }
		   
		  // Generate a random portfolio with lower bounds constraints
		  var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { minWeights: feasibleLowerBounds } });
		  
			// Check that the weights belong to the interval [l_i, 1] in case the assets have been selected
			var weightsBelongInterval = true;
			for (var k = 0; k < randomWeights.length; ++k) {
				if (randomWeights[k] != 0 && (randomWeights[k] > 1 || randomWeights[k] < feasibleLowerBounds[k])) {
					weightsBelongInterval = false;
					break;
				}
			}
			assert.equal(weightsBelongInterval, true, "Random weights portfolio with lower bounds constraints - Test " + i);


		  // Generate random feasible upper bounds with sum u_i = k with k > 1.
		  var integerBounds = new PortfolioAllocation.randomCompositionsIterator_(100, nbAssets).next();
		  var feasibleUpperBounds = new Array(nbAssets);
		  var feasibleUpperBoundsFactor = Math.random()*(1.9 - 1.1) + 1.1;
		  for (var j = 0; j < nbAssets; ++j) {
			  feasibleUpperBounds[j] = Math.min(1, integerBounds[j] * feasibleUpperBoundsFactor/100);
		  }

		  // Generate a random portfolio with upper bounds constraints
		  var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { maxWeights: feasibleUpperBounds } });

			// Check that the weights belong to the interval [0, u_i]
			var weightsBelongInterval = true;
			for (var k = 0; k < randomWeights.length; ++k) {
				if (randomWeights[k] != 0 && (randomWeights[k] > feasibleUpperBounds[k] || randomWeights[k] < 0)) {
					weightsBelongInterval = false;
					break;
				}
			}
			assert.equal(weightsBelongInterval, true, "Random weights portfolio with upper bounds constraints - Test " + i);

			
		  // Generate random feasible upper bounds with sum u_i = k with k > 1.
		  //
		  // Ensure u_i >= l_i.
		  var integerBounds = new PortfolioAllocation.randomCompositionsIterator_(100, nbAssets).next();
		  var feasibleUpperBounds = new Array(nbAssets);
		  var feasibleUpperBoundsFactor = Math.random()*(1.9 - 1.1) + 1.1;
		  for (var j = 0; j < nbAssets; ++j) {
			  feasibleUpperBounds[j] = Math.max(feasibleLowerBounds[j], Math.min(1, integerBounds[j] * feasibleUpperBoundsFactor/100));
	      }
		  		  
		  // Generate a random portfolio with lower and upper bounds constraints
		  var randomWeights = PortfolioAllocation.randomWeights(nbAssets, { constraints: { minWeights: feasibleLowerBounds, maxWeights: feasibleUpperBounds } });
		  
			// Check that the weights belong to the interval [l_i, u_i]
			var weightsBelongInterval = true;
			for (var k = 0; k < randomWeights.length; ++k) {
				if (randomWeights[k] != 0 && (randomWeights[k] > feasibleUpperBounds[k] || randomWeights[k] < feasibleLowerBounds[k])) {
					weightsBelongInterval = false;
					break;
				}
			}
			assert.equal(weightsBelongInterval, true, "Random weights portfolio with lower and upper bounds constraints - Test " + i);
	  }
	}
});


QUnit.test('Mean variance optimization - computation of maximum risk tolerance portfolio with the GSMO algorithm', function(assert) {    	
	// Test using static data
	// Test infeasible/unsupported cases:
	{
		// Lower bounds > Upper bounds
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	 { constraints: {minWeights: [0.6, 0.3], maxWeights: [0.2, 1]} }) },
			new Error('infeasible problem detected: lower bound strictly greater than upper bound'),
			"Error case - Lower bounds greater than upper bounds");
			
		//  Sum lower bounds > 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	 { constraints: {minWeights: [0.6, 0.5], maxWeights: [0.7, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Error case - Sum of lower bounds greater than one");
		
		//  Sum upper bounds < 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	{ constraints: {minWeights: [0.4, 0.4], maxWeights: [0.4, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Error case - Sum of upper bounds lower than one");
	}

	// Test using static data that assets with identical returns are properly managed
	{
		var covMat = new PortfolioAllocation.Matrix([[1, 0, 0],
													[0, 1, 0],
													[0,	0,	1]]);
		var returns = new PortfolioAllocation.Matrix([1, 1, 1]);

		var expectedPortfolio = [[0.3333333333333333, 0.3333333333333333, 0.3333333333333333], 1]; // True risk tolerance is 0

		var efficientFrontier = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat);
		var portfolio = efficientFrontier.getHighestRiskTolerancePortfolio();
		var lambda = efficientFrontier.getHighestRiskTolerance();
		
		assert.deepEqual(portfolio.toArray(), expectedPortfolio[0], "Equal assets returns, portfolio weights");
		assert.equal(lambda, expectedPortfolio[1], "Equal assets returns, risk tolerance");		
	}
	
	// Test using static data
	// Reference: https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm
	{
		var covMat = new PortfolioAllocation.Matrix([[1, 2.96, 2.31],
													[2.96, 54.76, 39.886],
													[2.31,	39.886,	237.16]]);
		var returns = new PortfolioAllocation.Matrix([2.8000, 6.3000, 10.8000]);
		
		var expectedPortfolio = [[0.2, 0.29999999999999993, 0.5], 32]; // True risk tolerance is 20.898844444444443

		var efficientFrontier = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat, { constraints: {minWeights: [0.2, 0.2, 0.2], maxWeights: [0.5, 0.5, 0.5]} });
		var portfolio = efficientFrontier.getHighestRiskTolerancePortfolio();
		var lambda = efficientFrontier.getHighestRiskTolerance();
			
		assert.deepEqual(portfolio.toArray(), expectedPortfolio[0], "Bounds constraints, portfolio weights");
		assert.equal(lambda, expectedPortfolio[1], "Bounds constraints, risk tolerance");		
	}
});

QUnit.test('Mean variance optimization - computation of minimum risk tolerance portfolio with the GSMO algorithm', function(assert) {    	
	// Test using static data
	// Test infeasible/unsupported cases for determining the critical line algorithm starting portfolio:
	{
		// Lower bounds > Upper bounds
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	 { constraints: {minWeights: [0.6, 0.3], maxWeights: [0.2, 1]} }) },
			new Error('infeasible problem detected: lower bound strictly greater than upper bound'),
			"Error case - Lower bounds greater than upper bounds");
			
		//  Sum lower bounds > 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	 new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	 { constraints: {minWeights: [0.6, 0.5], maxWeights: [0.7, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Error case - Sum of lower bounds greater than one");
		
		//  Sum upper bounds < 1
		assert.throws(function() { 
			new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(new PortfolioAllocation.Matrix([0.1, 0.2]), 
																	new PortfolioAllocation.Matrix([[1,0],[0,1]]), 
																	{ constraints: {minWeights: [0.4, 0.4], maxWeights: [0.4, 0.5]} }) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Error case - Sum of upper bounds lower than one");
	}

	// Test using static data that assets with identical returns are properly managed
	{
		var covMat = new PortfolioAllocation.Matrix([[1, 0, 0],
													[0, 1, 0],
													[0,	0,	1]]);
		var returns = new PortfolioAllocation.Matrix([1, 1, 1]);

		var expectedPortfolio = [[0.3333333333333333, 0.3333333333333333, 0.3333333333333333], 1]; // True risk tolerance is 0

		var efficientFrontier = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat);
		var portfolio = efficientFrontier.getLowestRiskTolerancePortfolio();
		var lambda = efficientFrontier.getLowestRiskTolerance();

		assert.deepEqual(portfolio.toArray(), expectedPortfolio[0], "Equal assets returns, portfolio weights");
		assert.equal(lambda, expectedPortfolio[1], "Equal assets returns, risk tolerance");		
	}
	
	// Test using static data
	// Reference: https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm
	{
		var covMat = new PortfolioAllocation.Matrix([[1, 2.96, 2.31],
													[2.96, 54.76, 39.886],
													[2.31,	39.886,	237.16]]);
		var returns = new PortfolioAllocation.Matrix([2.8000, 6.3000, 10.8000]);
		
		var expectedPortfolio = [[0.5, 0.29999999999999993, 0.2], 1]; // True risk tolerance is 0

		var efficientFrontier = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat, { constraints: {minWeights: [0.2, 0.2, 0.2], maxWeights: [0.5, 0.5, 0.5]} });
		var portfolio = efficientFrontier.getLowestRiskTolerancePortfolio();
		var lambda = efficientFrontier.getLowestRiskTolerance();
		
		assert.deepEqual(portfolio.toArray(), expectedPortfolio[0], "Bounds constraints, portfolio weights");
		assert.equal(lambda, expectedPortfolio[1], "Bounds constraints, risk tolerance");		
	}
	
	// Test using static data
	// Test that in case of semi-positive definite covariance matrix, the computed portfolio is not only of minimum volatility, 
	// but efficient (i.e., maximizes the return).
	//
	// Here, portfolio [2/3, 1/3, 0] would also be minimizing the volatility, but would not be efficient for instance.
	{
		var covMat = new PortfolioAllocation.Matrix([[1, 1, 1],
													 [1, 1, 1],
													 [1, 1, 3]]);
		var returns = new PortfolioAllocation.Matrix([1, 2, 3]);

		var expectedPortfolio = [[0, 0.9999990463256834, 9.536743164617612e-7], 0.0000019073486328125]; // True risk tolerance is 0

		var efficientFrontier = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat);
		var portfolio = efficientFrontier.getLowestRiskTolerancePortfolio();
		var lambda = efficientFrontier.getLowestRiskTolerance();
		
		assert.deepEqual(portfolio.toArray(), expectedPortfolio[0], "Semi-definite positive covariance matrix, portfolio weights");
		assert.equal(lambda, expectedPortfolio[1], "Semi-definite positive covariance matrix, risk tolerance");	
	}	
});

QUnit.test('Mean variance optimization - maximum Sharpe ratio internal computation', function(assert) {    	
	// Test using static data
	// Reference: An Open-Source Implementation of the Critical-Line Algorithm for Portfolio Optimization, David H. Bailey and Marcos Lopez de Prado
	{
		// Note: this portfolio has a Sharpe ratio of ~4.4535(3...), with a null risk free rate, which
		// is exactly the same value as in the reference.
		var expectedPortfolio = [[0.08397318948217423, 0.04890598613711377, 0, 0.21830925954049477, 0.0016773041709357821,
							     0.18120064671441183, 0, 0.031183038765169507, 0.007859012532850177, 0.42689156265685], 0.05105260106191406];
										
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
		
		var efficientFrontierCla = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat);
		var maxSharpeRatioPortfolioCla = efficientFrontierCla.computeMaximumSharpeRatioEfficientPortfolio(rf);
		assert.deepEqual(maxSharpeRatioPortfolioCla[0].toArray(), expectedPortfolio[0], "Test #1, portfolio weights, critical line");
		assert.equal(maxSharpeRatioPortfolioCla[1], expectedPortfolio[1], "Test #1, risk tolerance, critical line");
		
		var maxSharpeRatioPortfolioWeightsGsmo = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat);
		var maxSharpeRatioPortfolioGsmo = maxSharpeRatioPortfolioWeightsGsmo.computeMaximumSharpeRatioEfficientPortfolio(rf);

		var weightsOK = true;
		if (expectedPortfolio[0].length != maxSharpeRatioPortfolioGsmo[0].nbRows) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedPortfolio[0].length; ++i) {
				if (Math.abs(maxSharpeRatioPortfolioGsmo[0].data[i] - expectedPortfolio[0][i]) > 1e-5) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, "Test #1, portfolio weights, GSMO");
		assert.equal(Math.abs(maxSharpeRatioPortfolioGsmo[1] - expectedPortfolio[1]) <= 1e-6, true, "Test #1, risk tolerance, GSMO");
	}
	
	
	// Test using static data that the maximum Sharpe ratio is properly computed in a
	// case where it is constant on the right of the risk tolerance interval [0, 1].
	//
	// Reference: Properties of the most diversified portfolio
	{
		var covMat = [[0.0400, 0.0100], [0.0100, 0.0100]];
		var returns = [0.2, 0.1];
		var rf  = 0;
		
		var expectedPortfolio = [[0.33333333333333337, 0.6666666666666666], 0.1];		
		
		var efficientFrontierGsmo = new PortfolioAllocation.MeanVarianceEfficientFrontierGsmo(returns, covMat, {optimizationMethod: 'gsmo'});
		var maxSharpeRatioPortfolioGsmo = efficientFrontierGsmo.computeMaximumSharpeRatioEfficientPortfolio(rf);
		
		var efficientFrontierCla = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, {optimizationMethod: 'critical-line'});
		var maxSharpeRatioPortfolioCla = efficientFrontierCla.computeMaximumSharpeRatioEfficientPortfolio(rf);
		
		var weightsOK = true;
		if (expectedPortfolio[0].length != maxSharpeRatioPortfolioGsmo[0].nbRows || 
		    expectedPortfolio[0].length != maxSharpeRatioPortfolioCla[0].nbRows) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedPortfolio[0].length; ++i) {
				if (Math.abs(maxSharpeRatioPortfolioGsmo[0].data[i] - expectedPortfolio[0][i]) > 1e-6 ||
				    Math.abs(maxSharpeRatioPortfolioCla[0].data[i] - expectedPortfolio[0][i]) > 1e-12) {
					weightsOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, "Test #2, portfolio weights");
		assert.equal(Math.abs(maxSharpeRatioPortfolioGsmo[1] - expectedPortfolio[1]) <= 1e-6 && Math.abs(maxSharpeRatioPortfolioCla[1] - expectedPortfolio[1]) <= 1e-12, true, "Test #2, risk tolerance");
	}
});


