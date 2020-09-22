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

	
		// Test using data for equal returns, with 2 assets IN for the maximum return portfolio and associated KKT system not invertible,
		// but solvable by least squares
		{
			var covMat = [[0.0400, 0.0400], [0.0400, 0.0400]]; // determinant = 0
			var returns = [1, 1];
			assert.throws(function() { 
				var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																								new PortfolioAllocation.Matrix(covMat)); },
				new Error('internal error: impossible to solve the KKT system'),
				"Corner portfolios, non invertible KKT matrix");
		}
		
		// Test using data for equal returns, with 2 assets IN for the maximum return portfolio and associated KKT system not invertible,
		// (solvable by least squares, but cycling afterwards)
		{
			var covMat = [[0.1, 0.1, -0.2], [0.1, 0.1, -0.1], [-0.2, -0.1, 1]];
			var returns = [2, 2, 1];
		
			assert.throws(function() { 
				new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																		 new PortfolioAllocation.Matrix(covMat)) },
				new Error('internal error: impossible to solve the KKT system'),
				"Corner portfolios, non invertible KKT matrix, and cycling if done so");		 
		}
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
										[[0.2218073778034857, 0.5, 0.27819262219651436], 10.510888123471753],
										[[ 0.45191561095218635, 0.3480843890478137, 0.2], 7.551921700040866],
										[[0.5, 0.30000000000000016, 0.2], 6.867200000000002],
										[[0.5, 0.30000000000000004, 0.2], 0]];
		
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
										[[0, 0.22496808316614986, 0.7750319168338502], 0.14080643200194537], 
										[[0.8414051841746242, 0, 0.15859481582537574], 0.03332764893133249], 
										[[0.993103448275862, 0, 0.006896551724137945], 0]];
				
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
										[[0, 0.23479391919356143, 0.7652060808064386], 0.19184619136655556],
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
										[[0, 0.6485013623978203, 0.3514986376021797], 0.0547683923705722], 
										[[0.9754098360655737, 0, 0.024590163934426305], 0.0006557377049180432], 
										[[0.98, 0, 0.020000000000000018], 0]];
		
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
										[[0, 0.6485013623978203, 0.3514986376021797], 0.0547683923705722], 
										[[0.7, 0.18310626702997276, 0.11689373297002731], 0.015934604904632162], 
										[[0.7, 0.24380952380952386, 0.05619047619047618], 0]];
		
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
										[[0.649369407093181, 0.35063059290681897, 0, 0, 0, 0, 0, 0, 0, 0], 4.1742728458857545], 
										[[0.4339841341086233, 0.23124750065448724, 0, 0.3347683652368895, 0, 0, 0, 0, 0, 0], 1.945566141455883], 
										[[0.12688785385570886, 0.07234334721032556, 0, 0.2812537492633406, 0, 0, 0, 0, 0, 0.519515049670625], 0.1645811749447761],
										[[0.1232010040590673, 0.07044407130753662, 0, 0.2789935668090117, 0, 0, 0, 0.00643556436288719, 0, 0.5209257934614973], 0.1473887508934171],
										[[0.0869215492990579, 0.0504510422685584, 0, 0.22359401742288834, 0, 0.1738316150715647, 0, 0.030173015551356194, 0, 0.43502876038657434], 0.05617220400275158],
										[[0.08467094119962172, 0.04925385874111865, 0, 0.21963390336360714, 0, 0.18003923464176075, 0, 0.031029801855353523, 0.006485702415438527, 0.4288865577830998], 0.052048190674580205],
										[[0.07378925302280318, 0.043828660769718863, 0, 0.19897560805881462, 0.026158159857441944, 0.1981518722797053, 0, 0.033419586399198, 0.027902966026643647, 0.3977738935856744], 0.036521613747270656],
										[[0.06834400480527457, 0.0413870268206493, 0.015215259551836766, 0.1881344310704583, 0.0341624859927482, 0.2023194321474713, 0, 0.03392932355956689, 0.033632649591729444, 0.3828753864602651], 0.030971168861678753],
										[[0.036968581479215, 0.026900837800810484, 0.09494243056479873, 0.1257759521946725, 0.07674608103254753, 0.219355671316169, 0.029987096882220284, 0.03596328462138633, 0.06134983772972707, 0.2920102263784531], 0]];
										
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
	
	// Reference: Avoiding the Downside: A Practical Review of the Critical Line Algorithm for Mean-Semivariance Portfolio Optimization
	{
		var expectedCornerPortfolios = [[[0.1, 0.5, 0.4], 1.7567], 
										[[0.1000, 0.4000, 0.5000], 1.2203], 
										[[0.1000, 0.4000, 0.5000], 0.3142], 
										[[0.3764, 0.1236, 0.5000], 0.0973], 
										[[0.4644, 0.1000, 0.4356], 0.0853], 
										[[0.5000, 0.1000, 0.4000], 0.0770], 
										[[0.5000, 0.1000, 0.4000], 0]];

		var returns = [[-0.173,0.098,0.200,0.030,-0.183,0.067,0.300,0.103,0.216,-0.046,-0.071,0.056,0.038,0.089,0.090,0.083,0.035,0.176],
					   [-0.318,0.285,-0.047,0.104,-0.171,-0.039,0.149,0.260,0.419,-0.078,0.169,-0.035,0.133,0.732,0.021,0.131,0.006,0.908],
					   [-0.319,0.076,0.381,-0.051,0.087,0.262,0.341,0.227,0.352,0.153,-0.099,0.038,0.273,0.091,0.054,0.109,0.210,0.112]];
		var covMat = PortfolioAllocation.covarianceMatrix(returns, {method: "sample-covariance"});
		var returns = PortfolioAllocation.meanVector(returns);
		var lb = [0.1, 0.1, 0.1];
		var ub = [0.5, 0.5, 0.5];
		
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, {constraints: {minWeights: lb, maxWeights: ub}}).cornerPortfolios;
		
		var weightsOK = true;
		var lambdaOK = true;
		if (expectedCornerPortfolios.length != cornerPortfolios.length || 
		    expectedCornerPortfolios[0].length != cornerPortfolios[0].length) {
			weightsOK = false;
		}
		else {
			for (var i = 0; i < expectedCornerPortfolios.length; ++i) {
				var expectedWeights = expectedCornerPortfolios[i][0];
				var expectedLambda = expectedCornerPortfolios[i][1];
				
				var weights = cornerPortfolios[i][0].toArray();
				var lambda = cornerPortfolios[i][1];
				
				for (var j = 0; j < weights.length; ++j) {
					if (Math.abs(weights[j] - expectedWeights[j]) > 1e-4) {
						weightsOK = false;
						break;
					}
				}
				if (Math.abs(expectedLambda - lambda) > 1e-4) {
					lambdaOK = false;
					break;
				}
			}
		}
		assert.equal(weightsOK, true, "Mean variance portfolio - Corner portfolios #6, weights");
		assert.equal(lambdaOK, true, "Mean variance portfolio - Corner portfolios #6, lambdas");
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
		
		var expectedCornerPortfolios = [[[0.0012578867048532302, 0, 0, 0, 0.00025227580212938694, 0.015783537737661113, 0.007213947630812098, 0.9754923521245441], 0],									
										[[0.0025252262108544496, 0, 0, 0, 0.001887825588463257, 0.015508854157720443, 0.014692458694259876, 0.9653856353487019], 0.0016518985816731664], 
										[[0.0059458372051487424, 0, 0.011373637171205212, 0, 0.007925805020110858, 0, 0.05260108602780383, 0.9221536345757314], 0.0061468557113804265], 
										[[0.005539024486768344, 0, 0.01333767478516882, 0, 0.010453211204478393, 0, 0.06402244820190232, 0.9066476413216822], 0.007810872279412646],
										[[0, 0.005120342417396367, 0.01542305298430301, 0, 0.012201458745100379, 0, 0.07041370288764186, 0.8968414429655583], 0.009462582988609844],
										[[0, 0.028434275026287065, 0.11159560997120876, 0, 0.13077134305689622, 0, 0.729198771945608, 0],  0.11691392600773665],
										[[0, 0.25, 0.040598510405476254, 0, 0.009674272143084452, 0, 0.6997272174514393, 0], 0.26377183570091994],
										[[0, 0.25, 0.03158049051256295, 0, 0.02373398066197993, 0, 0.6946855288254571, 0], 0.39912101396530514],
										[[0.0462119683316885, 0.25, 0.00852378319201863, 0, 0, 0, 0.6952642484762929, 0], 0.44342871837702763],
										[[0.05663429563220417, 0.25, 0, 0, 0, 0, 0.6933657043677959, 0], 0.4576339467556026],
										[[0.25, 0.25, 0, 0, 0, 0, 0.5, 0], 1.8766716337311202],
										[[0.25,0.25,0,0,0,0,0.5,0], 8.899361909160637],
										[[0.25, 0.25, 0, 0, 0.25, 0, 0.25, 0], 18.562065370001065]].reverse();
						
		// Test that the algorithm is behaving properly
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																						new PortfolioAllocation.Matrix(covMat), 
																						{ constraints: { minWeights: minWeights, maxWeights: maxWeights }}).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #7');
	}

	// Test using data for equal returns, with only one IN asset for the maximum return portfolio
	// https://github.com/lequant40/portfolio_allocation_js/issues/8
	{
		var covMat = [[2.3/100, 2.5/100, 3.3/100, 1.7/100, 2.3/100, -0.4/100, -0.3/100, -0.3/100, 0.1/100, 1.7/100],
					[2.5/100, 2.9/100, 3.8/100, 1.9/100, 2.5/100, -0.4/100, -0.3/100, -0.3/100, 0.2/100, 1.9/100],
					[3.3/100, 3.8/100, 5.3/100, 2.3/100, 3.1/100, -0.5/100, -0.3/100, -0.4/100, 0.3/100, 2.9/100],
					[1.7/100, 1.9/100, 2.3/100, 1.4/100, 1.8/100, -0.3/100, -0.3/100, -0.3/100, 0.0/100, 1.1/100],
					[2.3/100, 2.5/100, 3.1/100, 1.8/100, 4.0/100, -0.1/100, 0.0/100, 0.2/100, 0.7/100, 1.8/100],
					[-0.4/100, -0.4/100, -0.5/100, -0.3/100, -0.1/100, 0.2/100, 0.2/100, 0.2/100, 0.2/100, -0.1/100],
					[-0.3/100, -0.3/100, -0.3/100, -0.3/100, 0.0/100, 0.2/100, 0.4/100, 0.2/100, 0.3/100, 0.0/100],
					[-0.3/100, -0.3/100, -0.4/100, -0.3/100, 0.2/100, 0.2/100, 0.2/100, 0.3/100, 0.2/100, 0.0/100],
					[0.1/100, 0.2/100, 0.3/100, 0.0/100, 0.7/100, 0.2/100, 0.3/100, 0.2/100, 0.5/100, 0.4/100],
					[1.7/100, 1.9/100, 2.9/100, 1.1/100, 1.8/100, -0.1/100, 0.0/100, 0.0/100, 0.4/100, 2.0/100]];
		var returns = [4.7/100, 5.0/100, 6.1/100, 4.0/100, 5.4/100, 2.0/100, 2.2/100, 2.4/100, 2.4/100, 3.6/100];
		var minWeights = [3/100, 3/100, 3/100, 0/100, 0/100, 0/100, 0/100, 3/100, 0/100, 0/100];
		var maxWeights = [35/100, 35/100, 35/100, 35/100, 35/100, 35/100, 35/100, 35/100, 35/100, 35/100];
		
		var expectedCornerPortfolios = [[[0.03,0.24,0.35,0,0.35,0,0,0.03,0,0],1.1957692307692307],
		                                [[0.03,0.04782383419689129,0.35,0,0.35,0,0,0.22217616580310873,0,0],0.9148963730569948],
										[[0.03,0.03695652173912947,0.35,0,0.3208695652173892,0,0,0.26217391304348137,0,0],0.8665217391304311],
										[[0.03567375886524782,0.030000000000000027,0.35,0,0.3147517730496442,0,0,0.269574468085108,0,0],0.8569503546099267],
										[[0.03,0.03,0.35,0,0.2926436781609211,0,0,0.2973563218390789,0,0],0.8231034482758638],
										[[0.03,0.03,0.35,0,0.24000000000000016,0,0,0.35,0,0],0.7546666666666666],
										[[0.03,0.03,0.35,0,0.24,0,0,0.35,0,0],0.7328571428571421],
										[[0.033888888888888635,0.03,0.35,0,0.23611111111111138,0,0,0.35,0,0],0.7234126984126983],
										[[0.030000000000000027,0.03,0.35,0.004999999999999449,0.23500000000000043,0,0,0.35,0,0],0.7214285714285718],
										[[0.03,0.03,0.35,0.007741935483870588,0.2322580645161294,0,0,0.35,0,0],0.7179032258064516],
										[[0.03,0.03,0.35,0.009003808378432322,0.18209861695730617,0,0.048897574664261534,0.35,0,0],0.6499198236119462],
										[[0.03,0.03,0.16565458422174467,0.22871073205401993,0.020810234541574546,0,0.17482444918266096,0.35,0,0],0.38088415067518994],
										[[0.03,0.03,0.14475897435897425,0.25005811965811964,0,0.025025641025640977,0.17015726495726524,0.35,0,0],0.34491623931623916],
										[[0.03,0.03,0.03,0.2657070707070706,0,0.18151515151515218,0.11277777777777728,0.35,0,0],0.1727777777777777],
										[[0.03,0.03,0.03,0.2086363636363638,0,0.3013636363636366,0.0499999999999997,0.35,0,0],0.11000000000000014],
										[[0.03,0.03,0.03,0.20115384615384624,0,0.35,0.04176923076923074,0.3170769230769231,0,0],0.1017692307692308],
										[[0.03,0.03,0.03,0.12632352941176475,0,0.35,0.13455882352941162,0.2991176470588237,0,0],0]];
		
		//
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																						new PortfolioAllocation.Matrix(covMat), 
																						{ constraints: { minWeights: minWeights, maxWeights: maxWeights }}).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #8, equal returns, one IN asset in maximum return portfolio');
	}

	// Test using data for equal returns, with 2 assets IN for the maximum return portfolio, with covariance matrix invertible
	{
		var covMat = [[1, 0.1, 0.2], [0.1, 1, -0.1], [0.2, -0.1, 1]];
		var returns = [2, 2, 1];
		
		var expectedCornerPortfolios = [[[0.5, 0.5, 0],0.5], [[0.261904761904762, 0.38095238095238093, 0.3571428571428571], 0]];
		
		//
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																						new PortfolioAllocation.Matrix(covMat)).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #9, equal returns, 2 IN asset in maximum return portfolio');
	}
	
	// Test using data for equal returns, with 2 assets IN for the maximum return portfolio, with covariance matrix not invertible but first KKT system invertible
	{
		var covMat = [[1, 0.1, 0, 0], [0.1, 1, 0, 0], [0, 0, 0.2, 0.2], [0, 0, 0.2, 0.2]]; // determinant = 0
		var returns = [2, 2, 1, 0.9];

		var expectedCornerPortfolios = [[[0.5, 0.5, 0, 0],0.55], [[0.13333333333333336, 0.13333333333333336, 0.7333333333333334, 0], 0]];
		
		//
		var cornerPortfolios = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(new PortfolioAllocation.Matrix(returns), 
																						new PortfolioAllocation.Matrix(covMat)).cornerPortfolios;
		assert.deepEqual(cornerPortfoliosToArray(cornerPortfolios), expectedCornerPortfolios, 'Mean variance portfolio - Corner portfolios #9, equal returns, 2 IN asset in maximum return portfolio, non invertible covariance matrix');
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
		var expectedPortfolio = [[0.08397318948217397, 0.04890598613711383, 0, 0.21830925954049438, 0.0016773041709360515, 0.18120064671441213, 0, 0.031183038765169584, 0.007859012532850748, 0.42689156265684947], 0.05105260106191383];
										
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
	
	// Test using static data
	// Test the case of a positive semi-definite covariance matrix with a 0 volatility point, as well as a 0 return point, on the efficient frontier
	// Test also that a partial investment constraint provides the same result.
	// The result has been validated with a grid search
	// Reference: Portfolio Selection, H. Markowitz example, chapter II "Illustrative portfolio analyses"
	//
	// Note that this test is sensitive to the value of epsVolatility, because for all values of epsVolatility,
	// a portfolio with the same Sharpe ratio as the optimal portfolio below (SR ~= 0.8520782020723375) can be constructed
	// with a proportion of cash (the last asset), so that there is no test on portfolio weights.
	{		
		var expectedSharpeRatio = 0.8520782020723374;
	
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

		var covMatNoZero = [[0.05338816358,0.02149069753,0.02865533642,0.04896485802,0.01624895062,0.03223945062,0.02425553395,0.03999812963,0.0361509784],
							[0.02149069753,0.01468446914,0.01878391358,0.02441658642,0.008041938272,0.01002193827,0.01448993827,0.02536259259,0.02083593827],
							[0.02865533642,0.01878391358,0.08550016358,0.06260714198,0.04439938272,0.01328671605,0.01043991049,0.06864603704,0.0420215216],
							[0.04896485802,0.02441658642,0.06260714198,0.09546446914,0.05153806173,0.02902461728,0.02077028395,0.09002012963,0.03664589506],
							[0.01624895062,0.008041938272,0.04439938272,0.05153806173,0.1278900988,0.0128384321,0.02091715432,0.1015344074,0.04497232099],
							[0.03223945062,0.01002193827,0.01328671605,0.02902461728,0.0128384321,0.04125832099,0.01127854321,0.02960762963,0.02165332099],
							[0.02425553395,0.01448993827,0.01043991049,0.02077028395,0.02091715432,0.01127854321,0.02883379321,0.02913762963,0.01739445988],
							[0.03999812963,0.02536259259,0.06864603704,0.09002012963,0.1015344074,0.02960762963,0.02913762963,0.1467278889,0.05284057407],
							[0.0361509784,0.02083593827,0.0420215216,0.03664589506,0.04497232099,0.02165332099,0.01739445988,0.05284057407,0.07926979321]];
		var returnsNoZero = [0.06594444444,0.06155555556,0.1460555556,0.1734444444,0.1981111111,0.05511111111,0.1276111111,0.1903333333,0.1156111111];
		var rf = 0;
		var epsVolatility = Math.random() * (0.15 - 1e-4) + 1e-4; // generate a random minimum volatility between 1e-4 and 0.15 (portfolios with target SR are found until ~0.16) 
		
		//
		var efficientFrontierCla = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returns, covMat, {optimizationMethod: 'critical-line'});	
		var maxSharpeRatioPortfolioWeightsCla = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf, {optimizationMethod: 'critical-line', epsVolatility: epsVolatility});
		var srCla = efficientFrontierCla.computePortfolioSharpeRatio(new PortfolioAllocation.Matrix(maxSharpeRatioPortfolioWeightsCla), rf);
		
		var maxSharpeRatioPortfolioWeightsGsmo = PortfolioAllocation.maximumSharpeRatioWeights(returns, covMat, rf, {epsVolatility: epsVolatility});
		var srGsmo = efficientFrontierCla.computePortfolioSharpeRatio(new PortfolioAllocation.Matrix(maxSharpeRatioPortfolioWeightsGsmo), rf);

		assert.equal(Math.abs(srCla - expectedSharpeRatio) <= 1e-14, true, "Test #5, CLA");
		assert.equal(Math.abs(srGsmo - expectedSharpeRatio) <= 1e-08, true, "Test #5, GSMO");

		//
		var efficientFrontierClaNoZero = new PortfolioAllocation.MeanVarianceEfficientFrontierCla(returnsNoZero, covMatNoZero, {optimizationMethod: 'critical-line', constraints: {fullInvestment: false}});	
		var maxSharpeRatioPortfolioWeightsClaNoZero = PortfolioAllocation.maximumSharpeRatioWeights(returnsNoZero, covMatNoZero, rf, {optimizationMethod: 'critical-line', constraints: {fullInvestment: false}});
		var srClaNoZero = efficientFrontierClaNoZero.computePortfolioSharpeRatio(new PortfolioAllocation.Matrix(maxSharpeRatioPortfolioWeightsClaNoZero), rf);
		
		var maxSharpeRatioPortfolioWeightsGsmoNoZero = PortfolioAllocation.maximumSharpeRatioWeights(returnsNoZero, covMatNoZero, rf, {constraints: {fullInvestment: false}});
		var srGsmo = efficientFrontierClaNoZero.computePortfolioSharpeRatio(new PortfolioAllocation.Matrix(maxSharpeRatioPortfolioWeightsGsmoNoZero), rf);
		
		assert.equal(Math.abs(srCla - expectedSharpeRatio) <= 1e-14, true, "Test #5, partial investment, CLA");
		assert.equal(Math.abs(srGsmo - expectedSharpeRatio) <= 1e-08, true, "Test #5, partial investment,  GSMO");
	}
});

