// ------------------------------------------------------------
QUnit.module('Computational geometry internal module', {
  before: function() {
    // 
  }
});



QUnit.test('Geometric center', function(assert) {    
	// Functions used to generate random problems
	function generateRandomDimension(min, max) {
		return Math.floor(Math.random()*(max-min+1) + min);
	}
	
	function generateRandomValue(minVal, maxVal) {
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test with static data
	// Reference: Geospatial Analysis Online http://www.spatialanalysisonline.com/
	{
		// Define the points
		var points = [new PortfolioAllocation.Matrix([3, 0, 0]), 
		              new PortfolioAllocation.Matrix([14, 3, 0]),
					  new PortfolioAllocation.Matrix([10, 4, 0]),
					  new PortfolioAllocation.Matrix([13, 11, 0]),
					  new PortfolioAllocation.Matrix([4, 13, 0]),
					  new PortfolioAllocation.Matrix([0, 8, 0])];

	    // Compute the geometric center associated to the point above
		var y = PortfolioAllocation.geometricCenter_(points);
		
		// Compare it to the expected solution
		var expectedGeometricCenter = new PortfolioAllocation.Matrix([7.33, 6.5, 0]);
		assert.equal(PortfolioAllocation.Matrix.areEqual(y, expectedGeometricCenter, 1e-2), true, 'Geometric center #1');
	}
	
    // Test with random data
	// Limit case: the geometric center of one point must be this point
	{
        // Generate the dimension
		var n = generateRandomDimension(1, 100); 

        // The point
		var x_1 = PortfolioAllocation.Matrix.fill(n, 1, 
                                                  function(i,j) { 
                                                      return generateRandomValue(-1, 1);
                                                  });
         
		// Compute the geometric center associated to the point above
		var y = PortfolioAllocation.geometricCenter_([x_1]);

        // Ensure both are equal
		assert.equal(PortfolioAllocation.Matrix.areEqual(y, x_1), true, 'Geometric center of one point, n = ' + n);
	}
	
	// Test with random data
    //
	// The geometric center of m points x_1, ..., x_m is the point y which 
    // minimizes the sum of the squared Euclidean distances between itself and each point:
    //
    // y = argmin_x in R^n f(x) = sum (||y - x_i||_2^2), i = 1..m
	//
	// So, a necessary and sufficient condition for a point y to be the geometric center
	// of m points x_1, ... , x_m is that gradf(y) = 0, with:
	//
	// gradf(y) = sum (2 * (y - x_i)), i = 1..m
	//          = 2*m*y - 2*sum (x_i), i = 1..m
	//
	// This allows testing using random data.
    {
        // Generate the dimension
		var n = generateRandomDimension(1, 100);

        // Generate the points
        var m = generateRandomDimension(10, 1000);
		var points = new Array(m);
		for (var k = 0; k < m; ++k) {
			points[k] = PortfolioAllocation.Matrix.fill(n, 1, 
                                                        function(i,j) { 
                                                            return generateRandomValue(-1, 1);
                                                        });
		}

        // Compute the geometric center associated to the points above
		var y = PortfolioAllocation.geometricCenter_(points);

        // Generate the gradient function
		var gradf = function(y) {
			return PortfolioAllocation.Matrix.fill(n, 1, 
                                                   function(i,j) { 
                                                       var sum_i = 0;
													   for (var k = 0; k < m; ++k) {
														   sum_i += points[k].getValue(i, 1);
													   }
													   return 2*m*y.getValue(i, 1) - 2*sum_i;
                                                   });
		}

        // Ensure the gradient is null at the computed geometric center
		assert.equal(gradf(y).vectorNorm('infinity') <= 1e-10, true, 'Geometric center of m points, m = ' + m + ', n = ' + n);
    }
});


QUnit.test('Geometric median', function(assert) {    
	// Functions used to validate random problems
	function f(y, points) {
		var m = points.length;
		
		var res = 0;		
		for (var k = 0; k < m; ++k) {
			res += PortfolioAllocation.Matrix.xmy(y, points[k]).vectorNorm('two');;
		}
		return res;
	}
	function modified_gradf(y, points) {
		var n = y.nbRows;
		var m = points.length;
		
		var res = PortfolioAllocation.Matrix.zeros(n, 1);
		
		for (var k = 0; k < m; ++k) {
			// Compute (y - x_k)/||y - x_k||_2 and add it
			// to the currently computed gradient value.
			//
			// In case ||y - x_k||_2 is numerically zero, 
			// replace it with a small positive number.
			var y_m_x_k = PortfolioAllocation.Matrix.xmy(y, points[k]);
			var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');
			res = PortfolioAllocation.Matrix.axpby(1, res, 1/Math.max(y_m_x_k_two_norm, 1e-12), y_m_x_k, res);
		}
		return res;
	}
	
	// Functions used to generate random problems
	function generateRandomDimension(min, max) {
		return Math.floor(Math.random()*(max-min+1) + min);
	}
	function generateRandomValue(minVal, maxVal) {
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	// Test with static data
	// Reference: Geospatial Analysis Online http://www.spatialanalysisonline.com/
	{
		// Define the points
		var points = [new PortfolioAllocation.Matrix([3, 0, 0]), 
		              new PortfolioAllocation.Matrix([14, 3, 0]),
					  new PortfolioAllocation.Matrix([10, 4, 0]),
					  new PortfolioAllocation.Matrix([13, 11, 0]),
					  new PortfolioAllocation.Matrix([4, 13, 0]),
					  new PortfolioAllocation.Matrix([0, 8, 0])];

	    // Compute the geometric center associated to the point above
		var y = PortfolioAllocation.geometricMedian_(points);
		
		// Compare it to the expected solution
		var expectedGeometricMedian = new PortfolioAllocation.Matrix([8.58, 5.61, 0]);
		assert.equal(PortfolioAllocation.Matrix.areEqual(y, expectedGeometricMedian, 1e-2), true, 'Geometric median #1');
	}
	
	// Test with static data
	// Reference: A NOTE ON FERMAT'S PROBLEM, Harold W. KUHN
	{
		// Define the points
		var points = [new PortfolioAllocation.Matrix([-2, 0]), 
		              new PortfolioAllocation.Matrix([-1, 0]),
					  new PortfolioAllocation.Matrix([1, 0]),
					  new PortfolioAllocation.Matrix([2, 0]),
					  new PortfolioAllocation.Matrix([0, 1]),
					  new PortfolioAllocation.Matrix([0, -1])];

	    // Compute the geometric center associated to the point above
		var y = PortfolioAllocation.geometricMedian_(points);
		
		// Compare it to the expected solution
		var expectedGeometricMedian = new PortfolioAllocation.Matrix([0, 0]);
		assert.equal(PortfolioAllocation.Matrix.areEqual(y, expectedGeometricMedian, 1e-12), true, 'Geometric median #2');
	}
	
    // Test with random data
	// Limit case: the geometric median of one point must be this point
	{
        // Generate the dimension
		var n = generateRandomDimension(1, 100); 

        // The point
		var x_1 = PortfolioAllocation.Matrix.fill(n, 1, 
                                                  function(i,j) { 
                                                      return generateRandomValue(-1, 1);
                                                  });
         
		// Compute the geometric median associated to the point above
		var y = PortfolioAllocation.geometricMedian_([x_1]);

        // Ensure both are equal
		assert.equal(PortfolioAllocation.Matrix.areEqual(y, x_1), true, 'Geometric median of one point, n = ' + n);
	}

    // Test with random data
	// Limit case: the geometric median of collinear points must be one of these points
	//
	// As this condition is numerically difficult to test, this test validates that
	// the optimal geometric median function value (which is attained at one of the input
	// points, per construction) is equal to the geometric median function value 
	// at the computed solution.
	{
        // Generate the dimension
		var n = generateRandomDimension(1, 100); 

        // Generate m collinear points: x_k = x + t_k * d, 
		// with x and d in R^n and t_p in R, k = 1..m 
        var m = generateRandomDimension(1, 100);
		var points = new Array(m);
		var x = PortfolioAllocation.Matrix.fill(n, 1, 
                                                function(i,j) { 
                                                    return generateRandomValue(-1, 1);
                                                });
		var d = PortfolioAllocation.Matrix.fill(n, 1, 
                                                function(i,j) { 
                                                    return generateRandomValue(-1, 1);
                                                });
		for (var k = 0; k < m; ++k) {
			var t_k = generateRandomValue(-1, 1);
			points[k] = PortfolioAllocation.Matrix.axpby(1, x, t_k, d);
		}
         
		// Compute the geometric median associated to the point above
		var y = PortfolioAllocation.geometricMedian_(points);
		
        // Ensure that the computed optimal geometric median function value is the same as 
		// the optimal geometric median function value attained at one of the input points.
		var minValF = Number.POSITIVE_INFINITY;
		for (var k = 0; k < m; ++k) {
			var valF = f(points[k], points);
			if (valF < minValF) {
				minValF = valF;
			}
		}		
		assert.equal(Math.abs(minValF - f(y, points)) <= 1e-6, true, 'Geometric median of m collinear points, m = ' + m + ', n = ' + n);
	}

	// Test with random data
	//
	// If the optimal point is not one of the input points, a necessary and sufficient condition 
	// for a point y to be the geometric median of m points x_1, ... , x_m is that the modified
	// gradient of f at y is (numerically) zero: modified_gradf(y) = 0, c.f. A NOTE ON FERMAT'S PROBLEM
	// Harold W. KUHN.
	//
	// This allows testing using random data.
	{
		// Generate the dimension
		var n = generateRandomDimension(10, 100);

        // Generate the points
        var m = generateRandomDimension(100, 10000);
		var points = new Array(m);
		for (var k = 0; k < m; ++k) {
			points[k] = PortfolioAllocation.Matrix.fill(n, 1, 
                                                        function(i,j) { 
                                                            return generateRandomValue(-1, 1);
                                                        });
		}


        // Compute the geometric median associated to the points above
		var y = PortfolioAllocation.geometricMedian_(points);
		
		// Ensure that the modified gradient of f is null at the
		// computed point.
		assert.equal(modified_gradf(y, points).vectorNorm("infinity") <= 1e-4, true, 'Geometric median of m points, m = ' + m + ', n = ' + n);
	}
	
	// Test with static data
	// Reference: TSPLIB P654 Drilling problem (Reinelt), https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/ 
	{
		// The static data
		var drillingPoints = [[1245.00,1255.00],[1807.50,2192.50],[1822.50,2207.50],[1867.50,2732.50],[2377.50,2612.50],[2467.50,3587.50],[2497.50,3587.50],[2572.50,3677.50],[2647.50,3782.50],[2662.50,3977.50],[3307.50,3377.50],[3337.50,3347.50],[3322.50,3632.50],[3712.50,3677.50],[4267.50,2372.50],[5655.00,1255.00],[4852.50,2507.50],[4882.50,3497.50],[4852.50,3617.50],[4402.50,4067.50],[4417.50,4052.50],[4237.50,4022.50],[4267.50,4007.50],[4402.50,4007.50],[4357.50,3977.50],[4417.50,3962.50],[4357.50,3947.50],[4312.50,3932.50],[4252.50,3932.50],[4222.50,3902.50],[4327.50,3887.50],[4297.50,3872.50],[4297.50,3842.50],[4327.50,3827.50],[4207.50,3827.50],[4222.50,3767.50],[4252.50,3767.50],[4357.50,3752.50],[4357.50,3722.50],[4237.50,3722.50],[4267.50,3707.50],[4312.50,3692.50],[4252.50,3662.50],[4312.50,3647.50],[4297.50,3632.50],[4342.50,3617.50],[4222.50,3602.50],[4327.50,3587.50],[4327.50,3527.50],[4387.50,3467.50],[4282.50,3467.50],[4207.50,3467.50],[4342.50,3437.50],[4327.50,3422.50],[4327.50,3392.50],[4342.50,3377.50],[4237.50,3362.50],[4387.50,3317.50],[4312.50,3272.50],[4222.50,3257.50],[4207.50,3227.50],[3697.50,3482.50],[2692.50,4667.50],[2647.50,4682.50],[3067.50,4922.50],[3127.50,5042.50],[3082.50,5087.50],[3007.50,5117.50],[3142.50,5222.50],[2707.50,5117.50],[2692.50,5057.50],[2572.50,5027.50],[2377.50,5132.50],[1477.50,4022.50],[1237.50,3992.50],[1357.50,3962.50],[1417.50,3932.50],[1252.50,3932.50],[1312.50,3902.50],[1282.50,3827.50],[1372.50,3827.50],[1507.50,3827.50],[1267.50,3812.50],[1282.50,3797.50],[1267.50,3767.50],[1432.50,3767.50],[1372.50,3737.50],[1327.50,3722.50],[1267.50,3707.50],[1492.50,3662.50],[1267.50,3632.50],[1357.50,3602.50],[1507.50,3587.50],[1267.50,3572.50],[1387.50,3527.50],[1342.50,3512.50],[1327.50,3452.50],[1342.50,3377.50],[1327.50,3317.50],[1432.50,3272.50],[1462.50,3257.50],[1402.50,3227.50],[1807.50,4712.50],[1822.50,4727.50],[1867.50,5252.50],[1245.00,5845.00],[1042.50,2462.50],[1057.50,2462.50],[1072.50,2462.50],[1087.50,2462.50],[1042.50,2477.50],[1057.50,2477.50],[1072.50,2477.50],[1087.50,2477.50],[1042.50,2492.50],[1057.50,2492.50],[1072.50,2492.50],[1087.50,2492.50],[1042.50,2507.50],[1057.50,2507.50],[1072.50,2507.50],[1087.50,2507.50],[1042.50,2522.50],[1057.50,2522.50],[1072.50,2522.50],[1087.50,2522.50],[1042.50,2537.50],[1057.50,2537.50],[1072.50,2537.50],[1087.50,2537.50],[1042.50,2552.50],[1057.50,2552.50],[1072.50,2552.50],[1087.50,2552.50],[1042.50,2567.50],[1057.50,2567.50],[1072.50,2567.50],[1087.50,2567.50],[1042.50,2582.50],[1057.50,2582.50],[1072.50,2582.50],[1087.50,2582.50],[1042.50,2597.50],[1057.50,2597.50],[1072.50,2597.50],[1087.50,2597.50],[1042.50,2612.50],[1057.50,2612.50],[1072.50,2612.50],[1087.50,2612.50],[1042.50,2627.50],[1057.50,2627.50],[1072.50,2627.50],[1087.50,2627.50],[1042.50,2642.50],[1057.50,2642.50],[1072.50,2642.50],[1087.50,2642.50],[1042.50,2657.50],[1057.50,2657.50],[1072.50,2657.50],[1087.50,2657.50],[1042.50,2672.50],[1057.50,2672.50],[1072.50,2672.50],[1087.50,2672.50],[1042.50,2687.50],[1057.50,2687.50],[1072.50,2687.50],[1087.50,2687.50],[1042.50,2702.50],[1057.50,2702.50],[1072.50,2702.50],[1087.50,2702.50],[1042.50,2717.50],[1057.50,2717.50],[1072.50,2717.50],[1087.50,2717.50],[1042.50,2732.50],[1057.50,2732.50],[1072.50,2732.50],[1087.50,2732.50],[1042.50,2747.50],[1057.50,2747.50],[1072.50,2747.50],[1087.50,2747.50],[1042.50,2762.50],[1057.50,2762.50],[1072.50,2762.50],[1087.50,2762.50],[1042.50,2777.50],[1057.50,2777.50],[1072.50,2777.50],[1087.50,2777.50],[1042.50,2792.50],[1057.50,2792.50],[1072.50,2792.50],[1087.50,2792.50],[1042.50,2807.50],[1057.50,2807.50],[1072.50,2807.50],[1087.50,2807.50],[1042.50,4952.50],[1057.50,4952.50],[1072.50,4952.50],[1087.50,4952.50],[1042.50,4967.50],[1057.50,4967.50],[1072.50,4967.50],[1087.50,4967.50],[1042.50,4982.50],[1057.50,4982.50],[1072.50,4982.50],[1087.50,4982.50],[1042.50,4997.50],[1057.50,4997.50],[1072.50,4997.50],[1087.50,4997.50],[1042.50,5012.50],[1057.50,5012.50],[1072.50,5012.50],[1087.50,5012.50],[1042.50,5027.50],[1057.50,5027.50],[1072.50,5027.50],[1087.50,5027.50],[1042.50,5042.50],[1057.50,5042.50],[1072.50,5042.50],[1087.50,5042.50],[1042.50,5057.50],[1057.50,5057.50],[1072.50,5057.50],[1087.50,5057.50],[1042.50,5072.50],[1057.50,5072.50],[1072.50,5072.50],[1087.50,5072.50],[1042.50,5087.50],[1057.50,5087.50],[1072.50,5087.50],[1087.50,5087.50],[1042.50,5102.50],[1057.50,5102.50],[1072.50,5102.50],[1087.50,5102.50],[1042.50,5117.50],[1057.50,5117.50],[1072.50,5117.50],[1087.50,5117.50],[1042.50,5132.50],[1057.50,5132.50],[1072.50,5132.50],[1087.50,5132.50],[1042.50,5147.50],[1057.50,5147.50],[1072.50,5147.50],[1087.50,5147.50],[1042.50,5162.50],[1057.50,5162.50],[1072.50,5162.50],[1087.50,5162.50],[1042.50,5177.50],[1057.50,5177.50],[1072.50,5177.50],[1087.50,5177.50],[1042.50,5192.50],[1057.50,5192.50],[1072.50,5192.50],[1087.50,5192.50],[1042.50,5207.50],[1057.50,5207.50],[1072.50,5207.50],[1087.50,5207.50],[1042.50,5222.50],[1057.50,5222.50],[1072.50,5222.50],[1087.50,5222.50],[1042.50,5237.50],[1057.50,5237.50],[1072.50,5237.50],[1087.50,5237.50],[1042.50,5252.50],[1057.50,5252.50],[1072.50,5252.50],[1087.50,5252.50],[1042.50,5267.50],[1057.50,5267.50],[1072.50,5267.50],[1087.50,5267.50],[1042.50,5282.50],[1057.50,5282.50],[1072.50,5282.50],[1087.50,5282.50],[1042.50,5297.50],[1057.50,5297.50],[1072.50,5297.50],[1087.50,5297.50],[5812.50,2462.50],[5827.50,2462.50],[5842.50,2462.50],[5857.50,2462.50],[5812.50,2477.50],[5827.50,2477.50],[5842.50,2477.50],[5857.50,2477.50],[5812.50,2492.50],[5827.50,2492.50],[5842.50,2492.50],[5857.50,2492.50],[5812.50,2507.50],[5827.50,2507.50],[5842.50,2507.50],[5857.50,2507.50],[5812.50,2522.50],[5827.50,2522.50],[5842.50,2522.50],[5857.50,2522.50],[5812.50,2537.50],[5827.50,2537.50],[5842.50,2537.50],[5857.50,2537.50],[5812.50,2552.50],[5827.50,2552.50],[5842.50,2552.50],[5857.50,2552.50],[5812.50,2567.50],[5827.50,2567.50],[5842.50,2567.50],[5857.50,2567.50],[5812.50,2582.50],[5827.50,2582.50],[5842.50,2582.50],[5857.50,2582.50],[5812.50,2597.50],[5827.50,2597.50],[5842.50,2597.50],[5857.50,2597.50],[5812.50,2612.50],[5827.50,2612.50],[5842.50,2612.50],[5857.50,2612.50],[5812.50,2627.50],[5827.50,2627.50],[5842.50,2627.50],[5857.50,2627.50],[5812.50,2642.50],[5827.50,2642.50],[5842.50,2642.50],[5857.50,2642.50],[5812.50,2657.50],[5827.50,2657.50],[5842.50,2657.50],[5857.50,2657.50],[5812.50,2672.50],[5827.50,2672.50],[5842.50,2672.50],[5857.50,2672.50],[5812.50,2687.50],[5827.50,2687.50],[5842.50,2687.50],[5857.50,2687.50],[5812.50,2702.50],[5827.50,2702.50],[5842.50,2702.50],[5857.50,2702.50],[5812.50,2717.50],[5827.50,2717.50],[5842.50,2717.50],[5857.50,2717.50],[5812.50,2732.50],[5827.50,2732.50],[5842.50,2732.50],[5857.50,2732.50],[5812.50,2747.50],[5827.50,2747.50],[5842.50,2747.50],[5857.50,2747.50],[5812.50,2762.50],[5827.50,2762.50],[5842.50,2762.50],[5857.50,2762.50],[5812.50,2777.50],[5827.50,2777.50],[5842.50,2777.50],[5857.50,2777.50],[5812.50,2792.50],[5827.50,2792.50],[5842.50,2792.50],[5857.50,2792.50],[5812.50,2807.50],[5827.50,2807.50],[5842.50,2807.50],[5857.50,2807.50],[5812.50,4952.50],[5827.50,4952.50],[5842.50,4952.50],[5857.50,4952.50],[5812.50,4967.50],[5827.50,4967.50],[5842.50,4967.50],[5857.50,4967.50],[5812.50,4982.50],[5827.50,4982.50],[5842.50,4982.50],[5857.50,4982.50],[5812.50,4997.50],[5827.50,4997.50],[5842.50,4997.50],[5857.50,4997.50],[5812.50,5012.50],[5827.50,5012.50],[5842.50,5012.50],[5857.50,5012.50],[5812.50,5027.50],[5827.50,5027.50],[5842.50,5027.50],[5857.50,5027.50],[5812.50,5042.50],[5827.50,5042.50],[5842.50,5042.50],[5857.50,5042.50],[5812.50,5057.50],[5827.50,5057.50],[5842.50,5057.50],[5857.50,5057.50],[5812.50,5072.50],[5827.50,5072.50],[5842.50,5072.50],[5857.50,5072.50],[5812.50,5087.50],[5827.50,5087.50],[5842.50,5087.50],[5857.50,5087.50],[5812.50,5102.50],[5827.50,5102.50],[5842.50,5102.50],[5857.50,5102.50],[5812.50,5117.50],[5827.50,5117.50],[5842.50,5117.50],[5857.50,5117.50],[5812.50,5132.50],[5827.50,5132.50],[5842.50,5132.50],[5857.50,5132.50],[5812.50,5147.50],[5827.50,5147.50],[5842.50,5147.50],[5857.50,5147.50],[5812.50,5162.50],[5827.50,5162.50],[5842.50,5162.50],[5857.50,5162.50],[5812.50,5177.50],[5827.50,5177.50],[5842.50,5177.50],[5857.50,5177.50],[5812.50,5192.50],[5827.50,5192.50],[5842.50,5192.50],[5857.50,5192.50],[5812.50,5207.50],[5827.50,5207.50],[5842.50,5207.50],[5857.50,5207.50],[5812.50,5222.50],[5827.50,5222.50],[5842.50,5222.50],[5857.50,5222.50],[5812.50,5237.50],[5827.50,5237.50],[5842.50,5237.50],[5857.50,5237.50],[5812.50,5252.50],[5827.50,5252.50],[5842.50,5252.50],[5857.50,5252.50],[5812.50,5267.50],[5827.50,5267.50],[5842.50,5267.50],[5857.50,5267.50],[5812.50,5282.50],[5827.50,5282.50],[5842.50,5282.50],[5857.50,5282.50],[5812.50,5297.50],[5827.50,5297.50],[5842.50,5297.50],[5857.50,5297.50],[1042.50,1802.50],[1042.50,1817.50],[1042.50,1832.50],[1042.50,1877.50],[1042.50,1892.50],[1042.50,1922.50],[1042.50,1952.50],[1042.50,1997.50],[1042.50,2072.50],[1042.50,2117.50],[1042.50,2162.50],[1042.50,2222.50],[1042.50,2252.50],[1042.50,2357.50],[1042.50,2402.50],[1042.50,4292.50],[1042.50,4307.50],[1042.50,4322.50],[1042.50,4367.50],[1042.50,4382.50],[1042.50,4412.50],[1042.50,4442.50],[1042.50,4487.50],[1042.50,4562.50],[1042.50,4607.50],[1042.50,4652.50],[1042.50,4712.50],[1042.50,4742.50],[1042.50,4847.50],[1042.50,4892.50],[1057.50,1847.50],[1057.50,1952.50],[1057.50,2027.50],[1057.50,2087.50],[1057.50,2147.50],[1057.50,2282.50],[1057.50,4337.50],[1057.50,4442.50],[1057.50,4517.50],[1057.50,4577.50],[1057.50,4637.50],[1057.50,4772.50],[1072.50,1817.50],[1072.50,1832.50],[1072.50,1892.50],[1072.50,1967.50],[1072.50,2162.50],[1072.50,2222.50],[1072.50,2282.50],[1072.50,2357.50],[1072.50,4307.50],[1072.50,4322.50],[1072.50,4382.50],[1072.50,4457.50],[1072.50,4652.50],[1072.50,4712.50],[1072.50,4772.50],[1072.50,4847.50],[1087.50,1787.50],[1087.50,1847.50],[1087.50,1877.50],[1087.50,1967.50],[1087.50,1997.50],[1087.50,2027.50],[1087.50,2087.50],[1087.50,2117.50],[1087.50,2147.50],[1087.50,2192.50],[1087.50,2252.50],[1087.50,2402.50],[1087.50,4277.50],[1087.50,4337.50],[1087.50,4367.50],[1087.50,4457.50],[1087.50,4487.50],[1087.50,4517.50],[1087.50,4577.50],[1087.50,4607.50],[1087.50,4637.50],[1087.50,4682.50],[1087.50,4742.50],[1087.50,4892.50],[5812.50,1802.50],[5812.50,1817.50],[5812.50,1832.50],[5812.50,1877.50],[5812.50,1892.50],[5812.50,1922.50],[5812.50,1952.50],[5812.50,1997.50],[5812.50,2072.50],[5812.50,2117.50],[5812.50,2162.50],[5812.50,2222.50],[5812.50,2252.50],[5812.50,2357.50],[5812.50,2402.50],[5812.50,4292.50],[5812.50,4307.50],[5812.50,4322.50],[5812.50,4367.50],[5812.50,4382.50],[5812.50,4412.50],[5812.50,4442.50],[5812.50,4487.50],[5812.50,4562.50],[5812.50,4607.50],[5812.50,4652.50],[5812.50,4712.50],[5812.50,4742.50],[5812.50,4847.50],[5812.50,4892.50],[5827.50,1847.50],[5827.50,1952.50],[5827.50,2027.50],[5827.50,2087.50],[5827.50,2147.50],[5827.50,2282.50],[5827.50,4337.50],[5827.50,4442.50],[5827.50,4517.50],[5827.50,4577.50],[5827.50,4637.50],[5827.50,4772.50],[5842.50,1817.50],[5842.50,1832.50],[5842.50,1892.50],[5842.50,1967.50],[5842.50,2162.50],[5842.50,2222.50],[5842.50,2282.50],[5842.50,2357.50],[5842.50,4307.50],[5842.50,4322.50],[5842.50,4382.50],[5842.50,4457.50],[5842.50,4652.50],[5842.50,4712.50],[5842.50,4772.50],[5842.50,4847.50],[5857.50,1787.50],[5857.50,1847.50],[5857.50,1877.50],[5857.50,1967.50],[5857.50,1997.50],[5857.50,2027.50],[5857.50,2087.50],[5857.50,2117.50],[5857.50,2147.50],[5857.50,2192.50],[5857.50,2252.50],[5857.50,2402.50],[5857.50,4277.50],[5857.50,4337.50],[5857.50,4367.50],[5857.50,4457.50],[5857.50,4487.50],[5857.50,4517.50],[5857.50,4577.50],[5857.50,4607.50],[5857.50,4637.50],[5857.50,4682.50],[5857.50,4742.50],[5857.50,4892.50]];
		
		// Generate the points
		var n = 2;	
		var m = drillingPoints.length;
		var points = new Array(m);
		for (var k = 0; k < m; ++k) {
			points[k] = new PortfolioAllocation.Matrix(drillingPoints[k]);
		}
		
        // Compute the geometric median associated to the points above
		var y = PortfolioAllocation.geometricMedian_(points);
		
		// Ensure that the modified gradient of f is null at the
		// computed point.
		assert.equal(modified_gradf(y, points).vectorNorm("infinity") <= 1e-4, true, 'Geometric median, TSPLIB P654');
	}
});
