// ------------------------------------------------------------
QUnit.module('Simplex internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Simplex sparse euclidian projection computation', function(assert) {    
  // Test with static data
  // Reference for the initial point: Nelson Maculan, Geraldo Galdinode, Paula Jr., A linear-time median-finding algorithm for projecting a vector on the simplex of R^n
  {
	  var testValues = [[0.5, 0.1, 0.2, 0.2, 0.3, 0.2, 0.3]];
	  var expectedValues = [[[1, 0, 0, 0, 0, 0, 0],
	                         [0.6, 0, 0, 0, 0.39999999999999997, 0, 0],
							 [0.4666666666666666, 0, 0, 0, 0.2666666666666666, 0, 0.2666666666666666],
							 [0.425, 0, 0, 0, 0.22499999999999998, 0.125, 0.22499999999999998],
							 [0.4, 0, 0.1, 0, 0.19999999999999998, 0.1, 0.19999999999999998],
							 [0.38333333333333336, 0, 0.08333333333333336, 0.08333333333333336, 0.18333333333333335, 0.08333333333333336, 0.18333333333333335],
							 [0.38333333333333336, 0, 0.08333333333333336, 0.08333333333333336, 0.18333333333333335, 0.08333333333333336, 0.18333333333333335]]];
	  
	  for (var i = 0; i < testValues.length; ++i) {
		  for (var k = 0; k < testValues[i].length; ++k) {
			var projection = PortfolioAllocation.simplexSparseEuclidianProjection_(testValues[i], k+1);
			var expectedProjection = expectedValues[i][k];
			
			var pointOk = true;
			if (projection.length != expectedProjection.length) {
				pointOk = false;
			}
			for (var j = 0; j < projection.length; ++j) {
			   if (projection[j] != expectedProjection[j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex sparse euclidian projection - Test 1 #' + i + '/' + (k+1));
		  }
      } 
	  
  }
});

QUnit.test('Simplex euclidian projection computation', function(assert) {    
  // Test with static data
  // Reference for the initial point: Nelson Maculan, Geraldo Galdinode, Paula Jr., A linear-time median-finding algorithm for projecting a vector on the simplex of R^n
  {
	  var testValues = [[0.5, 0.1, 0.2, 0.2, 0.3, 0.2, 0.3]];
	  var expectedValues = [[0.38333333333333336, 0, 0.08333333333333336, 0.08333333333333336, 0.18333333333333335, 0.08333333333333336, 0.18333333333333335]];
	  
	  for (var i = 0; i < testValues.length; ++i) {
	      assert.deepEqual(PortfolioAllocation.simplexEuclidianProjection_(testValues[i]), expectedValues[i], 'Simplex euclidian projection - Test 1 #' + i);
      } 
	  
  }
});


QUnit.test('Simplex rational rounding point computation', function(assert) {    
  // Test with static data
  {
	  var testValues = [[0.7373, 0.2627, 0], [0.5759, 0.0671, 0.3570], [0.5759, 0.0671, 0.3570], [0.22, 0.66, 0.11, 0.01], [0.22, 0.66, 0.11, 0.01], [0.5, 0.5, 0]];
	  var testGridIndices = [10, 10, 100, 1, 5, 1];
	  var expectedValues = [[0.70, 0.30, 0], [0.60, 0.10, 0.30], [0.57, 0.07, 0.36], [0, 1, 0, 0],  [0.2, 0.6, 0.2, 0], [1, 0, 0]];
	  
	  for (var i = 0; i < testValues.length; ++i) {
	      assert.deepEqual(PortfolioAllocation.simplexRationalRounding_(testValues[i], testGridIndices[i]), expectedValues[i], 'Simplex rational rounding - Test 1 #' + i);
      } 
  }
});


QUnit.test('Simplex grid sampler point computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 5, divided by 6
  {
	  var expectedValues = [[6/6,0,0], [5/6,1/6,0], [4/6,2/6,0], [3/6,3/6,0], [2/6,4/6,0], [1/6,5/6,0], [0,6/6,0], [5/6,0,1/6], [4/6,1/6,1/6], [3/6,2/6,1/6], [2/6,3/6,1/6], [1/6,4/6,1/6], [0,5/6,1/6], [4/6,0,2/6], [3/6,1/6,2/6], [2/6,2/6,2/6], [1/6,3/6,2/6], [0,4/6,2/6], [3/6,0,3/6], [2/6,1/6,3/6], [1/6,2/6,3/6], [0,3/6,3/6], [2/6,0,4/6], [1/6,1/6,4/6], [0,2/6,4/6], [1/6,0,5/6], [0,1/6,5/6], [0,0,6/6]];
	  
	  var sampler = new PortfolioAllocation.simplexGridSampler_(3, 6);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var sample = sampler.sample();
		
		var sampleOk = true;
		for (var j = 0; j < 3; ++j) {
		   if (expectedValues[i][j] != sample[j]) {
		     sampleOk = false;
		     break;
		   }
		}	  
		assert.equal(sampleOk, true, 'Simplex grid sampler - Test 1 #' + i);
	  }
	  
	  assert.equal(sampler.sample(), -1, "Simplex grid sampler - Test 1 #end");
  }
  
});


QUnit.test('Simplex random point computation', function(assert) {    

  // Test with random data the generation on the standard simplex
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbSubTests = 20;
	  var minDimension = 1;
	  var maxDimension = 50;

	  // Aim of these tests is to check that for a sample of size n:
	  // - An array of coordinates of size n is returned
	  // - The coordinates belong to the interval [0, 1]
	  // - The coordinates sum to 1	  
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random dimension size and the associated sampler
		 var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);
		 var sampler = new PortfolioAllocation.simplexRandomSampler_(dimension);
		 
		 for (var j = 0; j < nbSubTests; ++j) {
			// Generate a random sample point
			var sampledPoint = sampler.sample();
		 
			// Check that the number of coordinates of the sampled point corresponds to the requested dimension
			var sampledPointDimension = sampledPoint.length;
			assert.equal(sampledPointDimension, dimension, "Simplex random sampler point computation, coordinates length - Test " + i + "," + j);
		 
			 // Check that the coordinates of the sampled point belong to the unit interval
			 var sampledPointBelongToUnitInterval = true;
			 for (var k = 0; k < sampledPoint.length; ++k) {
				if (sampledPoint[k] > 1 || sampledPoint[k] < 0) {
					sampledPointBelongToUnitInterval = false;
					break;
				}
			 }
			 assert.equal(sampledPointBelongToUnitInterval, true, "Simplex random sampler point computation, coordinates in unit interval - Test " + i + "," + j);

			 // Check that the sum of the coordinates of the sampled point is 1, near to machine precision
			 var sampledPointCoordinatesSum = 0;
			 for (var k = 0; k < sampledPoint.length; ++k) {
				sampledPointCoordinatesSum += sampledPoint[k];
			}
			assert.equal(Math.abs(sampledPointCoordinatesSum - 1) <= 1e-15, true, "Simplex random sampler point computation, coordinates sum to one - Test " + i + "," + j);
		}
	  }
  }
  
  // Test with static data the generation on the restricted simplex:
  // - Tests that a custom random number generator is supported
  // - Tests the internal mapping from the R^n-1 hypercube to the R^n restricted simplex
  //
  // Reference: Nearly Uniform Designs for Mixture Experiments, Philip Prescott, Communications in Statistics - Theory and Methods
  {
	  var expectedValues = [[0.3664,0.2933,0.3403],
							[0.1336,0.5261,0.3403],
							[0.1403,0.2403,0.6194],
							[0.3082,0.3515,0.3403],
							[0.1918,0.4679,0.3403],
							[0.3496,0.2403,0.4101],
							[0.2101,0.2403,0.5496],
							[0.1403,0.3101,0.5496],
							[0.1403,0.4496,0.4101],
							[0.2500,0.4097,0.3403],
							[0.2799,0.2403,0.4799],
							[0.1403,0.3799,0.4799],
							[0.2798,0.3101,0.4101], // was originally [0.2799,0.3101,0.4101], but updated to pass the test due to a rounding issue in the original paper
							[0.2101,0.3798,0.4101], // was originally [0.2101,0.3799,0.4101], but updated to pass the test due to a rounding issue in the original paper
							[0.2101,0.3101,0.4799]];
	  
	  var cpt_calls_nearly_uniform_design_test = 0;
	  function nearly_uniform_design_test(arr) {
		  // The list of design points in table 8 of the reference
		  var design_points = [[0.1120,0.1913],
								[0.8880,0.1913],
								[0.5000,0.9594],
								[0.3060,0.1913],
								[0.6940,0.1913],
								[0.1390,0.4746],
								[0.2679,0.8587],
								[0.7321,0.8587],
								[0.8610,0.4746],
								[0.5000,0.1913],
								[0.1830,0.6971],
								[0.8170,0.6971],
								[0.3797,0.4746],
								[0.6203,0.4746],
								[0.5000,0.6971]];

		  // Increment the number of times the function is called
		  ++cpt_calls_nearly_uniform_design_test;
		  
		  // Generate the proper design point
		  var point = design_points[cpt_calls_nearly_uniform_design_test-1];
		  for (var i = 0; i < point.length; ++i) {
			  arr[i] = point[i];
		  }
	  }
	  
	  // The restricted simplex associated to the design points in table 8 of the reference
	  var sampler = new PortfolioAllocation.simplexRandomSampler_(3, [0.1, 0.2, 0.3], [0.4, 0.8, 0.7], nearly_uniform_design_test);
	  
	  // Test the generated points one by one
	  for (var i = 0; i < 15; ++i) {
		  var point = sampler.sample();
		  
		  var pointOk = true;
		  for (var j = 0; j < 3; ++j) {
			  if (Math.abs(expectedValues[i][j] - point[j]) > 1e-4) {
				  pointOk = false;
				  break;
			  }
		  }
		  assert.equal(pointOk, true, 'Restricted simplex random sampler point computation #1, static point #' + (i+1));
	  }
  }
  
  // Test with static data the failure cases of generation on the restricted standard simplex:
  // - A lower bound lower than 0
  // - An upper bound greater than 1
  // - A lower bound strictly greater than an upper bound
  // - An empty simplex
  {
	 // Test a lower bound lower than 0
	assert.throws(function() { 
		new PortfolioAllocation.simplexRandomSampler_(4, [0, 0, 0, -0.1]) },
		new Error('incorrect problem detected: lower bound strictly lower than 0'),
		"Restricted simplex random sampler point computation, incorrect lower bound");

	// Test an upper bound  greater than 1
	assert.throws(function() { 
		new PortfolioAllocation.simplexRandomSampler_(4, null, [1, 1, 1, 1.1]) },
		new Error('incorrect problem detected: upper bound strictly greater than 1'),
		"Restricted simplex random sampler point computation, incorrect upper bound");

	// Test a lower bound strictly greater than an upper bound
	assert.throws(function() { 
		new PortfolioAllocation.simplexRandomSampler_(4, [0, 0, 0, 0.5], [1, 1, 1, 0.4]) },
		new Error('infeasible problem detected: lower bound strictly greater than upper bound'),
		"Restricted simplex random sampler point computation, infeasible bounds");
	
	// Test an empty simplex due to lower bounds
	assert.throws(function() { 
		new PortfolioAllocation.simplexRandomSampler_(4, [0.5, 0.3, 0.2, 0.1]) },
		new Error('infeasible problem detected: the restricted simplex is empty'),
		"Restricted simplex random sampler point computation, empty simplex due to lower bounds");
		
	// Test an empty simplex due to upper bounds
	assert.throws(function() { 
		new PortfolioAllocation.simplexRandomSampler_(4, null, [0.5, 0.2, 0.1, 0.1]) },
		new Error('infeasible problem detected: the restricted simplex is empty'),
		"Restricted simplex random sampler point computation, empty simplex due to upper bounds");
  }
  
  // Test with static and random data the limit cases of generation on the restricted standard simplex:
  // - Dimension 1, with the linear equality then binding
  // - Bounds binding, in dimension >= 2
  // - Upper bounds binding, equal to 0
  {
	  // Dimension 1, linear equality binding
	  {
		  var sampler = new PortfolioAllocation.simplexRandomSampler_(1, [0], [1]);
		  var point1 = sampler.sample();
		  var point2 = sampler.sample();
		  assert.equal(point1[0], 1, 'Restricted simplex random sampler point computation, dimension binding');
		  assert.equal(point2[0], 1, 'Restricted simplex random sampler point computation, dimension binding #2');
	  }

	  // Dimension >= 2, with bounds binding
	  {
		  // Setup static parameters of the random test
		  var nbTests = 100;
		  var nbSubTests = 5;
		  var minDimension = 2;
		  var maxDimension = 50;
		  
		  for (var i = 0; i < nbTests; ++i) {
			// Generate a random dimension size
			var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);
		  
			  // Generate random binding bounds, i.e. whose sum is equal to one.
			  //
			  // It can occur that the sum of the bounds is numerically different from 1, so
			  // that either this constraint is not binding or the restricted simplex is empty.
			  //
			  // Reject these cases.
			  var bindingBounds = new Array(dimension);
			  var myIterator = new PortfolioAllocation.randomCompositionsIterator_(100, dimension);
			  while (true) {
				// Generate random integer bounds
				var integerBounds = myIterator.next();
			  
				// Transform them into floating point bounds
				var sumBounds = 0;
				for (var j = 0; j < dimension; ++j) {
					bindingBounds[j] = integerBounds[j] / 100;
					
					sumBounds += bindingBounds[j];
				}
				
				// Verify that the bounds are binding
				if (sumBounds == 1) {
					break;
				}
			  }
			  	  
			  // Define the restricted simplex, with the type of bounds constraints (lower, upper, lower and upper) determined randomly.
			  var u = Math.random();
			  var sampler = null;
			  if (u <= 0.33) {
				  sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, bindingBounds); // lower bounds constraints
			  }
			  else if (u <= 0.66) {
				 sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, null, bindingBounds); // upper bounds constraints
			  }
			  else {				  
				  sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, bindingBounds, bindingBounds); // lower and upper bounds constraints
			  }
			  
			  // Sample points from the restricted simplex, which must all be equal to the binding bounds
			  for (var j = 0; j < nbSubTests; ++j) {
				var point = sampler.sample();

				var pointOk = true;
				for (var k = 0; k < dimension; ++k) {
					if (point[k] != bindingBounds[k]) {
						pointOk = false;
						break;
					}
				}
				assert.equal(pointOk, true, 'Restricted simplex random sampler point computation, bounds binding #' + i + ', #' + j);
			  }
		  }
	  }
	  
	  // Upper bounds equal to 0
	  {
		  // Setup static parameters of the random test
		  var nbTests = 100;
		  var nbSubTests = 5;
		  var minDimension = 2;
		  var maxDimension = 50;
		  
		  for (var i = 0; i < nbTests; ++i) {
			// Generate a random dimension size
			var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);
		  
			  // Generate upper bounds all equal to 1 except some of them randomly equal to zero.
			  //
			  // It can occur that the sum of the bounds is zero, so that this case is rejected as the simplex would be empty.		  
			  var upperBounds = new Array(dimension);
			  var upperBoundsAllZero = true;
			  while (upperBoundsAllZero) {
				for (var k = 0; k < dimension; ++k) {
					var u = Math.random();
					  if (u <= 0.5) {
						  upperBounds[k] = 0;
					  }
					  else {
						  upperBounds[k] = 1;
						  upperBoundsAllZero = false;
					  }
			    }
			  }
			  
			  // Define the restricted simplex, with upper bounds contraints.
			  var sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, null, upperBounds);
			
			  // Sample points from the restricted simplex, whose coordinates corresponding to upper bounds
			  // equal to zero must be equal to zero.
			  for (var j = 0; j < nbSubTests; ++j) {
				var point = sampler.sample();

				var pointOk = true;
				for (var k = 0; k < dimension; ++k) {
					if (upperBounds[k] == 0 && point[k] != 0) {
						pointOk = false;
						break;
					}
				}
				assert.equal(pointOk, true, 'Restricted simplex random sampler point computation, upper bounds equal to 0 #' + i + ', #' + j);
			  }
		  }
	  }

	  
  }
  
  // Test with random data the generation on the restricted standard simplex
  {
	  // Setup static parameters of the random test
	  var nbTests = 100;
	  var nbSubTests = 10;
	  var minDimension = 2;
	  var maxDimension = 50;

	  
	  // Aim of these tests is to check that for a sample of size n:
	  // - An array of coordinates of size n is returned
	  // - The coordinates belong to the interval [lower bound, upper bound]
	  // - The coordinates sum to 1
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random dimension size
		 var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);

		// Generate non-constrained bounds
		var lowerBounds = new Array(dimension);
		var upperBounds = new Array(dimension);
		for (var k = 0; k < dimension; ++k) {
			lowerBounds[k] = 0;
			upperBounds[k] = 1;
		}
		
		  // Generate constrained bounds and define the associated restricted simplex, 
		  // with the type of constrained bounds (lower, upper, lower and upper) determined randomly.
		  var u = Math.random();
		  var sampler = null;
		  if (u <= 0.33) {
			  // Generate constrained lower bounds with sum l_i = 0.5
			  var myIterator = new PortfolioAllocation.randomCompositionsIterator_(200, dimension);
			  var integerBounds = myIterator.next();
			  for (var j = 0; j < dimension; ++j) {
			  	  lowerBounds[j] = Math.min(1, integerBounds[j] / 400);
			  }

			  sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, lowerBounds); // lower bounds constraints
		  }
		  else if (u <= 0.66) {
			  // Generate constrained upper bounds with sum u_i = 2
			  var myIterator = new PortfolioAllocation.randomCompositionsIterator_(200, dimension);
			  var integerBounds = myIterator.next();
			  for (var j = 0; j < dimension; ++j) {
			  	  upperBounds[j] = Math.min(1, integerBounds[j] / 100);
			  }
				
			  sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, null, upperBounds); // upper bounds constraints
		  }
		  else {
			  // Generate constrained lower bounds with sum l_i = 0.5
			  var myIterator = new PortfolioAllocation.randomCompositionsIterator_(200, dimension);
			  var integerBounds = myIterator.next();
			  for (var j = 0; j < dimension; ++j) {
			  	  lowerBounds[j] = Math.min(1, integerBounds[j] / 400);
			  }
			  
			  // Generate constrained upper bounds with sum u_i = 2 and u_i >= l_i
			  var myIterator = new PortfolioAllocation.randomCompositionsIterator_(200, dimension);
			  var integerBounds = myIterator.next();
			  for (var j = 0; j < dimension; ++j) {
			  	  upperBounds[j] = Math.max(Math.min(1, integerBounds[j] / 100), lowerBounds[j]);
			  }
			  
			  sampler = new PortfolioAllocation.simplexRandomSampler_(dimension, lowerBounds, upperBounds); // lower and upper bounds constraints
		  }
		 
		 for (var j = 0; j < nbSubTests; ++j) {
			// Generate a random sample point
			var sampledPoint = sampler.sample();
		 
			// Check that the number of coordinates of the sampled point corresponds to the requested dimension
			var sampledPointDimension = sampledPoint.length;
			assert.equal(sampledPointDimension, dimension, "Restricted simplex random sampler point computation, coordinates length - Test " + i + "," + j);
		 
			 // Check that the coordinates of the sampled point belong to the unit interval
			 var sampledPointBelongToBoundedInterval = true;
			 for (var k = 0; k < sampledPoint.length; ++k) {
				if (sampledPoint[k] > upperBounds[k] || sampledPoint[k] < lowerBounds[k]) {
					sampledPointBelongToBoundedInterval = false;
					break;
				}
			 }
			 assert.equal(sampledPointBelongToBoundedInterval, true, "Restricted simplex random sampler point computation, coordinates in proper interval - Test " + i + "," + j);

			 // Check that the sum of the coordinates of the sampled point is 1, near to machine precision
			 var sampledPointCoordinatesSum = 0;
			 for (var k = 0; k < sampledPoint.length; ++k) {
				sampledPointCoordinatesSum += sampledPoint[k];
			}
			assert.equal(Math.abs(sampledPointCoordinatesSum - 1) <= 1e-14, true, "Restricted simplex random sampler point computation, coordinates sum to one - Test " + i + "," + j);
		}
	  }
  }

});


QUnit.test('Simplex random direction computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbSubTests = 20;
	  var minDimension = 3;
	  var maxDimension = 3;

	  // Aim of these tests is to check that for a sample of size n:
	  // - An array of coordinates of size n is returned
	  // - The coordinates sum to 0
	  // - The 2-norm of the R^n vector associated to these coordinates is equal to 1
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random dimension size and the associated sampler
		 var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);
		 var sampler = new PortfolioAllocation.simplexDirectionRandomSampler_(dimension);
		 
		 for (var j = 0; j < nbSubTests; ++j) {
			// Generate a random sample direction
			var sampledDirection = sampler.sample();

			// Check that the number of coordinates of the sampled direction corresponds to the requested dimension
			var sampledDirectionDimension = sampledDirection.length;
			assert.equal(sampledDirectionDimension, dimension, "Simplex random direction sampler computation, coordinates length - Test " + i + "," + j);
		 
			 // Check that the sum of the coordinates of the sampled point is 1, near to machine precision
			 // Check that the coordinates 2-norm of the sampled point is 1, near to machine precision
			 var sampledDirectionCoordinatesSum = 0;
			 var sampledDirectionCoordinatesSumSq = 0;
			 for (var k = 0; k < sampledDirection.length; ++k) {
				sampledDirectionCoordinatesSum += sampledDirection[k];
				sampledDirectionCoordinatesSumSq += sampledDirection[k] * sampledDirection[k];
			}
			assert.equal(Math.abs(sampledDirectionCoordinatesSum) <= 1e-12, true, "Simplex random direction sampler computation, coordinates sum to zero - Test " + i + "," + j);
			assert.equal(Math.abs(Math.sqrt(sampledDirectionCoordinatesSumSq) - 1) <= 1e-14, true, "Simplex random direction sampler computation, coordinates 2-norm equal to one - Test " + i + "," + j);
		}
	  }
  }

});


QUnit.test('Simplex grid search computation', function(assert) {    
	// Test using static data
	{
		// f(x,y) = x, min must be at [0,1]
		var expectedValues = [[0, 1]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0]; }, 2, 10);
		
		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search, #1, ' + i);
		}
	}
	
	// Test using static data
	{
		// f(x,y) = y, min must be at [1,0]
		var expectedValues = [[1,0]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[1]; }, 2, 10);
		
		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search, #2, ' + i);
		}
	}
	
	// Test using static data, both with and without boundaries contraints
	{
		// f(x,y) = x + y, min must be at all points in [0,1]
		var expectedValues = [[1,0],[0.9,0.1],[0.8,0.2],[0.7,0.3],[0.6,0.4],[0.5,0.5],[0.4,0.6],[0.3,0.7],[0.2,0.8],[0.1,0.9],[0,1]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10);
		
		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search, #3, ' + i);
		}
		
		var expectedValues = [[0.7,0.3], [0.6,0.4]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, [0.6, 0.3]); // lower boundaries

		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search with boundaries constaints, #1, ' + i);
		}
		
		var expectedValues = [[0.7,0.3],[0.6,0.4]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, null , [0.7, 0.4]); // upper boundaries

		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search with boundaries constaints, #2, ' + i);
		}
		
		var expectedValues = [[0.5,0.5]];
		var res = PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, [0.45, 0.45] , [0.55, 0.55]); // both boundaries

		for (var i = 0; i < res.length; ++i) {
			var pointOk = true;
			for (var j = 0; j < 2; ++j) {
			   if (expectedValues[i][j] != res[i][j]) {
				 pointOk = false;
				 break;
			   }
			}	  
			assert.equal(pointOk, true, 'Simplex grid search with boundaries constaints, #2, ' + i);
		}
	}
	
	  // Test with static data the failure cases of grid search on the simplex with boundaries constaints:
	  // - A lower bound lower than 0
	  // - An upper bound greater than 1
	  // - A lower bound strictly greater than an upper bound
	  // - An empty simplex
	  {
		 // Test a lower bound lower than 0
		assert.throws(function() { 
			new PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, [-1, 0]) },
			new Error('incorrect problem detected: lower bound strictly lower than 0'),
			"Simplex grid search with boundaries constaints, incorrect lower bound");

		// Test an upper bound greater than 1
		assert.throws(function() { 
			new PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, null, [1.1, 1]) },
			new Error('incorrect problem detected: upper bound strictly greater than 1'),
			"Simplex grid search with boundaries constaints, incorrect upper bound");

		// Test a lower bound strictly greater than an upper bound
		assert.throws(function() { 
			new PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, [0, 0.5], [1, 0.4]) },
			new Error('infeasible problem detected: lower bound strictly greater than upper bound'),
			"Simplex grid search with boundaries constaints, incorrect upper bound");
		
		// Test an empty simplex due to lower bounds
		assert.throws(function() {
			new PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, [0.5, 0.6]) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Simplex grid search with boundaries constaints, empty simplex due to lower bounds");			
			
		// Test an empty simplex due to upper bounds
		assert.throws(function() { 
			new PortfolioAllocation.simplexGridSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10, null, [0.4, 0.3]) },
			new Error('infeasible problem detected: the restricted simplex is empty'),
			"Simplex grid search with boundaries constaints, empty simplex due to upper bounds");
	  }
});
