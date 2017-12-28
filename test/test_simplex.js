// ------------------------------------------------------------
QUnit.module('Simplex internal module', {
  before: function() {
    // 
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


QUnit.test('Simplex deterministic rational sampler point computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 5, divided by 6
  {
	  var expectedValues = [[6/6,0,0], [5/6,1/6,0], [4/6,2/6,0], [3/6,3/6,0], [2/6,4/6,0], [1/6,5/6,0], [0,6/6,0], [5/6,0,1/6], [4/6,1/6,1/6], [3/6,2/6,1/6], [2/6,3/6,1/6], [1/6,4/6,1/6], [0,5/6,1/6], [4/6,0,2/6], [3/6,1/6,2/6], [2/6,2/6,2/6], [1/6,3/6,2/6], [0,4/6,2/6], [3/6,0,3/6], [2/6,1/6,3/6], [1/6,2/6,3/6], [0,3/6,3/6], [2/6,0,4/6], [1/6,1/6,4/6], [0,2/6,4/6], [1/6,0,5/6], [0,1/6,5/6], [0,0,6/6]];
	  
	  var sampler = new PortfolioAllocation.simplexDeterministicRationalSampler_(3, 6);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var sample = sampler.sample();
		assert.deepEqual(sample, expectedValues[i], 'Simplex deterministic rational sampler - Test 1 #' + i);
	  }
	  
	  assert.equal(sampler.sample(), null, "Simplex deterministic rational sampler - Test 1 #end");
  }
});


QUnit.test('Simplex random sampler point computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbSubTests = 10;
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

});

QUnit.test('Simplex rational grid search computation', function(assert) {    
	// Test using static data
	{
		// f(x,y) = x, min must be at [0,1]
		assert.deepEqual(PortfolioAllocation.simplexRationalGirdSearch_(function(arr) { return arr[0]; }, 2, 10), [[0,1]], 'Simplex rational grid search Test #1');
	}
	
	// Test using static data
	{
		// f(x,y) = y, min must be at [1,0]
		assert.deepEqual(PortfolioAllocation.simplexRationalGirdSearch_(function(arr) { return arr[1]; }, 2, 10), [[1,0]], 'Simplex rational grid search Test #2');
	}
	
	// Test using static data
	{
		// f(x,y) = x + y, min must be at all points in [0,1]
		assert.deepEqual(PortfolioAllocation.simplexRationalGirdSearch_(function(arr) { return arr[0] + arr[1]; }, 2, 10), 
		                 [[1,0],[0.9,0.1],[0.8,0.2],[0.7,0.3],[0.6,0.4],[0.5,0.5],[0.4,0.6],[0.3,0.7],[0.2,0.8],[0.1,0.9],[0,1]], 
						 'Simplex rational grid search Test #3');
	}

	

});