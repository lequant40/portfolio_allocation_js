// ------------------------------------------------------------
QUnit.module('Optimisation internal module', {
  before: function() {
    // 
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