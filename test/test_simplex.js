// ------------------------------------------------------------
QUnit.module('Simplex internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Closest rational point computation', function(assert) {    
  // Test with static data
  {
	  var testValues = [[0.7373, 0.2627, 0], [0.5759, 0.0671, 0.3570], [0.22, 0.66, 0.11, 0.01], [0.22, 0.66, 0.11, 0.01], [0.5, 0.5]];
	  var testGridIndices = [10, 10, 1, 5, 1];
	  var expectedValues = [[0.70, 0.30, 0], [0.60, 0.10, 0.30], [0, 1, 0, 0], [0.2, 0.6, 0.2, 0], [1, 0]];
	  
	  for (var i = 0; i < testValues.length; ++i) {
	      assert.deepEqual(PortfolioAllocation.closestRationalPoint_(testValues[i], testGridIndices[i]), expectedValues[i], 'Closest rational point #' + i);
		  //assert.ok(Math.abs( (PortfolioAnalytics.normcdf_(x[i]) - p[i]) ) <= 8e-14, 'Normcdf comparison v.s. Wolfram Alpha');
      } 
  }  

});
