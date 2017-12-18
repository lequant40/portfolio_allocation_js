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



QUnit.test('Simplex rational sampler point computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 5, divided by 6
  {
	  var expectedValues = [[6/6,0,0], [5/6,1/6,0], [4/6,2/6,0], [3/6,3/6,0], [2/6,4/6,0], [1/6,5/6,0], [0,6/6,0], [5/6,0,1/6], [4/6,1/6,1/6], [3/6,2/6,1/6], [2/6,3/6,1/6], [1/6,4/6,1/6], [0,5/6,1/6], [4/6,0,2/6], [3/6,1/6,2/6], [2/6,2/6,2/6], [1/6,3/6,2/6], [0,4/6,2/6], [3/6,0,3/6], [2/6,1/6,3/6], [1/6,2/6,3/6], [0,3/6,3/6], [2/6,0,4/6], [1/6,1/6,4/6], [0,2/6,4/6], [1/6,0,5/6], [0,1/6,5/6], [0,0,6/6]];
	  
	  var sampler = new PortfolioAllocation.simplexRationalSampler(3, 6);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var sample = sampler.nextDeterministicSample();
		assert.deepEqual(sample, expectedValues[i], 'Simplex rational sampler - Test 1 #' + i);
	  }
	  
	  assert.equal(sampler.nextDeterministicSample(), null, "Simplex rational sampler - Test 1 #end");
  }
});
