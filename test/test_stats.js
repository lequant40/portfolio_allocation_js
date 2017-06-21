// ------------------------------------------------------------
QUnit.module('Statistics internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Rank computation', function(assert) {    
  // Test with static data
  {
	  var testValues = [12, 13, 15, 10, 12];
	  var expectedRanksDescending = [3, 2, 1, 5, 3];
	  var expectedRanksAscending = [2, 4, 5, 1, 2];
	  
	  assert.deepEqual(PortfolioAllocation.rank_(testValues, 0), expectedRanksDescending, 'Rank descending');
	  assert.deepEqual(PortfolioAllocation.rank_(testValues, 1), expectedRanksAscending, 'Rank ascending');
  }  
  
  // Second test with static data
  {
	var testValues = [0.13, 0.63, 0.79]; 
	var expectedRanks = [3, 2, 1];
	assert.deepEqual(PortfolioAllocation.rank_(testValues, 0), expectedRanks, 'Rank #2');
  }
});

