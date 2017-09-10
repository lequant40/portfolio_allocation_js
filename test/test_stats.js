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



QUnit.test('FTCA computation', function(assert) {    
  // FTCA test, using static data
  {
	  var corrMat = new PortfolioAllocation.Matrix([[1.0000000,  0.2378848,  0.2483431,  0.3163914,  0.1796639], [0.2378848,  1.0000000,  0.2487009, -0.1986677, -0.2165444], [0.2483431,  0.2487009,  1.0000000, -0.3179188,  0.3713964], [0.3163914, -0.1986677, -0.3179188, 1.0000000,  0.4131639],[0.1796639, -0.2165444,  0.3713964,  0.4131639,  1.0000000]]);
	
      // Compute FTCA with varying thresholds
	  // A threshold of 1 will always output as many clusters as initial elements
	  var thresholds = [-0.10, 0.0, 0.10, 0.25, 0.33, 0.5, 1];
	  var expectedClusters = [ [[1,2,3,4,5]], [[1,2,3,4],[5]], [[1,2,3],[5,4]], [[1,4],[2],[5,3]], [[1],[2],[5,3],[4]], [[1],[2],[5],[3],[4]], [[1],[2],[5],[3],[4]]];
	  for (var i = 0; i < thresholds.length; ++i) {
		var clusters = PortfolioAllocation.ftca_(corrMat, thresholds[i]);
		assert.deepEqual(clusters, expectedClusters[i], "FTCA - Test 1 #" + i);
	  }	  
  }
  
  // FTCA test, limit case with static data
  {
	var corrMat = new PortfolioAllocation.Matrix([[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]]);
	var clusters = PortfolioAllocation.ftca_(corrMat, 1);
	var expectedClusters = [[2],[4],[3],[1]];
	assert.deepEqual(clusters, expectedClusters, "FTCA - Test 2");
  }
  
});


  