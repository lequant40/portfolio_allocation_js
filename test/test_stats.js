// ------------------------------------------------------------
QUnit.module('Statistics internal module', {
  before: function() {
    // 
  }
});



QUnit.test('Median computation', function(assert) {    
  // Test with static data
  {
	  assert.equal(PortfolioAllocation.median_([2,4,1]), 2, 'Median computation #1');
	  assert.equal(PortfolioAllocation.median_([2,4,1,3]), 2, 'Median computation #2');
  }  
  
  //TODO: use random data
});


QUnit.test('Smallest k element computation', function(assert) {    
	function generateRandomDimension(min, max) { // used for n
		return Math.floor(Math.random()*(max-min+1) + min);
	}
	
	function generateRandomValue(minVal, maxVal) { // used for each array element		
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	function generateRandomArray(n) { // used for arr
		var minVal = -10;
		var maxVal = 10;
		
		var arr = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		
		for (var i = 0; i < n; ++i) {
			arr[i] = generateRandomValue(minVal, maxVal);
		}
		
		return arr;
	}
	
  // Test with random data
  {
	var nbTests = 100;
	for (var j = 0; j < nbTests; ++j) {
	  // Generate a random array of a random size between 1 and 1000
	  var n = generateRandomDimension(1, 1000);
	  var arr = generateRandomArray(n);
		
	  // Sort (ascending) a copy of the array, to be used for the expected values
	  var copy_arr = arr.slice().sort(function(a, b) { return a - b; });
	  
	  // Generate a random integer k between 1 and the size of the array
	  var k = generateRandomDimension(1, n);
	  
	  // Compute the k-th smallest element of the array arr using the function
	  // SELECT
	  var k_smallest_elem = PortfolioAllocation.select_(arr, k);
	  
	  // Test that the k-th smallest element of the array arr, computed
	  // from the function SELECT, is the same as arr[k-1].
	  assert.equal(k_smallest_elem, arr[k-1], 'Smallest k element computation #1 ' + (j+1));

	  // Test that the k-th smallest element of the array arr, computed
	  // from the function SELECT, is the same as the k-th element of 
	  // the sorted array copy_arr.
	  assert.equal(k_smallest_elem, copy_arr[k-1], 'Smallest k element computation #2 ' + (j+1));
	  
	  // Test that smallest k elements of arr are x[i], i=0..k-1
	  var correct_order_left_k = true;
	  for (var i = 0; i < k; ++i) {
		if (arr[i] > k_smallest_elem) {
			correct_order_left_k = false;
			break;
		}
	  }
	  assert.equal(correct_order_left_k, true, 'Smallest k element computation #3 ' + (j+1));
	  
	  // Test that largest n-k elements of arr are x[i], i=k+1..n-1
	  var correct_order_right_k = true;
	  for (var i = k; i < n; ++i) {
		if (arr[i] < k_smallest_elem) {
			correct_order_right_k = false;
			break;
		}
	  }
	  assert.equal(correct_order_right_k, true, 'Smallest k element computation #4 ' + (j+1));
    }
  }  
});

QUnit.test('Hypothenuse computation', function(assert) {    
  // Tests with static data
  {
	  assert.ok(Math.abs(PortfolioAllocation.hypot_(3e200, 4e200) - 5e+200)/5e+200 <= 1e-15, 'Hypothenuse computation with no overflow');
	  assert.equal(PortfolioAllocation.hypot_(3, 4), 5, 'Hypothenuse computation #1');
	  assert.equal(PortfolioAllocation.hypot_(-2, 0), 2, 'Hypothenuse computation one zero argument #1');
	  assert.equal(PortfolioAllocation.hypot_(0, -2), 2, 'Hypothenuse computation one zero argument #2');
	  assert.equal(PortfolioAllocation.hypot_(0, 0), 0, 'Hypothenuse computation two zero arguments');
  }  
  
  // Tests with the formula and random data
  {
	  var nbTests = 50;
	  for (var i = 0; i < nbTests; ++i) {
          var x = Math.random();
		  var y = Math.random();
		  var naiveHypothenuse = Math.sqrt(x*x + y*y);
		  assert.ok(Math.abs(PortfolioAllocation.hypot_(x, y) - naiveHypothenuse) <= 1e-14, 'Hypothenuse computation #' + (i + 2));
	  }
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
	  var corrMat = [[1.0000000,  0.2378848,  0.2483431,  0.3163914,  0.1796639], [0.2378848,  1.0000000,  0.2487009, -0.1986677, -0.2165444], [0.2483431,  0.2487009,  1.0000000, -0.3179188,  0.3713964], [0.3163914, -0.1986677, -0.3179188, 1.0000000,  0.4131639],[0.1796639, -0.2165444,  0.3713964,  0.4131639,  1.0000000]];
	
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
	var corrMat = [[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]];
	var clusters = PortfolioAllocation.ftca_(corrMat, 1);
	var expectedClusters = [[2],[4],[3],[1]];
	assert.deepEqual(clusters, expectedClusters, "FTCA - Test 2");
  }
  
});
