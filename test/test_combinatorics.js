// ------------------------------------------------------------
QUnit.module('Combinatorics internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Alias method sampler', function(assert) {    
  // Test with random data
  {
	  // Setup parameters of the random test
	  var nbTests = 10;
	  var minSizeProbabilityDistribution = 1;
	  var maxSizeProbabilityDistribution = 50;
	  var nbSamplesToDraw = 1000;
	  
	  // For each test, generate a probability distribution, and then sample
	  // from this probability distribution, ensuring the generated values are coherent.
	  for (var i = 0; i < nbTests; ++i) {
		// Generate a random discrete probability distribution of a random size
		var sizeProbabilityDistribution = Math.floor(Math.random()*(maxSizeProbabilityDistribution - minSizeProbabilityDistribution + 1) + minSizeProbabilityDistribution);
		var probabilityDistribution = new Array(sizeProbabilityDistribution);
		var normalizationConstant = 0;
		for (var j = 0; j < probabilityDistribution.length; ++j) { 
			probabilityDistribution[j] = Math.random(); 
			normalizationConstant += probabilityDistribution[j];
		}
		for (var j = 0; j < probabilityDistribution.length; ++j) { 
			probabilityDistribution[j] /= normalizationConstant;
		}
		
		// Construct an alias sampler associated to the generated probability distribution.
		var aliasSampler = new PortfolioAllocation.aliasMethodSampler_(probabilityDistribution);
		
		// Sample from the alias sampler.
		var sampleOk = true;
		for (var j = 0; j < nbSamplesToDraw; ++j) {
			var s = aliasSampler.sample();
			if (!(s >= 0 && s <= probabilityDistribution.length - 1)) {
				sampleOk = false;
				break;
			}
		}
		assert.equal(sampleOk, true, "Alias method sampler, elements in correct range - Test " + i);
	  }	  
  }
});


QUnit.test('Next composition computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 5
  {
	var expectedValues = [[6,0,0],[5,1,0],[4,2,0],[3,3,0],[2,4,0],[1,5,0], [0,6,0], [5,0,1], [4,1,1], [3,2,1], [2,3,1], [1,4,1], [0,5,1], [4,0,2], [3,1,2], [2,2,2], [1,3,2], [0,4,2], [3,0,3], [2,1,3], [1,2,3], [0,3,3], [2,0,4], [1,1,4], [0,2,4], [1,0,5], [0,1,5], [0,0,6], -1];

	var nextCompositionIterator = new PortfolioAllocation.compositionsIterator_(6, 3);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var nextComposition = nextCompositionIterator.next();
		
		var compositionOk = true;
		for (var j = 0; j < 3; ++j) {
		   if (expectedValues[i][j] != nextComposition[j]) {
		     compositionOk = false;
		     break;
		   }
		}	  
		assert.equal(compositionOk, true, 'Next composition - Test 1 #' + i);
	  }  
  }

});

QUnit.test('Next random composition computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbCompositions = 10;
	  var minN = 1;
	  var maxN = 1000;

	  // Aim of these tests is to check that for any generated k-composition of an integer, 
	  // an array of k integers summing to n is returned
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random integer, and an associated random k-composition iterator
		 var n = Math.floor(Math.random()*(maxN - minN + 1) + minN);
		 var k  = Math.floor(Math.random()*(n - 1 + 1) + 1);
		 var nextRandomCompositionIterator = new PortfolioAllocation.randomCompositionsIterator_(n, k);
		 
		 for (var j = 0; j < nbCompositions; ++j) {
			 // Generate a random k-composition
			 var generatedComposition = nextRandomCompositionIterator.next();
			 
			 // Compute the length of the generated composition
			 var generatedCompositionLength = generatedComposition.length;
			 assert.equal(generatedCompositionLength, k, "Next random composition computation, k-composition length - Test " + i + "," + j);
			 
			 // Check that the generated composition elements sum to n
			 var generatedCompositionSum = 0;
			 for (var k = 0; k < generatedComposition.length; ++k) {
				generatedCompositionSum += generatedComposition[k];
			 }
			 assert.equal(generatedCompositionSum, n, "Next random composition computation, k-composition elements sum to n - Test " + i + "," + j);
		  }
	  }
  }

});


QUnit.test('Next subset computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 1
  {
	  var expectedValues = [[],[1],[1,2],[2],[2,3],[1,2,3],[1,3],[3],[3,4],[1,3,4],[1,2,3,4],[2,3,4],[2,4],[1,2,4],[1,4],[4],[4,5],[1,4,5],[1,2,4,5],[2,4,5],[2,3,4,5],[1,2,3,4,5],[1,3,4,5],[3,4,5],[3,5],[1,3,5],[1,2,3,5],[2,3,5],[2,5],[1,2,5],[1,5], [5]];
	  
	  var nextSubsetIterator = new PortfolioAllocation.subsetsIterator_(5);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var nextSubset = nextSubsetIterator.next();
		
		var subsetOk = true;
		for (var j = 0; j < nextSubset.length; ++j) {
		   if (expectedValues[i][j] != nextSubset[j]) {
		     subsetOk = false;
		     break;
		   }
		}	  
		assert.equal(subsetOk, true, 'Next subset - Test 1 #' + i);
	  }
	  
	  var nextSubset = nextSubsetIterator.next();
	  assert.equal(nextSubset, -1, 'Next subset - End value');
  }

});



QUnit.test('Next k-subset computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 3
  {
	  var expectedValues = [[1,2,3], [1,2,4], [1,2,5], [1,3,4], [1,3,5], 
	                        [1,4,5], [2,3,4],[2,3,5], [2,4,5], [3,4,5]];
	  
	  var nextKSubsetIterator = new PortfolioAllocation.kSubsetsIterator_(5, 3);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var nextKSubset = nextKSubsetIterator.next();
		
		var subsetOk = true;
		for (var j = 0; j < nextKSubset.length; ++j) {
		   if (expectedValues[i][j] != nextKSubset[j]) {
		     subsetOk = false;
		     break;
		   }
		}	  
		assert.equal(subsetOk, true, 'Next k-subset - Test 1 #' + i);
	  }

	  var nextKSubset = nextKSubsetIterator.next();
	  assert.equal(nextKSubset, -1, 'Next k-subset - End value');
  }

});


QUnit.test('Next random k-subset computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbSubTests = 100;
	  var minN = 1;
	  var maxN = 1000;

	  // Aim of these tests is to check that for any generated k-subset of a n-set, 
	  // an array of k strictly increasing integers belonging to the n-set is returned
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random n-set, and an associated k-subset iterator
		 var n = Math.floor(Math.random()*(maxN - minN + 1) + minN);
		 var k  = Math.floor(Math.random()*(n - 1 + 1) + 1);
		 var nextRandomKSubsetIterator = new PortfolioAllocation.randomKSubsetIterator_(n, k);
		 
		 for (var j = 0; j < nbSubTests; ++j) {
			 // Generate a random k-subset
			 var generatedSubset = nextRandomKSubsetIterator.next();
			 
			 // Compute the length of the generatedSubset
			 var generatedSubsetLength = generatedSubset.length;
			 assert.equal(generatedSubsetLength, k, "Next random k-subset computation, k-subset length - Test " + i + "," + j);
			 
			 // Check that the subset elements belong to the n-set
			 var generatedSubsetElementsBelongToNSet = true;
			 for (var k = 0; k < generatedSubset.length; ++k) {
				if (generatedSubset[k] > n || generatedSubset[k] < 0) {
					generatedSubsetElementsBelongToNSet = false;
					break;
				}
			 }
			 assert.equal(generatedSubsetElementsBelongToNSet, true, "Next random k-subset computation, elements in n-set - Test " + i + "," + j);

			 // Check that the subset elements are in strictly increasing order
			 var generatedSubsetElementsInOrder = true;
			 for (var k = 1; k < generatedSubset.length; ++k) {
				if (generatedSubset[k] <= generatedSubset[k-1]) {
					generatedSubsetElementsInOrder = false;
					break;
				}
			 }
			assert.equal(generatedSubsetElementsInOrder, true, "Next random k-subset computation, elements in strictly increasing order - Test " + i + "," + j);
		  }
	  }
  }

});


QUnit.test('Binomial coefficient computation', function(assert) {    
	// Test with the static data examples
	{
		assert.equal(PortfolioAllocation.binomial_(0, 0), 1, 'Binomial coefficient - Test 1');
		assert.equal(PortfolioAllocation.binomial_(7, 5), 21, 'Binomial coefficient - Test 2');
		assert.equal(PortfolioAllocation.binomial_(20, 15), 15504, 'Binomial coefficient - Test 3');
		assert.equal(PortfolioAllocation.binomial_(63, 7), 553270671, 'Binomial coefficient - Test 4');
		assert.equal(PortfolioAllocation.binomial_(25, 6), 177100, 'Binomial coefficient - Test 5');
		assert.equal(PortfolioAllocation.binomial_(5+100-1, 5-1), 4598126, 'Binomial coefficient - Test 6');
	}
	
	// Test using random data and the recursive formula defining the binomial coefficients
	{
		// Generate a random integer n
		var n = Math.floor(Math.random() * (50 - 2 + 1)) + 2; // max 50 min 2
		
		// Generate a random integer 1 <= k <= n-1
		var k = Math.floor(Math.random() * (n-1 - 1 + 1)) + 1; // max n-1 min 1
		
		// Compute initial values
		assert.equal(PortfolioAllocation.binomial_(n, 0), 1, 'Binomial coefficient - Test 7 #1');
		assert.equal(PortfolioAllocation.binomial_(n, n), 1, 'Binomial coefficient - Test 7 #2');
		
		// Compute binomial(n,k) and compare it to binomial(n-1,k-1) + binomial(n-1,k)
		assert.equal(PortfolioAllocation.binomial_(n, k), 
		             PortfolioAllocation.binomial_(n-1, k-1) + PortfolioAllocation.binomial_(n-1, k), 
					 'Binomial coefficient - Test 7 #3');
	}
  
});