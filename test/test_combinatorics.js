// ------------------------------------------------------------
QUnit.module('Combinatorics internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Next composition computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 5
  {
	  var expectedValues = [[true,[6,0,0]],[true, [5,1,0]],[true, [4,2,0]],[true, [3,3,0]],[true, [2,4,0]],[true, [1,5,0]],[true, [0,6,0]],[true, [5,0,1]],[true, [4,1,1]],[true, [3,2,1]],[true, [2,3,1]],[true, [1,4,1]],[true, [0,5,1]],[true, [4,0,2]],[true, [3,1,2]],[true, [2,2,2]],[true, [1,3,2]],[true, [0,4,2]],[true, [3,0,3]],[true, [2,1,3]],[true, [1,2,3]],[true, [0,3,3]],[true, [2,0,4]],[true, [1,1,4]],[true, [0,2,4]],[true, [1,0,5]],[true, [0,1,5]],[false, [0,0,6]]];
	  
	  var nextCompositionIterator = new PortfolioAllocation.compositionsIterator_(6, 3);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var nextComposition = nextCompositionIterator.next();
		assert.deepEqual(nextComposition, expectedValues[i], 'Next composition - Test 1 #' + i);
	  }  
  }

});


QUnit.test('Next subset computation', function(assert) {    
  // Reference: Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
  // Test with the static data examples of section 1
  {
	  var expectedValues = [[true,[]], [true,[1]], [true,[1,2]], [true,[2]], [true,[2,3]], [true,[1,2,3]], [true,[1,3]],[true,[3]], [true,[3,4]], [true,[1,3,4]], [true,[1,2,3,4]], [true,[2,3,4]], [true,[2,4]], [true,[1,2,4]], [true,[1,4]], [true,[4]], [true,[4,5]], [true,[1,4,5]], [true,[1,2,4,5]], [true,[2,4,5]], [true,[2,3,4,5]], [true,[1,2,3,4,5]], [true,[1,3,4,5]], [true,[3,4,5]], [true,[3,5]], [true,[1,3,5]], [true,[1,2,3,5]], [true,[2,3,5]], [true,[2,5]], [true,[1,2,5]], [true,[1,5]], [false,[5]]];
	  
	  var nextSubsetIterator = new PortfolioAllocation.subsetsIterator_(5);
	  for (var i = 0; i < expectedValues.length; ++i) {
		var nextSubset = nextSubsetIterator.next();
		assert.deepEqual(nextSubset, expectedValues[i], 'Next subset - Test 1 #' + i);
	  }  
  }

});


QUnit.test('Next random k-subset computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 10;
	  var nbSubTests = 10;
	  var minN = 1;
	  var maxN = 50;

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