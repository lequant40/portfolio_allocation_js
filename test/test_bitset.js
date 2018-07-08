// ------------------------------------------------------------
QUnit.module('Bit set internal module', {
  before: function() {
    // 
  }
});



QUnit.test('Bit set data structure', function(assert) {    
  // Test with static data
  {
	  // Initialize the bit set
	  var bs = new PortfolioAllocation.BitSet_();
	  
	  // Test to string
	  assert.equal(bs.toString(), "", 'Bit set, static tests - initial base-2 string representation');

	  // Test cardinality computation
	  assert.equal(bs.nbSetBits(), 0, 'Bit set, static tests - initial cardinality');

	  // Test isEmpty computation
	  assert.equal(bs.isEmpty(), true, 'Bit set, static tests - initial isEmpty');
	  
	  // Test to array
	  assert.deepEqual(bs.toArray(), [], 'Bit set, static tests - initial array representation');
	  
	  // Test setting one bit to 1
	  assert.equal(bs.get(31), false, 'Bit set, static tests - set one bit to 1 1/2');
	  bs.set(31)
	  assert.equal(bs.get(31), true, 'Bit set, static tests - set one bit to 1 2/2');
	  
	  // Test to string
	  assert.equal(bs.toString(), "10000000000000000000000000000000", 'Bit set, static tests - base-2 string representation with one bit set to 1');

	  // Test cardinality computation
	  assert.equal(bs.nbSetBits(), 1, 'Bit set, static tests - cardinality with one bit set to 1');

	  // Test isEmpty computation
	  assert.equal(bs.isEmpty(), false, 'Bit set, static tests - isEmpty with one bit set to 1');
	  
	  // Test to array
	  assert.deepEqual(bs.toArray(), [31], 'Bit set, static tests - array representation with one bit set to 1');
	  
	  // Test unsetting one bit
	  bs.unset(31)
	  assert.equal(bs.get(31), false, 'Bit set, static tests - unset one bit to 1');
	  
	  // Test flipping one bit
	  bs.flip(31)
	  assert.equal(bs.get(31), true, 'Bit set, static tests - flip one bit to 1');
	  
	  // Test resize through get, string, cardinality, isEmpty and to array
	  // Note: the bit set is on purpose creating an empty word in the "middle" of the bit set underlying array
	  bs.set(94);
	  assert.equal(bs.get(94), true, 'Bit set, static tests - set one bit to 1 to resize 1/5');
	  assert.equal(bs.toString(), "100000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000", 'Bit set, static tests - set one bit to 1 to resize 2/5');
	  assert.equal(bs.nbSetBits(), 2, 'Bit set, static tests - set one bit to 1 to resize 3/5');
	  assert.equal(bs.isEmpty(), false, 'Bit set, static tests - set one bit to 1 to resize 4/5');
	  assert.deepEqual(bs.toArray(), [31, 94], 'Bit set, static tests - set one bit to 1 to resize 5/5');
	  
	  // Test setting a range of bits to 1 through iteration
	  bs.setRange(1,3);
	  var bs_it = new bs.iterator();
	  assert.deepEqual(bs_it.next(), 1, 'Bit set, static tests - iteration over two bits set to 1 1/6');
	  assert.deepEqual(bs_it.next(), 2, 'Bit set, static tests - iteration over two bits set to 1 2/6');
	  assert.deepEqual(bs_it.next(), 3, 'Bit set, static tests - iteration over two bits set to 1 3/6');
	  assert.deepEqual(bs_it.next(), 31, 'Bit set, static tests - iteration over two bits set to 13 4/6');
	  assert.deepEqual(bs_it.next(), 94, 'Bit set, static tests - iteration over two bits set to 1 5/6');
	  assert.deepEqual(bs_it.next(), 0, 'Bit set, static tests - iteration over two bits set to 1 6/6');
	  
	  // Test clear through isEmpty
	  bs.clear();
	  assert.equal(bs.isEmpty(), true, 'Bit set, static tests - bit set cleared');
	  
  }  
  
  //TODO: use random data
});
