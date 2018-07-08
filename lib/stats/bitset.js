/**
 * @file Functions related to bit set object.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
self.BitSet_ = BitSet_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function BitSet_
*
* @summary Construct a bit set.
*
* @description This function constructs an empty bit set (a.k.a. bit array, bit vector),
* which is a data structure taylored to storing (relatively small) integers.
*
* The internal way to handle bit sets has been fully adapted from: 
* - https://github.com/infusion/BitSet.js
* - https://github.com/lemire/FastBitSet.js
*
* @see <a href="https://en.wikipedia.org/wiki/Bit_array">Bit array</a>
*
* @return {this} the constructed bit set.
*
* @example
* var myBitSet = new BitSet_();
*/
function BitSet_() {
    // Catches incorrect usage of var b = BitSet_() instead of var b = new BitSet_()
	if (!(this instanceof BitSet_)) {
      return new BitSet_();
    }

	// The number of bits hold in a word
	this.WORD_LENGTH = 32;

	// The log base 2 of WORD_LENGTH
	this.WORD_LOG = 5;
	
	// The "list" of words hold by the BitSet
	this.words = typeof Uint32Array === 'function' ? new Uint32Array(0) : new Array(0);
	
	// For subsequent usage in initialization sub-functions
	var that = this;
	
	/**
	* @function iterator
	*
	* @summary Returns an iterator to compute the indexes
	* corresponding to the bits set to 1 in the bit set.
	*
	* @description This function constructs an iterator to compute the indexes
	* corresponding to the bits set to 1 in the bit set.
	*
	* To be noted that once an index corresponding to a bit set to 1 has been
	* iterated over by the iterator, this index and all the indexes lower than
	* this index can be altered (i.e., set to 0) without invalidating the iterator.
	*
	* @memberof BitSet_
	* @return {function} a function to be used as an iterator through its .next() method, computing  
	* the indexes corresponding to the bits set to 1 in the bit set.
	*
	* @example
	* var myBitSet = new BitSet_().add(2);
	* var myIterator = new myBitSet.iterator();
	* myIterator.next(); myIterator.next();
	* // 2; -1;
	*/
	this.iterator = function() {
		// Initialize the current words index and the current word to
		// the first non-null word, if existing.
		this.w_idx = -1;
		this.w = 0;
		while (this.w_idx < that.words.length && this.w == 0) {
			++this.w_idx;
			this.w = that.words[this.w_idx];
		}

		/**
		* @function next
		*
		* @summary Returns the index corresponding to the next bit set to 1 in the bit set.
		*
		* @description This function computes the index corresponding to the next bit set to 1
		* in the bit set.
		*
		* The initial index computed by the first call to this function is the index corresponding
		* to the first bit set to 1, and each subsequent call to this function will result in computing
		* the index corresponding to the next bit set to 1, in increasing order, until the index 
		* corresponding to the last bit set to 1 is reached.
		*
		* A subsequent call to this function when the index corresponding to the last bit set to 1
		* has been reached will result in 0 being returned.
		*
		* @memberof BitSet_.iterator
		* @return {number} a natural integer corresponding to the index of the computed bit set to 1 
		* or 0 in case all the bits set to 1 have already been iterated over.
		*/
		this.next = function() {
			// If the end of the bit set has already been reached, there is nothing more to do
			if (this.w_idx == that.words.length) {
				return 0;
			}			
			
			// Otherwise, extract the next bit set to 1 and compute its associated index, from the current
			// remaining word.
			var t = this.w & -this.w;
			var value = (this.w_idx << that.WORD_LOG) + BitSet_.populationCount((t - 1) | 0);
			this.w ^= t;
			
			// In case the current word has been exhausted, next iteration needs to
			// take place on the next non-null word, if existing.
			while (this.w_idx < that.words.length && this.w == 0) {
				++this.w_idx;
				this.w = that.words[this.w_idx];
			}
		
			// If the end of the bit set is reached, no subsequent call to .next are necessary
			return value;
		}
	}
}

/**
* @function populationCount_
*
* @summary Return the number of 1-bits in a 32-bit word.
*
* @description This function computes the number of 1-bits in a 32-bit word,
* using the formula 5-2 of the chapter 5 of the reference.
*
* @see Warren, H. (2009), Hacker`s Delight, New York, NY: Addison-Wesley
* 
* @param {number} x a 32-bit word.
* @return {number} the number of bits set to 1 in the binary representation of x.
*
* @example
* populationCount_(5);
* // 2
*/
BitSet_.populationCount = function(x) {
	x = x - ((x >>> 1) & 0x55555555);
	x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
	x = (x + (x >>> 4)) & 0x0F0F0F0F;
	x = x + (x >>> 8);
	x = x + (x >>> 16);
	return x & 0x0000003F;
};

BitSet_.prototype = {
    constructor: BitSet_,

	/**
	* @function resize
	*
	* @summary Resize the bit set.
	*
	* @description This function resizes the bit set so that
	* the bit corresponding to an index is stored within the bit set.
	*
	* @memberof BitSet_
	* @param {number} idx the index of the bit to be stored within the bit set, 
	* a natural integer.
	*
	* @example
	* resize_(123);
	* // 
	*/
	resize : function(idx) {
	    // Short circuit in case there is nothing to do
		var c = this.words.length;
		if ((c << this.WORD_LOG) > idx) {
			return;
		}
		
		// Compute the total number of words needed in order to
		// store the bit corresponding to the index provided in input.
	    var count = (idx + this.WORD_LENGTH) >>> this.WORD_LOG;
		
		// Depending on whether typed arrays are supported by the JavaScript
		// engine, resize the bit set differently.
		if (typeof Uint32Array === 'function') {
			var words_new = new Uint32Array(count); // the new typed array is (automatically) initialized with 0s
			words_new.set(this.words); // copy the previous words into the new typed array
			this.words = words_new;
		}
		else {  
			// Expand the array
			this.words[count-1] = 0; // copy (automatically) the previous words into the new array
			
			// Fill the expanded array with 0s
			for (var i = c; i < count-1; ++i) {
				this.words[i] = 0;
			}
		}
	},
	
	/**
	* @function toString
	*
	* @summary Return a string representation of the bit set as 0s and 1s.
	*
	* @description This function builds a base-2 string representation of
	* the bit set.
	* 
	* @memberof BitSet_
	* @return {string} a base-2 string representation of the bit set' content.
	*
	* @example
	* BitSet_().add(5).toString()
	* // 00000000000000000000000000000101
	*/
	toString: function () {
		// Initialize the final string
		var fullStr = '';

		// Concatenate the base-2 representation of each word hold by the bit set, 
		// possibly padded with leading WORD_LENGTH 0s.
		var c = this.words.length;
		for (var i = 0; i < c; ++i) {
			// Compute the base-2 string representation of the i-th word of the bit set.
			//
			// Note: if the underlying array of words is a standard array, words greater than
			// 2^(this.WORD_LENGTH-1)-1 will be considered as negative by the toString(2)
			// method below, so that the unsigned right shift bitwise operator (>>>) is
			// used to coerce the word to an unsigned integer.
			//
			// C.f. https://stackoverflow.com/questions/9939760/how-do-i-convert-an-integer-to-binary-in-javascript
			var str = "";
			if (typeof Uint32Array === 'function') {
				str = this.words[i].toString(2);
			}
			else {
				str = (this.words[i] >>> 0).toString(2);
			}

			// Concatenate the (possibly) padded string above with the other words
			// already built.
			fullStr += str.length >= this.WORD_LENGTH ? str : new Array(this.WORD_LENGTH - str.length + 1).join('0') + str;
		}

		// Return the computed string
		return fullStr;
	},

	/**
	* @function set
	*
	* @summary Set a single bit to 1.
	*
	* @description This function sets the bit corresponding to an index to 1.
	* 
	* @memberof BitSet_
	* @param {number} idx the index of the bit to set to 1, a natural integer.
	* @return {BitSet_} this.
	*
	* @example
	* BitSet_().set(5)
	* //
	*/
	set: function(idx) {
		// Logic to transparently grow the bit set
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
		
		// Set the proper bit to 1, in the proper word
		this.words[idx >>> this.WORD_LOG] |= (1 << idx);
		
		// Return the altered bit set
		return this;
	},

	/**
	* @function setRange
	*
	* @summary Set a continuous range of bits to 1.
	*
	* @description This function sets the bits within a continuous range of indexes to 1.
	* 
	* @memberof BitSet_
	* @param {number} idxFrom the index of the first bit to set to 1, a natural integer.
	* @param {number} idxTo the index of the last bit to set to 1, a natural integer greater than or equal to idxFrom.
	* @return {BitSet_} this.
	*
	* @example
	* BitSet_().setRange(5, 10)
	* //
	*/
	setRange: function (idxFrom, idxTo) {
		// Logic to transparently grow the bit set
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idxTo) {
			this.resize(idxTo);
		}
	  
		// Set the proper bits to 1, in the proper words
		for (var i = idxFrom; i <= idxTo; ++i) {
			this.words[i >>> this.WORD_LOG] |= (1 << i);
		}

		// Return the altered bit set
		return this;
	},
	
	/**
	* @function unset
	*
	* @summary Set a single bit to 0.
	*
	* @description This function sets the bit corresponding to an index to 0.
	* 
	* @memberof BitSet_
	* @param {number} idx the index of the bit to set to 0, a natural integer.
	* @return {BitSet_} this.
	*
	* @example
	* BitSet_().unset(5)
	* //
	*/
	unset: function(idx) {
		// Logic to transparently grow the bit set
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
		
		// Set the proper bit to 0, in the proper word
		this.words[idx >>> this.WORD_LOG] &= ~(1 << idx);
		
		// Return the altered bit set
		return this;
	},

	/**
	* @function get
	*
	* @summary Return the bit corresponding to an index.
	*
	* @description This function returns the bit corresponding to an index.
	* 
	* @memberof BitSet_
	* @param {number} idx the index of the bit to retrieve, a natural integer.
	* @return {boolean} true if the bit corresponding to the index idx is set to true,
	* false if the bit corresponding to the index idx does not exist in the bit set 
	* or is set to false.
	*
	* @example
	* BitSet_().set(5).get(5)
	* // true
	*/
	get: function(idx) {
		return (this.words[idx  >>> this.WORD_LOG] & (1 << idx)) !== 0;
	},	
	
	/**
	* @function clear
	*
	* @summary Clear the bit set.
	*
	* @description This function clears the bit set by resetting it to 0.
	* 
	* @memberof BitSet_
	* @return {BitSet_} this.
	*
	* @example
	* BitSet_().clear()
	* //
	*/
	clear: function() {
		// Re-initialize the bit set
		this.words = typeof Uint32Array === 'function' ? new Uint32Array(0) : new Array(0);
		
		// Return the altered bit set
		return this;
	},
	
	/**
	* @function flip
	*
	* @summary Flip/toggle the value of a single bit.
	*
	* @description This function flips/toggles the value of the bit corresponding to
	* an index.
	* 
	* @memberof BitSet_
	* @param {number} idx the index of the bit to flipt/toggle, a natural integer.
	* @return {BitSet_} this.
	*
	* @example
	* BitSet_().flip(5)
	* //
	*/
	flip: function(idx) {
		// Logic to transparently grow the bit set
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
		
	  // Set the proper bit to its boolean negation, in the proper word
	  this.words[idx >>> this.WORD_LOG] ^= 1 << idx;
	  
	  // Return the altered bit set
	  return this;
	},
	 
	/**
	* @function isEmpty
	*
	* @summary Determine whether the bit set is empty.
	*
	* @description This function determines whether the bit set is empty
	* (i.e., has no bit set to 1).
	* 
	* @memberof BitSet_
	* @return {boolean} true if the bit set is empty, false otherwise.
	*
	* @example
	* BitSet_().set(5).isEmpty()
	* // false
	*/
	isEmpty: function() {
	  // Loop over all the words help by the bit set to detetmine
	  // if there is a non-null word.
	  var c = this.words.length;
	  for (var  i = 0; i < c; ++i) {
		if (this.words[i] != 0) {
			return false;
		}
	  }
	  
	  // Arrived here, the bit set has only null words (if any), so that
	  // it is empty.
	  return true;
	},

	/**
	* @function nbSetBits
	*
	* @summary Compute the number of bits set to 1.
	*
	* @description This function computes the number of bits set to 1 in the bit set.
	* 
	* @memberof BitSet_
	* @return {number} the number of bits set to 1 in the bit set, a natural integer.
	*
	* @example
	* BitSet_().set(5).nbSetBits()
	* // 1
	*/
	nbSetBits: function() {
	  // Loop over all the words help by the bit set and compute their
	  // number of bits set to 1 thanks to the populationCount function.
	  var s = 0;
	  var c = this.words.length;
	  for (var i = 0; i < c; ++i) {
		s += BitSet_.populationCount(this.words[i] | 0);
	  }
	  
	  // Return the computed value
	  return s;
	},

	/**
	* @function toArray
	*
	* @summary Return an array representation of the bit set.
	*
	* @description This function builds an array representation of the bit set,
	* which contains the indexes of the bits set to 1 in increasing order.
	* 
	* @memberof BitSet_
	* @return {Uint32Array<number>} an array representation of the bit set' content, 
	* a newly allocated array of length the number of bits set to 1 in the
	* bit set.
	*
	* @example
	* BitSet_().add(5).add(10).toArray()
	* // [5, 10]
	*/
	toArray: function() {
	  // Initialize the output array, and its associated elements pointer
	  var arr = new Array(this.nbSetBits());
	  var pos = 0;
	  
	  // Loop over all the words help by the bit set 
	  var c = this.words.length;
	  for (var i = 0; i < c; ++i) {
		var w = this.words[i];
		
		// For each word help by the bit set, extract the bits set to 1
		// and compute their associated index.
		while (w != 0) {
		  var t = w & -w;
		  arr[pos++] = (i << this.WORD_LOG) + BitSet_.populationCount((t - 1) | 0);
		  w ^= t;
		}
	  }
	  
	  // Return the computed array
	  return arr;
	},

};
