/**
 * @file Misc. combinatorics functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.compositionsIterator_ = compositionsIterator_;
self.subsetsIterator_ = subsetsIterator_;
self.binomial_ = function(n, k) { return binomial_(n, k); }
/* End Wrapper private methods - Unit tests usage only */
 
 
/**
* @function compositionsIterator_
*
* @summary Returns an iterator to compute all the compositions of a non negative integer.
*
* @description This function constructs an iterator to compute all the k-compositions of a non-negative integer n, 
* using the algorithm NEXCOM described in section 5 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="*https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n a non-negative integer whose composition are desired.
* @param {number} k a non-negative integer indicating the number of parts of desired composition of n.
* @return {function} a function to be used as an iterator through its .next() method, sequentially computing all 
* the k-compositions of n.
*
* @example
* var myIterator = new compositionsIterator_(6, 3);
* myIterator.next(); myIterator.next();
* // [true, [6,0,0]]; [true, [5,1,0]];
*/
function compositionsIterator_(n, k) {
	// Initialize n and k
	this.n = n;
	this.k = k;
	
	// Variables required for NEXTCOM internal computations,
	// initialized so as to generate the first composition upon
	// the first call to .next() function.
	this.mtc = false;
	this.r = new Array(k);
	this.t = this.n;
	this.h = 0;

	/**
	* @function next
	*
	* @summary Returns the next composition of a non negative integer.
	*
	* @description This function computes the next k-composition of n.
	*
	* The initial k-composition computed by the first call to this function is n00...0, and each subsequent call to 
	* this function will result in a new k-composition until the final k-composition 00...0n is reached.
	*
	* A subsequent call to this function when the final k-composition has been reached will result in
	* the recomputation of all the k-compositions of n, starting from the initial k-composition.
	*
	* @memberof compositionsIterator_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one k-composition of n remains to be computed
	* and arr[1] an array of k elements containing the computed k-composition of n.
	*
	*/
	this.next = function() {
		if (this.mtc) { // There is still a composition to generate
			if (this.t > 1) {
				this.h = 0;
			}
			this.h++;
			this.t = this.r[this.h - 1];
			this.r[this.h - 1] = 0;
			this.r[0] = this.t - 1;
			++this.r[this.h];
		}
		else  { 
		    // No more composition to generate, so, (re) generation of the first composition, equals to n00...0
			this.r[0] = this.n;
			for (var i = 1; i <= this.k - 1; ++i) {
				this.r[i] = 0;
			}
		}
		
		// End logic
		this.mtc = (this.r[this.k - 1] != this.n);
		
		// Return a copy of the r array, so that callers can alter it
		return [this.mtc, this.r.slice()];
	}
}


/**
* @function randomKSubsetIterator_
*
* @summary Returns an infinite iterator to compute random k-subsets of a n-set.
*
* @description This function constructs an iterator to compute random k-subsets of the n-set {1,...,n}, 
* using both the algorithms RANKSB and RKS2 described in section 4 of the reference.
*
* From the discussion following the examples in the reference, the random k-subsets are probably generated
* uniformly, but this is not written in the reference.
*
* The algorithm used to compute the random k-subsets is either RANKSB when k < n/2,
* or RKS2 when k >= n/2, so that linear performances in O(k) are guaranteed.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
*
* @param {number} n the number of elements of the n-set {1,...,n} whose k-subsets are desired, a non-negative integer.
* @param {number} k a non-negative integer, with 0 <= k <= n.
* @return {function} a function to be used as an iterator through its .next() method, computing random 
* k-subsets of the n-set {1,...,n}.
*
* @example
* var myIterator = new randomKSubsetIterator_(6, 3);
* myIterator.next();
* // [1, 2, 5];
*/
function randomKSubsetIterator_(n, k) {
	// Initialize n and k
	this.n = n;
	this.k = k;
	
	// Initialize the array to hold the k-subsets
	this.a = new Array(k);

	/**
	* @function next_ranksb
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the algorithm RANKSB of the reference.
	*
	* @memberof randomKSubsetIterator_
	* @return {Array.<number>} a random k-subset of the n-set {1,...,n}, a sorted array of k increasing strictly positive integers.
	*
	*/
	this.next_ranksb = function() {
		// Step A - Initialization of a
		for (var i = 1; i <= this.k; ++i) {
			this.a[i-1] = Math.floor((i - 1) * this.n / this.k);
		}
		
		// Step B
		// Note: in the reference, the c variable is initialized to k and is decremented until 0
		// each time a generated x is accepted: this is a reverse for loop in disguise.
		var x;
		var l;
		for (var c = this.k; c > 0; --c) {
			do {
				var u = Math.random();
				x = 1 + Math.floor(u * this.n);
				l = 1 + Math.floor((x * this.k - 1) / this.n);
			} while (x <= this.a[l-1]);
			this.a[l-1] = this.a[l-1] + 1;
		}
		var p = 0;
		var s = this.k;
		
		// Step C
		// Note: in the reference, the i variable is initialized to 0 and is incremented
		// until k each time: this is a for loop in disguise.
		for (var i = 1; i <= this.k; ++i) {
			if (this.a[i-1] == Math.floor((i - 1) * this.n / this.k)) {
				this.a[i-1] = 0;
			}
			else {
				p = p + 1;
				var m = this.a[i-1];
				this.a[i-1] = 0;
				this.a[p-1] = m;
			}
		}
		
		// Step D
		// Note: in the reference, the p variable is initialized to whatever value it has, and is decremented
		// until 0 each time: this is a reverse for loop in disguise.
		for (; p > 0; --p) {
			l = 1 + Math.floor((this.a[p-1] * this.k - 1) / this.n);
			var delta_s = this.a[p-1] - Math.floor((l - 1) * this.n / this.k);
			this.a[p-1] = 0;
			this.a[s-1] = l;
			s = s - delta_s;			
		}
		l = k;
		
		// Steps E to H
		// Note: in the reference, the l variable is initialized at this step to k, and is decremented
		// until 0 each time: this is a reverse for loop in disguise.
		var r;
		for (; l > 0; --l) {
			// Step E
			var m_0;
			if (this.a[l-1] != 0) {
				r = l;
				m_0 = 1 + Math.floor((this.a[l-1] - 1) * this.n / this.k);
				m = Math.floor(this.a[l-1] * this.n / this.k) - m_0 + 1;
			}

			// Step F
			var u = Math.random();
			x = m_0 + Math.floor(u * m);
			i = l;
			
			// Step G
			++i;
			while (i <= r && x >= this.a[i-1]) {
				this.a[i-2] = this.a[i-1];
				x = x + 1;
				++i;
			}
			
			// Step H
			this.a[i-2] = x;
			m = m - 1;
		}
		
		// Return a copy of the computed array
		return this.a.slice();
	}
	
	/**
	* @function next_rks2
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the algorithm RKS2 of the reference.
	*
	* @memberof randomKSubsetIterator_
	* @return {Array.<number>} a random k-subset of the n-set {1,...,n}, a sorted array of k increasing strictly positive integers.
	*
	*/
	this.next_rks2 = function() {
		// Initializations
		var c_1 = this.k;
		var c_2 = this.n;
		var k_0 = 0;
		var i = 0;

		// Main loop of the RKS2 algorithm
		while (c_1 > 0) {
			++i;
			var u = Math.random();
			if (u <= c_1/c_2) {
				c_1 = c_1 - 1;
				this.a[k_0] = i;
				k_0 = k_0 + 1; // this line is inversed compared to the reference because of JavaScript arrays starting at index 0
			}
			c_2 = c_2 - 1;
		}
		
		// Return a copy of the computed array
		return this.a.slice();
	}
	
	// Initialize the appropriate iterator to keep the required labor to O(k) uniformly for 1 <= k <= n
	if (k < n/2) {
		this.next = this.next_ranksb;
	}
	else {
		this.next = this.next_rks2;
	}
}


/**
* @function subsetsIterator_
*
* @summary Returns an iterator to compute all the subsets of a n-set.
*
* @description This function constructs an iterator to compute all the subsets of the n-set {1,...,n}, 
* using the algorithm NEXSUB described in section 1 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="*https://en.wikipedia.org/wiki/Power_set">Power set</a>
*
* @param {number} n the number of elements of the n-set {1,...,n} whose subsets are desired, a non-negative integer.
* @return {function} a function to be used as an iterator through its .next() method, sequentially computing all 
* the subsets of the n-set {1,...,n}.
*
* @example
* var myIterator = new subsetsIterator_(5);
* myIterator.next(); myIterator.next();
* // [true, []]; [true, [1]];
*/
function subsetsIterator_(n) {
	// Initialize n
	this.n = n;
	
	// Variables required for NEXSUB internal computations,
	// initialized so as to generate the first subset upon
	// the first call to .next() function.
	this.mtc = false;
	this.iin = new Array(n);
	this.ncard = 0;
	
	/**
	* @function next
	*
	* @summary Returns the next subset of a set.
	*
	* @description This function computes the next subset of the n-set {1,...,n}.
	*
	* The initial subset computed by the first call to this function is the empty subset {}, and each subsequent call to 
	* this function will result in a new subset until the final subset {1,...,n} is reached.
	*
	* A subsequent call to this function when the final subset {1,...,n} has been reached will result in
	* the recomputation of all the subsets, re-starting from the initial subset.
	*
	* @memberof subsetsIterator_
	* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one subset
	* of the n-set {1,...,n} remains to be computed and arr[1] an array containing the computed sorted subset
	* of the n-set {1,...,n}.
	*
	*/
	this.next = function() {
		// The output array containing the computed subset
		var nextSubset = [];
		
		if (this.mtc) { // There is still a subset to generate
			var j = 0;
			if (this.ncard % 2 != 0) {
				++j;
				while (this.iin[j - 1] == 0) {
					++j;
				}
			}
			this.iin[j] = 1 - this.iin[j];
			this.ncard = this.ncard + 2*this.iin[j] - 1;

			// Build the output array
			nextSubset = new Array(this.ncard);
			var idx = 0;
			for (var i = 0; i <= this.n - 1; ++i) {
				if (this.iin[i] == 1) {
					nextSubset[idx++] = i + 1;
				}
			}
			
			// End logic
			this.mtc = (this.ncard != this.iin[this.n -1]);
		}
		else  { 
		    // No more subset to generate, so, (re) generation of the first subset, equals to {}
			for (var i = 0; i <= this.n - 1; ++i) {
				this.iin[i] = 0;
			}
			
			// The output array is already built in this case (empty)
			
			// Specific end logic
			this.mtc = true;
		}

		// Return the computed array, not used anymore by this function
		return [this.mtc, nextSubset];
	}
}


/**
* @function binomial_
*
* @summary Returns a binomial coefficient.
*
* @description This function computes the k-th binomial coefficient of order n, which is
* the coefficient of the x^k term in the polynomial expansion of the binomial power (1 + x)^n.
*
* This coefficient is also the number of ways to choose a subset of k elements,
* disregarding their order, from a set of n elements.
*
* The algorithm used is a multiplicative formula, c.f. the reference.
*
* @see <a href="*https://en.wikipedia.org/wiki/Binomial_coefficient">Binomial coefficient</a>
*
* @param {number} n a non-negative integer.
* @param {number} k a non-negative integer, with 0 <= k <= n.
* @return {number} the computed binomial coefficient.
*
* @example
* binomial_(7, 5);
* // 21
*/
function binomial_(n, k) {
    // Checks
    if (n < 0) {
        throw new Error('n must be a positive integer');
    }
    if (k < 0) {
        throw new Error('k must be a positive integer');
    }
    if (k > n) {
        throw new Error('k must be less than or equal to n');
    }
	  
	// Compute the binomial coefficient using the multiplicative formula of the reference.
	var val = 1;
    for (var i = 1; i <= k; ++i) {
        val *= (n + 1 - i);
		val /= i; // Note: separating the computation of val in two steps guarantees (unless compiler optimisations) that val is an integer.
    }
	
	// Return it
    return val;
}
     