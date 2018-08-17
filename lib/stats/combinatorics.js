/**
 * @file Misc. combinatorics functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.aliasMethodSampler_ = aliasMethodSampler_;
self.compositionsIterator_ = compositionsIterator_;
self.subsetsIterator_ = subsetsIterator_;
self.kSubsetsIterator_ = kSubsetsIterator_;
self.randomKSubsetIterator_ = randomKSubsetIterator_;
self.binomial_ = binomial_;
/* End Wrapper private methods - Unit tests usage only */
 

/**
* @function aliasMethodSampler_
*
* @summary Returns a function to generate random values sampled from a discrete
* finite probability distribution.
*
* @description This function constructs a function to generate random values from the set
* {0,...,n-1} sampled according to the provided discrete finite probability distribution
* {p_0,...,p_n-1}.
* 
* The algorithm used is the Vose's algorithm, which is a variation of the alias method
* allowing to sample random values from a finite discrete probability distribution in O(1) time
* after a O(n) time preprocessing step, c.f. the reference.
* 
* @see <a href="https://doi.org/10.1109/32.92917">M. D. Vose, A linear algorithm for generating random numbers 
* with a given distribution, IEEE Transactions on Software Engineering, vol. 17, no. 9, pp. 972-975, Sep 1991.</a>
*
* @param {Array.<number>} p, an array of n positive real numbers p_0,...,p_n-1 with sum_i p_i = 1.
* @return {function} a function to be used through its .sample() method, generating an integer i from the set
* {0,...,n-1} with probability p_i.
*
* @example
* var mySampler = new aliasMethodSampler_([0, 0.1, 0.4, 0.5]);
* mySampler.sample();
* // 3;
*/
function aliasMethodSampler_(p) {
	// ----
	// init function, c.f. paragraph B of section III of the reference.
	// ----

	// Initializations.
	this.prob = typeof Float64Array === 'function' ? new Float64Array(p.length) : new Array(p.length);
	this.alias = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);

	// TODO: Checks on probabilities (positive, sum to one)

	// Computation of the average probability.
    var avgProb = 1 / p.length;
		 
	// Initializations of the small and large stacks, together with their associated indexes.
	var small = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var s = 0;
	var large = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var l = 0;
		
	// Population of the small and large stacks with the probabilities indexes.
	for (var j = 0; j < p.length; ++j) {
		if (p[j] > avgProb) {
			large[l] = j;
			++l;
		}
		else {
			small[s] = j;
			++s;
		}
	}
		
	// Main loop of the algorithm, populating the prob and alias arrays.
	var p = p.slice(0); // local copy of the probabilities, as they are updated below
	while (s > 0 && l > 0) {
		// Get the index of the small and the large probabilities.
		--s;
		var j = small[s];
		
		--l;
		var k = large[l];
		
		// Update the prob and alias arrays.
		this.prob[j] = p.length * p[j];
		this.alias[j] = k;
		
		// Update the probabilities.
		p[k] = p[k] + (p[j] - avgProb);
		
		// Update the large and small stacks.
		if (p[k] > avgProb) {
			large[l] = k;
			++l;
		}
		else {
			small[s]= k;
			++s;
		}
	}
		
	// Process the remaining elements of the small stack.
	while (s > 0) {
		--s;
		this.prob[small[s]] = 1;
	}
	
	// Process the remaining elements of the large stack.
	//
	// Theoretically not needed, but due to round off errors, practicaly needed.
	while (l > 0) {
		--l;
		this.prob[large[l]] = 1;
	}
	
	
	// ----
	// rand function, c.f. paragraph A of section III of the reference.
	// ----
	
	/**
	* @function sample
	*
	* @summary Returns a random value sampled from the underlying probability distribution.
	*
	* @description This function computes a random value sampled from the underlying 
	* probability distribution using the method described in the reference.
	*
	* @memberof aliasMethodSampler_
	* @return {number} an integer i belonging to the set {0,...,n-1} with probability p_i.
	*/
    this.sample = function() {
		var u = Math.random() * this.prob.length; // Uniform real number belonging to [0..n[
		var j = Math.floor(u);
		if (u - j <= this.prob[j]) {
			return j;
		}
		else {
			return this.alias[j];
		}
    }
}

	
 /**
* @function compositionsIterator_
*
* @summary Returns an iterator to compute all the compositions of a non negative integer.
*
* @description This function constructs an iterator to compute all the k-compositions of a non-negative integer n, 
* using the algorithm NEXCOM described in section 5 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n a non-negative integer whose composition are desired.
* @param {number} k a non-negative integer indicating the number of parts of desired composition of n.
* @return {function} a function to be used as an iterator through its .next() method, computing all 
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
	* @description This function computes the next k-composition of a non negative integer n.
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
			++this.h;
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
* @function kSubsetsIterator_
*
* @summary Returns an iterator to compute all the k-subsets of a n-set.
*
* @description This function constructs an iterator to compute all the k-subsets of the n-set {1,...,n}, 
* using the algorithm NEXKSB described in section 3 of the first reference.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="https://en.wikipedia.org/wiki/Power_set">Power set</a>
*
* @param {number} n the number of elements of the n-set {1,...,n} whose k-subsets are desired, a non-negative integer.
* @param {number} k a non-negative integer, with 1 <= k <= n.
* @return {function} a function to be used as an iterator through its .next() method, computing all the 
* k-subsets of the n-set {1,...,n} in lexicographic order.
*
* @example
* var myIterator = new kSubsetsIterator_(5, 3);
* myIterator.next(); myIterator.next();
* // [1, 2, 3]; [1, 2, 4]; ...; -1
*/
function kSubsetsIterator_(n, k) {
	// Initialize n and k
	this.n = n;
	this.k = k;
	
	// Initialize the array to hold the k-subsets
	this.a = new Array(k);
	
	// Variables required for NEXKSB internal computations,
	// initialized so as to generate the first subset upon
	// the first call to .next() function.
	this.firstcall = true;
	this.mtc = false;
	this.m2 = 0;
	this.h = this.k;
	this.endval = this.n - this.k + 1;	
	
	/**
	* @function next
	*
	* @summary Returns the next k-subset of a set.
	*
	* @description This function computes the next k-subset of the n-set {1,...,n},
	* in lexicographic order.
	*
	* The initial k-subset computed by the first call to this function is the subset {1,...,k}, 
	* and each subsequent call to this function will result in a new k-subset until the final 
	* k-subset {n-k+1,...,n} is reached.
	*
	* A subsequent call to this function when the final k-subset has been reached will result in
	* the function returning -1.
	*
	* @memberof kSubsetsIterator_
	* @return {Array.<number>|number} either an array containing a newly generated sorted k-subset
	* of the n-set {1,...,n} or -1 to indicate that all the k-subsets have already been generated.
	*/
	this.next = function() {
		if (!this.firstcall && !this.mtc) {
			// No more k-subset to generate
			return -1;
		}
		
		if (this.firstcall) {		
			// The first call has now been made
			this.firstcall = false;
		}
		else {
			// There is still a k-subset to generate
			if (this.m2 < this.n - this.h) {
				this.h = 0;
			}
			
			++this.h;
			this.m2 = this.a[this.k - this.h];
		}
				
		// Fill the k-subset array with the initial values {1,...,k}
		for (var j = 1; j <= this.h; ++j) {
			this.a[this.k + j - this.h - 1] = this.m2 + j;
		}

		// End logic
		this.mtc = (this.a[0] != this.endval);

		// Return a copy of the computed array
		return this.a.slice();
	}
}

/**
* @function randomKSubsetIterator_
*
* @summary Returns an infinite iterator to compute random k-subsets of a n-set.
*
* @description This function constructs an iterator to compute random k-subsets of the n-set {1,...,n}, 
* as lists of k distinct increasing integers in {1,...,n}, using the method D of the references.
*
* From the references, the random k-subsets are generated uniformly at random, in O(k) time 
* and O(1) additional space.
*
* @see J.S. Vitter. Faster Methods for Random Sampling. Communications of the ACM, 27, (July 1984), 703-718
* @see J.S. Vitter. An efficient algorithm for sequential random sampling. RR-0624, INRIA. 1987. <inria-00075929>
*
* @param {number} n the number of elements of the n-set {1,...,n} whose k-subsets are desired, a non-negative integer.
* @param {number} k a non-negative integer, with 1 <= k <= n.
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
	
	// Initialize an array to hold the k-subsets
	this.a = null;

	// Initializations for the method D
	// - N, the number of records that have not yet been processed
	// - nn, the number of records remaining to be selected
	// - idx_a, the array index used to write the selected records in the array a 
	// - selected_record, the value of the selected record 
	this.N;
	this.nn;
	this.idx_a;
	this.selected_record;
	
	
	/**
	* @function uniformrv
	*
	* @summary Returns a number generated uniformly at random in interval ]0,1[.
	*
	* @description This function computes a number uniformly at random in interval ]0,1[.
	*
	* @memberof randomKSubsetIterator_
	* @return {number} a number generated uniformly at random in interval ]0,1[, a real number
	*
	*/
	function uniformrv() {
		// Generate a random number in the [0, 1[ interval
		var rnd = Math.random();
		
		// While the generated random number is (exactly) equal to 0,
		// reject it.
		while (rnd === 0) {
			rnd = Math.random();
		}
		
		// Return the generated random number, which is then 
		// generated uniformly at random in the ]0,1[ interval.
		return rnd;
	}

	/**
	* @function method_a
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the method A of the references,
	* and is used by the method D to avoid its worst-case behaviour, c.f. the references.
	*
	* @memberof randomKSubsetIterator_
	*
	*/
	this.method_a = function() {
		// Initializations
		var top = this.N - this.nn;
		
		// General case
		while (this.nn >= 2) {
			// Step A1
			var V = uniformrv();
			
			// Step A2
			var S = 0;
			var quot = top / this.N;
			
			while (quot > V) {
				++S;
				--top;
				--this.N;
				quot = (quot * top) / this.N;
			}
			
			// Step A3
			// Skip over the next S records and select the following one for the sample
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
			
			--this.N;
			--this.nn;
		}
		
		// Special case nn = 1
		var S = Math.floor(this.N * uniformrv());
		if (S === this.N) { // the out of range value S = N must never be generated (which could be in case of roundoff errors)
			S = this.N - 1;
		}
		
		// Skip over the next S records and select the following one for the sample
		this.selected_record += S + 1;
		this.a[this.idx_a++] = this.selected_record;
	}

	/**
	* @function method_d
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the method D of the references.
	*
	* @memberof randomKSubsetIterator_
	*
	*/
	this.method_d = function() {
		// Initializations
		var ninv = 1/this.nn;
		var Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
		var qu1 = -this.nn + 1 + this.N;
		var negalphainv = -13;
		var threshold = -negalphainv * this.nn;
		
		while (this.nn > 1 && threshold < this.N) {
			var nmin1inv = 1 / (-1 + this.nn);
			
			var X;
			var S;
			while (true) {
				// Step D2: Generate U and X
				while(true) {
					X = this.N 	* (-Vprime + 1);
					S = Math.floor(X);
					if (S < qu1) {
						break;
					}
					Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
				}
				var U = uniformrv();
				
				// Step D3: Accept ?
				var y1 = Math.pow(U * this.N / qu1, nmin1inv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
				Vprime = y1 * (-X/this.N + 1) * (qu1 / (-S + qu1));
				if (Vprime <= 1) {
					break;
				}
				
				// Step D4: Accept ?
				var y2 = 1;
				var top = -1 + this.N;
				
				var bottom;
				var limit;
				if (-1 + this.nn > S) {
					bottom = -this.nn + this.N;
					limit = -S + this.N;
				}
				else {
					bottom = -1 - S + this.N;
					limit = qu1;
				}
				
				for (var t = -1 + this.N; t >= limit; --t) {
					y2 = y2 * top / bottom;
					--top;
					--bottom;
				}
				
				if (this.N / (-X + this.N) >= y1 * Math.pow(y2, nmin1inv)) { // Math.exp(Math.log(a) * b) === Math.pow(a, b)
					// Accept
					Vprime = Math.pow(uniformrv(), nmin1inv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
					break;
				}
				Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
			}
			
			// Step D5: Select the (S + 1)st record
			// Skip over the next S records and select the following one for the sample
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
			
			// Prepare for the next iteration
			this.N = -S + (-1 + this.N);
			--this.nn;
			ninv = nmin1inv;
			qu1 = -S + qu1;
			threshold = threshold + negalphainv;
		}
		
		// If nn > 1 (i.e., threshold < N), use method A to finish the sampling,
		// otherwise, if nn == 1, deal with the special case nn = 1.
		if (this.nn > 1) {
			this.method_a();
		}
		else {
			var S = Math.floor(this.N * Vprime);
			if (S === this.N) { // the out of range value S = N must never be generated (which could be in case of roundoff errors)
				S = this.N - 1;
			}
		
			// Skip over the next S records and select the following one for the sample
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
		}		
	}
   
	/**
	* @function next
	*
	* @summary Returns a random k-subset of a n-set.
	*
	* @description This function computes a random k-subset of a n-set using the method D of the references.
	*
	* @memberof randomKSubsetIterator_
	* @return {Array.<number>|Uint32Array} a random k-subset of the n-set {1,...,n}, a sorted array of k increasing strictly positive integers.
	*
	*/
	this.next = function() {
		// Initializations
		// Array holding the k-subset
		if (typeof Uint32Array === 'function') {
			this.a = new Uint32Array(this.k);
		}
		else {
			this.a = new Array(this.k);
			for (var i = 0; i < this.k; ++i) {
				this.a[i] = 0;
			}
		}	

		// Misc. internal variables required by the method D
		this.N = this.n;
		this.nn = this.k;
		this.idx_a = 0;
		this.selected_record = 0;
			
		// Call the method D, which will proceed with the effective
		// generation of the k-subset.
		this.method_d();
		
		// Return the computed array, which is not used anymore
		return this.a;
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
* @return {function} a function to be used as an iterator through its .next() method, computing all 
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
		var nextSubset = new Array(0);
		
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
     