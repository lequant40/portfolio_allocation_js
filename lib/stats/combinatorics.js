/**
 * @file Misc. combinatorics functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.compositionsIterator_ = compositionsIterator_;
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
	this.k = k
	
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
	* A subsequent call to this function when the final k-composition 00...0n has been reached will result in
	* the recomputation of all the k-compositions of n, starting from the initial k-composition n00...0.
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
			this.r[this.h]++;		
		}
		else  { 
		    // No more composition to generate, so, (re) generation of the first composition, equals to n00...0
			this.r[0] = this.n;
			for (var i = 1; i <= this.k-1; ++i) {
				this.r[i] = 0;
			}
		}
		
		// End logic
		this.mtc = (this.r[this.k-1] != this.n);
		
		// Return a copy of the r array, so that callers can alter it
		return [this.mtc, this.r.slice()];
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
     