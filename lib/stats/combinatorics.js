/**
 * @file Misc. combinatorics functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.nexcom_ = function(n, k) { return nexcom_(n, k); }
self.binomial_ = function(n, k) { return binomial_(n, k); }
/* End Wrapper private methods - Unit tests usage only */
 
 
/**
* @function nexcom_
*
* @summary Returns the next composition of a non negative integer.
*
* @description This function computes the next k-composition of a non-negative integer n, 
* using the algorithm NEXCOM described in section 5 of the first reference.
*
* The initial k-composition computed by the function is n00...0, and each subsequent call to the function
* will result in a new k-composition until the final k-composition 00...0n is reached.
*
* A subsequent call to the function when the final k-composition 00...0n has been reached will result in
* the recomputation of all the k-compositions of n, starting from the initial k-composition n00...0.
*
* @see Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators. 2d ed. New York: Academic Press.
* @see <a href="*https://en.wikipedia.org/wiki/Composition_(combinatorics)">Composition (combinatorics)</a>
*
* @param {number} n a non-negative integer whose composition are desired.
* @param {number} k a non-negative integer indicating the number of parts of desired composition of n.
* @return {Array} an array arr of 2 elements, with arr[0] a boolean indicating whether at least one k-composition of n remains to be computed
* and arr[1] an array of k elements containing the computed k-composition of n.
*
* @example
* nexcom_(6, 3); nexcom_(6, 3);
* // [true, [6,0,0]]; [true, [5,1,0]];
*/
function nexcom_(n, k) {
	// Variable mtc has been statified for ease of use reasons
	if (nexcom_.mtc === undefined) {
		nexcom_.mtc = false;
	}
	
	// NEXTCOM computation logic
	if (nexcom_.mtc) { // There is still a composition to generate
		if (nexcom_.t > 1) {
			nexcom_.h = 0;
		}
		nexcom_.h++;
		nexcom_.t = nexcom_.r[nexcom_.h - 1];
		nexcom_.r[nexcom_.h - 1] = 0;
		nexcom_.r[0] = nexcom_.t - 1;
		nexcom_.r[nexcom_.h]++;		
	}
	else  { // No more composition to generate, so, (re) generation of the first composition, equals to n00...0
		nexcom_.r = new Array(k);
		nexcom_.r[0] = n;

		// Variables h and t need to be static 
		nexcom_.t = n;
		nexcom_.h = 0;
		
		for (var i = 1; i <= k-1; ++i) {
			nexcom_.r[i] = 0;
		}
	}
	
	// NEXTCOM end logic
	nexcom_.mtc = (nexcom_.r[k-1] != n);
		
	// Return the value of the mtc variable, so that calls to nexcom_ function
	// can be chained in a while loop
	// Return a copy of the r array, so that callers can alter it
	return [nexcom_.mtc, nexcom_.r.slice()];
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
     