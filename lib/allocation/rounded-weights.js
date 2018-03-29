/**
* @file Functions related to (rational) rounding of floating-point portfolio weights.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function roundedWeights
*
* @summary Compute the closest rational approximation of the weights of a portfolio.
*
* @description Given n (floating-point) weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, 
* this function returns n (rational) weights wr_1,...,wr_n associated to a fully invested and long-only portfolio of n assets
* satisfying:
* - k * wr_i is a natural integer, i=1..n
* - wr_1,...,wr_n are the closest weights to w_1,...,w_n, in the sense defined in the reference.
*
* To be noted that typical values of k are 10 (rounding to 10%), 20 (rounding to 5%) and 100 (rounding to 1%).
*
* @see <a href="https://doi.org/10.1007/s10898-013-0126-2">.M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258.</a>
* 
* @param {Array.<number>} originalWeights the weights w_1,...,w_n associated to a fully invested and long-only portfolio of n assets, array of n real numbers.
* @param {number} k the value to which the rounded weights will be a multiple of the inverse, natural integer greater than or equal to 1.
* @return {Array.<number>} the rounded weights wr_1,...,wr_n, array of n real numbers.
*
* @example
* roundedWeights([0.5759, 0.0671, 0.3570], 10);
* // [0.6, 0.1, 0.3]
* roundedWeights([0.5759, 0.0671, 0.3570], 20);
* // [0.6, 0.05, 0.35]
* roundedWeights([0.5759, 0.0671, 0.3570], 100);
* // [0.57, 0.07, 0.36]
*/
self.roundedWeights = function (originalWeights, k) {
	// ------
	
	// Call to the simplex rational rounding method
	var roundedWeights = simplexRationalRounding_(originalWeights, k);

	// Return the computed weights
	return roundedWeights;	
}
