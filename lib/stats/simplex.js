/**
* @file Functions related to computations on the unit simplex.
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/


/* Start Wrapper private methods - Unit tests usage only */
self.closestRationalPoint_ = function(x, r) { return closestRationalPoint_(x, r); }
/* End Wrapper private methods - Unit tests usage only */


/**
* @function closestRationalPoint
*
* @summary Compute a rational point on the unit simplex that is the closest to a point on the unit simplex.
*
* @description Given a point x = (x_1,...,x_n) on the standard simplex of R^n, which is the n-1 dimensional set of R^n containing 
* the points y = (y_1,...,y_n) statisfying sum y_i = 1 and y_i >= 0, i = 1..n, this function computes a proximal point xr = (xr_1,...,xr_n) on 
* the r-th rational grid of the unit simplex of R^n, 1/r * I_n(r), with I_n(r) the set of N^n containing the points m = (m_1,...,m_n) 
* statisfying sum m_i = r with r a strictly positive natural integer, so that the computed proximal point xr is one of the closest points to x 
* on this grid with respect to any norm in a large class, c.f. the reference.
*
* @see <a href="https://link.springer.com/article/10.1007/s10898-013-0126-2">.M. Bomze, S. Gollowitzer, and E.A. Yıldırım, Rounding on the standard simplex: Regular grids for global optimization, J. Global Optim. 59 (2014), pp. 243–258.</a>
* 
* @param {Array.<number>} x a point belonging to the standard simplex of R^n.
* @param {number} r the indice of the rational grid of the unit simplex of R^n, 1/r * I_n(r), on which to compute the closest point to x.
* @return {Array.<number>} the computed closest point to x on the r-th rational grid of the unit simplex of R^n.
*
* @example
* closestRationalPoint_([], 20);
* // XX
*/
function closestRationalPoint_(x, r) {
	// TODO: Checks, if enabled

	// Compute the integer and fractional parts of the coordinates of the input point multiplied by r, as described in paragraph 2 of the reference.
	// In parallel, compute k as also defined in paragraph 2 of the reference. 
	var k = r;
	var xPartsWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		//
		var rx = r * x[i];
		var integerPart = Math.floor(rx);
		var fractionalPart = rx - integerPart;
		
		//
		k -= integerPart;
		
		//
		xPartsWithIndexes[i] = [integerPart, fractionalPart, i];
	}

	// Re-order the fractional parts in non-increasing order, c.f. theorem 1 of the reference.
	xPartsWithIndexes.sort(function(a, b) {
		return b[1] - a[1];
	}); 

	// Rounding rule: round the k largest fractional parts up to one and all other fractional parts down to zero,
	// as described in paragraph 2 of the reference.
	var xr = new Array(x.length);
	for (var i = 0; i < k; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = (xPartsWithIndexes[i][0] + 1) / r;
	}
	for (var i = k; i < xPartsWithIndexes.length; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = xPartsWithIndexes[i][0] / r;
	}

	// Return the computed point
	return xr;
}
