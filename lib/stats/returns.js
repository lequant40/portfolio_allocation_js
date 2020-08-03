/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function returns
*
* @summary Compute the period-to-period returns of a series of values.
*
* @description This function returns the period-to-period returns of a series of values.
*
* The period-to-period returns of a series of values are defined as the series of the returns of the series of values
* over each of its valuation period, with the return associated to the first period being undefined.
*
* @see <a href="https://en.wikipedia.org/wiki/Rate_of_return">https://en.wikipedia.org/wiki/Rate_of_return</a>
* 
* @param {Array.<number>} x the series of values, an array of T real numbers, with T corresponding to the number of valuation periods.
* @param {object} opt optional parameters for the returns computation.
* @param {string} opt.method the method to use to compute the returns, a string either equals to:
* - "arithmetic", in order to compute the arithmetic returns, c.f. the reference
* - "logarithmic", in order to compute the logarithmic returns, c.f. the reference
; defaults to "arithmetic"

* @return {Array.<number>} the period-to-period returns of x.
*
* @example
* returns([1, 2, 1]); 
* // [1.0, -0.5], i.e. 100% arithmetic return from the first period to the second period, 
* // and -50% arithmetic return from the second period to the third period
*/
self.returns = function(x, opt) {
	// Initialize default parameters
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "arithmetic";
	}
	
	
	// Decode the parameters
	var method = opt.method;
	if (method != "arithmetic" &&
	    method != "logarithmic") {
			throw new Error('unsupported returns computation method');
	}
	
	
	// Compute the returns
	var returns = typeof Float64Array === 'function' ? new Float64Array(x.length - 1) : new Array(x.length - 1); 

	if (method == "arithmetic") {
		for (var i = 0; i < x.length - 1; ++i) {
			returns[i] = (x[i+1] - x[i])/x[i];
		}
	}
	else if (method == "logarithmic") {
		for (var i = 0; i < x.length - 1; ++i) {
			returns[i] = Math.log(x[i+1]/x[i]);
		}
	}
	else {
		throw new Error('internal error: unsupported returns computation method');
	}
	

	// Return the computed returns
	return returns;
}