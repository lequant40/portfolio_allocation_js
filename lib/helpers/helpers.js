/**
 * @file Functions related to helpers for portfolio allocation.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function Helpers
*
* @summary Construct a portfolio allocation helper object, unused.
*
* @description This function constructs a portfolio allocation helper object, unused.
* 
* @return {this} the constructed portfolio allocation helper, unused.
*
* @example
* Helpers();
*/
self.Helpers = function() {
	// Return
	return this;
}


/**
* @function covarianceMatrix
*
* @summary Returns the covariance matrix of a series of values.
*
* @description This function computes the covariance matrix of a series of values, provided as 
* a variable number of arrays of real numbers of the same length.
*
* @param {...Array.<number>} var_args, arrays of real numbers of the same length.
* @return {Array.<Array.<number>>} an array of array of real numbers representing the covariance matrix
* of the input series of values, organized by rows and then by columns.
*
* @example
* covarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
* // [[0.00036,  -0.00053], [-0.00053, 0.00107]]
*/
self.Helpers.covarianceMatrix = function(var_args) {
	// Construct the covariance matrix
	var arrs = arguments;	
	var cov = Matrix_.fillSymetric(arrs.length, function(i, j) { 
		return covariance_(arrs[i-1], arrs[j-1]); 
	});
	
	// Return it
	return cov.toDoubleArray();
}


/**
* @function sampleCovarianceMatrix
*
* @summary Returns the sample covariance matrix of a series of values.
*
* @description This function computes the sample covariance matrix of a series of values, provided as 
* a variable number of arrays of real numbers of the same length.
*
* @param {...Array.<number>} var_args, arrays of real numbers of the same length.
* @return {Array.<Array.<number>>} an array of array of real numbers representing the sample covariance matrix
* of the input series of values, organized by rows and then by columns.
*
* @example
* sampleCovarianceMatrix([0.05, 0.01, 0.01], [-0.05, 0.03, -0.01]);
* // [[0.00053, -0.0008], [-0.0008, 0.0016]]
*/
self.Helpers.sampleCovarianceMatrix = function(var_args) {
	// Construct the sample covariance matrix
	var arrs = arguments;
	var cov = Matrix_.fillSymetric(arrs.length, function(i, j) { 
		return sampleCovariance_(arrs[i-1], arrs[j-1]); 
	});
	
	// Return it
	return cov.toDoubleArray();
}