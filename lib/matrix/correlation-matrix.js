/**
* @author Roman Rubsamen <roman.rubsamen@gmail.com>
*/

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */



/**
* @function randomCorrelationMatrix
*
* @summary Returns a random correlation matrix.
*
* @description This function computes a random n by n correlation matrix, 
* using the algorithm described in the reference
* 
* @see <a href="https://link.springer.com/article/10.1023/A:1022384216930">Philip I. Davies and Nicholas J. Higham, Numerically Stable Generation of Correlation Matrices and Their Factors,BIT Numerical Mathematics volume 40, pages 640–651 (2000)</a>
*
* @param {number} n the row/column length of the matrix to construct, natural integer greater than or equal to 1.
* @param {object} opt optional parameters for the random correlation matrix generation algorithm.
* @param {number} opt.eps tolerance for the convergence of the algorithm, a strictly positive real number; defaults to 1e-14.
* @param {Array.<number>} opt.lambda the desired eigenvalues lambda_i, i=1..n of the generated correlation matrix, an array of n real numbers lambda_i,
* which must satisfy 0 <= lambda_i, i = 1..n and sum_i lambda_i = n; defaults to an array of random numbers belonging to [0,1] and summing to one.
* @param {number} opt.epsLambda tolerance for the condition that the provided eigenvalues must sum to one; defaults to 1e-8.
*
* @return {Matrix_} the computed matrix.
*
*/
self.randomCorrelationMatrix = function(n, opt) {
	// Generate the random correlation matrix
	return Matrix_.randomCorrelation(n, opt);
}

