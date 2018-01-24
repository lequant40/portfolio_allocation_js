// ------------------------------------------------------------
QUnit.module('Optimisation internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Linear programming solver - Chambolle-Pock', function(assert) {    
  // Reference:  Vanderbei, Robert J, Linear Programming Foundations and Extensions 4th edition
  // Problem, p. 11
/*
  // Infeasible problem, p. 7
  {
	var A = new PortfolioAllocation.Matrix([[1, 1], [-2, -2]]);
	var b = new PortfolioAllocation.Matrix([2, -9]);
	var c = new PortfolioAllocation.Matrix([-5, -4]);

	assert.throws(function() { PortfolioAllocation.lpsolveAffineScaling_(A, b, c, {problemForm: 'standard'}) },
		                         new Error('infeasible problem detected'),
		                         "");
	
  }
  
  // Infeasible problem, p. 62
  {
	var A = new PortfolioAllocation.Matrix([[1, -1], [-1, 1]]);
	var b = new PortfolioAllocation.Matrix([1, -2]);
	var c = new PortfolioAllocation.Matrix([-2, 1]);
	
	assert.throws(function() { PortfolioAllocation.lpsolveAffineScaling_(A, b, c, {problemForm: 'standard'}) },
		                         new Error('infeasible problem detected'),
		                         "");
	
  }
  
  // Unbounded problem, p. 7
  {
	var A = new PortfolioAllocation.Matrix([[-2, 1], [-1, -2]]);
	var b = new PortfolioAllocation.Matrix([-1, -2]);
	var c = new PortfolioAllocation.Matrix([-1, 4]);
	
	assert.throws(function() { PortfolioAllocation.lpsolveAffineScaling_(A, b, c, {problemForm: 'standard'}) },
	                         new Error('unbounded problem detected'),
	                         "");
  }
 */
 
  // Reference:  Vanderbei, Robert J, Linear Programming Foundations and Extensions 4th edition
  // Problem, p. 11
 {
	var Ai = new PortfolioAllocation.Matrix([[2, 3, 1], [4, 1, 2], [3, 4, 2]]);
	var bi = new PortfolioAllocation.Matrix([5, 11, 8]);
	var c = new PortfolioAllocation.Matrix([-5, -4, -3]);
	
	var sol = PortfolioAllocation.lpsolveChambollePock_(null, null, Ai, bi, c);

	var expectedX = new PortfolioAllocation.Matrix([2, 0, 1]);
	var expectedMinVal = -13;
	
	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Simple feasible #1 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Simple feasible #1 - 2/2');
  }
  
  // Problem, p. 19
  {
	var Ai = new PortfolioAllocation.Matrix([[-1, 3], [1, 1], [2, -1]]);
	var bi = new PortfolioAllocation.Matrix([12, 8, 10]);
	var c = new PortfolioAllocation.Matrix([-3, -2]);
	
	var sol = PortfolioAllocation.lpsolveChambollePock_(null, null, Ai, bi, c);
	
	var expectedX = new PortfolioAllocation.Matrix([6, 2]);
	var expectedMinVal = -22;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Simple feasible #2 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Simple feasible #2 - 2/2');
  }
  
  // Reference: Dantzig George, Linear programming and extensions
  // Illustrative example 1, p. 108
  {
	var Ae = new PortfolioAllocation.Matrix([[5, -4, 13, -2, 1], [1, -1, 5, -1, 1]]);
	var be = new PortfolioAllocation.Matrix([20, 8]);
	var c = new PortfolioAllocation.Matrix([1, 6, -7, 1, 5]);
	
	var sol = PortfolioAllocation.lpsolveChambollePock_(Ae, be, null, null, c);
	
	var expectedX = new PortfolioAllocation.Matrix([0, 4/7, 12/7, 0, 0]);
	var expectedMinVal = -60/7;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Simple feasible #3 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Simple feasible #3 - 2/2');
  }
  
  // Illustrative example 1, p. 110
  {
	var Ae = new PortfolioAllocation.Matrix([[-1, -1, -1, -1, -1, -1, -1, -1,-1], 
											[-0.1, -0.1, -0.4, -0.6, -0.3, -0.3, -0.3, -0.5, -0.2], 
											[-0.1, -0.3, -0.5, -0.3, -0.3, -0.4, -0.2, -0.4, -0.3], 
											[-0.8, -0.6, -0.1, -0.1, -0.4, -0.3, -0.5, -0.1, -0.5]]);
	var be = new PortfolioAllocation.Matrix([-1, -0.3, -0.3, -0.4]);
	var c = new PortfolioAllocation.Matrix([4.1, 4.3, 5.8, 6.0, 7.6, 7.5, 7.3, 6.9, 7.3]);
	 
	var sol = PortfolioAllocation.lpsolveChambollePock_(Ae, be, null, null, c);
	
	var expectedX = new PortfolioAllocation.Matrix([0, 0.6, 0, 0.4, 0, 0, 0, 0, 0]);
	var expectedMinVal = 4.98;
	
	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Simple feasible #4 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Simple feasible #4 - 2/2');
  }
 
  // Reference: https://math.stackexchange.com/questions/59429/berlin-airlift-linear-optimization-problem
  // Note: with 9000 and 5000 as 3rd constrains instead of 9 and 5, and with the associated limit 300 000 instead of 300, 
  // the A matrix and the b vector are badly scaled (ex. for b - min value: 44, max value 300 000), which slows down the algorithm.
  {
	var Ai = new PortfolioAllocation.Matrix([[1, 1], 
											[16, 8], 
											[9, 5]]); 
	var bi = new PortfolioAllocation.Matrix([44, 512, 300]);
	var c = new PortfolioAllocation.Matrix([-30000, -20000]);
	
	var sol = PortfolioAllocation.lpsolveChambollePock_(null, null, Ai, bi, c);
	
	var expectedX = new PortfolioAllocation.Matrix([20, 24]);
	var expectedMinVal = -1080000;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Simple feasible #5 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Simple feasible #5 - 2/2');
  }
});