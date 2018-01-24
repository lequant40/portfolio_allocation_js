// ------------------------------------------------------------
QUnit.module('Optimisation internal module', {
  before: function() {
    // 
  }
});


QUnit.test('Linear programming solver - Primal Dual Hybrid Gradient', function(assert) {    
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
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(null, null, Ai, bi, c);

	var expectedX = new PortfolioAllocation.Matrix([2, 0, 1]);
	var expectedMinVal = -13;
	
	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Feasible #1 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Feasible #1 - 2/2');
  }
  
  // Problem, p. 19
  {
	var Ai = new PortfolioAllocation.Matrix([[-1, 3], [1, 1], [2, -1]]);
	var bi = new PortfolioAllocation.Matrix([12, 8, 10]);
	var c = new PortfolioAllocation.Matrix([-3, -2]);
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(null, null, Ai, bi, c);
	
	var expectedX = new PortfolioAllocation.Matrix([6, 2]);
	var expectedMinVal = -22;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Feasible #2 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Feasible #2 - 2/2');
  }
  
  // Reference: Dantzig George, Linear programming and extensions
  // Illustrative example 1, p. 108
  {
	var Ae = new PortfolioAllocation.Matrix([[5, -4, 13, -2, 1], [1, -1, 5, -1, 1]]);
	var be = new PortfolioAllocation.Matrix([20, 8]);
	var c = new PortfolioAllocation.Matrix([1, 6, -7, 1, 5]);
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(Ae, be, null, null, c);
	
	var expectedX = new PortfolioAllocation.Matrix([0, 4/7, 12/7, 0, 0]);
	var expectedMinVal = -60/7;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Feasible #3 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Feasible #3 - 2/2');
  }
  
  // Illustrative example 1, p. 110
  {
	var Ae = new PortfolioAllocation.Matrix([[-1, -1, -1, -1, -1, -1, -1, -1,-1], 
											[-0.1, -0.1, -0.4, -0.6, -0.3, -0.3, -0.3, -0.5, -0.2], 
											[-0.1, -0.3, -0.5, -0.3, -0.3, -0.4, -0.2, -0.4, -0.3], 
											[-0.8, -0.6, -0.1, -0.1, -0.4, -0.3, -0.5, -0.1, -0.5]]);
	var be = new PortfolioAllocation.Matrix([-1, -0.3, -0.3, -0.4]);
	var c = new PortfolioAllocation.Matrix([4.1, 4.3, 5.8, 6.0, 7.6, 7.5, 7.3, 6.9, 7.3]);
	 
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(Ae, be, null, null, c);
	
	var expectedX = new PortfolioAllocation.Matrix([0, 0.6, 0, 0.4, 0, 0, 0, 0, 0]);
	var expectedMinVal = 4.98;
	
	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Feasible #4 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Feasible #4 - 2/2');
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
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(null, null, Ai, bi, c);
	
	var expectedX = new PortfolioAllocation.Matrix([20, 24]);
	var expectedMinVal = -1080000;

	assert.equal(PortfolioAllocation.Matrix.areEqual(sol[0], expectedX, 1e-04), true, 'Feasible #5 - 1/2');
	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'Feasible #5 - 2/2');
  }
  
  // Reference: http://www.numerical.rl.ac.uk/cute/netlib.html, AFIRO problem
  {
	var Ae = new PortfolioAllocation.Matrix([[-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
											  [-1.06, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ],
											  [0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, -1.06, -1.06, -0.96, -0.86, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.43, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.43, -0.43, -0.39, -0.37, 0, 0, 0, 0, 0, 0, 1, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 1, 0, 1]]);
	var be = new PortfolioAllocation.Matrix([0, 0, 0, 0, 0, 0, 0, 44]);
	var Ai = new PortfolioAllocation.Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.4, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 2.364, 2.386, 2.408, 2.429, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2.191, 2.219, 2.249, 2.279, 0, 0, 0, 0],
											  [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.109, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0.109, 0.108, 0.108, 0.107, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0.301, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0.301, 0.313, 0.313, 0.326, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
											  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
											  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]]); 
	var bi = new PortfolioAllocation.Matrix([80, 0, 80, 0, 0, 0, 500, 0, 500, 0, 0, 0, 0, 0, 0, 0, 0, 310, 300]);
	var c = new PortfolioAllocation.Matrix([0, -0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.32, 0, 0, 0, -0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.48, 0, 0, 10]);
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(Ae, be, Ai, bi, c);
	
	var expectedMinVal = -4.6475314286e02; // Direclty from Netlib README file

	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'AFIRO');
  }
  
  // Reference: http://www.numerical.rl.ac.uk/cute/netlib.html, SC50B problem
  {	
	var Ae = new PortfolioAllocation.Matrix([[0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,1.1000000000000001,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1000000000000001,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1000000000000001,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1000000000000001,0,0,0,0,0,0,0,0,0,0,-1]]);
	var be = new PortfolioAllocation.Matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
	var Ai = new PortfolioAllocation.Matrix([[3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], // 2 null rows removed here
											[0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,-1,0,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0.40000000000000002,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0.59999999999999998,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40000000000000002,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.59999999999999998,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40000000000000002,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.59999999999999998,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,3,3,3,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40000000000000002,0,0,0,0,-1,0,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.59999999999999998,0,0,0,0,0,-1,0,0,0,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.69999999999999996,0.29999999999999999,0.29999999999999999,0],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0.40000000000000002],
											[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0.59999999999999998]]);
	var bi = new PortfolioAllocation.Matrix([300, 0, 0, 0, 300, 0, 0, 0, 0, 0, 300, 0, 0, 0, 0, 0, 300, 0, 0, 0, 0, 0, 300, 0, 0, 0, 0, 0]);  // 2 null rows removed at indexes 2 and 3
	var c = new PortfolioAllocation.Matrix([0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
	
	var sol = PortfolioAllocation.lpsolvePrimalDualHybridGradient_(Ae, be, Ai, bi, c);
	
	var expectedMinVal = -7.0000000000e01; // Direclty from Netlib README file

	assert.equal(Math.abs(sol[1] - expectedMinVal) <= Math.abs(expectedMinVal) * 1e-06, true, 'SC50B');
  }
});