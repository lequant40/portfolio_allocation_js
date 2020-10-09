### 0.0.12 - xx/yy/2020

### 0.0.11 - 09/10/2020
- Cosmetic updates related correlation matrices fixing, and added a unit test to cover the infinite precision of the nearest correlation matrix computation
- Improved the numerical stability of the Constrained Risk Budgeting algorithm: added auto-expanding bisection interval, and early stopping criterion
- Added methods to swap two rows or two columns of a matrix
- Added an LU decomposition method for square matrices, with full pivoting
- Added a method to solve lower triangular linear systems
- Reworked the internals of the critical line algorithm, to solve the KKT linear systems with a LU decomposition with full pivoting instead, to guard against numerical issues
- Reworked the internals of the critical line algorithm to allow for equal returns; the only case in which the CLA algorithm will not work is if the covariance matrix is not invertible, if there are more than 1 asset IN in the maximum return portfolio, and if the associated KKT matrix is also non invertible
- Fixed a bug in the GSMO optimization method, which could seemingly lead to cycling (or at least required an unreasonably huge number of iterations before convergence)
- Added the critical line method as additional optimization method for the GMV portfolio
- Added the critical line method as additional optimization method for the MDP portfolio
- Improved the numerical precision of the corner portfolios, of the mean-variance efficient portfolios, and of the maximum Sharpe ratio portfolios
- Added a RSO-MSR optimization method
- Redesigned the asset allocation non-regression to factorize the tests on portfolio weights using an equality function
- Redesigned the internals of the mean-variance optimization method, to automatically selects the critical line algorithm, or the GSMO algorithm, as optimization routine; this impacts all portfolio computations using MVO
- Allowed for an infinite number of cycles in the risk budgeting computation
- Added lower bounds and upper bounds constraints in the equal risk bounding method
- Updated the covariance matrix computation method to be able to use Bessel's correction with Ledoit-Wolf shrinkage (= non-backward compatible change)
- Updated the mean vector computation method to be consistent with the covariance computation method w.r.t. the regularizationMethod parameter (= non-backward compatible change)
- Updated the returns computation method to generate an error in case of division by zero
- Added a method to test that a matrix is a covariance matrix
- Fixed a bug in Cholesky decomposition, which wrongly allowed for semi-definite positive matrices to be decomposed
- Added a method to compute a Cholesky decomposition for a positive semi-definite matrix
- Misc. checks added in the covariance matrix creation methods
- Improved the method to check a correlation matrix by allowing details on what property the matrix is missing
- Added a method to compute a covariance matrix from a supposed covariance matrix
- Updated the Threshold Accepting method to enable optimization in computations of successive values of the function f, and used this optimization with the postOptimizationWeights and minimumTrackingErrorWeights methods
- Removed partial investment constraints on MVO, Max Sharpe portfolios
- Improved the performances of the minimum tracking error portfolio in case of cardinality constraints 
- Removed the post optimization method "rationalizing" (= non backward compatible change), as this one is useless in practice
- Misc. performances and randomization improvements in the post optimization method internals
- Automatic computation of the portfolio target cash position, as well as management of levered portfolio, in the post optimization weights method
- Fixed a numerical bug in correlation matrix extraction from covariance matrix, which was not strictly symmetric and unit diagonal
- Polished the output of the nearest correlation matrix computation to have a strictly symmetric and unit diagonal matrix in output whatever the requested precision
- Misc. optimizations in the critical-line algorithm when assets returns are identical
- Misc. optimizations in the computation of efficient portfolios in case the GSMO algorithm is used
- Added the possibility to provide an initial point to the GSMO algorithm for warm start
- Fixed a bug in the projection on the efficient frontier in case the GSMO algorithm is used, which lead to the same efficient portfolios being computed on the efficient frontier, thus biaising the projection
- Removed partial investment constraint support from all MV-related algorithms (= non backward compatible changes), and introduced a proper portfolio exposure management in the RSO-MV related algorithms
- Removed partial investment constraint support from minimum tracking error algorithm, in case of no cardinality constraints (= non backward compatible changes)
- Fixed a bug in all MV-related algorithms where in case of tight min/max weight constraints, the algorithms would loop infinitely

### 0.0.10 - 14/09/2020
- Added the Jacobi method in order to compute eigenvalues and eigenvectors of real symmetric matrices
- Updated the method to generate random correlation matrices to output truly symmetric matrices, and not numerically symmetric matrices only
- Added a method to repair the (semi) definite-positiveness of correlation matrices, using spectral decomposition
- Change of method name for Ledoit-Wolf covariance matrix shrinkage to "linear-shrinkage" (= non-backward compatible change)
- Change of method name for DeMiguel mean vector shrinkage to "linear-shrinkage" (= non-backward compatible change)
- Added a method to repair the (semi) definite-positiveness of correlation matrices, using shrinkage towards a target matrix
- Added the possibility to provide standard deviations for random mean vectors and random variances vectors creation
- Added a method to randomly perturb a correlation matrix
- Added a method to compute the nearest correlation matrix in Frobenius norm, including a lower bound constraint on the smallest eigenvalue(s)
- Added a method to repair the (semi) definite-positiveness of correlation matrices, using the computation of the nearest correlation matrix
- Updated grunt version
- Added a method to compute a covariance matrix from a correlation matrix and a vector of either variances or standard deviations
- Updated the internals of the Minimum Tracking Error portfolio: now Threshold Accepting is the default method in case of cardinality constraints because of speed/precision/stability improvements in the neighbourhood generation function
- Updated the internals of the Minimum Tracking Error portfolio: the generation of the initial point is now customizable in terms of the number of iterations
- Updated the interface of the Minimum Tracking Error portfolio (= non-backward compatible change), the optimization method is now either "exact" or "heuristic" instead of "combinatorial" and "thresholdAccepting"
- Updated the Minimum Tracking Error portfolio, adding the possibility to allow for partial investment
- Fix for issue https://github.com/lequant40/portfolio_allocation_js/issues/7 (new way to instantiate Matrix objects and detect arrays)
- Added a method to symmetrize a square matrix
- Added a method to set 1's on the diagonal of a square matrix
- Added a method to test if a matrix is unit diagonal
- Added the Cholesky decomposition of a positive definite matrix
- Added a method to determine if a matrix is a correlation matrix
- Redesigned the meanVarianceEfficientFrontierPortfolios method to output portfolios uniformly distributed w.r.t. either return, 
volatility of risk tolerance parameter; default is now return to have the whole efficient frontier properly covered (= non-backward compatible change)
- Added a method to test for the emptyness of the restricted full unit simplex of R^n
- Added a method to compute the euclidean projection of a point of R^n in the restricted full unit simplex of R^n
- Added the characteristic function of the restricted full unit simplex of R^n
- Fixed a bug in meanVarianceEfficientFrontierPortfolios, in which partial investment constraint was not taken into account
- Fixed a bug in meanVarianceEfficientFrontierNearestWeights, in which partial investment constraint was not taken into account
- Misc. internal changes in the maximumSharpeRatioWeights method, plus fixed a bug in which partial investment constraint was not taken into account
- Updated the error message in case the critical line algorithm is failing due to equal returns, to better explain the issue behind https://github.com/lequant40/portfolio_allocation_js/issues/8 
- Removed the possibility to compute corner portfolios from the outside of the library, as preparatory step to use another algorithm for the mean variance efficient frontier computation
- Added a function to minimize a linear function over the restricted unit simplex of R^n
- Improved the precision of the method to project on the restricted unit simplex of R^n, which was not precise for big numbers, typically encountered in FISTA algorithm
- As a side effect, the minimumTrackingErrorWeights method now computes better portfolios
- Reworked some internals of the computeCornerPortfolios_ method related to the computation of the E-maximizing portfolio
- Reword the internals of the mean-variance optimization methods, so that other algorithms that the critical line algorithm can be used to compute MV-efficient portfolios
- Reworked the internals of the computeEfficientFrontierPortfolios, to remove the limit cases around the number of portfolios to be computed
- Reworked the internals of the minimumTrackingErrorWeights method to improve the way partial investment constraint is managed (better numerical precision and better speed)
- Reworked the mean-variance optimization to use the GSMO algorithm by default instead of the critical line algorithm, to ensure definiteness in all cases 
- Reworked the internals of the globalMinimumVarianceWeights method to handle partial investment constraint, optional returns, as well as to re-use the new internals of the mean-variance optimization methods
- Fixed misc. bugs in goldenSectionSearch optimization method, allowing the case of upper bound = lower bound, allowing the case where the minimum is found on one of the bounds, and allowing a constant function on the right of the interval
- Reworked the internals of the maximumSharpeRatioWeights method to re-use the new internals of the mean-variance optimization methods
- Reworked the internals of the mostDiversifiedWeights method to re-use the new internals of the mean-variance optimization methods, which allows for min/max weights and partial investment constraints
- Added a RSO-GMV optimization method
- Updated the randomSubspaceMeanVarianceOptimizatioWeights interface (non backward-compatible change), renamed the options subsetsOpt

### 0.0.9 - 03/08/2020
- Removed comments from the generated Google Sheet script using a new grunt plugin
- Added a method to extract the columns of a matrix
- Added a method to generate random normal numbers, using the inverse method
- Added a method to generate random matrices made of random normal numbers
- Added a method to generate random orthogonal matrices
- Added a method to generate random correlation matrices
- Reworked the internals of the risk budgeting portfolio algorithm, and implemented 3 new coordinates sampler algorithms
- Fixed bug https://github.com/lequant40/portfolio_allocation_js/issues/5 related to MVO optimization with maxTargetVolatility optimization method
- Added min/max weights constraints in the ERC and RB portfolio allocation methods
- Added a method to generate random normal numbers, with positive support (i.e., truncated to R^+)
- Added a method to generate randoms variances
- Added a method to randomly perturb a variances vector
- Added a method to test if a matrix is symmetric
- Fixed bug https://github.com/lequant40/portfolio_allocation_js/issues/6 related to MVO optimization (corner portfolios computation)
- Added a method to compute p-quantiles of a series of values
- Updated the method to generate random portfolios: providing no cardinality constraints now defaults to strictly taking into account the min/max weights constraints
- Added an heuristic stochastic optimization method: Threshold Accepting
- Updated the minimum tracking error portfolio computation: added the Threshold Accepting optimization method to compute an approximate solution in case of cardinality constraints
- Merged the computation of the covariance matrix and of the sample covariance matrix, breaking the public interface (= non-backward compatible interface change)
- Added 4 shrinkage methods to compute covariance matrices (all based on Ledoit-Wolf's linear shrinkage)
- Misc. refactoring of the mean-variance optimization methods, breaking the public interface (= non-backward compatible interface change)
- Removed the possibility to add soft inequality constraints to the minimum tracking error portfolio (= non-backward compatible change)
- Added a new constraint in the mean-variance optimization algorithm: risk tolerance
- Updated the internals of the mean-variance optimization, so that all corner portfolios are now retained (numerically close portfolios were merged, which is incompatible with the addition of a risk tolerance constraint due to kinks, c.f. the unit tests)
- Redesigned the meanVarianceEfficientFrontierPortfolios method to output portfolios uniformly distributed w.r.t. the risk tolerance parameter
- Renamed the roundedWeights optimization method into postOptimizationWeights to better reflect what this method is doing, plus added a method to convert numerical portfolio weights into integer number of shares (roundlotting)
- Added a method to compute the mean vector of series of values
- Added a method to randomly perturb a mean vector
- Added a method to generate random mean vectors
- Renammed equalRiskBudgetWeights to inverseVolatilityWeights for better consistency (= non-backward compatible change)
- Added a method to compute a shrinked mean vector of a series of values
- Added a method to compute the euclidean projection on a line segment in R^n
- Added a method to compute the nearest portfolio located on the mean-variance efficient frontier w.r.t. a given input portfolio
- Added a method to compute arithmetic and logarithmic returns

### 0.0.8 - 16/04/2020
- Introduced a function simplexEmptinessCheck_ to avoid copy/pasting feasibility checks on the restricted simplex everywhere
- Added min/max exposure constraints, and min/max weights constraints in the random portfolio generation method
- Proper management of infeasible subsets in the randomSubspaceMeanVarianceOptimizationWeights method
- Added a method to generate random permutations of integers
- Used an accurate algorithm to compute the l2 norm of vectors in the methods hypersphereRandomSampler_ and simplexDirectionRandomSampler_
- Added a method to generate all the permutations of a n-set (Heap's algorithm)
- Added a method to compute the permutation entropy of a time series
- Added an hyper rectangle regular grid sampler
- Added a unidimensional root finding method: the bisection method
- Added the computation of the best constantly rebalanced portfolio
- Improved the convergence of the geometric median computation
- Added a unidimensional optimization method method: the golden section search method
- Added a generating set search optimization method
- Added the computation of the minimum tracking error portfolio

### 0.0.7 - 25/05/2019
- Redesigned the meanVarianceEfficientFrontier method to compute a desired number of efficient portfolios with their returns and volatilities
- Added the meanVarianceCornerPortfolios method to compute the corner portfolios with their returns and volatilities
- Added the computation of the maximum Sharpe ratio portfolio, through the efficient frontier computation
- Added a max function
- Added toTypedArray method for bit sets, and updated the critical line algorithm to use typed arrays for variables indexes
- Reorganized the bit set structure files
- Misc. tests improvements
- Added a FISTA-like optimization method for composite convex problems
- Added the computation of geometric center and geometric median of m points in R^n, using the FISTA method
- Added the computation of all the k-subsets of a n-set
- Added a maximum volatility target constraint for the mean-variance optimization algorithm
- Added an helper method to compute portfolios: the random subspace method
- Applied the random subspace method to the mean variance optimization portfolio
- Replaced the randomKSubsetIterator_ internal method RANKSB from A. Nijenhuis and ‎H.S. Wilf by the method D from J.S. Vitter, which is provably faster and uniform
- Added a method to approximately compute the inverse of the standard normal cumulative distribution function
- Added a method to generate uniformly distributed vectors on the R^n unit hypersphere
- Added a method to generate uniformly distributed vectors on the intersection of the R^n unit hypersphere and of the R^n hyperplane defined by the equation <(1,1,...,1)/x> = 0,
in preparation for random optimization algorithms
- Fix for issue https://github.com/lequant40/portfolio_allocation_js/issues/3 (proper management of semi-definite positive covariance matrices in the ERC/RB algorithm)
- Added the computation of random compositions of an integer, using the algorithm RANCOM from A. Nijenhuis and ‎H.S. Wilf
- Added the support of bounds contraints to the sampling of points uniformly at random on the simplex of R^n
- Misc. rework of the simplex rational grid search function and associated functions
- Fix for issue https://github.com/lequant40/portfolio_allocation_js/issues/4 (updated README, to be consistent between version displayed/published on third-party site) 
- Reworked the generic numerical optimization portfolio algorithm, now called numericalOptimizationWeights, and added bound contraints
- Renamed the MCA and MVA algorithms to minimumCorrelationWeights and proportionalMinimumVarianceWeights to be consistent with the other algorithms
- Updated the grunt dependencies versions

### 0.0.6 - 10/06/2018
- Updated the Maximin portfolio method output signature to allow outputting only the portfolio weights
- Misc. matrix operations refactor
- Added the computation of the Markowitz efficient frontier, through the critical line algorithm
- Added the computation of mean-variance efficient portfolios with a given return or volatility target, through the critical line algorithm

### 0.0.5 - 02/04/2018
- Misc. Matrix/Vector operations additions/refactor (in particular, toDoubleArray method replaced by toRowArray)
- Added two iterative linear system solvers (least square sense): Kaczmarz extended method and randomized Kaczmarz extended method
- Added a feasible/bounded linear program solver: primal dual hybrid gradient algorithm
- New portfolio allocation methods: Random Portfolio, Minimax Portfolio
- Updated portfolio allocation helper method to round portfolio weights 
- Added a quadratic program solver: generalized sequential minimization optimization
- Added a continuous quadratic knapsack problem solver using a breakpoint searching algorithm
- Updated GMV and MDP portfolios to use the quadratic program solver
- Updated GMV to manage assets weights constraints
- Added a version of the select algorithm from Floyd and Rivest
- Added a median finding algorithm, using the select algorithm from Floyd and Rivest
- Optimized the ERC and ERB portfolio allocation algorithms, going from ~5 minutes to find the ERB portfolio for 20 assets to ~30 seconds

### 0.0.4 - 01/12/2017
- Misc. Matrix/Vector operations additions/refactor
- New portfolio allocation methods: Equal Risk Bounding, Global Minimum Variance, Grid Search
- New portfolio allocation helper method: "optimal" floating-point portfolio weights transformation to rational portfolio weights

### 0.0.3 - 10/09/2017
- Misc. Matrix/Vector operations additions/refactor
- New portfolio allocation method: Cluster Risk Parity, with a included clustering method: Fast Threshold Clustering Algorithm

### 0.0.2 - 24/06/2017
- New portfolio allocation methods: most diversified portfolio, minimum correlation portfolio, proportional minimum variance portfolio
- Included a covariance computation function using a corrected two pass formula
- New Matrix/Vector operations, plus code refactoring for future split

### 0.0.1 - 09/06/2017
- Initial commit
- Four portfolio allocation methods supported: equal weights, equal risk budget, equal risk contributions, risk budgeting
- Several Matrix/Vector basic operations implemented to support above algorithms
