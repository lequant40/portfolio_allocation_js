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
