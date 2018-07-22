### 0.0.7 - XX/YY/2018
- Redesigned the meanVarianceEfficientFrontier method to compute a desired number of efficient portfolios with their returns and volatilities.
- Added the meanVarianceCornerPortfolios to compute the corner portfolios with their returns and volatilities.
- Added the computation of the maximum Sharpe ratio portfolio, through the efficient frontier computation.
- Added a max function
- Added toTypedArray method for bit sets, and updated the critical line algorithm to use typed arrays for variables indexes
- Reorganized the bit set structure files
- Misc. tests improvements
- Added a FISTA-like optimization method for composite convex problems
- Added the computation of geometric center and geometric median
- Added a remark on the bias of the method randomKSubsetIterator_
- Added the computation of the random subspace mean variance optimization portfolio
- Added the computation of all the k-subsets of a n-set

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
