/**
 * @file Functions related to cluster risk parity portfolio.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

 
/* Start Wrapper private methods - Unit tests usage only */
/* End Wrapper private methods - Unit tests usage only */


/**
* @function clusterRiskParityWeights
*
* @summary Compute the weights of the cluster risk parity portfolio.
*
* @description This function returns the weights w_1,...,w_n associated to the fully invested and long-only 
* cluster risk parity portfolio of n assets, as computed by the cluster risk parity algorithm described in 
* the reference.
*
* This algorithm combines the use of a clustering method to isolate groups of assets and then allocate 
* both within and across these groups using equal risk contribution (ERC) weights.
*
* To be noted that the choice of the clustering method is not detailled in the reference, so that two possibilites
* are offered by this function:
* - (Default) Using automatically the Fast Threshold Clustering Algorithm (FTCA) from David Varadi, c.f. the ftca_ function.
* - Using a list of clusters provided by the user, typically constructed with the clustering algorithm of her choice
* 
* This portfolio is unique, provided the covariance matrix of the assets is definite positive.
* 
* @see <a href="https://cssanalytics.wordpress.com/2013/01/03/cluster-risk-parity/">Cluster Risk Parity</a>
* 
* @param {Matrix_|Array.<Array.<number>>} sigma the covariance matrix (sigma_ij),i,j=1..n of the n assets in the considered universe, square n by n Matrix or array of n arrays of n real numbers statisfying sigma[i-1][j-1] = sigma_ij.
* @param {object} opt the optional parameters for the algorithms used by the function.
* @param {number} opt.eps the tolerance parameter for the convergence of the ERC algorithms, a strictly positive real number; defaults to 1e-8.
* @param {number} opt.maxIter the maximum number of iterations of the ERC algorithms, a strictly positive natural integer; defaults to 10000.
* @param {number} opt.clusteringMode the method to use for the clusters computation, a string either equals to:
* - 'ftca': usage of the list of clusters automatically constructed by the FTCA algorithm from David Varadi
* - 'manual': usage of a list of clusters provided in input in the opt.clusters option
* ; defaults to 'ftca'.
* @param {number} opt.ftcaThreshold the correlation threshold to use in the FTCA algorithm in case opt.clusteringMode is equal to 'ftca', a real number beloning to [0,1]; defaults to 0.5.
* @param {Array.<Array.<number>>} opt.clusters the list of clusters to use in the algorithm in case opt.clusteringMode is equal to 'manual', an array of m arrays of strictly positive integers representing the indexes of the assets in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
* @return {Array.<number>} the weights corresponding to the cluster risk parity portfolio, array of n real numbers.
*
* @example
* clusterRiskParityWeights([[0.1,0], [0,0.2]], {clusteringMode: 'manual', clusters: [[1], [2]]});
*  
*/
self.clusterRiskParityWeights = function (sigma, opt) {
	// Decode options
	if (opt === undefined) {
		opt = {};
	}
	
	// Clustering options
	var clusteringMode = opt.clusteringMode || 'ftca';
	
	// Convert sigma to matrix format and convert it to acovariance matrix
	var sigma = new Matrix_(sigma).toCovarianceMatrix(sigma);
	var nbAssets = sigma.nbRows;

	
	// ------
	// Clustering logic
	var clusters = [];
	
	// In case clusters are provided, check that they form a partition of the set [1..nbAssets]
	if (clusteringMode === 'manual') {
		// Decode the manual clustering options
		clusters = opt.clusters || [];
		
		// Prepare a partition of the set of integers [1..nbAssets]
		var partition = new Array(nbAssets);
		for (var i = 0; i < partition.length; ++i) {
		    partition[i] = 0;
		}
		
		// Count the number of times each asset index appears in the list of clusters
	 	for (var j = 0; j < clusters.length; ++j) {
	        // Detect empty clusters
			if (clusters[j].length === 0) {
				throw new Error('empty cluster at index: ' + j);
			}
			
			// Detect out of bounds asset indexes
			// Count the number of times each asset index appears in the list of clusters
			for (var k = 0; k < clusters[j].length; ++k) {
	            var assetIdx = clusters[j][k];
	            if (assetIdx < 1 || assetIdx > nbAssets) {
	                throw new Error('asset index out of bounds: ' + assetIdx);
	            }
	            else {
	                partition[assetIdx-1]++;
	            }
	        }
	    }
		    
		// Check that each integer in the set [1..nbAssets] is appears once and only once 
		// in the list of clusters
		for (var i = 0; i < partition.length; ++i) {
		    if (partition[i] !== 1) {
		        if (partition[i] === 0) {
		            throw new Error('missing asset index: ' + (i+1));
		        }
		        else if (partition[i] > 1) {
		            throw new Error('duplicate asset index: ' + (i+1));
		        }
		        else {
		            throw new Error('unknown error');
		        }
		    }
		}

		// All checks passed
	}
	// Otherwise, compute the clusters using the Fast Threshold Clustering Algorithm from David Varadi
	else if (clusteringMode === 'ftca') {
		// Extract the correlation matrix from the covariance matrix, as a double array
		var corrMat = sigma.getCorrelationMatrix().toRowArray();
		
		// Compute the clusters using the FTCA algorithm
		clusters = ftca_(corrMat, opt.ftcaThreshold);
	}
	else {
		//
		throw new Error('unsupported clustering method');
	}
	

	// ------
	// The cluster risk parity portfolio is constructed in three steps:
	// 1 - Clusters creation (done above)
	// 2 - Assets allocation within each cluster using ERC weights
	// 3 - Portfolio allocation across all clusters using ERC weights and considering each cluster as a (synthetic) asset
	
	// 1 - N/A
	var nbClusters = clusters.length;
	
	// Instantiate the matrix mapping the initial assets to the clusters space (e.g., a change of base matrix),
	// to be populated with proper weights in step 2 below.
	var assetsToClustersWeights = Matrix_.zeros(nbClusters, nbAssets); 
	
	// 2 - For each cluster:
	for (var i = 0; i < nbClusters; ++i) {
		// Extract the covariance matrix associated to the assets it contains 
		var assetsIndexes = clusters[i].slice().sort(function(a, b) { return a - b; });
		
		var clusterSigma = sigma.submatrix(assetsIndexes, assetsIndexes);
		
		// Compute ERC weights for these assets
		var assetsWeights = self.equalRiskContributionWeights(clusterSigma, opt);
		
		// - Populate the change of base matrix for the current cluster with the computed ERC weights
		for (var j = 0; j < assetsWeights.length; ++j) {
			assetsToClustersWeights.setValueAt(i+1, assetsIndexes[j], assetsWeights[j]);
		}
	}
	
	// 2(tmp) - Compute the transpose of the matrix mapping the initial assets to the clusters space
	var assetsToClustersWeightsT = assetsToClustersWeights.transpose()
	
	// 3a - Compute the the covariance matrix associated to the weighted clusters space, using formula Var(A*X) = A * Var(X) * A'.
	var clustersSigma = Matrix_.xy(assetsToClustersWeights, Matrix_.xy(sigma, assetsToClustersWeightsT));
	
	// 3b - Compute ERC weights in the clusters space
	var clustersWeights = self.equalRiskContributionWeights(clustersSigma, opt);
	
	// 3c - Compute original assets weights, using formula A' * Y
	var weights = Matrix_.xy(assetsToClustersWeightsT, new Matrix_(clustersWeights));
	
	// Return them (already normalized)
	return weights.toArray();
}

