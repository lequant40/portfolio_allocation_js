/**
 * @file Misc. statistical functions.
 * @author Roman Rubsamen <roman.rubsamen@gmail.com>
 */

/* Start Wrapper private methods - Unit tests usage only */
self.max_ = max_;
self.median_ = median_;
self.select_ = select_;
self.hypot_ = hypot_;
self.rank_ = rank_;
self.ftca_ = ftca_;
self.normcdf_ = function(x) { return normcdf_(x); }
self.norminv_ = norminv_;
self.hypersphereRandomSampler_ = hypersphereRandomSampler_;
/* End Wrapper private methods - Unit tests usage only */


/**
* @function max_
*
* @summary Compute the maximum of a serie of values.
*
* @description This function returns the maximum of a serie of values [x_1,...,x_n],
* as well as its index.
*
* In case there are several identical maximum values, the one corresponding to the
* lowest indice in the array x is returned.
*
* @param {Array.<number>} x an array of real numbers.
* @param {function} compareFunction an optional sort function that defines the sort order, using the standard prototype for JavaScript sort functions (c.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/sort).
* @return {Array.<number>} an array arr of two elements:
* arr[0], the maximum of the values of the array x, a real number
* arr[1], the index of the maximum of the values of the array x, a positive integer
*
* @example
* max_([2,4,4,1]);
* // [4, 1]
*/
function max_(x, compareFunction) {
	// Initialisations.
	var defaultCompareFct = function (a, b) {
		return a - b;
	};
	var compareFunction = compareFunction || defaultCompareFct;
	
	var n = x.length;
	
	// Core loop
	var maxValue = x[0];
	var maxValueIdx = 0;
	for (var i = 1; i < n; ++i) {
		//if (x[i] > maxValue) {
		if (compareFunction(x[i], maxValue) > 0) {
			maxValue = x[i];
			maxValueIdx = i;
		}
	}
	
	// Return the computed maximum value and its index
	return [maxValue, maxValueIdx];
}


/**
* @function median_
*
* @summary Compute the median of a serie of values.
*
* @description This function returns the median of a serie of values [x_1,...,x_n], 
* which is defined as:
* - When n is odd, the (n+1)/2-th smallest element of the p values x_1,...,x_n
* - When n is even, the n/2-th smallest element of the p values x_1,...,x_n
*
* The algorithm used internally is based on the O(n) SELECT algorithm of the reference.
*
* @see <a href="https://www.sciencedirect.com/science/article/pii/S0304397505004081">Krzysztof C. Kiwiel, On Floyd and Rivest's SELECT algorithm, Theoretical Computer Science, Volume 347, Issues 1–2, 2005, Pages 214-238</a>
* 
* @param {Array.<number>} x an array of real numbers.
* @param {function} compareFunction an optional sort function that defines the sort order, using the standard prototype for JavaScript sort functions (c.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/sort).
* @return {number} the median of the values of the array x.
*
* @example
* median_([2,4,1]);
* // 2
*
* median_([2,4,1,3]);
* // 2.5
*/
function median_(x, compareFunction) {
	// Initialisations.
	var n = x.length;
	var xx = x.slice(); // to avoid altering the array x
	
	// Compute the smallest |-n/2-| element of the array, which corresponds to the median
	return select_(xx, Math.ceil(n/2), compareFunction);
}


/**
* @function select_
*
* @summary Compute the smallest k element of a serie of values.
*
* @description This function permutes a serie of values x = [x_1,...,x_n] so that:
* - The smallest k elements of x are x[i], i=0..k-1 (in an arbitrary order)
* - The k-th smallest element of x is x[k-1]
* - The n-k-th largest elements of x are x[i], i=k..n-1 (in an arbitrary order)
*
* The algorithm used internally is the O(n) algorithm of the reference.
*
* This code is a port to JavaScript by Roman Rubsamen of the Fortran 77 code
* written by K.C. Kiwiel, version of the 8 March 2006, kiwiel@ibspan.waw.pl.,
* except for the indices computation part which is new.
*
* The Fortran 77 version was a Fortran code for the Algol 68 procedure from
* the second reference, including some modifications suggested in the third 
* reference.
*
* @see <a href="https://www.sciencedirect.com/science/article/pii/S0304397505004081">Krzysztof C. Kiwiel, On Floyd and Rivest's SELECT algorithm, Theoretical Computer Science, Volume 347, Issues 1–2, 2005, Pages 214-238</a>
* @see <a href="https://dl.acm.org/citation.cfm?doid=360680.360694">R.W. Floyd and R.L. Rivest: "Algorithm 489: The Algorithm SELECT---for Finding the $i$th Smallest of $n$ Elements", Comm. ACM 18, 3 (1975) 173</a>
* @see <a href="https://dl.acm.org/citation.cfm?id=355704">T. Brown: "Remark on Algorithm 489", ACM Trans. Math. Software 3, 2 (1976), 301-304.</a>
* 
* @param {Array.<number>} x an array of real numbers.
* @param {number} k a strictly positive natural integer specifying which k-th smallest element of x is to be selected.
* @param {function} compareFunction an optional sort function that defines the sort order, using the standard prototype for JavaScript sort functions (c.f. https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/sort).
* @return {number} the k-th smallest element of x
*
* @example
* select_([2,4,1], 2);
* // 2
* // [2,4,1] is permuted into [1,2,4]
*/
function select_(x, k, compareFunction) {
	// ------
	
	// Initialisations.
	var defaultCompareFct = function (a, b) {
		return a - b;
	};
	var compareFunction = compareFunction || defaultCompareFct;
	
	var n = x.length;

	var cutoff = 600;
	var cs = 0.5; // Brown's version: cs = 0.5
	var csd = 0.5; // Brown's version: cs = 0.1

	// The arrays stack_1 and stack_2 of nstack elements permits up to
	// nstack levels of recursion.
    // For standard parameters cs <= 1 and cutoff >= 600,
    // nstack = 5 suffices for n up to 2**31-1 (maximum integer*4).	
	var nstack = 10;
	var stack_1 = typeof UInt32Array === 'function' ? new UInt32Array(nstack) : new Array(nstack);
	var stack_2 = typeof UInt32Array === 'function' ? new UInt32Array(nstack) : new Array(nstack);
	var jstack = 0; // number of elements in the stacks stack_1 and stack_2
	
	var l = 0;
    var r = n - 1; // -1 because Fortran convention is to start the arrays at index 1
    var k = k - 1; // same as above
	
	
	// ------
	
	// entry to SELECT( x, n, l, r, k)
	// SELECT will rearrange the values of the array segment x[l:r] so
	// that x(k) (for some given k; 0 <= k <= r-1) will contain the
	// (k-l+1)-th smallest value, l <= i <= k will imply x(i) <= x(k),
	// and k <= i <= r will imply x(k) <= x(i).
	while (true) {
		// Note: Rules of FORTRAN 77 rounding of real numbers to integers can
		// be found here -> https://gcc.gnu.org/onlinedocs/gfortran/INT.html
		
		// The additional test below prevents stack overflow.
		if (r - l > cutoff &&  jstack < nstack) {
			// Use SELECT recursively on a sample of size s to get an
			// estimate for the (k-l+1)-th smallest element into x(k),
			// biased slightly so that the (k-l+1)-th element is
			// expected to lie in the smaller set after partitioning.
			var m = r - l + 1;
			var i = k - l + 1;
			var dm = m;
			
			var z = Math.log(dm);
			var s = Math.floor(cs * Math.exp(2*z/3) + 0.5); // from the code, s is a positive integer
			var sign = i >= dm/2 ? 1 : - 1 // emulates sign(1,i-dm/2)
			var sd = csd * Math.sqrt(z*s*(1-s/dm)) * sign + 0.5; // sd is supposed to be an integer, and can be positive or negative, so, emulates FORTRAN rounding
			if (-1 < sd && sd < 1) {
				sd = 0;
			}
			else if (sd >= 1) {
				sd = Math.floor(sd);
			}
			else {
				sd = Math.ceil(sd);
			}
			// Brown's modification: sd = csd*Math.sqrt(z*s*(1-s/dm))*(2*i/dm-1)+0.5;
			if (i == m/2) {
				sd = 0;
			}
			
			// Push the current l and r on the stack.
			stack_1[jstack] = l;
			stack_2[jstack] = r;
			jstack++;
			
			// Find new l and r for the next recursion.
			var comp = k - i*(s/dm) + sd;
			if (l < comp) {
				l = Math.floor(comp + 0.5); // l is a positive integer
			}
			if (r > comp + s) {
				r = Math.floor(comp + s + 0.5); // r is a positive integer
			}
			// call SELECT( x, n, l, r, k)
		}
		else {
			if (l >= r) {
				// Exit if the stack is empty.
				if (jstack == 0) {
					return x[k];
				}
				
				// Pop l and r from the stack.
				--jstack;
				l = stack_1[jstack];
				r = stack_2[jstack];
				
				// Continue as if after a return from a recursive call.
			}
			
			// Partition x[l:r] about the pivot v := x(k).
			var v = x[k];
			
			// Initialize pointers for partitioning.
			i = l;
			j = r;
			
			// Swap x(l) and x(k).
			x[k] = x[l];
			x[l] = v;

			//if (v < x[r]) {
			if (compareFunction(v, x[r]) < 0) {
				// Swap x(l) and x(r).
				x[l] = x[r];
				x[r] = v;
			}
			
			while (i < j) {
	            //Swap x(i) and x(j).
				var tmp = x[j];
	            x[j] = x[i];
	            x[i] = tmp;
				
				++i;
				--j;
				
				// Scan up to find element >= v.
	            //while (x[i] < v) {
				while (compareFunction(x[i], v) < 0) {
					++i;
				}
				
				// Scan down to find element <= v.
				//while (x[j] > v) {
				while (compareFunction(x[j], v) > 0) {
					--j;
				}
			}
			
			//if (x[l] == v) {
			if (compareFunction(x[l], v) == 0) {
				// Swap x(l) and x(j).
				var tmp = x[l];
				x[l] = x[j];
				x[j] = tmp;
			} 
			else {
				++j;
				
				// Swap x(j) and x(r).
				var tmp = x[j];
				x[j] = x[r];
				x[r] = tmp;
			}
			
			// Now adjust l, r so that they surround the subset containing
			// the (k-l+1)-th smallest element.
			if (j <= k) {
				l = j + 1;
			}
			if (k <= j) {
				r = j - 1;
			}
		}
	}
}
 
  /**
* @function nextUp_
*
* @summary Returns the next double-precision number larger than a number.
*
* @description This function computes the next double-precision number
* larger than a number x.
*
* This function has been copied/pasted from https://gist.github.com/Yaffle/4654250,
* with no adaptation.
*
* @param {number} x a real number.
* @return {number} the next double-precision number larger than x, a real number.
*
* @example
* nextUp_(1.0000000000000002);
* // 
*/
function nextUp_(x) {
	var EPSILON = Math.pow(2, -52);
	var MAX_VALUE = (2 - EPSILON) * Math.pow(2, 1023);
	var MIN_VALUE = Math.pow(2, -1022);

	if (x !== x) {
	  return x;
	}
	if (x === -1 / 0) {
	  return -MAX_VALUE;
	}
	if (x === +1 / 0) {
	  return +1 / 0;
	}
	if (x === +MAX_VALUE) {
	  return +1 / 0;
	}
	var y = x * (x < 0 ? 1 - EPSILON / 2 : 1 + EPSILON);
	if (y === x) {
	  y = MIN_VALUE * EPSILON > 0 ? x + MIN_VALUE * EPSILON : x + MIN_VALUE;
	}
	if (y === +1 / 0) {
	  y = +MAX_VALUE;
	}
	var b = x + (y - x) / 2;
	if (x < b && b < y) {
	  y = b;
	}
	var c = (y + x) / 2;
	if (x < c && c < y) {
	  y = c;
	}
	return y === 0 ? -0 : y;
}


 /**
* @function hypot_
*
* @summary Returns the square root of the sum of the squares of two numbers (i.e., the hypotenuse).
*
* @description This function computes the value of sqrt(abs(x)^2 + abs(y)^2) in a way
* to avoid as much as possible underflow and overflow.
*
* @see <a href="https://en.wikipedia.org/wiki/Hypot#Implementation">Hypot</a>
*
* @param {number} x a real number.
* @param {number} y a real number.
* @return {number} the value of sqrt(abs(x)^2 + abs(y)^2), a real number.
*
* @example
* hypot_(3, 4);
* // 5
*/
function hypot_(x, y) {
    // Initialization
	var r = 0;
    
	// Main algorithm
	var absX = Math.abs(x);
	var absY = Math.abs(y);
	if (absX > absY) {
	   r = y/x;
	   r = absX * Math.sqrt(1 + r*r);
    } 
	else if (y != 0) {
	   r = x/y;
	   r = absY * Math.sqrt(1 + r*r);
    }
	else {
	   r = 0;
    }
    
	// Return the computed value
	return r;
}

 
/**
* @function rank_
*
* @summary Returns the rank of each value in a serie of values.
*
* @description This function computes the rank of each value in a serie of values, which is computed 
* by first sorting the serie of values, either by ascending or descending order, and then by computing 
* the position of each value in the sorted serie.
*
* Duplicate values in the serie of values all have the same rank, defined as the bottom rank of these duplicate values.
*
* This function mimics the Excel function RANK.EQ.
*
* @param {Array.<number>} x an array of real numbers.
* @param {number} order an integer equals to 0 to sort the serie of values in descending order, or equals to 1 to sort the serie of values in ascending order.
* @return {Array.<number>} an array of real numbers of the same size as x, containing the rank of each value of x.
*
* @example
* rank_([12, 13, 15, 10, 12], 1);
* // [2, 4, 5, 1, 2]
*/
function rank_(x, order) {
	// Transform the input array into an array with indexes
	var xWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		xWithIndexes[i] = [x[i], i];
	}
	
	// Sort the transformed array
	if (order == 0) {
		xWithIndexes.sort(function(a, b) {
			return a[0] > b[0] ? -1 : 1;
		}); 
	}
	else if (order == 1) {
		xWithIndexes.sort(function(a, b) {
			return a[0] < b[0] ? -1 : 1;
		}); 
	}
	
	// Compute the ranks of the values, setting an equal rank for all identical values
	// and skipping the next ranks values
	var xRanks = new Array(x.length);
	xRanks[xWithIndexes[0][1]] = 1; // first rank is always 1
	for (var i = 1; i < x.length; ++i) {
	    if (xWithIndexes[i][0] == xWithIndexes[i-1][0]) {
	  	  xRanks[xWithIndexes[i][1]] = xRanks[xWithIndexes[i-1][1]];
	    }
	    else {
		  xRanks[xWithIndexes[i][1]] = i + 1;
	    }
	}
	
	// Returnt the computed ranks
	return xRanks;
}


/**
* @function ftca_
*
* @summary Compute a clustering representation of a universe of elements using the Fast Threshold Clustering Algorithm (FTCA).
*
* @description This function returns the clustering representation of a universe of n elements based on their correlation
* and a correlation threshold, as computed by the FTCA algorithm of David Varadi described in the reference, initially created
* to deal with financial assets.
*
* This algorithm has many desirable properties that traditional clustering algorithms do not:
* 1) it produces fairly stable clusters
* 2) it is fast and deterministic 
* 3) it is easy to understand. 
*
* By default, David Varadi used a correlation threshold of 0.5 (approximately the level of statistical significance) to separate similar from dissimilar elements (assets).
* The choice of the threshold will change the number and stability of the clusters, with higher thresholds showing more clusters and a greater change in membership than lower thresholds. 
*
* To be noted that the FTCA works similar to the Minimum Correlation Algorithm from the same author in that it uses the average correlation of each element (asset)
* to all other elements (asset)s as a means of determining how closely or distantly related an element (asset) is to the universe of elements (assets) chosen.
*
* @see <a href="https://cssanalytics.wordpress.com/2013/11/26/fast-threshold-clustering-algorithm-ftca/">Fast Threshold Clustering Algorithm (FTCA)</a>
* 
* @param {Array.<Array.<number>>} correlationMatrix the correlation matrix (rho_ij),i,j=1..n, an array arr of n arrays of n real numbers 
* satisfying arr[i-1][j-1] = rho_ij, i,j=1..n, where n is a strictly positive natural integer.
* @param {number} threshold the correlation threshold to use in the FTCA algorithm, a real number typically belonging to interval [-1, 1].
* @return {Array.<Array.<number>>} the list of clusters as computed by the FTCA algorithm, array of m arrays of strictly positive integers representing the indexes of the elements in the considered universe, where m is the number of clusters, with the m arrays forming a partition of the set [1..n].
*
* @example
* ftca_([[1, 0], [0,1]]), 0.5);
*  // [[2],[1]]
*/
function ftca_(correlationMatrix, threshold) {
	// Decode the optional threshold
	var threshold = threshold;
	if (threshold === undefined) {
		threshold = 0.5;
	}
	
	// Convert the correlation matrix to matrix format
	var correlationMatrix = new Matrix_(correlationMatrix);
	
	// The list of output clusters, to be populated
	var clusters = [];

	// The list of elements indexes not assigned to any cluster, initialized with all elements indexes (initially, no clusters are existing)
	var nbElements = correlationMatrix.nbRows;
	var unassignedElementsIdx = new Array(nbElements);
	for (var i = 0; i < unassignedElementsIdx.length; ++i) {
		unassignedElementsIdx[i] = i + 1;
	}

	// While there are elements that have not been assigned to a cluster
	while (unassignedElementsIdx.length != 0) {
		// If only one element remaining then
		if (unassignedElementsIdx.length === 1) {
			// Add a new cluster
			// Only member is the remaining element, set as not unassigned anymore
			var newCluster = [unassignedElementsIdx[0]];
			unassignedElementsIdx[0] = null;
			
			// Effectively add the new cluster into the list of clusters
			clusters.push(newCluster);
		}	
		else {
			// Get the (sub)correlation matix of the unassigned elements
			var subCorrMat = correlationMatrix.submatrix(unassignedElementsIdx, unassignedElementsIdx);
			
			// Compute the average correlation of each unassigned element to all the other unassigned elements
			// Computation is done for each row
			var subCorrMatRows = subCorrMat.toRowArray(function(i, j, val) {
				return i != j;
			});
			var avgCorrelation = new Array(unassignedElementsIdx);
			for (var i = 0; i < unassignedElementsIdx.length; ++i) {
				avgCorrelation[i] = mean_(subCorrMatRows[i]);
			}
				
			// Find the element with the Highest Average Correlation (HC) to all elements not yet been assigned to a Cluster
			// Find the element with the Lowest Average Correlation (LC) to all elements not yet assigned to a Cluster
			// Note: When only 2 elements are remaining, HC will be equal to LC
			var hc = 0;
			var hcIdx = -1;
			var highestAvgCorrelation = -1;
			var lc = 0;
			var lcIdx = -1;
			var lowestAvgCorrelation = 1;		
			for (var i = 0; i < unassignedElementsIdx.length; ++i) {
				if (avgCorrelation[i] >= highestAvgCorrelation) {
					hc = unassignedElementsIdx[i];
					hcIdx = i;
					highestAvgCorrelation = avgCorrelation[i];
				}
				if (avgCorrelation[i] <= lowestAvgCorrelation) {
					lc = unassignedElementsIdx[i];
					lcIdx = i;
					lowestAvgCorrelation = avgCorrelation[i];
				}
			}
			
			// If Correlation between HC and LC > Threshold
			if (correlationMatrix.getValueAt(hc, lc) > threshold) {
				// Add a new Cluster made of HC and LC and set these two elements as not unassigned anymore
				// (Unless HC == LC, which can happen, for instance when there are only two elements remaining)
				var newClusterHcLc = (hc === lc ? [hc] : [hc, lc]);
				unassignedElementsIdx[hcIdx] = null;
				unassignedElementsIdx[lcIdx] = null;
				
				// Add to Cluster all other elements that have yet been assigned to a Cluster and have an Average Correlation to HC and LC > Threshold
				// Note: In Systematic Investor R code, all remaining elements are put inthe HcLc cluster, disregarding the condition on the correlation above.
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned elements (HC and LC)
						var avgHcLcAssetCorrelation = (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) + correlationMatrix.getValueAt(unassignedElementsIdx[i], lc)) / 2;
						if (avgHcLcAssetCorrelation  > threshold) {
							newClusterHcLc.push(unassignedElementsIdx[i]);
							
							// Set the element as not unassigned anymore
							unassignedElementsIdx[i] = null;
						}
					}
				}
			   
				// Effectively add the new cluster into the list of clusters				
				clusters.push(newClusterHcLc);
			}
			// Else
			else {
				// Add a Cluster made of HC and set this element as not unassigned anymore
				var newClusterHc = [hc];
				unassignedElementsIdx[hcIdx] = null;
				
				// Add to Cluster all other assets that have yet been assigned to a Cluster and have a Correlation to HC > Threshold
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC)
						if (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) > threshold) {
							newClusterHc.push(unassignedElementsIdx[i]);
							
							// Set the element as not unassigned anymore
							unassignedElementsIdx[i] = null;
						}
					}
				}
				
				// Effectively add the new cluster into the list of clusters				
				clusters.push(newClusterHc);

				// Add a Cluster made of LC and set this element as not unassigned anymore
				// (Unless HC == LC, which can happen, for instance when there are only two elements remaining)
				if (hc !== lc) {
					// Note: At this stage, the LC element cannot have been assigned to the Hc cluster above if LC <> HC, since
					// otherwise, it would mean corr(lc, hc) > threshold, which is incompatible with the "else" branch in which
					// the code currently is; Lc cluster is thus always non empty.
					var newClusterLc = [lc];
					unassignedElementsIdx[lcIdx] = null;
					
					// Add to Cluster all other assets that have yet been assigned to a Cluster and have Correlation to LC > Threshold
					for (var i = 0; i < unassignedElementsIdx.length; ++i) {
						if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC with its correlated assets, and LC)
							if (correlationMatrix.getValueAt(unassignedElementsIdx[i], lc) > threshold) {
								newClusterLc.push(unassignedElementsIdx[i]);
								
								// Set the element as not unassigned anymore
								unassignedElementsIdx[i] = null;
							}
						}
					}

					// Effectively add the new cluster into the list of clusters				
					clusters.push(newClusterLc);
				}
				
				// Note: In Systematic Investor R code, it is possible for an element to belong to the two clusters Hc and Lc, in which case Hc is the final cluster,
				// which is conform to the description of the reference.
			}
		}
		
		// Effectively remove the assigned elements indexes (now pointing to null)  from the list of unassigned elements
		var newUnassignedElementsIdx = [];
		for (var i = 0; i < unassignedElementsIdx.length; ++i) {
			if (unassignedElementsIdx[i] !== null) {
				newUnassignedElementsIdx.push(unassignedElementsIdx[i]);
			}
		}
		unassignedElementsIdx = newUnassignedElementsIdx;
	}

	// Return the computed list of clusters
	return clusters;
}



/**
* @function mean_
*
* @summary Compute the arithmetic mean of a serie of values.
*
* @description This function returns the arithmetic mean of a serie of values [x_1,...,x_p], 
* which is defined as the sum of the p values x_1,...,x_p, divided by p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
* 
* @param {Array.<number>} x an array of real numbers.
* @return {number} the arithmetic mean of the values of the array x.
*
* @example
* mean_([2,4]); 
* // 3
*/
function mean_(x) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the values of the input numeric array (first pass)
	var tmpMean = 0.0;
	var sum = 0.0;
	for (var i=0; i<nn; ++i) {
		sum += x[i];
	}
	tmpMean = sum/nn;

	// Compute the correction factor (second pass)
	// C.f. M_3 formula of the reference
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		sumDiff += (x[i] - tmpMean);
	}

	// Return the corrected mean
	return (sum + sumDiff)/nn;
}


/**
* @function variance_
*
* @summary Compute the variance of a serie of values.
*
* @description This function returns the variance of a serie of values [x_1,...,x_p], 
* which is defined as the arithmetic mean of the p values (x_1-m)^2,...,(x_p-m)^2, where m is the arithmetic mean
* of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the variance of the values of the array x.
*
* @example
* variance_([4, 7, 13, 16]); 
* // 22.5
*/
function variance_(x) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the input numeric array (first pass)
	var meanX = mean_(x);

	// Compute the squared deviations plus the correction factor (second pass)
	// C.f. S_4 formula of the reference
	var sumSquareDiff = 0.0;
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		var diff = (x[i] - meanX);
		sumSquareDiff += diff * diff;
		sumDiff += diff;
	}

	// Compute the corrected sum of squares of the deviations from the mean
	var S = sumSquareDiff - ((sumDiff * sumDiff) / nn);

	// Return the corrected variance
	return S/nn;
}


/**
* @function sampleVariance_
*
* @summary Compute the sample variance of a serie of values.
*
* @description This function returns the sample variance of a serie of values [x_1,...,x_p], 
* which is defined as the variance of the p values x_1,...,x_p multiplied by p/(p-1).
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function variance_.
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the variance of the values of the array x.
*
* @example
* sampleVariance_([4, 7, 13, 16]); 
* // 30
*/
function sampleVariance_(x) {
	var nn = x.length;
	return variance_(x) * nn/(nn - 1);
}


/**
* @function stddev_
*
* @description Compute the standard deviation of a serie of values.
*
* @description This function returns the standard deviation of a serie of values [x_1,...,x_p], 
* which is defined as the square root of the variance of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function variance_.
*
* @see <a href="https://en.wikipedia.org/wiki/Standard_deviation">https://en.wikipedia.org/wiki/Standard_deviation</a>
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the standard deviation of the values of the array x.
*
* @example
* stddev_([1, 2, 3, 4]); 
* // ~1.12
*/
function stddev_(x) {
	return Math.sqrt(variance_(x));
}


/**
* @function sampleStddev_
*
* @description Compute the sample standard deviation of a serie of values.
*
* @description This function returns the sample standard deviation of a serie of values [x_1,...,x_p], 
* which is defined as the square root of the sample variance of the p values x_1,...,x_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function sampleVariance_.
*
* @param {Array.<number>} x an array of real numbers.
* @return {number} the standard deviation of the values of the array x.
*
* @example
* sampleStddev_([1, 2, 3, 4]); 
* // ~1.29
*/
function sampleStddev_(x) {
	return Math.sqrt(sampleVariance_(x));
}



/**
* @function normcdf_
*
* @summary Compute the the standard normal cumulative distribution function.
*
* @description This function returns an approximation of the standard normal cumulative distribution function, i.e.
* given x a real number, it returns an approximation to p = Pr{Z <= x} where Z is a
* random variable following a standard normal distribution law.
*
* This function is also called Phi in the statistical litterature.
*
* The algorithm uses a Taylor expansion around 0 of a well chosen function of Phi,
* and has a theorical absolute error of less than 8e−16.
*
* @author George Marsaglia
*
* @see <a href="https://www.jstatsoft.org/article/view/v011i04/v11i04.pdf"> G. Marsaglia. Evaluating the normal distribution. Journal of Statistical Software, 11(4):1–11, 2004.</a>
* 
* @param {number} x a real number.
* @return {number} an approximation to the p value satisfying p = Pr{Z <= x} where Z is a random variable following a standard normal distribution law.
*
* @example
* normcdf_(0);
* // 0.5
*/
function normcdf_(x) {
	// Initialisations
	var sum = x;
	var term = 0;
	var next_term = x;
	var power = x*x;
	var i = 1;

	// Limit cases, as described in the reference.
	if (x < -8.0) {
		return 0.0;
	}
	else if (x > 8.0) {
		return 1.0;
	}
	
	// The main loop corresponds to the computation of the Taylor serie of the function B around 0, 
	// c.f. page 5 of the reference.
	//
	// In a nutshell, the Taylor expansion is computed term by term until the addition of a new term 
	// stops to produce a change (from a numerical accuracy perspective).
	while (sum != term) {
		sum = (term = sum) + (next_term *= power/(i += 2));
	}

	// The formula linking Phi and the Taylor expansion above if Phi = 1/2 + normal density * B, c.f. page 5 of the reference.
	return 0.5 + sum * Math.exp(-0.5 * power - 0.91893853320467274178)
}


/**
* @function norminv_
*
* @summary Compute the inverse of the standard normal cumulative distribution function.
*
* @description This function returns an approximation of the inverse standard normal cumulative distribution function, i.e.
* given p in [0,1] it returns an approximation to the x value satisfying p = Pr{Z <= x} where Z is a
* random variable following a standard normal distribution law.
*
* x is also called a z-score.
*
* The algorithm uses three separate rational minimax approximations: one rational approximation is used for the central region and the outer region is split into two sub-regions.
* The algorithm has a relative error whose absolute value is in the order of 1e-15.
*
* @see <a href="https://www.jstor.org/stable/2347330">Michael Wichura, Algorithm AS 241: The Percentage Points of the Normal Distribution., Applied Statistics, Volume 37, Number 3, pages 477-484, 1988.</a>
* 
* @param {number} p a probability value, real number belonging to interval [0,1].
* @return {number} an approximation to the x value satisfying p = Pr{Z <= x} where Z is a random variable following a standard normal distribution law.
*
* @example
* norminv_(0.5);
* // 0
*/
function norminv_(p) {
    // Checks
	if (p <= 0.0 || p >= 1.0) {
		if (p == 0.0) {
			return -Infinity;
		}
		else if (p == 1.0) {
			return Infinity;
		}
		else {
			throw "The probality p must be bigger than 0 and smaller than 1";
		}
	}

    var q = p - 0.5;

	var ppnd16;
	
    if (Math.abs(q) <= 0.425) { // P CLOSE TO 1/2
        var r = 0.180625 - q * q;
        var num_ppnd16 = ((((((r * 2.5090809287301226727E3 + 3.3430575583588128105E4) * r + 6.7265770927008700853E4) * r + 4.5921953931549871457E4) * r + 1.3731693765509461125E4) * r + 1.9715909503065514427E3) * r + 1.3314166789178437745E2) * r + 3.3871328727963666080E0;
		var denom_ppnd16 = ((((((r * 5.2264952788528545610E3 + 2.8729085735721942674E4) * r + 3.9307895800092710610E4) * r + 2.1213794301586595867E4) * r + 5.3941960214247511077E3) * r + 6.8718700749205790830E2) * r + 4.2313330701600911252E1) * r + 1.0;
		ppnd16 = q * num_ppnd16 / denom_ppnd16;
    }
    else {
		var r;
		if ( q < 0.0 ) {
		  r = p;
		}
		else {
		  r = 1.0 - p;
		}

		//if ( r <= 0.0 ) {
		// No need for this check, as it has already been done at the beginning of the function
		//}
		
        r = Math.sqrt(-Math.log(r));

        if (r <= 5) { // P NEITHER CLOSE TO 1/2 NOR 0 OR 1
            r = r - 1.6;
			var num_ppnd16 = ((((((r * 7.74545014278341407640E-4 + 2.27238449892691845833E-2) * r + 2.41780725177450611770E-1) * r + 1.27045825245236838258E0) * r + 3.64784832476320460504E0) * r + 5.76949722146069140550E0) * r + 4.63033784615654529590E0) * r + 1.42343711074968357734E0;
			var denom_ppnd16 = ((((((r * 1.05075007164441684324E-9 + 5.47593808499534494600E-4) * r + 1.51986665636164571966E-2) * r + 1.48103976427480074590E-1) * r + 6.89767334985100004550E-1) * r + 1.67638483018380384940E0) * r + 2.05319162663775882187E0) * r + 1.0;
			ppnd16 = num_ppnd16 / denom_ppnd16;
        }
        else { // COEFFICIENTS FOR P NEAR 0 OR 1
            r = r - 5.0;
			var num_ppnd16 = (((((((r * 2.01033439929228813265E-7 + 2.71155556874348757815E-5) * r + 1.24266094738807843860E-3) * r + 2.65321895265761230930E-2) * r + 2.96560571828504891230E-1) * r + 1.78482653991729133580E0) * r + 5.46378491116411436990E0) * r + 6.65790464350110377720E0);
			var denom_ppnd16 = (((((((r * 2.04426310338993978564E-15 + 1.42151175831644588870E-7) * r + 1.84631831751005468180E-5) * r + 7.86869131145613259100E-4) * r + 1.48753612908506148525E-2) * r + 1.36929880922735805310E-1) * r + 5.99832206555887937690E-1) * r + 1.0);
			ppnd16 = num_ppnd16 / denom_ppnd16;
        }

        if (q < 0.0) {
            ppnd16 = -ppnd16;
        }
    }

	// Return the computed value
    return ppnd16;
}


/**
* @function hypersphereRandomSampler_
*
* @summary Returns a function to compute random points on the unit hypersphere of R^n.
*
* @description This function constructs a function to compute random points uniformly distributed on
* the unit hypersphere of R^n, using the algorithm of the reference.
* 
* @see <a href="https://dl.acm.org/citation.cfm?id=377946">	Mervin E. Muller, A note on a method for generating points uniformly on n-dimensional spheres, Communications of the ACM CACM Homepage archive
Volume 2 Issue 4, April 1959 Pages 19-20 </a>
*
* @param {number} n the dimension of the unit hypersphere of R^n, natural integer superior or equal to 1.
* @return {function} a function to be used through its .sample() method, computing random  
* points on the unit hypersphere of R^n.
*
* @example
* var mySampler = new hypersphereRandomSampler_(3);
* mySampler.sample();
* // [1, 0, 0]
*/
function hypersphereRandomSampler_(n) {
	// Initializations
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	
	/**
	* @function sample
	*
	* @summary Returns a random point on the unit hypersphere of R^n.
	*
	* @description This function computes a point choosen uniformly at random on the unit hypersphere of R^n,
	* using the algorithm of the reference.
	*
	* @memberof hypersphereRandomSampler_
	* @return {Array.<number>|Float64Array} an array of n real numbers corresponding to the coordinates of the computed point in R^n.
	*
	*/
	this.sample = function() {
		// Computation of n independent random variables from N(0,1), which will form the basis
		// of the coordinates of the point being sampled.
		var sum_sq = 0;
		for (var i = 0; i < this.n; ++i) {
			// Generate a random variable from N(0,1), using the inverse method
			var u = Math.random(); // u ~ U[0,1[
			while (u === 0.0) {
				u = Math.random();
			} // u ~ U]0,1[
			var r = norminv_(u); // r ~ N(0,1)
			
			// Set the i-th coordinate of the point being sampled.
			this.x[i] = r;
			
			// Compute the running sum of the squares of the normal variables, for the subsequent normalization step.
			sum_sq += r * r;
		}

		// Normalization of the computed coordinates of the point being sampled, so that
		// the 2-norm of the associated vector in R^n is equal to 1.
		var x_two_norm = Math.sqrt(sum_sq);
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/x_two_norm;
		}
		
		// Return a copy of the point being sampled, so that callers can alter it.
		return this.x.slice();
	}
}


/**
* @function covariance_
*
* @summary Compute the covariance of two serie of values.
*
* @description This function returns the covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p], 
* which is defined as the arithmetic mean of the p values (x_1-m_x)*(y_1-m_y),...,(x_p-m_x)*(y_p-m_y), 
* where m_x is the arithmetic mean of the p values x_1,...,x_p and m_y is the arithmetic mean of the p values y_1,...,y_p.
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the reference.
*
* @see <a href="http://dl.acm.org/citation.cfm?doid=365719.365958">Peter M. Neely (1966) Comparison of several algorithms for computation of means, standard deviations and correlation coefficients. Commun ACM 9(7):496–499.</a>
*
* @param {Array.<number>} x an array of real numbers.
* @param {Array.<number>} y an array of real numbers of the same length as x.
* @return {number} the covariance of the values of the arrays x and y.
*
* @example
* covariance_([4, 7, 13, 16], [4, 7, 13, 16]); 
* // 22.5
*/
function covariance_(x, y) {
	// Initialisations
	var nn = x.length;

	// Compute the mean of the input numeric arrays (first pass)
	var meanX = mean_(x);
	var meanY = mean_(y);

	// Compute the sum of the product of the deviations plus the correction factor (second pass)
	// C.f. P_4 formula of the reference
	var sumProdDiff = 0.0;
	var sumDiffX = 0.0;
	var sumDiffY = 0.0;
	for (var i=0; i<nn; ++i) {
		var diffX = (x[i] - meanX);
		var diffY = (y[i] - meanY);
		sumProdDiff += diffX * diffY;
		sumDiffX += diffX;
		sumDiffY += diffY;
	}

	// Compute the corrected sum of the product of the deviations from the means
	var C = sumProdDiff - ((sumDiffX * sumDiffY) / nn);

	// Return the corrected covariance
	return C/nn;
}


/**
* @function sampleCovariance_
*
* @summary Compute the sample covariance of two serie of values.
*
* @description This function returns the sample covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p], 
* which is defined as the covariance of two series of values [x_1,...,x_p] and [y_1,...,y_p] multiplied by p/(p-1).
*
* The algorithm implemented uses a two pass formula in order to reduce the computation error, c.f. the function covariance_.
*
* @param {Array.<number>} x an array of real numbers.
* @param {Array.<number>} y an array of real numbers of the same length as x.
* @return {number} the covariance of the values of the arrays x and y.
*
* @example
* sampleCovariance_([4, 7, 13, 16], [4, 7, 13, 16]); 
* // 30
*/
function sampleCovariance_(x, y) {
	var nn = x.length;
	return covariance_(x,y) * nn/(nn - 1);
}

