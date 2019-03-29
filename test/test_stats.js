// ------------------------------------------------------------
QUnit.module('Statistics internal module', {
  before: function() {
    // 
  }
});



QUnit.test('Median computation', function(assert) {    
  // Test with static data
  {
	  assert.equal(PortfolioAllocation.median_([2,4,1]), 2, 'Median computation #1');
	  assert.equal(PortfolioAllocation.median_([2,4,1,3]), 2, 'Median computation #2');
  }  
  
  //TODO: use random data
});



QUnit.test('Max computation', function(assert) {    
  // Test with static data
  {
	  assert.deepEqual(PortfolioAllocation.max_([2,4,4,1]), [4,1], 'Max computation #1');
  }  
  
  //TODO: use random data
});

QUnit.test('Smallest k element computation', function(assert) {    
	function generateRandomDimension(min, max) { // used for n
		return Math.floor(Math.random()*(max-min+1) + min);
	}
	
	function generateRandomValue(minVal, maxVal) { // used for each array element		
		return Math.random() * (maxVal - minVal) + minVal;
	}
	
	function generateRandomArray(n) { // used for arr
		var minVal = -10;
		var maxVal = 10;
		
		var arr = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		
		for (var i = 0; i < n; ++i) {
			arr[i] = generateRandomValue(minVal, maxVal);
		}
		
		return arr;
	}
	
  // Test with static data, to make sure duplicate values are properly handled
  {
	  // Array with duplicate values
	  var val = [0.5, 0.1, 0.2, 0.2, 0.3, 0.2, 0.3];
	  var val_idx = [0, 1, 2, 3, 4, 5, 6];
	  var expected_val_idx = [0, 4, 6, 5, 2, 1, 3];
	  
	  // Call the SELECT algorithm to partition the indexes of the above array.
	  var compareIndexes = function (a, b) {
		  return val[b] - val[a];
	  };
	  PortfolioAllocation.select_(val_idx, 3, compareIndexes);

	  // In case of incorrect duplicate values management, the SELECT algorithm
	  // will output [0, 6, 6, 5, 2, 1, 3]
	  assert.deepEqual(val_idx, expected_val_idx, 'Smallest k element computation #0');
  }
  
  // Test with random data, no indices requested in output
  {
	var nbTests = 100;
	for (var j = 0; j < nbTests; ++j) {
	  // Generate a random array of a random size between 1 and 1000
	  var n = generateRandomDimension(1, 1000);
	  var arr = generateRandomArray(n);
		
	  // Sort (ascending) a copy of the array, to be used for the expected values
	  var copy_arr = arr.slice().sort(function(a, b) { return a - b; });
	  
	  // Generate a random integer k between 1 and the size of the array
	  var k = generateRandomDimension(1, n);
	  
	  // Compute the k-th smallest element of the array arr using the function
	  // SELECT
	  var k_smallest_elem = PortfolioAllocation.select_(arr, k);
	  
	  // Test that the k-th smallest element of the array arr, computed
	  // from the function SELECT, is the same as arr[k-1].
	  assert.equal(k_smallest_elem, arr[k-1], 'Smallest k element computation #1 ' + (j+1));

	  // Test that the k-th smallest element of the array arr, computed
	  // from the function SELECT, is the same as the k-th element of 
	  // the sorted array copy_arr.
	  assert.equal(k_smallest_elem, copy_arr[k-1], 'Smallest k element computation #2 ' + (j+1));
	  
	  // Test that smallest k elements of arr are x[i], i=0..k-1
	  var correct_order_left_k = true;
	  for (var i = 0; i < k; ++i) {
		if (arr[i] > k_smallest_elem) {
			correct_order_left_k = false;
			break;
		}
	  }
	  assert.equal(correct_order_left_k, true, 'Smallest k element computation #3 ' + (j+1));
	  
	  // Test that largest n-k elements of arr are x[i], i=k+1..n-1
	  var correct_order_right_k = true;
	  for (var i = k; i < n; ++i) {
		if (arr[i] < k_smallest_elem) {
			correct_order_right_k = false;
			break;
		}
	  }
	  assert.equal(correct_order_right_k, true, 'Smallest k element computation #4 ' + (j+1));
    }
  }
  
  // Test with random data, using a custom comparator to extract indexes
  {
	var nbTests = 100;
	for (var j = 0; j < nbTests; ++j) {
	  // Generate a random array of a random size between 1 and 1000
	  var n = generateRandomDimension(1, 1000);
	  var arr = generateRandomArray(n);
	  var arr_copy = arr.slice();
	  
	  // Generate the associated indexes
	  var arr_idx = typeof UInt32Array === 'function' ? new UInt32Array(n) : new Array(n);
	  for (var i = 0; i < n; ++i) {
		arr_idx[i] = i;
	  }
	   
	  // Generate a random integer k between 1 and the size of the array
	  var k = generateRandomDimension(1, n);
	  
	  // Compute the index of the k-th smallest element of the array arr using the function
	  // SELECT and a custom comparator
      var compareIndexes = function (a, b) {
		  return arr[a] - arr[b];
	  };
	  var k_smallest_elem_index = PortfolioAllocation.select_(arr_idx, k, compareIndexes);
	  
	  // Compute the k-th smallest element of the array arr using the function
	  // SELECT
	  var k_smallest_elem = PortfolioAllocation.select_(arr, k);
	  
	  // Test that the k-th smallest element index provided in output is correct
	  assert.equal(k_smallest_elem == arr_copy[k_smallest_elem_index], true, 'Smallest k element computation #5 ' + (j+1));
	  
	  // Test that the indexes provided in output of SELECT 
	  // corresponds to the original array elements
	  var correct_indexes = true;
	  for (var i = 0; i < n; ++i) {
		if (arr[i] != arr_copy[arr_idx[i]]) {
			correct_indexes = false;
			break;
		}
	  }
	  assert.equal(correct_indexes, true, 'Smallest k element computation #6 ' + (j+1));
    }
  }
});

QUnit.test('Hypothenuse computation', function(assert) {    
  // Tests with static data
  {
	  assert.ok(Math.abs(PortfolioAllocation.hypot_(3e200, 4e200) - 5e+200)/5e+200 <= 1e-15, 'Hypothenuse computation with no overflow');
	  assert.equal(PortfolioAllocation.hypot_(3, 4), 5, 'Hypothenuse computation #1');
	  assert.equal(PortfolioAllocation.hypot_(-2, 0), 2, 'Hypothenuse computation one zero argument #1');
	  assert.equal(PortfolioAllocation.hypot_(0, -2), 2, 'Hypothenuse computation one zero argument #2');
	  assert.equal(PortfolioAllocation.hypot_(0, 0), 0, 'Hypothenuse computation two zero arguments');
  }  
  
  // Tests with the formula and random data
  {
	  var nbTests = 50;
	  for (var i = 0; i < nbTests; ++i) {
          var x = Math.random();
		  var y = Math.random();
		  var naiveHypothenuse = Math.sqrt(x*x + y*y);
		  assert.ok(Math.abs(PortfolioAllocation.hypot_(x, y) - naiveHypothenuse) <= 1e-14, 'Hypothenuse computation #' + (i + 2));
	  }
  }
});

QUnit.test('Rank computation', function(assert) {    
  // Test with static data
  {
	  var testValues = [12, 13, 15, 10, 12];
	  var expectedRanksDescending = [3, 2, 1, 5, 3];
	  var expectedRanksAscending = [2, 4, 5, 1, 2];
	  
	  assert.deepEqual(PortfolioAllocation.rank_(testValues, 0), expectedRanksDescending, 'Rank descending');
	  assert.deepEqual(PortfolioAllocation.rank_(testValues, 1), expectedRanksAscending, 'Rank ascending');
  }  
  
  // Second test with static data
  {
	var testValues = [0.13, 0.63, 0.79]; 
	var expectedRanks = [3, 2, 1];
	assert.deepEqual(PortfolioAllocation.rank_(testValues, 0), expectedRanks, 'Rank #2');
  }
});



QUnit.test('FTCA computation', function(assert) {    
  // FTCA test, using static data
  {
	  var corrMat = [[1.0000000,  0.2378848,  0.2483431,  0.3163914,  0.1796639], [0.2378848,  1.0000000,  0.2487009, -0.1986677, -0.2165444], [0.2483431,  0.2487009,  1.0000000, -0.3179188,  0.3713964], [0.3163914, -0.1986677, -0.3179188, 1.0000000,  0.4131639],[0.1796639, -0.2165444,  0.3713964,  0.4131639,  1.0000000]];
	
      // Compute FTCA with varying thresholds
	  // A threshold of 1 will always output as many clusters as initial elements
	  var thresholds = [-0.10, 0.0, 0.10, 0.25, 0.33, 0.5, 1];
	  var expectedClusters = [ [[1,2,3,4,5]], [[1,2,3,4],[5]], [[1,2,3],[5,4]], [[1,4],[2],[5,3]], [[1],[2],[5,3],[4]], [[1],[2],[5],[3],[4]], [[1],[2],[5],[3],[4]]];
	  for (var i = 0; i < thresholds.length; ++i) {
		var clusters = PortfolioAllocation.ftca_(corrMat, thresholds[i]);
		assert.deepEqual(clusters, expectedClusters[i], "FTCA - Test 1 #" + i);
	  }	  
  }
  
  // FTCA test, limit case with static data
  {
	var corrMat = [[0.01, 0.016, 0, 0], [0.016, 0.04, 0, 0], [0, 0, 0.09, -0.06], [0, 0, -0.06, 0.16]];
	var clusters = PortfolioAllocation.ftca_(corrMat, 1);
	var expectedClusters = [[2],[4],[3],[1]];
	assert.deepEqual(clusters, expectedClusters, "FTCA - Test 2");
  }
  
});


QUnit.test('Normcdf computation', function(assert) {    
  // Boundaries
  assert.equal(PortfolioAllocation.normcdf_(0), 0.5, 'Normcdf 0'); 


  // Values taken from Wolfram Alpha, with CDF[NormalDistribution[0, 1], i/10] command, 22 digits precision requested.
  var x = [-1/10,2/10,-3/10,-4/10,-5/10,-6/10,-7/10,-8/10,-9/10,-10/10,-11/10,-12/10,-13/10,-14/10,-15/10,-16/10,-17/10,-18/10,-19/10,-20/10,-21/10,-22/10,-23/10,-24/10,-25/10,-26/10,-27/10,-28/10,-29/10,-30/10,-31/10,-32/10,-33/10,-34/10,-35/10,-36/10,-37/10,-38/10,-39/10,-40/10,-41/10,-42/10,-43/10,-44/10,-45/10,-46/10,-47/10,-48/10,-49/10,-50/10,-51/10,-52/10,-53/10,-54/10,-55/10,-56/10,-57/10,-58/10,-59/10,-60/10,-61/10,-62/10,-63/10,-64/10,-65/10,-66/10,-67/10,-68/10,-69/10,-70/10,-71/10,-72/10,-73/10,-74/10,-75/10,-76/10,-77/10,-78/10,-79/10,-80/10,-81/10,-82/10,-83/10,-84/10,-85/10,-86/10,-87/10,-88/10,-89/10,-90/10,-91/10,-92/10,-93/10,-94/10,-95/10,-96/10,-97/10,-98/10,-99/10,-100/10, 1/10,2/10,3/10,4/10,5/10,6/10,7/10,8/10,9/10,10/10,11/10,12/10,13/10,14/10,15/10,16/10,17/10,18/10,19/10,20/10,21/10,22/10,23/10,24/10,25/10,26/10,27/10,28/10,29/10,30/10,31/10,32/10,33/10,34/10,35/10,36/10,37/10,38/10,39/10,40/10,41/10,42/10,43/10,44/10,45/10,46/10,47/10,48/10,49/10,50/10,51/10,52/10,53/10,54/10,55/10,56/10,57/10,58/10,59/10,60/10,61/10,62/10,63/10,64/10,65/10,66/10,67/10,68/10,69/10,70/10,71/10,72/10,73/10,74/10,75/10,76/10,77/10,78/10,79/10,80/10,81/10,82/10,83/10,84/10,85/10,86/10,87/10,88/10,89/10,90/10,91/10,92/10,93/10,94/10,95/10,96/10,97/10,98/10,99/10,100/10];
  var p =[0.4601721627229710185346,0.5792597094391030230424,0.3820885778110473626935,0.3445782583896758332631,0.3085375387259868963623,0.2742531177500735802944,0.2419636522230730147494,0.2118553985833966855755,0.1840601253467594885542,0.1586552539314570514148,0.1356660609463826751731,0.1150696702217082680222,0.09680048458561033315201,0.08075665923377104649619,0.06680720126885806600449,0.05479929169955799396047,0.04456546275854303948743,0.03593031911292580396033,0.02871655981600179940134,0.02275013194817920720028,0.01786442056281655678392,0.01390344751349861061475,0.01072411002167580539236,0.008197535924596129444387,0.006209665325776135166978,0.004661188023718750250993,0.003466973803040668495942,0.002555130330427932801531,0.001865813300384037950310,0.001349898031630094526652,0.0009676032132183568921157,0.0006871379379158484551177,0.0004834241423837772011101,0.0003369292656768809394098,0.0002326290790355250363499,0.0001591085901575338796651,0.0001077997334773883369375,0.00007234804392511997399341,0.00004809634401760271714671,0.00003167124183311992125377,0.00002065750691254673879525,0.00001334574901590633835309,8.539905470991804195354e-6,5.412543907703859841921e-6,3.397673124730060401687e-6,2.112454702502849769124e-6,1.300807453917282059602e-6,7.933281519755946161470e-7,4.791832765903198532984e-7,2.866515718791939116738e-7,1.698267407147598273938e-7,9.964426316933481269842e-8,5.790134039964588482705e-8,3.332044848542857284776e-8,1.898956246588771938385e-8,1.071759025831090735496e-8,5.990371401063534429834e-9,3.315745978326161338027e-9,1.817507863099432371362e-9,9.865876450376981407009e-10,5.303423262948829737329e-10,2.823158037043274469652e-10,1.488228221762310961320e-10,7.768847581709830408558e-11,4.016000583859117808346e-11,2.055788909399517967613e-11,1.042097698796519370792e-11,5.230957544144587508727e-12,2.600126965638172838025e-12,1.279812543885835004384e-12,6.237844463331575105220e-13,3.010627981117437487237e-13,1.438838638157585748676e-13,6.809224890620033184947e-14,3.190891672910896227767e-14,1.480653749004804708609e-14,6.803311540773970775939e-15,3.095358771958695450558e-15,1.394517146659268252841e-15,6.220960574271784123516e-16,2.747959392398220468854e-16,1.201935154273578710970e-16,5.205569744890285157996e-17,2.232393197288050341136e-17,9.479534822203318354151e-18,3.985804962848169561416e-18,1.659420869964773833918e-18,6.840807685935589292706e-19,2.792334374939655566956e-19,1.128588405953840647736e-19,4.516591491435442087441e-20,1.789748812014040789999e-20,7.022284240441672958777e-21,2.728153571346126105605e-21,1.049451507536260749283e-21,3.997221205726226709960e-22,1.507493168810194374794e-22,5.629282311376572223463e-23,2.081375219493213518484e-23,7.619853024160526065973e-24, 0.5398278372770289814654,0.5792597094391030230424,0.6179114221889526373065,0.6554217416103241667369,0.6914624612740131036377,0.7257468822499264197056,0.7580363477769269852506,0.7881446014166033144245,0.8159398746532405114458,0.8413447460685429485852,0.8643339390536173248269,0.8849303297782917319778,0.9031995154143896668480,0.9192433407662289535038,0.9331927987311419339955,0.9452007083004420060395,0.9554345372414569605126,0.9640696808870741960397,0.9712834401839982005987,0.9772498680518207927997,0.9821355794371834432161,0.9860965524865013893852,0.9892758899783241946076,0.9918024640754038705556,0.9937903346742238648330,0.9953388119762812497490,0.9965330261969593315041,0.9974448696695720671985,0.9981341866996159620497,0.9986501019683699054733,0.9990323967867816431079,0.9993128620620841515449,0.9995165758576162227989,0.9996630707343231190606,0.9997673709209644749637,0.9998408914098424661203,0.9998922002665226116631,0.9999276519560748800260,0.9999519036559823972829,0.9999683287581668800787,0.9999793424930874532612,0.9999866542509840936616,0.9999914600945290081958,0.9999945874560922961402,0.9999966023268752699396,0.9999978875452974971502,0.9999986991925460827179,0.9999992066718480244054,0.9999995208167234096801,0.9999997133484281208061,0.9999998301732592852402,0.9999999003557368306652,0.9999999420986596003541,0.9999999666795515145714,0.9999999810104375341123,0.9999999892824097416891,0.9999999940096285989365,0.9999999966842540216738,0.9999999981824921369006,0.9999999990134123549623,0.9999999994696576737051,0.9999999997176841962957,0.9999999998511771778238,0.9999999999223115241829,0.9999999999598399941614,0.9999999999794421109060,0.9999999999895790230120,0.9999999999947690424559,0.9999999999973998730344,0.9999999999987201874561,0.9999999999993762155537,0.9999999999996989372019,0.9999999999998561161362,0.9999999999999319077511,0.9999999999999680910833,0.9999999999999851934625,0.9999999999999931966885,0.9999999999999969046412,0.9999999999999986054829,0.9999999999999993779039,0.9999999999999997252041,0.9999999999999998798065,0.9999999999999999479443,0.9999999999999999776761,0.9999999999999999905205,0.9999999999999999960142,0.9999999999999999983406,0.9999999999999999993159,0.9999999999999999997208,0.9999999999999999998871,0.9999999999999999999548,0.9999999999999999999821,0.9999999999999999999930,0.9999999999999999999973,0.9999999999999999999990,0.9999999999999999999996,0.9999999999999999999998,0.9999999999999999999999,1.000000000000000000000,1.000000000000000000000];
  
  // The algorithm should have an absolute error of less than 8e−16 on the points sampled, but
  // this also depends on the precision of the exponential function used.
  //
  // Here, the absolute error is computed to be less than 2e-15,
  // and the relative error is computed to be less than 8e-16.
  for (var i=0; i<x.length; ++i) {
	assert.ok(Math.abs( (PortfolioAllocation.normcdf_(x[i]) - p[i]) ) <= 2e-15, 'Normcdf comparison v.s. Wolfram Alpha, absolute error');
	assert.ok(Math.abs( (PortfolioAllocation.normcdf_(x[i]) - p[i]) ) <= 8e-16 * Math.abs(x[i]), 'Normcdf comparison v.s. Wolfram Alpha, relative error');
  } 
});


QUnit.test('Inverse of the standard normal cumulative distribution function computation', function(assert) {    
  // Boundaries
  assert.equal(PortfolioAllocation.norminv_(0), Number.NEGATIVE_INFINITY, 'Norminv -inf');  
  assert.equal(PortfolioAllocation.norminv_(1), Number.POSITIVE_INFINITY, 'Norminv +inf');  
  
  // Values taken from Algorithm AS 241: The Percentage Points of the Normal Distribution, Michael J. Wichura
  // Each of them test a different internal code branch.
  //
  // Ideally, absolute error here should be 0.
  var p = [0.25, 0.001, 1e-20];
  var z = [-0.6744897501960817, -3.090232306167814, -9.262340089798408];
  for (var i =0 ; i < p.length; ++i) {
	assert.ok(Math.abs( (PortfolioAllocation.norminv_(p[i]) - z[i]) ) <= 5e-14, 'Norminv comparison v.s. Wichura');
  }
  
  // Values taken from Wolfram Alpha, with InverseCDF[NormalDistribution[0, 1], i/100] command, 22 digits precision requested.
  var p = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99];
  var x = [-2.326347874040841100886,-2.053748910631823052937,-1.880793608151250938868,-1.750686071252169979435,-1.644853626951472714864,-1.554773594596853541090,-1.475791028179170735221,-1.405071560309632555951,-1.340755033690216379613,-1.281551565544600466965,-1.226528120036610080449,-1.174986792066090005868,-1.126391129038800589205,-1.080319340814956118530,-1.036433389493789579713,-0.9944578832097531677397,-0.9541652531461944091517,-0.9153650878428140497859,-0.8778962950512285953771,-0.8416212335729142051787,-0.8064212470182402084887,-0.7721932141886846986886,-0.7388468491852136293212,-0.7063025628400874558801,-0.6744897501960817432022,-0.6433454053929169647476,-0.6128129910166272255783,-0.5828415072712162186919,-0.5533847195556728193145,-0.5244005127080407840383,-0.4958503473474533265668,-0.4676987991145082144086,-0.4399131656732338077535,-0.4124631294414047958027,-0.3853204664075676238108,-0.3584587932511937384668,-0.3318533464368165782307,-0.3054807880993973393709,-0.2793190344474541653227,-0.2533471031357997987982,-0.2275449766411494098089,-0.2018934791418508509514,-0.1763741647808613218049,-0.1509692154967772588685,-0.1256613468550740342102,-0.1004337205114697931399,-0.07526986209982982978492,-0.05015358346473361602091,-0.02506890825871103576236,0,0.02506890825871103576236,0.05015358346473361602091,0.07526986209982982978492,0.1004337205114697931399,0.1256613468550740342102,0.1509692154967772588685,0.1763741647808613218049,0.2018934791418508509514,0.2275449766411494098089,0.2533471031357997987982,0.2793190344474541653227,0.3054807880993973393709,0.3318533464368165782307,0.3584587932511937384668,0.3853204664075676238108,0.4124631294414047958027,0.4399131656732338077535,0.4676987991145082144086,0.4958503473474533265668,0.5244005127080407840383,0.5533847195556728193145,0.5828415072712162186919,0.6128129910166272255783,0.6433454053929169647476,0.6744897501960817432022,0.7063025628400874558801,0.7388468491852136293212,0.7721932141886846986886,0.8064212470182402084887,0.8416212335729142051787,0.8778962950512285953771,0.9153650878428140497859,0.9541652531461944091517,0.9944578832097531677397,1.036433389493789579713,1.080319340814956118530,1.126391129038800589205,1.174986792066090005868,1.226528120036610080449,1.281551565544600466965,1.340755033690216379613,1.405071560309632555951,1.475791028179170735221,1.554773594596853541090,1.644853626951472714864,1.750686071252169979435,1.880793608151250938868,2.053748910631823052937,2.326347874040841100886];

  // Check the absolute and relative error of the algorithm on [0.01, 0.99] interval
  for (var i=0; i<p.length; ++i) {
	assert.ok(Math.abs( (PortfolioAllocation.norminv_(p[i]) - x[i]) ) <= 2e-15 * Math.abs(x[i]), 'Norminv comparison v.s. Wolfram Alpha, relative');
	assert.ok(Math.abs( (PortfolioAllocation.norminv_(p[i]) - x[i]) ) <= 2e-15, 'Norminv comparison v.s. Wolfram Alpha, absolute');
  }  
});


QUnit.test('Hypersphere random sampler point computation', function(assert) {    
  // Test with random data
  {
	  // Setup static parameters of the random test
	  var nbTests = 100;
	  var nbSubTests = 10;
	  var minDimension = 1;
	  var maxDimension = 100;

	  // Aim of these tests is to check that for a sample of size n:
	  // - An array of coordinates of size n is returned
	  // - The 2-norm of the R^n vector associated to these coordinates is equal to 1
	  for (var i = 0; i < nbTests; ++i) {
		 // Generate a random dimension size and the associated sampler
		 var dimension = Math.floor(Math.random()*(maxDimension - minDimension + 1) + minDimension);
		 var sampler = new PortfolioAllocation.hypersphereRandomSampler_(dimension);
		 
		 for (var j = 0; j < nbSubTests; ++j) {
			// Generate a random sample point
			var sampledPoint = sampler.sample();
		 
			// Check that the number of coordinates of the sampled point corresponds to the requested dimension
			var sampledPointDimension = sampledPoint.length;
			assert.equal(sampledPointDimension, dimension, "Hypersphere random sampler point computation, coordinates length - Test " + i + "," + j);
		 
			 // Check that the coordinates 2-norm of the sampled point is 1, near to machine precision
			 var sampledPointCoordinatesSumSq = 0;
			 for (var k = 0; k < sampledPoint.length; ++k) {
				sampledPointCoordinatesSumSq += sampledPoint[k] * sampledPoint[k];
			}
			assert.equal(Math.abs(Math.sqrt(sampledPointCoordinatesSumSq) - 1) <= 5e-15, true, "Hypersphere random sampler point computation, coordinates 2-norm equal to one - Test " + i + "," + j);
		}
	  }
  }

});
