
function randomCorrelationMatrix(n, opt) {
	return Matrix_.randomCorrelation(n, opt);
}
function nearestCorrelationMatrix(mat, opt) {
	function unitDiagonalMatricesFrobeniusProjector(mat) {
		if (!mat.isSymmetric(epsSymmetric)) {
			throw new Error('internal error: input matrix must be symmetric');
		}
		var p = mat.unitDiagonalize();
		return p;
	}
	function pdMatricesFrobeniusProjector(mat, delta) {
		if (!mat.isSymmetric(epsSymmetric)) {
			throw new Error('internal error: input matrix must be symmetric');
		}
		var jacobi = Matrix_.eig(mat, {epsSymmetric: epsSymmetric});
		var q = jacobi[0];
		var lambdas = jacobi[1];
		var lambdas_p = lambdas.elemMap(function(i,j,val) { return Math.max(delta, val); });
		var p = Matrix_.axty(1, Matrix_.elementwiseProduct(q, lambdas_p.transpose()), q);
		p = p.symmetrize();
		return p;
	}
	var mat = new Matrix_(mat);
	if (opt === undefined) {
		opt = {};
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var eps = opt.eps
	if (eps == undefined) {
		eps = 1e-6;
	}
	var maxIterations = opt.maxIter
	if (maxIterations == undefined) {
		maxIterations = 100;
	}
	var minEigenvalue = opt.minEigenvalue
	if (minEigenvalue == undefined) {
		minEigenvalue = 1e-12;
	}
	if (!mat.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	var n = mat.nbRows;
	var x_k = Matrix_.copy(mat);
	var y_k = Matrix_.copy(mat);
	var delta_s_k = Matrix_.zeros(n, n);

	var diff_x_rel = Number.POSITIVE_INFINITY;
	var diff_y_rel = Number.POSITIVE_INFINITY;
	var diff_yx_rel = Number.POSITIVE_INFINITY;
	var iter = 0;
	while (true) {
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		var r = Matrix_.xmy(y_k, delta_s_k);
		var x_kp = pdMatricesFrobeniusProjector(r, minEigenvalue);
		var delta_s_kp = Matrix_.xmy(x_kp, r);
		var y_kp = unitDiagonalMatricesFrobeniusProjector(x_kp);
		diff_x_rel = Matrix_.xmy(x_kp, x_k).vectorNorm('infinity') / x_k.vectorNorm('infinity');
		diff_y_rel = Matrix_.xmy(y_kp, y_k).vectorNorm('infinity') / y_k.vectorNorm('infinity');
		diff_yx_rel = Matrix_.xmy(y_kp, x_kp).vectorNorm('infinity') / y_kp.vectorNorm('infinity');
		x_k = x_kp;
		y_k = y_kp;
		delta_s_k = delta_s_kp;
		if (eps == 0.0) {
			if (y_k.isCorrelationMatrix()) {
				break;
			}
		}
		else {
			if (Math.max(diff_x_rel, diff_y_rel, diff_yx_rel) <= eps) {
				break;
			}
		}
	}
	if (eps == 0.0) {
		return y_k;
	}
	else {
		return x_k;
	}
}
function repairCorrelationMatrix(corrMat, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var method = opt.method;
	if (method === undefined) {
		method = "nearest-correlation-matrix";
	}
	var minEigenvalue = opt.minEigenvalue
	if (minEigenvalue == undefined) {
		minEigenvalue = 0;
	}
	var epsNcm = opt.epsNcm
	if (epsNcm == undefined) {
		epsNcm = 1e-6;
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var epsUnitDiagonal = opt.epsUnitDiagonal
	if (epsUnitDiagonal == undefined) {
		epsUnitDiagonal = 1e-12;
	}
	var shrinkageTarget = opt.shrinkageTarget;
	if (method != "spectral" && method != "linear-shrinkage" && method != "nearest-correlation-matrix") {
		throw new Error('unsupported correlation matrix repair method');
	}
	var corrMat = new Matrix_(corrMat);
	if (!corrMat.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	if (!corrMat.isUnitDiagonal(epsUnitDiagonal)) {
		throw new Error('input matrix must be unit diagonal');
	}
	if (minEigenvalue < 0 || minEigenvalue > 1) {
		throw new Error('lower bound on the smallest eigenvalue(s) of the correlation matrix must belong to interval [0,1]');
	}
	var n = corrMat.nbRows;
	var jacobi = Matrix_.eig(corrMat, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
	if (jacobi[1].data[n-1] >= minEigenvalue) {
		return corrMat;
	}
	var c;
	if (method == "spectral") {
		var s = jacobi[0];
		var lambda = jacobi[1];
		var lambdap_sqrt = lambda.elemMap(function(i,j,val) { return Math.sqrt(Math.max(0, val)); });
		var bp = Matrix_.elementwiseProduct(s, lambdap_sqrt.transpose());
		var t_sqrt = Matrix_.fill(n, 1, function(i,j) { return 1 / bp.vectorNorm('two', 'row', i); });
		var b = Matrix_.elementwiseProduct(bp, t_sqrt);
		c = Matrix_.axty(1, b, b);
		c = c.unitDiagonalize();
	}
	else if (method == "linear-shrinkage") {
		if (shrinkageTarget == undefined) {
			shrinkageTarget = Matrix_.identity(n);
		}
		else {
			shrinkageTarget = new Matrix_(shrinkageTarget);
			if (!shrinkageTarget.isSymmetric(epsSymmetric)) {
				throw new Error('shrinkage target matrix must be symmetric');
			}
			if (!shrinkageTarget.isUnitDiagonal(epsUnitDiagonal)) {
				throw new Error('shrinkage target matrix must be unit diagonal');
			}
			var jacobi = Matrix_.eig(shrinkageTarget, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
			if (jacobi[1].data[n-1] < minEigenvalue) {
				throw new Error('smallest eigenvalue of the shrinkage target matrix strictly lower than the desired lower bound on the smallest eigenvalue(s) of the correlation matrix');
			}
		}
		var m0 = corrMat;
		var m1 = shrinkageTarget;
		
		var alpha_star_interval = bisection_(function (alpha) { 
												var s_alpha = Matrix_.axpby(alpha, m1, 1-alpha, m0);
												var jacobi = Matrix_.eig(s_alpha, {epsSymmetric: epsSymmetric, sortedEigenvalues: true});
												var lambdas = jacobi[1];
												var lambda_min = lambdas.data[n-1];
												return lambda_min - minEigenvalue; 
											}, 
											0, 1, {outputInterval: true});
		var alpha_star = alpha_star_interval[1];
		c = Matrix_.axpby(alpha_star, m1, 1-alpha_star, m0);
		c = c.symmetrize().unitDiagonalize();
 	}
	else if (method == "nearest-correlation-matrix") {
		c = nearestCorrelationMatrix(corrMat, {maxIter: -1, eps: epsNcm, minEigenvalue: minEigenvalue});
	}
	else {
		throw new Error('internal error: unsupported repair method');
	}
	return c;
}
function perturbedCorrelationMatrix(corrMat, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "additive-noise";
	}
	var epsCorrMat = opt.epsCorrMat
	if (epsCorrMat == undefined) {
		epsCorrMat = 1e-12;
	}
	var epsNcm = opt.epsNcm
	if (epsNcm == undefined) {
		epsNcm = 1e-6;
	}
	var maxIterNcm = opt.maxIterNcm
	if (maxIterNcm == undefined) {
		maxIterNcm = 100;
	}
	
	var sigma = opt.sigma;
	if (sigma == undefined) {
		sigma = 0.05;
	}
	
	var noiseLevelMax = opt.noiseLevelMax;
	if (noiseLevelMax == undefined) {
		noiseLevelMax = 0.01;
	}
	var noiseSpaceDimension = opt.noiseSpaceDimension;
	if (noiseSpaceDimension == undefined) {
		noiseSpaceDimension = 3;
	}
	var method = opt.method;
	if (method != "spectral-noise" && method != "multiplicative-noise" && method != "additive-noise") {
			throw new Error('unsupported perturbation method');
	}
	var corrMat = new Matrix_(corrMat);
	if (!corrMat.isCorrelationMatrix(epsCorrMat)) {
		throw new Error('input matrix must be symmetric, unit diagonal and positive semidefinite');
	}
	var n = corrMat.nbRows;
	var c;
	if (method == "spectral-noise"){		
		var jacobi = Matrix_.eig(corrMat, {epsSymmetric: epsCorrMat, sortedEigenvalues: true});
		var s = jacobi[0];
		var lambdas = jacobi[1];
		var lambdas_p = lambdas.elemMap(function(i,j,val) { return val * Math.exp(sigma * normrnd_(0, 1)); });
		var cp = Matrix_.axty(1, Matrix_.elementwiseProduct(s, lambdas_p.transpose()), s);
		c = cp.elemMap(function(i,j,val) { return val / Math.sqrt(cp.data[(i-1) * cp.nbColumns + (i-1)] * cp.data[(j-1) * cp.nbColumns + (j-1)]); });
		c = c.symmetrize().unitDiagonalize();
	}
	else if (method == "multiplicative-noise") {
		c = Matrix_.fillSymmetric(n, function(i,j) { 
										if (i == j) { 
											return 1; 
										} 
										else { 
											return Math.max(-1, Math.min(1, corrMat.data[(i-1) * corrMat.nbColumns + (j-1)] * (1 + sigma*normrnd_(0,1))));
										}
									});
		if (!c.isCorrelationMatrix()) {
			c = nearestCorrelationMatrix(c, {maxIter: maxIterNcm, eps: epsNcm});
		}
	}
	else if (method == "additive-noise") {
		var hypersphereRandomSampler = new hypersphereRandomSampler_(noiseSpaceDimension, true);
		var us = new Array(n);
		for (var i = 0; i < n; ++i) {
			us[i] = new Matrix_(hypersphereRandomSampler.sample());
		}
		c = Matrix_.fillSymmetric(n, function(i,j) { 
										if (i == j) { 
											return 1; 
										} 
										else { 
											var corrMat_ij = corrMat.data[(i-1) * corrMat.nbColumns + (j-1)];
											var noise_ij = noiseLevelMax * Matrix_.vectorDotProduct(us[i-1], us[j-1]);
											
											return Math.max(-1, Math.min(1, corrMat_ij + noise_ij));
										}
									});
		if (!c.isCorrelationMatrix()) {
			c = nearestCorrelationMatrix(c, {maxIter: maxIterNcm, eps: epsNcm});
		}									
	}
	else {
		throw new Error('internal error: unsupported perturbation method');
	}
	return c;
}
function covarianceMatrixFromCorrelationMatrix(corrMat, diagonalVec, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.diagonalVectorType === undefined) {
		opt.diagonalVectorType = "variances";
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var epsUnitDiagonal = opt.epsUnitDiagonal
	if (epsUnitDiagonal == undefined) {
		epsUnitDiagonal = 1e-12;
	}
	var diagonalVectorType = opt.diagonalVectorType;
	if (diagonalVectorType != "variances" && diagonalVectorType != "standard-deviations") {
		throw new Error('unsupported diagonal vector type');
	}
	var corrMat = new Matrix_(corrMat);
	var diagonalVec = new Matrix_(diagonalVec);
	if (!corrMat.isSymmetric(epsSymmetric)) {
		throw new Error('input correlation matrix not symmetric');
	}
	if (corrMat.nbRows != diagonalVec.nbRows) {
		throw new Error('input correlation matrix dimensions not compatible with input variances vector dimension');
	}
	if ( Matrix_.xmy(corrMat.diagonal(), Matrix_.ones(corrMat.nbRows, 1)).vectorNorm('infinity') > epsUnitDiagonal ) {
		throw new Error('input correlation matrix not unit diagonal');
	}
	var n = corrMat.nbRows;
	
	var stddev;
	if (diagonalVectorType == "variances") {
		stddev = diagonalVec.elemMap(function(i,j,val) { return Math.sqrt(val); });
	}
	else if (diagonalVectorType == "standard-deviations") {
		stddev = diagonalVec;
	}
	else {
		throw new Error('internal error: unsupported diagonal vector type');
	}
	var c = Matrix_.elementwiseProduct(Matrix_.elementwiseProduct(corrMat, stddev), stddev.transpose());
	c = c.symmetrize();
	return c;
}
function covarianceMatrix(arr, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "covariance";
	}
	var covarianceMethod = opt.method;
	if (covarianceMethod != "covariance" &&
	    covarianceMethod != "sample-covariance" &&
		covarianceMethod != "linear-shrinkage") {
			throw new Error('unsupported covariance matrix computation method');
	}

	var shrinkageTarget = opt.shrinkageTarget;
	if ( covarianceMethod == "linear-shrinkage" &&
	    ["constant-variance-null-correlation", 
		 "constant-variance-correlation", 
		 "null-correlation", 
		 "constant-correlation"].indexOf(shrinkageTarget) == -1 ) {
		 throw new Error('unsupported covariance matrix shrinkage target');
	}
	var nbSeries = arr.length;
	var nbObservations = arr[0].length;
	var obj;
	if (covarianceMethod == "covariance" || covarianceMethod == "sample-covariance") {
		var covarianceFunction = covariance_;
		if (covarianceMethod == "sample-covariance") {
			covarianceFunction = sampleCovariance_;
		}
		obj = allocateMatrix_(nbSeries, nbSeries);
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < i; ++j) {
				obj.data[i * obj.nbColumns + j] = obj.data[j * obj.nbColumns + i];
			}
			for (var j = i; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = covarianceFunction(arr[i], arr[j]);
			}
		}
	}
	else if (covarianceMethod == "linear-shrinkage") {
		var means = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]) });
		var x = Matrix_.fill(nbSeries, nbObservations, function(i,j) { return arr[i-1][j-1] - means.data[i-1] });
		var covMat = Matrix_.axty(1/nbObservations, x, x).toCovarianceMatrix();
		var stdvarVec = covMat.getStandardDeviations();	
		var corrMat = covMat.getCorrelationMatrix();
		var prior;
		var rBar;
		if (shrinkageTarget == "constant-variance-null-correlation") {
			var mu = covMat.trace() / nbSeries;
			prior = Matrix_.fill(nbSeries, nbSeries, function(i,j) { if (i == j) { return mu; } else { return 0; }  });
		}
		else if (shrinkageTarget == "constant-variance-correlation") {
			var phi = covMat.trace() / nbSeries;
			var nu = 0;
			for (var i = 0; i < nbSeries - 1; ++i) {
				for (var j = i+1; j < nbSeries; ++j) {
					nu += covMat.data[i * covMat.nbColumns + j];
				}
			}
			nu *= 2/(nbSeries*(nbSeries-1)); 
			prior = Matrix_.fill(nbSeries, nbSeries, function(i,j) { if (i == j) { return phi; } else { return nu; }  });
		}	
		else if (shrinkageTarget == "null-correlation") {
			prior = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									 if (i == j) { 
										 return stdvarVec.data[i-1]*stdvarVec.data[i-1]; 
									 } 
									 else { 
										 return 0;
									 } 
								 });
		}
		else if (shrinkageTarget == "constant-correlation") {
			rBar = 0;
			for (var i = 0; i < nbSeries - 1; ++i) {
				for (var j = i+1; j < nbSeries; ++j) {
					rBar += corrMat.data[i * corrMat.nbColumns + j];
				}
			}
			rBar *= 2/(nbSeries*(nbSeries-1)); 
			prior = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									 if (i == j) { 
										 return stdvarVec.data[i-1]*stdvarVec.data[i-1]; 
									 } 
									 else { 
										return rBar*stdvarVec.data[i-1]*stdvarVec.data[j-1] ;
									} 
								 });
		}
		var piMat = Matrix_.fill(nbSeries, nbSeries, 
								 function(i,j) { 
									var pi_ij = 0;
									for (var k = 0; k < nbObservations; ++k) {
										 pi_ij += Math.pow((arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)], 2);
									}
									pi_ij /= nbObservations;
									return pi_ij;
								 });
		var pi = piMat.sum();
		var rho;
		if (shrinkageTarget == "constant-variance-null-correlation") {
			rho = 0;
		}
		else if (shrinkageTarget == "constant-variance-correlation") {
			rho = 0;
		}	
		else if (shrinkageTarget == "null-correlation") {
			rho = piMat.trace();
		}
		else if (shrinkageTarget == "constant-correlation") {
			rho = 0;
			for (var i = 1; i <= nbSeries; ++i) {
				for (var j = 1; j <= nbSeries; ++j) {
					if (i == j) {
						continue;
					}
					var theta_ii__ij = 0;
					for (var k = 0; k < nbObservations; ++k) {				
						 theta_ii__ij += ( Math.pow(arr[i-1][k] - means.data[i-1], 2) - covMat.data[(i-1) * covMat.nbColumns + (i-1)] )
										  *
										  ( (arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)] );
					}
					theta_ii__ij /= nbObservations;
					
					var theta_jj__ij = 0;
					for (var k = 0; k < nbObservations; ++k) {
						 theta_jj__ij += ( Math.pow(arr[j-1][k] - means.data[j-1], 2) - covMat.data[(j-1) * covMat.nbColumns + (j-1)] )
										  *
										  ( (arr[i-1][k] - means.data[i-1])*(arr[j-1][k] - means.data[j-1]) - covMat.data[(i-1) * covMat.nbColumns + (j-1)] );
					}
					theta_jj__ij /= nbObservations;
					rho += (stdvarVec.data[j-1]/stdvarVec.data[i-1]) * theta_ii__ij  
						   +
						   (stdvarVec.data[i-1]/stdvarVec.data[j-1]) * theta_jj__ij;
				}
			}
			rho *= rBar/2;
			rho += piMat.trace();
		}
		var gamma = Math.pow(Matrix_.xmy(covMat, prior).matrixNorm('frobenius'), 2);
		var kappa = (pi - rho)/gamma;
		var shrinkage = Math.max(0, Math.min(1, kappa/nbObservations));
		obj = Matrix_.axpby(shrinkage, prior, 1-shrinkage, covMat);
	}
	else {
		throw new Error('internal error: unsupported covariance matrix computation method');
	}
	addCovarianceMatrixMethods_(obj);
	return obj;
}
function addCovarianceMatrixMethods_(matrix) {
	var methods = {
		'getCorrelationMatrix': function(out) { 
			var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
			for (var i = 0; i < obj.nbRows; ++i) {
				var stdDevI = Math.sqrt(this.data[i * this.nbColumns + i]);
				for (var j = 0; j < i; ++j) {
					obj.data[i * obj.nbColumns + j] = obj.data[j * this.nbColumns + i];
				}
				for (var j = i; j < obj.nbColumns; ++j) {
					var stdDevJ = Math.sqrt(this.data[j * this.nbColumns + j]);
				
					obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j] / ( stdDevI * stdDevJ );
				}
			}
			return obj;
		},
		'getVariances': function(out) { 
			return this.diagonal(out);
		},
		'getStandardDeviations': function(out) { 
			return this.diagonal(out).elemMap(function(i,j,val) { return Math.sqrt(val);}, out);
		},
    	'standardizedGeneralizedVariance': function () {
        	var det = this.determinant();
    		if (det < 0) {
    		    throw new Error('covariance matrix is not positive semi-definite');
    		}
    		var sgv = Math.pow(det, 1/this.nbRows);
    		return sgv;
    	},
	};
  for (var name in methods) {
    matrix[name] = methods[name];
  }
};
function Matrix_(input) {
	function isArray(arr) {
		var o = Object.prototype.toString.call(arr);
		
		if (o === '[object Array]' ||
		    o === '[object Int8Array]' || o === '[object Uint8Array]' || o === '[object Uint8ClampedArray]' ||
			o === '[object Int16Array]' || o === '[object Uint16Array]' ||
			o === '[object Int32Array]' || o === '[object Uint32Array]' ||
			o === '[object Float32Array]' || o === '[object Float64Array]') {
			return true;
		}
		else {
			return false;
		}		
	}
	
	function fromDoubleArray(dblarr) {
		that = allocateMatrix_(dblarr.length, dblarr[0].length);
		for (var i = 0; i < that.nbRows; ++i) {
			if (!isArray(dblarr[i])) {
				throw new Error('unsupported input type');
			}
			if (dblarr[i].length !== that.nbColumns) {
				throw new Error('unsupported input type');
			}
			for (var j = 0; j < that.nbColumns; ++j) {
				that.data[i * that.nbColumns + j] = dblarr[i][j];
			}
		}
		return that;
	}
	
	function fromArray(arr) {
		that = allocateMatrix_(arr.length, 1);
		for (var i = 0; i < that.nbRows; ++i) {
			that.data[i * that.nbColumns] = arr[i];
		}
		return that;
	}
	
	function fromMatrix(mat) {
		that = allocateMatrix_(mat.nbRows, mat.nbColumns);
		var nbElements = mat.nbRows * mat.nbColumns;
		for (var k = 0; k < nbElements; ++k) {
			that.data[k] = mat.data[k];
		}
		return that;
	}
	if (!(this instanceof Matrix_)) {
      return new Matrix_(input);
    }
	var that = this;
	if (isArray(input) && isArray(input[0])) { // Standard matrix
		return fromDoubleArray(input);
	}
	else if (isArray(input)) { // Simplified constructor for a column matrix (i.e., a vector)
		return fromArray(input);
	}
	else if (input instanceof Matrix_) { // Equivalent to a "clone" operation on the Matrix
		return fromMatrix(input);
	}
	else {
		throw new Error('unsupported input type');
	}
}
function allocateMatrix_(n, m, out) {
	var obj;
	if (out === undefined) {
    	obj = Object.create(Matrix_.prototype);
    	obj.nbRows = n;
    	obj.nbColumns = m;
    	obj.data = typeof Float64Array === 'function' ? new Float64Array(obj.nbRows * obj.nbColumns) : new Array(obj.nbRows * obj.nbColumns);
	}
	else if (!(out instanceof Matrix_)){
		throw new Error('provided output must be a matrix');
	}
	else if (out.nbRows !== n || out.nbColumns !== m) {
		throw new Error('provided output matrix size does not match expected matrix size: ' + '(' + out.nbRows + ',' + out.nbColumns + 
		') v.s. ' + '(' + n + ',' + m + ')');
	}
	else {
		obj = out;
	}
	return obj;
}


Matrix_.prototype = {
    constructor: Matrix_,
	toString: function() {
		var maxLen = 0;
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j=0; j < this.nbColumns; ++j) {
				var len = String(this.data[i * this.nbColumns + j]).length;
				if (len > maxLen) {
					maxLen = len;
				}
			}
		}
		var strMat = [];
		for (var i = 0; i < this.nbRows; ++i) {
			strMat.push('[ ');
			
			for (var j=0; j < this.nbColumns; ++j) {
				var strVal = String(this.data[i * this.nbColumns + j]);
				strMat.push(new Array(maxLen - strVal.length + 1).join(' ') + strVal + ' ');
			}
			
			strMat.push(']\n');
		}
		return strMat.join('');
	},
	toRowArray: function(fct) {
		fct = fct || function(i, j, val) { return true; };
		var arr = new Array(this.nbRows);
		for (var i = 0; i < this.nbRows; ++i) {
			arr[i] = new Array(this.nbColumns); // potential size, to be restrained after
			
			var k =0;
			for (var j=0; j < this.nbColumns; ++j) {
				var val = this.data[i * this.nbColumns + j];
				if (fct(i+1, j+1, val)) {
					arr[i][k++] = val;
				}
			}			
			arr[i].length = k; // Restrain the inner array size to effectively matching elements
		}
		return arr;
	},
	toArray: function(fct) {
		fct = fct || function(i, j, val) { return true; };
		var arr = new Array(this.nbRows * this.nbColumns);
		var k =0;
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j=0; j < this.nbColumns; ++j) {
				var val = this.data[i * this.nbColumns + j];
				if (fct(i+1, j+1, val)) {
					arr[k++] = val;
				}
			}			
		}
		arr.length = k; // Restrain the array size to effectively matching elements
		return arr;
	},
	setValueAt: function(i, j, val) {
		if (i < 1 || j < 1 || i > this.nbRows || j > this.nbColumns) {
			throw new Error(
			'index out of bounds when setting matrix value, (' + i + ',' + j +
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		this.data[(i-1) * this.nbColumns + (j-1)] = val;
		return this;
	},
	setValue: function(i, j, val) {
		this.data[(i-1) * this.nbColumns + (j-1)] = val;
		return this;
	},
	getValueAt: function(i, j) {
		if (i < 1 || j < 1 || i > this.nbRows || j > this.nbColumns) {
			throw new Error(
			'index out of bounds when getting matrix value, (' + i + ',' + j +
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		return this.data[(i-1) * this.nbColumns + (j-1)];
	},
	getValue: function(i, j) {
		return this.data[(i-1) * this.nbColumns + (j-1)];
	},	
	isSquare: function() {
	    return this.nbRows === this.nbColumns;
	},
	isSymmetric: function(eps) {
		if (eps == undefined) {
			eps = 0;
		}	
		if (!this.isSquare()) {
			return false;
		}
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (Math.abs(this.data[i * this.nbColumns + j] - this.data[j * this.nbColumns + i]) > eps) {
					return false;
				}
			}
		}
	    return true;
	},
	isUnitDiagonal: function(eps) {
		if (eps == undefined) {
			eps = 0;
		}	
		if (!this.isSquare()) {
			return false;
		}
		var max_a_ii = this.data[0];
		for (var i = 0; i < this.nbRows; ++i) {
			var a_ii = this.data[i * this.nbColumns + i];
			if (a_ii > max_a_ii) {
				max_a_ii = a_ii;
			}
		}
		if (Math.abs(max_a_ii - 1) > eps) {
			return false;
		}
	    return true;
	},
	isCorrelationMatrix: function(eps) {
		if (eps == undefined) {
			eps = 0;
		}	
		if (!this.isSquare()) {
			return false;
		}
		if (!this.isSymmetric(eps)) {
			return false;
		}
		if (!this.isUnitDiagonal(eps)) {
			return false;
		}
		try {
			var g = Matrix_.choleskyDecomposition(this);
		}
		catch (e) {
			if (e.message === "input matrix must be positive definite") {
				return false;
			}
			else {
				throw(e);
			}
		}
	    return true;
	},
	isVector: function() {
	    return this.nbColumns === 1;
	},
	isNonNegative: function() {
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] < 0) {
					return false;
				}
			}
		}
		return true;
	},
	isPositive: function() {
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] <= 0) {
					return false;
				}
			}
		}
		return true;
	},
	isNonPositive: function() {	
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] > 0) {
					return false;
				}
			}
		}
		return true;
	},
	isNegative: function() {
		for (var i = 0; i < this.nbRows; ++i) {
			for (var j = 0; j < this.nbColumns; ++j) {
				if (this.data[i * this.nbColumns + j] >= 0) {
					return false;
				}
			}
		}
		return true;
	},
	sum: function() {
		var sum = 0;
		for (var k = 0; k < this.data.length; ++k) {
			sum += this.data[k];
		}
		return sum;
	},
	min: function() {
		var minVal = Number.POSITIVE_INFINITY;
		for (var k = 0; k < this.data.length; ++k) {
		  if (this.data[k] < minVal) {
			minVal = this.data[k];
		  }
		}
		return minVal;
	},
	max: function() {
		var maxVal = Number.NEGATIVE_INFINITY;
		for (var k = 0; k < this.data.length; ++k) {
		  if (this.data[k] > maxVal) {
			maxVal = this.data[k];
		  }
		}
		return maxVal;
	},
	normalize: function(out) {
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		var sum = this.sum();
		if (sum === 0.0) {
			throw new Error('sum of coefficients of matrix is null');
		}
		for (var k = 0; k < this.data.length; ++k) {
			obj.data[k] = this.data[k] / sum; 
		}
		return obj;
	},
    diagonal: function (out) {
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		var obj = allocateMatrix_(this.nbRows, 1, out);
    	for (var i = 0; i < this.nbRows; ++i) {
    		obj.data[i] = this.data[i * (this.nbColumns + 1)] ;
    	}
        return obj;
    },
    row: function (i, out) {
		if (i < 1 || i > this.nbRows) {
			throw new Error(
			'index out of bounds when getting matrix row, (' + i + 
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		var obj = allocateMatrix_(this.nbColumns, 1, out);
    	for (var j = 0; j < this.nbColumns; ++j) {
    		obj.data[j] = this.data[(i-1) * this.nbColumns + j] ;
    	}
        return obj;
    },
    column: function (j, out) {
		if (j < 1 || j > this.nbColumns) {
			throw new Error(
			'index out of bounds when getting matrix column, (' + j + 
			') in size (' + this.nbRows + ',' + this.nbColumns + ')');
		}
		var obj = allocateMatrix_(this.nbRows, 1, out);
    	for (var i = 0; i < this.nbRows; ++i) {
    		obj.data[i] = this.data[i * this.nbColumns + (j-1)] ;
    	}
        return obj;
    },
	elemMap: function (fct, out) {
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = fct(i+1, j+1, this.data[i * this.nbColumns + j]);
			}
		}
		return obj;
	},
    submatrix : function(rindexes, cindexes, out) {
    	if (!(rindexes instanceof Array || rindexes instanceof Uint32Array) || rindexes.length == 0) {
    		throw new Error('first parameter must be a non empty array');
    	}
    	if (!(cindexes instanceof Array || cindexes instanceof Uint32Array) || cindexes.length == 0) {
    		throw new Error('second parameter must be a non empty array');
    	}
    	for (var i = 1; i < rindexes.length; ++i) {
    	    if (rindexes[i-1] >= rindexes[i]) {
    	        throw new Error('first parameter must be a sorted array');
    	    }
    	    
    	}
        for (var j = 1; j < cindexes.length; ++j) {
    	    if (cindexes[j-1] >= rindexes[j]) {
    	        throw new Error('second parameter must be a sorted array');
    	    }
    	    
    	}
		var obj = allocateMatrix_(rindexes.length, cindexes.length, out);
    	for (var i = 0; i < obj.nbRows; ++i) {
    		var rindex = rindexes[i] - 1; // row index of the original matrix
    		
    		for (var j = 0; j < obj.nbColumns; ++j) {
    		    var cindex = cindexes[j] - 1; // column index of the original matrix
    			obj.data[i * obj.nbColumns + j] = this.data[rindex * this.nbColumns + cindex];
    		}
    	}
        return obj;
    },
	transpose: function (out) {
		var obj = allocateMatrix_(this.nbColumns, this.nbRows, out);
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[j * this.nbColumns + i];
			}
		}
		return obj;
	},
	determinant: function () {
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		var R = Matrix_.qrDecomposition(this, {qLess: true});
		var det = 1.0;
		for (var k = 0; k < R.data.length; k += R.nbColumns + 1) {
		        det *= R.data[k];
		}
		return det;
	},
	trace: function () {
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		var tr = 0;
		for (var i = 0; i < this.nbRows; ++i) {
			tr += this.data[i * this.nbColumns + i];
		}
		return tr;
	},
    symmetrize: function (out) {
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		for (var i = 1; i <= obj.nbRows; ++i) {
			obj.data[(i-1) * obj.nbColumns + (i-1)] = this.data[(i-1) * this.nbColumns + (i-1)];
			
			for (var j = i+1; j <= obj.nbColumns; ++j) {
				var a_ij = this.data[(i-1) * this.nbColumns + (j-1)];
				var a_ji = this.data[(j-1) * this.nbColumns + (i-1)];
				var a_ij_fixed = (a_ij + a_ji) / 2;
				obj.data[(i-1) * obj.nbColumns + (j-1)] = a_ij_fixed;
				obj.data[(j-1) * obj.nbColumns + (i-1)] = a_ij_fixed;
			}
		}
        return obj;
    },
    unitDiagonalize: function (out) {
    	if (!this.isSquare()) {
    		throw new Error('matrix is not square: ' + '(' + this.nbRows + ',' + this.nbColumns + ')');
    	}
		var obj = allocateMatrix_(this.nbRows, this.nbColumns, out);
		for (var i = 0; i < obj.nbRows; ++i) {
			for (var j = 0; j < obj.nbColumns; ++j) {
				obj.data[i * obj.nbColumns + j] = this.data[i * this.nbColumns + j];
			}
			
			obj.data[i * obj.nbColumns + i] = 1;
		}
        return obj;
    },
	matrixNorm: function(p) {
		if (p == 'one') {
			var maxColSum = 0;
            for (var j = 0; j < this.nbColumns; ++j) {
                var colSum = 0;
			    for (var i = 0; i < this.nbRows; ++i) {
				    colSum += Math.abs(this.data[i * this.nbColumns + j]);
			     }
			     maxColSum = Math.max(maxColSum, colSum);
            }
            return maxColSum;
		}
		else if (p == 'infinity') {
		    var maxRowSum = 0;
		    for (var i = 0; i < this.nbRows; ++i) {
			    var rowSum = 0;
			    for (var j = 0; j < this.nbColumns; ++j) {
				    rowSum += Math.abs(this.data[i * this.nbColumns + j]);
			    }
			    maxRowSum = Math.max(maxRowSum, rowSum);
		    }
		    return maxRowSum;
		}
		else if (p == 'frobenius') {
			return this.vectorNorm('two');
		}
		else {
			throw new Error('unsupported matrix norm: ' + p);
		}
	},
	vectorNorm: function(p, subset, idx) {
		var idxStartRow;
		var idxEndRow;
		var idxStartColumn;
		var idxEndColumn;
		if (subset === undefined || subset == 'matrix') {
			idxStartRow = 0;
			idxEndRow = this.nbRows;
			idxStartColumn = 0;
			idxEndColumn = this.nbColumns;
		}
		else if (subset == 'row') {
			if (idx !== undefined) {
				if (idx < 1 || idx > this.nbRows) { 
					throw new Error('row index out of bounds: ' + idx);
				}
				idxStartRow = idx-1;
				idxEndRow = idx;
				idxStartColumn = 0;
				idxEndColumn = this.nbColumns;
			}
			else {
				throw new Error('undefined row index');
			}
		}
		else if (subset == 'column') {
			if (idx !== undefined) {
				if (idx < 1 || idx > this.nbColumns) { 
					throw new Error('column index out of bounds: ' + idx);
				}
				idxStartRow = 0;
				idxEndRow = this.nbRows;
				idxStartColumn = idx-1;
				idxEndColumn = idx;
			}
			else {
				throw new Error('undefined row index');
			}
		}
		else if (subset !== undefined) {
			throw new Error('unsupported matrix subset: ' + subset);
		}
		if (p == 'one') {
			var absSum = 0;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					absSum += Math.abs(this.data[ij_idx]);
					ij_idx++;
				}
				i_m += this.nbColumns;
			}
			return absSum;
		}
		else if (p == 'two') {
			var t = 0;
			var s = 1;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
			    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
				    var val = this.data[ij_idx];
					var absVal = Math.abs(val);
					if (absVal != 0) {
						if (absVal > t) {
							s = 1 + s * (t/val) * (t/val);
							t = absVal;
						}
						else  {
							s = s + (val/t) * (val/t);
						}
					}
					ij_idx++;
				}
				i_m += this.nbColumns;
			}
			return t * Math.sqrt(s);
		}
		else if (p == 'infinity') {
		    var maxAbsVal = 0;
			var i_m = idxStartRow * this.nbColumns; // i * this.nbColumns
			for (var i = idxStartRow; i < idxEndRow; ++i) {
    		    var ij_idx = i_m + idxStartColumn;
				for (var j = idxStartColumn; j < idxEndColumn; ++j) {
					 maxAbsVal = Math.max(maxAbsVal, Math.abs(this.data[ij_idx]));
					 ij_idx++;
				}
				i_m += this.nbColumns;
			}
		    return maxAbsVal;
		}
		else {
			throw new Error('unsupported vector norm: ' + p);
		}
	},
	toCovarianceMatrix: function() {
		addCovarianceMatrixMethods_(this);
		return this;
	},
	
	
};
Matrix_.areEqual = function (a, b, eps) {
	if (!(a instanceof Matrix_)) {
		throw new Error('a must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('b must be a matrix');
	}
	if (a.nbRows !== b.nbRows) {
		return false;
	}
	if (a.nbColumns !== b.nbColumns) {
		return false;
	}
	var nbRows = a.nbRows;
	var nbColumns = a.nbColumns;
	var nbElements = nbRows * nbColumns;
	var tol = eps || 0;
	for (var k = 0; k < nbElements; ++k) {
		if (Math.abs(a.data[k] - b.data[k]) > tol) {
			return false;
		}
	}
	return true;
};
Matrix_.xpy = function(X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = X.data[k] + Y.data[k];
	}
    return obj;
};
Matrix_.xmy = function(X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = X.data[k] - Y.data[k];
	}
    return obj;
};
Matrix_.axpby = function(a, X, b, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (X.nbColumns !== Y.nbColumns || X.nbRows !== Y.nbRows) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	var nbElements = X.nbRows * X.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = a*X.data[k] + b*Y.data[k];
	}
    return obj;
};
Matrix_.copy = function(A, out) {
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	var obj = allocateMatrix_(A.nbRows, A.nbColumns, out);
	var nbElements = A.nbRows * A.nbColumns;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = A.data[k];
	}
	return obj;
};
Matrix_.elementwiseProduct = function(X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(X.nbRows === Y.nbRows && X.nbColumns === Y.nbColumns || 
		  Y.nbRows === X.nbRows && Y.nbColumns === 1 ||
		  Y.nbRows === 1 && Y.nbColumns === X.nbColumns)) {
		throw new Error('input matrices sizes do not match: ' + '(' + X.nbRows + ',' + X.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	var nbRows = X.nbRows;
	var nbColumns = X.nbColumns;
	if (X.nbRows === Y.nbRows && X.nbColumns === Y.nbColumns) { // full matrix elementwise product	
		var nbElements = nbRows * nbColumns;
		for (var k = 0; k < nbElements; ++k) {
			obj.data[k] = X.data[k] * Y.data[k];
		}
	}
	else if (Y.nbRows === X.nbRows && Y.nbColumns === 1) { // row matrix elementwise product
		var x_ij_idx = 0;
		for (var i = 0; i < nbRows; ++i) {
			var y_i = Y.data[i];
			for (var j = 0; j < nbColumns; ++j) {
				obj.data[x_ij_idx] = X.data[x_ij_idx] * y_i;
				x_ij_idx++;
			}
		}
	}
	else if (Y.nbRows === 1 && Y.nbColumns === X.nbColumns) { // column matrix elementwise product
		var x_ij_idx = 0;
		for (var i = 0; i < nbRows; ++i) {
			for (var j = 0; j < nbColumns; ++j) {
				obj.data[x_ij_idx] = X.data[x_ij_idx] * Y.data[j];
				x_ij_idx++
			}
		}
	}
    return obj;
};
Matrix_.xy = function(X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (X.nbColumns !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, Y.nbColumns, out);
	var n = X.nbRows;
	var m = X.nbColumns;
	var p = Y.nbColumns;
	
	var x_ik_idx = 0;
	var i_p = 0; // == i*p
	for (var i = 0; i < n; ++i) {
		var obj_ij_idx = i_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = 0; // obj(i,j) = 0
			obj_ij_idx++;
		}
	
		var y_kj_idx = 0;
		for (var k = 0; k < m; ++k) {
			var x_ik = X.data[x_ik_idx]; // x(i,k)
			
			obj_ij_idx = i_p;
			for (var j = 0; j < p; ++j) {
				obj.data[obj_ij_idx] += x_ik * Y.data[y_kj_idx]; // obj(i,j) += x(i,k) * y(k,j);
				y_kj_idx++;
				obj_ij_idx++;
			}
			x_ik_idx++;
		}
		i_p += p;
	}
    return obj;
};
Matrix_.txy = function(X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	if (X.nbRows !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbColumns, Y.nbColumns, out);
	var n = X.nbColumns;
	var m = X.nbRows;
	var p = Y.nbColumns;
	
	var x_ki_idx = 0;
	var obj_ij_idx = 0;
	var k = 0;
	for (var i = 0; i < n; ++i) {
		var x_ki = X.data[x_ki_idx]; // x(k,i)
		var y_kj_idx = k_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = x_ki * Y.data[j]; // obj(i,j) = x(k,i) *y(k,j)
			obj_ij_idx++;
		}
		x_ki_idx++;
	}
	
	var k_p = p; // == k*p
    for (var k = 1; k < m; ++k) {
        obj_ij_idx = 0;
		for (var i = 0; i < n; ++i) {
            var x_ki = X.data[x_ki_idx];  // x(k,i)
			var y_kj_idx = k_p;
			for (var j = 0; j < p; ++j) {
                obj.data[obj_ij_idx] += x_ki * Y.data[y_kj_idx]; // obj(i,j) = x(k,i) *y(k,j)
				obj_ij_idx++;
				y_kj_idx++;
            }
			x_ki_idx++;
        }
		k_p += p;
    }
    return obj;
};
Matrix_.axy = function(a, X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	if (X.nbColumns !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, Y.nbColumns, out);
	var n = X.nbRows;
	var m = X.nbColumns;
	var p = Y.nbColumns;
	
	var x_ik_idx = 0;
	var i_p = 0; // == i*p
	for (var i = 0; i < n; ++i) {
		var obj_ij_idx = i_p;
		for (var j = 0; j < p; ++j) {
			obj.data[obj_ij_idx] = 0; // obj(i,j) = 0
			obj_ij_idx++;
		}
	
		var y_kj_idx = 0;
		for (var k = 0; k < m; ++k) {
			var x_ik = a * X.data[x_ik_idx]; // a * x(i,k)
			
			obj_ij_idx = i_p;
			for (var j = 0; j < p; ++j) {
				obj.data[obj_ij_idx] += x_ik * Y.data[y_kj_idx]; // obj(i,j) += x(i,k) * y(k,j);
				y_kj_idx++;
				obj_ij_idx++;
			}
			x_ik_idx++;
		}
		i_p += p;
	}
    return obj;
};
Matrix_.ax = function(a, X, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	var obj = allocateMatrix_(X.nbRows, X.nbColumns, out);
	var n = X.nbRows;
	var m = X.nbColumns;
	
	for (var i = 1; i <= n; ++i) {
		for (var j = 1; j <= m; ++j) {
			obj.setValue(i, j,
			             a * X.getValue(i, j));
		}
	}
    return obj;
};
Matrix_.axty = function(a, X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	if (X.nbColumns !== Y.nbColumns) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbRows, Y.nbRows, out);
    for (var i = 0; i < X.nbRows; ++i) {
        for (var j = 0; j < Y.nbRows; ++j) {
            obj.data[i * obj.nbColumns + j] = 0;
            
            for (var k = 0; k < X.nbColumns; ++k) {
                obj.data[i * obj.nbColumns + j] += X.data[i * X.nbColumns + k] * Y.data[j * Y.nbColumns + k];
            }
            
            obj.data[i * obj.nbColumns + j] = a * obj.data[i * obj.nbColumns + j];
        }
    }
    return obj;
};
Matrix_.atxy = function(a, X, Y, out) {
	if (!(X instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	if (!(Y instanceof Matrix_)) {
		throw new Error('third input must be a matrix');
	}
	if (X.nbRows !== Y.nbRows) {
		throw new Error('matrices sizes do not match: ' + '(' + X.nbRows + ',' + Y.nbColumns + 
		') - ' + '(' + Y.nbRows + ',' + Y.nbColumns + ')');
	}
	var obj = allocateMatrix_(X.nbColumns, Y.nbColumns, out);
	var k = 0;
	for (var i = 0; i < X.nbColumns; ++i) {
		var a_ik = a * X.data[k * X.nbColumns + i];
		
		for (var j = 0; j < Y.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] = 0;
		}

		for (var j = 0; j < Y.nbColumns; ++j) {
			obj.data[i * obj.nbColumns + j] += a_ik * Y.data[k * Y.nbColumns + j];
		} 
	}
    for (var k = 1; k < X.nbRows; ++k) {
        for (var i = 0; i < X.nbColumns; ++i) {
            var a_ik = a * X.data[k * X.nbColumns + i];
			
            for (var j = 0; j < Y.nbColumns; ++j) {
                obj.data[i * obj.nbColumns + j] += a_ik * Y.data[k * Y.nbColumns + j];
            } 
        }
    }
    return obj;
};
Matrix_.diagonal = function(x, out) {
	if (!(x instanceof Matrix_) || !x.isVector()) {
		throw new Error('first input must be a vector');
	}
	var n = x.nbRows;
	var obj = allocateMatrix_(n, n, out);
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < n; ++j) {
			obj.data[i * n + j] = 0;
		}
		obj.data[i * n + i] = x.data[i];
	}
    return obj;
};
Matrix_.fillSymmetric = function(n, fct, out) {
	if (n < 1) {
		throw new Error('input number of rows and columns out of bounds: ' + n);
	}
	var obj = allocateMatrix_(n, n, out);
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < i; ++j) {
			obj.data[i * n + j] = obj.data[j * n + i];
		}
		for (var j = i; j < obj.nbColumns; ++j) {
			obj.data[i * n + j] = fct(i + 1, j + 1);
		}	
	}
    return obj;
};
Matrix_.fill = function(n, m, fct, out) {
	if (n < 1) {
		throw new Error('input number of rows out of bounds: ' + n);
	}
	if (m < 1) {
		throw new Error('input number of columns out of bounds: ' + m);
	}
	var obj = allocateMatrix_(n, m, out);
	for (var i = 0; i < n; ++i) {
		for (var j = 0; j < m; ++j) {
			obj.data[i * m + j] = fct(i + 1 , j + 1);
		}
	}
    return obj;
};
Matrix_.zeros = function(n, m, out) {
	var obj = allocateMatrix_(n, m, out);
	var nbElements = n * m;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 0;
	}
    return obj;	

}
Matrix_.ones = function(n, m) {
	var obj = allocateMatrix_(n, m);
	var nbElements = n * m;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 1;
	}
    return obj;
}
Matrix_.identity = function(n) {
	var obj = allocateMatrix_(n, n);
	var nbElements = n * n;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = 0;
		if (k % (n + 1) == 0) {
			obj.data[k] = 1;
		}
	}
    return obj;
}
Matrix_.normrnd = function(n, m, mu, sigma) {
	var obj = allocateMatrix_(n, m);
	var nbElements = n * m;
	for (var k = 0; k < nbElements; ++k) {
		obj.data[k] = normrnd_(mu, sigma);
	}
    return obj;
}
Matrix_.vectorHadamardProduct = function(x, y) {
	return Matrix_.elementwiseProduct(x, y);
}
Matrix_.vectorDotProduct = function(x, y) {
	if (!x.isVector()) {
		throw new Error('first argument must be a vector');
	}
	if (!y.isVector()) {
		throw new Error('second argument must be a vector');
	}
	if (x.nbRows !== y.nbRows) {
		throw new Error('Vectors sizes do not match: ' + '(' + x.nbRows + ') - ' + '(' + y.nbRows + ')');
	}
	var dotProd = 0;
	var nbElements = x.nbRows;
	for (var i = 0; i < nbElements; ++i) {
		dotProd += x.data[i] * y.data[i]; 
	}
	return dotProd;
}
Matrix_.choleskyDecomposition = function(A, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	if (!(A instanceof Matrix_)) {
		throw new Error('input must be a matrix');
	}
	if (!A.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	var n = A.nbRows;
	var G = Matrix_.zeros(n,n);
	var v = Matrix_.zeros(n, 1);
	for (var j = 1; j <= n; ++j) {
		for (var i = j; i <= n; ++i) {
			v.data[i-1] = A.data[(i-1) * A.nbColumns + (j-1)];
		}
		for (var k = 1; k <= j-1; ++k) {
			for (var i = j; i <= n; ++i) {
				v.data[i-1] -= G.data[(j-1) * G.nbColumns + (k-1)] * G.data[(i-1) * G.nbColumns + (k-1)]
			}
		}
		if (v.data[j-1] < 0) {
			throw new Error('input matrix must be positive definite');
		}
		for (var i = j; i <= n; ++i) {
			G.data[(i-1) * G.nbColumns + (j-1)] = v.data[i-1] / Math.sqrt(v.data[j-1]);
		}		
	}
	return G;
}
Matrix_.qrDecomposition = function(A, opt) {
	function givens(a, b) {
		var c;
		var s;
		var r;
		
		if (b == 0) {
			c = (a >= 0) ? 1 : -1; // emulates sign(a)
			s = 0;
			r = Math.abs(a);
		}
		else if (a == 0) {
			c = 0;
			s = ((b >= 0) ? 1 : -1); // emulates sign(b)
			r = Math.abs(b);
		}
		else if (Math.abs(a) > Math.abs(b)) {
			var t = b/a;
			var u = ((a >= 0) ? 1 : -1) * Math.sqrt(1 + t*t);
			c = 1/u;
			s = t*c;
			r = a*u;
		}
		else {
			var t = a/b;
			var u = ((b >= 0) ? 1 : -1) * Math.sqrt(1 + t*t);
			s = 1/u;
			c = t*s;
			r = b*u;
		}
		
		return [c, -s, r]; // Compared to the second reference, the sign of s is altered so that it is made compatible with the notations of the first reference.
	}
	if (opt === undefined) {
		opt = {};
	}
	var isQLessDecomposition = opt.qLess || false;
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (A.nbRows < A.nbColumns) {
		throw new Error('matrix has more columns than rows: ' + '(' + A.nbRows + ') v.s. ' + '(' + A.nbColumns + ')');
	}
	var m = A.nbRows;
	var n = A.nbColumns;
	var R = new Matrix_(A); // represents R
	var Q = null;
	if (!isQLessDecomposition) {
		Q = Matrix_.identity(m); // represents Q
	}
	for (var j = 1; j <= n; ++j) {
		for (var i = m; i >= j+1; --i) {
			var R_im_j_idx = (i-2) * R.nbColumns + (j-1);
			var R_i_j_idx = (i-1) * R.nbColumns + (j-1);
			
			var a = R.data[R_im_j_idx]; // R(i-1, j)
			var b = R.data[R_i_j_idx]; // R(i, j)
			var givensArr = givens(a, b);
			var c = givensArr[0];
			var s = givensArr[1];
			var r = givensArr[2];
			R.data[R_im_j_idx] = r; // R(i-1, j) = r
			R.data[R_i_j_idx] = 0; // R(i, j) = 0
			for (var k = j+1; k <= n; ++k) {
				var R_im_k_idx = (i-2) * R.nbColumns + (k-1);
				var R_i_k_idx = (i-1) * R.nbColumns + (k-1);
			
				var t1 = R.data[R_im_k_idx]; // t1 = R(i-1, k)
				var t2 = R.data[R_i_k_idx]; // t2 = R(i, k)
				R.data[R_im_k_idx] = c*t1 - s*t2; // R(i-1, k) = ...
				R.data[R_i_k_idx] = s*t1 + c*t2; // R(i, k) = ...
			}
			if (!isQLessDecomposition) {
				for (var k = 1; k <= m; ++k) {
					var Q_k_im_idx = (k-1) * Q.nbColumns + (i-2);
					var Q_k_i_idx = (k-1) * Q.nbColumns + (i-1);
				
					var t1 = Q.data[Q_k_im_idx] // t1 = Q(k,i-1)
					var t2 = Q.data[Q_k_i_idx] // t2 = Q(k,i)
					Q.data[Q_k_im_idx] = c*t1 - s*t2; // Q(k,i-1) = ...
					Q.data[Q_k_i_idx] = s*t1 + c*t2; // Q(k,i) = ...
				}
			}
		}
	}
	if (!isQLessDecomposition) {
		return [Q, R];
	}
	else {
		return R;
	}
}
Matrix_.svdDecomposition = function(A, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-16;
	var maxIterations = opt.maxIter || 100;
	var svdForm = opt.svdForm || 'thin';
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (A.nbRows < A.nbColumns) {
		throw new Error('matrix has more columns than rows: ' + A.nbColumns + ' v.s. ' + A.nbRows);
	}
	var m = A.nbRows;
	var n = A.nbColumns;
	var uu = new Matrix_(A); // represents U 
	var u_frob_norm = uu.matrixNorm('frobenius');
	var vv = Matrix_.identity(n); // represents V
	var iter = 0;
	var u_columns_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	while (true) {
		++iter;
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
    	for (var j = 1; j <= n; ++j) {
    		var u_j_two_norm = uu.vectorNorm('two', 'column', j);
			u_columns_two_norm_sq[j-1] = u_j_two_norm * u_j_two_norm;
    	}
		var converged = true; // the convergence status is updated during the sweep
		for (var j = 2; j <= n; ++j) {
			for (var i = 1; i <= j-1; ++i) {
				var aa = u_columns_two_norm_sq[i-1]; // aa = sum U(k,i)^2, k=1..m = ||U(:,i)||_2^2
				var bb = u_columns_two_norm_sq[j-1]; // bb = sum U(k,j)^2, k=1..m = ||U(:,j)||_2^2
				var cc = 0; // cc = sum U(k,i)*U(k,j), k=1..m = <U(:,i)/U(:.j)>
			    for (var k = 1; k <= m; ++k) {
					cc += uu.data[(k-1) * uu.nbColumns + (i-1)] * uu.data[(k-1) * uu.nbColumns + (j-1)]; 
				}
				if (Math.abs(cc) <= eps * Math.sqrt(m * aa * bb)) { // this condition also covers the case cc = 0, which would make zeta below indefinite
					continue;
				}
				if (Math.sqrt(aa) <= eps * u_frob_norm || Math.sqrt(bb) <= eps * u_frob_norm) {
				    continue;
				}
				converged = false;
				var zeta = (aa - bb)/(2 * cc);
				var t = ((zeta >= 0) ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta)); // first part emulates sign(zeta)
				var cs = 1 / Math.sqrt(1 + t*t);
				var sn = cs * t;
				for (var k = 1; k <= m; ++k) {
					var t1 = uu.data[(k-1) * uu.nbColumns + (i-1)] // t1 = U(k,i)
					var t2 = uu.data[(k-1) * uu.nbColumns + (j-1)] // t2 = U(k,j)
					uu.data[(k-1) * uu.nbColumns + (i-1)] = cs*t1 + sn*t2 // U(k,i) = ...
					uu.data[(k-1) * uu.nbColumns + (j-1)] = -sn*t1 + cs*t2 // U(k,j) = ...
				}
                u_columns_two_norm_sq[i-1] += t * cc;
                u_columns_two_norm_sq[j-1] -= t * cc;
				for (var k = 1; k <= n; ++k) {
					var t1 = vv.data[(k-1) * vv.nbColumns + (i-1)] // t1 = V(k,i)
					var t2 = vv.data[(k-1) * vv.nbColumns + (j-1)] // t2 = V(k,j)
					vv.data[(k-1) * vv.nbColumns + (i-1)] = cs*t1 + sn*t2 // V(k,i) = ...
					vv.data[(k-1) * vv.nbColumns + (j-1)] = -sn*t1 + cs*t2 // V(k,j) = ...
                }
            }
		}
		if (converged == true) {
			break;
		}
	}
	var sigmas = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	var sigmas_idx = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
	for (var j = 1; j <= n; ++j) {
		var singVal_j = Math.sqrt(u_columns_two_norm_sq[j-1]); // == uu.vectorNorm('two', 'column', j);
		sigmas[j-1] = singVal_j; 
		sigmas_idx[j-1] = j;
	}
	sigmas_idx.sort(function(a, b) { return sigmas[b-1] - sigmas[a-1]; });
	var uuu = Matrix_.zeros(m, n); 
	var sss = Matrix_.zeros(n, n);
	var vvv = Matrix_.zeros(n, n);
	for (var j = 1; j <= n; ++j) {
		var sigma_j_col_idx = sigmas_idx[j-1];
		var sigma_j = sigmas[sigma_j_col_idx-1];
		sss.data[(j-1) * sss.nbColumns + (j-1)] = sigma_j;
		if (sigma_j != 0) { 
			for (var i = 1; i <= m; ++i) {
				uuu.data[(i-1) * uuu.nbColumns + (j-1)] = uu.data[(i-1) * uu.nbColumns + (sigma_j_col_idx-1)] / sigma_j;
			}
		}
		for (var i = 1; i <= n; ++i) {
			vvv.data[(i-1) * vvv.nbColumns + (j-1)] = vv.data[(i-1) * vv.nbColumns + (sigma_j_col_idx-1)];
		}		
	}
	if (svdForm == 'thin') {
		return [uuu, sss, vvv];
	}
	var qr = Matrix_.qrDecomposition(uuu);
	var q = qr[0];
	
	var uuuu = Matrix_.zeros(m, m);
	var ssss = Matrix_.zeros(m, n);
	for (var j = 1; j <= n; ++j) {
		var sigma_j_col_idx = sigmas_idx[j-1];
		var sigma_j = sigmas[sigma_j_col_idx-1];
		ssss.data[(j-1) * ssss.nbColumns + (j-1)] = sigma_j;
		if (sigma_j != 0) {
			for (var i = 1; i <= m; ++i) {
				uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = uuu.data[(i-1) * uuu.nbColumns + (j-1)];
			}
		}
		else {
			for (var i = 1; i <= m; ++i) {
				uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = q.data[(i-1) * q.nbColumns + (j-1)];
			}
		}
	}
	for (var j = n+1; j <= m; ++j) {
		for (var i = 1; i <= m; ++i) {
			uuuu.data[(i-1) * uuuu.nbColumns + (j-1)] = q.data[(i-1) * q.nbColumns + (j-1)];
		}
	}
	if (svdForm == 'full') {
		return [uuuu, ssss, vvv];
	}
}
Matrix_.eig = function(A, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var epsSymmetric = opt.epsSymmetric
	if (epsSymmetric == undefined) {
		epsSymmetric = 1e-12;
	}
	var maxIterations = opt.maxIter;
	if (maxIterations == undefined) {
		maxIterations = 100;
	}
	var sortedEigenvalues = opt.sortedEigenvalues;
	if (sortedEigenvalues == undefined) {
		sortedEigenvalues = false;
	}
	if (!(A instanceof Matrix_)) {
		throw new Error('input must be a matrix');
	}
	if (!A.isSymmetric(epsSymmetric)) {
		throw new Error('input matrix must be symmetric');
	}
	var n = A.nbRows;
	var aa = new Matrix_(A);
	var v = Matrix_.identity(n); // represents V
	var d = aa.diagonal(); // represents the vector D
	var b = new Matrix_(d);
	var z = Matrix_.zeros(n, 1);
	var iter = 0;
	while (true) {
		++iter;
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		var sm = 0;
		for (var p = 1; p <= n-1; ++p) {
			for (var q = p+1; q <= n; ++q) {
				sm += Math.abs( aa.data[(p-1) * aa.nbColumns + (q-1)] ); // a[p,q]
			}
		}
		if (sm == 0) {
			break;
		}
		var thresh = iter < 4 ? 0.2*sm/(n*n) : 0.0;
		for (var p = 1; p <= n-1; ++p) {
			for (var q = p+1; q <= n; ++q) {
				var g = 100 * Math.abs( aa.data[(p-1) * aa.nbColumns + (q-1)] ); // a[p,q]
				if ( iter > 4 && 
				     Math.abs(d.data[p-1]) + g == Math.abs(d.data[p-1]) && 
					 Math.abs(d.data[q-1]) + g == Math.abs(d.data[q-1]) ) {
					aa.data[(p-1) * aa.nbColumns + (q-1)] = 0;
				}
				else {
					if (Math.abs( aa.data[(p-1) * aa.nbColumns + (q-1)] ) > thresh) {
						var h = d.data[q-1] - d.data[p-1];
						var t;
						if (Math.abs(h) + g == Math.abs(h)) {
							t = aa.data[(p-1) * aa.nbColumns + (q-1)] / h;
						}
						else {
							var theta = 0.5 * h / aa.data[(p-1) * aa.nbColumns + (q-1)];
							t = 1 / ( Math.abs(theta) + Math.sqrt(1 + theta * theta) );
							if (theta < 0) {
								t = -t;
							}
						}
						var c = 1 / Math.sqrt( 1 + t * t );
						var s = t * c;
						var tau = s / ( 1 + c );
						var h = t * aa.data[(p-1) * aa.nbColumns + (q-1)];
						z.data[p-1] -= h; 
						z.data[q-1] += h;
						d.data[p-1] -= h; 
						d.data[q-1] += h;
						aa.data[(p-1) * aa.nbColumns + (q-1)] = 0;
						for (j = 1; j <= p-1; ++j) {
							var g = aa.data[(j-1) * aa.nbColumns + (p-1)];
							var h = aa.data[(j-1) * aa.nbColumns + (q-1)];
							aa.data[(j-1) * aa.nbColumns + (p-1)] = g - s * ( h + g * tau );
							aa.data[(j-1) * aa.nbColumns + (q-1)] = h + s * ( g - h * tau );
						}
						for (j = p+1; j <= q-1; ++j) {
							var g = aa.data[(p-1) * aa.nbColumns + (j-1)];
							var h = aa.data[(j-1) * aa.nbColumns + (q-1)];
							aa.data[(p-1) * aa.nbColumns + (j-1)] = g - s * ( h + g * tau );
							aa.data[(j-1) * aa.nbColumns + (q-1)] = h + s * ( g - h * tau );							
						}
						for (j = q+1; j <= n; ++j) {
							var g = aa.data[(p-1) * aa.nbColumns + (j-1)];
							var h = aa.data[(q-1) * aa.nbColumns + (j-1)];
							aa.data[(p-1) * aa.nbColumns + (j-1)] = g - s * ( h + g * tau );
							aa.data[(q-1) * aa.nbColumns + (j-1)] = h + s * ( g - h * tau );							
						}
						for (var j = 1; j <= n; ++j) {
							var g = v.data[(j-1) * v.nbColumns + (p-1)];
							var h = v.data[(j-1) * v.nbColumns + (q-1)];
							v.data[(j-1) * v.nbColumns + (p-1)] = g - s * ( h + g * tau );
							v.data[(j-1) * v.nbColumns + (q-1)] = h + s * ( g - h * tau );
						}
					}
				}
			}
		}
		for (var p = 1; p <= n; ++p) {
			b.data[p-1] += z.data[p-1];
			d.data[p-1] = b.data[p-1];
			z.data[p-1] = 0;
		}
	}
	var dd = d;
	var vv = v ;
	if (sortedEigenvalues) {
		var r = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
		for (var k = 1; k <= n; ++k) {
			r[k-1] = k;
		}
		r.sort(function(a, b) { return d.data[b-1] - d.data[a-1]; });
		dd = Matrix_.zeros(n, 1);
		vv = Matrix_.zeros(n, n);
		for (var k = 1; k <= n; ++k) {
			var idx = r[k-1];
			dd.data[k-1] = d.data[idx-1];
			for (var j = 1; j <= n; ++j) {
				vv.data[(j-1) * vv.nbColumns + (k-1)] = v.data[(j-1) * v.nbColumns + (idx-1)];
			}
		}
	}		
	return [vv, dd];
}
Matrix_.nullSpace = function(A, opt) {	
    if (opt === undefined) {
        opt = {};
    }
    var eps = opt.eps || undefined;
    if (!(A instanceof Matrix_)) {
        throw new Error('first input must be a matrix');
    }   
    var m = A.nbRows;
    var n = A.nbColumns;
    var a_ns = null;
    if (m >= n) {
		var svd = Matrix_.svdDecomposition(A, {maxIter: -1});		
        var u = svd[0];
        var s = svd[1];
        var v = svd[2];
		if (eps === undefined) {
			eps = m * (nextUp_(s.data[0]) - s.data[0]);
		}
        var r;
        for (r = n; r >= 1; --r) {
            if (s.data[(r-1)*s.nbColumns + (r-1)] > eps) {
                break;
            }
        }
        var p = r;
        if (p == n) {
            a_ns = Matrix_.zeros(n, 1);
        }
        else {
            a_ns = new Matrix_.zeros(n, n - p);
            for (var j = p + 1; j <= n; ++j) {
                for (var i = 1; i <= n; ++i) {
                    a_ns.data[(i-1)*a_ns.nbColumns + ((j-(p+1)+1)-1)] = v.data[(i-1)*v.nbColumns + (j-1)];
                }
            }
        }       
    }
    else {
        var svd = Matrix_.svdDecomposition(A.transpose(), {maxIter: -1, svdForm: 'full'});
        var u = svd[0];
        var s = svd[1];
        var v = svd[2];
		if (eps === undefined) {
			eps = n * (nextUp_(s.data[0]) - s.data[0]);
		}
        var r;
        for (r = m; r >= 1; --r) {
            if (s.data[(r-1)*s.nbColumns + (r-1)] > eps) {
                break;
            }
        }
        var p = r;
        if (p == m) {
            a_ns = Matrix_.zeros(n, 1);
        }
        else {
            var ta_rg = new Matrix_.zeros(n, p);
            for (var j = 1; j <= p; ++j) {
                for (var i = 1; i <= n; ++i) {
                    ta_rg.data[(i-1)*ta_rg.nbColumns + (j-1)] = u.data[(i-1)*u.nbColumns + (j-1)];
                }
            }
            var qr = Matrix_.qrDecomposition(ta_rg);
            var q = qr[0];
            var a_ns = new Matrix_.zeros(n, n-p);
            for (var j = p+1; j <= n; ++j) {
                for (var i = 1; i <= n; ++i) {
                    a_ns.data[(i-1)*a_ns.nbColumns + ((j-(p+1)+1)-1)] = q.data[(i-1)*q.nbColumns + (j-1)];
                }
            }           
        }
    }
    return a_ns;
}
Matrix_.linsolveBackSubstitution = function(A, b, out) {
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	if (!A.isSquare()) {
		throw new Error('matrix is not square: ' + '(' + A.nbRows + ',' + A.nbColumns + ')');
	}
	if (A.nbRows !== b.nbRows) {
		throw new Error('matrix and second member sizes do not match: ' + '(' + A.nbRows + ',' + A.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns !== 1) {
		throw new Error('b is not a vector: ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	var n = A.nbColumns;
	var x = allocateMatrix_(n, 1, out); // the solution vector
	for (var i = n; i >= 1; --i) {
		x.data[i-1] = b.data[i-1];
		for (var j = i + 1; j <= n; ++j) {
			x.data[i-1] = x.data[i-1] - A.data[(i-1) * A.nbColumns + (j-1)] * x.data[j-1];
		}
		var a_ii = A.data[(i-1) * A.nbColumns + (i-1)];
		if (a_ii == 0) {
			throw new Error('input matrix is not invertible: zero diagonal coefficient at index ' + i);
		}
		else {
			x.data[i-1] = x.data[i-1] / a_ii;
		}
	}
	return x;
}
Matrix_.randomOrthogonal = function(n) {
	var z = Matrix_.normrnd(n, n); // Unlike the first reference, and like the second reference, no division by SQRT(2) is made
	var qr = Matrix_.qrDecomposition(z);
	var q = qr[0];
	var r = qr[1];
	var lambda = Matrix_.fill(1, n, function(i,j) { return r.data[(j-1) * r.nbColumns + (j-1)] >= 0 ? 1 : -1; });
	return Matrix_.elementwiseProduct(q, lambda);
}
Matrix_.randomCorrelation = function(n, opt) {
	if (opt === undefined) {
		opt = { };
	}
	if (opt.eps === undefined) {
		opt.eps = 1e-14;
	}
	if (opt.epsLambda === undefined) {
		opt.epsLambda = 1e-8;
	}
	if (opt.lambda === undefined) {
		opt.lambda = new simplexRandomSampler_(n).sample();
		for (var i = 0; i < n; ++i) {
			opt.lambda[i] *= n;
		}
	}
	if (opt.lambda.length != n) {
		throw new Error('number of input eigenvalues not equal to ' + n + ', but to ' + opt.lambda.length);
	}
	var sum_lambda = 0;
	for (var i = 0; i < n; ++i) {
		sum_lambda += opt.lambda[i];
	}
	if (Math.abs(sum_lambda - n) > opt.epsLambda) {
		throw new Error('input eigenvalues not summing to ' + n);
	}
	var eps = opt.eps;
	var lambda = Matrix_.fill(1, n, function(i,j) { return opt.lambda[j-1]; });
	var u = Matrix_.randomOrthogonal(n);
	var a = Matrix_.axty(1, Matrix_.elementwiseProduct(u, lambda), u);
	for (var k = 1; k < n; ++k) {
		var converged = true;
		for (var l = 1; l <= n; ++l) {
			var a_ll = a.data[(l-1) * a.nbColumns + (l-1)];
			
			if (Math.abs(a_ll - 1) > eps) {
				converged = false;
				break;
			}
			else {
				if (a_ll != 1) {
					a.data[(l-1) * a.nbColumns + (l-1)] = 1;
				}
			}
		}
		if (converged) {
			break;
		}
		var i = -1; // first index for which a_ii < 1
		var j = -1; // last index for which a_jj >1
		for (var l = 1; l <= n; ++l) {
			if (a.data[(l-1) * a.nbColumns + (l-1)] < 1 && i === -1) {
				i = l;
			}
			if (a.data[(l-1) * a.nbColumns + (l-1)] > 1) {
				j = l;
			}
		}
		if (i > j) {
			i = -1; // first index for which a_ii > 1
			j = -1; // last index for which a_jj < 1
			for (var l = 1; l <= n; ++l) {
				if (a.data[(l-1) * a.nbColumns + (l-1)] > 1 && i === -1) {
					i = l;
				}
				if (a.data[(l-1) * a.nbColumns + (l-1)] < 1) {
					j = l;
				}
			}
		}
		var a_ii = a.data[(i-1) * a.nbColumns + (i-1)];
		var a_ij = a.data[(i-1) * a.nbColumns + (j-1)];
		var a_jj = a.data[(j-1) * a.nbColumns + (j-1)];
		var t = (a_ij + Math.sqrt(a_ij*a_ij - (a_ii-1)*(a_jj-1) )) / (a_jj - 1);
		var c = 1/Math.sqrt(1 + t*t);
		var s = c*t;
		for (var l = 1; l <= n; ++l) {
			var t1 = a.data[(i-1) * a.nbColumns + (l-1)]; // a(i,l)
			var t2 = a.data[(j-1) * a.nbColumns + (l-1)]; // a(j,l)
			
			a.data[(i-1) * a.nbColumns + (l-1)] = c*t1 - s*t2;
			a.data[(j-1) * a.nbColumns + (l-1)] = s*t1 + c*t2;
		}
		for (var l = 1; l <= n; ++l) {
			var t1 = a.data[(l-1) * a.nbColumns + (i-1)]; // a(l,i)
			var t2 = a.data[(l-1) * a.nbColumns + (j-1)]; // a(l,j)
			
			a.data[(l-1) * a.nbColumns + (i-1)] = c*t1 - s*t2;
			a.data[(l-1) * a.nbColumns + (j-1)] = s*t1 + c*t2;
		}
		a.data[(i-1) * a.nbColumns + (i-1)] = 1;
	}
	for (var l = 1; l <= n; ++l) {
		a.data[(l-1) * a.nbColumns + (l-1)] = 1;
	}
	for (var i = 1; i <= n; ++i) {
		for (var j = i+1; j <= n; ++j) {
			var a_ij = a.data[(i-1) * a.nbColumns + (j-1)];
			var a_ji = a.data[(j-1) * a.nbColumns + (i-1)];
			var a_ij_fixed = Math.min(1, Math.max(-1, (a_ij + a_ji) / 2));
			a.data[(i-1) * a.nbColumns + (j-1)] = a_ij_fixed;
			a.data[(j-1) * a.nbColumns + (i-1)] = a_ij_fixed;
		}
	}
	return a;
}
Matrix_.linsolveExtendedKaczmarz = function(A, b, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-12;
	var maxIterations = opt.maxIter || 100000;
	var randomized = false;
	if (opt.randomized !== undefined) {
		randomized = opt.randomized;
	}
	if (!(A instanceof Matrix_)) {
		throw new Error('first input must be a matrix');
	}
	if (!(b instanceof Matrix_)) {
		throw new Error('second input must be a matrix');
	}
	
	if (A.nbRows !== b.nbRows) {
		throw new Error('matrix and second member sizes do not match: ' + '(' + A.nbRows + ',' + A.nbColumns + 
		') - ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	if (b.nbColumns !== 1) {
		throw new Error('b is not a vector: ' + '(' + b.nbRows + ',' + b.nbColumns + ')');
	}
	var m = A.nbRows;
	var n = A.nbColumns;
	var x_k = Matrix_.zeros(n, 1); // the current solution
	var z_k = new Matrix_(b); // the current "adjusted" b
	
	var x_res = Matrix_.zeros(m, 1); // the x residuals vector
	var b_res = Matrix_.zeros(m, 1); // the b residuals vector
	var a_x_k = Matrix_.zeros(m, 1); // A*x_k
	var a_frob_norm = A.matrixNorm('frobenius');
	var a_frob_norm_sq = a_frob_norm * a_frob_norm;	
	if (a_frob_norm === 0) {
		return x_k;
	}
	var a_rows_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(m) : new Array(m);
	for (var i = 1; i <= m; ++i) {
		var a_i_two_norm = A.vectorNorm('two', 'row', i);
		a_rows_two_norm_sq[i-1] = a_i_two_norm * a_i_two_norm;
	}
	var a_columns_two_norm_sq = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
	for (var j = 1; j <= n; ++j) {
		var alpha_j_two_norm = A.vectorNorm('two', 'column', j);
		a_columns_two_norm_sq[j-1] = alpha_j_two_norm * alpha_j_two_norm;
	}
	if (randomized == false) {
		var iter = 0;
		while (true) {
			++iter;
			if (maxIterations !== -1 && iter >= maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			for (var i = 1; i <= m; ++i) {
				if (a_rows_two_norm_sq[i-1] == 0) {
					continue;
				}
				var a_i_x_k = 0;
				for (var j = 1; j <= n; ++j) {
					a_i_x_k += A.data[(i-1) * A.nbColumns + (j-1)] * x_k.data[(j-1) * x_k.nbColumns];
				}
				var r_k = a_i_x_k - (b.data[(i-1) * b.nbColumns + 0] - z_k.data[(i-1) * z_k.nbColumns + 0]);
				for (var j = 1; j <= n; ++j) {
					x_k.data[(j-1) * x_k.nbColumns] -= r_k / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
				}
			}
			for (var j = 1; j <= n; ++j) {
				if (a_columns_two_norm_sq[j-1] == 0) {
					continue;
				}
				var a_j_z_k = 0;
				for (var i = 1; i <= m; ++i) {
					a_j_z_k += A.data[(i-1) * A.nbColumns + (j-1)] * z_k.data[(i-1) * z_k.nbColumns + 0];
				}
				for (var i = 1; i <= m; ++i) {
					z_k.data[(i-1) * z_k.nbColumns + 0] -= a_j_z_k / a_columns_two_norm_sq[j-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
				}
			}
			var x_k_two_norm = x_k.vectorNorm('two');
			a_x_k =  Matrix_.xy(A, x_k, a_x_k);
			b_res = Matrix_.xmy(b, z_k, b_res);
			x_res = Matrix_.xmy(a_x_k, b_res, x_res);
			var x_res_two_norm = x_res.vectorNorm('two');
			if (x_res_two_norm <= eps * a_frob_norm * x_k_two_norm) {
				break;
			}
		}
	}
	else {
		var ta_z_k = Matrix_.zeros(n, 1); // A^t*z_k
		var q = typeof Float64Array === 'function' ? new Float64Array(m) : new Array(m);
		for (var i = 1; i <= m; ++i) {
			q[i-1] = a_rows_two_norm_sq[i-1]/a_frob_norm_sq;
		}

		var qSampler = new aliasMethodSampler_(q);
		var p = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var j = 1; j <= n; ++j) {
			p[j-1] = a_columns_two_norm_sq[j-1]/a_frob_norm_sq;
		}
		var pSampler = new aliasMethodSampler_(p);
		var iter = 0;
		while (true) {
			++iter;
			if (maxIterations !== -1 && iter >= maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			var i = qSampler.sample() + 1;
			var a_i_x_k = 0;
			for (var j = 1; j <= n; ++j) {
				a_i_x_k += A.data[(i-1) * A.nbColumns + (j-1)] * x_k.data[(j-1) * x_k.nbColumns];
			}
			var r_k = a_i_x_k - (b.data[(i-1) * b.nbColumns + 0] - z_k.data[(i-1) * z_k.nbColumns + 0]);
			for (var j = 1; j <= n; ++j) {
				x_k.data[(j-1) * x_k.nbColumns] -= r_k / a_rows_two_norm_sq[i-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
			var j = pSampler.sample() + 1;
			var a_j_z_k = 0;
			for (var i = 1; i <= m; ++i) {
				a_j_z_k += A.data[(i-1) * A.nbColumns + (j-1)] * z_k.data[(i-1) * z_k.nbColumns + 0];
			}
			for (var i = 1; i <= m; ++i) {
				z_k.data[(i-1) * z_k.nbColumns + 0] -= a_j_z_k / a_columns_two_norm_sq[j-1] * A.data[(i-1) * A.nbColumns + (j-1)]; 
			}
			if (iter % 8 * Math.min(m, n) === 0) {
				var x_k_two_norm = x_k.vectorNorm('two');
				a_x_k =  Matrix_.xy(A, x_k, a_x_k);
				b_res = Matrix_.xmy(b, z_k, b_res);
				x_res = Matrix_.xmy(a_x_k, b_res, x_res);
				var x_res_two_norm = x_res.vectorNorm('two');
				ta_z_k = Matrix_.txy(A, z_k, ta_z_k);
				var ta_z_k_two_norm = ta_z_k.vectorNorm('two');
				
				if (x_res_two_norm <= eps * a_frob_norm * x_k_two_norm && 
					ta_z_k_two_norm <= eps * a_frob_norm_sq * x_k_two_norm) {
					break;
				}
			}
		}
	}
	return x_k;
}
function meanVector(arr, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "sample-mean";
	}
	var meanMethod = opt.method;
	if (meanMethod != "sample-mean" && meanMethod != "linear-shrinkage" ) {
			throw new Error('unsupported mean vector computation method');
	}
	var nbSeries = arr.length;
	var nbObservations = arr[0].length;
	var means;
	if (meanMethod == "sample-mean") {
		var meanVector = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]); });
		means = meanVector;
	}
	else if (meanMethod == "linear-shrinkage") {
		var meanVector = Matrix_.fill(nbSeries, 1, function(i,j) { return mean_(arr[i-1]); });
		var nu = mean_(meanVector.toArray());
		var prior = Matrix_.fill(nbSeries, 1, function(i,j) { return nu; });
		var variances = Matrix_.fill(nbSeries, 1, function(i,j) { return covariance_(arr[i-1], arr[i-1]); });
		var sigma_sq_b = variances.sum()/nbSeries;
		var alpha = nbSeries/nbObservations * sigma_sq_b / ( nbSeries/nbObservations * sigma_sq_b + Math.pow(Matrix_.xmy(prior, meanVector).vectorNorm('two'), 2) );
		means = Matrix_.axpby(alpha, prior, 1-alpha, meanVector);
	}
	else {
		throw new Error('internal error: unsupported mean vector computation method');
	}
	return means;
}
function randomMeanVector(n, sigma) {
	var s = sigma == undefined ? 1 : sigma;
	var meanVect = Matrix_.fill(n, 1, function(i,j) { return normrnd_(0, s); });
	return meanVect;
}
function perturbedMeanVector(meanVect, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "multiplicative-noise";
	}
	var method = opt.method;
	if (method != "multiplicative-noise") {
			throw new Error('unsupported perturbation method');
	}
	
	var sigma = opt.sigma;
	if (opt.sigma === undefined) {
		sigma = 0.05;
	}
	var meanVect = new Matrix_(meanVect);
	if (!meanVect.isVector()) {
		throw new Error('input must be a vector');
	}
	var perturbedMeanVect = meanVect.elemMap(function(i,j,val) { return val * (1 + sigma*normrnd_(0,1));})
	return perturbedMeanVect;
}
function randomVariances(n, sigma) {
	var s = sigma == undefined ? 1 : sigma;
	return Matrix_.fill(n, 1, function(i,j) { return pnormrnd_(0, s); });	
}
function perturbedVariances(variancesVect, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "multiplicative-noise";
	}
	var method = opt.method;
	if (method != "multiplicative-noise") {
			throw new Error('unsupported perturbation method');
	}
	
	var sigma = opt.sigma;
	if (opt.sigma === undefined) {
		sigma = 0.05;
	}
	var variancesVect = new Matrix_(variancesVect);
	if (!variancesVect.isVector()) {
		throw new Error('input must be a vector');
	}
	var perturbedVariancesVect = variancesVect.elemMap(function(i,j,val) { 
	                                                       var pvar = val * (1 + sigma*normrnd_(0,1)); 
	                                                       while (pvar <= 0) { 
														       pvar = val * (1 + sigma*normrnd_(0,1));
														   }
														   return pvar;
													   })
	return perturbedVariancesVect;
}
function BitSet_(arr) {
	function fromArray(arr) {
		for (var i = 0; i < arr.length; ++i) {
			that.set(arr[i]);
		}
	}
	if (!(this instanceof BitSet_)) {
      return new BitSet_();
    }
	this.WORD_LENGTH = 32;
	this.WORD_LOG = 5;
	this.words = typeof Uint32Array === 'function' ? new Uint32Array(0) : new Array(0);
	var that = this;
	if (arr instanceof Array || (typeof Uint32Array === 'function' && arr instanceof Uint32Array)) {
		fromArray(arr);
	}
	this.iterator = function() {
		this.w_idx = -1;
		this.w = 0;
		while (this.w_idx < that.words.length && this.w == 0) {
			++this.w_idx;
			this.w = that.words[this.w_idx];
		}
		this.next = function() {
			if (this.w_idx == that.words.length) {
				return 0;
			}			
			var t = this.w & -this.w;
			var value = (this.w_idx << that.WORD_LOG) + BitSet_.populationCount((t - 1) | 0);
			this.w ^= t;
			while (this.w_idx < that.words.length && this.w == 0) {
				++this.w_idx;
				this.w = that.words[this.w_idx];
			}
			return value;
		}
	}
}
BitSet_.populationCount = function(x) {
	x = x - ((x >>> 1) & 0x55555555);
	x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
	x = (x + (x >>> 4)) & 0x0F0F0F0F;
	x = x + (x >>> 8);
	x = x + (x >>> 16);
	return x & 0x0000003F;
};

BitSet_.prototype = {
    constructor: BitSet_,
	resize : function(idx) {
		var c = this.words.length;
		if ((c << this.WORD_LOG) > idx) {
			return;
		}
	    var count = (idx + this.WORD_LENGTH) >>> this.WORD_LOG;
		if (typeof Uint32Array === 'function') {
			var words_new = new Uint32Array(count); // the new typed array is (automatically) initialized with 0s
			words_new.set(this.words); // copy the previous words into the new typed array
			this.words = words_new;
		}
		else {  
			this.words[count-1] = 0; // copy (automatically) the previous words into the new array
			for (var i = c; i < count-1; ++i) {
				this.words[i] = 0;
			}
		}
	},
	toString: function () {
		var fullStr = '';
		var c = this.words.length;
		for (var i = 0; i < c; ++i) {
			var str = "";
			if (typeof Uint32Array === 'function') {
				str = this.words[i].toString(2);
			}
			else {
				str = (this.words[i] >>> 0).toString(2);
			}
			fullStr += str.length >= this.WORD_LENGTH ? str : new Array(this.WORD_LENGTH - str.length + 1).join('0') + str;
		}
		return fullStr;
	},
	set: function(idx) {
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
		this.words[idx >>> this.WORD_LOG] |= (1 << idx);
		return this;
	},
	setRange: function (idxFrom, idxTo) {
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idxTo) {
			this.resize(idxTo);
		}
		for (var i = idxFrom; i <= idxTo; ++i) {
			this.words[i >>> this.WORD_LOG] |= (1 << i);
		}
		return this;
	},
	unset: function(idx) {
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
		this.words[idx >>> this.WORD_LOG] &= ~(1 << idx);
		return this;
	},
	get: function(idx) {
		return (this.words[idx  >>> this.WORD_LOG] & (1 << idx)) !== 0;
	},	
	clear: function() {
		this.words = typeof Uint32Array === 'function' ? new Uint32Array(0) : new Array(0);
		return this;
	},
	flip: function(idx) {
		var c = this.words.length;
		if ((c << this.WORD_LOG) <= idx) {
			this.resize(idx);
		}
	  this.words[idx >>> this.WORD_LOG] ^= 1 << idx;
	  return this;
	},
	isEmpty: function() {
	  var c = this.words.length;
	  for (var  i = 0; i < c; ++i) {
		if (this.words[i] != 0) {
			return false;
		}
	  }
	  return true;
	},
	nbSetBits: function() {
	  var s = 0;
	  var c = this.words.length;
	  for (var i = 0; i < c; ++i) {
		s += BitSet_.populationCount(this.words[i] | 0);
	  }
	  return s;
	},
	toArray: function() {
	  var arr = new Array(this.nbSetBits());
	  var pos = 0;
	  var c = this.words.length;
	  for (var i = 0; i < c; ++i) {
		var w = this.words[i];
		while (w != 0) {
		  var t = w & -w;
		  arr[pos++] = (i << this.WORD_LOG) + BitSet_.populationCount((t - 1) | 0);
		  w ^= t;
		}
	  }
	  return arr;
	},

};
function aliasMethodSampler_(p) {
	this.prob = typeof Float64Array === 'function' ? new Float64Array(p.length) : new Array(p.length);
	this.alias = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
    var avgProb = 1 / p.length;
	var small = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var s = 0;
	var large = typeof Uint32Array === 'function' ? new Uint32Array(p.length) : new Array(p.length);
	var l = 0;
	for (var j = 0; j < p.length; ++j) {
		if (p[j] > avgProb) {
			large[l] = j;
			++l;
		}
		else {
			small[s] = j;
			++s;
		}
	}
	var p = p.slice(0); // local copy of the probabilities, as they are updated below
	while (s > 0 && l > 0) {
		--s;
		var j = small[s];
		
		--l;
		var k = large[l];
		this.prob[j] = p.length * p[j];
		this.alias[j] = k;
		p[k] = p[k] + (p[j] - avgProb);
		if (p[k] > avgProb) {
			large[l] = k;
			++l;
		}
		else {
			small[s]= k;
			++s;
		}
	}
	while (s > 0) {
		--s;
		this.prob[small[s]] = 1;
	}
	while (l > 0) {
		--l;
		this.prob[large[l]] = 1;
	}
    this.sample = function() {
		var u = Math.random() * this.prob.length; // Uniform real number belonging to [0..n[
		var j = Math.floor(u);
		if (u - j <= this.prob[j]) {
			return j;
		}
		else {
			return this.alias[j];
		}
    }
}
function compositionsIterator_(n, k, reuseOutputArray) {
	this.n = n;
	this.k = k;
	this.reuseOutputArray = reuseOutputArray;
	this.firstcall = true;
	this.mtc = false;
	this.r = typeof Uint32Array === 'function' ? new Uint32Array(k) : new Array(k);
	this.t = this.n;
	this.h = 0;
	this.next = function() {
		if (!this.firstcall && !this.mtc) {
			return -1;
		}
		
		if (this.firstcall) {		
			this.firstcall = false;
			this.r[0] = this.n;
			for (var i = 1; i <= this.k - 1; ++i) {
				this.r[i] = 0;
			}	
		}
		else {
			if (this.t > 1) {
				this.h = 0;
			}
			++this.h;
			this.t = this.r[this.h - 1];
			this.r[this.h - 1] = 0;
			this.r[0] = this.t - 1;
			++this.r[this.h];
		}
		this.mtc = (this.r[this.k - 1] != this.n);
		if (this.reuseOutputArray) {
			return this.r;
		}
		else {
			return this.r.slice(0);
		}
	}
}
function permutationsIterator_(n, reuseOutputArray) {
	this.n = n;
	this.reuseOutputArray = reuseOutputArray;
	this.p = typeof Uint32Array === 'function' ? new Uint32Array(this.n) : new Array(this.n);
	for (var i = 0; i < this.n; ++i) {
		this.p[i] = i+1;
	}
	this.c = typeof Uint32Array === 'function' ? new Uint32Array(this.n) : new Array(this.n);
	for (var i = 0; i < this.n; ++i) {
		this.c[i] = 0;
	}
	this.j = -1;
	this.next = function() {	
		if (this.j >= this.n) {
			return -1;
		}
		if (this.j === -1) {
			++this.j;
			if (this.reuseOutputArray) {
				return this.p;
			}
			else {
				return this.p.slice(0);
			}
		}
		else {
			while (this.j < this.n) {
				if (this.c[this.j] < this.j) {
					if (this.j % 2 === 0) { // j is even, swap P(j) and P(0)
						var tmp = this.p[this.j];
						this.p[this.j] = this.p[0];
						this.p[0] = tmp;
					}
					else { // j is odd, swap P(j) and P(c(j))
						var tmp = this.p[this.j];
						this.p[this.j] = this.p[this.c[this.j]];
						this.p[this.c[this.j]] = tmp;
					}
					++this.c[this.j];
					this.j = 0;
					if (this.reuseOutputArray) {
						return this.p;
					}
					else {
						return this.p.slice(0);
					}
				}
				else {
					this.c[this.j] = 0;
					++this.j;
				}
			}
			return -1;
		}
	}
}
function randomCompositionsIterator_(n, k) {
	this.n = n;
	this.k = k;
	this.ranksb = new randomKSubsetsIterator_(n+k-1, k-1);
	this.r = typeof Uint32Array === 'function' ? new Uint32Array(k) : new Array(k);
	this.next = function() {
		var rr = this.ranksb.next();
		for (var i = 0; i < this.k-1; ++i) {
			this.r[i] = rr[i];
		}
		this.r[this.k-1] = this.n + this.k;
		var l = 0;
		for (var i = 1; i <= this.k; ++i) {
			var m = this.r[i-1];
			this.r[i-1] = m - l - 1;
			l = m;
		}
		return this.r.slice(0);
	}
}
function randomPermutationsIterator_(n, outputArray, reuseOutputArray) {
	if (outputArray === undefined) {
		this.r = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.r[i] = i+1;
		}
	}
	else {
		this.r = outputArray;
	}
	this.rrlen = this.r.length;
	this.reuseOutputArray = reuseOutputArray;
	function uniformrv() {
		var rnd = Math.random();
		while (rnd === 0) {
			rnd = Math.random();
		}
		return rnd;
	}
	this.next = function() {
		if (this.reuseOutputArray) {
			var rr = this.r;
		}
		else {
			var rr = this.r.slice(0);
		}
		for (var m = 1; m <= this.rrlen; ++m) {
			var l = m + Math.floor( uniformrv()*(this.rrlen + 1 - m) );
			
			var tmp = rr[l-1];
			rr[l-1] = rr[m-1];
			rr[m-1] = tmp;
		}
		return rr;
	}
}
function kSubsetsIterator_(n, k, useArrayCopy) {
	this.n = n;
	this.k = k;
	this.useArrayCopy = useArrayCopy;
	this.a = typeof Uint32Array === 'function' ? new Uint32Array(k) : new Array(k);
	this.firstcall = true;
	this.mtc = false;
	this.m2 = 0;
	this.h = this.k;
	this.endval = this.n - this.k + 1;	
	this.next = function() {
		if (!this.firstcall && (!this.mtc || this.k === 0)) {
			return -1;
		}
		
		if (this.firstcall) {		
			this.firstcall = false;
		}
		else {
			if (this.m2 < this.n - this.h) {
				this.h = 0;
			}
			
			++this.h;
			this.m2 = this.a[this.k - this.h];
		}
		for (var j = 1; j <= this.h; ++j) {
			this.a[this.k + j - this.h - 1] = this.m2 + j;
		}
		this.mtc = (this.a[0] != this.endval);
		if (this.useArrayCopy) {
			return this.a.slice(0);
		}
		else {
			return this.a;
		}
	}
}
function randomKSubsetsIterator_(n, k, useArrayCopy) {
	this.n = n;
	this.k = k;
	this.useArrayCopy = useArrayCopy;
	this.a = typeof Uint32Array === 'function' ? new Uint32Array(k) : new Array(k);
	this.N;
	this.nn;
	this.idx_a;
	this.selected_record;
	function uniformrv() {
		var rnd = Math.random();
		while (rnd === 0) {
			rnd = Math.random();
		}
		return rnd;
	}
	this.method_a = function() {
		var top = this.N - this.nn;
		while (this.nn >= 2) {
			var V = uniformrv();
			var S = 0;
			var quot = top / this.N;
			
			while (quot > V) {
				++S;
				--top;
				--this.N;
				quot = (quot * top) / this.N;
			}
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
			
			--this.N;
			--this.nn;
		}
		var S = Math.floor(this.N * uniformrv());
		if (S === this.N) { // the out of range value S = N must never be generated (which could be in case of roundoff errors)
			S = this.N - 1;
		}
		this.selected_record += S + 1;
		this.a[this.idx_a++] = this.selected_record;
	}
	this.method_d = function() {
		var ninv = 1/this.nn;
		var Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
		var qu1 = -this.nn + 1 + this.N;
		var negalphainv = -13;
		var threshold = -negalphainv * this.nn;
		
		while (this.nn > 1 && threshold < this.N) {
			var nmin1inv = 1 / (-1 + this.nn);
			
			var X;
			var S;
			while (true) {
				while(true) {
					X = this.N 	* (-Vprime + 1);
					S = Math.floor(X);
					if (S < qu1) {
						break;
					}
					Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
				}
				var U = uniformrv();
				var y1 = Math.pow(U * this.N / qu1, nmin1inv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
				Vprime = y1 * (-X/this.N + 1) * (qu1 / (-S + qu1));
				if (Vprime <= 1) {
					break;
				}
				var y2 = 1;
				var top = -1 + this.N;
				
				var bottom;
				var limit;
				if (-1 + this.nn > S) {
					bottom = -this.nn + this.N;
					limit = -S + this.N;
				}
				else {
					bottom = -1 - S + this.N;
					limit = qu1;
				}
				
				for (var t = -1 + this.N; t >= limit; --t) {
					y2 = y2 * top / bottom;
					--top;
					--bottom;
				}
				
				if (this.N / (-X + this.N) >= y1 * Math.pow(y2, nmin1inv)) { // Math.exp(Math.log(a) * b) === Math.pow(a, b)
					Vprime = Math.pow(uniformrv(), nmin1inv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
					break;
				}
				Vprime = Math.pow(uniformrv(), ninv); // Math.exp(Math.log(a) * b) === Math.pow(a, b)
			}
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
			this.N = -S + (-1 + this.N);
			--this.nn;
			ninv = nmin1inv;
			qu1 = -S + qu1;
			threshold = threshold + negalphainv;
		}
		if (this.nn > 1) {
			this.method_a();
		}
		else {
			var S = Math.floor(this.N * Vprime);
			if (S === this.N) { // the out of range value S = N must never be generated (which could be in case of roundoff errors)
				S = this.N - 1;
			}
			this.selected_record += S + 1;
			this.a[this.idx_a++] = this.selected_record;
		}		
	}
	this.next = function() {
		for (var i = 0; i < this.k; ++i) {
			this.a[i] = 0;
		}
		this.N = this.n;
		this.nn = this.k;
		this.idx_a = 0;
		this.selected_record = 0;
		this.method_d();
		if (this.useArrayCopy) {
			return this.a.slice(0);
		}
		else {
			return this.a;
		}
	}
}
function subsetsIterator_(n) {
	this.n = n;
	this.firstcall = true;
	this.mtc = false;
	this.iin = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
	this.ncard = 0;
	this.next = function() {
		var nextSubset = [];
		
		if (!this.firstcall && !this.mtc) {
			return -1;
		}
		
		if (this.firstcall) {
			this.firstcall = false;
			for (var i = 0; i <= this.n - 1; ++i) {
				this.iin[i] = 0;
			}
			this.mtc = true;
		}
		else {
			var j = 0;
			if (this.ncard % 2 != 0) {
				++j;
				while (this.iin[j - 1] == 0) {
					++j;
				}
			}
			this.iin[j] = 1 - this.iin[j];
			this.ncard = this.ncard + 2*this.iin[j] - 1;
			nextSubset = typeof Uint32Array === 'function' ? new Uint32Array(this.ncard) : new Array(this.ncard);
			var idx = 0;
			for (var i = 0; i <= this.n - 1; ++i) {
				if (this.iin[i] == 1) {
					nextSubset[idx++] = i + 1;
				}
			}
			this.mtc = (this.ncard != this.iin[this.n -1]);
		}
		return nextSubset;
	}
}
function binomial_(n, k) {
    if (n < 0) {
        throw new Error('n must be a positive integer');
    }
    if (k < 0) {
        throw new Error('k must be a positive integer');
    }
    if (k > n) {
        throw new Error('k must be less than or equal to n');
    }
	var val = 1;
    for (var i = 1; i <= k; ++i) {
        val *= (n + 1 - i);
		val /= i; // Note: separating the computation of val in two steps guarantees (unless compiler optimisations) that val is an integer.
    }
    return val;
}
function geometricCenter_(x) {
	var m = x.length;
	var n = x[0].nbRows;
	var y = Matrix_.zeros(n, 1);
	for (var i = 1; i <= n; ++i) {
		var sum_i = 0.0;
		for (var k = 0; k < m; ++k) {
			sum_i += x[k].getValue(i, 1);
		}
		var tmpMean_i = sum_i/m;
		var sumDiff_i = 0.0;
		for (var k = 0; k < m; ++k) {
			sumDiff_i += (x[k].getValue(i, 1) - tmpMean_i);
		}
		y.setValue(i, 1,
		           (sum_i + sumDiff_i)/m);
	}
	return y;
}
function geometricMedian_(x) {
	function f_eps(y) {
		var sum = 0.0;

		for (var k = 0; k < m; ++k) {
			var y_m_x_k = Matrix_.xmy(y, x[k], tmp_vec_n);
			var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');
			
			sum += hypot_(y_m_x_k_two_norm, eps_f);
		}
		
		return sum;
	}
	function gradf_eps(y) {
		var res = Matrix_.zeros(n, 1);
		
		for (var k = 0; k < m; ++k) {
			var y_m_x_k = Matrix_.xmy(y, x[k], tmp_vec_n);
			var y_m_x_k_two_norm = y_m_x_k.vectorNorm('two');

			res = Matrix_.axpby(1, res, 1/hypot_(y_m_x_k_two_norm, eps_f), y_m_x_k, res);
		}
		
		return res;
	}
	var m = x.length; // the number of points provided in input
	var n = x[0].nbRows; // the dimension of each point provided in input

	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var x0 = geometricCenter_(x);
	var g = function(x) {
		return 0;
	}
	var proxg = function(x, mu) {
		return x;
	}
	var eps_f = 1 / m;
	var sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, x0, {eps: 1e-4, maxIter: -1, maxLine: -1});

	eps_f = 1e-1 / m;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-4, maxIter: -1, maxLine: -1});

	eps_f = 1e-2 / m;
	sol = ccpsolveFISTA_(f_eps, gradf_eps, g, proxg, sol[0], {eps: 1e-4, maxIter: -1, maxLine: -1});
	return sol[0];
}
function lineSegmentEuclidianProjection_(x, a, b) {
	var eps = 1e-8; // the numerical zero
	var x = new Matrix_(x);
	
	var a = new Matrix_(a);
	var b = new Matrix_(b);
	if (a.nbRows != x.nbRows || b.nbRows != x.nbRows) {
		throw new Error('incompatible dimensions: ' + x.nbRows + ', ' + a.nbRows + ', ' + b.nbRows);
	}
	
	var u = Matrix_.xmy(b, a); // b-a
	var v = Matrix_.xmy(x, a); // x-a
	var u_two_norm = u.vectorNorm('two');
	if (u_two_norm <= eps) {
		return a.toArray();
	}
	var t = Math.max(0, Math.min(1, Matrix_.vectorDotProduct(u, v) / (u_two_norm * u_two_norm)));
	var proj_x = Matrix_.axpby(1, a, t, u);
	return proj_x.toArray();
}
function quantile_(x, p, inputArraySorted) {
	var xx = x;
	if (!inputArraySorted) {
		xx = x.slice(0).sort(function(a, b) { return a-b });
	}
    var i = p * (x.length - 1);
	var i0 = Math.floor(i);
	var res;
    if (i0 == i) {
    	res = xx[i0];
    } else {
        res = xx[i0] + (xx[i0 + 1] - xx[i0]) * (i - i0);
    }
    return res;
}
function boxRandomSampler_(n, l, u, reuseOutputArray) {
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.reuseOutputArray = reuseOutputArray;
	this.l = l;
	if (!this.l) {
		this.l = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n); 
		for (var i = 0; i < this.n; ++i) {
			this.l[i] = 0; // default lower bounds values
		}
	}
	this.u = u;
	if (!this.u) {
		this.u = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n); 
		for (var i = 0; i < this.n; ++i) {
			this.u[i] = 1; // default upper bounds values
		}
	}
	for (var i = 0; i < this.n; ++i) {
		if (this.l[i] > this.u[i]) {
			throw new Error('empty box detected: lower bound strictly greater than upper bound');
		}
	}
	this.sample = function() {
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = Math.random() * (this.u[i] - this.l[i]) + this.l[i];
		}
		if (this.reuseOutputArray) {
			return this.x;
		}
		else {
			return this.x.slice(0);
		}
	}
}
function boxGridSampler_(n, k, l, u, reuseOutputArray) {
	this.n = n;
	this.k = k;
	this.reuseOutputArray = reuseOutputArray;
	this.arrs = typeof UInt32Array === 'function' ? new UInt32Array(this.n) : new Array(this.n);
	for (var i = 0; i < this.n; ++i) {
		this.arrs[i] = 0;
	}
	this.x = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n);
	this.l = l;
	if (!this.l) {
		this.l = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n); 
		for (var i = 0; i < this.n; ++i) {
			this.l[i] = 0; // default lower bounds values
		}
	}
	this.u = u;
	if (!this.u) {
		this.u = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n); 
		for (var i = 0; i < this.n; ++i) {
			this.u[i] = 1; // default upper bounds values
		}
	}
	this.nbGridPoints = typeof UInt32Array === 'function' ? new UInt32Array(this.n) : new Array(this.n);
	this.gridSize = typeof Float64Array === 'function' ? new Float64Array(this.n) : new Array(this.n);
	for (var i = 0; i < this.n; ++i) {
		var k_i;
		if (this.k[i]) {
			 k_i = this.k[i];
		}
		else {
			k_i = this.k;
		}
		this.nbGridPoints[i] = k_i;
		if (k_i === 1) {
			this.gridSize[i] = 0;
		}
		else {
			this.gridSize[i] = (this.u[i] - this.l[i]) / (this.nbGridPoints[i] - 1);
		}
	}
	for (var i = 0; i < this.n; ++i) {
		var lowerBound = this.l[i];
		var upperBound = this.u[i];
		if (lowerBound > upperBound) {
			throw new Error('infeasible problem detected: lower bound ' + lowerBound + ' strictly greater than upper bound ' + upperBound);
		}
		else if (lowerBound < upperBound) {
			if (this.gridSize[i] <= 0) {
				throw new Error('incorrect number of grid points to generate on the interval [' + lowerBound + ', ' + upperBound + ']: ' + this.nbGridPoints[i]);
			}
		}
	}
	this.sample = function() {
		if (this.arrs[0] === this.nbGridPoints[0]) {
			return -1;
		}
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.l[i] + this.arrs[i] * this.gridSize[i];
		}
		var change = true;
		var r = this.n - 1;
		while (change && r >= 0) {
			++this.arrs[r];
            if (this.arrs[r] === this.nbGridPoints[r]) {        
				if (r === 0) {
					break;
				}
				this.arrs[r] = 0;
                change = true;
            }
            else {
				change = false;
			}
			--r;
        }
		if (this.reuseOutputArray) {
			return this.x;
		}
		else {
			return this.x.slice(0);
		}
	}
}
function permutationEntropy_(x, m) {
	var n = x.length;
	var permutationsCounter = {};
	var nbPermutations = 0;
	var idxStart = 0;
	var idxEnd = n - m + 1;	
	var xxPermutation = typeof UInt32Array === 'function' ? new UInt32Array(m) : new Array(m);
	for (var i = idxStart; i < idxEnd; ++i) {
		var xx = x.slice(i, i + m);
		for (var j = 0; j < m; ++j) {
			xxPermutation[j] = j + 1;
		}
		xxPermutation.sort(function(a, b) {
			return xx[a-1] - xx[b-1];
		});
		var permutationIndex = xxPermutation.toString();
		permutationsCounter[permutationIndex] = (permutationsCounter[permutationIndex] || 0) + 1;
		++nbPermutations;
	}
	var pE = 0;
	for (var key in permutationsCounter) {
		if (permutationsCounter.hasOwnProperty(key)) {           
			var permutationCounter = permutationsCounter[key];
			var permutationFrequency = permutationCounter / nbPermutations;
			
			pE += permutationFrequency * Math.log(permutationFrequency)
		}
	}
	return -pE / Math.log(2);
}
function max_(x, compareFunction) {
	var defaultCompareFct = function (a, b) {
		return a - b;
	};
	var compareFunction = compareFunction || defaultCompareFct;
	
	var n = x.length;
	var maxValue = x[0];
	var maxValueIdx = 0;
	for (var i = 1; i < n; ++i) {
		if (compareFunction(x[i], maxValue) > 0) {
			maxValue = x[i];
			maxValueIdx = i;
		}
	}
	return [maxValue, maxValueIdx];
}
function median_(x, compareFunction) {
	var n = x.length;
	var xx = x.slice(); // to avoid altering the array x
	return select_(xx, Math.ceil(n/2), compareFunction);
}
function select_(x, k, compareFunction) {
	var defaultCompareFct = function (a, b) {
		return a - b;
	};
	var compareFunction = compareFunction || defaultCompareFct;
	
	var n = x.length;

	var cutoff = 600;
	var cs = 0.5; // Brown's version: cs = 0.5
	var csd = 0.5; // Brown's version: cs = 0.1
	var nstack = 10;
	var stack_1 = typeof UInt32Array === 'function' ? new UInt32Array(nstack) : new Array(nstack);
	var stack_2 = typeof UInt32Array === 'function' ? new UInt32Array(nstack) : new Array(nstack);
	var jstack = 0; // number of elements in the stacks stack_1 and stack_2
	
	var l = 0;
    var r = n - 1; // -1 because Fortran convention is to start the arrays at index 1
    var k = k - 1; // same as above
	while (true) {
		if (r - l > cutoff &&  jstack < nstack) {
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
			if (i == m/2) {
				sd = 0;
			}
			stack_1[jstack] = l;
			stack_2[jstack] = r;
			jstack++;
			var comp = k - i*(s/dm) + sd;
			if (l < comp) {
				l = Math.floor(comp + 0.5); // l is a positive integer
			}
			if (r > comp + s) {
				r = Math.floor(comp + s + 0.5); // r is a positive integer
			}
		}
		else {
			if (l >= r) {
				if (jstack == 0) {
					return x[k];
				}
				--jstack;
				l = stack_1[jstack];
				r = stack_2[jstack];
			}
			var v = x[k];
			i = l;
			j = r;
			x[k] = x[l];
			x[l] = v;
			if (compareFunction(v, x[r]) < 0) {
				x[l] = x[r];
				x[r] = v;
			}
			
			while (i < j) {
				var tmp = x[j];
	            x[j] = x[i];
	            x[i] = tmp;
				
				++i;
				--j;
				while (compareFunction(x[i], v) < 0) {
					++i;
				}
				while (compareFunction(x[j], v) > 0) {
					--j;
				}
			}
			if (compareFunction(x[l], v) == 0) {
				var tmp = x[l];
				x[l] = x[j];
				x[j] = tmp;
			} 
			else {
				++j;
				var tmp = x[j];
				x[j] = x[r];
				x[r] = tmp;
			}
			if (j <= k) {
				l = j + 1;
			}
			if (k <= j) {
				r = j - 1;
			}
		}
	}
}
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
function hypot_(x, y) {
	var r = 0;
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
	return r;
}
function rank_(x, order) {
	var xWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		xWithIndexes[i] = [x[i], i];
	}
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
	return xRanks;
}
function ftca_(correlationMatrix, threshold) {
	var threshold = threshold;
	if (threshold === undefined) {
		threshold = 0.5;
	}
	var correlationMatrix = new Matrix_(correlationMatrix);
	var clusters = [];
	var nbElements = correlationMatrix.nbRows;
	var unassignedElementsIdx = new Array(nbElements);
	for (var i = 0; i < unassignedElementsIdx.length; ++i) {
		unassignedElementsIdx[i] = i + 1;
	}
	while (unassignedElementsIdx.length != 0) {
		if (unassignedElementsIdx.length === 1) {
			var newCluster = [unassignedElementsIdx[0]];
			unassignedElementsIdx[0] = null;
			clusters.push(newCluster);
		}	
		else {
			var subCorrMat = correlationMatrix.submatrix(unassignedElementsIdx, unassignedElementsIdx);
			var subCorrMatRows = subCorrMat.toRowArray(function(i, j, val) {
				return i != j;
			});
			var avgCorrelation = new Array(unassignedElementsIdx);
			for (var i = 0; i < unassignedElementsIdx.length; ++i) {
				avgCorrelation[i] = mean_(subCorrMatRows[i]);
			}
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
			if (correlationMatrix.getValueAt(hc, lc) > threshold) {
				var newClusterHcLc = (hc === lc ? [hc] : [hc, lc]);
				unassignedElementsIdx[hcIdx] = null;
				unassignedElementsIdx[lcIdx] = null;
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned elements (HC and LC)
						var avgHcLcAssetCorrelation = (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) + correlationMatrix.getValueAt(unassignedElementsIdx[i], lc)) / 2;
						if (avgHcLcAssetCorrelation  > threshold) {
							newClusterHcLc.push(unassignedElementsIdx[i]);
							unassignedElementsIdx[i] = null;
						}
					}
				}
				clusters.push(newClusterHcLc);
			}
			else {
				var newClusterHc = [hc];
				unassignedElementsIdx[hcIdx] = null;
				for (var i = 0; i < unassignedElementsIdx.length; ++i) {
					if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC)
						if (correlationMatrix.getValueAt(unassignedElementsIdx[i], hc) > threshold) {
							newClusterHc.push(unassignedElementsIdx[i]);
							unassignedElementsIdx[i] = null;
						}
					}
				}
				clusters.push(newClusterHc);
				if (hc !== lc) {
					var newClusterLc = [lc];
					unassignedElementsIdx[lcIdx] = null;
					for (var i = 0; i < unassignedElementsIdx.length; ++i) {
						if (unassignedElementsIdx[i] !== null) { // Skip assigned assets (HC with its correlated assets, and LC)
							if (correlationMatrix.getValueAt(unassignedElementsIdx[i], lc) > threshold) {
								newClusterLc.push(unassignedElementsIdx[i]);
								unassignedElementsIdx[i] = null;
							}
						}
					}
					clusters.push(newClusterLc);
				}
			}
		}
		var newUnassignedElementsIdx = [];
		for (var i = 0; i < unassignedElementsIdx.length; ++i) {
			if (unassignedElementsIdx[i] !== null) {
				newUnassignedElementsIdx.push(unassignedElementsIdx[i]);
			}
		}
		unassignedElementsIdx = newUnassignedElementsIdx;
	}
	return clusters;
}
function mean_(x) {
	var nn = x.length;
	var tmpMean = 0.0;
	var sum = 0.0;
	for (var i=0; i<nn; ++i) {
		sum += x[i];
	}
	tmpMean = sum/nn;
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		sumDiff += (x[i] - tmpMean);
	}
	return (sum + sumDiff)/nn;
}
function variance_(x) {
	var nn = x.length;
	var meanX = mean_(x);
	var sumSquareDiff = 0.0;
	var sumDiff = 0.0;
	for (var i=0; i<nn; ++i) {
		var diff = (x[i] - meanX);
		sumSquareDiff += diff * diff;
		sumDiff += diff;
	}
	var S = sumSquareDiff - ((sumDiff * sumDiff) / nn);
	return S/nn;
}
function sampleVariance_(x) {
	var nn = x.length;
	return variance_(x) * nn/(nn - 1);
}
function stddev_(x) {
	return Math.sqrt(variance_(x));
}
function sampleStddev_(x) {
	return Math.sqrt(sampleVariance_(x));
}
function normcdf_(x) {
	var sum = x;
	var term = 0;
	var next_term = x;
	var power = x*x;
	var i = 1;
	if (x < -8.0) {
		return 0.0;
	}
	else if (x > 8.0) {
		return 1.0;
	}
	while (sum != term) {
		sum = (term = sum) + (next_term *= power/(i += 2));
	}
	return 0.5 + sum * Math.exp(-0.5 * power - 0.91893853320467274178)
}
function normrnd_(mu, sigma) {
	if (mu == undefined) {
		mu = 0;
	}
	if (sigma == undefined) {
		sigma = 1;
	}
	var u = Math.random(); // u ~ U[0,1[
	while (u === 0.0) {
		u = Math.random();
	} // u ~ U]0,1[
	var r = norminv_(u); // r ~ N(0,1)
	return mu + sigma*r;
}
function pnormrnd_(mu, sigma) {
	if (mu == undefined) {
		mu = 0;
	}
	if (sigma == undefined) {
		sigma = 1;
	}
	var sigma_sq = sigma*sigma;
	var A = 1.136717791056118;
	var mu_a = (1 - A*A)/A*sigma;
	var mu_b = 0;
	var mu_c = sigma * Math.sqrt(Math.PI/2);
	var z;
	while (true) {
		var rho;
		if (mu < mu_a) {
			var a = (-mu + Math.sqrt(mu*mu + 4*sigma_sq)) / 2 / sigma_sq;
			
			z = -Math.log(1 - Math.random())/a;
			rho = Math.exp( -(z-mu)*(z-mu) / 2 / sigma_sq - a*(mu-z+a*sigma_sq/2) );
		}
		else if (mu <= mu_b) { 
			z = Math.abs(normrnd_())*sigma + mu;
			rho = (z >= 0) ? 1 : 0;
		}
		else if (mu < mu_c) { 
			var r = (Math.random() < mu/(mu + Math.sqrt(Math.PI/2)*sigma)) ? 1 : 0;
			var u = Math.random()*mu;
			var g = Math.abs(normrnd_()*sigma) + mu;
			
			z = r*u + (1-r)*g;
			rho = r*Math.exp(-(z-mu)*(z-mu) / 2 /sigma_sq) + (1-r);
		}
		else { 
			z = normrnd_()*sigma + mu;
			rho = (z >= 0) ? 1 : 0;
		}
		if (Math.random() <= rho){
			break;
		}		
	}
	return z;
}
function norminv_(p) {
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
    return ppnd16;
}
function hypersphereRandomSampler_(n, reuseOutputArray) {
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.reuseOutputArray = reuseOutputArray;
	this.sample = function() {
		var t = 0;
		var s = 1;
		for (var i = 0; i < this.n; ++i) {
			var r = normrnd_();
			this.x[i] = r;
			var absR = Math.abs(r);
			if (absR != 0) {
				if (absR > t) {
					s = 1 + s * (t/r) * (t/r);
					t = absR;
				}
				else  {
					s = s + (r/t) * (r/t);
				}
			}
		}
		var x_two_norm = t * Math.sqrt(s);
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/x_two_norm;
		}
		if (this.reuseOutputArray) {
			return this.x;
		}
		else {
			return this.x.slice(0);
		}
	}
}
function covariance_(x, y) {
	var nn = x.length;
	var meanX = mean_(x);
	var meanY = mean_(y);
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
	var C = sumProdDiff - ((sumDiffX * sumDiffY) / nn);
	return C/nn;
}
function sampleCovariance_(x, y) {
	var nn = x.length;
	return covariance_(x,y) * nn/(nn - 1);
}
function thresholdAcceptingSolve_(f, x0, opt) { 
	function defaultNeighbourGenerator(x, neighbourGeneratorParameters) {
		var alpha = neighbourGeneratorParameters.alpha;
		var l = neighbourGeneratorParameters.lowerBounds;
		var u = neighbourGeneratorParameters.upperBounds;
		var n = x.length;
		var xl = x.slice();
		var xu = x.slice();
		for (var i = 0; i < n; ++i) {
			xl[i] = x[i] - alpha/2;
			xu[i] = x[i] + alpha/2;
			if (l) {
				xl[i] = Math.max(xl[i], l[i]);
			}
			if (u) {
				xu[i] = Math.min(xu[i], u[i]);
			}
		}
		var xx = new boxRandomSampler_(n, xl, xu).sample();
		return xx;
	}
	function computeThresholds(f, x_c, neighbourGeneratorParameters, nDeltas, nRounds) {	
		var f_x_c = f(x_c);
		var deltas = typeof Float64Array === 'function' ? new Float64Array(nDeltas) : new Array(nDeltas);
		for (var i = 0; i < nDeltas; ++i) {
			var x_i = neighbourGenerator(x_c, neighbourGeneratorParameters);
			var f_x_i = f(x_i);
			deltas[i] = Math.abs(f_x_c - f_x_i);
			x_c = x_i;
			f_x_c = f_x_i;
		}
		deltas.sort(function(a, b) { return a-b });
		var taus = typeof Float64Array === 'function' ? new Float64Array(nRounds) : new Array(nRounds);
		for (var i = 0; i < nRounds; ++i) {
			taus[i] = quantile_(deltas, (nRounds - (i+1))/nRounds, true);
		}
		taus[nRounds-1] = 0;
		return taus;
	}
	var n = x0.nbRows;
	if (opt === undefined) {
		opt = {};
	}
	var nRounds = opt.nRounds;
	if (nRounds === undefined) {
		nRounds = 10;
	}
	var nSteps = opt.nSteps;
	if (nSteps === undefined) {
		nSteps = 5000;
	}
	var nDeltas = opt.nDeltas;
	if (nDeltas === undefined) {
		nDeltas = nSteps;
	}
	var neighbourGenerator = opt.neighbourGenerator;
	var neighbourGeneratorParameters = opt.neighbourGeneratorParameters;
	if (neighbourGenerator === undefined) {
		neighbourGenerator = defaultNeighbourGenerator;
		
		if (neighbourGeneratorParameters === undefined || 
		    neighbourGeneratorParameters.alpha === undefined) {
			neighbourGeneratorParameters = { alpha: 1e-3 };
		}
	}
	var x_c = x0.slice();
	var taus = computeThresholds(f, x_c, neighbourGeneratorParameters, nDeltas, nRounds);
	var x_best = x0.slice();
	var f_x_best = f(x_best);
	var x_c = x0.slice();
	var f_x_c = f(x_c);
	for (var r = 0; r < nRounds; ++r) {
		for (var i = 0; i < nSteps; ++i) {
			var x_c_tmp = x_c.slice();
			var x_ri = neighbourGenerator(x_c, neighbourGeneratorParameters);
			var f_x_ri = f(x_ri);
			var delta = f_x_ri - f_x_c;
			if (delta < taus[r]) {
				x_c = x_ri.slice();
				f_x_c = f_x_ri;
			}
			else {
				x_c = x_c_tmp;
			}
			if (f_x_ri < f_x_best) {
				f_x_best = f_x_ri;
				x_best = x_ri.slice();
			}
		}
	}
	return [x_best, f_x_best];
}
 function gssSolve_(f, x0, l, u, opt) {
	function coordinateDirectionsPollingSetGenerator(x, alpha, Ip, Im) {
		this.n = x.nbRows;
		this.alpha = alpha;
		this.nbGeneratedPollingDirections = 0;
		this.d = Matrix_.zeros(n, 1); // placeholder for the coordinate direction to generate at each call of the .next() method
		this.nbPollingDirections;
		this.pollingDirectionsIndices;
		if (Ip && Im) {	
			this.nbPollingDirections = Ip.length + Im.length;
			this.pollingDirectionsIndices = typeof Int32Array === 'function' ? new Int32Array(this.nbPollingDirections) : new Array(this.nbPollingDirections);
			for (var i = 0; i < Ip.length; ++i) {
				this.pollingDirectionsIndices[i] = Ip[i];
			}
			for (var i = Ip.length, j= 0; i < this.nbPollingDirections; ++i, ++j) {
				this.pollingDirectionsIndices[i] = -Im[j];
			}
		}
		else {
			this.nbPollingDirections = 2*this.n;
			this.pollingDirectionsIndices = typeof Int32Array === 'function' ? new Int32Array(this.nbPollingDirections) : new Array(this.nbPollingDirections);
			for (var i = 0; i < this.n; ++i) {
				this.pollingDirectionsIndices[i] = i + 1;
			}
			for (var i = this.n, j = 0; i < 2*this.n; ++i, ++j) {
				this.pollingDirectionsIndices[i] = - (j + 1);
			}
		}
		this.next = function() {		
			if (this.nbGeneratedPollingDirections >= this.nbPollingDirections) {
				return -1;
			}
			if (this.nbGeneratedPollingDirections >= 1) {
				var previousPollingDirectionIndice = this.pollingDirectionsIndices[this.nbGeneratedPollingDirections - 1];
				if (previousPollingDirectionIndice < 0) {
					this.d.data[-previousPollingDirectionIndice - 1] = 0;
				}
				else if (previousPollingDirectionIndice > 0) {
					this.d.data[previousPollingDirectionIndice - 1] = 0;
				}
				else {
					throw new Error('internal error: 0 polling direction indice detected');
				}
			}
			var pollingDirectionIndice = this.pollingDirectionsIndices[this.nbGeneratedPollingDirections];
			if (pollingDirectionIndice < 0) {
				this.d.data[-pollingDirectionIndice - 1] = -1;
			}
			else if (pollingDirectionIndice > 0) {
				this.d.data[pollingDirectionIndice - 1] = 1;
			}
			else {
				throw new Error('internal error: 0 polling direction indice detected');
			}
			++this.nbGeneratedPollingDirections;
			return this.d;
		}
	}
	function randomUnitSpherePollingSetGenerator(x, alpha) {
		this.n = x.nbRows;
		this.nbGeneratedPollingDirections = 0;
		this.d = Matrix_.zeros(n, 1); // placeholder for the direction to generate at each call	
		this.nbPollingDirections = Math.floor(Math.log( 1 - Math.log(theta)/Math.log(gamma) ) / Math.log(2)) + 1; // = m in section 5.4 of the second reference, the minimal number of random directions to generate
		this.hypersphereRandomSampler = new hypersphereRandomSampler_(this.n, true); // uniform sampler of directions on the unit sphere, with output array re-usage for improved performances
		this.next = function() {
			if (this.nbGeneratedPollingDirections >= this.nbPollingDirections) {
				return -1;
			}
			if (this.nbPollingDirections === 2 && this.nbGeneratedPollingDirections === 1) {
				this.d = Matrix_.ax(-1, this.d, this.d);
			}
			else {
				var vect = this.hypersphereRandomSampler.sample();
				this.d = Matrix_.fill(this.n, 1, function(i,j) { return vect[i-1]; }, this.d);
			}
			++this.nbGeneratedPollingDirections;
			return this.d;
		}
	}
	var n = x0.nbRows;
	var eps_tol = 1e-12; // used to numerically determine some conditions
	if (opt === undefined) {
		opt = {};
	}
	var maxIterations = opt.maxIter;
	if (maxIterations === undefined) {
		maxIterations = 10000;
	}
	var eps = opt.eps;
	if (eps === undefined) {
		eps = 1e-6;
	}
	if (l && u) {
		if (l.length != u.length) {
			throw new Error("incompatible number of lower bounds and upper bounds constraints: " + l.length + " v.s. " + u.length);
		}
		else {
			l = new Matrix_(l);
			u = new Matrix_(u);
		}
	}
	if (l && !u) {
		throw new Error('missing upper bounds constraints');
	}
	else if (!l && u) {
		throw new Error('missing lower bounds constraints');
	}
	var unconstrainedPollingSet = opt.unconstrainedPollingSet;
	if (unconstrainedPollingSet === undefined) {
		unconstrainedPollingSet = 'coordinateDirections';
	}
	if (unconstrainedPollingSet !== 'coordinateDirections' 
	    && unconstrainedPollingSet !== 'probabilisticDescentDirections'
		&& unconstrainedPollingSet !== 'custom') {
		throw new Error('unsupported unconstrained polling set');
	}
	var unconstrainedPollingSetGenerator;
	if (unconstrainedPollingSet === 'coordinateDirections') {
		unconstrainedPollingSetGenerator = coordinateDirectionsPollingSetGenerator;
	}
	else if (unconstrainedPollingSet === 'probabilisticDescentDirections') {
		unconstrainedPollingSetGenerator = randomUnitSpherePollingSetGenerator;
	}
	else if (unconstrainedPollingSet === 'custom') {
		unconstrainedPollingSetGenerator = opt.customUnconstrainedPollingSetGenerator;
	}
	else {
		throw new Error('internal error: unsupported unconstrained polling set detected');
	}
	var constrainedPollingSet = opt.constrainedPollingSet;
	if (constrainedPollingSet === undefined) {
		constrainedPollingSet = 'coordinateDirections';
	}
	if (constrainedPollingSet !== 'coordinateDirections' 
		&& constrainedPollingSet !== 'custom') {
		throw new Error('unsupported constrained polling set');
	}
	var constrainedPollingSetGenerator;
	if (constrainedPollingSet === 'coordinateDirections') {
		constrainedPollingSetGenerator = coordinateDirectionsPollingSetGenerator;
	}
	else if (constrainedPollingSet === 'custom') {
		constrainedPollingSetGenerator = opt.customConstraintedPollingSetGenerator;
	}
	else {
		throw new Error('internal error: unsupported constrained polling set detected');
	}
	var alphaZero = opt.alphaZero;
	if (alphaZero === undefined) {
		alphaZero = 1;
	}
	var alphaMax = opt.alphaMax;
	if (alphaMax === undefined) {
		alphaMax = 1e10;
	}
	var gamma = opt.gamma;
	if (gamma === undefined) {
		if (unconstrainedPollingSet === 'coordinateDirections') {
			gamma = 1;
		}
		else if (unconstrainedPollingSet === 'probabilisticDescentDirections') {
			gamma = 2;
		}
		else if (unconstrainedPollingSet === 'custom') {
			gamma = 1;
		}
		else {
			throw new Error('internal error: unsupported unconstrained polling set detected');
		}
	}
	var theta = opt.theta;
	if (theta === undefined) {
		theta = 0.5;
	}
	var pollingType = opt.pollingType;
	if (pollingType === undefined) {
		pollingType = 'opportunistic';
	}
	if (pollingType !== 'opportunistic' && pollingType !== 'complete') {
		throw new Error('unsupported polling type');
	}
    var rho = opt.rho; 
	if (rho === undefined) {
		rho = function (alpha) {
			return alpha * alpha / 2;
		}
	}
	var x_k = new Matrix_(x0); // placeholder for the x_k vector
	var x_kp = new Matrix_(x0); // placeholder for the x_kp vector
	var x_kp_best = new Matrix_(x0); // placeholder for the best x_kp vector, mostly used in case of complete polling
	var f_x_k = f(x_k); // placeholder for the f(x_k) value
	var f_x_kp = f(x_kp); // placeholder for the f(x_kp) value
	var f_x_kp_best = undefined; // placeholder for the best f(x_kp) value, mostly used in case of complete polling
	var alpha_k = alphaZero;
	var iter = 0;
	while (true) {
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		var rho_alpha_k = rho(alpha_k);
		var searchSuccessful = false;
		var pollSuccessful = false;

		if (searchSuccessful === false) {		
			if (l && u) {
				var Ip = new Array(0);
				var Im = new Array(0);
				
				for (var i = 0; i < n; ++i) {
					if (x_k.data[i] + alpha_k <= u.data[i] - eps_tol) {
						Ip.push(i + 1);
					}
					
					if (l.data[i] + eps_tol <= x_k.data[i] - alpha_k) {
						Im.push(i + 1);
					}
				}
			}
			var D_k;
			if ((l && u) && Ip.length + Im.length < 2*n) {
				D_k = new constrainedPollingSetGenerator(x_k, alpha_k, Ip, Im);
			}
			else {
				D_k = new unconstrainedPollingSetGenerator(x_k, alpha_k);
			}		
			var d_k = D_k.next();
			while (d_k != -1) {
				x_kp = Matrix_.axpby(1, x_k, alpha_k, d_k, x_kp);
				f_x_kp = f(x_kp);
				
				if ( f_x_kp < f_x_k - rho_alpha_k ) {
					pollSuccessful = true;
					if (pollingType === 'opportunistic') {
						x_kp_best = Matrix_.copy(x_kp, x_kp_best);
						f_x_kp_best = f_x_kp;
						
						break;
					}
					else if (pollingType === 'complete') {
						if (f_x_kp_best === undefined) {
							x_kp_best = Matrix_.copy(x_kp, x_kp_best);
							f_x_kp_best = f_x_kp;						
						}
						else {
							if (f_x_kp < f_x_kp_best) {
								x_kp_best = Matrix_.copy(x_kp, x_kp_best);
								f_x_kp_best = f_x_kp;
							}
						}
					}
				}
				
				d_k = D_k.next(); 
			}
		}
		if (searchSuccessful === true || pollSuccessful === true) {
			x_k = Matrix_.copy(x_kp_best, x_k);
			f_x_k = f_x_kp_best;
			f_x_kp_best = undefined;
		}
		if (pollSuccessful === true || searchSuccessful === true) {
			alpha_k = Math.min(gamma * alpha_k, alphaMax);
		}
		else {
			alpha_k = theta * alpha_k;			
		}
		if (alpha_k <= eps) {
			break;
		}
	}
	return [x_k, f_x_k];
}
function goldenSectionSearch_(f, x_min, x_max, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var maxIterations = opt.maxIter;
	if (maxIterations === undefined) {
		maxIterations = -1;
	}
	var eps = opt.eps;
	if (eps === undefined) {
		eps = 1e-6;
	}
	var inv_phi = (Math.sqrt(5) - 1) / 2; // ~0.618
	
	var a = x_min;
	var b = x_max;
	var c = b - inv_phi * (b - a);
	var d = a + inv_phi * (b - a);
	
	var f_x_min = f(x_min);
	var f_x_max = f(x_max);
	var f_c = f(c);
	var f_d = f(d);
	if (x_min > x_max) {
		throw new Error('bracketing interval lower bound ' + x_min + 
		                ' greater than bracketing interval upper bound ' + x_max);
	}
	var iter = 0;	
	while (true) {
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		if (f_c <= f_d) { // x^* belongs to [a, d]
			b = d;
			d = c;
			c = b - inv_phi * (b - a);
	
			f_d = f_c;
			f_c = f(c);
		}
		else if (f_c > f_d) { // x^* belongs to [c, b]
			a = c;
			c = d;
			d = a + inv_phi * (b - a);
			
			f_c = f_d;
			f_d = f(d);
		}
		if (Math.abs(b - a) <= eps ) {
			if (f_c < f_d) {
				if (f_x_min < f_c) {
					return [x_min, f_x_min];
				}
				else {
					return [c, f_c];
				}
			}
			else {
				if (f_x_max < f_d) {
					return [x_max, f_x_max];
				}
				else {
					return [d, f_d];
				}
			}
		}
	}
}
function bisection_(f, x_min, x_max, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var maxIterations = opt.maxIter;
	if (maxIterations === undefined) {
		maxIterations = 45;
	}
	var eps = opt.eps;
	if (eps === undefined) {
		eps = 1e-6;
	}
	var outputInterval = opt.outputInterval;
	if (outputInterval === undefined) {
		outputInterval = false;
	}
	var f_x_min = f(x_min);
	var f_x_max = f(x_max);
	var r; // a root of the function f in the interval [x_min, x_max]
	var dx; // the signed length of the interval [x_min, x_max], which will be halved in each iteration of the bisection algorithm
	if (f_x_min <= 0) { // by convention taken from the second reference, the bisection root search is oriented so that f(r) <= 0 (i.e., f > 0 is at r + dx)
		r = x_min;
		dx = x_max - x_min;
	}
	else {
		r = x_max;
		dx = x_min - x_max;
	}
	if (x_min >= x_max) {
		throw new Error('bracketing interval lower bound ' + x_min + 
		                ' greater than bracketing interval upper bound ' + x_max);
	}
	if (f_x_min == 0) {
		if (outputInterval) {
			return [x_min, x_min + eps];
		}
		else {
			return x_min; // a root has been found !
		}
	}
	if (f_x_max == 0) {
		if (outputInterval) {
			return [x_max - eps, x_max];
		}
		else {
			return x_max; // a root has been found !
		}
	}
	if (f_x_min*f_x_max > 0.0) {
		throw new Error('interval [' + x_min + ',' + x_max + '] is not a bracketing interval');
	}
	var iter = 0;	
	while (true) {
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		dx = dx / 2;
		var x_mid = r + dx;
		var f_mid = f(x_mid);
		if (f_mid <= 0) {
			r = x_mid;
		}
		if (Math.abs(dx) <= eps || f_mid == 0) {
			if (outputInterval) {
				return [Math.max(x_min, r - eps), Math.min(x_max, r + eps)];
			}
			else {
				return r;
			}
		}
	}
}
function ccpsolveFISTA_(f, gradf, g, proxg, x0, opt) {
	function F(x) {
		return f(x) + g(x);
	}
	function Q(mu, u, v, gradf_v) {
		var f_v = f(v);
		var u_m_v = Matrix_.xmy(u, v);
		var u_m_v_two_norm = u_m_v.vectorNorm('two');
		var g_u = g(u);
		var Q_mu_u_v = f_v + Matrix_.vectorDotProduct(u_m_v, gradf_v) + 1/(2 * mu) * u_m_v_two_norm * u_m_v_two_norm + g_u;
		return Q_mu_u_v;
	}
	function p(mu, v, gradf_v) {
		var v_m_mu_gradf_v = Matrix_.axpby(1, v, -mu, gradf_v);
		var p_mu_v = Matrix_.copy(proxg(v_m_mu_gradf_v, mu));
		return [v_m_mu_gradf_v, p_mu_v];
	}
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-04;
	var maxIterations = opt.maxIter || 10000;
	var maxLineSearches = opt.maxLine || 100;
	var beta = opt.beta || 0.5;
	var alphaMin = opt.alphaMin || 1e-10;
	var alphaMax = opt.alphaMax || 1e10;
	var restartPeriod = opt.restartPeriod || 1000;
	var n = x0.nbRows;
	var eps_tol = 1e-12; // used to numerically determine some conditions (backtrack, adaptative restart, stepsize)
	var t_km; 
	var t_k;
	var theta_km;
	var theta_k;
	var tau_k = 0.5;
	var mu_k_0 = alphaMin;
	var x_k; // placeholder for the x_k vector
	var x_km; // placeholder for the x_k-1 vector
	var x_kmm; // placeholder for the x_k-2 vector
	var x_km_m_x_kmm; //  placeholder for the x_k-1 - x_k-2 vector
	var gradf_x_k = Matrix_.zeros(n, 1); // placeholder for the gradf(x_k) vector
	var gradf_x_km = Matrix_.zeros(n, 1);; // placeholder for the gradf(x_k-1) vector
	var gradf_x_kmm = Matrix_.zeros(n, 1);; // placeholder for the gradf(x_k-2) vector
	var gradf_x_km_m_gradf_x_kmm; // // placeholder for the gradf(x_k-1) - gradf(x_k-2) vector
	var y_k = Matrix_.zeros(n, 1); // placeholder for the y_k vector
	var gradf_y_k = Matrix_.zeros(n, 1); // placeholder for the gradf(y_k) vector	
	var restart = true; // a first initialization is needed at the first iteration 

	var iter = 0;	
	while (true) {
		if (maxIterations !== -1 && iter > maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		if (iter % restartPeriod === 0) {
			x0 = x_k;
			restart = true;
		}
		if (restart === true) {
			t_km = 0; 
			t_k = 1;
			theta_km = 1;
			theta_k = null;
			x_k = new Matrix_(x0);
			x_km = new Matrix_(x_k);
			x_kmm = new Matrix_(x_km);
			x_km_m_x_kmm = Matrix_.zeros(n, 1);
			gradf_x_k = Matrix_.copy(gradf(x_k), gradf_x_k);
			gradf_x_km = Matrix_.copy(gradf_x_k, gradf_x_km);
			gradf_x_kmm = Matrix_.copy(gradf_x_km, gradf_x_kmm);
			gradf_x_km_m_gradf_x_kmm = Matrix_.zeros(n, 1);
			y_k = Matrix_.copy(x_k, y_k); 
			gradf_y_k = Matrix_.copy(gradf_x_k, gradf_y_k);
			restart = false;
		}
		var mu_k = mu_k_0;
		var p_mu = p(mu_k, y_k, gradf_y_k);
		var y_k_m_mu_k_gradf_y_k = p_mu[0];
		var p_mu_k_y_k = p_mu[1];
		
		var iter_ls = 0;
		while ( F(p_mu_k_y_k) > Q(mu_k, p_mu_k_y_k, y_k, gradf_y_k) + eps_tol ) {
			if (maxLineSearches !== -1 && iter_ls > maxLineSearches) {
				throw new Error('maximum number of line searches reached: ' + maxLineSearches + ' at iteration: ' + iter);
			}
			++iter_ls;
			mu_k = beta * mu_k;
			theta_km = theta_km/beta;
			t_k = ( 1 + Math.sqrt(1 + 4*theta_km*t_km*t_km) ) / 2;
			y_k = Matrix_.axpby(1, x_km, (t_km - 1)/t_k, x_km_m_x_kmm, y_k);
			gradf_y_k = Matrix_.axpby(1, gradf_x_km, (t_km - 1)/t_k, gradf_x_km_m_gradf_x_kmm, gradf_y_k);
			p_mu = p(mu_k, y_k, gradf_y_k);
			y_k_m_mu_k_gradf_y_k = p_mu[0];
			p_mu_k_y_k = p_mu[1];
		}
		x_k = p_mu_k_y_k;
		gradf_x_k = Matrix_.copy(gradf(x_k), gradf_x_k); // gradf(x_k)
		
		var x_k_m_x_km = Matrix_.xmy(x_k, x_km); // x_k - x_k-1
		var gradf_x_k_m_gradf_x_km = Matrix_.xmy(gradf_x_k, gradf_x_km); // gradf(x_k) - gradf(x_k-1)
		var s_k_d_s_k = Matrix_.vectorDotProduct(x_k_m_x_km, x_k_m_x_km); // <x_k - x_k-1/x_k - x_k-1>
		var z_k_d_z_k = Math.max(Matrix_.vectorDotProduct(gradf_x_k_m_gradf_x_km, gradf_x_k_m_gradf_x_km), // <gradf(x_k) - gradf(x_k-1)/gradf(x_k) - gradf(x_k-1)>
		                         eps_tol); // to avoid numerical issues in the division below
		var s_k_d_z_k = Matrix_.vectorDotProduct(x_k_m_x_km, gradf_x_k_m_gradf_x_km); // <x_k - x_k-1/gradf(x_k) - gradf(x_k-1)>

		if (s_k_d_z_k <= eps_tol) {
			mu_kp_0 = alphaMax;
		}
		else {
			var alpha_k_1 = Math.max(alphaMin, Math.min(s_k_d_s_k/s_k_d_z_k, alphaMax));
			var alpha_k_2 = Math.max(alphaMin, Math.min(s_k_d_z_k/z_k_d_z_k, alphaMax));
			
			if (alpha_k_2/alpha_k_1 <= tau_k) {
				mu_kp_0 = alpha_k_2;
				tau_k = 0.9 * tau_k;
			}
			else {
				mu_kp_0 = alpha_k_1;
				tau_k = 1.1 * tau_k;
			}
		}
		theta_k = mu_k/mu_kp_0;
		var t_kp = ( 1 + Math.sqrt(1 + 4*theta_k*t_k*t_k) ) / 2; // t_k+1
		var y_kp = Matrix_.axpby(1, x_k, (t_k - 1)/t_kp, x_k_m_x_km); // y_k+1	
		
		var gradf_y_kp = Matrix_.axpby(1, gradf_x_k, (t_k - 1)/t_kp, gradf_x_k_m_gradf_x_km); // gradf(y_k+1)
		var subgradg_x_k = Matrix_.axpby(1/mu_k, y_k_m_mu_k_gradf_y_k, -1/mu_k, x_k);
		var r_k = Matrix_.xpy(gradf_x_k, subgradg_x_k);
		if (r_k.vectorNorm('infinity') <= eps) {
			break;
		}
		var gs_k = Matrix_.vectorDotProduct(Matrix_.xmy(y_k, x_k), x_k_m_x_km);
		if (gs_k >= eps_tol) {
			x0 = x_k;
			restart = true;
		}
		mu_k_0 = mu_kp_0;
		x_kmm = x_km;
		x_km = x_k;
		x_km_m_x_kmm = x_k_m_x_km;
		gradf_x_kmm = Matrix_.copy(gradf_x_km, gradf_x_kmm);
		gradf_x_km = Matrix_.copy(gradf_x_k, gradf_x_km);
		gradf_x_km_m_gradf_x_kmm = gradf_x_k_m_gradf_x_km;
		y_k = Matrix_.copy(y_kp, y_k);
		gradf_y_k = Matrix_.copy(gradf_y_kp, gradf_y_k);
		theta_km = theta_k;
		t_km = t_k;
		t_k = t_kp;
	}
	return [x_k, F(x_k)];
}
function qksolveBS_(d, a, b, r, l, u, opt) {
	function resizeArray(arr, n) {
		if (arr instanceof Array) { // this restrict the size of the array to the first n elements
			arr.length = n; 
			return arr;
		}
		else if (arr instanceof Float64Array || arr instanceof Int32Array) { // this constructs a view on the first n array elements
			return arr.subarray(0, n);
		}
	}
	function g(t) {
		var g_t = (p - t * q) + s;
		var I_it = new I.iterator();
		var el = I_it.next();
		while (el != 0) {
			var i = el - 1;

			
			var x_i;
			if (t <= T_u[i]) {
				x_i = u.data[i];
			}
			else if (T_u[i] <= t && t <= T_l[i]) {
				x_i = (a.data[i] - t*b.data[i]) / d.data[i];
			}
			else if (T_l[i] <= t) {
				x_i = l.data[i];
			}
			
			g_t += b.data[i] * x_i;
			el = I_it.next();
		}
		return g_t;
	}
	function x(t) {
		var x_t = Matrix_.zeros(n, 1);
		for (var i = 0; i < n; ++i) {
			var x_i;
			if (t_star <= T_u[i]) {
				x_i = u.data[i];
			}
			else if (T_u[i] <= t_star && t_star <= T_l[i]) {
				x_i = (a.data[i] - t_star*b.data[i]) / d.data[i];
			}
			else if (T_l[i] <= t_star) {
				x_i = l.data[i];
			}
			x_t.data[i] = x_i;
		}
		return x_t;
	}
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-16;
	var outputLagrangeMultiplier = false || opt.outputLagrangeMultiplier;
	if (!(d instanceof Matrix_) || !d.isVector()) {
		throw new Error('first input must be a vector');
	}
	if (!(a instanceof Matrix_) || !a.isVector()) {
		throw new Error('second input must be a vector');
	}
	if (!(b instanceof Matrix_) || !b.isVector()) {
		throw new Error('third input must be a vector');
	}
	if (!(l instanceof Matrix_) || !l.isVector()) {
		throw new Error('fifth input must be a vector');
	}
	if (!(u instanceof Matrix_) || !u.isVector()) {
		throw new Error('sixth input must be a vector');
	}
	
	if (d.nbRows !== a.nbRows) {
		throw new Error('first and second inputs number of rows do not match: ' + d.nbRows + '-' + a.nbRows);
	}
	if (d.nbRows !== b.nbRows) {
		throw new Error('first and third inputs number of rows do not match: ' + d.nbRows + '-' + b.nbRows);
	}
	if (d.nbRows !== l.nbRows) {
		throw new Error('first and fifth inputs number of rows do not match: ' + d.nbRows + '-' + l.nbRows);
	}
	if (d.nbRows !== u.nbRows) {
		throw new Error('first and sixth inputs number of rows do not match: ' + d.nbRows + '-' + u.nbRows);
	}
	var n = b.nbRows;
	
	var abs_r = Math.abs(r);
	
	var T = typeof Float64Array === 'function' ? new Float64Array(2*n) : new Array(2*n); // the set of breakpoints T
	var T_l = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the list of breakpoints t_l_i, i=1..n
	var T_u = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the list of breakpoints t_u_i, i=1...n
	
	var I = new BitSet_().setRange(1, n); // the set of indices I
	var p = 0;
	var q = 0;
	var s = 0;
	var t_1 = Infinity;
	var t_r = -Infinity;
	for (var i = 0, j = 0; i < n; ++i) {
		if (l.data[i] > u.data[i]) {
			throw new Error('infeasible problem detected');
		}
		if (b.data[i] <= 0) {
			throw new Error('negative element detected in b');
		}
		if (d.data[i] <= 0) {
			throw new Error('negative element detected in d');
		}		
		var t_l_i = (a.data[i] - l.data[i]*d.data[i]) / b.data[i];
		T_l[i] = t_l_i;
		T[j++] = t_l_i;
		var t_u_i = (a.data[i] - u.data[i]*d.data[i]) / b.data[i];
		T_u[i] = t_u_i;
		T[j++] = t_u_i;
		if (t_l_i > t_r) {
			t_r = t_l_i;
		}
		if (t_u_i < t_1) {
			t_1 = t_u_i;
		}	
	}
	var g_t_1 = g(t_1);
	var g_t_r = g(t_r);
	if (g_t_1 < r || g_t_r > r) {
		throw new Error('infeasible problem detected');
	}
	var t_star = null;
	if (Math.abs(g_t_1 - r) <= eps * abs_r) {
		t_star = t_1;
	}
	else if (Math.abs(g_t_1 - r) <= eps * abs_r) {
		t_star = t_r;
	}
	else {
		var t_l = t_1; // t_1 was already computed to check feasibility, so there is no need to use t_0
		var t_u = t_r; // t_r was already computed to check feasibility, so there is no need to use t_rp

		while (T.length != 0) {
			var median_indice = Math.ceil(T.length/2);
			var t_hat = select_(T, median_indice);
			var g_t_hat = g(t_hat);
			if (Math.abs(g_t_hat - r) <= eps * abs_r) {
				t_star = t_hat;
				break;
			}
			else if (g_t_hat > r) {
				t_l = t_hat;
				var j = 0;
				for (var i = median_indice; i < T.length; ++i) {
					T[j++] = T[i];
				}
				T = resizeArray(T, j);
			}
			else if (g_t_hat < r) {
				t_u = t_hat;
				T = resizeArray(T, median_indice - 1);

			}
			var I_it = new I.iterator();
			var el = I_it.next();
			while (el != 0) {
				var i = el - 1;
				var remove_i = false;
				
				if (T_l[i] <= t_l) {
					s += b.data[i] * l.data[i];

					remove_i = true;
				}
				if (t_u <= T_u[i]) {
					s += b.data[i] * u.data[i];

					remove_i = true;
				}
				if (T_u[i] <= t_l && t_u <= T_l[i]) {
					var b_d = b.data[i] / d.data[i];
					p += a.data[i] * b_d;
					q += b.data[i] * b_d;

					remove_i = true;
				}
				
				if (remove_i === true) {
					I.unset(i+1);
				}
				el = I_it.next();
			}
		}
		if (T.length == 0) {
			t_star = (p + s - r) / q;
		}
	}
	var x_star = x(t_star);
	var fctVal = 1/2 * Matrix_.vectorDotProduct(x_star, Matrix_.elementwiseProduct(x_star, d)) - Matrix_.vectorDotProduct(a, x_star);
	if (outputLagrangeMultiplier === true) {
		return [x_star, fctVal, t_star];
	}
	else {
		return [x_star, fctVal];
	}
}
 function qpsolveGSMO_(Q, p, b, r, l, u, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-04;
	var maxIterations = opt.maxIter || 10000;
	if (!(Q instanceof Matrix_) || !Q.isSquare()) {
		throw new Error('first input must be a square matrix');
	}
	if (!(p instanceof Matrix_) || !p.isVector()) {
		throw new Error('second input must be a vector');
	}
	if (!(b instanceof Matrix_) || !b.isVector()) {
		throw new Error('third input must be a vector');
	}
	if (!(l instanceof Matrix_) || !l.isVector()) {
		throw new Error('fifth input must be a vector');
	}
	if (!(u instanceof Matrix_) || !u.isVector()) {
		throw new Error('sixth input must be a vector');
	}
	
	if (Q.nbRows !== p.nbRows) {
		throw new Error('first and second inputs number of rows do not match: ' + Q.nbRows + '-' + a.nbRows);
	}
	if (Q.nbRows !== b.nbRows) {
		throw new Error('first and third inputs number of rows do not match: ' + Q.nbRows + '-' + b.nbRows);
	}
	if (Q.nbRows !== l.nbRows) {
		throw new Error('first and fifth inputs number of rows do not match: ' + Q.nbRows + '-' + l.nbRows);
	}
	if (Q.nbRows !== u.nbRows) {
		throw new Error('first and sixth inputs number of rows do not match: ' + Q.nbRows + '-' + u.nbRows);
	}
	var n = Q.nbRows;
	var centroid = Matrix_.fill(n, 1, 
								function(i,j) { 
									return r/n* 1/b.data[i-1];
								});
	var p_centroid = qksolveBS_(Matrix_.ones(n, 1), centroid, b, r, l, u);
	var x = p_centroid[0];
	var grad_f_x = Matrix_.xpy(Matrix_.xy(Q, x), p);
	var iter = 0;
	while (true) {
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		var i_low = -1;
		var F_i_low = -Infinity;
		var i_up = -1;
		var F_i_up = Infinity;
		for (var i = 0; i < n; ++i) {
			F_i = grad_f_x.data[i] / b.data[i];
			if (l.data[i] < x.data[i] && x.data[i] < u.data[i]) {
				if (F_i > F_i_low) {
					F_i_low = F_i;
					i_low = i;
				}
				if (F_i < F_i_up) {
					F_i_up = F_i;
					i_up = i;
				}
			}
			else if (x.data[i] == l.data[i]) {
				if (F_i < F_i_up) {
					F_i_up = F_i;
					i_up = i;
				}		
			}
			else if (x.data[i] == u.data[i]) {
				if (F_i > F_i_low) {
					F_i_low = F_i;
					i_low = i;
				}		
			}
		}
		if (F_i_low - F_i_up <= eps) {
			break;
		}
		var i = i_up;
		var j = i_low;
		var t_min_l_i = (l.data[i] - x.data[i])*b.data[i];
		var t_min_u_j = (x.data[j] - u.data[j])*b.data[j];
		var t_min = t_min_l_i >= t_min_u_j ? t_min_l_i : t_min_u_j;
		var t_max_u_i = (u.data[i] - x.data[i])*b.data[i];
		var t_max_l_j = (x.data[j] - l.data[j])*b.data[j];
		var t_max = t_max_u_i <= t_max_l_j ? t_max_u_i : t_max_l_j;
		var dphi_0 = F_i_up - F_i_low;
		var ddphi_0 = Q.data[i*Q.nbColumns + i]/(b.data[i] * b.data[i]) + Q.data[j*Q.nbColumns + j]/(b.data[j] * b.data[j]) - 2*Q.data[i * Q.nbColumns + j]/(b.data[i] * b.data[j]);
		var t_star;
		if (ddphi_0 > 0) { // phi''(0) > 0
			t_star = -dphi_0/ddphi_0;
			if (t_star > t_max) {
				t_star = t_max;
			}
			else if (t_star < t_min) {
				t_star = t_min;
			}
		}
		else if (ddphi_0 == 0) { // phi''(0) = 0
			if (dphi_0 > 0) {
				t_star = t_min;
			}
			else {
				t_star = t_max;
			}
		}
		else { // phi''(0) < 0 implies that Q is not positive semi-definite
			throw new Error("internal error: the input matrix might not be positive semi-definite");
		}
		var old_x_i = x.data[i];
		var old_x_j = x.data[j];
		var delta_x_i = t_star/b.data[i];
		var delta_x_j = -t_star/b.data[j];		
		x.data[i] += delta_x_i;
		x.data[j] += delta_x_j;
		if (t_star == t_min) {
			if (t_min == t_min_l_i) {
				x.data[i] = l.data[i];
				delta_x_i = x.data[i] - old_x_i;
			}
			if (t_min == t_min_u_j) {
				x.data[j] = u.data[j];
				delta_x_j = x.data[j] - old_x_j;
			}
		}
		if (t_star == t_max) {
			if (t_max == t_max_u_i) {
				x.data[i] = u.data[i];
				delta_x_i = x.data[i] - old_x_i;
			}
			if (t_max == t_max_l_j) {
				x.data[j] = l.data[j];
				delta_x_j = x.data[j] - old_x_j;
			}		
		}		
		for (var k = 0; k < n; ++k) {
			grad_f_x.data[k] = grad_f_x.data[k] + Q.data[k*Q.nbColumns + i]*delta_x_i + Q.data[k*Q.nbColumns + j]*delta_x_j;
		}
	}
	var fctVal = 1/2 * Matrix_.vectorDotProduct(x, Matrix_.xy(Q, x)) + Matrix_.vectorDotProduct(p, x);
	return [x, fctVal];
}
 function lpsolvePDHG_(Ae, be, Ai, bi, c, lb, ub, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var eps = opt.eps || 1e-08;
	var maxIterations = opt.maxIter || 100000;
	var eqContraints = false;
	var ineqContraints = false;
	var boundContraints = false;
	
	if (Ae !== null && be !== null) {
		eqContraints = true;
	}
	else if (Ae !== null && be === null) {
		throw new Error('equality constraints vector is missing');
	}
	else if (Ae === null && be !== null) {
		throw new Error('equality constraints matrix is missing');
	}

	if (Ai !== null && bi !== null) {
		ineqContraints = true;
	}
	else if (Ai !== null && bi === null) {
		throw new Error('inequality constraints vector is missing');
	}
	else if (Ai === null && bi !== null) {
		throw new Error('inequality constraints matrix is missing');
	}
	
	if (lb !== null && ub !== null) {
	    boundContraints = true;
	}
	else if (lb !== null && ub === null ) {
    	throw new Error('upper bounds constraints vector is missing');
	}
	else if (lb === null && ub !== null ) {
    	throw new Error('lower bounds constraints vector is missing');
	}
	
	if (!(c instanceof Matrix_)) {
		throw new Error('fifth input must be a matrix');
	}
	if (c.nbColumns !== 1) {
		throw new Error('fifth input is not a vector: ' + c.nbColumns + '-' + c.nbRows);
	}
	
	if (eqContraints) {
		if (!(Ae instanceof Matrix_)) {
			throw new Error('first input must be a matrix');
		}
		if (!(be instanceof Matrix_)) {
			throw new Error('second input must be a matrix');
		}
		if (Ae.nbRows !== be.nbRows) {
			throw new Error('first and second inputs number of rows do not match: ' + Ae.nbRows + '-' + be.nbRows);
		}
		if (Ae.nbColumns !== c.nbRows) {
			throw new Error('first input number of columns and fifth input number of rows do not match: ' + Ae.nbColumns + '-' + c.nbRows);
		}
		if (be.nbColumns !== 1) {
			throw new Error('second input is not a vector: ' + be.nbColumns + '-' + be.nbRows);
		}
	}
	
	if (ineqContraints) {
		if (!(Ai instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (!(bi instanceof Matrix_)) {
			throw new Error('third input must be a matrix');
		}
		if (Ai.nbRows !== bi.nbRows) {
			throw new Error('third and fourth inputs number of rows do not match: ' + Ai.nbRows + '-' + bi.nbRows);
		}
		if (Ai.nbColumns !== c.nbRows) {
			throw new Error('third input number of columns and fifth input number of rows do not match: ' + Ai.nbColumns + '-' + c.nbRows);
		}
		if (bi.nbColumns !== 1) {
			throw new Error('fourth input is not a vector: ' + bi.nbColumns + '-' + bi.nbRows);
		}
	}
	
	if (boundContraints) {
		if (lb.nbRows !== null) {
			if (!(lb instanceof Matrix_)) {
				throw new Error('sixth input must be a matrix');
			}
			if (lb.nbRows !== c.nbRows) {
				throw new Error('sixth input number of rows and fifth input number of rows do not match: ' + lb.nbRows + '-' + c.nbRows);
			}
		}
        if (ub.nbRows !== null) {
			if (!(ub instanceof Matrix_)) {
				throw new Error('seventh input must be a matrix');
			}
			if (ub.nbRows !== c.nbRows) {
				throw new Error('seventh input number of rows and fifth input number of rows do not match: ' + ub.nbRows + '-' + c.nbRows);
			}
		}
	}
	var me = 0;
	var ye_k = null;
	var ye_kp = null;
	var res_ye_kp_ye_k = null;
	if (eqContraints) {
	    me = Ae.nbRows; // the number of equality constaints
    	ye_k = Matrix_.zeros(me, 1); // the dual iterate ye_k
    	ye_kp = Matrix_.zeros(me, 1); // the dual iterate ye_k+1
    	res_ye_kp_ye_k = Matrix_.zeros(me, 1); // the residual ye_kp - ye_k
	}
	var mi = 0;
	var yi_k = null;
	var yi_kp = null;
	var res_yi_kp_yi_k = null;
	if (ineqContraints) {
	    mi = Ai.nbRows; // the number of inequality (<=) constaints
    	yi_k = Matrix_.zeros(mi, 1); // the dual iterate yi_k
    	yi_kp = Matrix_.zeros(mi, 1); // the dual iterate yi_k+1
    	res_yi_kp_yi_k = Matrix_.zeros(mi, 1); // the residual yi_kp - yi_k
	}
	var m = me + mi; // the total number of constraints
    var n = c.nbRows;
	var x_k = Matrix_.ones(n, 1); // the primal iterate x_k
	var x_kp = Matrix_.ones(n, 1); // the primal iterate x_k+1
	var z_k = Matrix_.ones(n, 1); // the relaxed iterate z_k = 2*x_k+1 - x_k
	var res_x_kp_x_k = Matrix_.zeros(n, 1); // the residual x_kp - x_k
	var tmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var ttmp_vec_n = Matrix_.zeros(n, 1); // a temporary placeholder vector of dimension n
	var tmp_vec_me = null;
	if (eqContraints) {
		tmp_vec_me = Matrix_.zeros(me, 1); // a temporary placeholder vector of dimension me
	}
	var tmp_vec_mi = null;
	if (ineqContraints) {
		tmp_vec_mi = Matrix_.zeros(mi, 1); // a temporary placeholder vector of dimension mi
	}
	var mu = 0.9995;
	var nu = 0.9995;
	var T = Matrix_.fill(n, 1, 
				function(i,j) { 
						var columnNorm = 0;
						if (eqContraints) {
						    var aeColNorm = Ae.vectorNorm('one', 'column', i); 
							columnNorm += aeColNorm;
						}
						if (ineqContraints) {
						    var aiColNorm = Ai.vectorNorm('one', 'column', i);
							columnNorm += aiColNorm;
						}
						if (columnNorm == 0) {
							columnNorm = 1;
						}
						return mu * 1/columnNorm;
				});
	
	var Se = null;
	if (eqContraints) {
	    Se = Matrix_.fill(me, 1, 
	                      function(i,j) { 
						    var aeRowNorm = Ae.vectorNorm('one', 'row', i);
							if (aeRowNorm == 0) {
								aeRowNorm = 1;
							}
							return nu * 1/aeRowNorm;
				          });
	}
	var Si = null;
	if (ineqContraints) {
	    Si = Matrix_.fill(mi, 1, 
	                      function(i,j) { 
						    var aiRowNorm = Ai.vectorNorm('one', 'row', i);
							if (aiRowNorm == 0) {
								aiRowNorm = 1;
							}
			        	    return nu * 1/aiRowNorm;
				          });
	}
	var iter = 0;
	while (true) {
		if (maxIterations !== -1 && iter >= maxIterations) {
			throw new Error('maximum number of iterations reached: ' + maxIterations);
		}
		++iter;
		if (eqContraints && ineqContraints) {
			x_kp = Matrix_.xpy(Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), Matrix_.txy(Ai, yi_k, ttmp_vec_n), tmp_vec_n), c, x_kp);
		}
		else if (eqContraints) {
			x_kp = Matrix_.xpy(Matrix_.txy(Ae, ye_k, tmp_vec_n), c, x_kp);
		}
		else if (ineqContraints) {
			x_kp = Matrix_.xpy( Matrix_.txy(Ai, yi_k, tmp_vec_n), c, x_kp);
		}
		x_kp = Matrix_.xmy(x_k, Matrix_.elementwiseProduct(x_kp, T, tmp_vec_n), x_kp);
		if (boundContraints) {
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < lb.data[i-1]) {
					x_kp.data[i-1] = lb.data[i-1];
				}
				else if (x_kp.data[i-1] > ub.data[i-1]) {
					x_kp.data[i-1] = ub.data[i-1];
				}
			}
		}
		else {			
			for (var i = 1; i <= n; ++i) {
				if (x_kp.data[i-1] < 0) {
					x_kp.data[i-1] = 0;
				}
			}
		}
		z_k = Matrix_.axpby(2, x_kp, -1, x_k, z_k);
		if (eqContraints) {
			ye_kp = Matrix_.xpy(ye_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ae, z_k, tmp_vec_me), be, tmp_vec_me), Se, tmp_vec_me), ye_kp);
		}
		if (ineqContraints) {
			yi_kp = Matrix_.xpy(yi_k, Matrix_.elementwiseProduct(Matrix_.xmy(Matrix_.xy(Ai, z_k, tmp_vec_mi), bi, tmp_vec_mi), Si, tmp_vec_mi), yi_kp);
			for (var i = 1; i <= mi; ++i) {
				if (yi_kp.data[i-1] < 0) {
					yi_kp.data[i-1] = 0;
				}
			}
		}
		res_x_kp_x_k = Matrix_.xmy(x_kp, x_k, res_x_kp_x_k);
		var res_x_kp_x_k_inf_norm = res_x_kp_x_k.vectorNorm('infinity');
		var x_kp_inf_norm = x_kp.vectorNorm('infinity');
		
		var res_ye_kp_ye_k_inf_norm = 0;
		var ye_kp_inf_norm = 0;
		if (eqContraints) {
		    res_ye_kp_ye_k = Matrix_.xmy(ye_kp, ye_k, res_ye_kp_ye_k);
		    res_ye_kp_ye_k_inf_norm = res_ye_kp_ye_k.vectorNorm('infinity');
											 
			ye_kp_inf_norm = ye_kp.vectorNorm('infinity');
		}
		
		var res_yi_kp_yi_k_inf_norm = 0;
		var yi_kp_inf_norm = 0;
		if (ineqContraints) {
		    res_yi_kp_yi_k = Matrix_.xmy(yi_kp, yi_k, res_yi_kp_yi_k);
		    res_yi_kp_yi_k_inf_norm = res_yi_kp_yi_k.vectorNorm('infinity');
			
			yi_kp_inf_norm = yi_kp.vectorNorm('infinity');
		}
		
		if (res_x_kp_x_k_inf_norm <= eps * x_kp_inf_norm  && 
		    res_ye_kp_ye_k_inf_norm <= eps * ye_kp_inf_norm &&
			res_yi_kp_yi_k_inf_norm <= eps * yi_kp_inf_norm) {
			break;
		}
		x_k = Matrix_.copy(x_kp, x_k);
		if (eqContraints) {
		    ye_k = Matrix_.copy(ye_kp, ye_k);
		}
		if (ineqContraints) {
		    yi_k = Matrix_.copy(yi_kp, yi_k);
		}
	}
	var fctVal = Matrix_.vectorDotProduct(c, x_kp);
	return [x_kp, fctVal];		
 }
function returns(x, opt) {
	var opt = opt;
	if (opt === undefined) {
		opt = { };
	}
	if (opt.method === undefined) {
		opt.method = "arithmetic";
	}
	var method = opt.method;
	if (method != "arithmetic" &&
	    method != "logarithmic") {
			throw new Error('unsupported returns computation method');
	}
	var returns = typeof Float64Array === 'function' ? new Float64Array(x.length - 1) : new Array(x.length - 1); 

	if (method == "arithmetic") {
		for (var i = 0; i < x.length - 1; ++i) {
			returns[i] = (x[i+1] - x[i])/x[i];
		}
	}
	else if (method == "logarithmic") {
		for (var i = 0; i < x.length - 1; ++i) {
			returns[i] = Math.log(x[i+1]/x[i]);
		}
	}
	else {
		throw new Error('internal error: unsupported returns computation method');
	}
	return returns;
}
function simplexLpSolve_(c, l, u) {
	var eps = 1e-16; // the numerical zero
	
	var c = new Matrix_(c);
	var n = c.nbRows;
	var l = l ? new Matrix_(l) : Matrix_.zeros(n, 1);
	var u = u ? new Matrix_(u) : Matrix_.ones(n, 1);
	simplexEmptinessCheck_(n, l.toArray(), u.toArray());
	var idx_c = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
	for (var j = 0; j < n; ++j) {		
		idx_c[j] = j + 1;
	}
	idx_c.sort(function(a, b) { 
		return c.getValue(a, 1) - c.getValue(b, 1);
	});
	var x = new Matrix_(l);
	var delta_sum_x = 1 - x.sum();
	var idx_i = -1;
	for (var i = 0; i < n; ++i) {	
		if (Math.abs(delta_sum_x) <= eps) {
			break;
		}
		idx_i = idx_c[i];
		var x_idx_i = x.getValue(idx_i, 1);
		if (delta_sum_x >= u.getValue(idx_i, 1) - x_idx_i) {
			x.setValue(idx_i, 1, u.getValue(idx_i, 1));
			delta_sum_x -= u.getValue(idx_i, 1) - x_idx_i;
		}
		else {
			x.setValue(idx_i, 1, x_idx_i + delta_sum_x);
			delta_sum_x = 0;
		}
	}	
	return [x.toArray(), Matrix_.vectorDotProduct(c, x)];
}
function simplexCharacteristicFunction_(x, l, u) {
	var n = x.length;
	var eps_tol = 1e-12; // used to numerically determine some conditions
	simplexEmptinessCheck_(n, l, u);
	var sum_xi = 0;
	for (var i = 0; i < n; ++i) {	
		var lb_i = 0;
		if (l) {
			lb_i = l[i];
		}
		var up_i = 1;
		if (u) {
			up_i = u[i];
		}
		var x_i = x[i];
		if (x_i < lb_i || x_i > up_i) {
			return Number.POSITIVE_INFINITY;
		}
		sum_xi += x_i;
	}
	if (Math.abs(sum_xi - 1) > eps_tol) {
		return Number.POSITIVE_INFINITY;
	}
	return 0;
}
function fullSimplexCharacteristicFunction_(x, l, u) {
	var n = x.length;
	var eps_tol = 1e-12; // used to numerically determine some conditions
	fullSimplexEmptinessCheck_(n, l, u);
	var sum_xi = 0;
	for (var i = 0; i < n; ++i) {	
		var lb_i = 0;
		if (l) {
			lb_i = l[i];
		}
		var up_i = 1;
		if (u) {
			up_i = u[i];
		}
		var x_i = x[i];
		if (x_i < lb_i || x_i > up_i) {
			return Number.POSITIVE_INFINITY;
		}
		sum_xi += x_i;
	}
	if (sum_xi - 1 > eps_tol) {
		return Number.POSITIVE_INFINITY;
	}
	return 0;
}
function fullSimplexEmptinessCheck_(n, l, u) {
	var sumLowerBounds = 0;
	var sumUpperBounds = 0;
	for (var i = 0; i < n; ++i) {
		var lowerBound = 0;
		if (l) {
			lowerBound = l[i];
		}

		var upperBound = 1;
		if (u) {
			upperBound = u[i];
		}
		if (lowerBound > upperBound) {
			throw new Error('infeasible problem detected: lower bound strictly greater than upper bound');
		}
		if (lowerBound < 0) {
			throw new Error('incorrect problem detected: lower bound strictly lower than 0');
		}
		if (upperBound > 1) {
			throw new Error('incorrect problem detected: upper bound strictly greater than 1');
		}
		sumLowerBounds += lowerBound;
		sumUpperBounds += upperBound;		
	}
	
	if (sumLowerBounds > 1) {
		throw new Error('infeasible problem detected: the restricted simplex is empty');
	}
	return [sumLowerBounds, sumUpperBounds];
}
function simplexEmptinessCheck_(n, l, u) {
	var sumLowerBounds = 0;
	var sumUpperBounds = 0;
	for (var i = 0; i < n; ++i) {
		var lowerBound = 0;
		if (l) {
			lowerBound = l[i];
		}

		var upperBound = 1;
		if (u) {
			upperBound = u[i];
		}
		if (lowerBound > upperBound) {
			throw new Error('infeasible problem detected: lower bound strictly greater than upper bound');
		}
		if (lowerBound < 0) {
			throw new Error('incorrect problem detected: lower bound strictly lower than 0');
		}
		if (upperBound > 1) {
			throw new Error('incorrect problem detected: upper bound strictly greater than 1');
		}
		sumLowerBounds += lowerBound;
		sumUpperBounds += upperBound;		
	}
	
	if (sumLowerBounds > 1 || sumUpperBounds < 1) {
		throw new Error('infeasible problem detected: the restricted simplex is empty');
	}
	return [sumLowerBounds, sumUpperBounds];
}
function simplexSparseEuclidianProjection_(x, k, l, u) {
	var n = x.length;
	if (k === n) {
		return simplexEuclidianProjection_(x, l, u);
	}
	if (!l && !u) {
		var idx = typeof UInt32Array === 'function' ? new UInt32Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			idx[i] = i;
		}
		var compareIndexes = function (a, b) {
			return x[b] - x[a];
		};
		select_(idx, k, compareIndexes);
		var x_k = x.slice(0, k);
		for (var i = 0; i < k; ++i) {
			x_k[i] = x[idx[i]];
		}
		var proj_x_k = simplexEuclidianProjection_(x_k);
		var y = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < k;  ++i) {
			y[idx[i]] = proj_x_k[i];
		}
		for (var i = k; i < n;  ++i) {
			y[idx[i]] = 0;
		}
		return y;
	}
	else {
		var x = new Matrix_(x);
		var minDistanceValue = Infinity;
		var minDistanceIndexes = [];
		var minDistanceProjX = null;
		for (var K = 1; K <= k; ++K) {
			var nextKSubsetIterator = new kSubsetsIterator_(n, K, false);
			var nextKSubset = nextKSubsetIterator.next();
			
			while (nextKSubset != -1) {
				var subsetNbIndexes = nextKSubset.length;
				var subsetIndexes = typeof UInt32Array === 'function' ? new UInt32Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				for (var i = 0; i < subsetNbIndexes; ++i) {
					subsetIndexes[i] = nextKSubset[i];
				}
				var subsetX = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				for (var i = 0; i < subsetNbIndexes; ++i) {
					subsetX[i] = x.data[subsetIndexes[i]-1];
				}
				var subsetL = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				if (l) {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetL[i] = l[subsetIndexes[i]-1];
					}
				}
				else {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetL[i] = 0;
					}
				}
				var subsetU = typeof Float64Array === 'function' ? new Float64Array(subsetNbIndexes) : new Array(subsetNbIndexes);
				if (u) {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetU[i] = u[subsetIndexes[i]-1];
					}
				}
				else {
					for (var i = 0; i < subsetNbIndexes; ++i) {
						subsetU[i] = 1;
					}
				}
				try {
					var proj_subsetX = simplexEuclidianProjection_(subsetX, subsetL, subsetU);
					var proj_x = Matrix_.zeros(n, 1);
					for (var i = 0; i < subsetNbIndexes; ++i) {
						proj_x.data[subsetIndexes[i] - 1] = proj_subsetX[i];
					}
					var d_x_proj_x = Matrix_.axpby(1, x, -1, proj_x).vectorNorm('two');
					if (d_x_proj_x < minDistanceValue) {
						minDistanceValue = d_x_proj_x;
						minDistanceIndexes = subsetIndexes;
						minDistanceProjX = proj_x;
					}
					
				}
				catch (e) {
					if (e.message !== "infeasible problem detected: the restricted simplex is empty") {
						throw(e);
					}
				}
				var nextKSubset = nextKSubsetIterator.next();
			}
		}
		if (minDistanceValue != Infinity) {
			return minDistanceProjX.toArray();
		}
		else {
			throw new Error('infeasible problem detected');
		}
	}
}
function fullSimplexEuclidianProjection_(x, l, u) {
	var n = x.length;
	fullSimplexEmptinessCheck_(n, l, u);
	var p_box = new Matrix_(x).elemMap(function(i,j,val) { 
	                                      var l_i = l ? l[i-1] : 0;
										  var u_i = u ? u[i-1] : 1;
	                                      
										  return Math.max(-l_i, Math.min(u_i, val));
									   });
	var a = Matrix_.ones(n, 1);
	var ps = Matrix_.vectorDotProduct(a, p_box);
	
	var p;
	if (ps <= 1) {
		p = p_box.toArray();
	}
	else {
		p = simplexEuclidianProjection_(x, l, u);
	}
	return p;
}
function simplexEuclidianProjection_(x, l, u) {
	var n = x.length;
	simplexEmptinessCheck_(n, l, u);
	var d = Matrix_.ones(n, 1);
	var a = new Matrix_(x);
	var b = Matrix_.ones(n, 1);
	var r = 1;
	var lb = l ? new Matrix_(l) : Matrix_.zeros(n, 1);
	var ub = u ? new Matrix_(u) : Matrix_.ones(n, 1);
	var sol = qksolveBS_(d, a, b, r, lb, ub);
	var y = sol[0].toArray();
	while (simplexCharacteristicFunction_(y, l, u) == Number.POSITIVE_INFINITY) {
		a = new Matrix_(y);
		sol = qksolveBS_(d, a, b, r, lb, ub);
		y = sol[0].toArray();
	}
	return y;
}
function simplexGridSampler_(n, k, reuseOutputArray) {
	this.n = n;
	this.k = k;
	
	this.reuseOutputArray = reuseOutputArray;
	
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.compositionIterator = new compositionsIterator_(k, n, true); // reuse the ouput array for better performances
	this.sample = function() {
		var comp = this.compositionIterator.next();
		if (comp == -1) {
			return -1;
		}
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = comp[i] / this.k;
		}
		if (this.reuseOutputArray) {
			return this.x;
		}
		else {
			return this.x.slice(0);
		}
	}
}
function simplexRandomSampler_(n, l, u, rnd) {
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.u = typeof Float64Array === 'function' ? new Float64Array(n-1) : new Array(n-1); // the coordinates of a random point in the unit hypercube of R^n-1
	this.randomGenerator = function(arr) {
		for (var i = 0; i < this.n; ++i) {
			arr[i] = Math.random();
		}
	}
	if (rnd) {
		if (typeof rnd !== "function") {
			throw new Error('the random number generator must be a function');
		}
		this.randomGenerator = rnd;
	}
	if (l) {
		this.lowerBounds = l; // no input checks
		
		if (!u) {
			this.upperBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			for (var i = 0; i < this.n; ++i) {
				this.upperBounds[i] = 1;
			}
		}
	}
	if (u) {
		this.upperBounds = u; // no input checks
		
		if (!l) {
			this.lowerBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			for (var i = 0; i < this.n; ++i) {
				this.lowerBounds[i] = 0;
			}
		}
	}
	if (this.lowerBounds || this.upperBounds) {
		var sumBounds = simplexEmptinessCheck_(n, this.lowerBounds, this.upperBounds);
		this.sumLowerBounds = sumBounds[0];
		this.sumUpperBounds = sumBounds[1];
		if (this.sumLowerBounds == 1 || this.sumUpperBounds == 1) {
		}
		else {
			var updatedSumLowerBounds = 0;
			var updatedSumUpperBounds = 0;
			for (var i = 0; i < this.n; ++i) {
				var lowerBound = this.lowerBounds[i];
				var upperBound = this.upperBounds[i];
				
				var updatedLowerBound = Math.max(lowerBound, upperBound + 1 - this.sumUpperBounds);
				var updatedUpperBound = Math.min(upperBound, lowerBound + 1 - this.sumLowerBounds);
				
				this.lowerBounds[i] = updatedLowerBound;
				this.upperBounds[i] = updatedUpperBound;

				updatedSumLowerBounds += updatedLowerBound;
				updatedSumUpperBounds += updatedUpperBound;
			}
			this.sumLowerBounds = updatedSumLowerBounds;
			this.sumUpperBounds = updatedSumUpperBounds;
			this.upperStarBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			this.runningSumUpperStarBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			
			var sumUpperStarBounds = 0;
			for (var i = 0; i < this.n; ++i) {
				this.upperStarBounds[i] = (this.upperBounds[i] - this.lowerBounds[i]) / (1-this.sumLowerBounds);
				
				sumUpperStarBounds += this.upperStarBounds[i];
				this.runningSumUpperStarBounds[i] = sumUpperStarBounds;
			}
		}
	}
	this.sample_no_bounds = function() {
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			var u = 1 - Math.random(); // belongs to ]0,1], so that the logarithm below is properly defined
			var e = Math.log(u); // e ~ EXP(1)
			this.x[i] = e;
			sum += e;
		}
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/sum;
		}
		return this.x.slice(0);
	}
	this.sample_binding_bounds = function() {
		var bindingBounds = null;
		if (this.sumLowerBounds == 1) {
			bindingBounds = this.lowerBounds;
		}
		else if (this.sumUpperBounds == 1) {
			bindingBounds = this.upperBounds;
		}
		else {
			throw new Error('internal error');
		}
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = bindingBounds[i];
		}
		return this.x.slice(0);
	}
	this.sample_with_bounds = function() {
		this.randomGenerator(this.u);
		var u = this.u;
		var delta_k = 1;
		for (var k = n; k >= 2; --k) {
			if (Math.abs(delta_k) <= 1e-14) {
				for (var kk = k; kk >= 1; --kk) {
					this.x[kk-1] = 0;
				}
				break;
			}
			var u_k = u[k-2];
			var d_k = Math.max(0, 1 - this.runningSumUpperStarBounds[k-2]/delta_k);			
			var phi_k = Math.min(1, this.upperStarBounds[k-1]/delta_k);			
			var y_k = delta_k * (1 - Math.pow(u_k*Math.pow(1-phi_k, k-1) + (1-u_k)*Math.pow(1-d_k, k-1), 1/(k-1)));
			this.x[k-1] = y_k;
			delta_k -= y_k;			
		}
		if (k == 1) { // In case the main loop above exited with delta not numerically null
			this.x[0] = delta_k; // Update the 1st coordinate of the generated random point
		}
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = (1-this.sumLowerBounds)*this.x[i] + this.lowerBounds[i];
		}
		return this.x.slice(0);
	}
	if (this.lowerBounds || this.upperBounds) {
		if (this.sumLowerBounds == 1 || this.sumUpperBounds == 1) {
			this.sample = this.sample_binding_bounds;
		}
		else {
			this.sample = this.sample_with_bounds;
		}
	}
	else {
		this.sample = this.sample_no_bounds;
	}
}
function simplexDirectionRandomSampler_(n) {
	this.n = n;
	this.x = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n); // the coordinates of a point being sampled
	this.sample = function() {
		var sum = 0;
		for (var i = 0; i < this.n; ++i) {
			var r = normrnd_(0, 1);
			this.x[i] = r;
			sum += r;
		}
		var sum_d_n = sum / this.n;
		var t = 0;
		var s = 1;
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i] - sum_d_n;
			var absX = Math.abs(this.x[i]);
			if (absX != 0) {
				if (absX > t) {
					s = 1 + s * (t/this.x[i]) * (t/this.x[i]);
					t = absX;
				}
				else  {
					s = s + (this.x[i]/t) * (this.x[i]/t);
				}
			}
		}
		var x_two_norm = t * Math.sqrt(s);
		for (var i = 0; i < this.n; ++i) {
			this.x[i] = this.x[i]/x_two_norm;
		}
		return this.x.slice(0);
	}
}
function simplexRationalRounding_(x, r) {
	var k = r;
	var xPartsWithIndexes = new Array(x.length);
	for (var i = 0; i < x.length; ++i) {
		var rx = r * x[i];
		var integerPart = Math.floor(rx);
		var fractionalPart = rx - integerPart;
		k -= integerPart;
		xPartsWithIndexes[i] = [integerPart, fractionalPart, i];
	}
	xPartsWithIndexes.sort(function(a, b) {
		if (b[1] < a[1]) {
			return -1;
		}
		else if (b[1] > a[1]) {
			return 1;
		}
		else { // coordinates have equal fractional parts
			return b[0] - a[0];
		}
	}); 
	var xr = new Array(x.length);
	for (var i = 0; i < k; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = (xPartsWithIndexes[i][0] + 1) / r;
	}
	for (var i = k; i < xPartsWithIndexes.length; ++i) {
		var index = xPartsWithIndexes[i][2];
		xr[index] = xPartsWithIndexes[i][0] / r;
	}
	return xr;
}
function simplexGridSearch_(f, n, k, l, u) {
	var lowerBounds = null;
	var upperBounds = null;
	if (l) {
		lowerBounds = l; // no input checks
		
		if (!u) {
			upperBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			for (var i = 0; i < n; ++i) {
				upperBounds[i] = 1;
			}
		}
	}
	if (u) {
		upperBounds = u; // no input checks
		
		if (!l) {
			lowerBounds = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
			for (var i = 0; i < n; ++i) {
				lowerBounds[i] = 0;
			}
		}
	}
	if (lowerBounds || upperBounds) {
		var sumBounds = simplexEmptinessCheck_(n, lowerBounds, upperBounds);
	}
	var minValue = Number.POSITIVE_INFINITY;
	var minValueGridPoints = [];
	var sampler = new simplexGridSampler_(n, k, true); // use no array copy in the simplex grid sampler to improve performances
	var weights = sampler.sample();
	while (weights !== -1) {  
		var withinBounds = true;
		if (lowerBounds) {
			for (var i = 0; i < n; ++i) {
				if (lowerBounds[i] > weights[i]) {
					withinBounds = false;
					break;
				}
			}
		}
		if (upperBounds) {
			for (var i = 0; i < n; ++i) {
				if (weights[i] > upperBounds[i]) {
					withinBounds = false;
					break;
				}
			}
		}
		if (!withinBounds) {
			weights = sampler.sample();
			continue;
		}
		var fctValue = f(weights);
		if (fctValue < minValue) {
			minValue = fctValue;
			minValueGridPoints = [weights.slice(0)];
		}
		else if (fctValue == minValue) {
			minValueGridPoints.push(weights.slice(0));
		}
		weights = sampler.sample();
	}
	return minValueGridPoints;
}
function bestConstantlyRebalancedWeights (priceRelatives, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	if (opt.eps === undefined) {
		opt.eps = 1e-04;
	}
	if (opt.maxIter === undefined) {
		opt.maxIter = 10000;
	}
	var eps = opt.eps;
	var maxIterations = opt.maxIter;
	var lowerBounds = opt.constraints.minWeights;
	var upperBounds = opt.constraints.maxWeights;
	var nbAssets = priceRelatives.length; // m in the reference
	var nbPeriods = priceRelatives[0].length; // n in the reference
	function f(b) {
		var x_k = Matrix_.zeros(nbAssets, 1);
		var prodPortfolioRelatives = 1.0;
		for (var k = 0; k < nbPeriods; ++k) {
			x_k = Matrix_.fill(nbAssets, 1, function(i,j) { return priceRelatives[i-1][k]; }, x_k);
			prodPortfolioRelatives *= Matrix_.vectorDotProduct(b, x_k);
		}
		return -prodPortfolioRelatives;
	}
	function gradf(b) {		
		var x_k = Matrix_.zeros(nbAssets, 1);
		var portfolioRelatives = typeof Float64Array === 'function' ? new Float64Array(nbPeriods) : new Array(nbPeriods);
		var prodPortfolioRelatives = 1.0;
		var nullPortfolioRelatives = false;
		for (var k = 0; k < nbPeriods; ++k) {
			x_k = Matrix_.fill(nbAssets, 1, function(i,j) { return priceRelatives[i-1][k]; }, x_k);
			var b_d_x_k = Matrix_.vectorDotProduct(b, x_k);
			portfolioRelatives[k] = b_d_x_k;
			if (nullPortfolioRelatives === false) {
				prodPortfolioRelatives *= b_d_x_k;
				if (Math.abs(b_d_x_k) <= 1e-12) {
					nullPortfolioRelatives = true;
				}
			}
		}
		var partialPortfolioRelatives;
		if (nullPortfolioRelatives === false) {
			partialPortfolioRelatives = Matrix_.fill(nbPeriods, 1, 
			                                         function(i,j) { return prodPortfolioRelatives / portfolioRelatives[i-1]; });
		}
		else {
			throw new Error('null portfolio relative detected, unsuported case');
		}
		var xx_k = Matrix_.zeros(nbPeriods, 1);
		var res = Matrix_.zeros(nbAssets, 1);
		for (var k = 0; k < nbAssets; ++k) {
			xx_k = Matrix_.fill(nbPeriods, 1, function(i,j) { return priceRelatives[k][i-1]; }, xx_k);
			res.data[k] = -Matrix_.vectorDotProduct(xx_k, partialPortfolioRelatives);
		}
		return res;
	}
	function g(b) {
		return simplexCharacteristicFunction_(b.data, lowerBounds, upperBounds);
	}
	function proxg(b) {
		return new Matrix_(simplexEuclidianProjection_(b.data, lowerBounds, upperBounds));
	}
	var x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(nbAssets, 1).data, lowerBounds, upperBounds));
	var sol = ccpsolveFISTA_(f, gradf, g, proxg, x0, {eps: eps, maxIter: maxIterations, maxLine: maxIterations});
	var weights = sol[0];
	return weights.toArray();
}
function clusterRiskParityWeights (sigma, opt) {
	if (opt === undefined) {
		opt = {};
	}
	var clusteringMode = opt.clusteringMode || 'ftca';
	var sigma = new Matrix_(sigma).toCovarianceMatrix(sigma);
	var nbAssets = sigma.nbRows;
	var clusters = [];
	if (clusteringMode === 'manual') {
		clusters = opt.clusters || [];
		var partition = new Array(nbAssets);
		for (var i = 0; i < partition.length; ++i) {
		    partition[i] = 0;
		}
	 	for (var j = 0; j < clusters.length; ++j) {
			if (clusters[j].length === 0) {
				throw new Error('empty cluster at index: ' + j);
			}
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
	}
	else if (clusteringMode === 'ftca') {
		var corrMat = sigma.getCorrelationMatrix().toRowArray();
		clusters = ftca_(corrMat, opt.ftcaThreshold);
	}
	else {
		throw new Error('unsupported clustering method');
	}
	var nbClusters = clusters.length;
	var assetsToClustersWeights = Matrix_.zeros(nbClusters, nbAssets); 
	for (var i = 0; i < nbClusters; ++i) {
		var assetsIndexes = clusters[i].slice().sort(function(a, b) { return a - b; });
		
		var clusterSigma = sigma.submatrix(assetsIndexes, assetsIndexes);
		var assetsWeights = equalRiskContributionWeights(clusterSigma, opt);
		for (var j = 0; j < assetsWeights.length; ++j) {
			assetsToClustersWeights.setValueAt(i+1, assetsIndexes[j], assetsWeights[j]);
		}
	}
	var clustersSigma = Matrix_.xy(assetsToClustersWeights, Matrix_.axty(1, sigma, assetsToClustersWeights));
	var clustersWeights = equalRiskContributionWeights(clustersSigma, opt);
	var weights = Matrix_.txy(assetsToClustersWeights, new Matrix_(clustersWeights));
	return weights.toArray();
}
function equalRiskBoundingWeights (sigma, opt) {	
	if (opt === undefined) {
		opt = {};
	}
	opt.outputPortfolioVolatility = true;
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var minRCValue = Infinity;
	var minRCAssetsIndexes = [];
	var minRCAssetsWeights = [];
	var nextSubsetIterator = new subsetsIterator_(nbAssets);
	var nextSubset = nextSubsetIterator.next(); // empty set
	var nextSubset = nextSubsetIterator.next(); // "true" first set
	while (nextSubset != -1) {	
		var subsetAssetsIdx = nextSubset;
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);
		var sol = equalRiskContributionWeights(subsetSigma, opt);
		var subsetAssetsWeights = sol[0];
		var subsetPortfolioVolatility = sol[1];
		var rcValue = subsetPortfolioVolatility * subsetPortfolioVolatility / subsetAssetsIdx.length;
		if (rcValue < minRCValue) {
			minRCValue = rcValue;
			minRCAssetsIndexes = subsetAssetsIdx;
			minRCAssetsWeights = subsetAssetsWeights;
		}
		var nextSubset = nextSubsetIterator.next();
	}
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < minRCAssetsIndexes.length; ++i) {
		weights.setValueAt(minRCAssetsIndexes[i], 1, 
		                   minRCAssetsWeights[i]);
	}
	return weights.toArray();
}
function equalRiskContributionWeights (sigma, opt) {	
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var rb = Matrix_.fill(nbAssets, 1, function (i,j) { return 1/nbAssets; });
	return riskBudgetingWeights(sigma, rb, opt);
}
function equalWeights (nbAssets, opt) {
	var weights = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
	return weights.toArray();
}
function globalMinimumVarianceWeights (sigma, opt) {
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var mu;
	if (opt.mu == undefined) {
		mu = Matrix_.zeros(nbAssets, 1);
	}
	else {
		mu = new Matrix_(opt.mu);
	}
	var efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	var lowestRiskTolerancePortfolioWeights = efficientFrontier.getLowestRiskTolerancePortfolio();
	return lowestRiskTolerancePortfolioWeights.toArray();
}
function inverseVolatilityWeights (sigma, opt) {
	var weights = new Matrix_(sigma).elemMap(function(i,j,val) { return 1/Math.sqrt(val); })
	weights = weights.normalize(weights);
	return weights.toArray();
}
function maximumSharpeRatioWeights(mu, sigma, rf, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'gsmo';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo') {
		throw new Error('unsupported optimisation method');
	}
	var epsVolatility = opt.epsVolatility;
	if (epsVolatility === undefined) {
		epsVolatility = 1e-4;
	}
	var efficientFrontier;
	if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	efficientFrontier.restrict("minVolatility", epsVolatility);
	efficientFrontier.restrict("minReturn", rf);
	var portfolio = efficientFrontier.computeMaximumSharpeRatioEfficientPortfolio(rf);
	return portfolio[0].toArray();
}
function MeanVarianceEfficientFrontierCla(mu, sigma, opt) {
	MeanVarianceEfficientFrontier.call(this, mu, sigma, opt);
	if (this.fullInvestment == true) {
		this.cornerPortfolios = computeCornerPortfolios(this.mu, this.sigma, this.lowerBounds, this.upperBounds, opt);
	}
	else {
		var that = this;
		
		this.alteredCornerPortfolios = computeCornerPortfolios(this.alteredMu, this.alteredSigma, this.alteredLowerBounds, this.alteredUpperBounds, opt);
		
		this.cornerPortfolios = new Array(this.alteredCornerPortfolios.length);
		for (var k = 0; k < this.alteredCornerPortfolios.length; ++k) {
			var alteredCornerPortfolio_k = this.alteredCornerPortfolios[k][0];
			var riskTolerance_k = this.alteredCornerPortfolios[k][1];
			
			var cornerPortfolio_k = Matrix_.fill(this.nbAssets, 1, 
												function(i, j) { 
													if (i <= that.nbAssets) {
														return alteredCornerPortfolio_k.getValue(i, 1);
													}
												});
			
			this.cornerPortfolios[k] = [cornerPortfolio_k, riskTolerance_k];			
		}
	}
	if (this.cornerPortfolios.length == 0) {
		throw new Error('internal error: no corner portfolio could be computed');
	}
	function computeCornerPortfolios(mu, sigma, lowerBounds, upperBounds, opt) {	
		var eps = 1e-8; // the numerical zero
		function variablesStatusManager_(nbAssets, nbEqualityConstraints) {
			this.STATUS_UNDEF = -1;
			this.STATUS_IN = 0;
			this.STATUS_LOW = 1;
			this.STATUS_UP = 2;
			this.nbAssets = nbAssets;
			this.nbEqualityConstraints = nbEqualityConstraints;
			
			this.varIn = new BitSet_();
			this.varIn.resize(nbAssets + nbEqualityConstraints);
			this.varOut = new BitSet_();
			this.varOut.resize(nbAssets + nbEqualityConstraints);
			this.varLow = new BitSet_();
			this.varLow.resize(nbAssets + nbEqualityConstraints);
			this.varUp = new BitSet_();
			this.varUp.resize(nbAssets + nbEqualityConstraints);
			this.setIn = function(idx) {
				this.varIn.set(idx);
				this.varOut.unset(idx);
				this.varLow.unset(idx);
				this.varUp.unset(idx);
			}
			this.setOnLowerBound = function(idx) {
				this.varLow.set(idx);
				this.varOut.set(idx);
				this.varIn.unset(idx);
				this.varUp.unset(idx);
			}
			this.setOnUpperBound = function(idx) {
				this.varUp.set(idx);
				this.varOut.set(idx);
				this.varIn.unset(idx);
				this.varLow.unset(idx);
			}
			this.setLambdasIn = function() {
				for (var i = this.nbAssets + 1; i <= this.nbAssets + this.nbEqualityConstraints; ++i) {
					this.varIn.set(i);
					this.varOut.unset(i);
					this.varLow.unset(i);
					this.varUp.unset(i);
				}
			}
			this.isAsset = function(idx) {
				return (idx >= 1) && (idx <= this.nbAssets);
			}
			this.isLambda = function(idx) {
				return (idx >= this.nbAssets + 1) && (idx <= this.nbAssets + this.nbEqualityConstraints);
			}
			this.isIn = function(idx) {
				return this.varIn.get(idx);
			}
			this.isOnLowerBound = function(idx) {
				return this.varLow.get(idx);
			}
			this.isOnUpperBound = function(idx) {
				return this.varUp.get(idx);
			}
			this.isOut = function(idx) {
				return this.varOut.get(idx);
			}
			this.getInIndexes = function() {
				return this.varIn.toArray();
			}
			this.getOutIndexes = function() {
				return this.varOut.toArray();
			}
		}
		if (opt === undefined) {
			opt = { constraints: {} };
		}
		if (opt.optimizationMethodParams === undefined) {
			opt.optimizationMethodParams = {};
		}
		var maxIterations = opt.optimizationMethodParams.maxIterCriticalLine;
		if (maxIterations == undefined) {
			maxIterations = 1000;
		}
		var nbAssets = sigma.nbColumns;
		var A = Matrix_.ones(1, nbAssets); // the matrix holding the equality constraints
		var nbEqualityConstraints = A.nbRows; // the number of equality constraints 
		var b = Matrix_.ones(1, 1); // the vector holding the right member of the equality constraints
		
		var lb = lowerBounds;
		var ub = upperBounds;
		
		var cornerPortfoliosWeights = [];
		var currentCornerPortfolioWeights = null;
		
		var variablesStatusManager = new variablesStatusManager_(nbAssets, nbEqualityConstraints);
		var mu_idx = typeof Uint32Array === 'function' ? new Uint32Array(nbAssets) : new Array(nbAssets);
		for (var j = 0; j < nbAssets; ++j) {		
			mu_idx[j] = j + 1;
		}
		mu_idx.sort(function(a, b) {  // Order the assets in descending order w.r.t. their returns
			return mu.getValue(b, 1) - mu.getValue(a, 1);
		});

		for (var i = 1; i < nbAssets; ++i) {
			if (mu.getValue(mu_idx[i], 1) == mu.getValue(mu_idx[i-1], 1)) {
				throw new Error('unsupported problem: equal returns are not supported by the critical line algorithm');
			}
		}
		var maxReturnSolution = simplexLpSolve_(mu.elemMap(function(i,j,val) { return -val;}), lb, ub);
		currentCornerPortfolioWeights = new Matrix_(maxReturnSolution[0]);
		var epsBounds = 1e-14;
		var cptIn = 0;
		var cptLow = 0;
		var cptUp = 0;
		for (var i = 1; i <= nbAssets; ++i) {
			if (Math.abs(currentCornerPortfolioWeights.data[i-1] - lb.data[i-1]) <= epsBounds) {
				variablesStatusManager.setOnLowerBound(i);
				++cptLow;
			}
			else if (Math.abs(currentCornerPortfolioWeights.data[i-1] - ub.data[i-1]) <= epsBounds) {
				variablesStatusManager.setOnUpperBound(i);
				++cptUp;
			}
			else {
				variablesStatusManager.setIn(i);
				++cptIn;
			}
		}
		if (cptLow == nbAssets || cptUp == nbAssets) {
			var weights = new Matrix_(currentCornerPortfolioWeights);
			cornerPortfoliosWeights.push([weights, 0]);

			return cornerPortfoliosWeights;
		}
		else if (cptIn == 0) {
			for (var i = 1; i <= nbAssets - 1; ++i) {
				var idx_i = mu_idx[i-1];
				var idx_ip = mu_idx[i];
				if (variablesStatusManager.isOnUpperBound(idx_i) && variablesStatusManager.isOnLowerBound(idx_ip)) {
					variablesStatusManager.setIn(idx_i);
				}
			}
		}
		var assetsInIdx = variablesStatusManager.getInIndexes();
		var assetsOutIdx = variablesStatusManager.getOutIndexes();
		var xi = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
		var alpha = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
		for (var i = 1; i <= assetsOutIdx.length; ++i) {
			var out_idx_i = assetsOutIdx[i-1];
			
			alpha.setValue(out_idx_i, 1, 
						   currentCornerPortfolioWeights.getValue(out_idx_i, 1));
		}
		var beta = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
		var M = Matrix_.fill(nbAssets + nbEqualityConstraints, nbAssets + nbEqualityConstraints, 
									function(i,j) { 
										if (i <= nbAssets && j <= nbAssets) {
											return sigma.data[(i-1)*sigma.nbColumns + (j-1)]; // Sigma(i, j)
										}
										else if (i >= nbAssets + 1 && j <= nbAssets) {
											return A.data[(i-nbAssets-1)*A.nbColumns + (j-1)]; // A(i-nbAssets, j)
										}
										else if (i <= nbAssets && j >= nbAssets + 1) {
											return A.data[(j-nbAssets-1)*A.nbColumns + (i-1)]; // A(j-nbAssets, i) == A(i-nbAssets, j)^t
										}
										else {
											return 0;
										}
									});	
		var Ai = Matrix_.ones(1, 1);
		var Mi = Matrix_.zeros(nbAssets + nbEqualityConstraints, nbAssets + nbEqualityConstraints);
		var T = Matrix_.atxy(-1, Ai, Matrix_.xy(sigma.submatrix(assetsInIdx, assetsInIdx), Ai));
		for (var j = 1; j <= assetsInIdx.length; ++j) {
			for (var k = 1; k <= assetsInIdx.length; ++k) {
				var var_in_idx_j = assetsInIdx[j-1];
				
				var Ai_j_k = Ai.getValue(j, k);
				Mi.setValue(nbAssets + k, var_in_idx_j, 
							Ai_j_k);		
				Mi.setValue(var_in_idx_j, nbAssets + k, 
							Ai_j_k);
				
				Mi.setValue(nbAssets + j, nbAssets + k, 
							T.getValue(j, k)); 
			}
		}
		variablesStatusManager.setLambdasIn();
		var b_bar = Matrix_.zeros(nbAssets + nbEqualityConstraints, 1);
		for (var i = 1; i <= assetsInIdx.length; ++i) {
			var in_idx_i = assetsInIdx[i-1];
			var b_bar_in_idx_i = 0;
			for (var j = 1; j <= assetsOutIdx.length; ++j) {
				var out_idx_j = assetsOutIdx[j-1];
				
				b_bar_in_idx_i -= M.getValue(in_idx_i, out_idx_j) * currentCornerPortfolioWeights.getValue(out_idx_j, 1);
			}
			b_bar.setValue(in_idx_i, 1, 
						   b_bar_in_idx_i);
		}
		for (var i = nbAssets + 1; i <= nbAssets + nbEqualityConstraints; ++i) {
			var in_idx_i = i;
			var b_bar_in_idx_i = b.getValue(in_idx_i-nbAssets, 1);
			for (var j = 1; j <= assetsOutIdx.length; ++j) {
				var out_idx_j = assetsOutIdx[j-1];
				
				b_bar_in_idx_i -= M.getValue(in_idx_i, out_idx_j) * currentCornerPortfolioWeights.getValue(out_idx_j, 1);
			}
			b_bar.setValue(in_idx_i, 1, 
						   b_bar_in_idx_i);
		}	
		var iter = 0;
		var lambda_e = 0;	
		var idx_out = -1;
		var lambda_out = 0;
		var status_out = variablesStatusManager.STATUS_UNDEF;
		var idx_in = -1;
		var lambda_in = 0;
		while (true) {
			if (maxIterations !== -1 && iter >= maxIterations) {
				throw new Error('maximum number of iterations reached: ' + maxIterations);
			}
			++iter;
			if (iter >= 2) {
				if (lambda_out >= lambda_in) { // an asset IN goes OUT
					alpha.setValue(idx_out, 1, 
								   currentCornerPortfolioWeights.getValue(idx_out, 1));
					beta.setValue(idx_out, 1, 
								  0);
					if (status_out == variablesStatusManager.STATUS_LOW) {
						variablesStatusManager.setOnLowerBound(idx_out);
					}
					else {
						variablesStatusManager.setOnUpperBound(idx_out);
					}
					var variablesInIdx = variablesStatusManager.getInIndexes();
					var Mi_out_idx_out_idx = Mi.getValue(idx_out, idx_out);				
					if (Mi_out_idx_out_idx <= eps) {
						Mi_out_idx_out_idx = eps;
					}
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						for (var j = 1; j <= variablesInIdx.length; ++j) {
							var in_idx_j = variablesInIdx[j-1];
							
							Mi.setValue(in_idx_i, in_idx_j, 
										Mi.getValue(in_idx_i, in_idx_j) - Mi.getValue(in_idx_i, idx_out) * Mi.getValue(idx_out, in_idx_j) / Mi_out_idx_out_idx);
						}
						
					}
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						b_bar.setValue(in_idx_i, 1, 
									   b_bar.getValue(in_idx_i, 1) - M.getValue(in_idx_i, idx_out) * currentCornerPortfolioWeights.getValue(idx_out, 1));
					}
				}
				else { // an asset OUT goes IN				
					var variablesInIdx = variablesStatusManager.getInIndexes();
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						var xi_in_idx_i = 0;
						for (var j = 1; j <= variablesInIdx.length; ++j) {
							var in_idx_j = variablesInIdx[j-1];
							
							xi_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * M.getValue(in_idx_j, idx_in);
						}	
						xi.setValue(in_idx_i, 1, 
									xi_in_idx_i);
					}
					var xi_j = M.getValue(idx_in, idx_in);
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						xi_j -= M.getValue(idx_in, in_idx_i) * xi.getValue(in_idx_i, 1);
					}
					if (xi_j <= eps) {
						xi_j = eps;
					}
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						for (var j = 1; j <= variablesInIdx.length; ++j) {
							var in_idx_j = variablesInIdx[j-1];
							
							Mi.setValue(in_idx_i, in_idx_j, 
										Mi.getValue(in_idx_i, in_idx_j) + xi.getValue(in_idx_i, 1) * xi.getValue(in_idx_j, 1) / xi_j);
						}
						Mi.setValue(in_idx_i, idx_in, 
									-xi.getValue(in_idx_i, 1)/xi_j);
						Mi.setValue(idx_in, in_idx_i, 
									-xi.getValue(in_idx_i, 1)/xi_j);
					}
					Mi.setValue(idx_in, idx_in, 
								1/xi_j);
					for (var i = 1; i <= variablesInIdx.length; ++i) {
						var in_idx_i = variablesInIdx[i-1];
						
						b_bar.setValue(in_idx_i, 1, 
									   b_bar.getValue(in_idx_i, 1) + M.getValue(in_idx_i, idx_in) * currentCornerPortfolioWeights.getValue(idx_in, 1));
					}
					variablesStatusManager.setIn(idx_in);
					var assetsOutIdx = variablesStatusManager.getOutIndexes();
					var b_bar_in_idx_i = 0;
					for (var i = 1; i <= assetsOutIdx.length; ++i) {
						var out_idx_i = assetsOutIdx[i-1];
						
						b_bar_in_idx_i -= M.getValue(idx_in, out_idx_i) * currentCornerPortfolioWeights.getValue(out_idx_i, 1);
					}
					b_bar.setValue(idx_in, 1, 
								   b_bar_in_idx_i);
					
				}
			}
			var variablesInIdx = variablesStatusManager.getInIndexes();
			var assetsOutIdx = variablesStatusManager.getOutIndexes(); // only assets indexes per construction
			idx_out = -1;
			lambda_out = 0;
			status_out = variablesStatusManager.STATUS_UNDEF;
			for (var i = 1; i <= variablesInIdx.length; ++i) {
				var in_idx_i = variablesInIdx[i-1];
				var alpha_in_idx_i = 0;
				var beta_in_idx_i = 0;
				for (var j = 1; j <= variablesInIdx.length; ++j) {
					var in_idx_j = variablesInIdx[j-1];
					
					alpha_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * b_bar.getValue(in_idx_j, 1);
				
					if (in_idx_j <= nbAssets) {
						beta_in_idx_i += Mi.getValue(in_idx_i, in_idx_j) * mu.getValue(in_idx_j, 1);
					}
				}
				alpha.setValue(in_idx_i, 1, 
							   alpha_in_idx_i);
				beta.setValue(in_idx_i, 1, 
							  beta_in_idx_i);
				if (variablesStatusManager.isAsset(in_idx_i)) {
					if (beta_in_idx_i > eps) {
						var lb_idx_in_i = lb.getValue(in_idx_i, 1);
						
						var tmp_lambda_out = (lb_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
						if (tmp_lambda_out >= lambda_out) {
							idx_out = in_idx_i;
							lambda_out = tmp_lambda_out;
							status_out = variablesStatusManager.STATUS_LOW;
						}
					}
					else if (beta_in_idx_i < -eps) {
						var ub_idx_in_i = ub.getValue(in_idx_i, 1);
						
						var tmp_lambda_out = (ub_idx_in_i - alpha_in_idx_i)/beta_in_idx_i;
						if (tmp_lambda_out >= lambda_out) {
							idx_out = in_idx_i;
							lambda_out = tmp_lambda_out;
							status_out = variablesStatusManager.STATUS_UP;
						}
					}
				}
			}
			idx_in = -1;
			lambda_in = 0;
			for (var i = 1; i <= assetsOutIdx.length; ++i) {
				var out_idx_i = assetsOutIdx[i-1];
				var gamma_out_idx_i = 0;
				var delta_out_idx_i = -mu.getValue(out_idx_i, 1);
				for (var j = 1; j <= nbAssets + nbEqualityConstraints; ++j) {
					var M_out_idx_i_j = M.getValue(out_idx_i, j);
					
					gamma_out_idx_i += M_out_idx_i_j * alpha.getValue(j, 1);
					delta_out_idx_i += M_out_idx_i_j * beta.getValue(j, 1);
				}
				if (variablesStatusManager.isOnLowerBound(out_idx_i)) {
					if (delta_out_idx_i > eps) {
						var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
					
						if (tmp_lambda_in >= lambda_in) {
							idx_in = out_idx_i;
							lambda_in = tmp_lambda_in;
						}
					}
				}
				else {
					if (delta_out_idx_i < -eps) {
						var tmp_lambda_in = -gamma_out_idx_i/delta_out_idx_i;
						
						if (tmp_lambda_in >= lambda_in) {
							idx_in = out_idx_i;
							lambda_in = tmp_lambda_in;
						}
					}
				}
			}
			lambda_e = Math.max(lambda_out, lambda_in, 0);
			for (var i = 1; i <= variablesInIdx.length; ++i) {
				var in_idx_i = variablesInIdx[i-1];
				if (variablesStatusManager.isAsset(in_idx_i)) {
					currentCornerPortfolioWeights.setValue(in_idx_i, 1, 
														   alpha.getValue(in_idx_i, 1) + lambda_e * beta.getValue(in_idx_i, 1));
				}
			}
			var weights = new Matrix_(currentCornerPortfolioWeights);
			cornerPortfoliosWeights.push([weights, lambda_e]);
			if (lambda_e < eps) {
				break;
			}
		}
		return cornerPortfoliosWeights;	
	}
};
MeanVarianceEfficientFrontierCla.prototype = Object.create(MeanVarianceEfficientFrontier.prototype);
MeanVarianceEfficientFrontierCla.prototype.constructor = MeanVarianceEfficientFrontierCla;

MeanVarianceEfficientFrontierCla.prototype.getHighestRiskTolerancePortfolio = function(x) {
	return this.cornerPortfolios[0][0];
};
MeanVarianceEfficientFrontierCla.prototype.getHighestRiskTolerance = function(x) {
	return this.cornerPortfolios[0][1];
};
MeanVarianceEfficientFrontierCla.prototype.getLowestRiskTolerancePortfolio = function(x) {
	return this.cornerPortfolios[this.cornerPortfolios.length-1][0];
};
MeanVarianceEfficientFrontierCla.prototype.getLowestRiskTolerance = function(x) {
	return this.cornerPortfolios[this.cornerPortfolios.length-1][1];
};
MeanVarianceEfficientFrontierCla.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}

	var cornerPortfolioConstraintFct;
	var that = this;
	if (constraintType == "return") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return that.computePortfolioReturn(cornerPortfolio[0]); };
	}
	else if (constraintType == "volatility") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return that.computePortfolioVolatility(cornerPortfolio[0]); };
	}
	else if (constraintType == "riskTolerance") {
		cornerPortfolioConstraintFct = function (cornerPortfolio) { return cornerPortfolio[1]; };
	}
	else {
		throw new Error('unknown constraint type');
	}
	
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}
	var enclosingCornerPortfolios = computeEnclosingCornerPortfolios(cornerPortfolioConstraintFct, 
	                                                                 constraintValue, 
																	 this.cornerPortfolios,
																	 this.epsEfficientPortfolioComputation);
	if (enclosingCornerPortfolios.length == 0) {
		return [];
	}
	else if (enclosingCornerPortfolios.length == 1) {
		var idx_min = enclosingCornerPortfolios[0];
		var weights = this.cornerPortfolios[idx_min][0];
		var lambda = this.cornerPortfolios[idx_min][1];
		return [weights, lambda];
	}
	else {
		var idx_min = enclosingCornerPortfolios[0];
		var weights_min = this.cornerPortfolios[idx_min][0];
		var fct_min = cornerPortfolioConstraintFct(this.cornerPortfolios[idx_min]);

		var idx_max = enclosingCornerPortfolios[1];
		var weights_max = this.cornerPortfolios[idx_max][0];
		var fct_max = cornerPortfolioConstraintFct(this.cornerPortfolios[idx_max]);
		var t;
		if (constraintType === "return") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}		
		else if (constraintType === "volatility") {			
			var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
			var a = fct_min * fct_min + fct_max * fct_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
			var b = -2 * (fct_max * fct_max - variance_cross); // 
			var c = fct_max * fct_max - constraintValue*constraintValue; //always > 0
			var b_p = b/2; // reduced discriminant
			var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
			var disc = b_p*b_p - a*c;
			if (disc < 0) {
				throw new Error('internal error; the covariance matrix might not be semi-definite positive');
			}
			var q = -(b_p + sign_b_p * Math.sqrt(disc));
			var r1 = q/a;
			var r2 = c/q;
			
			if (r1 > 0 && r1 < 1) {
				t = r1;
			}
			else if (r2 > 0 && r2 < 1) {
				t = r2;
			}
			else {
				throw new Error('internal error: the covariance matrix might not be semi-definite positive');
			}
		}
		else if (constraintType === "riskTolerance") {
			t = (fct_max - constraintValue)/(fct_max - fct_min);
		}
		else {
			throw new Error('internal error: unknown constraint type');
		}
		var weights = Matrix_.fill(weights_min.nbRows, 1, 
							   	   function(i,j) { 
									   return t*weights_min.getValue(i, 1) + (1-t)*weights_max.getValue(i, 1); 
								   });
		var lambda = t*this.cornerPortfolios[idx_min][1] + (1-t)*this.cornerPortfolios[idx_max][1];
		return [weights, lambda];	
	}
	function computeEnclosingCornerPortfolios(fct, fctValue, cornerPortfolios, eps) {		
		var idx_min = cornerPortfolios.length - 1;
		var idx_max = 0

		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];
		
		var fct_min = fct(cornerPortfolios[idx_min]);
		var fct_max = fct(cornerPortfolios[idx_max]);
		if (fctValue > fct_max + eps || fctValue < fct_min - eps) {
			return [];
		}
		if (Math.abs(fctValue - fct_min) <= eps) {
			return [idx_min];
		}
		else if (Math.abs(fctValue - fct_max) <= eps) {
			return [idx_max];
		}
		while (idx_min - idx_max != 1) { 
			var idx_middle = Math.floor((idx_min + idx_max)/2);
			var fct_middle = fct(cornerPortfolios[idx_middle]);
			if (fct_middle > fctValue) {
				idx_max = idx_middle;
			}
			else if (fct_middle < fctValue) {
				idx_min = idx_middle;
			}
			else { // the target constraint function value is exactly attained on the idx_middle-th corner portfolio
				return [idx_middle];
			}
		}
		return [idx_min, idx_max];
	}
};
MeanVarianceEfficientFrontierCla.prototype.getCornerPortfolios = function() {
	var portfolios = new Array(this.cornerPortfolios.length);
	for (var i = 0; i < portfolios.length; ++i) {
		portfolios[i] = new Matrix_(this.cornerPortfolios[i][0]);
	}
	return portfolios;
};
MeanVarianceEfficientFrontierCla.prototype.restrict = function(constraintType, constraintValue) {
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}
	if (constraintType == "minReturn") {	
		var lowestReturn = this.getLowestReturn();
		if (lowestReturn < constraintValue) {
			var highestReturn = this.getHighestReturn();
			if (highestReturn < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum return constraint is not feasible');
			}
			var efficientPortfolio = this.computeEfficientPortfolio("return", constraintValue + this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a return greater than ' + constraintValue);
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			var updatedCornerPortfolios = new Array();
			for (i = 0; i < this.cornerPortfolios.length - 1; ++i) {	
				if (this.cornerPortfolios[i][1] > efficientPortfolioRiskTolerance) {
					updatedCornerPortfolios.push(this.cornerPortfolios[i]);
				}
				else {
					break;
				}
			}
			updatedCornerPortfolios.push([efficientPortfolioWeights, efficientPortfolioRiskTolerance]);
			this.cornerPortfolios = updatedCornerPortfolios;
		}
		
	}
	else if (constraintType == "minVolatility") {
		var lowestVolatility = this.getLowestVolatility();
		if (lowestVolatility < constraintValue) {
			var highestVolatility = this.getHighestVolatility();
			if (highestVolatility < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum volatility constraint is not feasible');
			}
			var efficientPortfolio = this.computeEfficientPortfolio("volatility", constraintValue + this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a volatility greater than ' + constraintValue + ': the covariance matrix might not be semi-definite positive');
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			var updatedCornerPortfolios = new Array();
			for (i = 0; i < this.cornerPortfolios.length - 1; ++i) {	
				if (this.cornerPortfolios[i][1] > efficientPortfolioRiskTolerance) {
					updatedCornerPortfolios.push(this.cornerPortfolios[i]);
				}
				else {
					break;
				}
			}
			updatedCornerPortfolios.push([efficientPortfolioWeights, efficientPortfolioRiskTolerance]);
			this.cornerPortfolios = updatedCornerPortfolios;
		}
	}
	else {
		throw new Error('unknown constraint type');
	}
};
MeanVarianceEfficientFrontierCla.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	var idx = computeMaximumSharpeRatioCornerPortfolio.call(this, rf);
	var weights_idx = this.cornerPortfolios[idx][0];
	var sr_idx = this.computePortfolioSharpeRatio(weights_idx, rf);
	var candidatePortfolios = [[weights_idx, sr_idx, this.cornerPortfolios[idx][1]]];
	if (idx <= this.cornerPortfolios.length - 2) {
		candidatePortfolios.push( computeMaximumSharpeRatioEfficientSegmentPortfolio.call(this, rf, idx + 1, idx) );
	}
	if (idx >= 1) {
        candidatePortfolios.push( computeMaximumSharpeRatioEfficientSegmentPortfolio.call(this, rf, idx, idx - 1) );
	}
	var compareSharpeRatios = function (a, b) {
		return a[1] - b[1];
	};
	var maxSharpeRatioPortfolio = max_(candidatePortfolios, compareSharpeRatios)[0];
	return [maxSharpeRatioPortfolio[0], maxSharpeRatioPortfolio[2]];
	function computeMaximumSharpeRatioCornerPortfolio(rf) {
		var idx_min = this.cornerPortfolios.length - 1;
		var idx_max = 0;
		if (idx_min == idx_max) {
			return idx_min;
		}
		while (idx_min - idx_max != 1) { 
			var idx_middle = Math.floor((idx_min + idx_max)/2);
			var weights_middle = this.cornerPortfolios[idx_middle][0];
			var sharpeRatio_middle = this.computePortfolioSharpeRatio(weights_middle, rf);

			var idx_middle_p = idx_middle + 1; 
			var weights_middle_p = this.cornerPortfolios[idx_middle_p][0];
			var sharpeRatio_middle_p = this.computePortfolioSharpeRatio(weights_middle_p, rf);
			if (sharpeRatio_middle > sharpeRatio_middle_p) {
				idx_min = idx_middle;
			}
			else if (sharpeRatio_middle	< sharpeRatio_middle_p) {
				idx_max = idx_middle;
			}
			else {
				idx_min = idx_middle_p;		
				idx_max = idx_middle;

				break;
			}
		}
		return idx_min;
	}
	function computeMaximumSharpeRatioEfficientSegmentPortfolio(rf, idx_min, idx_max) {
		var cornerPortfolios = this.cornerPortfolios;
		var weights_min = cornerPortfolios[idx_min][0];
		var weights_max = cornerPortfolios[idx_max][0];
		var sr_min = this.computePortfolioSharpeRatio(weights_min, rf);
		var return_min = this.computePortfolioReturn(weights_min);
		var volatility_min = this.computePortfolioVolatility(weights_min);
		var variance_min = volatility_min * volatility_min;
		
		var sr_max = this.computePortfolioSharpeRatio(weights_max, rf);
		var return_max = this.computePortfolioReturn(weights_max);
		var volatility_max = this.computePortfolioVolatility(weights_max);
		var variance_max = volatility_max * volatility_max;
		var return_min_m_max = return_min - return_max;
		var return_max_m_rf = return_max - rf;
		var a = return_min_m_max * return_min_m_max;
		var b = 2 * return_min_m_max * return_max_m_rf;
		var c = return_max_m_rf * return_max_m_rf;
		
		var variance_cross = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, weights_min), weights_max); // <Sigma*w_min/w_max>
		var d = variance_min + variance_max - 2 * variance_cross; // always >= 0, by semi-definite positivity of the covariance matrix
		var e = -2 * (variance_max - variance_cross); // 
		var f = variance_max; //always > 0
		var aa = a*e - b*d;
		var bb = 2*(a*f - c*d);
		var cc = b*f - c*e;
		var bb_p = bb/2; // reduced discriminant
		var sign_bb_p = (bb_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for bb_p == 0 this returns 1
		var disc = bb_p*bb_p - aa*cc;
		if (disc < 0) {
			throw new Error('internal error, the covariance matrix might not be semi-definite positive');
		}
		var qq = -(bb_p + sign_bb_p * Math.sqrt(disc));
		var t1 = qq/aa;
		var t2 = cc/qq;
		var candidateSharpeRatios = [[weights_min, sr_min, cornerPortfolios[idx_min][1]], [weights_max, sr_max, cornerPortfolios[idx_max][1]]]; // t = 0 and t = 1 portfolios are always present
		
		if (t1 > 0 && t1 < 1) { // t1 belongs to ]0,1[
			var weights_t1 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t1*weights_min.getValue(i, 1) + (1-t1)*weights_max.getValue(i, 1); 
										})
			var sr_t1 = this.computePortfolioSharpeRatio(weights_t1, rf);
			var lambda_t1 = t1*cornerPortfolios[idx_min][1] + (1-t1)*cornerPortfolios[idx_max][1];

			candidateSharpeRatios.push([weights_t1, sr_t1, lambda_t1]);
		}

		if (t2 > 0 && t2 < 1) { // t2 belongs to ]0,1[
			var weights_t2 = Matrix_.fill(weights_min.nbRows, 1, 
										function(i,j) { 
											return t2*weights_min.getValue(i, 1) + (1-t2)*weights_max.getValue(i, 1); 
										})
			var sr_t2 = this.computePortfolioSharpeRatio(weights_t2, rf);
			var lambda_t2 = t2*cornerPortfolios[idx_min][1] + (1-t2)*cornerPortfolios[idx_max][1];
			
			candidateSharpeRatios.push([weights_t2, sr_t2, lambda_t2]);
		}
		var compareSharpeRatios = function (a, b) {
			return a[1] - b[1];
		};
		return max_(candidateSharpeRatios, compareSharpeRatios)[0];
	}
};
function MeanVarianceEfficientFrontierGsmo(mu, sigma, opt) {
	MeanVarianceEfficientFrontier.call(this, mu, sigma, opt);
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	this.epsGsmo = opt.optimizationMethodParams.epsGsmo;
	if (this.epsGsmo == undefined) {
		this.epsGsmo = 1e-6;
	}
	this.maxIterationsGsmo = opt.optimizationMethodParams.maxIterGsmo;
	if (this.maxIterationsGsmo == undefined) {
		this.maxIterationsGsmo = 10000;
	}
	var l = computeMinimumRiskTolerancePortfolio.call(this);
	this.lowestRiskTolerance = l[1];
	this.lowestRiskTolerancePortfolio = l[0];

	var h = computeMaximumRiskTolerancePortfolio.call(this);
	this.highestRiskTolerance = h[1];
	this.highestRiskTolerancePortfolio = h[0];
	function computeMaximumRiskTolerancePortfolio() {
		var nbAssets;
		var mu;
		var sigma;
		var lowerBounds;
		var upperBounds;
		if (this.fullInvestment == true) {
			nbAssets = this.nbAssets;
			mu = this.mu;
			sigma = this.sigma;
			lowerBounds = this.lowerBounds;
			upperBounds = this.upperBounds;
		}
		else {
			nbAssets = this.alteredNbAssets;
			mu = this.alteredMu;
			sigma = this.alteredSigma;
			lowerBounds = this.alteredLowerBounds;
			upperBounds = this.alteredUpperBounds;
		}	
		var maxReturnSolution = simplexLpSolve_(mu.elemMap(function(i,j,val) { return -val;}), lowerBounds, upperBounds);
		var maxReturn = -maxReturnSolution[1];
		var Q = sigma;
		var b = Matrix_.ones(nbAssets, 1);
		var r = 1;
		var l = lowerBounds;
		var u = upperBounds;
		var epsSearch =  1e-12;
		var nbIterSearch = 0;
		var maxNbIterSearch = 54; // The default is taken to be 54, because 2^54 is quite close to an already unreasonable high value
		var rt = 0.5;
		do {
			++nbIterSearch;
			if (nbIterSearch > maxNbIterSearch) {
				throw new Error('internal error: maximum number of iterations reached when searching for the efficient portfolio with maximum return');
			}
			rt = 2 * rt;
			var p = Matrix_.ax(-rt, mu);
			var efficientSol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
			var ret = this.computePortfolioReturn(efficientSol[0], mu, sigma);
		}
		while ( Math.abs(ret - maxReturn) > epsSearch );
		var portfolio = efficientSol[0];
		if (this.fullInvestment == false) {
			portfolio = Matrix_.fill(nbAssets, 1, 
									function(i, j) { 
										if (i <= nbAssets) {
											return portfolio.getValue(i, 1);
										}
									});
		}
		return [portfolio, rt];
	}
	function computeMinimumRiskTolerancePortfolio() {
		var nbAssets;
		var mu;
		var sigma;
		var lowerBounds;
		var upperBounds;
		if (this.fullInvestment == true) {
			nbAssets = this.nbAssets;
			mu = this.mu;
			sigma = this.sigma;
			lowerBounds = this.lowerBounds;
			upperBounds = this.upperBounds;
		}
		else {
			nbAssets = this.alteredNbAssets;
			mu = this.alteredMu;
			sigma = this.alteredSigma;
			lowerBounds = this.alteredLowerBounds;
			upperBounds = this.alteredUpperBounds;
		}
		var Q = sigma;
		var p = Matrix_.zeros(nbAssets, 1);
		var b = Matrix_.ones(nbAssets, 1);
		var r = 1;
		var l = lowerBounds;
		var u = upperBounds;
		var minVolatilitySolution = qpsolveGSMO_(Q, p, b, r, l, u, {eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
		var minVolatility = this.computePortfolioVolatility(minVolatilitySolution[0]);
		var epsSearch =  1e-12;
		var rt = 2;
		do {
			rt = rt / 2;
			var p = Matrix_.ax(-rt, mu);
			var efficientSol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
			var vol = this.computePortfolioVolatility(efficientSol[0]);
		}
		while ( Math.abs(vol - minVolatility) > epsSearch );
		var portfolio = efficientSol[0];
		if (this.fullInvestment == false) {
			portfolio = Matrix_.fill(nbAssets, 1, 
									function(i, j) { 
										if (i <= nbAssets) {
											return portfolio.getValue(i, 1);
										}
									});
		}
		return [portfolio, rt];
	}

}
MeanVarianceEfficientFrontierGsmo.prototype = Object.create(MeanVarianceEfficientFrontier.prototype);
MeanVarianceEfficientFrontierGsmo.prototype.constructor = MeanVarianceEfficientFrontierGsmo;

MeanVarianceEfficientFrontierGsmo.prototype.getHighestRiskTolerancePortfolio = function(x) {
	var portfolio = this.highestRiskTolerancePortfolio;
	if (this.fullInvestment == false) {
		var that = this;
		portfolio = Matrix_.fill(this.nbAssets, 1, 
								function(i, j) { 
									if (i <= that.nbAssets) {
										return portfolio.getValue(i, 1);
									}
								});
	}
	return portfolio;
};
MeanVarianceEfficientFrontierGsmo.prototype.getHighestRiskTolerance = function(x) {
	return this.highestRiskTolerance;
};
MeanVarianceEfficientFrontierGsmo.prototype.getLowestRiskTolerancePortfolio = function(x) {
	var portfolio = this.lowestRiskTolerancePortfolio;
	if (this.fullInvestment == false) {
		var that = this;
		portfolio = Matrix_.fill(this.nbAssets, 1, 
								function(i, j) { 
									if (i <= that.nbAssets) {
										return portfolio.getValue(i, 1);
									}
								});
	}
	return portfolio;
};
MeanVarianceEfficientFrontierGsmo.prototype.getLowestRiskTolerance = function(x) {
	return this.lowestRiskTolerance;
};
MeanVarianceEfficientFrontierGsmo.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	if (constraintType === undefined || constraintType === null) {
		throw new Error('internal error: missing constraint type');
	}

	var constraintFct;
	var that = this;
	if (constraintType == "return") {
		constraintFct = function (portfolio) { return that.computePortfolioReturn(portfolio); };
	}
	else if (constraintType == "volatility") {
		constraintFct = function (portfolio) { return that.computePortfolioVolatility(portfolio); };
	}
	else if (constraintType == "riskTolerance") {
	}
	else {
		throw new Error('internal error: unknown constraint type');
	}
	
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('internal error: missing constraint value');
	}
	var nbAssets;
	var mu;
	var sigma;
	var lowerBounds;
	var upperBounds;
	var eps = this.epsEfficientPortfolioComputation;
	if (this.fullInvestment == true) {
		nbAssets = this.nbAssets;
		mu = this.mu;
		sigma = this.sigma;
		lowerBounds = this.lowerBounds;
		upperBounds = this.upperBounds;
	}
	else {
		nbAssets = this.alteredNbAssets;
		mu = this.alteredMu;
		sigma = this.alteredSigma;
		lowerBounds = this.alteredLowerBounds;
		upperBounds = this.alteredUpperBounds;
	}
	var Q = sigma;
	var b = Matrix_.ones(nbAssets, 1);
	var r = 1;
	var l = lowerBounds;
	var u = upperBounds;
	var portfolioWeights;
	var riskTolerance;
	if (constraintType === "return" || constraintType === "volatility") {
		var weights_min = this.getLowestRiskTolerancePortfolio();
		var constraintFct_min = constraintFct(weights_min);
		
		var weights_max = this.getHighestRiskTolerancePortfolio()
		var constraintFct_max = constraintFct(weights_max);
		if (constraintValue > constraintFct_max + eps || constraintValue < constraintFct_min - eps) {
			return [];
		}
		if (Math.abs(constraintValue - constraintFct_min) <= eps) {
			portfolioWeights = weights_min;
			riskTolerance = this.getLowestRiskTolerance();
		}
		else if (Math.abs(constraintValue - constraintFct_max) <= eps) {
			portfolioWeights = weights_max;
			riskTolerance = this.getHighestRiskTolerance();
		}
		else {
			var that = this;
			riskTolerance = bisection_(function (rt) { 
											var p = Matrix_.ax(-rt, mu);
											var sol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: that.epsGsmo, maxIter: that.maxIterationsGsmo});
											
											portfolioWeights = sol[0];
											return constraintFct(portfolioWeights) - constraintValue; 
										}, 
										this.getLowestRiskTolerance(), this.getHighestRiskTolerance());
		}
	}		
	else if (constraintType === "riskTolerance") {
		var rt = constraintValue;
		var p = Matrix_.ax(-rt, mu);
		var sol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: this.epsGsmo, maxIter: this.maxIterationsGsmo});
		
		portfolioWeights = sol[0];
		riskTolerance = rt;
	}
	else {
		throw new Error('internal error: unknown constraint type');
	}
	if (this.fullInvestment == false) {
		var that = this;
		portfolioWeights = Matrix_.fill(this.nbAssets, 1, 
										function(i, j) { 
											if (i <= that.nbAssets) {
												return portfolioWeights.getValue(i, 1);
											}
										});
	}	
	
	return [portfolioWeights, riskTolerance];
};
MeanVarianceEfficientFrontierGsmo.prototype.restrict = function(constraintType, constraintValue) {
	if (constraintType === undefined || constraintType === null) {
		throw new Error('missing constraint type');
	}
	if (constraintValue === undefined || constraintValue === null) {
		throw new Error('missing constraint value');
	}
	if (constraintType == "minReturn") {	
		var lowestReturn = this.getLowestReturn();
		if (lowestReturn < constraintValue) {
			var highestReturn = this.getHighestReturn();
			if (highestReturn < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum return constraint is not feasible');
			}
			var efficientPortfolio = this.computeEfficientPortfolio("return", constraintValue + this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a return greater than ' + constraintValue);
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			this.lowestRiskTolerancePortfolio = efficientPortfolioWeights;
			this.lowestRiskTolerance = efficientPortfolioRiskTolerance;
			if (this.lowestRiskTolerance >= this.highestRiskTolerance) {
				this.highestRiskTolerancePortfolio = efficientPortfolioWeights;
				this.highestRiskTolerance = efficientPortfolioRiskTolerance;
			}			
		}
	}
	else if (constraintType == "minVolatility") {	
		var lowestVolatility = this.getLowestVolatility();
		if (lowestVolatility < constraintValue) {
			var highestVolatility = this.getHighestVolatility();
			if (highestVolatility < constraintValue) {
				throw new Error('impossible to restrict the efficient frontier: the minimum volatility constraint is not feasible');
			}
			var efficientPortfolio = this.computeEfficientPortfolio("volatility", constraintValue + this.epsEfficientPortfolioComputation);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no efficient portfolio with a volatility greater than ' + constraintValue + ': the covariance matrix might not be semi-definite positive');
			}
			var efficientPortfolioWeights = efficientPortfolio[0];
			var efficientPortfolioRiskTolerance = efficientPortfolio[1];
			this.lowestRiskTolerancePortfolio = efficientPortfolioWeights;
			this.lowestRiskTolerance = efficientPortfolioRiskTolerance;
			if (this.lowestRiskTolerance >= this.highestRiskTolerance) {
				this.highestRiskTolerancePortfolio = efficientPortfolioWeights;
				this.highestRiskTolerance = efficientPortfolioRiskTolerance;
			}			
		}
	}
	else {
		throw new Error('unknown constraint type');
	}
};
MeanVarianceEfficientFrontierGsmo.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	var nbAssets;
	var mu;
	var sigma;
	var lowerBounds;
	var upperBounds;
	var eps = this.epsEfficientPortfolioComputation;
	if (this.fullInvestment == true) {
		nbAssets = this.nbAssets;
		mu = this.mu;
		sigma = this.sigma;
		lowerBounds = this.lowerBounds;
		upperBounds = this.upperBounds;
	}
	else {
		nbAssets = this.alteredNbAssets;
		mu = this.alteredMu;
		sigma = this.alteredSigma;
		lowerBounds = this.alteredLowerBounds;
		upperBounds = this.alteredUpperBounds;
	}
	var Q = sigma;
	var b = Matrix_.ones(nbAssets, 1);
	var r = 1;
	var l = lowerBounds;
	var u = upperBounds;
	var portfolioWeights;
	var that = this;
	var riskTolerance = goldenSectionSearch_(function(rt) { 
												  var p = Matrix_.ax(-rt, mu);
												  var sol = qpsolveGSMO_(Q, p, b, r, l, u, {eps: that.epsGsmo, maxIter: that.maxIterationsGsmo});
														
												  portfolioWeights = sol[0];
												  return -that.computePortfolioSharpeRatio(portfolioWeights, rf); 
											 }, 
											 this.getLowestRiskTolerance(), this.getHighestRiskTolerance());
	if (this.fullInvestment == false) {
		var that = this;
		portfolioWeights = Matrix_.fill(this.nbAssets, 1, 
										function(i, j) { 
											if (i <= that.nbAssets) {
												return portfolioWeights.getValue(i, 1);
											}
										});
	}	
	return [portfolioWeights, riskTolerance[0]];
};
function MeanVarianceEfficientFrontier(mu, sigma, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	this.epsPortfolioVolatility = 1e-14; // the tolerance for numerically zero volatility
	this.epsEfficientPortfolioComputation = 1e-10; // the tolerance for numerically searching for an efficient portfolio with a given constraint value
	
	this.mu = new Matrix_(mu);
	this.sigma = new Matrix_(sigma);
	this.nbAssets = this.mu.nbRows;
	
	this.fullInvestment = opt.constraints.fullInvestment == undefined ? true : opt.constraints.fullInvestment;
	this.lowerBounds = opt.constraints.minWeights == undefined ? Matrix_.zeros(this.nbAssets, 1) : new Matrix_(opt.constraints.minWeights);
	this.upperBounds = opt.constraints.maxWeights == undefined  ? Matrix_.ones(this.nbAssets, 1) : new Matrix_(opt.constraints.maxWeights);
	if (this.fullInvestment == false) {
		var that = this;
		this.alteredMu = Matrix_.fill(this.nbAssets + 1, 1, 
								function(i, j) { 
									if (i <= that.nbAssets) {
										return that.mu.getValue(i, 1);
									}
									else {
										return 0; // risk free asset, no return
									}
								});
		this.alteredSigma = Matrix_.fill(this.nbAssets + 1, this.nbAssets + 1, 
								function(i, j) { 
									if (i <= that.nbAssets && j <= that.nbAssets) {
										return that.sigma.getValue(i, j);
									}
									else {
										return 0; // risk free asset, no variance/covariance
									}
								});
								
		this.alteredLowerBounds = Matrix_.fill(this.nbAssets + 1, 1, 
										function(i, j) { 
											if (i <= that.nbAssets) {
												return that.lowerBounds.getValue(i, 1);
											}
											else {
												return 0; // risk free asset, no min weight
											}
										});	

		this.alteredUpperBounds = Matrix_.fill(this.nbAssets + 1, 1, 
										function(i, j) { 
											if (i <= that.nbAssets) {
												return that.upperBounds.getValue(i, 1);
											}
											else {
												return 1; // risk free asset, no max weight
											}
										});	
		this.alteredNbAssets = this.nbAssets + 1;
	}
	if (this.fullInvestment == true) {
		simplexEmptinessCheck_(this.nbAssets, this.lowerBounds.toArray(), this.upperBounds.toArray());	
	}
	else {
		if (this.fullInvestment == false) {
			simplexEmptinessCheck_(this.alteredNbAssets, this.alteredLowerBounds.toArray(), this.alteredUpperBounds.toArray());	
		}
	}
};

MeanVarianceEfficientFrontier.prototype.getHighestReturnPortfolio = function(x) {
	return this.getHighestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getHighestReturn = function(x) {
	return this.computePortfolioReturn(this.getHighestReturnPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getLowestReturnPortfolio = function(x) {
	return this.getLowestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getLowestReturn = function(x) {
	return this.computePortfolioReturn(this.getLowestReturnPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getHighestVolatilityPortfolio = function(x) {
	return this.getHighestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getHighestVolatility = function(x) {
	return this.computePortfolioVolatility(this.getHighestVolatilityPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getLowestVolatilityPortfolio = function(x) {
	return this.getLowestRiskTolerancePortfolio();
};
MeanVarianceEfficientFrontier.prototype.getLowestVolatility = function(x) {
	return this.computePortfolioVolatility(this.getLowestVolatilityPortfolio());
};
MeanVarianceEfficientFrontier.prototype.getHighestRiskTolerancePortfolio = function(x) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getHighestRiskTolerance = function(x) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getLowestRiskTolerancePortfolio = function(x) {
	throw new Error('internal error: unction is not implemented');
};
MeanVarianceEfficientFrontier.prototype.getLowestRiskTolerance = function(x) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.computePortfolioReturn = function(x) {
	if (this.fullInvestment == true) {
		return Matrix_.vectorDotProduct(this.mu, x);
	}
	else {
		if (x.nbRows == this.alteredNbAssets) {
			return Matrix_.vectorDotProduct(this.alteredMu, x);
		}
		else if (x.nbRows == this.nbAssets) {
			return Matrix_.vectorDotProduct(this.mu, x);
		}
		else {
			throw new Error("internal error: unexpected portfolio number of weights: " + x.nbRows);
		}
	}
};
MeanVarianceEfficientFrontier.prototype.computePortfolioVolatility = function(x) {
	var sigma_x_x;
	if (this.fullInvestment == true) {
		sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, x), x);		
	}
	else {
		if (x.nbRows == this.alteredNbAssets) {
			sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(this.alteredSigma, x), x);	
		}
		else if (x.nbRows == this.nbAssets) {
			sigma_x_x = Matrix_.vectorDotProduct(Matrix_.xy(this.sigma, x), x);		
		}
		else {
			throw new Error("internal error: unexpected portfolio number of weights: " + x.nbRows);
		}
	} 
	if (Math.abs(sigma_x_x) <= this.epsPortfolioVolatility) {
		sigma_x_x = 0;
	}
	else if (sigma_x_x < 0 && sigma_x_x < -this.epsPortfolioVolatility) {
		throw new Error('internal error: negative volatility, covariance matrix might not be semi-definite positive');
	}
	var s_x = Math.sqrt(sigma_x_x);
	return s_x;
};
MeanVarianceEfficientFrontier.prototype.computePortfolioSharpeRatio = function(x, rf) {
	if (rf == undefined || rf == null) {
		throw new Error('internal error: missing risk free rate');
	}
	var ret = this.computePortfolioReturn(x);
	var excessRet = ret - rf;
	var vol = Math.max(this.epsPortfolioVolatility, this.computePortfolioVolatility(x));
	if (vol == 0) {
		throw new Error('internal error: null volatility when computing the Sharpe ratio');
	}
	var sharpeRatio = excessRet/vol;
	return sharpeRatio;
};
MeanVarianceEfficientFrontier.prototype.computeEfficientPortfolio = function(constraintType, constraintValue) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.restrict = function(constraintType, constraintValue) {
	throw new Error('internal error: function is not implemented');
};
MeanVarianceEfficientFrontier.prototype.computeMaximumSharpeRatioEfficientPortfolio = function(rf) {
	throw new Error('internal error: function is not implemented');
};
function meanVarianceOptimizationWeights(mu, sigma, opt) {	
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'gsmo';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo') {
		throw new Error('unsupported optimisation method');
	}
	var returnConstraint = opt.constraints["return"]; // .return does not work within Google Script
	var volatilityConstraint = opt.constraints.volatility;
	var maxVolatilityConstraint = opt.constraints.maxVolatility;
	var riskToleranceConstraint = opt.constraints.riskTolerance;
	if (returnConstraint === undefined && volatilityConstraint === undefined &&
		maxVolatilityConstraint === undefined && riskToleranceConstraint == undefined) {
		throw new Error('missing return, volatility or risk tolerance constraints');
	}
	if ( (returnConstraint !== undefined && volatilityConstraint !== undefined) ||
		 (returnConstraint !== undefined && maxVolatilityConstraint !== undefined) ) {
		throw new Error('simultaneous return and volatility constraints');
	}
	if ( riskToleranceConstraint !== undefined && 
	     (returnConstraint !== undefined || maxVolatilityConstraint !== undefined || volatilityConstraint !== undefined) ) {
		throw new Error('simultaneous risk tolerance and return or volatility constraints');
	}
	var efficientFrontier;
	if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	var efficientPortfolioWeights;
	if (returnConstraint !== undefined) {
		var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("return", returnConstraint);
		if (efficientPortfolio.length == 0) {
			throw new Error('no matching efficient portfolio with a return equal to ' + returnConstraint);
		}
		else {
			efficientPortfolioWeights = efficientPortfolio[0];
		}
	}
	else if (volatilityConstraint !== undefined) {
		var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("volatility", volatilityConstraint);
		if (efficientPortfolio.length == 0) {
			throw new Error('no matching efficient portfolio with a volatility equal to ' + volatilityConstraint);
		}
		else {
			efficientPortfolioWeights = efficientPortfolio[0];
		}
	}
	else if (riskToleranceConstraint !== undefined) {
		var highestRiskTolerancePortfolioWeights = efficientFrontier.getHighestRiskTolerancePortfolio();
		var highestRiskTolerance = efficientFrontier.getHighestRiskTolerance();
		var lowestRiskTolerancePortfolioWeights = efficientFrontier.getLowestRiskTolerancePortfolio();
		var lowestRiskTolerance = efficientFrontier.getLowestRiskTolerance();
		var eps = efficientFrontier.epsEfficientPortfolioComputation;
		if (riskToleranceConstraint > highestRiskTolerance - eps) {
			efficientPortfolioWeights = highestRiskTolerancePortfolioWeights;
		}
		else if (riskToleranceConstraint < lowestRiskTolerance + eps) {
			efficientPortfolioWeights = lowestRiskTolerancePortfolioWeights;
		}
		else {
			var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("riskTolerance", riskToleranceConstraint);
			if (efficientPortfolio.length == 0) {
				throw new Error('internal error: no matching efficient portfolio with a risk tolerance equal to ' + riskToleranceConstraint);
			}
			else {
				efficientPortfolioWeights = efficientPortfolio[0];
			}
		}																						  
	}
	else if (maxVolatilityConstraint !== undefined) {
		var highestVolatilityPortfolioWeights = efficientFrontier.getHighestVolatilityPortfolio();
		var highestVolatility = efficientFrontier.getHighestVolatility();
		var eps = efficientFrontier.epsEfficientPortfolioComputation;
		if (maxVolatilityConstraint > highestVolatility - eps) {
			efficientPortfolioWeights = highestVolatilityPortfolioWeights;
		}
		else {
			var efficientPortfolio = efficientFrontier.computeEfficientPortfolio("volatility", maxVolatilityConstraint);
			if (efficientPortfolio.length == 0) {
				throw new Error('no matching efficient portfolio with a volatility lower than or equal to ' + maxVolatilityConstraint);
			}
			else {
				efficientPortfolioWeights = efficientPortfolio[0];
			}
		}
	}	
	return efficientPortfolioWeights.toArray();
}
function meanVarianceEfficientFrontierNearestWeights(inputWeights, mu, sigma, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'gsmo';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo') {
		throw new Error('unsupported optimisation method');
	}
	var inputWeights = new Matrix_(inputWeights);
	var efficientFrontier;
	var efficientFrontierPortfolios;
	if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
		 efficientFrontierPortfolios = efficientFrontier.getCornerPortfolios();
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
		var nbEfficientPortfolios = opt.optimizationMethodParams.nbPortfoliosGsmo;
		if (nbEfficientPortfolios == undefined) {
			nbEfficientPortfolios = 100;
		}
		efficientFrontierPortfolios = new Array(nbEfficientPortfolios);
		
		var rt_min = efficientFrontier.getLowestRiskTolerance();
		var rt_max = efficientFrontier.getHighestRiskTolerance();	
		var delta_rt = (rt_max - rt_min)/(nbEfficientPortfolios - 1);	
		for (var i = 0; i < nbEfficientPortfolios; ++i) {
			var rt_i = rt_min + i * delta_rt;
			var weights = efficientFrontier.computeEfficientPortfolio("riskTolerance", rt_i);
			if (weights.length == 0) {
				throw new Error('internal error: no matching efficient portfolio with a risk tolerance of ' + rt_i);
			}
			else {
				efficientFrontierPortfolios[i] = weights[0];
			}
		}
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	var weights;
	var minDist = Number.POSITIVE_INFINITY;
	for (var i = 0; i < efficientFrontierPortfolios.length - 1; ++i) {
		var w_e = efficientFrontierPortfolios[i];
		var w_b = efficientFrontierPortfolios[i+1];
		var proj = lineSegmentEuclidianProjection_(inputWeights, w_b, w_e);
		var dist = Matrix_.xmy(inputWeights, new Matrix_(proj)).vectorNorm('two');
		if (dist <= minDist) {
			weights = proj;
			minDist = dist;
		}
	}
	return weights;
}
function meanVarianceEfficientFrontierPortfolios(mu, sigma, opt) {	
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'gsmo';
	}
	if (optimizationMethod != 'critical-line' && optimizationMethod != 'gsmo') {
		throw new Error('unsupported optimisation method');
	}

	var nbPortfolios = opt.nbPortfolios || 100;

	var discretizationType = opt.discretizationType;
	if (discretizationType == undefined) {
		discretizationType = "return";
	}
	if (discretizationType != "return" && discretizationType != "volatility" && discretizationType != "riskTolerance") {
		throw new Error('unsupported discretization type');
	}
	var efficientFrontier;
	if (optimizationMethod == "critical-line") {
		 efficientFrontier = new MeanVarianceEfficientFrontierCla(mu, sigma, opt);
	}
	else if (optimizationMethod == "gsmo") {
		efficientFrontier = new MeanVarianceEfficientFrontierGsmo(mu, sigma, opt);
	}
	else {
		throw new Error('internal error: unsupported optimisation method');
	}
	var efficientFrontierPortfolios = new Array(nbPortfolios);
	var t_min;
	var t_max;
	if (discretizationType == "return") {
		t_min = efficientFrontier.getLowestReturn();
		t_max = efficientFrontier.getHighestReturn();
	}
	else if (discretizationType == "volatility") {
		t_min = efficientFrontier.getLowestVolatility();
		t_max = efficientFrontier.getHighestVolatility();
	}
	else if (discretizationType == "riskTolerance") {
		t_min = efficientFrontier.getLowestRiskTolerance();
		t_max = efficientFrontier.getHighestRiskTolerance();
	}
	else {
		throw new Error('internal error: unsupported discretization type');
	}
	var delta_t = nbPortfolios == 1 ? -1 : (t_max - t_min)/(nbPortfolios - 1);
	var t_start = nbPortfolios == 1 ? (t_min + t_max)/2 : t_min;
	for (var i = 0; i < nbPortfolios; ++i) {
		var t_i = t_start + i * delta_t;
		var weights = efficientFrontier.computeEfficientPortfolio(discretizationType, t_i);
		if (weights.length == 0) {
			throw new Error('internal error: no matching efficient portfolio with a constraint value ' + discretizationType + ' equal to ' + t_i);
		}
		else {
			weights = weights[0];
			var ret = efficientFrontier.computePortfolioReturn(weights);
			var vol = efficientFrontier.computePortfolioVolatility(weights);
			efficientFrontierPortfolios[i] = [weights.toArray(), ret, vol];
		}
	}
	return efficientFrontierPortfolios;
}
function minimaxWeights (assetsReturns, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	var fullInvestmentContraint = true;
	if (opt.constraints.fullInvestment !== undefined) {
		fullInvestmentContraint = opt.constraints.fullInvestment;
	}
	var outputMinimumPortfolioReturn = false || opt.outputMinimumPortfolioReturn; 
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;
	var c = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 0 : -1; }); // c = [0,...,0,-1]
	var Ae = null;
	var be = null;
	if (fullInvestmentContraint) {
		Ae = Matrix_.fill(1, nbAssets + 1, function(i,j) { return j <= nbAssets ? 1 : 0; }); // Ae = [1,...,1,0]
		be = Matrix_.ones(1,1); // be = [1]
	}
	var Ai = Matrix_.fill(nbPeriods + (fullInvestmentContraint ? 0 : 1), nbAssets + 1, 
						  function(i,j) { 
								if (i <= nbPeriods) { return j <= nbAssets ? -assetsReturns[j-1][i-1] : 1; }
								else if (i == nbPeriods + 1) { return j <= nbAssets ? 1 : 0; }
						  }); // Ai = [[-ret11,...,-retN1,1], ..., [-ret1T,...,-retNT,1] (optional: , [1,...,1,0])]
	var bi = Matrix_.fill(nbPeriods + (fullInvestmentContraint ? 0 : 1), 1, 
						  function(i,j) { 
							  if (i <= nbPeriods) { return 0; }
							  else if (i == nbPeriods + 1) { return 1; }
						  }); // bi = [0, ..., 0 (optional: , 1)]
	var lb = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 0 : -Infinity; }); // lb = [0,...,0, -Infinity]
	var ub = Matrix_.fill(nbAssets + 1, 1, function(i,j) { return i <= nbAssets ? 1 : Infinity; });  // ub = [1,...,1, Infinity]
	var lpSolution = lpsolvePDHG_(Ae, be, Ai, bi, c, lb, ub, {maxIter: -1});
	var sol = lpSolution[0];
	var weights = sol.toArray(function(i, j, val) { return i != nbAssets + 1; });
	var minPortfolioReturn = sol.data[nbAssets];
	if (outputMinimumPortfolioReturn === true) {
		return [weights, minPortfolioReturn];
	}
	else {
		return weights;
	}
}
function minimumCorrelationWeights (sigma, opt) {
	var sigma = new Matrix_(sigma).toCovarianceMatrix();
	var variances = sigma.getVariances();
	var invStddevs = variances.elemMap(function(i,j,val) { return 1/Math.sqrt(val); });
	var rho = sigma.getCorrelationMatrix();
	var nbAssets = rho.nbRows;
	if (nbAssets <= 2) {
		return inverseVolatilityWeights(variances, opt);
	}
	var elemRho = rho.toArray(function(i, j, val) {
		return j > i;
	});
	var elementsMean = mean_(elemRho);
	var elementsStddev = sampleStddev_(elemRho);
	var adjustedRho = rho.elemMap(function(i, j, val) { 
			if (i == j) {
				return 0; // Necessary for step 6
			}
			else {
				return 1 - normcdf_((val - elementsMean)/elementsStddev);
			}
		});
	var rowsElements = adjustedRho.toRowArray(function(i, j, val) {
		return i != j;
	});
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	var ranks = rank_(rowsAverages, 0);
	var weights = new Matrix_(ranks, 1).normalize();
	weights = Matrix_.xy(adjustedRho, weights).normalize();
	weights = Matrix_.vectorHadamardProduct(weights, invStddevs).normalize();
	return weights.toArray();
}
function minimumTrackingErrorWeights (assetsReturns, benchmarkReturns, opt) {
	if (opt === undefined) {
		opt = { constraints: {} };
	}
	if (opt.constraints  === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	if (opt.optimizationMethodParams.eps === undefined) {
		opt.optimizationMethodParams.eps = 1e-04;
	}
	if (opt.optimizationMethodParams.maxIter === undefined) {
		opt.optimizationMethodParams.maxIter = 10000;
	}
	if (opt.constraints.minNbAssets === undefined && opt.constraints.maxNbAssets) {
		opt.constraints.minNbAssets = 1;
	}
	if (opt.constraints.maxNbAssets === undefined && opt.constraints.minNbAssets) {
		opt.constraints.maxNbAssets = assetsReturns.length;
	}
	var eps = opt.optimizationMethodParams.eps;
	var maxIterations = opt.optimizationMethodParams.maxIter;

	var fullInvestmentConstraint = true;
	if (opt.constraints.fullInvestment !== undefined) {
		fullInvestmentConstraint = opt.constraints.fullInvestment;
	}
	
	var minNbAssets = opt.constraints.minNbAssets;
	var maxNbAssets = opt.constraints.maxNbAssets;
	var cardinalityConstraints = minNbAssets || maxNbAssets;
	var optimizationMethod;
	if (cardinalityConstraints) {
		optimizationMethod = opt.optimizationMethod;
		
		if (optimizationMethod === undefined) {
			optimizationMethod = 'heuristic';
		}
		
		if (optimizationMethod != 'heuristic' && optimizationMethod != 'exact') {
			throw new Error('unsupported optimisation method');
		}
	}
	var nRounds;
	var nSteps;
	var nDeltas;
	var maxIterationsInitPoint;
	if (optimizationMethod == 'heuristic') {
		nRounds = opt.optimizationMethodParams.nRounds;
		if (nRounds === undefined) {
			nRounds = 3;
		}
		
		nSteps = opt.optimizationMethodParams.nSteps;
		if (nSteps === undefined) {
			nSteps = 5000;
		}

		nDeltas = opt.optimizationMethodParams.nDeltas;
		if (nDeltas === undefined) {
			nDeltas = nSteps;
		}
		
		maxIterationsInitPoint = opt.optimizationMethodParams.maxIterationsInitPoint;
		if (maxIterationsInitPoint === undefined) {
			maxIterationsInitPoint = 10000;
		}
	}
	var nbAssets = assetsReturns.length;
	var nbPeriods = assetsReturns[0].length;
	var benchmarkReturns = new Matrix_(benchmarkReturns);
	var assetsReturns = Matrix_.fill(nbPeriods, nbAssets, function(i,j) { return assetsReturns[j-1][i-1]; });
	var lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
	var upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbAssets) : new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		lowerBounds[i] = opt.constraints.minWeights ? opt.constraints.minWeights[i] : 0;
		upperBounds[i] = opt.constraints.maxWeights ? opt.constraints.maxWeights[i] : 1;
	}
	function computeMinimumTrackingErrorVolatilityPortfolio(assetsReturns, benchmarkReturns, lowerBounds, upperBounds, fullInvestmentConstraint) {	
		function f(w) {
			var te = Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns).vectorNorm('two');
			
			return 0.5 * te * te;
		}
		function gradf(w) {
			var gte = Matrix_.txy(assetsReturns, Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns));
			
			return gte;
		}
		function g(w) {
			if (fullInvestmentConstraint == true) {
				return simplexCharacteristicFunction_(w.data, lowerBounds, upperBounds);
			}
			else {
				return fullSimplexCharacteristicFunction_(w.data, lowerBounds, upperBounds);
			}
		}
		function proxg(w) {
			if (fullInvestmentConstraint == true) {
				return new Matrix_(simplexEuclidianProjection_(w.data, lowerBounds, upperBounds));
			}
			else {
				return new Matrix_(fullSimplexEuclidianProjection_(w.data, lowerBounds, upperBounds));
			}
		}
		var x0;
		if (fullInvestmentConstraint == true) {
			x0 = new Matrix_(simplexEuclidianProjection_(Matrix_.ones(assetsReturns.nbColumns, 1).data, lowerBounds, upperBounds));
		}
		else {
			x0 = new Matrix_(fullSimplexEuclidianProjection_(Matrix_.ones(assetsReturns.nbColumns, 1).data, lowerBounds, upperBounds));
		}
		var sol = ccpsolveFISTA_(f, gradf, g, proxg, x0, {eps: eps, maxIter: maxIterations, maxLine: maxIterations});
		return sol;
	} 
	var weights;
	if (!cardinalityConstraints) {
		weights = computeMinimumTrackingErrorVolatilityPortfolio(assetsReturns, benchmarkReturns, 
																 lowerBounds, upperBounds,
																 fullInvestmentConstraint)[0];
	}
	else {
		if (optimizationMethod == 'heuristic') {
			function f(w) {
				var w = new Matrix_(w);
				
				var te = Matrix_.xmy(Matrix_.xy(assetsReturns, w), benchmarkReturns).vectorNorm('two');
				
				return 0.5 * te * te;
			}
			function neighbourGenerator(x, neighbourGeneratorParameters) {		
				function randomIndex(min, max) { // 0, n --> 0, n-1
					return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
				}
				function random(min, max) { // [min, max[
					return Math.random() * (max - min) + min;
				}
				var l = neighbourGeneratorParameters.lowerBounds;
				var u = neighbourGeneratorParameters.upperBounds;
				var xx = x;
				var nbNonNullAssets = 0;
				var cashPosition = 1;
				for (var i = 0; i < nbAssets; ++i) {
					if ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) {
						++nbNonNullAssets;
						cashPosition -= xx[i];
						
						if (xx[i] < l[i]) {
							throw new Error('internal error: asset present in the portfolio, but strictly below its lower bound');
						}
						
						if (xx[i] > u[i]) {
							throw new Error('internal error: asset present in the portfolio, but strictly above its upper bound');
						}
					}
				}
				cashPosition = Math.max(0, cashPosition);
				
				if (nbNonNullAssets > maxNbAssets) {
					throw new Error('internal error: number of assets strictly greater than the allowed maximum number of assets');
				}
				else if (nbNonNullAssets < minNbAssets) {
					throw new Error('internal error: number of assets strictly lower than the allowed minimum number of assets');
				}				
				var assetsSwitches = []
				for (var i = 0; i < nbAssets; ++i) {
					if ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) { // The asset i is existing in the portfolio and can be partially of fully sold
						var alphaPartialSellMax = xx[i] - l[i]; // -eps is in order to prevent null positions in case there is no lower bound for the asset i, otherwise, the cardinality constraint would be violated
						var alphaPartialSellMin = Math.min(0, alphaPartialSellMax);
						var alphaFullSell = xx[i];
						if (!fullInvestmentConstraint) {
							assetsSwitches.push([i, -1, alphaPartialSellMin, alphaPartialSellMax]);
							
							if (nbNonNullAssets > minNbAssets) {
								assetsSwitches.push([i, -1, alphaFullSell]);
							}
						}						
						for (var j = 0; j < nbAssets; ++j) {
							if ((xx[j] < l[j] || xx[j] == 0) && xx[j] < u[j]) {
								var alphaBuyMin = l[j];
								var alphaBuyMax = Math.max(alphaBuyMin, u[j] - xx[j]);
								if (nbNonNullAssets < maxNbAssets) {
									if (alphaPartialSellMin <= alphaBuyMin && alphaBuyMin <= alphaPartialSellMax) {
										var alphaMin = alphaBuyMin;
										var alphaMax = Math.min(alphaBuyMax, alphaPartialSellMax);
										assetsSwitches.push([i, j, alphaMin, alphaMax]);
									}
								}
								if (alphaBuyMin <= alphaFullSell && alphaFullSell <= alphaBuyMax) {
									assetsSwitches.push([i, j, alphaFullSell]);
								}
							}
							if (((xx[j] > 0 && l[j] == 0) || (xx[j] >= l[j] && l[j] > 0)) && xx[j] < u[j]) {
								var alphaBuyMax = u[j] - xx[j];
								var alphaBuyMin = Math.min(0, alphaBuyMax);
								if (nbNonNullAssets > minNbAssets) {
									if (alphaBuyMin <= alphaFullSell && alphaFullSell <= alphaBuyMax) {
										assetsSwitches.push([i, j, alphaFullSell]);
									}
								}
								if (alphaPartialSellMin <= alphaBuyMin && alphaBuyMin <= alphaPartialSellMax) {
									var alphaMin = alphaBuyMin;
									var alphaMax = Math.min(alphaBuyMax, alphaPartialSellMax);
									assetsSwitches.push([i, j, alphaMin, alphaMax]);
								}
							}
						}
					}
				}
				if (!fullInvestmentConstraint) {
					for (var i = 0; i < nbAssets; ++i) {
						if ((xx[i] > 0 && l[i] == 0) || (xx[i] >= l[i] && l[i] > 0)) { // The asset i is existing in the portfolio
							var alphaBuyMax = u[i] - xx[i];
							var alphaBuyMin = Math.min(0, alphaBuyMax);
							if (cashPosition >= alphaBuyMin) {
								var alphaMin = alphaBuyMin;
								var alphaMax = Math.min(alphaBuyMax, cashPosition);
								assetsSwitches.push([-1, i, alphaMin, alphaMax]);
							}
						}

						if ((xx[i] < l[i] || xx[i] == 0) && xx[i] < u[i]) { // The asset i is NOT existing in the portfolio and it might be possible to buy it.
							if (nbNonNullAssets < maxNbAssets) {
								var alphaBuyMin = l[i];
								var alphaBuyMax = Math.max(alphaBuyMin, u[i] - xx[i]);
								if (cashPosition >= alphaBuyMin) {
									var alphaMin = alphaBuyMin;
									var alphaMax = Math.min(alphaBuyMax, cashPosition);
									assetsSwitches.push([-1, i, alphaMin, alphaMax]);
								}
							}
						}
					}
				}
				if (assetsSwitches.length == 0) {
					throw new Error('internal error: no feasible switch');
				}
				var assetSwitch = assetsSwitches[randomIndex(0, assetsSwitches.length)];
				var sellAssetIdx = assetSwitch[0];
				var buyAssetIdx = assetSwitch[1];
				var quantity;
				if (assetSwitch.length == 3) {
					quantity = assetSwitch[2];
					if (sellAssetIdx != -1) {
						xx[sellAssetIdx] = 0; // full sell case
					}
					if (buyAssetIdx != -1) {
						xx[buyAssetIdx] += quantity;
					}
				}
				else if (assetSwitch.length == 4) {
					quantity = random(assetSwitch[2], assetSwitch[3]);
					if (sellAssetIdx != -1) {
						xx[sellAssetIdx] -= quantity;
					}
					if (buyAssetIdx != -1) {
						xx[buyAssetIdx] += quantity;
					}
				}
				else {
					throw new Error('internal error: unexpected lenght of the switches structure');
				}
				return xx;
			}
			var x0;
			try {
				var minExposure = 1;
				if (!fullInvestmentConstraint) {
					minExposure = Math.random();
				}
				
				x0 = randomWeights(nbAssets, { maxIter: maxIterationsInitPoint,
			                                        constraints: {
											            minExposure: minExposure, maxExposure: 1,
														minNbAssets: minNbAssets, maxNbAssets: maxNbAssets, 
											            minWeights: lowerBounds, maxWeights: upperBounds
											        }
												  });
			}
			catch (e) {
				if (e.message === "maximum number of iterations reached") {
					throw new Error("infeasible problem detected");
				}
				else {
					throw(e);
				}
			}
			weights = thresholdAcceptingSolve_(f, x0,
											   {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
											    neighbourGenerator: neighbourGenerator, 
												neighbourGeneratorParameters: { lowerBounds: lowerBounds, 
		                                                                        upperBounds: upperBounds }})[0];
			weights = new Matrix_(weights);
		}
		else if (optimizationMethod == 'exact') {
			var minTEValue = Infinity;
			var minTEAssetsIndexes = [];
			var minTEAssetsWeights = [];
			for (var K = minNbAssets; K <= maxNbAssets; ++K) {
				var nextKSubsetIterator = new kSubsetsIterator_(nbAssets, K, false);
				var nextKSubset = nextKSubsetIterator.next();
				
				while (nextKSubset != -1) {
					var subsetNbAssets = nextKSubset.length;
					var subsetAssetsIdx = typeof UInt32Array === 'function' ? new UInt32Array(subsetNbAssets) : new Array(subsetNbAssets);
					for (var i = 0; i < nextKSubset.length; ++i) {
						subsetAssetsIdx[i] = nextKSubset[i];
					}
					var subsetAssetsReturns = Matrix_.fill(nbPeriods, subsetNbAssets, 
															function(i,j) { 
																return assetsReturns.getValueAt(i, subsetAssetsIdx[j-1]); 
															});
					var subsetLowerBounds = typeof Float64Array === 'function' ? new Float64Array(subsetNbAssets) : new Array(subsetNbAssets);
					var subsetUpperBounds = typeof Float64Array === 'function' ? new Float64Array(subsetNbAssets) : new Array(subsetNbAssets);
					for (var i = 0; i < subsetNbAssets; ++i) {
						subsetLowerBounds[i] = lowerBounds[subsetAssetsIdx[i]-1];
						subsetUpperBounds[i] = upperBounds[subsetAssetsIdx[i]-1];
					}
					var subsetAssetsWeights;
					var subsetPortfolioTrackingError = Infinity;
					try {
						var subsetSol = computeMinimumTrackingErrorVolatilityPortfolio(subsetAssetsReturns, benchmarkReturns, 
																					   subsetLowerBounds, subsetUpperBounds,
																					   fullInvestmentConstraint);
																					   
						subsetAssetsWeights = subsetSol[0];
						subsetPortfolioTrackingError = subsetSol[1];
					}
					catch (e) {
						if (e.message !== "infeasible problem detected: the restricted simplex is empty") {
							throw(e);
						}
					}
					if (subsetPortfolioTrackingError < minTEValue) {
						minTEValue = subsetPortfolioTrackingError;
						minTEAssetsIndexes = subsetAssetsIdx;
						minTEAssetsWeights = subsetAssetsWeights;
					}
					var nextKSubset = nextKSubsetIterator.next();
				}
			}
			if (minTEValue != Infinity) {
				weights = Matrix_.zeros(nbAssets, 1);
				for (var i = 0; i < minTEAssetsIndexes.length; ++i) {
					weights.data[minTEAssetsIndexes[i] - 1] = minTEAssetsWeights.data[i];
				}
			}
			else {
				throw new Error('infeasible problem detected');
			}
		}
		else {
			throw new Error('internal error: unsupported optimisation method');
		}
	}
	return weights.toArray();
}
function mostDiversifiedWeights (sigma, opt) {
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	var epsVolatility = opt.epsVolatility;
	if (epsVolatility === undefined) {
		epsVolatility = 1e-4;
	}
	var vol = new Matrix_(sigma).toCovarianceMatrix().getStandardDeviations(); // assets volatilities
	var efficientFrontier = new MeanVarianceEfficientFrontierGsmo(vol, sigma, { optimizationMethodParams: { epsGsmo: opt.eps, 
									                                                                        maxIterGsmo: opt.maxIter },
																				constraints: { fullInvestment: opt.constraints.fullInvestment, 
																				               minWeights: opt.constraints.minWeights, 
																							   maxWeights: opt.constraints.maxWeights}
																			  });
	efficientFrontier.restrict("minVolatility", epsVolatility);
	var rf = 0;
	var portfolio = efficientFrontier.computeMaximumSharpeRatioEfficientPortfolio(rf);
	return portfolio[0].toArray();
}
function numericalOptimizationWeights (nbAssets, fct, opt) {
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.optimizationMethodParams === undefined) {
		opt.optimizationMethodParams = {};
	}
	var optimizationMethod = opt.optimizationMethod;
	if (optimizationMethod === undefined) {
		optimizationMethod = 'grid-search';
	}
	if (optimizationMethod === 'grid-search') {
		if (opt.optimizationMethodParams.k === undefined) {
			opt.optimizationMethodParams.k = nbAssets;
		}
	}
	if (optimizationMethod === 'grid-search') {
		return simplexGridSearch_(fct, nbAssets, opt.optimizationMethodParams.k, opt.constraints.minWeights, opt.constraints.maxWeights);
	}
	else {
	    throw new Error('unsupported optimisation method');
	}
}
function postOptimizationWeights (originalWeights, opt) {
	if (opt === undefined) {
		opt = { roundingMethodParams: {}, 
		        roundingOptimizationMethodParams: {} };
	}
	if (opt.roundingOptimizationMethodParams === undefined) {
		opt.roundingOptimizationMethodParams = {};
	}
	if (opt.roundingMethodParams === undefined) {
		opt.roundingMethodParams = {};
	}
	if (opt.roundingMethod === undefined) {
		opt.roundingMethod = "rationalizing";
	}
	if (opt.roundingMethod == "rationalizing" && opt.roundingMethodParams.k == undefined) {
		opt.roundingMethodParams.k = 20;
	}
	var roundingMethod = opt.roundingMethod;
	if (roundingMethod != 'rationalizing' && roundingMethod != 'roundlotting') {
		throw new Error('unsupported rounding method');
	}
	var k;
	if (roundingMethod == "rationalizing") {
		k = opt.roundingMethodParams.k;
	}
	var portfolioValue = opt.roundingMethodParams.portfolioValue;
	var assetsPrices = opt.roundingMethodParams.assetsPrices;
	var sizeLots = opt.roundingMethodParams.sizeLots;
	if (roundingMethod == 'roundlotting') {
		if (portfolioValue == undefined || portfolioValue == 0) {
			throw new Error('missing portfolio value');
		}
		if (assetsPrices == undefined) {
			throw new Error('missing assets prices');
		}
		if (sizeLots === undefined) {
			sizeLots = typeof Float64Array === 'function' ? new Float64Array(originalWeights.length) : new Array(originalWeights.length);		
			for (var i = 0; i < sizeLots.length; ++i) {
				sizeLots[i] = 1;
			}
		}
	}
	
	var nRounds;
	var nSteps;
	var nDeltas;
	if (roundingMethod == 'roundlotting') {
		nRounds = opt.roundingOptimizationMethodParams.nRounds;
		if (nRounds === undefined) {
			nRounds = 3;
		}
		
		nSteps = opt.roundingOptimizationMethodParams.nSteps;
		if (nSteps === undefined) {
			nSteps = 5000;
		}

		nDeltas = opt.roundingOptimizationMethodParams.nDeltas;
		if (nDeltas === undefined) {
			nDeltas = nSteps;
		}
	}
	if (roundingMethod == "rationalizing") {
		var roundedWeights = simplexRationalRounding_(originalWeights, k);
		return roundedWeights;
	}
	else if  (roundingMethod == "roundlotting") {
		var eps = 1e-12;
		var nbAssets = originalWeights.length;
		var originalWeights = new Matrix_.fill(nbAssets + 1, 1, function(i,j) { if (i <= nbAssets) { return originalWeights[i-1]; } else { return 0;} });
		function f(q) {
			var w = new Matrix_.fill(nbAssets + 1, 1, 
			                         function(i,j) { 
									    if (i <= nbAssets) { // true asset
											return q[i-1] * sizeLots[i-1] * assetsPrices[i-1] / portfolioValue;
										}
                                        else { // cash
											return q[nbAssets] / portfolioValue;
									    }
									  });
			var rmse = Matrix_.xmy(originalWeights, w).vectorNorm('two');
			return rmse;
		}
		function neighbourGenerator(q, neighbourGeneratorParameters) {					
			function randomIndex(min, max) { // 0, n --> 0, n-1
				return Math.floor(Math.random() * (max - min)) + min; //The maximum is exclusive and the minimum is inclusive
			}
			var nbSellableAssets = 0;
			for (var i = 0; i < nbAssets; ++i) {
				if (q[i] > 0) {
					++nbSellableAssets;
				}
			}
			var sellableAssets = typeof UInt32Array === 'function' ? new UInt32Array(nbSellableAssets) : new Array(nbSellableAssets);
			for (var i = 0, j = 0; i < nbAssets; ++i) {
				if (q[i] > 0) {
					sellableAssets[j] = i;
					++j;
				}
			}
			if (sellableAssets.length != 0) {
				var toSellIndex = sellableAssets[randomIndex(0, nbSellableAssets)];
				var toSellNbLots = randomIndex(1, q[toSellIndex] + 1);
				q[toSellIndex] -= toSellNbLots;
				var cashAddition = assetsPrices[toSellIndex] * sizeLots[toSellIndex] * toSellNbLots;
				q[nbAssets] += cashAddition;
			}
			var nbBuyableAssets = 1; // cash
			for (var i = 0; i < nbAssets; ++i) {
				if (originalWeights.getValue(i + 1, 1) > eps) {
					++nbBuyableAssets;
				}
			}
			var buyableAssets = typeof UInt32Array === 'function' ? new UInt32Array(nbBuyableAssets) : new Array(nbBuyableAssets);
			for (var i = 0, j = 0; i < nbAssets; ++i) {
				if (originalWeights.getValue(i + 1, 1) > eps) {
					buyableAssets[j] = i;
					++j;
				}
			}
			buyableAssets[nbBuyableAssets - 1] = nbAssets;
			var toBuyIndex = buyableAssets[randomIndex(0, nbBuyableAssets)];
			if (toBuyIndex != nbAssets) {
				var toBuyMaxNbLots = Math.floor( q[nbAssets] / (assetsPrices[toBuyIndex] * sizeLots[toBuyIndex]) );
				var toBuyNbLots = randomIndex(0, toBuyMaxNbLots + 1);
				q[toBuyIndex] += toBuyNbLots;
				var cashRemoval = assetsPrices[toBuyIndex] * sizeLots[toBuyIndex] * toBuyNbLots;
				q[nbAssets] -= cashRemoval;
			}	
			return q;
		}
		var q0 = typeof Float64Array === 'function' ? new Float64Array(nbAssets + 1) : new Array(nbAssets + 1);
		var val = 0;
		for (var i = 0; i < nbAssets; ++i) {
			var lotPrice = assetsPrices[i] * sizeLots[i];
			q0[i] = Math.floor( originalWeights.getValue(i + 1, 1) * portfolioValue / lotPrice );

			val += q0[i] * lotPrice;
		}
		var cash = portfolioValue - val;
		if (cash < 0) {
			throw new Error('internal error: negative cash during intitial feasible point generation');
		}
		q0[nbAssets] = cash;
		while (true) {
			var assetIdx = -1;
			var fMin = f(q0);
			for (var i = 0; i < nbAssets; ++i) {
				var lotPrice = assetsPrices[i] * sizeLots[i];
				if (lotPrice < cash) {
					q0[i] += 1;
					q0[nbAssets] = cash - lotPrice;
					
					var fNew = f(q0);
					if (fNew <= fMin) {
						assetIdx = i;
						fMin = fNew;
					}
					q0[i] -= 1;
					q0[nbAssets] = cash;
				}
			}
			if (assetIdx == -1) {
				break;
			}
			q0[assetIdx] += 1;
			cash -= assetsPrices[assetIdx] * sizeLots[assetIdx];
		}
		if (cash < 0) {
			throw new Error('internal error: negative cash during intitial feasible point generation');
		}
		q0[nbAssets] = cash;
		var q = thresholdAcceptingSolve_(f, q0,
					  					 {nSteps: nSteps, nRounds: nRounds, nDeltas: nDeltas,
										  neighbourGenerator: neighbourGenerator})[0];
		var roundedWeights = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1] * sizeLots[i-1] * assetsPrices[i-1] / portfolioValue; }).toArray();
		var qLots = new Matrix_.fill(nbAssets, 1, function(i,j) { return q[i-1]; }).toArray();
		var cash = q[nbAssets];
		return [qLots, roundedWeights, cash];
	}
	else {
		throw new Error('internal error: unsupported rounding method');
	}
}
function proportionalMinimumVarianceWeights (sigma, opt) {
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbRows;
	var rowsElements = sigma.toRowArray();
	var rowsAverages = new Array(nbAssets);
	for (var i = 0; i < nbAssets; ++i) {
		rowsAverages[i] = mean_(rowsElements[i]);
	}
	var elementsMean = mean_(rowsAverages);
	var elementsStddev = sampleStddev_(rowsAverages);
	var weights = new Matrix_(rowsAverages, 1).elemMap(function(i, j, val) { 
		return 1 - normcdf_((val - elementsMean)/elementsStddev);
	});
	weights = weights.normalize(weights);
	var invVariancesWeights = sigma.diagonal().elemMap(function(i,j,val) { return 1/val; }).normalize();
	weights = Matrix_.vectorHadamardProduct(weights, invVariancesWeights).normalize();
	return weights.toArray();
}
function randomSubspaceOptimizationWeights(nbAssets, subsetOptimizationFct, opt) {			
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetOptimizationFctOpt === undefined) {
		opt.subsetOptimizationFctOpt = {};
	}
	if (opt.subsetOptimizationFctOpt.constraints === undefined) {
		opt.subsetOptimizationFctOpt.constraints = {};
	}
	if (nbAssets <= 1) {
		return [1];
	}
	var sizeSubsets = opt.sizeSubsets;
	if (sizeSubsets === undefined) {
		sizeSubsets = Math.max(2, Math.floor(Math.sqrt(nbAssets)));
	}
	if (sizeSubsets <= 1 || sizeSubsets >= nbAssets + 1) {
		throw new Error('the size of the subsets of assets must be between 2 and ' + nbAssets.toString());
	}
	var subsetsGenerationMethod = opt.subsetsGenerationMethod;
	if (subsetsGenerationMethod === undefined) {
		subsetsGenerationMethod = 'random';
	}
	var nbRandomSubsets = opt.nbRandomSubsets;
	if (nbRandomSubsets === undefined) {
		nbRandomSubsets = 128;		
	}
	if (sizeSubsets === nbAssets) {
		nbRandomSubsets = 1;
	}
	var nbSubsets;
	if (subsetsGenerationMethod === 'random') {
		nbSubsets = nbRandomSubsets;
	}
	else if (subsetsGenerationMethod === 'deterministic') {
		nbSubsets = binomial_(nbAssets, sizeSubsets);
	}
	else {
		throw new Error('unsupported subsets of assets generation method');
	}
	var subsetsAggregationMethod = opt.subsetsAggregationMethod;
	if (subsetsAggregationMethod === undefined) {
		subsetsAggregationMethod =  'average';
	}
	if (subsetsAggregationMethod !== 'average' && 
	    subsetsAggregationMethod !== 'median') {
		throw new Error('unsupported aggregation method');
	}
	var maxInfeasibleSubsetsRatio = opt.maxInfeasibleSubsetsRatio;
	if (maxInfeasibleSubsetsRatio === undefined) {
		maxInfeasibleSubsetsRatio =  0;
	}
	var subsetAssetsIdxIterator; 
	if (subsetsGenerationMethod === 'random') {
		subsetAssetsIdxIterator = new randomKSubsetsIterator_(nbAssets, sizeSubsets, false); // use no array copy in the subsets generation to improve performances
	}
	else if (subsetsGenerationMethod === 'deterministic') {
		subsetAssetsIdxIterator = new kSubsetsIterator_(nbAssets, sizeSubsets, false); // use no array copy in the subsets generation to improve performances
	}
	else {
		throw new Error('unsupported subsets generation method');
	}
	var subsetOptimizationFctOpt = opt.subsetOptimizationFctOpt;
	var nbFeasibleGeneratedPortfolios = 0;
	var generatedPortfoliosWeights = new Array(nbSubsets); 
	if (opt.constraints.minWeights) {
		subsetOptimizationFctOpt.constraints.minWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubsets) : new Array(sizeSubsets);
	}
	if (opt.constraints.maxWeights) {
		subsetOptimizationFctOpt.constraints.maxWeights = typeof Float64Array === 'function' ? new Float64Array(sizeSubsets) : new Array(sizeSubsets);
	}
	for (var k = 0; k < nbSubsets; ++k) {
		var subsetAssetsIdx = subsetAssetsIdxIterator.next();
		if (opt.constraints.minWeights) {
			for (var i = 0; i < sizeSubsets; ++i) {
				subsetOptimizationFctOpt.constraints.minWeights[i] = opt.constraints.minWeights[subsetAssetsIdx[i]-1];
			}
		}
		if (opt.constraints.maxWeights) {
			for (var i = 0; i < sizeSubsets; ++i) {
				subsetOptimizationFctOpt.constraints.maxWeights[i] = opt.constraints.maxWeights[subsetAssetsIdx[i]-1]
			}
		}
		try {
			var subsetWeights = subsetOptimizationFct(subsetAssetsIdx, subsetOptimizationFctOpt);
			
			if (subsetWeights.length != sizeSubsets) {
				throw new Error('internal error: the portfolio optimization method did not return the expected number of weights');
			}
		}
		catch (e) {
			if (e.message === "infeasible portfolio optimization problem") {
				continue;
			}
			else {
				throw(e);
			}
		}
		var weights = Matrix_.zeros(nbAssets, 1);
		for (var i = 0; i < sizeSubsets; ++i) {
			weights.setValueAt(subsetAssetsIdx[i], 1, 
							   subsetWeights[i]);
		}
		generatedPortfoliosWeights[nbFeasibleGeneratedPortfolios++] = weights;
	}
	generatedPortfoliosWeights.length = nbFeasibleGeneratedPortfolios;
	if (nbFeasibleGeneratedPortfolios === 0) {
		throw new Error('no feasible portfolio generated');
	}
	else {
		var nbInfeasibleGeneratedPortfolios = nbSubsets - nbFeasibleGeneratedPortfolios;
		
		if (nbInfeasibleGeneratedPortfolios >= 1 && nbInfeasibleGeneratedPortfolios/nbSubsets >= maxInfeasibleSubsetsRatio) {
			throw new Error('too many infeasible portfolios generated');
		}
	}
	var weights = null;
	if (subsetsAggregationMethod == 'average') {
		weights = geometricCenter_(generatedPortfoliosWeights);
	}
	else if (subsetsAggregationMethod == 'median') {
		weights = geometricMedian_(generatedPortfoliosWeights);
	}
	else  {
		throw new Error('internal error');
	}
	return weights.toArray();
}
function randomSubspaceMeanVarianceOptimizationWeights(mu, sigma, opt) {	
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetsOpt === undefined) {
		opt.subsetsOpt = {};
	}
	if (opt.subsetsOpt.constraints === undefined) {
		opt.subsetsOpt.constraints = {};
	}
	var mu = new Matrix_(mu);
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbColumns;
	if (opt.sizeSubsets === undefined) {
		var a = 1;
		var b = 3;
		var c = -2*Math.sqrt(nbAssets*(nbAssets+3)/2);
		var b_p = b/2; // reduced discriminant
		var sign_b_p = 1;
		var disc = b_p*b_p - a*c; // > 0 because a,b are positive and c is negative
		var q = -(b_p + Math.sqrt(disc));
		var r2 = c/q; // positive because c and q are negative
		
		opt.sizeSubsets = Math.max(2, Math.floor(r2));
	}
	opt.subsetOptimizationFctOpt = opt.subsetsOpt;
	function subsetMeanVarianceOptimization(subsetAssetsIdx, subsetOptimizationFctOpt) {
		var subsetMu = mu.submatrix(subsetAssetsIdx, [1]);
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);
		try {
			return meanVarianceOptimizationWeights(subsetMu, subsetSigma, subsetOptimizationFctOpt);
		}
		catch (e) {
			if (e.message.includes('no matching efficient portfolio') ||
				e.message === 'infeasible problem detected: the restricted simplex is empty') {
				throw new Error('infeasible portfolio optimization problem');
			}
			else {
				throw(e);
			}
		}
		
	}
	return randomSubspaceOptimizationWeights(nbAssets, subsetMeanVarianceOptimization, opt);
}
function randomSubspaceGlobalMinimumVarianceOptimizationWeights(sigma, opt) {	
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.subsetsOpt === undefined) {
		opt.subsetsOpt = {};
	}
	if (opt.subsetsOpt.optimizationMethodParams === undefined) {
		opt.subsetsOpt.optimizationMethodParams = {};
	}
	if (opt.subsetsOpt.constraints === undefined) {
		opt.subsetsOpt.constraints = {};
	}
	var mu = opt.mu; // the optional assets returns
	if (mu != undefined) {
		mu = new Matrix_(mu);
	}
	var sigma = new Matrix_(sigma);
	var nbAssets = sigma.nbColumns;
	if (opt.sizeSubsets === undefined) {
		var a = 1;
		var b = mu == undefined ? 1 : 3;
		var c = mu == undefined ? -2*Math.sqrt(nbAssets*(nbAssets+1)/2) : -2*Math.sqrt(nbAssets*(nbAssets+3)/2);
		var b_p = b/2; // reduced discriminant
		var sign_b_p = 1;
		var disc = b_p*b_p - a*c; // > 0 because a,b are positive and c is negative
		var q = -(b_p + Math.sqrt(disc));
		var r2 = c/q; // positive because c and q are negative
		
		opt.sizeSubsets = Math.max(2, Math.floor(r2));
	}
	opt.subsetOptimizationFctOpt = opt.subsetsOpt;
	function subsetGlobalMinimumVarianceOptimization(subsetAssetsIdx, subsetOptimizationFctOpt) {
		var subsetSigma = sigma.submatrix(subsetAssetsIdx, subsetAssetsIdx);
		var subsetMu;
		if (mu != undefined) {
			subsetMu = mu.submatrix(subsetAssetsIdx, [1]);
			subsetOptimizationFctOpt.mu = subsetMu;
		}
		try {
			return globalMinimumVarianceWeights(subsetSigma, subsetOptimizationFctOpt);
		}
		catch (e) {
			if (e.message === 'infeasible problem detected: the restricted simplex is empty') {
				throw new Error('infeasible portfolio optimization problem');
			}
			else {
				throw(e);
			}
		}
		
	}
	return randomSubspaceOptimizationWeights(nbAssets, subsetGlobalMinimumVarianceOptimization, opt);
}
function randomWeights (nbAssets, opt) {
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	var cardinalityConstraints = opt.constraints.minNbAssets || opt.constraints.maxNbAssets;
	var nbMinAssets = opt.constraints.minNbAssets;
	if (nbMinAssets === undefined) {
		nbMinAssets = 1;
	}
	var nbMaxAssets = opt.constraints.maxNbAssets;
	if (nbMaxAssets === undefined) {
		nbMaxAssets = nbAssets;
	}
	var minExposure = opt.constraints.minExposure;
	if (minExposure === undefined) {
		minExposure = 1;
	}
	var maxExposure = opt.constraints.maxExposure;
	if (maxExposure === undefined) {
		maxExposure = 1;
	}
	var maxIterations = opt.maxIter;
	if (maxIterations === undefined) {
		maxIterations = 10000;
	}
	var nbIter = -1;
	while (true) {
		++nbIter;
		if (maxIterations !== -1 && nbIter > maxIterations) {
			throw new Error('maximum number of iterations reached');
		}
		var nbSelectedAssets = Math.floor(Math.random() * (nbMaxAssets - nbMinAssets +1)) + nbMinAssets;
		var selectedAssetsIdx = new randomKSubsetsIterator_(nbAssets, nbSelectedAssets, false).next();
		var portfolioExposure = Math.random() * (maxExposure - minExposure) + minExposure;
		var nbSlackAssets = 0;
		if (portfolioExposure !== 1) {
			nbSlackAssets = 1;
			var lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbSelectedAssets + nbSlackAssets) : new Array(nbSelectedAssets + nbSlackAssets);
			var upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbSelectedAssets + nbSlackAssets) : new Array(nbSelectedAssets + nbSlackAssets);
			for (var i = 0; i < nbSelectedAssets; ++i) {
				lowerBounds[i] = 0;
				upperBounds[i] = 1;
			}
			var portfolioExposureWeightConstraint = 1 - portfolioExposure;	
			lowerBounds[nbSelectedAssets] = portfolioExposureWeightConstraint;
			upperBounds[nbSelectedAssets] = portfolioExposureWeightConstraint;
		}
		if (opt.constraints.minWeights) {
			if (portfolioExposure === 1) {
				var lowerBounds = typeof Float64Array === 'function' ? new Float64Array(nbSelectedAssets) : new Array(nbSelectedAssets);
			}
			for (var i = 0; i < nbSelectedAssets; ++i) {
				lowerBounds[i] = opt.constraints.minWeights[selectedAssetsIdx[i]-1];
			}
		}
		if (opt.constraints.maxWeights) {
			if (portfolioExposure === 1) {
				var upperBounds = typeof Float64Array === 'function' ? new Float64Array(nbSelectedAssets) : new Array(nbSelectedAssets);
			}
			for (var i = 0; i < nbSelectedAssets; ++i) {
				upperBounds[i] = opt.constraints.maxWeights[selectedAssetsIdx[i]-1];
			}
		}		
		try {
			var sumBounds = simplexEmptinessCheck_(nbSelectedAssets + nbSlackAssets, lowerBounds, upperBounds);
		}
		catch (e) {
			continue;
		}
		var selectedAssetsWeights = new simplexRandomSampler_(nbSelectedAssets + nbSlackAssets, lowerBounds, upperBounds).sample();
		try {
			if (cardinalityConstraints) {
				for (var i = 0; i < nbSelectedAssets; ++i) {
					if (selectedAssetsWeights[i] == 0) { 
						throw new Error('generated weights not compatible with cardinality constraints');
					}
				}
			}
		}
		catch (e) {
			continue;
		}
		break;
	}
	var weights = Matrix_.zeros(nbAssets, 1);
	for (var i = 0; i < nbSelectedAssets; ++i) {
		var assetIdx = selectedAssetsIdx[i];
		var assetWeight = selectedAssetsWeights[i];
		weights.setValueAt(assetIdx, 1, assetWeight);
	}
	return weights.toArray();
}
function riskBudgetingWeights (sigma, rb, opt) {
	function computeVolatility(x, sigma_x) {
		var sigma_x_x = Matrix_.vectorDotProduct(sigma_x, x);		
		if (Math.abs(sigma_x_x) <= epsSdp) {
			sigma_x_x = 0;
		}
		else if (sigma_x_x < 0 && sigma_x_x < -epsSdp) {
			throw new Error('negative volatility during iteration ' + iter + ', covariance matrix might not be semi-definite positive');
		}
		var s_x = Math.sqrt(sigma_x_x);
		return s_x;
	}
	function coordinatesDescentSolve(lambda, lowerBounds, upperBounds) {
		var x = Matrix_.fill(nbAssets, 1, function(i,j) { return 1/nbAssets; });
		var log_x = x.elemMap(function(i,j,val) { return Math.log(val);});
		var x_old = Matrix_.copy(x);
		var sigma_x = Matrix_.xy(sigma, x); // SIGMA*x
		var s_x = computeVolatility(x, sigma_x); // sigma(x)
		
		var rb_d_log_x = Matrix_.vectorDotProduct(rb, log_x); // Matrix_.vectorDotProduct(rb, log_x)
		var obj_old = 0.5 * s_x - rb_d_log_x; // the objective function to be minimized, c.f. formula 3 of the second reference
		var cs = new coordinatesSampler(nbAssets);
		var cycle = 0;
		while (true) {
			if	(nbCycles != -1 && cycle == nbCycles) {
				break;
			}
			if (nbCycles == -1 && cycle >= maxCycles) {
				throw new Error('maximum number of cycles reached: ' + maxCycles);
			}
			++cycle;
			cs.generateCoordinates();
			var obj_old_i
			var obj_new_i;
			while (true) {
				var i = cs.sampleCoordinate();
				if (i == -1) {
					break;
				}
				if (opt.coordinatesSampler === 'acf') {
					obj_old_i = 0.5 * s_x - rb_d_log_x;
				}
				var xi_old = x.data[i-1];
				var log_xi_old = log_x.data[i-1];
				var a = sigma.data[(i-1)*sigma.nbColumns + (i-1)]; // sigma_i^2, always > 0
				var b = sigma_x.data[i-1] - x.data[i-1] * a; // (SIGMA*x)_i - x_i*sigma_i^2, might be any sign
				var c = -lambda * rb.data[i-1] * s_x; // -lambda * b_i * sigma(x), always <= 0 (== 0 iff the covariance matrix is semi-definite positive and sigma(x) == 0)
				var b_p = b/2; // reduced discriminant
				var sign_b_p = (b_p >= 0) ? 1 : -1; // Math.sign is not supported everywhere plus it is mandatory that for b_p == 0 this returns 1
				var disc = b_p*b_p - a*c;
				if (disc < 0) {
				   throw new Error('internal error: negative discriminant detected, the covariance matrix might not be semi-definite positive');
				}
				var q = -(b_p + sign_b_p * Math.sqrt(disc));
				var r1 = q/a;
				var r2 = c/q;
				
				var xi_star;
				if (r1 > 0) {
					xi_star = r1;
				}
				else if (r2 > 0) {
					xi_star = r2;
				}
				else {
					throw new Error('internal error: no strictly positive root detected, the covariance matrix might not be semi-definite positive');
				}
				if (lowerBounds && xi_star < lowerBounds[i-1]) {
					xi_star = lowerBounds[i-1];
				}
				if (upperBounds && xi_star > upperBounds[i-1]) {
					xi_star = upperBounds[i-1];
				}
				x.data[i-1] = xi_star;
				log_x.data[i-1] = Math.log(xi_star);
				
				rb_d_log_x -= log_xi_old * rb.data[i-1];
				rb_d_log_x += log_x.data[i-1] * rb.data[i-1];
				for (var j = 1; j <= nbAssets; ++j) {
					sigma_x.data[j-1] += sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_star;
					sigma_x.data[j-1] -= sigma.data[(j-1)*sigma.nbColumns + (i-1)] * xi_old;
				}
				s_x = computeVolatility(x, sigma_x);
				if (opt.coordinatesSampler === 'acf') {
					obj_new_i = 0.5 * s_x - rb_d_log_x;
					var gain = obj_old_i - obj_new_i;
					if (cycle == 1) {
						cs.updateAverageGain(gain);
					}
					else {
						cs.updateSchedulingPreference(i, gain);
					}
				}	
			}
			var obj_new = 0.5 * s_x - rb_d_log_x;
			if (nbCycles == -1 && Math.abs(obj_new - obj_old) <= eps) {
				break;
			}
			obj_old = obj_new;
		}
		return x;
	}
	function cyclicCoordinatesSampler(n) {
		this.n = n;
		this.k = n;
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				return ++this.k;
			}        
		}
		this.generateCoordinates = function() {
			this.k = 0;
		}
	}
	function shuffledCyclicCoordinatesSampler(n) {
		this.n = n;
		this.k = n;
		this.r = typeof Uint32Array === 'function' ? new Uint32Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.r[i] = i+1;
		}
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				return this.r[this.k++];
			}        
		}
		this.generateCoordinates = function() {
			this.r = new randomPermutationsIterator_(undefined, this.r, true).next();
			this.k = 0;
		}
	}
	function randomizedCoordinatesSampler(n) {
		this.n = n;
		this.k = n;
		
		var p = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			p[i] = 1/n;
		}
		this.s = new aliasMethodSampler_(p);
		this.sampleCoordinate = function() {
			if (this.k == this.n) {
				return -1;
			}
			else {
				++this.k;
				return this.s.sample() + 1;
			}        
		}
		this.generateCoordinates = function() {
			this.k = 0;
		}
	}
	function adaptiveCoordinateFrequenciesSampler(n) {
		this.n = n;
		this.k = 0;
		this.change_rate = 1/5; // c in the fifth reference
		this.p_min = 1/20;
		this.p_max = 20;	
		this.idx_set_tmp = typeof Uint32Array === 'function' ? new Uint32Array(2*n) : new Array(2*n);
		this.idx_set = null;
		this.idx_set_len = 0;
		this.pref = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.pref[i] = 1;
		}
		this.prefsum = n; // sum p_i, i=1..n
		this.acc = typeof Float64Array === 'function' ? new Float64Array(n) : new Array(n);
		for (var i = 0; i < n; ++i) {
			this.acc[i] = 0.5;
		}
		this.gain_learning_rate = 1/n;
		this.average_gain = 0;
		this.sampleCoordinate = function() {
			if (this.k == this.idx_set_len) {
				return -1;
			}
			else {
				return this.idx_set[this.k++];
			}        
		}
		this.generateCoordinates = function() {
			this.k = 0;
			this.idx_set_len = 0;
			var q = this.n / this.prefsum;
			for (var i = 0; i < this.n; ++i) {
				var a_i = this.acc[i] + q * this.pref[i];
				var a_i_f = Math.floor(a_i);
				
				for (var j = 0; j < a_i_f; ++j) {
					this.idx_set_tmp[this.idx_set_len] = i + 1;
					this.idx_set_len++;
				}
				
				this.acc[i] = a_i - a_i_f;
			}
			this.idx_set = new randomPermutationsIterator_(undefined, this.idx_set_tmp.slice(0, this.idx_set_len), true).next();
		}
		this.updateSchedulingPreference = function(i, gain) {
			var p_new = this.pref[i-1] * Math.exp(this.change_rate * (gain / this.average_gain - 1));
			if (p_new < this.p_min) {
				p_new = this.p_min;
			}
			else if (p_new > this.p_max) {
				p_new = this.p_max;
			}
			this.prefsum = this.prefsum + p_new - this.pref[i-1];
			this.pref[i-1] = p_new;
			this.average_gain = (1 - this.gain_learning_rate) * this.average_gain + this.gain_learning_rate * gain;
		}
		this.updateAverageGain = function(gain) {
			this.average_gain += gain / this.n;
		}
	}
	if (opt === undefined) {
		opt = {};
	}
	if (opt.constraints === undefined) {
		opt.constraints = {};
	}
	if (opt.eps === undefined) {
		opt.eps = 1e-10;
	}
	if (opt.epsSdp === undefined) {
		opt.epsSdp = 1e-12;
	}
	if (opt.maxCycles === undefined) {
		opt.maxCycles = 10000;
	}
	if (opt.nbCycles === undefined) {
		opt.nbCycles = -1;
	}
	if (opt.outputPortfolioVolatility === undefined) {
		opt.outputPortfolioVolatility = false;
	}
	if (opt.coordinatesSampler === undefined) {
		opt.coordinatesSampler = 'cyclic';
	}
	var eps = opt.eps;
	var epsSdp = opt.epsSdp;
	var maxCycles = opt.maxCycles;
	var nbCycles = opt.nbCycles;
	
	var outputPortfolioVolatility = opt.outputPortfolioVolatility; 
	
	var coordinatesSampler;
	if (opt.coordinatesSampler === 'cyclic') {
		coordinatesSampler = cyclicCoordinatesSampler;
	}
	else if (opt.coordinatesSampler === 'shuffledCyclic') {
		coordinatesSampler = shuffledCyclicCoordinatesSampler;
	}	
	else if (opt.coordinatesSampler === 'randomized') {
		coordinatesSampler = randomizedCoordinatesSampler;
	}
	else if (opt.coordinatesSampler === 'acf') {
		coordinatesSampler = adaptiveCoordinateFrequenciesSampler;
	}
	else {
		throw new Error('unsupported coordinates sampler');
	}
	var lowerBounds = null;
	var upperBounds = null;
	if (opt.constraints.minWeights) {
		lowerBounds = opt.constraints.minWeights; // no input checks
		
		if (!opt.constraints.maxWeights) {
			upperBounds = typeof Float64Array === 'function' ? new Float64Array(opt.constraints.minWeights.length) : new Array(opt.constraints.minWeights.length);
			for (var i = 0; i < upperBounds.length; ++i) {
				upperBounds[i] = 1;
			}
		}
	}
	if (opt.constraints.maxWeights) {
		upperBounds = opt.constraints.maxWeights; // no input checks
		
		if (!opt.constraints.minWeights) {
			lowerBounds = typeof Float64Array === 'function' ? new Float64Array(opt.constraints.maxWeights.length) : new Array(opt.constraints.maxWeights.length);
			for (var i = 0; i < lowerBounds.length; ++i) {
				lowerBounds[i] = 0;
			}
		}
	}
	var sigma = new Matrix_(sigma);
	var rb = new Matrix_(rb);
	var nbAssets = sigma.nbRows;
	var lambda = 1;
	var weightsURb = coordinatesDescentSolve(lambda);
	weightsURb = weightsURb.normalize();
	var x;
	if (lowerBounds || upperBounds) {
		var sumBounds = simplexEmptinessCheck_(nbAssets, lowerBounds, upperBounds);
		var volURb = computeVolatility(weightsURb, Matrix_.xy(sigma, weightsURb));
		var a = 0.5 * volURb;
		var b = 2 * volURb;
		var lambda_star = bisection_(function (lambda) { 
						 		         var weights = coordinatesDescentSolve(lambda, lowerBounds, upperBounds); 
								         return weights.sum() - 1; 
							         }, 
							         a, b);
		x = coordinatesDescentSolve(lambda_star, lowerBounds, upperBounds);
		
	}
	else {
		x = weightsURb;
	}
	if (outputPortfolioVolatility === true) {
		return [x.toArray(), computeVolatility(x, Matrix_.xy(sigma, x))];
	}
	else {
		return x.toArray();
	}
}

