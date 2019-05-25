# PortfolioAllocation v0.0.7 ([Changelog](changelog.md))

[![Travis Build Status](https://travis-ci.org/lequant40/portfolio_allocation_js.svg?style=flat)](https://travis-ci.org/lequant40/portfolio_allocation_js)

PortfolioAllocation is a JavaScript library designed to help solving the mathematical problem of portfolio allocation and optimization.

In layman's terms, imagine you are faced with the problem of deciding how to invest your available funds into different financial instruments: stocks, bonds, mutual funds, exchange traded funds, cryptocurrencies...

A library such as PortfolioAllocation allows you to easily determine which proportion of which instrument you need to hold so that your available funds are allocated the best possible way.

Do not hesitate to report any bug / request additional features !


## Features

- Compatible with Google Sheets
- Compatible with any browser supporting ECMAScript 5 for front-end development
- Compatible with [Node.js](https://nodejs.org/) for back-end development
- Code continuously tested and integrated by [Travis CI](https://travis-ci.org/)
- Code documented for internal developers using [JSDoc](http://usejsdoc.org/)
- Code documented for users using [GitHub Pages](https://lequant40.github.io/portfolio_allocation_js/)


## Included algorithms

### Portfolio allocation and optimization algorithms

- Equal weights (EW)  
  Analyzed by Victor DeMiguel and al. in their research paper [Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?](https://doi.org/10.1093/rfs/hhm075).

- Equal risk budgets, a.k.a. naive risk parity  
  Described by Raul Leote de Carvalho and al. in the research paper [Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description](https://doi.org/10.3905/jpm.2012.38.3.056).

- Equal risk contributions (ERC) and risk budgeting (RB)  
  Extensively studied by [Thierry Roncalli](http://www.thierry-roncalli.com/) and al in misc. research papers ([The properties of equally weighted risk contribution portfolios](https://doi.org/10.3905/jpm.2010.36.4.060), [Managing Risk Exposures Using the Risk Budgeting Approach](https://ssrn.com/abstract=2009778)...).

- Equal risk bounding (ERB)  
  Described in the research paper [Equal Risk Bounding is better than Risk Parity for portfolio selection](https://doi.org/10.1007/s10898-016-0477-6) by Francesco Cesarone and Fabio Tardella, the ERB portfolio is an ERC portfolio possibly not containing all the assets in the considered universe.

- Cluster risk parity (CRP)  
  Discovered by [David Varadi](https://cssanalytics.wordpress.com/) and [Michael Kapler](http://systematicinvestor.wordpress.com/), the CRP portfolio combines the usage of a clustering algorithm (for instance, the Fast Threshold Clustering Algorithm - FTCA - of David Varadi) with the ERC portfolio.

- Most diversified portfolio (MDP)  
  Introduced in the research paper [Toward Maximum Diversification](https://doi.org/10.3905/JPM.2008.35.1.40) by [Yves Choueifaty](http://www.tobam.fr/yves-choueifaty/) and al., it maximizes what the authors call the diversification ratio, which is the weighted average of the assets volatilities divided by the portfolio total volatility.

- Minimum correlation algorithm (MCA)  
  Discovered by [David Varadi](https://cssanalytics.wordpress.com/), the MCA portfolio is meant to be an approximation of the MDP portfolio.

- Mean-variance optimization (MVO)  
  Based on the modern portfolio theory introduced by Harry M. Markowitz in numerous articles and books ([Portfolio Selection: Efficient Diversification of Investments](https://www.jstor.org/stable/j.ctt1bh4c8h)...), the portfolio obtained through mean-variance optimization is mean-variance efficient, that is, for a given level of return, it possesses the lowest attainable volatility and for a given level of volatility, it possesses the highest attainable return.

- Global minimum variance (GMV)  
  The leftmost portfolio on the mean-variance efficient frontier, the GMV portfolio possesses the smallest attainable volatility among all the mean-variance efficient portfolios.

- Proportional minimum variance algorithm (MVA)  
  Discovered by [David Varadi](https://cssanalytics.wordpress.com/), the MVA portfolio is meant to be an approximation of the GMV portfolio.

- Minimax portfolio  
  Introduced by Martin Young in the research paper [A Minimax Portfolio Selection Rule with Linear Programming Solution](http://www.jstor.org/stable/2634472), the minimax portfolio uses the minimum return as a measure of risk instead of the variance as in the Markowitz framework.

- Random portfolio  
  Random portfolios are generally used to benchmark the performances of portfolio allocation and optimization algorithms, as pioneered by Ronald J. Surz in the article [Portfolio Opportunity Distributions](https://doi.org/10.3905/joi.3.2.36) and latter complemented by Patrick Burns in the article [Random Portfolios for Performance Measurement](https://doi.org/10.1007/3-540-36626-1_11).

- Maximum Sharpe ratio (MSR), a.k.a. (Lintner) tangency portfolio  
  Introduced by John Lintner in the research paper [The Valuation of Risk Assets and the Selection of Risky Investments in Stock Portfolios and Capital Budgets](https://www.jstor.org/stable/1924119), the MSR portfolio possesses the highest Sharpe ratio among all the mean-variance efficient portfolios.

- Random subspace mean-variance optimization portfolio (RSO-MVO)  
  Discovered by [David Varadi](https://cssanalytics.wordpress.com/) and first formally studied by [Benjamin J. Gillen](http://www.its.caltech.edu/~bgillen/index.htm) in the research paper [Subset Optimization for Asset Allocation](http://www.its.caltech.edu/~bgillen/papers/Subsets.pdf), the RSO-MVO portfolio combines the usage of a [random subspace optimization method](https://en.wikipedia.org/wiki/Random_subspace_method) with a mean-variance optimization method.


### Misc. helper algorithms
- Portfolio weights rounding  
  The portfolio weights obtained through a portfolio optimization algorithm are usually provided with several digits after the decimal point. For any practical usage, these theoretical weights need to be rounded off, which can be done thanks to the algorithm described in the research paper [Rounding on the standard simplex: Regular grids for global optimization](https://doi.org/10.1007/s10898-013-0126-2) from Immanuel M. Bomze and al..

- Mean-variance efficient frontier and corner portfolios computation  
  The set of all mean-variance efficient portfolios (the mean-variance efficient frontier), as well its generating discrete set (the set of corner portfolios) can be efficiently computed thanks to a specialized algorithm developped by Harry M. Markowitz: [the critical line method](https://web.stanford.edu/~wfsharpe/mia/opt/mia_opt3.htm).

- Generic random subspace optimization method  
  As a direct extension of the RSO-MVO method, a generic random subspace optimization method has been implemented and is usable with any portfolio optimization method.

- Generic numerical optimization algorithms  
  When no specialized numerical algorithm exist to solve a particular portfolio optimization problem, a slow-but-always-working solution is to use generic numerical optimization algorithms instead (e.g., grid search on the simplex).


## Usage

### Usage in Google Sheets

If you would like to use PortfolioAllocation in Google Sheets, you can either:

- (Recommended) [Import the external Google Apps Script library](https://developers.google.com/apps-script/guide_libraries) with Script ID 1cTTAt3VZRZKptyhXtjF3EK-jKagdTUl6t29Pc7YidrC5m5ABR6LUy8sC into your spreadsheet script

or:

- Import the JavaScript files from the [dist/gs directory](https://github.com/lequant40/portfolio_allocation_js/tree/master/dist/gs) into your spreadsheet script

In both cases, calling the PortfolioAllocation functions is then accomplished your preferred way:

- Using a wrapper function in your spreadsheet script, directly accessible from your spreadsheet, to which you can provide standard data ranges (A1:B3...), e.g.:

```js
function computeERCPortfolioWeights(covarianceMatrix) {
  // Note: The input range coming from the spreadsheet and representing a covariance matrix
  // is directly usable.
    
  // Compute the ERC portfolio weights
  var ercWeights = PortfolioAllocation.equalRiskContributionWeights(covarianceMatrix);
  
  // Return them to the spreadsheet
  return ercWeights;
}
```

- Using pure Google Apps Script functions, for which computations are allowed by Google to last longer, e.g.:

```js
function computeERCPortfolioWeights() {
  // Adapted from https://developers.google.com/apps-script/reference/spreadsheet/sheet#getrangerow-column-numrows
 var ss = SpreadsheetApp.getActiveSpreadsheet();
 var sheet = ss.getSheets()[0];
 var range = sheet.getRange(1, 1, 3, 3); // A1:B3
 var values = range.getValues();

 // Convert the above range into an array
 var covarianceMatrix = [];
 for (var row in values) {
   for (var col in values[row]) {
     covarianceMatrix.push(values[row][col]);
   }
 }
 
  // Compute the ERC portfolio weights
  var ercWeights = PortfolioAllocation.equalRiskContributionWeights(covarianceMatrix);
  
  // Do something with them (use them in a computation, write them back to the spreadsheet, etc.)
  ...
}
```

You can find examples of PortfolioAllocation usage in [this spreadsheet](https://docs.google.com/spreadsheets/d/1ScrwSjr9EgwXfRyPN4IaqVxZvDnqw-hWvVQcJ9Ak590). 


### Usage inside a browser

If you would like to use PortfolioAllocation inside a browser you can download [its source code](http://raw.github.com/lequant40/portfolio_allocation_js/master/dist/portfolio_allocation.dist.js) and/or [its minified source code](http://raw.github.com/lequant40/portfolio_allocation_js/master/dist/portfolio_allocation.dist.min.js).

You then just need to include this code in an HTML page to use it, e.g.:
```html
<script src="portfolio_allocation.dist.min.js" type="text/javascript"></script>
<script type="text/javascript">
  var w = PortfolioAllocation.riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75]);
  // w = [0.44948974243459944, 0.5505102575654006]
</script>
```


### Usage with Node.js

If you would like to use PortfolioAllocation with [Node.js](https://nodejs.org/en/), you simply need to declare it as a dependency of your project in your `package.json` file.

Then, this is standard Node.js code, e.g.:

```js
var PortfolioAllocation = require('portfolio-allocation');
...
var w = PortfolioAllocation.riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75]);
// w = [0.44948974243459944, 0.5505102575654006]
```


## Documentation

A complete documentation, including code examples, can be found [on the GitHub Pages associated to this repository](https://lequant40.github.io/portfolio_allocation_js/).



## How to contribute ?

### Fork the project from [Github](https://github.com/)...


### Instal the [Grunt](http://gruntjs.com/) dependencies and command line

```
npm install
npm install -g grunt-cli
```

### Develop...

### Compile

- The following command generates the files to be used inside a browser or with Node.js in the `dist` directory:

```
grunt deliver
```

- The following command generates the files to be used in Google Sheets in the `dist\gs` directory:

```
grunt deliver-gs
```

### Test

Any of the following two commands run the [QUnit](https://qunitjs.com/) unit tests contained in the `test` directory on the generated file `dist\portfolio_allocation.dev.min.js`:

```
npm test
```

```
grunt test
```

### Submit a pull-request...


## License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)

