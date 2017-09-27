# PortfolioAllocation v0.0.3 ([Changelog](changelog.md))

[![Travis Build Status](https://travis-ci.org/lequant40/portfolio_allocation_js.svg?style=flat)](https://travis-ci.org/lequant40/portfolio_allocation_js)

PortfolioAllocation is a JavaScript library to allocate portfolios of finantical assets, i.e. to compute the proportions of a given set of financial instruments (securities, bonds, exchange traded funds, mutual funds...) to hold in a portfolio so as to optimize specific quantitative criteria related to this portfolio.

PortfolioAllocation is complementary to the JavaScript library [PortfolioAnalytics](https://github.com/lequant40/portfolio_analytics_js), and has been developed in JavaScript for the same reason: I heavily use [Google Sheets](https://www.google.com/sheets/about/) to analyse trading strategies on my blog [Le Quant 40](http://www.lequant40.com/), as well as in my personal trading, and Google Sheets is easily extensible thanks to [Google Apps Script](https://developers.google.com/apps-script/), a JavaScript-based language.


## Features

- Compatible with Google Sheets
- Compatible with any browser supporting ECMAScript 5 (i.e., front-end development)
- Compatible with [Node.js](https://nodejs.org/) (i.e., back-end development)
- (Performances) Automatically uses JavaScript Typed Arrays if available
- Code continuously tested and integrated by [Travis CI](https://travis-ci.org/)
- Code heavily documented using [JSDoc](http://usejsdoc.org/)


## Supported portfolio allocation algorithms

- Equal weights (EW)
  The 1/N portfolio allocation has been popularized by DeMiguel and al. in their research paper *Optimal Versus Naive Diversification: How Inefficient is the 1/N Portfolio Strategy?*.

- Equal risk budgets (ERB)
  Also known as naive risk parity, this portfolio allocation has been analysed by Carvalho and al. in the research paper *Demystifying Equity Risk-Based Strategies: A Simple Alpha Plus Beta Description*.

- Equal risk contributions (ERC) and Risk budgeting (RB)
  Both portfolio allocations have been extensively studied by [Thierry Roncalli](http://www.thierry-roncalli.com/) and al in misc. research papers (*The properties of equally weighted risk contribution portfolios*, *Managing Risk Exposures Using the Risk Budgeting Approach*...), the former being a specific case of the latter.

- Most diversified portfolio (MDP)
  This portfolio allocation is described in the research paper *Toward Maximum Diversification* by Choueifaty and al., and uses the authors' diversification ratio as the criterion for portfolio construction.

- Minimum correlation algorithm (MCA) and Proportional minimum variance algorithm (MVA)
  Both algorithms were discovered by [David Varadi](https://cssanalytics.wordpress.com/), the former being meant as an approximation of the MDP portfolio and the latter as an approximation of the GMV portfolio.

- Cluster risk parity portfolio (CRP)
  This portfolio allocation algorithm was discovered by [David Varadi](https://cssanalytics.wordpress.com/) and [Michael Kapler](http://systematicinvestor.wordpress.com/), and combines the usage of clustering algorithms (for instance, the Fast Threshold Clustering Algorithm - FTCA - of David Varadi) with ERC.

- Global minimum variance portfolio (GMV)
  The global minimum variance portfolio is the portfolio located at the leftmost part of the Markowitz mean-variance efficient frontier, which makes it the portfolio with the smallest attainable volatility among all the mean-variance efficient portfolios.

  
## Usage

### Usage in Google Sheets

If you would like to use PortfolioAllocation in Google Sheets, you can either:

- (Recommended) [Import the external Google Apps Script library](https://developers.google.com/apps-script/guide_libraries) with Script ID 1cTTAt3VZRZKptyhXtjF3EK-jKagdTUl6t29Pc7YidrC5m5ABR6LUy8sC into your spreadsheet script

or:

- Import the JavaScript files from the [dist/gs directory](https://github.com/lequant40/portfolio_allocation_js/tree/master/dist/gs) into your spreadsheet script

In both cases, providing data to the PortfolioAllocation functions is then accomplished your preferred way (c.f. the [PortfolioAnalytics JavaScript library README](https://github.com/lequant40/portfolio_analytics_js) for detailled examples):

- Using a wrapper function in your spreadsheet script, directly accessible from your spreadsheet, to which you can provide a standard data range (A1:A99...)
- Using pure Google Apps Script functions - typically the getRange(...) familly of functions -, optimized for speed
- ...

You can find examples of PortfolioAllocation usage in [this spreadsheet](https://docs.google.com/spreadsheets/d/1ScrwSjr9EgwXfRyPN4IaqVxZvDnqw-hWvVQcJ9Ak590). 


### Usage inside a browser

If you would like to use PortfolioAllocation inside a browser you can download [its source code](http://raw.github.com/lequant40/portfolio_allocation_js/master/dist/portfolio_allocation.dist.js) and/or [its minified source code](http://raw.github.com/lequant40/portfolio_allocation_js/master/dist/portfolio_allocation.dist.min.js).

You then just need to include this code in an HTML page to use it, e.g.:
```html
<script src="portfolio_allocation.dist.min.js" type="text/javascript"></script>
```


### Usage with Node.js

If you would like to use PortfolioAllocation with [Node.js](https://nodejs.org/en/), you simply need to declare it as a dependency of your project 
in your `package.json` file.

Then, this is standard Node.js code, e.g.:

```js
var PortfolioAllocation = require('portfolio-allocation');
...
var w = PortfolioAllocation.riskBudgetingWeights([[0.1,0], [0,0.2]], [0.25, 0.75]);
// w = [0.44948974243459944, 0.5505102575654006]
```


### Examples

#### Risk-based portfolio allocations

```js
PortfolioAllocation.equalWeights(5); 
// EW portfolio

PortfolioAllocation.equalRiskContributionWeights([[0.1, 0], [0, 0.2]]); 
// ERC portfolio

PortfolioAllocation.mostDiversifiedWeights([[0.1, 0], [0, 0.2]], {eps: 1e-10, maxIter: 10000});
// MDP portfolio

PortfolioAllocation.minCorrWeights([[0.1, 0], [0, 0.2]]);
// MCA portfolio

PortfolioAllocation.clusterRiskParityWeights([[0.1,0], [0,0.2]], {clusteringMode: 'ftca'});
// CRP portfolio
```


## How to contribute ?

### Fork the projet from [Github](https://github.com/)...


### Instal the [Grunt](http://gruntjs.com/) dependencies

```
npm install
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

