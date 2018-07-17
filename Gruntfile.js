module.exports = function(grunt) {
  var path = require('path');

  // Load plugins
  grunt.loadNpmTasks('grunt-contrib-qunit');
  grunt.loadNpmTasks('grunt-contrib-uglify');
  grunt.loadNpmTasks('grunt-strip-code');
  grunt.loadNpmTasks('grunt-contrib-concat');
  grunt.loadNpmTasks('grunt-replace');

  grunt.initConfig({
    pkg: grunt.file.readJSON('package.json'),
	
	strip_code: {
      portfolio_allocation_gs: {
		options: {
			blocks: [
			{ start_block: "/* Start Not to be used as is in Google Sheets */",
			  end_block: "/* End Not to be used as is in Google Sheets */"
			}
			]
		},
        files: [
		  { src: ['dist/portfolio_allocation.dist.js'], dest: 'dist/gs/portfolio_allocation.gs'}
		]
      },
      portfolio_allocation_dist: {
		options: {
			blocks: [
			{ start_block: "/* Start Wrapper private methods - Unit tests usage only */",
			  end_block: "/* End Wrapper private methods - Unit tests usage only */"
			}
			]
		},
        files: [
		  { src: ['dist/portfolio_allocation.dev.js'], dest: 'dist/portfolio_allocation.dist.js'}
		]
      },
    },

    replace: {
      portfolio_allocation_gs: {
        options: {
		  preserveOrder: true,
          patterns: [
            {
              match: /self\.(\w+)\s=\sfunction/g,
              replacement: 'function $1'
            },
            {
              match: /self\./g,
              replacement: ''
            }			
          ]
        },
        files: [
          { src: ['dist/gs/portfolio_allocation.gs'], dest: 'dist/gs/portfolio_allocation.gs'}
        ]
      },
	  portfolio_allocation : {
        options: {
		  preserveOrder: true,
          patterns: [
            {
              match: /NONE\((.+?)\);/g,
              replacement: "NONE($1);"
            }		
          ]
        },
        files: [
          { src: ['dist/portfolio_allocation.dev.js'], dest: 'dist/portfolio_allocation.dev.js'}
        ]
	  }
    },

	concat: {
	  portfolio_allocation: {
	    src: ['lib/header.js', 'lib/matrix/matrix.js', 'lib/matrix/covmatrix.js', 'lib/stats/*.js', 'lib/allocation/*.js', 'lib/footer.js'],
	    dest: 'dist/portfolio_allocation.dev.js',
	  }
	},
  
    uglify: {
	  options: {
	    banner: '/*! <%= pkg.name %> - v<%= pkg.version %> (<%= pkg.date %>), <%= pkg.author %> */'
      },
      portfolio_allocation_dev: {
        files: {
          'dist/portfolio_allocation.dev.min.js': ['dist/portfolio_allocation.dev.js']
        }
      },
      portfolio_allocation_dist: {
        files: {
          'dist/portfolio_allocation.dist.min.js': ['dist/portfolio_allocation.dist.js']
        }
      }	  
    },
	
	qunit: { 
      dev: ['test/index_optimisation.html', 'test/index_bitset.html', 'test/index_matrix.html', 'test/index_simplex.html', 'test/index_vector.html', 'test/index_stats.html', 'test/index_combinatorics.html', 'test/index_covmatrix_dev.html', 'test/index_allocation_dev.html', 'test/index_computational_geometry_dev.html'],
	  dist: ['test/index_allocation_dist.html', 'test/index_covmatrix_dist.html']
    }

  });


  //
  grunt.registerTask('test', 'Tests the app.', function() {
	// Build the app in dev mode, and run all dev unit tests, 
	// which are a superset of dist unit tests
	grunt.task.run('concat:portfolio_allocation');
	grunt.task.run('replace:portfolio_allocation')
	grunt.task.run('uglify:portfolio_allocation_dev');
	grunt.task.run('qunit:dev');
  });

  //
  grunt.registerTask('deliver', 'Builds the app into a distributable package.', function() {
    // Build the app in dist mode
	grunt.task.run('concat:portfolio_allocation');
	grunt.task.run('replace:portfolio_allocation')
    grunt.task.run('strip_code:portfolio_allocation_dist');
	grunt.task.run('uglify:portfolio_allocation_dist');

	// Run the non-dev unit tests
	grunt.task.run('qunit:dist');	
  });
  
  //
  grunt.registerTask('deliver-gs', 'Builds the app into a distributable package for Google Sheets.', function() {
    // Starting from dist mode, remove header/footer
	grunt.task.run('strip_code:portfolio_allocation_gs');
	
	// Then change function definitions
	grunt.task.run('replace:portfolio_allocation_gs')
  });
};