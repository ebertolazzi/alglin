Alglin linear algebra package repository
========================================

To use the library you must have installed on your system
a C++11 compiler, rake/ruby and cmake.

This command fetch all the submodules to the last commit

~~~~
ruby setup.rb --latest
~~~~

Use rake to compile all the tests and the library

~~~~
rake
~~~~

Run all the tests using rake

~~~
rake run
~~~

tested on OSX, Ubuntu linux and windows with visual studio 2017.
