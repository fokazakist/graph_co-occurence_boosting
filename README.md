# graph_co-occurence_boosting

This package is graph lpboost with co-occurence patterns and wildcard patterns

1. Prepare
First compile libraries GLPK, then

need to use:
c++ library, boost
-->install on Ubuntu
     sudo apt-get install libboost-dev

2. How to compile
How to make:
$ cd src
$ make
$ cd ../eval_wild
$ make
$ cd ..

3. USAGE

$./lpboost [-m minsup] [-x maxpat] [-w wildcard] [-n v] [-e conv_epsilon] [-c coocitr] [-o] TrainingFile
$./eval_wild  model TestFile
