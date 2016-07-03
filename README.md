# graph_co-occurence_boosting

This package is gBoost with wildcard

- Prepare
need to use:
c++ library, boost
-->install on Ubuntu
     sudo apt-get install libboost-dev

- How to compile
How to make:
$ cd src
$ make
$ cd ../eval_wild
$ make
$ cd ..

- How to use
$./lpboost [-m minsup] [-x maxpat] [-w wildcard] [-n v] [-e conv_epsilon] [-c coocitr] [-o] TrainingFile

$./eval_wild  model TestFile
