
#data_list=('cpdb' 'mutag')
data_list=('mutag' '1nci1000' '47nci1000')
#maxpat_list=('7' '8' '9' '10' '11' '12')
maxpat_list=('6' '5' '4' '7' '8' '9' '10' '11' '12')
wildcard_list=('0' '1')
for data in ${data_list[@]}; do
    echo start ${data}
    python split_cross_validation.py data/${data}.gsp split/${data}
    for maxpat in ${maxpat_list[@]}; do
	for wildcard in ${wildcard_list[@]}; do
	    echo   w${wildcard} x$maxpat
	    for i in 0 1 2 3 4 5 6 7 8 9 ;do
		echo data $i
		time ./src/lpboost -w $wildcard -x $maxpat split/${data}train$i.gsp > log.txt
		./eval_wild/evaluator -f eval_wild/eval/${data}w${wildcard}x$maxpat$i.ev model split/${data}test$i.gsp > log
	    done
	    cd eval_wild/eval
	    python auc_eval.py "${data}w${wildcard}x${maxpat}*" > ${data}/w${wildcard}x$maxpat.txt
	    cd ../../
	done
    done
done
