
data_list=('cpdb' 'mutag')
#data_list=('cpdb_seed100_')
#data_list=('mutag' '1nci1000' '47nci1000')
#maxpat_list=('7' '8' '9' '10' '11' '12')
maxpat_list=('8' '9' '10')
c_list=(' ' '-o' '-c 1')
#nu_list=('01' '1' '2' '3' '4' '5' '6')
nu_list=('7' '8' '9')
#nu_list=('4')
#wildcard_list=('0' '1')

echo cooc after conv # -o
for data in ${data_list[@]}; do
    echo start ${data}
    for maxpat in ${maxpat_list[@]}; do
	for nu in ${nu_list[@]};do
	    echo $data nu $nu x$maxpat   "-o"
	    accs=0.0
	    aucs=0.0
	    for i in 0 1 2 3 4 5 6 7 8 9 ;do
		echo -n $i
	        ./../src/lpboost -x $maxpat -n 0.$nu -o ../split/${data}train$i.gsp > log.txt
		./../eval_wild/evaluator -f ../eval_wild/eval/${data}on${nu}x$maxpat$i.ev model ../split/${data}test$i.gsp > log
	    done
	    cd ../eval_wild/eval
	    for i in 0 1 2 3 4 5 6 7 8 9 ;do
		#echo $i
		acc=`python auc_eval.py ${data}on${nu}x${maxpat}$i.ev  0`
		auc=`python auc_eval.py ${data}on${nu}x${maxpat}$i.ev  1`
		#echo $acc $auc
		#echo "scale=7;$accs + $acc" | bc
		accs=`echo "scale=14;$accs + $acc" | bc`
		aucs=`echo "scale=14;$aucs + $auc" | bc`
	    done
	    echo
	    #echo $accs $aucs
	    accs=`echo "scale=14;${accs}/10"|bc`
	    aucs=`echo "scale=14;${aucs}/10"|bc`
	    echo 0$accs 0$aucs #>../../Result_cooc/${data}/on${nu}x$maxpat.txt
	    cd ../../tools	    
	done
    done
done
<<EOF
echo cooc  # -c 1
for data in ${data_list[@]}; do
    echo start ${data}
    for maxpat in ${maxpat_list[@]}; do
	for nu in ${nu_list[@]};do
	    for i in 0 1 2 3 4 5 6 7 8 9 ;do
		echo $data nu $nu x$maxpat  $i "-c 1 "
		time ./../src/lpboost -x $maxpat -n 0.$nu -c 1 ../split/${data}train$i.gsp > log.txt
		./../eval_wild/evaluator -f ../eval_wild/eval/${data}c1n${nu}x$maxpat$i.ev model ../split/${data}test$i.gsp > log
	    done
	    cd ../eval_wild/eval
	    python auc_eval.py "${data}c1n${nu}x$maxpat*" > ../../Result_cooc/${data}/c1n${nu}x$maxpat.txt
	    cd ../../tools     
	done
    done
done
EOF
