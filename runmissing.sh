for i in $(ls -d run_* | awk -F'_' '{ print $2 }' | awk -F'.sh' '{ print $1 }'); do
    dirname="run_"$i"_dirichlet5"
    num=$(ls $dirname/*.out | wc -l)
    cd $dirname
    for i in $(seq $num 24); do
        qsub run.sh
    done
    cd ..
    # if test -f $outputname; then 
    #     echo "$outputname exists"; 
    # else 
    #     echo "$outputname does not exist. re-qsubbing..."; 
    #     qsub $runname
    # fi ;
done
