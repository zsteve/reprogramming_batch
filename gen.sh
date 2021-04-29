lamda_vals=(0.0025 0.005 0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28)
for lamda in ${lamda_vals[@]}; do 
    dirname="run_"$lamda"_dirichlet5"
    mkdir $dirname
    cd $dirname
    cp ../run.sh .
    ln -s ../batch.py .
    sed -i "s/__LAMDA__/$lamda/g" run.sh
    for i in $(seq 1 25); do 
        qsub run.sh
    done
    cd ..
done
