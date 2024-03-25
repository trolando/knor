#!/bin/bash

nuXmv=nuXmv

if ! [ -x "nuXmv" ]; then
    nuXmv=../test/nuXmv
fi

if ! [ -x "$nuXmv" ]; then
    echo "nuXmv not found! please place nuXmv binary on the path or in the test folder!"
    exit
fi

for e in `ls -XrS ../examples/*hoa`; do
    f=$e
    #for s in sym tl; do
        #for o in "" "--onehot" "--bisim" "--isop" "--bisim --isop"; do
    for s in sym; do
        #for o in "--isop --onehot"; do
        for o in "--sop"; do
            echo -e "\033[32m$f, solver $s, options $o\033[m"
            ../build/knor $f --$s $o -v -a > $f-output-$s 2> $f-log
            res=$?
            if [[ $res == 10 ]]; then
                # echo "preparing model checking..."
                cat $f | ../test/hoa2aig > $f-monitor.aag
                # tail -n +2 $f-output-$s > $f-sol.aag
                ../test/combine-aiger $f-monitor.aag $f-output-$s > $f-combined.aag
                # echo "running model checker..."
                # cat "read_aiger_model -i $f-combined.aag; encode_variables; build_boolean_model; check_ltlspec_ic3; quit" 
                echo "read_aiger_model -i $f-combined.aag; encode_variables; build_boolean_model; check_ltlspec_ic3; quit" | timeout -k 10 10 $nuXmv -int > $f-res 2>&1
                # sed -n '/^nuXmv >/,$p' $f-res
                cat $f-res | grep 'specification .* is '
            else
                echo "knor reports unreal"
            fi
        done
    done
    echo
done
