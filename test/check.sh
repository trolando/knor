#!/bin/bash

nuXmv=nuXmv

if ! [ -x "nuXmv" ]; then
    nuXmv=../test/nuXmv
fi

if ! [ -x "$nuXmv" ]; then
    echo "nuXmv not found! please place nuXmv binary on the path or in the test folder!"
    exit
fi

f=$1
shift 1
o=$*

echo "solving $f using options $o"
../build/knor $f $o -v -a 1> sol.aag
res=$?
if [[ $res == 10 ]]; then
    echo "preparing model checking..."
    cat $f | ../test/hoa2aig > monitor.aag
    ../test/combine-aiger monitor.aag sol.aag > combined.aag
    echo "running model checker..."
    echo "read_aiger_model -i combined.aag; encode_variables; build_boolean_model; check_ltlspec_ic3; quit" | timeout -k 10 10 $nuXmv -int > res 2>&1
    cat res | grep 'specification .* is '
else
    echo "knor reports unreal"
fi
