# Knor

Knor, a simple synthesis tool for HOA parity automata

**Currently implemented features**:
- parser for eHOA based on [HOA-TOOLS](https://github.com/gaperez64/hoa-tools/) by [Guillermo Perez](https://github.com/gaperez64/)
- BDD-based splitting of the parity automaton into a symbolic parity game using [Sylvan](https://github.com/trolando/sylvan/)
- solving the symbolic parity game explicitly using [Oink](https://github.com/trolando/oink/)
- solving the symbolic parity game symbolically using an internal solver based on [Sylvan](https://github.com/trolando/sylvan/)
- using symbolic bisimulation minimization to reduce the number of states, similar to [SigrefMC](https://github.com/trolando/sigrefmc/)
- synthesizing the controller as AIG from the strategy BDD either using ITE, or ISOP, or one-hot encoding
- post-synthesis compression using [ABC](https://github.com/berkeley-abc/abc)

## Usage

Compile Knor using the following command line:

```
mkdir build && cd build
cmake ..
make
```

Dependencies include `flex` and `bison` to generate the eHOA parser, and `boost`.

Use `knor --help` to get a list of parameters.
With `-v` or `--verbose`, Knor will write some information to standard error, such as the time it takes for each step of the process.
If invoked without a filename, Knor will attempt to read an eHOA file from `stdin`.

### Game construction
- Use `--naive` for the (old) naive splitting procedure, which is not recommended but possible for the purpose of comparing it with the BDD-based procedure.
- Use `--explicit` for the (also old) explicit parity game construction, also not recommended.
- **Recommendation**: no option

### Solving
- You can select a solver for the parity game solving, e.g., `--npp`, `--fpi`, `--fpj` are recommended solvers.  
  See further the documentation of [Oink](https://github.com/trolando/oink/).
- Use `--sym` to run a symbolic parity game solver that is tailor-made for Knor.
- Use `--print-game` to print the parity game to standard output instead of solving it.
- **Recommendation**: `--sym` unless the game is an artificial hard game, then `--tl`

### Synthesis
- Use `--real` to not do synthesis, only realizability.
- Use `--bisim` to minimize the state space using bisimulation minimization prior to synthesis.
- Use `--onehot` to encode the states using one-hot encoding instead of logarithmic encoding.
- Use `--isop` to use ZDD covers for the conversion to AIG.
- Use `--best` to find the smallest result of combinations of `--bisim`, `--isop` and `--onehot`.
- Use `--compress` to compress the AIG after construction using ABC.
- **Recommendation**: `--bisim --onehot` to get a fast result, or `--best --compress` to get the best result; `--bisim --onehot --compress` gives decent results too.

### Output
- Use `-a` to write the result to ascii AIGER format.
- Use `-b` to write the result to binary AIGER format.
- **Recommendation**: use `-a` for the SYNTCOMP competition and `-b` for smaller output size.

## Extended HOA

HOA is the [Hanoi Omega Automata](http://adl.github.io/hoaf/) file format that supports finite automata for languages on infinite words.
eHOA is HOA extended to record which atomic propositions are controllable, thus supporting synthesis.
See:
1. Guillermo A. Pérez (2019). The Extended HOA Format for Synthesis. In CoRR, [http://arxiv.org/abs/1912.05793](http://arxiv.org/abs/1912.05793).

## Splitting of the automaton

We use splitting of the parity automaton to a game using BDDs.
We use a variable ordering with uncontrollable APs before controllable APs.
For each state p, we do the following:
1. create an MTBDD encoding all transitions from p, where we encode the target state q and a priority in the MTBDD leaves.
2. create a set containing the unique MTBDDs for every valuation to uncontrollable variables (using straightforward enumeration in Sylvan)
3. for every unique MTBDD, compute all unique (priority, target state) pairs.
4. create a vertex owned by player 0 for every unique (priority, target state) pair, with the given priority, and the target state as successor state
5. create a vertex owned by player 0 for every MTBDD of step 3; from here, player 0 (controller) can choose any of the successor vertices created in step 4.
6. create a vertex owned by player 1 corresponding to p; from here, player 1 (environment) can choose any of the successor vertices created in step 5.

All vertices have priority 0 except for full-valuation vertices which have the priority (i.e. acceptance set) of the automaton's transition + 2. (The cases for min and odd parity automata are similar.)

## StarExec support

For the [SYNTCOMP competition](http://www.syntcomp.org/), the `prepare_se.sh` script prepares a gzipped tarball `knor.tar.gz` which can be uploaded to the StarExec environment.

## License

HOA-TOOLS is licensed under GPL v3 and Oink is licensed under Apache 2.0.
Thus Knor is licensed under GPL v3.
