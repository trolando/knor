# Knor

Knor, a symbolic synthesis tool for HOA parity automata.

## Features

| Feature | Description |
|---------|-------------|
| eHOA parser | Reading eHOA files, parser based on [HOA-TOOLS](https://github.com/gaperez64/hoa-tools/) by [Guillermo Perez](https://github.com/gaperez64/) |
| Conversion to game | Three different methods to convert the parity automaton into a parity game: (1) a naive explicit construction, (2) a BDD-supported explicit construction, (3) a fully symbolic construction |
| Explicit solving | Solving the symbolic parity game explicitly using [Oink](https://github.com/trolando/oink/) |
| Symbolic solving | Solving the symbolic parity game symbolically using an internal solver based on [Sylvan](https://github.com/trolando/sylvan/) |
| Bisimulation minimization | Using symbolic bisimulation minimization to reduce the number of states after solving the parity game, similar to [SigrefMC](https://github.com/trolando/sigrefmc/) |
| Translation to controller | Synthesizing the controller as AIG from the strategy BDD either using ITE or ISOP to compute the logic network, and either a logarithmic or a one-hot encoding of the state |
| Post-synthesis compression | Compressing the AIG after synthesis using [ABC](https://github.com/berkeley-abc/abc) |

## Building Knor

Knor depends on `flex` and `bison` to generate the eHOA parser, as well as on the `boost` libraries for C++.

Knor also downloads and builds [Oink](https://github.com/trolando/oink), [Sylvan](https://github.com/trolando/sylvan) and [ABC](https://github.com/berkeley-abc/abc) automatically.

Compile Knor using the following command line:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Using Knor

Use `knor --help` to get a list of parameters.
With `-v` or `--verbose`, Knor will write some information to standard error, such as the time it takes for each step of the process.
If invoked without a filename, Knor will attempt to read an eHOA file from `stdin`.
Knor returns value 10 if the specification is realizable, and value 20 if the specification is unrealizable.

Example command line:
```bash
./knor ../examples/amba_decomposed_arbiter_10.tlsf.ehoa --bisim --onehot --compress -a -v > controller.aag
```

### Game construction
- Use `--naive` for the naive explicit splitting procedure, which is not recommended but possible for the purpose of comparing it with the BDD-based procedures.
- Use `--explicit` for the BDD-supported explicit parity game construction, which is not recommended but possible for the purpose of comparing it with the fully symbolic construction.
- By default, the fully symbolic construction is used.
- **Recommendation**: no option, so the fully symbolic construction is used.

### Printing the game
- Use `--print-game` to print the constructed parity game to standard output, and terminate instead of continuing with solving the parity game.

### Solving the game
- Use `--sym` to run the internal symbolic parity game solver. This is the default option; using `--sym` is optional.
- You can select an explicit solver for the parity game solving, e.g., `--npp`, `--fpi`, `--fpj` are recommended solvers for most games, but `--tl`, `--rtl` are suitable for artificial hard games. See further the documentation of [Oink](https://github.com/trolando/oink/).
- **Recommendation**: `--sym` unless the game is an artificial hard game, then `--tl` or `--rtl`.

### Controller synthesis
- Use `--real` to not do synthesis, only check realizability.
- Use `--bisim` to minimize the state space using bisimulation minimization prior to synthesis.
- Use `--onehot` to encode the states using one-hot encoding instead of logarithmic encoding.
- Use `--isop` to use ZDD covers for the conversion to AIG. By default, a direct conversion based on ITE (Shannon expansion) is used.
- Use `--best` to automatically find the smallest result of combinations of `--bisim`, `--isop` and `--onehot`.
- Use `--compress` to compress the AIG after construction using ABC.
- **Recommendation**: `--bisim --onehot` to get a fast result, or `--best --compress` to get the best result; `--bisim --onehot --compress` gives reasonably good results.

### Output of AIG controller
- Use `-a` to write the result to ascii AIGER format.
- Use `-b` to write the result to binary AIGER format.
- **Recommendation**: use `-a` for the SYNTCOMP competition and `-b` for smaller output size.

## Extended HOA

HOA is the [Hanoi Omega Automata](http://adl.github.io/hoaf/) file format that supports finite automata for languages on infinite words.
eHOA is HOA extended to record which atomic propositions are controllable, thus supporting synthesis.
See:
1. Guillermo A. Pérez (2019). The Extended HOA Format for Synthesis. In CoRR, [http://arxiv.org/abs/1912.05793](http://arxiv.org/abs/1912.05793).

## StarExec support

For the [SYNTCOMP competition](http://www.syntcomp.org/), the `prepare_se.sh` script prepares a gzipped tarball `knor.tar.gz` which can be uploaded to the StarExec environment.

## License

HOA-TOOLS is licensed under GPL v3. Sylvan and Oink are licensed under Apache 2.0.
Thus Knor is licensed under GPL v3.
