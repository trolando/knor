# Knor

Knor, a symbolic synthesis tool for HOA parity automata.

## Features

| Feature | Description |
|---------|-------------|
| eHOA parser | Reading eHOA files: the [HOA format](https://adl.github.io/hoaf/) extended with `controllable-AP`, using a fast hand-written parser |
| Conversion to parity | Deterministic automata with other acceptance conditions are converted to parity automata: conditions that are parity conditions in disguise (Büchi, co-Büchi, Rabin/Streett with one pair, ...) automatically, and all others (generalized Büchi, Rabin, Streett, Muller, ...) with `--paritize`, via a product with the Zielonka tree of the condition |
| Conversion to game | Three different methods to convert the parity automaton into a parity game: (1) a naive explicit construction, (2) a BDD-supported explicit construction, (3) a fully symbolic construction |
| Explicit solving | Solving the symbolic parity game explicitly using [Oink](https://github.com/trolando/oink/) |
| Symbolic solving | Solving the symbolic parity game symbolically using an internal solver based on [Sylvan](https://github.com/trolando/sylvan/) |
| Bisimulation minimization | Using symbolic bisimulation minimization to reduce the number of states after solving the parity game, similar to [SigrefMC](https://github.com/trolando/sigrefmc/) |
| Translation to controller | Synthesizing the controller as AIG from the strategy BDD either using ITE or ISOP to compute the logic network, and either a logarithmic or a one-hot encoding of the state |
| Post-synthesis compression | Compressing the AIG after synthesis using [ABC](https://github.com/berkeley-abc/abc) |

## Building Knor

Knor depends on the `boost` libraries for C++.

Knor also downloads and builds [Oink](https://github.com/trolando/oink), [Sylvan](https://github.com/trolando/sylvan) and [ABC](https://github.com/berkeley-abc/abc) automatically.

Compile Knor using the following command line:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Run the regression tests (verdicts of all examples, parser, acceptance conditions) with `make test` in the build directory.

## Using Knor

Use `knor --help` to get a list of parameters.
With `-v` or `--verbose`, Knor will write some information to standard error, such as the time it takes for each step of the process.
If invoked without a filename, Knor will attempt to read an eHOA file from `stdin`.
Knor returns value 10 if the specification is realizable, and value 20 if the specification is unrealizable.

### Input automata
- The input must be a complete deterministic automaton with a single initial state.
- Acceptance conditions that a priority per edge can express on the same automaton are converted to `parity max even` automatically, for example Büchi, co-Büchi, `t`, `f`, Rabin and Streett with one pair, parity conditions in any notation, and edges in zero or several acceptance sets.
- Use `--paritize` to also convert all other acceptance conditions, such as generalized Büchi or Rabin and Streett with several pairs. This multiplies the automaton with the Zielonka tree automaton of the acceptance condition. Without `--paritize`, such inputs are rejected with exit code 104.

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
- By default, the game is solved by an explicit solver of [Oink](https://github.com/trolando/oink/), tangle learning (`--tl`). You can select another one, e.g., `--npp`, `--fpi`, `--fpj`; see further the documentation of Oink.
- Use `--sym` to run the internal symbolic parity game solver instead.
- **Recommendation**: the default.

### Controller synthesis
- Use `--real` to not do synthesis, only check realizability.
- Use `--bisim` to minimize the state space using bisimulation minimization prior to synthesis.
- Use `--onehot` to encode the states using one-hot encoding instead of logarithmic encoding.
- Use `--isop` to use ZDD covers for the conversion to AIG. By default, a direct conversion based on ITE (Shannon expansion) is used.
- Use `--auto-encode` for a small controller: the multi-strategy of fixpoint iteration (`--fpim`), minimisation of the controller as a Mealy machine with a SAT solver (`--mealy-sat`), encoding with don't cares (`--dc`), and for small controllers a short search over all winning moves. It works with an explicit solver; `--fpjm` may replace `--fpim`.
- Use `--compress` to compress the AIG after construction using ABC.
- **Recommendation**: `--auto-encode` for the smallest result, optionally with `--compress`; `--bisim --onehot` for a fast result.

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

Knor is licensed under the [Apache License 2.0](LICENSE).
Its dependencies Lace, Sylvan and Oink are licensed under Apache 2.0 as well, ABC under a BSD-style license, and the bundled AIGER library (`src/aiger.c`) and cxxopts (`src/cxxopts.hpp`) under MIT-style licenses.
The third-party test tools in `test/` keep their own licenses, see [test/README.md](test/README.md).
