# Knor

Knor, a simple synthesis tool for HOA parity automata

**Currently implemented features**:
- parser for eHOA based on [HOA-TOOLS](https://github.com/gaperez64/hoa-tools/) by [Guillermo Perez](https://github.com/gaperez64/)
- naive splitting of the parity automaton into a parity game
- solving the parity game using [Oink](https://github.com/trolando/oink/)

Knor does not yet support synthesizing the controller, and only computes whether or not the controller is realizable.

## Usage

Compile Knor using the following command line:

```
mkdir build && cd build
cmake ..
make
```

Dependencies include `flex` and `bison` to generate the eHOA parser, and `boost` for Oink.

Knor takes two parameters: `knor <filename> [solver]` where `<filename>` points at an eHOA file and the optional parameter `[solver]` selects a parity game solver implemented in Oink, e.g., `npp`, `fpi`, `fpj`, `zlk`, `psi`, etc.
See further the documentation of [Oink](https://github.com/trolando/oink/).
If invoked without parameters, Knor will attempt to read an eHOA file from `stdin`.

## Extended HOA

HOA is the [Hanoi Omega Automata](http://adl.github.io/hoaf/) file format that supports finite automata for languages on infinite words.
eHOA is HOA extended to record which atomic propositions are controllable, thus supporting synthesis.
See:
1. Guillermo A. Pérez (2019). The Extended HOA Format for Synthesis. In CoRR, [http://arxiv.org/abs/1912.05793](http://arxiv.org/abs/1912.05793).

## Splitting of the automaton

We use a naive splitting of the parity automaton to a game, as in HOA-TOOLS.
For each transition from p to q, via a valuation Val, we add vertices to the game:
1. a vertex owned by player 1 corresponding to p; from here, player 1 (environment) can choose
   a valuation of the uncontrollable atomic propositions;
2. a vertex owned by player 0 for each such partial valuation; from here,
   player 0 (controller) can choose a valuation of the controllable APs;
3. a vertex owned by player 0 for each such full valuation Val; from here there is no choice but to advance to q.

All vertices have priority 0 except for full-valuation vertices which have the
priority (i.e. acceptance set) of the automaton's transition + 2. (The cases
for min and odd parity automata are similar.)

This procedure is the most naive method to split an automaton into a game and likely much smarter methods exist. Nevertheless, the construction results in a game that increases the size of the automaton by at most 2 vertices and 2 edges for every edge in the original automaton. Hence a quadratic size increase.

## StarExec support

For the [SYNTCOMP competition](http://www.syntcomp.org/), the `prepare_se.sh` script prepares a gzipped tarball `knor.tar.gz` which can be uploaded to the StarExec environment.

## License

HOA-TOOLS is licensed under GPL v3 and Oink is licensed under Apache 2.0.
Thus Knor is licensed under GPL v3.
