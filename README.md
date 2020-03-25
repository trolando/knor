# HOA tools for synthesis
This is a set of utilities to work with the extended HOA format for reactive
synthesis. Most of the code consists in a _Flex_ + _Bison_ parser for the
format. The parser is used to provide:
1. A translator from extended HOA to PGSolver format
2. A translation to an Aiger file for model checking

## Dependencies
* If you are modifying the parser, you will need _Flex_ and _Bison_ or similar
  tools.

# TODO
[ ] Implement the translation to Aiger
[ ] Plug in HOA to PGSolver translation from hoa2pg
