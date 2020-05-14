# HOA tools for synthesis
This is a set of utilities to work with the extended HOA format for reactive
synthesis. Part of the code consists in a _Flex_ + _Bison_ parser for the
format. The parser is used to provide:
1. A translator from extended HOA to PGSolver format
2. A translation to an Aiger file for model checking

Additionally, some sample extended-HOA files are included. Most were generated
using scripts from https://github.com/gaperez64/tlsf2gpg and dumping an
extended HOA file instead of going to a parity game.

## HOA2PG
TODO: port https://github.com/gaperez64/hoa2pg to here in order to use the new
HOA parser

## HOA2AIG
The translation from a HOA file, specifying a *complete* and *deterministic*
parity automaton, to an Aiger file is based on [1]. Note that NuSMV includes
a similar utility although both the input and output formats differ from what
we need in this case: LTL vs HOA and SMV vs Aiger.

The main idea is to encode
the transitions of the automaton in the Aiger in a natural way while also
adding fairness constraints to make sure the parity objective is correctly
observed. More precisely, for each possible maximal/minimal priority which
is witnessed infinitely often, we add justice constraints to the Aiger file --
these are conjunctions of fairness constraints, see [2]. To make sure that
each justice constraint (respectively, priority) occurring infinitely often,
is not also "trumped" infinitely often by others, we add an input and a
"bad state" (i.e. a safety) constraint. The latter allow us to require that
the "trumping" occurs only finitely often. (The reset input essentially marks
the point in time after which the system promises eventual safety is in fact
safety.)

1. Clarke, E., Grumberg, O., & Hamaguchi, K. (1994, June). Another look at LTL
   model checking. In _International Conference on Computer Aided
   Verification_
   (pp. 415-427). Springer, Berlin, Heidelberg.
2. Biere, A., Heljanko, K., & Wieringa, S. (2011). AIGER 1.9 and beyond. 
   _Available at fmv.jku.at/hwmcc11/beyond1.pdf._

## Dependencies
* If you are modifying the parser, you will need _Flex_ and _Bison_ or similar
  tools.

# Citing
If you use HOA tools for your academic work, please cite the report regarding
the extended HOA format. For the translation to PGSolver format from HOA,
please cite our _Reachability Problems_ paper on generalized parity games,
where a first version of the translation was featured.

```
@article{DBLP:journals/corr/abs-1912-05793,
  author    = {Guillermo A. P{\'{e}}rez},
  title     = {The Extended {HOA} Format for Synthesis},
  journal   = {CoRR},
  volume    = {abs/1912.05793},
  year      = {2019},
  url       = {http://arxiv.org/abs/1912.05793},
  archivePrefix = {arXiv},
  eprint    = {1912.05793},
  timestamp = {Tue, 03 Mar 2020 16:02:55 +0100},
  biburl    = {https://dblp.org/rec/journals/corr/abs-1912-05793.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}

@inproceedings{DBLP:conf/rp/BruyerePRT19,
  author    = {V{\'{e}}ronique Bruy{\`{e}}re and
               Guillermo A. P{\'{e}}rez and
               Jean{-}Fran{\c{c}}ois Raskin and
               Cl{\'{e}}ment Tamines},
  editor    = {Emmanuel Filiot and
               Rapha{\"{e}}l M. Jungers and
               Igor Potapov},
  title     = {Partial Solvers for Generalized Parity Games},
  booktitle = {Reachability Problems - 13th International Conference, {RP} 2019,
               Brussels, Belgium, September 11-13, 2019, Proceedings},
  series    = {Lecture Notes in Computer Science},
  volume    = {11674},
  pages     = {63--78},
  publisher = {Springer},
  year      = {2019},
  url       = {https://doi.org/10.1007/978-3-030-30806-3\_6},
  doi       = {10.1007/978-3-030-30806-3\_6},
  timestamp = {Mon, 09 Sep 2019 15:37:02 +0200},
  biburl    = {https://dblp.org/rec/conf/rp/BruyerePRT19.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
