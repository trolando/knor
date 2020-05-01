# HOA tools for synthesis
This is a set of utilities to work with the extended HOA format for reactive
synthesis. Part of the code consists in a _Flex_ + _Bison_ parser for the
format. The parser is used to provide:
1. A translator from extended HOA to PGSolver format
2. A translation to an Aiger file for model checking

## Dependencies
* If you are modifying the parser, you will need _Flex_ and _Bison_ or similar
  tools.

# TODO
- [ ] Test the translation to Aiger
- [ ] Plug in HOA to PGSolver translation from hoa2pg

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
