# Test tools

The model-checking scripts `check.sh` and `test.sh` use two prebuilt
(Linux x86-64) third-party tools. They are not part of knor and keep their
own licenses:

| Tool | Source | License |
|------|--------|---------|
| `hoa2aig` | `hoa2aig.c` from [HOA-TOOLS](https://github.com/SYNTCOMP/hoa-tools), Copyright (c) 2019- Guillermo A. Pérez | GPL v3 or later, see [hoa2aig.LICENSE](hoa2aig.LICENSE) |
| `combine-aiger` | [combine-aiger](https://github.com/SYNTCOMP/combine-aiger), Copyright (c) 2016 Leander Tentrup, Saarland University | ISC, see [combine-aiger.LICENSE](combine-aiger.LICENSE) |

Both include the AIGER library by Armin Biere, whose license is in the
header of `src/aiger.h`. The complete source code of both tools is
available at the links above.
