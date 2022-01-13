# Welcome
core package, as an inter-lingua to work with sim, FTS, and potentially other simulation packages.

can see some more [documentation](http://equalapriori.github.io/sim-utils/).
- (which should be up to date for the `main` branch)

Current planned sketch:
IO:
- Files ←→ dictionaries
- Force field
- Topology
- Mapping

Core: cgfts
- System object: holds FF, Topology, Misc settings (e.g. CGMD or PFTS)
Moleculetype dictionary
- FF object (Charles)
  - infer bead types from ff definition
- Topology module? (My, Kevin)
  - Beadtype names, bonds, maybe: alternative comb definitions/helper functions
- Tosim
- Topfts
- *parameters* for bead types should be in force field. Topology simply directly uses bead names

Core Tools?
- mapping
- Creating Combs

Non-core tools? (maybe another package, mdfts_utils)
- Running FTS?
- Computations?



