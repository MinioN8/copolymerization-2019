These codes form the basis of data generation that was done in the published article:

Lynd, Nathaniel A., Robert C. Ferrier, and Bryan S. Beckingham. “Recommendation for Accurate Experimental Determination of Reactivity Ratios in Chain Copolymerization.” Macromolecules 52, no. 6 (March 26, 2019): 2277–85. https://doi.org/10.1021/acs.macromol.8b01752.

If you use these codes or modify them, please cite the original article. There are three codes in this repository. ScriptGen.cpp creates command-line input scripts for AB-terminal.cpp (in a *nix environment) to numerically integrate the coupled system of differential equations and create output files containing the data. There are Mathematica notebooks that read and interpret the data. AB-Terminal.init.cpp was an unpublished (not peer-reviewed) study where I investigated if the reason for the copolymer equation's inaccuracy had to do with my data explicitely modeling initiation whereas the copolymer equation should only apply to the steady-state regime of a copolymerization. Turns out it was a negative result. Initiation was not the problem. The problem is the copolymer equation and that it becomes increasingly invalid as a copolymerization proceeds.

Mayo, Frank R., and Frederick M. Lewis. “Copolymerization. I. A Basis for Comparing the Behavior of Monomers in Copolymerization; The Copolymerization of Styrene and Methyl Methacrylate.” Journal of the American Chemical Society 66, no. 9 (September 1944): 1594–1601. https://doi.org/10.1021/ja01237a052.

If you have any further requests, please send them to lynd@che.utexas.edu and I'd be happy to accommodate reasonable requests.
