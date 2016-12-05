# Single and Multi-trait_simulations (Polygenic/Pedigree approach)

### other packages the script depends on (please install them)  
    - MASS  
    - pedigree  

Simple R-codes (2 scripts) for multi-trait simulations

### simulating a base population
    -function : 'makebasepop'

### simulation the offspring population (desired number of generations)
    -function : 'makeoff'

## Options
- User defined number of traits  
- User defined number of generations (Numgen=)  
- selection of parents based on (sd=)   
    A) random - no selection (rnd)  
    B) Phenotypic selection (Phen/h or Phen/low)
    C) TBV ('tbv/h' or 'tbv/l')  
    D) selection index based on TBV ('index/h' or 'index/l')  
    D) selection index based on phenotypes ('phenindex/h' or 'phenindex/l')  
- mating design for selected animals (md=)  
    A) random mating with random union of gamete (rnd_ug) [leads to pseudo-overlapping generations: mostly used in cattle breeding] 
    B) nested mating designs, i.e. specific mating ratio (e.g. 1:2 ; 1:25) [purely non-overlapping generations: mostly used in fish breeding]  
- litter size per dam (ls=)  
- which of the traits you want to select when selection is not random (trsel=)  

### Additional fixed arguments  
- Sex ratio is fixed at 50% male and 50% female  
- G and R covariance structure are needed in matrix format (even for a single triat simulation)  

### Output columns 
    - Generation
    - Pedigree structure (ID, Sire, Dam)
    - Sex
    - Inbreeding (Fped)
    - True breeding values (number of columns equals number of traits)
    - Residuals (number of columns equals number of traits)
    - Phenotypes ((number of columns equals number of traits))


##### This code was written in close discussions with  
	- Panya Sae-Lim  (panya.sae-lim@nofima.no)
	- Binyam Dagnachew  
	- Kristine Hov Martinsen  
	- Bjarne Gjerde  

