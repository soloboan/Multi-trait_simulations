# Multi-trait_simulations

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
    B) TBV ('tbv/h' or 'tbv/l')  
    C) selection index ('index/h' or 'index/l')
- mating design for selected animals (md=)  
    A) random mating with random union of gamete (rnd_ug)  
    B) nested mating designs, i.e. specific mating ratio (e.g. 1:2 ; 1:25)  
- litter size per dam (ls=)  
- which of the traits you want to select when selection is not random (trsel=)  

### Additional fixed arguments  
- Sex ratio is fixed at 50% male and 50% female  
- G-struture and R-struture covariance structure are need in matrix format  





##### This code was written in close discussions with  
	- Panya Sae-Lim  
	- Binyam Dagnachew  
	- Kristine Hov Martinsen  
	- Bjarne Gjerde  

