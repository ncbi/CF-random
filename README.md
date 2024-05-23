# Data and code for "Sequence clustering confounds AlphaFold2"

Each directory contains all data and code used to generate each figure.

CF-random structure alternative predictions were generated with ColabFold1.5.5 as follows:

->  KaiB proteins: --num-seeds 33 --num-models 5 --model-type alphafold2 --max-seq 1 --max-extra-seq 2       
->  Mad2:          --num-seeds 19 --num-models 5 --model-type alphafold2 --max-seq 8 --max-extra-seq 16     
->  RfaH:          --num-seeds 25 --num-models 5 --model-type alphafold2 --max-seq 64 --max-extra-seq 128     

Programs were run on Linux and Mac OS 14.3.1.

Code for analysis written in python3 using the following libraries:

-> pandas     
-> numpy     
-> mdtraj     
-> matplotlib  
-> seaborn  
-> biopython


