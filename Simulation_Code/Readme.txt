Installation:
--------------
	1) Install python 3
	2) Install required packages: requirment.txt 
		-to install reqirment.txt run 'pip install -r reqirment.txt'

File list
----------
	main.py
	sim.py
	addNode.py
	Readme.txt 
	requirment.txt

How to run
-----------
	----------------------
	Command Line Arguments
	----------------------
	python main.py [-h] -start_d START_DEVICE [-end_d END_DEVICE]
               [-ncellusr N_CELL_USERS] [-pcell PWR_CELL_DB]
               [-pd2d PWR_D2D_DB] [-neve N_EAVESDROPPER] -ncpu NCPU -nsubslot
               NSUBSLOT -seed SEED [-rdth RDTHRESHOLD] -exp EXPERIMENT

    required arguments:
    --------------------
    -start_d 	:	 K starting D2D device range
    -ncpu       :	 number of cpu for parallel processing
    -seed       :	 end range of numpy random seed to average
    -exp        :	 experiment-> 1,2 or 3(1:Amal project, 2:Senario 1, 3:Senario 2)
		                    
    optional arguments:
    --------------------

	-end_d 	   :	ending D2D device range (PS: end device excluded)
	-ncellusr  :    number of cellular user (default: 5)
	-pcell 	   : 	Power of cellular user (default: 10)
	-pd2d 	   : 	Power of D2D user (default: 5)
	-neve      :  	number of eavesdroppers (default: 2)
	-rdth 	   : 	Percentage of RD Max threshold

Example 1: 
-----------
    python main.py -start_d=6 -nsubslots=2 -ncpu=2 -seed=10 -exp=3

    (K=6, number of sub-slot=2, average of 10 seeds)

Example 2: 
-----------
	python main.py -start_d=6 -nsubslots=2 -ncpu=2 -seed=10 -exp=3 -rdth=20

    (K=6, number of sub-slot=2, average of 10 seeds, MaxRD threshold 20%)

Example 3: 
-----------
	python main.py -start_d=6 -end_d=9 -nsubslots=2 -ncpu=2 -seed=10 -exp=3 -rdth=20

    (for K=6,7,8  with number of sub-slot=2, average of 10 seeds, MaxRD threshold=20)

How to read output file
------------------------

The results will be saved as multiple csv files which could be categorized as follows:
	1. csv files starting with A: per seed results.
	2. csv files starting with E: expection across the input seeds.
Results include - K, KS, KJ, R_throughput, RC_sm, Rcsm_noD2D, R_th, Time taken.

