# epLSAP-Aligner
The implementation of "epLSAP-Align: A non-sequential protein structural alignment solver with entropy-regularized partial linear sum assignment
problem formulation".     
            
## Python Dependencies:
For the Python implementation, we highly recommond using Anaconda enviornment. Suppose that the conda enviornment name is *eplsap*, then:

`conda create --name eplsap python=3.9`

After constructing conda enviroment, 

`pip install biopython`

`pip install pymican`

`pip install tmtools`

`pip install scipy`

or directly

`pip install -r requirements.txt`
            
By leveraging the concise grammar style of Python programming, the Python scripts briefly shows how the proposed non-sequential works with the superposition produced by TMalign and MICAN.     

### Usage 
To test with query protein "a.pdb" and the target protein "b.pdb" with TMalign/MICAN superposition, one can type   
`python score_test.py --q_pdb a.pdb --t_pdb b.pdb --mode TM/mican`

For example, you can try the pdb example included in the file as follows:
`python score_test.py --q_pdb ./USalign-epLSAP/PDB1.pdb --t_pdb ./USalign-epLSAP/PDB2.pdb --mode TM`

The output of the main alignment function, epLSAP_Main in epLSAP_fun.py, is Number of align (OT_Nali), RMSD (OT_RMSD), TM-score (OT_TMscore), structure overlap (OT_SO). If parameter fast_opt is set true (default false), the residue-to-residue correspondence will be the last output.

#### ***Notice:***   We have tested the Python version on MacOS and Ubuntu system and it works well, the function refer to system difference is the Python link of C/C++ codes "TMalign_tools.cpp" for TM-score calculation. The .so links are generated in the Github Actions part. If you meet bugs with other systems, please raise issues. 

## C/C++ implementation:
The C/C++ version of epLSAP-TM is also implemented in the File USalign-epLSAP, where the codes are built based on USalign2 (https://github.com/pylelab/USalign).    
In this C/C++ implementation, in the "SOalign.h" file, we replace the EGS solver of USalign2 with our epLSAP solver, and use SP-score instead of TM-score. The usage of epLSAP-TM is the same as USalign2 NS mode.  


