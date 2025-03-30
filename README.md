# epLSAP-Aligner
The implementation of "epLSAP-Align: A non-sequential protein structural alignment solver with entropy-regularized partial linear sum assignment
problem formulation".     
            
## ***Python Dependencies:***
For the Python implementation, we highly recommond using Anaconda enviornment. Suppose that the conda enviornment name is *eplsap*, then:
`conda create --name eplsap python=3.9`

            
By leveraging the concise grammar style of Python programming, the Python scripts briefly shows how the proposed non-sequential works with the superposition produced by TMalign and MICAN.     

To test with query protein "a.pdb" and the target protein "b.pdb" with TMalign/MICAN superposition, one can type   
`python score_test.py --q_pdb a.pdb --t_pdb b.pdb --mode TM/mican`

The C/C++ version of epLSAP-TM is also implemented in the File USalign-epLSAP, where the codes are built based on USalign2 (https://github.com/pylelab/USalign).    
In this C/C++ implementation, in the "SOalign.h" file, we replace the EGS solver of USalign2 with our epLSAP solver, and use SP-score instead of TM-score. The usage of epLSAP-TM is the same as USalign2 NS mode.  

