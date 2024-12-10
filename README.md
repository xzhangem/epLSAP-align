# epLSAP-Aligner
The implementation of "Accurate and fast non-sequential protein structural alignment solver with entropy-regularized partial linear sum assignment problem formulation" (epLSAP-Aligner).     
                
The dependencies and libraries needed for the Python script is listed in requirement.yaml.      
            
By leveraging the concise grammar style of Python programming, the Python scripts briefly shows how the proposed non-sequential works with the superposition produced by TMalign and MICAN.     

To test with query protein "a.pdb" and the target protein "b.pdb" with TMalign/MICAN superposition, one can type   
`python score_test.py --q_pdb a.pdb --t_pdb b.pdb --mode TM/mican`

As the example of integrating epLSAP-aligner into other the existing NS alignment in a computationally efficient way, currently we use USalign2.    
We implement epLSAP-aligner to replace EGS in original USalign2. The implementations are mainly contained in "SOIalign.h" file. 

