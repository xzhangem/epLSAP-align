import sys
from Bio import PDB 
from Bio.PDB.Polypeptide import is_aa
import epLSAP_fun
import csv
import os 
import numpy as np
import math
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--q_pdb', type=str, help='filename of query pdb')
    parser.add_argument('--t_pdb', type=str, help='filename of target pdb')
    parser.add_argument('--mode', type=str, default='TM', help='using TM or mican for superposition')
    args = parser.parse_args()

    OT_Nali, OT_RMSD, OT_tmscore, OT_SO = epLSAP_fun.epLSAP_Main(args.q_pdb, args.t_pdb, False, args.mode)
    
    print("OT Nali: {}; OT RMSD: {}; OT tmscore: {}; OT SO: {}".format(OT_Nali, OT_RMSD, OT_tmscore, OT_SO))
