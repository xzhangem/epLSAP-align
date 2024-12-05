import numpy as np
import math
import os
import TMalign_tools
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from Bio.SVDSuperimposer import SVDSuperimposer
from pymican import mican
from tmtools import tm_align
from scipy.spatial.transform import Rotation as R

from Bio import PDB, SeqIO
from scipy import io


def parameter_set4search(xlen, ylen):
    D0_MIN = 0.5
    dcu0 = 4.25
    Lnorm = min(xlen, ylen)
    if (Lnorm <= 19):
        d0 = 0.168
    else:
        d0 = 1.24 * (Lnorm - 15) ** (1.0 / 3.0) - 1.8
    D0_MIN = d0 + 0.8
    d0_search = D0_MIN
    if (d0_search > 8.0):
        d0_search = 8.0
    if (d0_search < 4.5):
        d0_search = 4.5
    score_d8 = (1.5 * (Lnorm * 1.0) ** (0.3)) + 3.5

    return D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0


def fasta_align(fasta_file, pdb1, pdb2):
    coord1, _ = CA_coord(pdb1)
    coord2, _ = CA_coord(pdb2)
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    fasta_info = np.loadtxt(fasta_file, dtype=str)
    fasta_seq_1 = fasta_info[0]
    fasta_seq_2 = fasta_info[1]
    match_indx_1 = []
    match_indx_2 = []
    tag_1 = 0
    tag_2 = 0
    for i in range(len(fasta_seq_1)):
        if fasta_seq_1[i] != '-' and fasta_seq_2[i] != '-':
            match_indx_1.append(tag_1)
            match_indx_2.append(tag_2)
            tag_1 += 1
            tag_2 += 1
        elif fasta_seq_1[i] == '-' and fasta_seq_2[i] != '-':
            tag_2 += 1
        elif fasta_seq_1[i] != '-' and fasta_seq_2[i] == '-':
            tag_1 += 1
    return superpose_set(coord1[match_indx_1, :], coord2[match_indx_2, :]), match_indx_1, match_indx_2



def lsap2ot(score_mat, gap, epsilon, lmbd, T):
    gap = -abs(gap)
    c = np.max(score_mat) + epsilon
    m, n = score_mat.shape
    S_mat = np.zeros((m+1, n+1))
    S_mat[0:m, 0:n] = score_mat
    S_mat[m, :] = gap
    S_mat[:, n] = gap
    S_mat[m, n] = 0
    D_mat = np.zeros((m+1, n+1))
    
    D_mat[0:m, 0:n] = 2 * c
    D_mat[m, :] = c
    D_mat[:, n] = c
    D_mat[m,n] = 0

    C_mat = D_mat - S_mat

    # Initialize 
    K_mat = np.exp(-lmbd * C_mat)
    y = np.ones((n+1, 1))
    conv = False
    i = 0
    x = np.ones((m+1,1))

    while (conv == False and i < T):
        xp = 1 / (np.dot(K_mat, y))
        xp[m,0] = 1
        if i > 0:
            err_x = ((xp / x - np.ones((m+1,1))) * (xp / x - np.ones((m+1,1)))).sum()
        x = xp
        yp = 1 / (np.dot(K_mat.T, x))
        yp[n,0] = 1
        err_y = ((yp / y - np.ones((n+1,1))) * (yp / y - np.ones((n+1,1)))).sum()
        y = yp
        if i > 0 and min(err_x, err_y) < epsilon:
            conv = True
        i += 1
    #print(i)

    x = x.squeeze()
    y = y.squeeze()
    result = np.dot(np.diag(x), K_mat)
    result = np.dot(result, np.diag(y))

    return result

def parameter_set4search(xlen, ylen):
    D0_MIN = 0.5
    dcu0 = 4.25
    Lnorm = min(xlen, ylen)
    if (Lnorm <= 19):
        d0 = 0.168
    else:
        d0 = 1.24 * (Lnorm - 15) ** (1.0 / 3.0) - 1.8
    D0_MIN = d0 + 0.8
    d0_search = D0_MIN
    if (d0_search > 8.0):
        d0_search = 8.0
    if (d0_search < 4.5):
        d0_search = 4.5
    score_d8 = (1.5 * (Lnorm * 1.0) ** (0.3)) + 3.5

    return D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0

def parameter_set4final(length):
    D0_MIN = 0.5
    Lnorm = length
    if (Lnorm <= 21):
        d0 = 0.5
    else:
        d0 = 1.24 * (Lnorm - 15) ** (1.0 / 3.0) - 1.8
    if (d0 < D0_MIN):
        d0 = D0_MIN
    d0_search = d0
    if (d0_search > 8):
        d0_search = 9
    if (d0_search < 4.5):
        d0_search = 4.5

    return D0_MIN, Lnorm, d0, d0_search

def superpose_set(query_set, target_set):
    sup = SVDSuperimposer()
    sup.set(target_set, query_set)
    sup.run()
    rot, trans = sup.get_rotran()
    return rot, trans

def sec_str(dis13, dis14, dis15, dis24, dis25, dis35, turn=True):
    s = 'C'
    delta = 2.1
    fabs = math.fabs
    if (fabs(dis15-6.37)<delta and fabs(dis14-5.18)<delta and
        fabs(dis25-5.18)<delta and fabs(dis13-5.45)<delta and
        fabs(dis24-5.45)<delta and fabs(dis35-5.45)<delta):
        s = 'H'
    
    delta = 1.42
    if (fabs(dis15-13  )<delta and fabs(dis14-10.4)<delta and
        fabs(dis25-10.4)<delta and fabs(dis13-6.1 )<delta and
        fabs(dis24-6.1 )<delta and fabs(dis35-6.1 )<delta):
        s = 'E'

    if turn == True:
        if (dis15 < 8.0):
            s = 'T'
    #if (dis15 < 8.0):
    #    s = 'T'

    return s

def CA_coord(pdb_id, chain_id='A'):
    aa3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
   'MSE':'M', 'PCA':'A',}

    CA_coord = []
    CA_seq = []
    pdb_name = 'input_pdb'
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_name, pdb_id)
    model = structure[0]
    if pdb_id.endswith('pdb') == False:
        chain = model
    else:
        chain = model
    for res in chain.get_residues():
        for atom in res:
            if atom.get_name() == 'CA':
                CA_coord.append(atom.get_coord())
                CA_seq.append(aa3to1[res.resname])

    CA_seq = ''.join(CA_seq)
    return CA_coord, CA_seq


def output_superposition(seq_x, seq_y, align_x, align_y):
    align_len = len(align_x)
    x_len = len(seq_x)
    y_len = len(seq_y)
    common_len = x_len + y_len - align_len
    out_x = []
    out_y = []
    result_x = []
    result_y = []
    x_point = 0
    y_point = 0
    for i in range(x_len):
        if i < align_x[x_point]:
            out_x.append('-')
        elif i == align_x[x_point]:
            out_x.append(align_y[x_point])
            x_point = min(x_point+1, align_len-1)
        elif i > align_x[x_point]:
            out_x.append('-')

    for j in range(y_len):
        if j < align_y[y_point]:
            out_y.append('-')
        elif j == align_y[y_point]:
            out_y.append(align_x[y_point])
            y_point = min(y_point+1, align_len-1)
        elif j > align_y[y_point]:
            out_y.append('-')

    xp = 0
    yp = 0
    #print(len(out_x), x_len)
    #print(len(out_y), y_len)

    for k in range(common_len):
        if (out_x[xp] == '-') and (out_y[yp] == '-'):
            result_x.append(seq_x[xp])
            result_y.append('-')
            xp = min(xp+1, len(out_x)-1)
        elif (out_x[xp] != '-') and (out_y[yp] == '-'):
            result_x.append('-')
            result_y.append(seq_y[yp])
            yp = min(yp+1, len(out_y)-1)
        elif (out_x[xp] == '-') and (out_y[yp] != '-'):
            result_x.append(seq_x[xp])
            result_y.append('-')
            xp = min(xp+1, len(out_x)-1)
        elif (out_x[xp] != '-') and (out_y[yp] != '-'):
            result_x.append(seq_x[xp])
            result_y.append(seq_y[yp])
            xp = min(xp+1, len(out_x)-1)
            yp = min(yp+1, len(out_y)-1)

    result_x = ''.join(result_x)
    result_y = ''.join(result_y)
    return result_x, result_y



def approx_sse(CA_coord, turn=True):
    CA_len = len(CA_coord)
    sse = []
    for i in range(CA_len):
        s = 'C'
        j1 = i - 2
        j2 = i - 1
        j3 = i
        j4 = i + 1
        j5 = i + 2
        if (j1 >= 0 and j5 < CA_len):
            d13 = np.linalg.norm(CA_coord[j1] - CA_coord[j3])
            d14 = np.linalg.norm(CA_coord[j1] - CA_coord[j4])
            d15 = np.linalg.norm(CA_coord[j1] - CA_coord[j5])
            d24 = np.linalg.norm(CA_coord[j2] - CA_coord[j4])
            d25 = np.linalg.norm(CA_coord[j2] - CA_coord[j5])
            d35 = np.linalg.norm(CA_coord[j3] - CA_coord[j5])
            s = sec_str(d13, d14, d15, d24, d25, d35, turn)

        sse.append(s)

    return sse

def TMmatrix(coord_1, coord_2):
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)
    letter_len_1 = coord_1.shape[0]
    letter_len_2 = coord_2.shape[0]
    d_min = min(letter_len_1, letter_len_2)
    d0 = 1.24 * (d_min - 15) ** (1.0 / 3.0) - 1.8
    #d0 = d0 + 0.8
    d02 = d0 * d0
    score_d8 = 1.5 * (d_min**0.3) + 3.5
    score_d82 = score_d8 * score_d8
    TMmat = np.zeros((letter_len_1, letter_len_2))
    #print("d0: {}".format(d0))
    #print("score_d82: {}".format(score_d8))

    for i in range(letter_len_1):
        for j in range(letter_len_2):
            s = coord_1[i,:] - coord_2[j,:]
            s = (s * s).sum()
            if (s > score_d82):
                s = 0
            else:
                s = 1.0 / (1.0 + s / d02)
            TMmat[i, j] = TMmat[i, j] + s

    return TMmat

def SPmatrix(coord_1, coord_2, Effo=1.3):
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)
    letter_len_1 = coord_1.shape[0]
    letter_len_2 = coord_2.shape[0]
    d_min = min(letter_len_1, letter_len_2)
    #d0 = 1.24 * (d_min - 15) ** (1.0 / 3.0) - 1.8
    d0 = 4
    d02 = d0 * d0
    SPMat = np.zeros((letter_len_1, letter_len_2))
    for i in range(letter_len_1):
        for j in range(letter_len_2):
            s = ((coord_1[i,:] - coord_2[j,:]) * (coord_1[i,:] - coord_2[j,:])).sum()
            if s < 4 * Effo * d02:
                SPMat[i,j] = (1.0 / (1.0 + s * Effo / d02) - 1.0 / (1.0 + 4 * Effo))
            else:
                SPMat[i,j] = 0
            if SPMat[i, j] < 0:
                SPMat[i, j] = 0
    #print(SPMat)
    return SPMat
            

def mican_matrix(coord_1, coord_2, sse_1, sse_2, d_R):
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)
    letter_len_1 = coord_1.shape[0]
    letter_len_2 = coord_2.shape[0]
    d_min = min(letter_len_1, letter_len_2)
    d0 = 1.24 * (d_min - 15) ** (1.0 / 3.0) - 1.8
    d02 = d0 * d0
    d_R_2 = d_R * d_R
    mican_mat = np.zeros((letter_len_1, letter_len_2))
    
    for i in range(letter_len_1):
        for j in range(letter_len_2):
            s = ((coord_1[i,:] - coord_2[j,:]) * (coord_1[i,:] - coord_2[j,:])).sum()
            s = 1.0 / (1.0 + s / d02) * ((float(sse_1[i] == sse_2[j]) + 1) / (2)) * float(s < d_R_2)
            mican_mat[i, j] = mican_mat[i, j] + s

    return mican_mat

def ft_matrix(coord_1, coord_2, sse_1, sse_2, d_R):
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)
    letter_len_1 = coord_1.shape[0]
    letter_len_2 = coord_2.shape[0]
    d_min = min(letter_len_1, letter_len_2)
    d0 = 1.24 * (d_min - 15) ** (1.0 / 3.0) - 1.8
    d02 = d0 * d0
    d_R_2 = d_R * d_R
    ft_mat = np.zeros((letter_len_1, letter_len_2))
    delta = 0.0

    for i in range(letter_len_1):
        for j in range(letter_len_2):
            s = ((coord_1[i,:] - coord_2[j,:]) * (coord_1[i,:] - coord_2[j,:])).sum()
            if sse_1[i] == sse_2[j]:
                if sse_1[i] == 'C':
                    delta = 1.0
                else:
                    delta = 2.0
            elif (sse_1[i] == 'C' and sse_2[j] == 'H') or (sse_1[i] == 'H' and sse_2[j] == 'C'):
                delta = 0.0
            elif (sse_1[i] == 'C' and sse_2[j] == 'E') or (sse_1[i] == 'E' and sse_2[j] == 'C'):
                delta = 0.5
            elif (sse_1[i] == 'H' and sse_2[j] == 'E') or (sse_1[i] == 'E' and sse_2[j] == 'H'):
                delta = -0.5
            s = delta / (1.0 + s / d02) * float(s < d_R_2)
            ft_mat[i, j] = ft_mat[i, j] + s
    return ft_mat


def IndxFromPro(Pmat):
    RowIndx = np.argmax(Pmat, axis=1)
    ColIndx = np.argmax(Pmat, axis=0)
    row_indx_max = len(RowIndx) - 1
    col_indx_max = len(ColIndx) - 1

    row_align_indx = []
    col_align_indx = []
    #print("row_indx_max: {}".format(row_indx_max))
    for i in range(Pmat.shape[0]-1):
        if (RowIndx[i] != col_indx_max):
            row_align_indx.append(RowIndx[i])
        else:
            row_align_indx.append(-1)

    for j in range(Pmat.shape[1]-1):
        if (ColIndx[j] != row_indx_max):
            col_align_indx.append(ColIndx[j])
        else:
            col_align_indx.append(-1)

    return row_align_indx, col_align_indx

def OT_align_refine(coord_1, coord_2, gap, epsilon, lmbd, T, mat):
    TMmat = mat
    Pmat = lsap2ot(TMmat, gap, epsilon, lmbd, T)
    row_align_indx, col_align_indx = IndxFromPro(Pmat)
    if True:
        row_indx_unique = list(set(row_align_indx))
        for i in range(len(row_indx_unique)):
            if row_indx_unique[i] != -1:
                unique_indx = (row_align_indx == row_indx_unique[i]).astype(float)
                if (np.sum(unique_indx)) > 1:
                    indx_list = [j for j,val in enumerate(row_align_indx) if val==row_indx_unique[i]]
                    score_list = [TMmat[indx_list[k], row_indx_unique[i]] for k in range(len(indx_list))]
                    max_indx = np.argmax(score_list)
                    indx_list.pop(max_indx)
                    for m, mm in enumerate(indx_list):
                        row_align_indx[mm] = -1

        col_indx_unique = list(set(col_align_indx))
        for i in range(len(col_indx_unique)):
            if col_indx_unique[i] != -1:
                unique_indx = (col_align_indx == col_indx_unique[i]).astype(float)
                if (np.sum(unique_indx)) > 1:
                    indx_list = [j for j,val in enumerate(col_align_indx) if val==col_indx_unique[i]]
                    score_list = [TMmat[col_indx_unique[i], indx_list[k]] for k in range(len(indx_list))]
                    max_indx = np.argmax(score_list)
                    indx_list.pop(max_indx)
                    for m, mm in enumerate(indx_list):
                        col_align_indx[mm] = -1

    return row_align_indx, col_align_indx



def OT_iter(coord_1, coord_2, g1, g2, gap, epsilon, lmbd, T, iteration_max, local_d0_search, d0, score_sum_method, score_d8):
    gap_open = [-0.6, 0]
    d02 = d0 * d0
    simplify_step = 40
    opt_align_1 = []
    opt_align_2 = []
    tmscore_max = -1

    for g in range(g1, g2):
        for iteration in range(iteration_max):
            TMmat = TMmatrix(coord_1, coord_2)
            Pmat = lsap2ot(TMmat, gap, epsilon, lmbd, T)

    return Pmat




def epLSAP_Main(pdb1, pdb2, fast_opt=False, mode='TM', lmbd=50):
    SO_th = 3.5
    coord1, seq1 = CA_coord(pdb1)
    coord2, seq2 = CA_coord(pdb2)
    coord1_array = np.array(coord1)
    coord2_array = np.array(coord2)
    sse1 = approx_sse(coord1)
    sse2 = approx_sse(coord2)
    SO_th_2 = SO_th * SO_th

    D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0 = parameter_set4search(len(coord1), len(coord2))
    local_d0_search = d0_search
    simplify_step = 1#40
    score_sum_method = 8

    if mode == 'TM':
        tmtool_score = tm_align(coord1_array, coord2_array, seq1, seq2)
        coord1_array_t1 = np.dot(coord1_array, (tmtool_score.u).T) + tmtool_score.t
        u0 = (tmtool_score.u).T
        t0 = tmtool_score.t
    else:
        mic = mican()
        aln = mic.align(pdb1, pdb2)
        mic_rot, mic_trans = (aln.translation_rot).T, aln.translation_vec
        coord1_array_t1 = np.dot(coord1_array, mic_rot) + mic_trans
        u0 = mic_rot
        t0 = mic_trans
    TMmat = SPmatrix(coord1_array_t1, coord2_array)
    gap = -0.5 # SPacore 0.5, TMscore 0.6
    epsilon = 0.05 #0.05
    T = 5000 #50
    OT_Nali = 0
    OT_RMSD = 10000
    OT_TMscore = 0
    OT_SO = -1
    SO_m = min(coord1_array.shape[0], coord2_array.shape[0])
    minlen = SO_m
    coord1_array_list = [coord1_array_t1]#, coord1_array_t2]
    row_len, col_len = TMmat.shape[0], TMmat.shape[1]
    D0_MIN, Lnorm, d0, d0_search = parameter_set4final(minlen)
    if row_len > col_len:
        which_maj = 'ROW'
    else:
        which_maj = 'COL'
    print(coord1_array.shape[0], coord2_array.shape[0])
    output_fwdmap = []
    for idx in range(len(coord1_array_list)):
        coord1_array_t = coord1_array_list[idx]
        row_align_indx, col_align_indx = OT_align_refine(coord1_array_t, coord2_array, gap, epsilon, lmbd, T, mat=TMmat)
        print(row_align_indx)
        print(col_align_indx)
        #np.savetxt('row_info.txt', list(map(int, row_align_indx)))
        #np.savetxt('col_info.txt', list(map(int, col_align_indx)))
        fwdmap = row_align_indx
        invmap = [-1 for x in range(len(row_align_indx))]
        for i in range(len(col_align_indx)):
            if col_align_indx[i] != -1:
                invmap[col_align_indx[i]] = i
        print("row align info")
        row_align_len = 0
        row_align_RMSD = 0
        row_align_SO = 0
        xtm = []
        ytm = []
        xt = []
        tm_row_len = 0
        row_dst = []
        for i in range(len(row_align_indx)):
            if row_align_indx[i] != -1:
                dst2 = ((coord1_array_t[i,:] - coord2_array[row_align_indx[i],:]) * (coord1_array_t[i,:] - coord2_array[row_align_indx[i],:])).sum()
                row_dst.append(np.sqrt(dst2))
                row_align_len += 1
                row_align_RMSD += dst2 #((coord1_array_t[i,:] - coord2_array[row_align_indx[i],:]) * (coord1_array_t[i,:] - coord2_array[row_align_indx[i],:])).sum()
                if dst2 < SO_th_2:#((coord1_array_t[i,:] - coord2_array[row_align_indx[i],:]) * (coord1_array_t[i,:] - coord2_array[row_align_indx[i],:])).sum() < SO_th_2:
                    row_align_SO += 1
                if dst2 < score_d8:
                    xtm.append(coord1_array[i,:])
                    xt.append(coord1_array_t[i,:])
                    ytm.append(coord2_array[row_align_indx[i],:])
                    tm_row_len += 1
        row_align_RMSD = (row_align_RMSD / (row_align_len + 1e-7)) ** 0.5
        #np.savetxt('row_dst.txt', row_dst)
        xtm = np.array(xtm)
        ytm = np.array(ytm)
        xt = np.array(xt)
        print(local_d0_search, score_d8, d0)
        TM_row = TMalign_tools.tmscore(xtm, ytm, xt, minlen, tm_row_len, t0, u0.T, simplify_step, score_sum_method, local_d0_search, Lnorm, score_d8, d0)
        print("{} OT row alignment Nali: {}, RMSD: {}, TMscore: {}, SO: {}".format(idx, row_align_len, row_align_RMSD, TM_row, row_align_SO/SO_m))
        if row_align_RMSD < OT_RMSD: #row_align_SO/SO_m > OT_SO:
            OT_Nali = row_align_len
            OT_RMSD = row_align_RMSD
            OT_TMscore = TM_row
            OT_SO = row_align_SO/SO_m
            which_SO = 'ROW'
            output_fwdmap = fwdmap

        #print(col_align_indx)
        #print(len(col_align_indx))
        col_align_len = 0
        col_align_RMSD = 0
        col_align_SO = 0
        xtm = []
        ytm = []
        xt = []
        tm_col_len = 0
        col_dst = []
        for i in range(len(col_align_indx)):
            if col_align_indx[i] != -1:
                col_align_len += 1
                dst2 = ((coord1_array_t[col_align_indx[i],:] - coord2_array[i,:]) * (coord1_array_t[col_align_indx[i],:] - coord2_array[i,:])).sum()
                col_dst.append(np.sqrt(dst2))
                col_align_RMSD += dst2 #((coord1_array_t[col_align_indx[i],:] - coord2_array[i,:]) * (coord1_array_t[col_align_indx[i],:] - coord2_array[i,:])).sum()
                if dst2 < SO_th_2:#((coord1_array_t[col_align_indx[i],:] - coord2_array[i,:]) * (coord1_array_t[col_align_indx[i],:] - coord2_array[i,:])).sum() < SO_th_2:
                    col_align_SO += 1
                dst2 = ((coord1_array_t[col_align_indx[i],:] - coord2_array[i,:]) * (coord1_array_t[col_align_indx[i],:] - coord2_array[i,:])).sum()
                if dst2 < score_d8:
                    xtm.append(coord1_array[col_align_indx[i],:])
                    xt.append(coord1_array_t[col_align_indx[i],:])
                    ytm.append(coord2_array[i,:])
                    tm_col_len += 1
        col_align_RMSD = (col_align_RMSD / (col_align_len + 1e-7)) ** 0.5
        #np.savetxt('col_dst.txt', row_dst)
        xtm = np.array(xtm)
        ytm = np.array(ytm)
        xt = np.array(xt)
        TM_col = TMalign_tools.tmscore(xtm, ytm, xt, minlen, tm_col_len, t0, u0.T, simplify_step, score_sum_method, local_d0_search, Lnorm, score_d8, d0)
        print("{} OT col alignment Nali: {}, RMSD: {}, TMscore: {}, SO: {}".format(idx, col_align_len, col_align_RMSD, TM_col, col_align_SO/SO_m))
        if col_align_RMSD < OT_RMSD:#col_align_SO/SO_m > OT_SO:
            OT_Nali = col_align_len
            OT_RMSD = col_align_RMSD
            OT_TMscore = TM_col
            OT_SO = col_align_SO/SO_m
            which_SO = 'COL'
            output_fwdmap = invmap

    print(TM_row, TM_col)

    print("Final OT alignment Nali: {}, RMSD: {}, TMscore: {}, SO: {}\n".format(OT_Nali, OT_RMSD, OT_TMscore, OT_SO))
    if fast_opt == False:
        return OT_Nali, OT_RMSD, OT_TMscore, OT_SO
    else:
        return OT_Nali, OT_RMSD, OT_TMscore, OT_SO, output_fwdmap




    

    


            

    
    

    

