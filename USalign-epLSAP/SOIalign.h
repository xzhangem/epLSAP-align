#ifndef SOIalign_h
#define SOIalign_h 1

#include "TMalign.h"

void print_invmap(int *invmap, const int ylen)
{
    int i,j;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i>=0) cout<<" ("<<i<<","<<j<<")";
    }
    cout<<endl;
}

void assign_sec_bond(int **secx_bond, const char *secx, const int xlen)
{
    int i,j;
    int starti=-1;
    int endi=-1;
    char ss;
    char prev_ss=0;
    for (i=0; i<xlen; i++)
    {
        ss=secx[i];
        secx_bond[i][0]=secx_bond[i][1]=-1;
        if (ss!=prev_ss && !(ss=='C' && prev_ss=='T') 
                        && !(ss=='T' && prev_ss=='C'))
        {
            if (starti>=0) // previous SSE end
            {
                endi=i;
                for (j=starti;j<endi;j++)
                {
                    secx_bond[j][0]=starti;
                    secx_bond[j][1]=endi;
                }
            }
            if (ss=='H' || ss=='E' || ss=='<' || ss=='>') starti=i;
            else starti=-1;
        }
        prev_ss=secx[i];
    }
    if (starti>=0) // previous SSE end
    {
        endi=i;
        for (j=starti;j<endi;j++)
        {
            secx_bond[j][0]=starti;
            secx_bond[j][1]=endi;
        }
    }
    for (i=0;i<xlen;i++) if (secx_bond[i][1]-secx_bond[i][0]==1)
        secx_bond[i][0]=secx_bond[i][1]=-1;
}

void getCloseK(double **xa, const int xlen, const int closeK_opt, double **xk)
{
    double **score;
    NewArray(&score, xlen+1, xlen+1);
    vector<pair<double,int> > close_idx_vec(xlen, make_pair(0,0));
    int i,j,k;
    for (i=0;i<xlen;i++)
    {
        score[i+1][i+1]=0;
        for (j=i+1;j<xlen;j++) score[j+1][i+1]=score[i+1][j+1]=dist(xa[i], xa[j]);
    }
    for (i=0;i<xlen;i++)
    {
        for (j=0;j<xlen;j++)
        {
            close_idx_vec[j].first=score[i+1][j+1];
            close_idx_vec[j].second=j;
        }
        sort(close_idx_vec.begin(), close_idx_vec.end());
        for (k=0;k<closeK_opt;k++)
        {
            j=close_idx_vec[k % xlen].second;
            xk[i*closeK_opt+k][0]=xa[j][0];
            xk[i*closeK_opt+k][1]=xa[j][1];
            xk[i*closeK_opt+k][2]=xa[j][2];
        }
    }

    /* clean up */
    vector<pair<double,int> >().swap(close_idx_vec);
    DeleteArray(&score, xlen+1);
}

/* check if pairing i to j conform to sequantiality within the SSE */
inline bool sec2sq(const int i, const int j,
    int **secx_bond, int **secy_bond, int *fwdmap, int *invmap)
{
    if (i<0 || j<0) return true;
    int ii,jj;
    if (secx_bond[i][0]>=0)
    {
        for (ii=secx_bond[i][0];ii<secx_bond[i][1];ii++)
        {
            jj=fwdmap[ii];
            if (jj>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    if (secy_bond[j][0]>=0)
    {
        for (jj=secy_bond[j][0];jj<secy_bond[j][1];jj++)
        {
            ii=invmap[jj];
            if (ii>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    return true;
}



void ot_gl_soi_egs(double **TMave_mat, const int chain1_num, const int chain2_num, int *fwdmap, int *invmap, int **secx_bond, int **secy_bond, const int mm_opt){
        const int ITER = 35;
        const double EPSILON = 0.05;
        const double LMBD = -15.0;
        const double ERR = 0.5;

        double total_score = 0.0;

        /* initialize parameters */
	for (int i = 0; i < chain1_num; i++) fwdmap[i] = -1;
        for (int j = 0; j < chain2_num; j++) invmap[j] = -1;

        /* initialize cost function from TM_mat */
        double C = -100.0;
        for (int i = 0; i < chain1_num; ++i){
                for (int j = 0; j < chain2_num; ++j){
                        if (TMave_mat[i][j] > C) C = TMave_mat[i][j];
                }
        }
        C += EPSILON;
        double C_2 = C * 2;
        vector<vector<double> > otMat, expMat, stochasticMat, stochasticMat_T;
        build_2D_array(otMat, chain1_num, chain2_num);
        build_2D_array(expMat, chain1_num, chain2_num);
        build_2D_array(stochasticMat, chain1_num, chain2_num);
        build_2D_array(stochasticMat_T, chain2_num, chain1_num);

        for (int i = 0; i < chain1_num; i++){
                for (int j = 0; j < chain2_num; j++){
			otMat[i][j] = C_2 - TMave_mat[i][j];
			expMat[i][j] = exp(LMBD * otMat[i][j]);
		}
	}
	vector<double > y_vec(chain2_num, 1);
        vector<double > yp(chain2_num, 1);
        bool conv = false;
        int tag = 0;
        vector<double > x_vec(chain1_num, 1);
        vector<double > xp(chain1_num, 1);
	while (conv == false && tag < ITER){
                for (int i = 0; i < chain1_num; ++i){
                        double xp_temp = 0.0;
                        for (int j = 0; j < chain2_num; ++j){
                                xp_temp += expMat[i][j] * y_vec[j];
                        }
                        xp_temp = 1.0 / xp_temp;
                        xp[i] = xp_temp;
                }
                double err_x = 0.0;
                double err_y = 0.0;
                if (tag > 0){
                        for (int i = 0; i < chain1_num; ++i){
                                err_x += (xp[i] / x_vec[i] - 1) * (xp[i] / x_vec[i] - 1);
                        }
                }
                //else{
                for (int i = 0; i < chain1_num; ++i){
                        x_vec[i] = xp[i];
                }
                for (int j = 0; j < chain2_num; ++j){
                        double yp_temp = 0.0;
                        for (int i = 0; i < chain1_num; ++i){
                                yp_temp += expMat[i][j] * x_vec[i];
                        }
                        yp_temp = 1.0 / yp_temp;
                        yp[j] = yp_temp;
                }
                for (int j = 0; j < chain2_num; ++j){
                        err_y += (yp[j] / y_vec[j] - 1) * (yp[j] / y_vec[j] - 1);
                }
                for (int j = 0; j < chain2_num; ++j){
                        y_vec[j] = yp[j];
                }
                //}
                if (tag > 0 && err_x < ERR && err_y < ERR){
                        conv = true;
                }
                tag += 1;
        }
	double mat_temp;
        for (int i = 0; i < chain1_num; ++i){
                for (int j = 0; j < chain2_num; ++j){
                        mat_temp = x_vec[i] * expMat[i][j] * y_vec[j];
                        stochasticMat[i][j] = mat_temp;
                        stochasticMat_T[j][i] = mat_temp;
                }
        }

        /* roughly extract the positions of maximum value of each row and col */
        vector<int> rowIndx(chain1_num);
        vector<int> colIndx(chain2_num);
        int indx_temp;
        for (int i = 0; i < chain1_num; i++){
                indx_temp = maxloc(chain2_num, stochasticMat[i]);
                rowIndx[i] = indx_temp;
        }
        for (int j = 0; j < chain2_num; j++){
                indx_temp = maxloc(chain1_num, stochasticMat_T[j]);
                colIndx[j] = indx_temp;
        }
        /* delete the repeated elements */
        set<int> rowSet(rowIndx.begin(), rowIndx.end());
        set<int> colSet(colIndx.begin(), colIndx.end());
        int rowSetLen = rowSet.size();
        int colSetlen = colSet.size();
        vector<int> indxTemp;
        vector<double> valueTemp;
	/* Align according to row information */
        for (auto it = rowSet.begin(); it != rowSet.end(); ++it){
                if (*it != -1){
                        for (int le = 0; le < chain1_num; le++){
                                if (rowIndx[le] == *it){
                                        indxTemp.push_back(le);
                                        valueTemp.push_back(TMave_mat[le][rowIndx[le]]);
                                }
                        }
                        if (valueTemp.size() > 1){
                                indx_temp = maxloc(valueTemp.size(), valueTemp);
                                for (int li = 0; li < indxTemp.size(); ++li){
                                        if (indxTemp[li] != indx_temp) rowIndx[indxTemp[li]] = -1;
                                }
                        }
                        if (indxTemp.empty() == false) indxTemp.clear();
                        if (valueTemp.empty() == false) valueTemp.clear();
                }
        }
	/* Align according to col information */
        for (auto jt = colSet.begin(); jt != colSet.end(); ++jt){
                if (*jt != -1){
                        for (int le = 0; le < chain2_num; le++){
                                if (colIndx[le] == *jt){
                                        indxTemp.push_back(le);
                                        valueTemp.push_back(TMave_mat[colIndx[le]][le]);
                                }
                        }
                        if (valueTemp.size() > 1){
                                indx_temp = maxloc(valueTemp.size(), valueTemp);
                                for (int lj = 0; lj < indxTemp.size(); lj++){
                                        if (indxTemp[lj] != indx_temp) colIndx[indxTemp[lj]] = -1;
                                }
                        }
                        if (indxTemp.empty() == false) indxTemp.clear();
                        if (valueTemp.empty() == false) valueTemp.clear();
                }
        }
	for (int i = 0; i < chain1_num; i++){
                fwdmap[i] = rowIndx[i];
        }
        for (int j = 0; j < chain2_num; j++){
                invmap[j] = colIndx[j];
        }
}



void ot_gp_soi_egs(double **TMave_mat, const int chain1_num, const int chain2_num, int *fwdmap, int *invmap, int **secx_bond, int **secy_bond, const int mm_opt){
        const int ITER = 5000;
        const double EPSILON = 0.05;
        const double GAP = -0.50;
        const double LMBD = -100.0;
        const double ERR = 0.05;

        double total_score = 0.0;

        /* initialize parameters */
	for (int i = 0; i < chain1_num; i++) fwdmap[i] = -1;
	for (int j = 0; j < chain2_num; j++) invmap[j] = -1;

        /* initialize cost function from TM_mat */
        double C = -100.0;
        for (int i = 0; i < chain1_num; ++i){
                for (int j = 0; j < chain2_num; ++j){
                        if (TMave_mat[i][j] > C) C = TMave_mat[i][j];
                }
        }
        C += EPSILON;
        double C_2 = C * 2;
        vector<vector<double> > otMat, expMat, stochasticMat, stochasticMat_T;
        build_2D_array(otMat, chain1_num+1, chain2_num+1);
        build_2D_array(expMat, chain1_num+1, chain2_num+1);
        build_2D_array(stochasticMat, chain1_num+1, chain2_num+1);
        build_2D_array(stochasticMat_T, chain2_num+1, chain1_num+1);

        for (int i = 0; i < chain1_num+1; i++){
                for (int j = 0; j < chain2_num+1; j++){
                        if (i < chain1_num && j < chain2_num){
                                otMat[i][j] = C_2 - TMave_mat[i][j];
                                expMat[i][j] = exp(LMBD * otMat[i][j]);
                        } else{
                                otMat[i][j] = C - GAP;
                                expMat[i][j] = exp(LMBD * otMat[i][j]);
                        }
                }
        }
        otMat[chain1_num][chain2_num] = 0.0;
        expMat[chain1_num][chain2_num] = 1.0;

        /* Obtain the stochastic matrix by Sinkhorn algorithm */
        vector<double > y_vec(chain2_num+1, 1);
        vector<double > yp(chain2_num+1, 1);
	bool conv = false;
        int tag = 0;
        vector<double > x_vec(chain1_num+1, 1);
        vector<double > xp(chain1_num+1);
	while (conv == false && tag < ITER){
                for (int i = 0; i < chain1_num+1; ++i){
                        double xp_temp = 0.0;
                        for (int j = 0; j < chain2_num+1; ++j){
                                xp_temp += expMat[i][j] * y_vec[j];
                        }
                        xp_temp = 1.0 / xp_temp;
                        xp[i] = xp_temp;
                }
                xp[chain1_num] = 1.0;
                double err_x = 0.0;
                double err_y = 0.0;
                if (tag > 0){
                        for (int i = 0; i < chain1_num+1; ++i){
                                err_x += (xp[i] / x_vec[i] - 1) * (xp[i] / x_vec[i] - 1);
                        }
                }
                //else{
                for (int i = 0; i < chain1_num+1; ++i){
                        x_vec[i] = xp[i];
                }
                for (int j = 0; j < chain2_num+1; ++j){
                        double yp_temp = 0.0;
                        for (int i = 0; i < chain1_num+1; ++i){
                                yp_temp += expMat[i][j] * x_vec[i];
                        }
                        yp_temp = 1.0 / yp_temp;
                        yp[j] = yp_temp;
                }
                yp[chain2_num] = 1.0;
                for (int j = 0; j < chain2_num+1; ++j){
                        err_y += (yp[j] / y_vec[j] - 1) * (yp[j] / y_vec[j] - 1);
                }
                for (int j = 0; j < chain2_num+1; ++j){
                        y_vec[j] = yp[j];
                }
                //}
                if (tag > 0 && err_x < ERR && err_y < ERR){
                        conv = true;
                }
                tag += 1;
        }
	double mat_temp;
        for (int i = 0; i < chain1_num+1; ++i){
                for (int j = 0; j < chain2_num+1; ++j){
                        mat_temp = x_vec[i] * expMat[i][j] * y_vec[j];
                        stochasticMat[i][j] = mat_temp;
                        stochasticMat_T[j][i] = mat_temp;
                }
        }

        /* roughly extract the positions of maximum value of each row and col */
        vector<int> rowIndx(chain1_num);
        vector<int> colIndx(chain2_num);
        int indx_temp;
        for (int i = 0; i < chain1_num; i++){
                indx_temp = maxloc(chain2_num+1, stochasticMat[i]);
                if (indx_temp == chain2_num) indx_temp = -1;
                rowIndx[i] = indx_temp;
        }
        for (int j = 0; j < chain2_num; j++){
                indx_temp = maxloc(chain1_num+1, stochasticMat_T[j]);
                if (indx_temp == chain1_num) indx_temp = -1;
                colIndx[j] = indx_temp;
        }
        /* delete the repeated elements */
        set<int> rowSet(rowIndx.begin(), rowIndx.end());
        set<int> colSet(colIndx.begin(), colIndx.end());
        int rowSetLen = rowSet.size();
        int colSetlen = colSet.size();
        vector<int> indxTemp;
        vector<double> valueTemp;
	/* Align according to row information */
	for (auto it = rowSet.begin(); it != rowSet.end(); ++it){
                if (*it != -1){
                        for (int le = 0; le < chain1_num; le++){
                                if (rowIndx[le] == *it){
                                        indxTemp.push_back(le);
                                        valueTemp.push_back(TMave_mat[le][rowIndx[le]]);
                                }
                        }
                        if (valueTemp.size() > 1){
                                indx_temp = maxloc(valueTemp.size(), valueTemp);
                                for (int li = 0; li < indxTemp.size(); ++li){
                                        if (indxTemp[li] != indx_temp) rowIndx[indxTemp[li]] = -1;
                                }
                        }
                        if (indxTemp.empty() == false) indxTemp.clear();
                        if (valueTemp.empty() == false) valueTemp.clear();
                }
        }
	/* Align according to col information */
	for (auto jt = colSet.begin(); jt != colSet.end(); ++jt){
                if (*jt != -1){
                        for (int le = 0; le < chain2_num; le++){
                                if (colIndx[le] == *jt){
                                        indxTemp.push_back(le);
                                        valueTemp.push_back(TMave_mat[colIndx[le]][le]);
                                }
                        }
                        if (valueTemp.size() > 1){
                                indx_temp = maxloc(valueTemp.size(), valueTemp);
                                for (int lj = 0; lj < indxTemp.size(); lj++){
                                        if (indxTemp[lj] != indx_temp) colIndx[indxTemp[lj]] = -1;
                                }
                        }
                        if (indxTemp.empty() == false) indxTemp.clear();
                        if (valueTemp.empty() == false) valueTemp.clear();
                }
        }
	for (int i = 0; i < chain1_num; i++){
		fwdmap[i] = rowIndx[i];
	}
	for (int j = 0; j < chain2_num; j++){
		invmap[j] = colIndx[j];
	}
}


void soi_egs(double **score, const int xlen, const int ylen, int *invmap,
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    int i,j;
    int *fwdmap=new int[xlen]; // j=fwdmap[i];
    for (i=0; i<xlen; i++) fwdmap[i]=-1;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i>=0) fwdmap[i]=j;
    }

    /* stage 1 - make initial assignment, starting from the highest score pair */
    double max_score;
    int maxi,maxj;
    while(1)
    {
        max_score=0;
        maxi=maxj=-1;
        for (i=0;i<xlen;i++)
        {
            if (fwdmap[i]>=0) continue;
            for (j=0;j<ylen;j++)
            {
                if (invmap[j]>=0 || score[i+1][j+1]<=max_score) continue;
                if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond,
                    fwdmap,invmap)) continue;
                maxi=i;
                maxj=j;
                max_score=score[i+1][j+1];
            }
        }
        if (maxi<0) break; // no assignment;
        invmap[maxj]=maxi;
        fwdmap[maxi]=maxj;
    }

    double total_score=0;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i>=0) total_score+=score[i+1][j+1];
    }

    /* stage 2 - swap assignment until total score cannot be improved */
    int iter;
    int oldi,oldj;
    double delta_score;
    for (iter=0; iter<getmin(xlen,ylen)*5; iter++)
    {
        //cout<<"total_score="<<total_score<<".iter="<<iter<<endl;
        //print_invmap(invmap,ylen);
        delta_score=-1;
        for (i=0;i<xlen;i++)
        {
            oldj=fwdmap[i];
            for (j=0;j<ylen;j++)
            {
                oldi=invmap[j];
                if (score[i+1][j+1]<=0 || oldi==i) continue;
                if (mm_opt==6 && (!sec2sq(i,j,secx_bond,secy_bond,fwdmap,invmap) ||
                            !sec2sq(oldi,oldj,secx_bond,secy_bond,fwdmap,invmap)))
                    continue;
                delta_score=score[i+1][j+1];
                if (oldi>=0 && oldj>=0) delta_score+=score[oldi+1][oldj+1];
                if (oldi>=0) delta_score-=score[oldi+1][j+1];
                if (oldj>=0) delta_score-=score[i+1][oldj+1];

                if (delta_score>0) // successful swap
                {
                    fwdmap[i]=j;
                    if (oldi>=0) fwdmap[oldi]=oldj;
                    invmap[j]=i;
                    if (oldj>=0) invmap[oldj]=oldi;
                    total_score+=delta_score;
                    break;
                }
            }
        }
        if (delta_score<=0) break; // cannot make further swap
    }

    /* clean up */
    delete[]fwdmap;
}

/* entry function for se
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * */
int soi_se_main(
    double **xa, double **ya, const char *seqx, const char *seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen, 
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool a_opt, const int u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, int *invmap, double *dist_list,
    int **secx_bond, int **secy_bond, const int mm_opt, double &L_SO)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double **score;       // score for aligning a residue pair
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  

    int *m1=NULL;
    int *m2=NULL;
    int i,j;
    double d;
    if (outfmt_opt<2)
    {
        m1=new int[xlen]; //alignd index in x
        m2=new int[ylen]; //alignd index in y
    }

    /***********************/
    /* allocate memory     */
    /***********************/
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path,  xlen+1, ylen+1);
    NewArray(&val,   xlen+1, ylen+1);
    //int *invmap          = new int[ylen+1];

    /* set d0 */
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
        score_d8, d0, d0_search, dcu0); // set score_d8
    parameter_set4final(xlen, D0_MIN, Lnorm,
        d0B, d0_search, mol_type); // set d0B
    parameter_set4final(ylen, D0_MIN, Lnorm,
        d0A, d0_search, mol_type); // set d0A
    if (a_opt)
        parameter_set4final((xlen+ylen)*0.5, D0_MIN, Lnorm,
            d0a, d0_search, mol_type); // set d0a
    if (u_opt)
    {
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0u, d0_search, mol_type); // set d0u
        if (u_opt==2)
        {
            parameter_set4search(Lnorm_ass, Lnorm_ass, D0_MIN, Lnorm,
                score_d8, d0, d0_search, dcu0); // set score_d8
        }
    }

    /* perform alignment */
    for(j=0; j<ylen; j++) invmap[j]=-1;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for(i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xa[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        dist_list[j]=-1;
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(&xa[i][0], &ya[j][0]));
            dist_list[j]=d;
            if (score[i+1][j+1]>0)
            {
                if (outfmt_opt<2)
                {
                    m1[k]=i;
                    m2[k]=j;
                }
                k++;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);

    if (outfmt_opt>=2)
    {
        DeleteArray(&score, xlen+1);
        return 0;
    }

    /* extract aligned sequence */
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    int *fwdmap = new int [xlen+1];
    for (i=0;i<xlen;i++) fwdmap[i]=-1;
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap[j];
        if (i<0) continue;
        d=sqrt(dist(xa[i], ya[j]));
        if (d<d0_out) seqM[j]=':';
	if (d<3.5) L_SO += 1;
        else seqM[j]='.';
        fwdmap[i]=j;
        seqxA[j]=seqx[i];
        Liden+=(seqxA[k]==seqyA[k]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }

    /* free memory */
    delete [] fwdmap;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);
    return 0; // zero for no exception
}

inline void SOI_SPsuper2score(double **xt, double **ya, const int xlen,
		const int ylen, double **sp_score)
{
	int i,j;
	double d2;
	int d_min = xlen;
	if (xlen > ylen) d_min = ylen;
	double Effo = 1.3;
	double d0 = 4.0;
	double d02 = 16.0;
	double th = 4 * Effo * d02;
	for (i=0;i<xlen;i++){
		for (j=0;j<ylen;j++){
			d2 = dist(xt[i], ya[j]);
			if (d2 > th) sp_score[i+1][j+1] = 0.0;
			else sp_score[i+1][j+1] = (1.0 / (1.0 + d2 * Effo / d02) - 1.0 / (1.0 + 4 * Effo));
			if (sp_score[i+1][j+1] < 0) sp_score[i+1][j+1] = 0;
		}
	}
}


inline void SOI_super2rmsd(double **xt, double **ya, const int xlen, const int ylen, double **rmsd_mat){
	int i, j;
	for (int i = 0; i < xlen; i++){
		for (int j = 0; j < ylen; j++){
			rmsd_mat[i][j] = dist(xt[i], ya[j]);
		}
	}
}

inline void SOI_super2score(double **xt, double **ya, const int xlen,
    const int ylen, double **score, double d0, double score_d8)
{
    int i,j;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for (i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xt[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1] = 0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
}

inline void trans_OTscore(double **tmscore, double **otscore, const int xlen, const int ylen){
	for (int i = 0; i < xlen; i++){
		for (int j = 0; j < ylen; j++){
			otscore[i][j] = tmscore[i+1][j+1];
		}
	}
}

double OT_TMiter(double **r1, double **r2, double **xtm, double **ytm, 
	double **xt, double **tmscore, bool **path, double **val, double **xa, double **ya, 
	int xlen, int ylen, double t[3], double u[3][3], int *invmap0, int *fwdmap0,
	int iteration_max, double local_d0_search,
	double Lnorm, double d0, double score_d8,
	int **secx_bond, int **secy_bond, const int mm_opt, double **rmsd_matrix, double & TM1, double & TM2, const bool init_invmap=false){

	double **xtm_temp, **ytm_temp, **otscore, **xtm_score_d8, **ytm_score_d8;
	int minlen = min(xlen, ylen);
	NewArray(&xtm_temp, minlen, 3);
	NewArray(&ytm_temp, minlen, 3);
	NewArray(&xtm_score_d8, minlen, 3);
	NewArray(&ytm_score_d8, minlen, 3);
	NewArray(&otscore, xlen, ylen);

	double rmsd;
	int *fwdmap = new int[xlen];
	int *invmap = new int[ylen];
	
	int iteration, row_k, col_k, k;
	double tmscore_row, tmscore_col, tmscore_row_1, tmscore_row_2, tmscore_col_1, tmscore_col_2;
	double tmscore_max = -1;
	double out_rmsd, out_SO;
        int out_Nali;

	double d02 = d0 * d0;
	double score_d82=score_d8*score_d8;
	double d2;
	for (iteration = 0; iteration < 1; iteration++){
		trans_OTscore(tmscore, otscore, xlen, ylen);
		//for (int ii = 0; ii < xlen; ii++){
		//	for (int jj = 0; jj < ylen; jj++){
		//		cout << otscore[ii][jj] << " ";
		//		}
		//	cout << endl;
		//}
		ot_gp_soi_egs(otscore, xlen, ylen, fwdmap, invmap, secx_bond, secy_bond, mm_opt);
		// row-based TMscore
		int row_k = 0;
		int row_score_k = 0;
		double row_rmsd = 0.0;
		double row_SO = 0.0;

		for (int i = 0; i < xlen; i++){
			int j = fwdmap[i];
			if (j < 0) continue;

			xtm_temp[row_k][0] = xa[i][0];
			xtm_temp[row_k][1] = xa[i][1];
			xtm_temp[row_k][2] = xa[i][2];

			ytm_temp[row_k][0] = ya[j][0];
			ytm_temp[row_k][1] = ya[j][1];
			ytm_temp[row_k][2] = ya[j][2];
			row_k++;
			row_rmsd += rmsd_matrix[i][j];
			if (rmsd_matrix[i][j] < score_d82){
				xtm_score_d8[row_score_k][0] = xa[i][0];
				xtm_score_d8[row_score_k][1] = xa[i][1];
				xtm_score_d8[row_score_k][2] = xa[i][2];

				ytm_score_d8[row_score_k][0] = ya[j][0];
				ytm_score_d8[row_score_k][1] = ya[j][1];
				ytm_score_d8[row_score_k][2] = ya[j][2];
				row_score_k++;
			}
			if (rmsd_matrix[i][j] < 12.25) row_SO = row_SO + 1.0;
		}
		row_rmsd = sqrt(row_rmsd / double(row_k));
		row_SO = row_SO / (double(minlen));

		tmscore_row_2 = TMscore8_search(r1, r2, xtm_score_d8, ytm_score_d8, xt, row_score_k, t, u, 40, 8, &rmsd, local_d0_search, xlen, score_d8, d0);
		tmscore_row_1 = TMscore8_search(r1, r2, xtm_score_d8, ytm_score_d8, xt, row_score_k, t, u, 40, 8, &rmsd, local_d0_search, ylen, score_d8, d0);
		tmscore_row = (tmscore_row_1 > tmscore_row_2)?(tmscore_row_1):(tmscore_row_2);
		cout << "Row Nali: " << row_k << "; " << "Row RMSD: " << row_rmsd << "; " << "Row TMscore: " << tmscore_row << "; " << "Row SO: " << row_SO << endl;
		cout << local_d0_search << ";" << Lnorm << ";" << score_d8 << ";" << d0 << endl;
		// col-based TMscore
		int col_k = 0;
		int col_score_k = 0;
		double col_rmsd = 0.0;
		double col_SO = 0.0;

		for (int j = 0; j < ylen; j++){
			int i = invmap[j];
			if (i < 0) continue;

			xtm_temp[col_k][0] = xa[i][0];
			xtm_temp[col_k][1] = xa[i][1];
			xtm_temp[col_k][2] = xa[i][2];

			ytm_temp[col_k][0] = ya[j][0];
			ytm_temp[col_k][1] = ya[j][1];
			ytm_temp[col_k][2] = ya[j][2];
			col_k++;
			col_rmsd += rmsd_matrix[i][j];
			if (rmsd_matrix[i][j] < score_d82){
                                xtm_score_d8[col_score_k][0] = xa[i][0];
                                xtm_score_d8[col_score_k][1] = xa[i][1];
                                xtm_score_d8[col_score_k][2] = xa[i][2];

                                ytm_score_d8[col_score_k][0] = ya[j][0];
                                ytm_score_d8[col_score_k][1] = ya[j][1];
                                ytm_score_d8[col_score_k][2] = ya[j][2];
                                col_score_k++;
                        }
			if (rmsd_matrix[i][j] < 12.25) col_SO = col_SO + 1.0;
		}
		col_rmsd = sqrt(col_rmsd / double(col_k));
		col_SO = col_SO / (double(minlen));

		tmscore_col_2 = TMscore8_search(r1, r2, xtm_score_d8, ytm_score_d8, xt, col_score_k, t, u, 40, 8, &rmsd, local_d0_search, xlen, score_d8, d0);
		tmscore_col_1 = TMscore8_search(r1, r2, xtm_score_d8, ytm_score_d8, xt, col_score_k, t, u, 40, 8, &rmsd, local_d0_search, ylen, score_d8, d0);
		tmscore_col = (tmscore_col_1 > tmscore_col_2)?(tmscore_col_1):(tmscore_col_2);
		cout << "Col Nali: " << col_k << "; " << "Col RMSD: " << col_rmsd << "; " << "Col TMscore: " << tmscore_col << "; " << "Col SO: " << col_SO << endl;
		if (tmscore_row > tmscore_col){
			k = 0;
			for (int j = 0; j < ylen; j++) invmap0[j] = -1;
			for (int i = 0; i < xlen; i++){
				fwdmap0[i] = fwdmap[i]; 
				int j = fwdmap[i];
				invmap0[j] = i;
                        	if (j < 0) continue;

                        	xtm[k][0] = xa[i][0];
                        	xtm[k][1] = xa[i][1];
                        	xtm[k][2] = xa[i][2];

                        	ytm[k][0] = ya[j][0];
                        	ytm[k][1] = ya[j][1];
                        	ytm[k][2] = ya[j][2];
                        	k++;
			}
			out_Nali = row_k;
			out_rmsd = row_rmsd;
			out_SO = row_SO;
			tmscore_max = tmscore_row;
			TM1 = tmscore_row_1;
			TM2 = tmscore_row_2;
		}else{
			k = 0;
			for (int i = 0; i < xlen; i++) fwdmap0[i] = -1;
			for (int j = 0; j < ylen; j++){
				invmap0[j] = invmap[j];
                                int i = invmap[j];
				fwdmap0[j] = i;
                                if (i < 0) continue;

                                xtm[k][0] = xa[i][0];
                                xtm[k][1] = xa[i][1];
                                xtm[k][2] = xa[i][2];

                                ytm[k][0] = ya[j][0];
                                ytm[k][1] = ya[j][1];
                                ytm[k][2] = ya[j][2];
                                k++;
                        }
			out_Nali = col_k;
			out_rmsd = col_rmsd;
			out_SO = col_SO;
			tmscore_max = tmscore_col;
			TM1 = tmscore_col_1;
			TM2 = tmscore_col_2;
		}
	}
	cout << "TMscore: " << tmscore_max << endl;
	DeleteArray(&xtm_temp, minlen);
	DeleteArray(&ytm_temp, minlen);
	DeleteArray(&otscore, xlen);
	DeleteArray(&xtm_score_d8, minlen);
	DeleteArray(&ytm_score_d8, minlen);

	fstream output_f;
        output_f.open("./output.txt", ios::app);
        output_f << out_Nali << " " << out_rmsd << " " << tmscore_max << " " << out_SO << endl;
        output_f.close();

	delete []invmap;
	delete []fwdmap;

	return tmscore_max;
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double SOI_iter(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **score, bool **path, double **val, double **xa, double **ya,
    int xlen, int ylen, double t[3], double u[3][3], int *invmap0,
    int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8,
    int **secx_bond, int **secy_bond, const int mm_opt, const bool init_invmap=false)
{
    double rmsd; 
    int *invmap=new int[ylen+1];
    
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;    
    tmscore_max=-1;

    //double d01=d0+1.5;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for (iteration=0; iteration<iteration_max; iteration++)
    {
        if (iteration==0 && init_invmap) 
            for (j=0;j<ylen;j++) invmap[j]=invmap0[j];
        else
        {
            for (j=0; j<ylen; j++) invmap[j]=-1;
            if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
        }
        soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);
    
        k=0;
        for (j=0; j<ylen; j++) 
        {
            i=invmap[j];
            if (i<0) continue;

            xtm[k][0]=xa[i][0];
            xtm[k][1]=xa[i][1];
            xtm[k][2]=xa[i][2];
            
            ytm[k][0]=ya[j][0];
            ytm[k][1]=ya[j][1];
            ytm[k][2]=ya[j][2];
            k++;
        }

        tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
            40, 8, &rmsd, local_d0_search, Lnorm, score_d8, d0);

        if (tmscore>tmscore_max)
        {
            tmscore_max=tmscore;
            for (j=0; j<ylen; j++) invmap0[j]=invmap[j];
        }
    
        if (iteration>0 && fabs(tmscore_old-tmscore)<0.000001) break;       
        tmscore_old=tmscore;
        do_rotation(xa, xt, xlen, t, u);
        SOI_super2score(xt, ya, xlen, ylen, score, d0, score_d8);
    }// for iteration
    
    delete []invmap;
    return tmscore_max;
}

void get_SOI_initial_assign(double **xk, double **yk, const int closeK_opt,
    double **score, bool **path, double **val, const int xlen, const int ylen,
    double t[3], double u[3][3], int invmap[], 
    double local_d0_search, double d0, double score_d8,
    int **secx_bond, int **secy_bond, const int mm_opt)
{
    int i,j,k;
    int minlen = min(xlen, ylen);
    double **xfrag;
    double **xtran;
    double **yfrag;
    double **otscore;
    NewArray(&xfrag, closeK_opt, 3);
    NewArray(&xtran, closeK_opt, 3);
    NewArray(&yfrag, closeK_opt, 3);
    NewArray(&otscore, xlen, ylen);
    int *fwdmap = new int[xlen];
    double rmsd;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;

    /* fill in score */
    for (i=0;i<xlen;i++)
    {
        for (k=0;k<closeK_opt;k++)
        {
            xfrag[k][0]=xk[i*closeK_opt+k][0];
            xfrag[k][1]=xk[i*closeK_opt+k][1];
            xfrag[k][2]=xk[i*closeK_opt+k][2];
        }

        for (j=0;j<ylen;j++)
        {
            for (k=0;k<closeK_opt;k++)
            {
                yfrag[k][0]=yk[j*closeK_opt+k][0];
                yfrag[k][1]=yk[j*closeK_opt+k][1];
                yfrag[k][2]=yk[j*closeK_opt+k][2];
            }
            Kabsch(xfrag, yfrag, closeK_opt, 1, &rmsd, t, u);
            do_rotation(xfrag, xtran, closeK_opt, t, u);
            
            //for (k=0; k<closeK_opt; k++)
            //{
                //d2=dist(xtran[k], yfrag[k]);
                //if (d2>score_d82) score[i+1][j+1]=0;
                //else score[i+1][j+1]=1./(1+d2/d02);
            //}
            k=closeK_opt-1;
            d2=dist(xtran[k], yfrag[k]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+d2/d02);
        }
    }

    /* initial assignment */
    for (j=0;j<ylen;j++) invmap[j]=-1;
    //if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
    //for (j=0; j<ylen;j++) i=invmap[j];
    //trans_OTscore(score, otscore, xlen, ylen);
    //ot_gp_soi_egs(score, xlen, ylen, fwdmap, invmap, secx_bond, secy_bond, mm_opt);

    /* clean up */
    DeleteArray(&xfrag, closeK_opt);
    DeleteArray(&xtran, closeK_opt);
    DeleteArray(&yfrag, closeK_opt);
    DeleteArray(&otscore, xlen);
    delete []fwdmap;
}

void SOI_assign2super(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **xa, double **ya,
    const int xlen, const int ylen, double t[3], double u[3][3], int invmap[], 
    double local_d0_search, double Lnorm, double d0, double score_d8)
{
    int i,j,k;
    double rmsd;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;

    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i<0) continue;
        xtm[k][0]=xa[i][0];
        xtm[k][1]=xa[i][1];
        xtm[k][2]=xa[i][2];

        ytm[k][0]=ya[j][0];
        ytm[k][1]=ya[j][1];
        ytm[k][2]=ya[j][2];
        k++;
    }
    TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
        40, 8, &rmsd, local_d0_search, Lnorm, score_d8, d0);
    do_rotation(xa, xt, xlen, t, u);
}

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */
int SOIalign_main(double **xa, double **ya,
    double **xk, double **yk, const int closeK_opt,
    const char *seqx, const char *seqy, const char *secx, const char *secy,
    double t0[3], double u0[3][3],
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA,
    int *invmap, double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen,
    const vector<string> sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, double *dist_list, 
    int **secx_bond, int **secy_bond, const int mm_opt, double &L_SO)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8,d0,d0_search,dcu0;//for TMscore search
    double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    double **score;       // Input score table for enhanced greedy search
    double **scoret;      // Transposed score table for enhanced greedy search
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  
    double **xtm, **ytm;  // for TMscore search engine
    double **xt;          //for saving the superposed version of r_1 or xtm
    double **yt;          //for saving the superposed version of r_2 or ytm
    double **r1, **r2;    // for Kabsch rotation
    double **match_score;
    double **rmsd_mat;

    /***********************/
    /* allocate memory     */
    /***********************/
    int minlen = min(xlen, ylen);
    int maxlen = (xlen>ylen)?xlen:ylen;
    NewArray(&score,  xlen+1, ylen+1);
    NewArray(&scoret, ylen+1, xlen+1);
    NewArray(&match_score, xlen+1, ylen+1);
    NewArray(&rmsd_mat, xlen, ylen);
    NewArray(&path, maxlen+1, maxlen+1);
    NewArray(&val,  maxlen+1, maxlen+1);
    NewArray(&xtm, maxlen, 3);
    NewArray(&ytm, maxlen, 3);
    NewArray(&xt, xlen, 3);
    NewArray(&yt, ylen, 3);
    NewArray(&r1, maxlen, 3);
    NewArray(&r2, maxlen, 3);

    /***********************/
    /*    parameter set    */
    /***********************/
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm, 
        score_d8, d0, d0_search, dcu0);
    //score_d8 = 3.5;
    int simplify_step    = 40; //for simplified search engine
    int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    int i,j;
    int *fwdmap0         = new int[xlen];
    int *invmap0         = new int[ylen];
    
    double TMmax=-1, TM=-1;
    for(i=0; i<xlen; i++) fwdmap0[i]=-1;
    for(j=0; j<ylen; j++) invmap0[j]=-1;
    double local_d0_search = d0_search;
    int iteration_max=(fast_opt)?2:30;
    //iteration_max=1;

    /*************************************************************/
    /* initial alignment with sequence order dependent alignment */
    /*************************************************************/
    vector<double> do_vec;
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt,
        mol_type,-1);
    do_vec.clear();
    if (mm_opt==6)
    {
        i=0;
        j=0;
        for (int r=0;r<seqxA.size();r++)
        {
            if (seqxA[r]=='*') // circular permutation point
            {
                for (int jj=0;jj<j;jj++) if (invmap0[jj]>=0)
                    invmap0[jj]+=xlen - i;
                i=0;
                continue;
            }
            if (seqyA[r]!='-')
            {
                if (seqxA[r]!='-') invmap0[j]=i;
                j++;
            }
            if (seqxA[r]!='-') i++;
        }
        for (j=0;j<ylen;j++)
        {
            i=invmap0[j];
            if (i>=0) fwdmap0[i]=j;
        }
    }
    do_rotation(xa, xt, xlen, t0, u0);

    SOI_super2rmsd(xt, ya, xlen, ylen, rmsd_mat);
    SOI_super2score(xt, ya, xlen, ylen, score, d0, score_d8);
    SOI_SPsuper2score(xt, ya, xlen, ylen, match_score);
    /*
    */
    TMmax = OT_TMiter(r1, r2, xtm, ytm, xt, match_score, path, val, xa, ya,
        xlen, ylen, t0, u0, invmap0, fwdmap0, iteration_max, 
    	local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt, rmsd_mat, TM1, TM2, true);

    for (int cc = 0; cc < ylen; ++cc) invmap[cc] = invmap0[cc];



    /* derive alignment from superposition */
    double d;
    int k;
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap0[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    //do_rotation(xa, xt, xlen, t, u);
    do_rotation(xa, xt, xlen, t0, u0);

    Liden=0;
    //double SO=0;
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap0[j];
        dist_list[j]=-1;
        if (i<0) continue;
        d=sqrt(dist(xt[i], ya[j]));
	if (d<3.5) L_SO += 1;
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        dist_list[j]=d;
        //SO+=(d<3.5);
        seqxA[j]=seqx[i];
        Liden+=(seqx[i]==seqy[j]);
    }
    //SO/=getmin(xlen,ylen);
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap0[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }
    //cout<<n_ali8<<'\t'
        //<<rmsd0<<'\t'
        //<<100.*SO<<endl;


    /* clean up */
    //*
    DeleteArray(&score, xlen+1);
    DeleteArray(&scoret,ylen+1);
    DeleteArray(&match_score, xlen+1);
    DeleteArray(&rmsd_mat, xlen);
    DeleteArray(&path,maxlen+1);
    DeleteArray(&val, maxlen+1);
    DeleteArray(&xtm, maxlen);
    DeleteArray(&ytm, maxlen);
    DeleteArray(&xt,xlen);
    DeleteArray(&yt,ylen);
    DeleteArray(&r1, maxlen);
    DeleteArray(&r2, maxlen);
    cout << "Finished" << endl;
    //delete[]invmap0;
    //delete[]fwdmap0;
    //delete[]m1;
    //delete[]m2;
    //* /
    
    return 0;
}
#endif

