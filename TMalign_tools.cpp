#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;

namespace py = pybind11;

//using c_array = py::array_t<double, py::array::c_style | py::array::forcecast>;

using c_array = py::array_t<double>;

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}



double dist(double x[3], double y[3])
{
    double d1=x[0]-y[0];
    double d2=x[1]-y[1];
    double d3=x[2]-y[2];
 
    return (d1*d1 + d2*d2 + d3*d3);
}


//PYBIND11_MODULE(TMalign, m) {
//    m.doc() = "pybind11 example plugin"; // optional module docstring

//    m.def("dist", &dist, "A function that adds two numbers");
//}


double dot(double *a, double *b)
{
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

//PYBIND11_MODULE(TMalign_tools, m) {
//    m.doc() = "pybind11 example plugin"; // optional module docstring

//    m.def("dot", &dot, "A function that adds two numbers");
//}


bool Kabsch(double **x, double **y, int n, int mode, double *rms,
    double t[3], double u[3][3])
{
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    double xc[3], yc[3];
    double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    double sqrt3 = 1.73205080756888, tol = 0.01;
    int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
    int ip2312[] = { 1, 2, 0, 1 };

    int a_failed = 0, b_failed = 0;
    double epsilon = 0.00000001;

    //initializtation
    *rms = 0;
    rms1 = 0;
    e0 = 0;
    double c1[3], c2[3];
    double s1[3], s2[3];
    double sx[3], sy[3], sz[3];
    for (i = 0; i < 3; i++)
    {
        s1[i] = 0.0;
        s2[i] = 0.0;

        sx[i] = 0.0;
        sy[i] = 0.0;
        sz[i] = 0.0;
    }

    for (i = 0; i<3; i++)
    {
        xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;
        for (j = 0; j<3; j++)
        {
            u[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            if (i == j)
            {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    if (n<1) return false;

    //compute centers for vector sets x, y
    for (i = 0; i<n; i++)
    {
        for (j = 0; j < 3; j++)
        {
            c1[j] = x[i][j];
            c2[j] = y[i][j];

            s1[j] += c1[j];
            s2[j] += c2[j];
        }

        for (j = 0; j < 3; j++)
        {
            sx[j] += c1[0] * c2[j];
            sy[j] += c1[1] * c2[j];
            sz[j] += c1[2] * c2[j];
        }
    }
    for (i = 0; i < 3; i++)
    {
        xc[i] = s1[i] / n;
        yc[i] = s2[i] / n;
    }
    if (mode == 2 || mode == 0)
        for (int mm = 0; mm < n; mm++)
            for (int nn = 0; nn < 3; nn++)
                e0 += (x[mm][nn] - xc[nn]) * (x[mm][nn] - xc[nn]) +
                      (y[mm][nn] - yc[nn]) * (y[mm][nn] - yc[nn]);
    for (j = 0; j < 3; j++)
    {
        r[j][0] = sx[j] - s1[0] * s2[j] / n;
        r[j][1] = sy[j] - s1[1] * s2[j] / n;
        r[j][2] = sz[j] - s1[2] * s2[j] / n;
    }

    //compute determinat of matrix r
    det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det;

    //compute tras(r)*r
    m = 0;
    for (j = 0; j<3; j++)
    {
        for (i = 0; i <= j; i++)
        {
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m++;
        }
    }

    double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
        - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
    det = det*det;

    for (i = 0; i<3; i++) e[i] = spur;

    if (spur>0)
    {
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det) / 2.0 - spur*h;

        if (h>0)
        {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if (d<0.0) d = 0.0;
            d = atan2(sqrt(d), -g) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0] = (spur + cth) + cth;
            e[1] = (spur - cth) + sth;
            e[2] = (spur - cth) - sth;

            if (mode != 0)
            {//compute a
                for (l = 0; l<3; l = l + 2)
                {
                    d = e[l];
                    ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                    if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
                    if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
                    if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
                    if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
                    if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
                    if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

                    if (fabs(ss[0]) >= fabs(ss[2]))
                    {
                        j = 0;
                        if (fabs(ss[0]) < fabs(ss[5])) j = 2;
                    }
                    else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
                    else j = 2;

                    d = 0.0;
                    j = 3 * j;
                    for (i = 0; i<3; i++)
                    {
                        k = ip[i + j];
                        a[i][l] = ss[k];
                        d = d + ss[k] * ss[k];
                    }


                    //if( d > 0.0 ) d = 1.0 / sqrt(d);
                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d = 0.0;
                    for (i = 0; i<3; i++) a[i][l] = a[i][l] * d;
                }//for l

                d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
                if ((e[0] - e[1]) >(e[1] - e[2]))
                {
                    m1 = 2;
                    m = 0;
                }
                else
                {
                    m1 = 0;
                    m = 2;
                }
                p = 0;
                for (i = 0; i<3; i++)
                {
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1] * a[i][m1];
                }
                if (p <= tol)
                {
                    p = 1.0;
                    for (i = 0; i<3; i++)
                    {
                        if (p < fabs(a[i][m])) continue;
                        p = fabs(a[i][m]);
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j + 1];
                    p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
                    if (p > tol)
                    {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    }
                    else a_failed = 1;
                }//if p<=tol
                else
                {
                    p = 1.0 / sqrt(p);
                    for (i = 0; i<3; i++) a[i][m1] = a[i][m1] * p;
                }//else p<=tol
                if (a_failed != 1)
                {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }
            }//if(mode!=0)
        }//h>0

        //compute b anyway
        if (mode != 0 && a_failed != 1)//a is computed correctly
        {
            //compute b
            for (l = 0; l<2; l++)
            {
                d = 0.0;
                for (i = 0; i<3; i++)
                {
                    b[i][l] = r[i][0] * a[0][l] +
                              r[i][1] * a[1][l] + r[i][2] * a[2][l];
                    d = d + b[i][l] * b[i][l];
                }
                //if( d > 0 ) d = 1.0 / sqrt(d);
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d = 0.0;
                for (i = 0; i<3; i++) b[i][l] = b[i][l] * d;
            }
            d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
            p = 0.0;

            for (i = 0; i<3; i++)
            {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1] * b[i][1];
            }

            if (p <= tol)
            {
                p = 1.0;
                for (i = 0; i<3; i++)
                {
                    if (p<fabs(b[i][0])) continue;
                    p = fabs(b[i][0]);
                    j = i;
                }
                k = ip2312[j];
                l = ip2312[j + 1];
                p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
                if (p > tol)
                {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0] / p;
                    b[l][1] = b[k][0] / p;
                }
                else b_failed = 1;
            }//if( p <= tol )
            else
            {
                p = 1.0 / sqrt(p);
                for (i = 0; i<3; i++) b[i][1] = b[i][1] * p;
            }
            if (b_failed != 1)
            {
                b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
                //compute u
                for (i = 0; i<3; i++)
                    for (j = 0; j<3; j++)
                        u[i][j] = b[i][0] * a[j][0] +
                                  b[i][1] * a[j][1] + b[i][2] * a[j][2];
            }

            //compute t
            for (i = 0; i<3; i++)
                t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                                                    u[i][2] * xc[2];
        }//if(mode!=0 && a_failed!=1)
    }//spur>0
    else //just compute t and errors
    {
        //compute t
        for (i = 0; i<3; i++)
            t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
                                                u[i][2] * xc[2];
    }//else spur>0

    //compute rms
    for (i = 0; i<3; i++)
    {
        if (e[i] < 0) e[i] = 0;
        e[i] = sqrt(e[i]);
    }
    d = e[2];
    if (sigma < 0.0) d = -d;
    d = (d + e[1]) + e[0];

    if (mode == 2 || mode == 0)
    {
        rms1 = (e0 - d) - d;
        if (rms1 < 0.0) rms1 = 0.0;
    }

    *rms = rms1;
    return true;
}


void transform(double t[3], double u[3][3], double *x, double *x1)
{
    x1[0]=t[0]+dot(&u[0][0], x);
    x1[1]=t[1]+dot(&u[1][0], x);
    x1[2]=t[2]+dot(&u[2][0], x);
}

//PYBIND11_MODULE(TMalign_tools, m) {
//   m.doc() = "pybind11 example plugin"; // optional module docstring
//    m.def("transform", *transform, "A function that adds two numbers");
//}


void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3])
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, &x[i][0], &x1[i][0]);
    }    
}

int score_fun8( double **xa, double **ya, int n_ali, double d, int i_ali[],
    double *score1, int score_sum_method, const double Lnorm,
    const double score_d8, const double d0)
{
    double score_sum=0, di;
    double d_tmp=d*d;
    double d02=d0*d0;
    double score_d8_cut = score_d8*score_d8;

    int i, n_cut, inc=0;

    while(1)
    {
        n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {
                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
            }
            else score_sum += 1/(1+di/d02);
        }
        //there are not enough feasible pairs, reliefe the threshold
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    *score1=score_sum/Lnorm;
    return n_cut;
}

double TMscore8_search(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, int Lali, double* t0, double  **u0, int simplify_step,
    int score_sum_method, double *Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0)
{
    int i, m;
    double score_max, score, rmsd;    
    const int kmax=Lali;    
    int k_ali[kmax], ka, k;
    double t[3];
    double u[3][3];
    double d;
    

    //iterative parameters
    int n_it=20;            //maximum number of iterations
    int n_init_max=6; //maximum number of different fragment length 
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;   

    int n_init=0, i_init;      
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }
    
    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment
    
    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;
      
        i=0;   
        while(1)
        {
            //extract the fragment starting from position i 
            ka=0;
            for(k=0; k<L_frag; k++)
            {
                int kk=k+i;
                r1[k][0]=xtm[kk][0];  
                r1[k][1]=xtm[kk][1]; 
                r1[k][2]=xtm[kk][2];   
                
                r2[k][0]=ytm[kk][0];  
                r2[k][1]=ytm[kk][1]; 
                r2[k][2]=ytm[kk][2];
                
                k_ali[ka]=kk;
                ka++;
            }
            
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if (simplify_step != 1)
                *Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);
            
            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, 
                score_sum_method, Lnorm, score_d8, d0);
            if(score>score_max)
            {
                score_max=score;
                
                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }
            
            //try to extend the alignment iteratively            
            d = local_d0_search + 1;
            for(int it=0; it<n_it; it++)            
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k][0]=xtm[m][0];  
                    r1[k][1]=xtm[m][1]; 
                    r1[k][2]=xtm[m][2];
                    
                    r2[k][0]=ytm[m][0];  
                    r2[k][1]=ytm[m][1]; 
                    r2[k][2]=ytm[m][2];
                    
                    k_ali[ka]=m;
                    ka++;
                } 
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, 
                    score_sum_method, Lnorm, score_d8, d0);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }                     
                }
                
                //check if it converges            
                if(n_cut==ka)
                {                
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k]) break;
                    }
                    if(k==n_cut) break;
                }                                                               
            } //for iteration            

            if(i<iL_max)
            {
                i=i+simplify_step; //shift the fragment        
                if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
            }
            else if(i>=iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}

static double* _vec_to_raw(const c_array &arr){
	py::buffer_info buf = arr.request();
	double *data = static_cast<double *>(buf.ptr);
	return data;
}

static std::vector<double *> _to_raw(const c_array &arr) {
  py::buffer_info buf = arr.request();
  double *ptr = static_cast<double *>(buf.ptr);

  std::vector<double *> data(arr.shape(0));

  for (int i = 0; i < arr.shape(0); i++) {
    data[i] = ptr;
    ptr += arr.shape(1);
  }

  return data;
}

double tmscore(c_array x, c_array y, c_array xt, int minlen, int Lali, c_array t0, c_array u0, int simplify_step, int score_sum_method, double local_d0_search, double Lnorm, double score_d8, double d0){
	double **r1, **r2;
	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);
	auto raw_t0 = _vec_to_raw(t0);
	auto raw_u0 = _to_raw(u0);
	//auto raw_r1 = _to_raw(r1);
	//auto raw_r2 = _to_raw(r2);
	auto raw_x = _to_raw(x);
	auto raw_y = _to_raw(y);
	auto raw_xt = _to_raw(xt);
	double Rcomm = 0;
	double TM_score;
	TM_score = TMscore8_search(r1, r2, raw_x.data(), raw_y.data(), raw_xt.data(), Lali, raw_t0, raw_u0.data(), simplify_step, score_sum_method, &Rcomm, local_d0_search, Lnorm, score_d8, d0);
	DeleteArray(&r1, minlen);
        DeleteArray(&r2, minlen);
	return TM_score;
}

PYBIND11_MODULE(TMalign_tools, m){
	m.doc() = "TMscore calculation";
	m.def("tmscore", &tmscore, "TMscore calcualtion");
}







//PYBIND11_MODULE(TMalign, m) {
//    m.doc() = "pybind11 example plugin"; // optional module docstring

//    m.def("TMscore8_search", &TMscore8_search, "A function that adds two numbers");
//}
