#include "mex.h"
#include "waveCuda.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    if (nrhs != 41)
    {
        mexErrMsgTxt("Invalid number of input arguments ");
        return;
    }
    if (nlhs != 4)
    {
        mexErrMsgTxt("Invalid number of outputs ");
        return;
    }

    double *my_real_w = (double *)mxGetPr(prhs[0]);
    double *my_img_w = (double *)mxGetPr(prhs[1]);

    double *aper = (double *)mxGetPr(prhs[2]);
    double *atom_slice = (double *)mxGetPr(prhs[3]);
    double *atom_slice_i = (double *)mxGetPr(prhs[5]);

    int absorp_flag = (int)mxGetScalar(prhs[30]);
    int absorp_i_flag = (int)mxGetScalar(prhs[31]);
    double *absorp_n = NULL;
    double *absorp_n_i = NULL;
    if (absorp_flag == 1)
        absorp_n = (double *)mxGetPr(prhs[4]);
    if (absorp_i_flag == 1)
        absorp_n_i = (double *)mxGetPr(prhs[6]);

    double *proj_coff_mat = (double *)mxGetPr(prhs[7]);
    double *proj_coff_mat_i = (double *)mxGetPr(prhs[8]);

    double *Kx = (double *)mxGetPr(prhs[9]);
    double *Ky = (double *)mxGetPr(prhs[10]);
    double *S_2 = (double *)mxGetPr(prhs[11]);

    double *my_real_const = (double *)mxGetPr(prhs[12]);
    double *my_img_const = (double *)mxGetPr(prhs[13]);

    double *gfsf = (double *)mxGetPr(prhs[14]);

    int Height = (int)mxGetScalar(prhs[15]);
    int Width = (int)mxGetScalar(prhs[16]);
    int Depth = (int)mxGetScalar(prhs[17]);

    int P_Height = (int)mxGetScalar(prhs[18]);
    int P_Width = (int)mxGetScalar(prhs[19]);

    int Slices = (int)mxGetScalar(prhs[20]);
    int Step = (int)mxGetScalar(prhs[21]);

    int beginRow = (int)mxGetScalar(prhs[22]);
    int beginCol = (int)mxGetScalar(prhs[23]);
    int width_red = (int)mxGetScalar(prhs[24]);
    double sigma = (double)mxGetScalar(prhs[25]);
    int w_num = (int)mxGetScalar(prhs[26]);
    int aper_num = (int)mxGetScalar(prhs[27]);
    double parameter = (double)mxGetScalar(prhs[28]);
    double *aperTrue = (double *)mxGetPr(prhs[29]);

    // 20201109 create by ypj
    int *mid_ceng = (int *)mxGetPr(prhs[32]);
    int ceng_len = (int)mxGetScalar(prhs[33]);
    double *aper2 = (double *)mxGetPr(prhs[34]);
    // corr
    double* seriesn_corr = (double *)mxGetPr(prhs[35]);
    double* seriesn_i_corr = (double *)mxGetPr(prhs[36]);
    double* elen_corr = (double *)mxGetPr(prhs[37]);
    double* elen_i_corr = (double *)mxGetPr(prhs[38]);
    double* corrinfo_matrix = NULL;
    int cimlen = (int)mxGetScalar(prhs[40]);
    if(cimlen != 0)
    {
        corrinfo_matrix = (double *)mxGetPr(prhs[39]);
    }

    
    
    
    //int flag = (int)mxGetScalar(prhs[35]);
    plhs[0] = mxCreateNumericMatrix(Depth, aper_num, mxDOUBLE_CLASS, mxREAL);
    double *return_result = (double *)mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(Depth * aper_num, ceng_len, mxDOUBLE_CLASS, mxREAL);
    double *mid_ceng_mat = (double *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericMatrix(P_Height * P_Width, Slices, mxDOUBLE_CLASS, mxREAL);
    double *potentialx = (double *)mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericMatrix(P_Height * P_Width, Slices, mxDOUBLE_CLASS, mxREAL);
    double *potentialy = (double *)mxGetPr(plhs[3]);
    //plhs[0] = mxCreateNumericMatrix(aper_num * Depth,w_num,mxSINGLE_CLASS,mxCOMPLEX);
    //float *return_real = (float *)mxGetPr(plhs[0]);
    //float *return_img = (float *)mxGetPi(plhs[0]);
    //double *result = (double *)malloc(sizeof(double) * aper_num * Depth);
    /*for(int i=0;i<Height * Width * w_num;i++)
            mexPrintf("%lf ",my_real_w[i]);*/
    myProjectfun(my_real_w, my_img_w, aper, atom_slice, absorp_n, atom_slice_i, absorp_n_i, proj_coff_mat, proj_coff_mat_i, Kx, Ky, S_2, my_real_const, my_img_const, gfsf, Height, Width, Depth, P_Height, P_Width, Slices, Step, beginRow - 1, beginCol - 1, width_red, sigma, w_num, aper_num, parameter, aperTrue, mid_ceng, ceng_len, aper2, seriesn_corr, seriesn_i_corr, elen_corr, elen_i_corr, corrinfo_matrix, cimlen, return_result, mid_ceng_mat, potentialx, potentialy);
    return;
}