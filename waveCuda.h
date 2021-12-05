#ifndef __CUDA_H__
#define __CUDA_H__

/*extern void myProjectfun(double *my_real_w,double *my_imag_w,double *atom_slice,double *proj_coff_mat
		,double *kx,double *ky,double *atom_desc_mat,double *ion_desc_mat,double *my_real_const,double *my_imag_const,int Height,int Width,int Depth,int Slices
		,int Kinds,int step,float *return_real,float *return_img);*/
/*extern void myProjectfun(double *my_real_w,double *my_imag_w,double *atom_slice,double *proj_coff_mat
		,double *kx,double *ky,double *s_2,double *my_real_const,double *my_imag_const,int Height,int Width,int Depth,int Slices
		,int step,float *return_real,float *return_img);*/
extern void myProjectfun(double *my_real_w, double *my_imag_w, double *aper, 
        double *atom_slice, double *absorp_n, double *atom_slice_i, 
        double *absorp_n_i, double *proj_coff_mat, double *proj_coff_mat_i, 
        double *kx, double *ky, double *s_2, double *my_real_const, 
        double *my_img_const, double *gfsf, int Height, int Width, int Depth, 
        int P_Height, int P_Width, int Slices, int step, int beginRow, int beginCol, 
        int width_red, double sigma, int w_num, int aper_num, double parameter, 
        double *AperTrue, int *mid_ceng, int ceng_len, double *aper2,
        double* series_n_corr, double* series_n_i_corr,
		double* ele_n_corr, double* ele_n_i_corr,
		double* corr_info_matrix, int cim_len, 
        double *return_result, double *mid_ceng_mat, double *potentialx, double *potentialy);
#endif
