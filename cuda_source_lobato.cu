#include "cufft.h"
#include <cstdio>
#include <cuda_runtime.h>
#include <cstdlib>
#include <iostream>
#include "device_launch_parameters.h"
#include <time.h>
#include <math.h>

// 计时
clock_t start, stop; 
double duration;
size_t avail, total;

#ifdef __linux__
 #define  CLK_TCK CLOCKS_PER_SEC
#endif

dim3 threadsPerBlock(32, 32);
size_t maxThreads;
int max_gpu_index = 0;


//初始化置零函数
__global__ void set_Zero(cufftComplex* cuda_result, int Height, int Width)
{
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < Height && col < Width)
	{
		cuda_result[row * Width + col].x = 0;
		cuda_result[row * Width + col].y = 0;
	}
}


//对w矩阵初始化赋值函数
__global__ void initKernel_W(cufftComplex* devPitchedPtr, size_t pitch, int Height, int Width, int Depth, int VOL,
                             double* init_real, double* init_img)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;

	if (swi_id < Depth)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)devPitchedPtr + swi_id * pitch);
		if (col < Width && row < Height)
		{
			rowHead[row * Width + col].x = (float)init_real[row * Width + col];
			rowHead[row * Width + col].y = (float)init_img[row * Width + col];
		}
	}
}

//组织cal_atom_ion_part1、exp_kxX_kyY1、sum_exp_prod执行顺序函数
__global__ void cal_atomORion_fun(cufftComplex* cuda_result, float* coff, double* s_2, int atom_nums
                                  , double* kx, double* ky, int Height, int Width, int VOL)
{
	//计算离子性或原子性
	//cal_atom_ion_part1计算后的结果存放在sum中
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int atom = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if(atom < atom_nums)
	{
		if (row < Height && col < Width)
		{
			float part1 = coff[atom * 14 + 3] * (2 + coff[atom * 14 + 4] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 14 + 4] * s_2[row * Width + col] * 4, 2);
			float part2 = coff[atom * 14 + 5] * (2 + coff[atom * 14 + 6] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 14 + 6] * s_2[row * Width + col] * 4, 2);
			float part3 = coff[atom * 14 + 7] * (2 + coff[atom * 14 + 8] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 14 + 8] * s_2[row * Width + col] * 4, 2);
			float part4 = coff[atom * 14 + 9] * (2 + coff[atom * 14 + 10] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 14 + 10] * s_2[row * Width + col] * 4, 2);
			float part5 = coff[atom * 14 + 11] * (2 + coff[atom * 14 + 12] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 14 + 12] * s_2[row * Width + col] * 4, 2);
			float sum = coff[atom * 14 + 2] * (part1 + part2 + part3 + part4 + part5) * expf((-1) * coff[atom * 14 + 13] * s_2[row * Width + col]);
			float temp = (-2) * (coff[atom * 14 + 0] * kx[row * Width + col] + coff[atom * 14 + 1] * ky[row * Width + col]
			) * 3.14159265;
			cuda_result[row * Width + col].x = cuda_result[row * Width + col].x + cosf(temp) * sum;
			cuda_result[row * Width + col].y = cuda_result[row * Width + col].y + sinf(temp) * sum;
		}
	}
}

__global__ void cal_absorb_fun(cufftComplex* cuda_result, float* coff, float* coff_absob, double* s_2, int atom_nums
                               , double* kx, double* ky, int Height, int Width, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int atom = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if(atom < atom_nums)
	{
		if (row < Height && col < Width)
		{
			
			float part1 = coff_absob[atom * 10 + 0] * expf((-1) * coff_absob[atom * 10 + 1] * s_2[row * Width + col]);
			float part2 = coff_absob[atom * 10 + 2] * expf((-1) * coff_absob[atom * 10 + 3] * s_2[row * Width + col]);
			float part3 = coff_absob[atom * 10 + 4] * expf((-1) * coff_absob[atom * 10 + 5] * s_2[row * Width + col]);
			float part4 = coff_absob[atom * 10 + 6] * expf((-1) * coff_absob[atom * 10 + 7] * s_2[row * Width + col]);
			float part5 = coff_absob[atom * 10 + 8] * expf((-1) * coff_absob[atom * 10 + 9] * s_2[row * Width + col]);
			float sum = coff[atom * 14 + 2] * (part1 + part2 + part3 + part4 + part5);
			float temp = (-2) * (coff[atom * 14 + 0] * kx[row * Width + col] + coff[atom * 14 + 1] * ky[row * Width + col]
			) * 3.14159265;
			cuda_result[row * Width + col].x = cuda_result[row * Width + col].x + (-1) * sinf(temp) * sum;
			cuda_result[row * Width + col].y = cuda_result[row * Width + col].y + cosf(temp) * sum;
		
		}
	}
}


//一层所有元素原子性或离子性加上对应的吸收值，结果保存在cuda_result中
__global__ void add_fun_absob(cufftComplex* cuda_result, cufftComplex* cuda_absorb, int Height, int Width)
{
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < Height && col < Width)
	{
		cuda_result[row * Width + col].x = cuda_result[row * Width + col].x + cuda_absorb[row * Width + col].x;
		cuda_result[row * Width + col].y = cuda_result[row * Width + col].y + cuda_absorb[row * Width + col].y;
	}
}


__global__ void cal_atomAndion_fun(cufftComplex* cuda_result,float* coff, double* s_2, int atomi_nums
                                   , double* kx, double* ky
                                   , int Height, int Width, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int atom = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if(atom < atomi_nums)
	{
		if (row < Height && col < Width)
		{
			float atom_part1 = coff[atom * 26 + 4] * (2 + coff[atom * 26 + 5] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 5] * s_2[row * Width + col] * 4, 2);
			float atom_part2 = coff[atom * 26 + 6] * (2 + coff[atom * 26 + 7] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 7] * s_2[row * Width + col] * 4, 2);
			float atom_part3 = coff[atom * 26 + 8] * (2 + coff[atom * 26 + 9] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 9] * s_2[row * Width + col] * 4, 2);
			float atom_part4 = coff[atom * 26 + 10] * (2 + coff[atom * 26 + 11] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 11] * s_2[row * Width + col] * 4, 2);
			float atom_part5 = coff[atom * 26 + 12] * (2 + coff[atom * 26 + 13] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 13] * s_2[row * Width + col] * 4, 2);
			float atom_sum = coff[atom * 26 + 3] * (atom_part1 + atom_part2 + atom_part3 + atom_part4 + atom_part5);

			float ion_part1 = coff[atom * 26 + 15] * (2 + coff[atom * 26 + 16] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 16] * s_2[row * Width + col] * 4, 2);
			float ion_part2 = coff[atom * 26 + 17] * (2 + coff[atom * 26 + 18] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 18] * s_2[row * Width + col] * 4, 2);
			float ion_part3 = coff[atom * 26 + 19] * (2 + coff[atom * 26 + 20] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 20] * s_2[row * Width + col] * 4, 2);
			float ion_part4 = coff[atom * 26 + 21] * (2 + coff[atom * 26 + 22] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 22] * s_2[row * Width + col] * 4, 2);
			float ion_part5 = coff[atom * 26 + 23] * (2 + coff[atom * 26 + 24] * s_2[row * Width + col] * 4) / powf(1 + coff[atom * 26 + 24] * s_2[row * Width + col] * 4, 2);
			float ion_sum = coff[atom * 26 + 14] * (ion_part1 + ion_part2 + ion_part3 + ion_part4 + ion_part5);

			float sum = coff[atom * 26 + 2] * (atom_sum + ion_sum) * expf((-1) * coff[atom * 26 + 25] * s_2[row * Width + col]);

			float temp = (-2) * (coff[atom * 26 + 0] * kx[row * Width + col] + coff[atom * 26 + 1] * ky[row * Width + col]
			) * 3.14159265;
			cuda_result[row * Width + col].x = cuda_result[row * Width + col].x + cosf(temp) * sum;
			cuda_result[row * Width + col].y = cuda_result[row * Width + col].y + sinf(temp) * sum;
			
		}
	}	
}

__global__ void cal_absorb_fun2(cufftComplex* cuda_result, float* coff, float* coff_absob, double* s_2, int atomi_nums
                                , double* kx, double* ky
                                , int Height, int Width, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int atom = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if(atom < atomi_nums)
	{
		if (row < Height && col < Width)
		{	
			float part1 = coff_absob[atom * 10 + 0] * expf((-1) * coff_absob[atom * 10 + 1] * s_2[row * Width + col]);
			float part2 = coff_absob[atom * 10 + 2] * expf((-1) * coff_absob[atom * 10 + 3] * s_2[row * Width + col]);
			float part3 = coff_absob[atom * 10 + 4] * expf((-1) * coff_absob[atom * 10 + 5] * s_2[row * Width + col]);
			float part4 = coff_absob[atom * 10 + 6] * expf((-1) * coff_absob[atom * 10 + 7] * s_2[row * Width + col]);
			float part5 = coff_absob[atom * 10 + 8] * expf((-1) * coff_absob[atom * 10 + 9] * s_2[row * Width + col]);
			float sum = coff[atom * 26 + 2] * (part1 + part2 + part3 + part4 + part5);
			float temp = (-2) * (coff[atom * 26 + 0] * kx[row * Width + col] + coff[atom * 26 + 1] * ky[row * Width + col]
			) * 3.14159265;
			cuda_result[row * Width + col].x = cuda_result[row * Width + col].x + (-1) * sinf(temp) * sum;
			cuda_result[row * Width + col].y = cuda_result[row * Width + col].y + cosf(temp) * sum;
		}
	}
}


__global__ void copy_result_to_p(cufftComplex* result, cufftComplex* P, int Height, int Width, float parameter)
{
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < Height && col < Width)
	{
		P[row * Width + col].x = result[row * Width + col].x * parameter + P[row * Width + col].x;
		P[row * Width + col].y = result[row * Width + col].y * parameter + P[row * Width + col].y;
	}
}

//对P进行ifftShift操作
__global__ void ifftShift(cufftComplex* P_PitchedPtr, size_t pitch, int P_Height, int P_Width, int P_Slices, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int d_id = id / (VOL / 2);
	int b_id = id % (VOL / 2) / (VOL / 4);
	int e_id = id % (VOL / 2) % (VOL / 4);

	int row = e_id / (P_Width / 2);
	int col = e_id % (P_Width / 2);

	int dest_row, dest_col;

	if (d_id < P_Slices)
	{
		if (b_id == 0)
		{
			// if(row == 0 && col == 0 && d_id == 0)
			// 	printf("invoke");
			dest_row = row + (P_Height / 2);
			dest_col = col + (P_Width / 2);
		}
		else if (b_id == 1)
		{
			col = col + (P_Width / 2);
			dest_row = row + (P_Height / 2);
			dest_col = col - (P_Width / 2);
		}
		// if(row == 0 && col == 2 && d_id == 0)
		// 	printf("(%d,%d,%d),(%d,%d),(%d,%d)\n",b_id,row,col,dest_row,dest_col,P_Height,P_Width);
		cufftComplex* rowHead = (cufftComplex *)((char *)P_PitchedPtr + d_id * pitch);
		float e1_real = rowHead[row * P_Width + col].x;
		float e1_img = rowHead[row * P_Width + col].y;

		float e2_real = rowHead[dest_row * P_Width + dest_col].x;
		float e2_img = rowHead[dest_row * P_Width + dest_col].y;

		rowHead[row * P_Width + col].x = e2_real;
		rowHead[row * P_Width + col].y = e2_img;

		rowHead[dest_row * P_Width + dest_col].x = e1_real;
		rowHead[dest_row * P_Width + dest_col].y = e1_img;
	}
}

//p与AperTrue矩阵点乘
__global__ void p_aperTrue_pointMul(cufftComplex* P_PitchedPtr, size_t pitch, double* aperTrue, int P_Height,
                                    int P_Width, int Slices, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int p_id = id / VOL;
	int row = (id % VOL) / P_Width;
	int col = (id % VOL) % P_Width;

	if (p_id < Slices)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)P_PitchedPtr + p_id * pitch);

		if (col < P_Width && row < P_Height)
		{
			rowHead[col + row * P_Width].x = rowHead[col + row * P_Width].x * aperTrue[col + row * P_Width];
			rowHead[col + row * P_Width].y = rowHead[col + row * P_Width].y * aperTrue[col + row * P_Width];
		}
	}
}

//对进行ifft2后的矩阵做归一化操作函数
__global__ void scaler(cufftComplex* devPitchedPtr, size_t pitch, int Height, int Width, int Depth, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;

	if (swi_id < Depth)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)devPitchedPtr + swi_id * pitch);

		if (col < Width && row < Height)
		{
			rowHead[col + row * Width].x = rowHead[col + row * Width].x / (Height * Width);
			rowHead[col + row * Width].y = rowHead[col + row * Width].y / (Height * Width);
		}
	}
}


//进行exp(i * simga * P)
__global__ void p_exp_sigma(cufftComplex* P_PitchedPtr, size_t pitch, int P_Height, int P_Width, int Slices, int VOL,
                            float sigma)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int p_id = id / VOL;
	int row = (id % VOL) / P_Width;
	int col = (id % VOL) % P_Width;

	if (p_id < Slices)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)P_PitchedPtr + p_id * pitch);

		if (col < P_Width && row < P_Height)
		{
			float real = (-1) * rowHead[col + row * P_Width].y * sigma;
			float img = rowHead[col + row * P_Width].x * sigma;

			rowHead[col + row * P_Width].x = expf(real) * cosf(img);
			rowHead[col + row * P_Width].y = expf(real) * sinf(img);
		}
	}
}

//与proj矩阵点乘操作函数
__global__ void w_p_pointMul(cufftComplex* W, cufftComplex* P, size_t pitch1
                             , int Height, int Width, int w_Depth, int VOL, int P_Width
                             , int beginRow, int beginCol
                             , int width_red, int step, int start)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	//定位到线程是在第几个w内
	int swi_id = id / VOL + start;
	//在每个w内部定位到是第几行
	int row = (id % VOL) / Width;
	//在每个w内部定位到是第几列
	int col = (id % VOL) % Width;
	int a = swi_id - start;
	if (a < w_Depth)
	{
		//获取到每一个点的W二维矩阵的首地址
		cufftComplex* rowHead1 = (cufftComplex *)((char *)W + a * pitch1);
		// 映射到大p里面第swi_id个小p的首地址
		int x_mat = beginRow + (swi_id / width_red) * step;
		int y_mat = beginCol + (swi_id % width_red) * step;
		//进行点乘操作
		float temp;
		if (col < Width && row < Height)
		{
			temp = rowHead1[row * Width + col].x;
			//小p与w做计算的元素在P的坐标
			int p_e_x = x_mat + row;
			int p_e_y = y_mat + col;
			//小p与w做计算的元素在P以一维存储时的下标
			int p_e = p_e_x * P_Width + p_e_y;
			int w_e = row * Width + col;

			rowHead1[w_e].x = temp * P[p_e].x - rowHead1[w_e].y * P[p_e].y;
			rowHead1[w_e].y = temp * P[p_e].y + rowHead1[w_e].y * P[p_e].x;
		}
	}
}

//与真空传播矩阵点乘操作函数
__global__ void w_constM_pointMul(cufftComplex* devPitchedPtr, size_t pitch, cufftComplex* constMat, int Height,
                                  int Width, int Depth, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if (swi_id < Depth)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)devPitchedPtr + swi_id * pitch);
		float temp;
		if (col < Width && row < Height)
		{
			// if (row==0 && col==0 && swi_id == 0)
			// {
			// 	printf("row:%d - col:%d\n",row,col);
			// 	printf("rowHead.x:%0.5lf + rowHead.y:%0.5f\n",rowHead[row * Width + col].x,rowHead[row * Width + col].y);
			// 	printf("constMat.x:%0.5lf + constMat.y:%0.5f\n",constMat[col + row * Width].x,constMat[row * Width + col].y);
			// }
			temp = rowHead[col + row * Width].x;
			rowHead[col + row * Width].x = temp * constMat[col + row * Width].x - rowHead[col + row * Width].y *
				constMat[col + row * Width].y;
			rowHead[col + row * Width].y = temp * constMat[col + row * Width].y + rowHead[col + row * Width].y *
				constMat[col + row * Width].x;
			//printf("rowHead.x:%0.5lf + rowHead.y:%0.5f\n",rowHead[row * Width + col].x,rowHead[row * Width + col].y);
		}
	}
}

__global__ void abs_w_final(cufftComplex* devPitchedPtr, double* aper_result, size_t pitch, size_t p, int Height,
                            int Width, int Depth, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;

	if (swi_id < Depth)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)devPitchedPtr + swi_id * pitch);
		double* rHead = (double *)((char *)aper_result + swi_id * p);
		if (col < Width && row < Height)
		{
			double temp1 = rowHead[col + row * Width].x;
			double temp2 = rowHead[col + row * Width].y;
			//hypotf为计算两个浮点数平方和的平方根
			rHead[col + row * Width] =  powf( hypotf(temp1, temp2), 2);
			// rowHead[col + row * Width].x = hypotf(temp1, temp2);
			// rowHead[col + row * Width].y = 0;
		}
	}
}


// w_final矩阵与aper矩阵进行点乘
__global__ void w_aper_pointMul(double* result, double* aper, size_t pitch, int Height, int Width, int Depth, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;
	if (swi_id < Depth)
	{
		double* rowHead = (double *)((char *)result + swi_id * pitch);
		double temp;
		if (col < Width && row < Height)
		{
			temp = rowHead[col + row * Width];
			rowHead[col + row * Width] = temp * aper[col + row * Width];
		}
	}
}


__global__ void absw_and_aperpm(cufftComplex* devPitchedPtr, double* aper_result, double* aper, size_t pitch, size_t p, int Height,
                            int Width, int Depth, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int swi_id = id / VOL;
	int row = (id % VOL) / Width;
	int col = (id % VOL) % Width;

	if (swi_id < Depth)
	{
		cufftComplex* rowHead = (cufftComplex *)((char *)devPitchedPtr + swi_id * pitch);
		double* rHead = (double *)((char *)aper_result + swi_id * p);
		if (col < Width && row < Height)
		{
			double temp1 = rowHead[col + row * Width].x;
			double temp2 = rowHead[col + row * Width].y;
			//hypotf为计算两个浮点数平方和的平方根
			rHead[col + row * Width] =  powf( hypotf(temp1, temp2), 2) * aper[col + row * Width];
		}
	}
}


//W_final与aper点乘以后的矩阵求和结果存放在result_sum_dev_head中
__global__ void sum_w_kernel(double* aper_result, size_t pitch
							 , int result_size, int thread_size, 
							 double* result_sum_dev_head, int Depth)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + 
	(threadIdx.y * blockDim.x + threadIdx.x);
	int result_id = id / thread_size; //开的线程数是thread_size*Depth，thread_size是处理一个result矩阵规约操作的线程数
	int r_e = id % thread_size; //r_e代表当前线程要处理的第swi_id个w的元素
	if (r_e + thread_size > result_size)
	{
		//代表已经超过w矩阵元素的个数了
		return;
	}
	double* result_rowHead = (double *)((char *)aper_result + result_id * pitch);
	if (result_id < Depth)
	{
		result_rowHead[r_e] = result_rowHead[r_e] + result_rowHead[r_e + thread_size];
		//printf("%d,%d,%lf\n",r_e,id,result_rowHead[r_e]);
	}
	if (thread_size == 1)
	{
		//代表每个result矩阵都只剩下最后两个元素的求和，所以只需要一个线程（即thread_size = 1）
		//但是实际核函数开了Depth个线程去处理Depth个result矩阵，此时id对应的就是第几个result
		// printf("%d,%d,%lf\n",r_e,id,result_rowHead[0]);
		result_sum_dev_head[id] = result_rowHead[0];
	}
}


void printf_2d(cufftComplex* w_PitchedPtr, int width, int height, int depth, size_t pitch)
{
	//w_s是所有组成的三维矩阵的主机端内存的首地址指针
	cufftComplex* w_p = (cufftComplex *)malloc(sizeof(cufftComplex) * width * height * depth);
	cudaMemcpy2D(w_p, width * height * sizeof(cufftComplex), w_PitchedPtr, pitch,
	             width * height * sizeof(cufftComplex), depth, cudaMemcpyDeviceToHost);
	//printf("%.12lf + (%.12lf)i\n", w_p[127*256+127].x, w_p[127*256+127].y);
	int k = 0;
	//for (int i = 496 * 235 + 273; i < 496 * 235 + 273 + 5; i++)
	for (int i = 0; i < 5; i++)
	{
		printf("%d================%.12lf + (%.12lf)i \n", i, w_p[i].x, w_p[i].y);
		k++;
		if ((k % (height * width)) == 0)
		{
			printf("\n=======================\n");
			printf("\n");
		}
	}
	if (w_p!=NULL)
	{
		free(w_p);
	}
}

void printf_1d(cufftComplex* cuda_result, int height, int width)
{
	cufftComplex* test = (cufftComplex *)malloc(sizeof(cufftComplex) * height * width);
	cudaMemcpy(test, cuda_result, height * sizeof(cufftComplex) * width, cudaMemcpyDeviceToHost);
	int k = 0;
	for (int i = 0; i < height * width; i++)
	{
		k++;
		printf("%lf + (%lf)i ", test[i].x, test[i].y);
		if (k % width == 0)
		{
			printf("\n");
		}
	}
	if (test!=NULL)
	{
		free(test);
	}
}

void printf_dd(double* cuda_result, int start, int len)
{
	double* test = (double *)malloc(sizeof(double) * len);
	cudaMemcpy(test, cuda_result+start, sizeof(double) * len, cudaMemcpyDeviceToHost);
	for (int i = 0; i < len; i++)
	{
		if (fabs(test[i]-0) > 1E-6){
			printf("printf_dd[%d]: %lf\n", start+i, test[i]);
		}
	}
	if (test!=NULL)
	{
		free(test);
	}
}



void setThreadGrid(int nums, int VOL, int* blockX, int* blockY)
{
	int blockThread = threadsPerBlock.x * threadsPerBlock.y;
	int blockSum = (nums * VOL  + blockThread - 1) / blockThread;
	if (blockSum < blockThread)
	{
		*blockX = 1;
		*blockY = blockSum;
	}
	else
	{
		*blockX = (blockSum + blockThread - 1) / blockThread;
		*blockY = (blockSum + *blockX - 1) / *blockX;
	}
}

size_t get_memory(int gpu_index)
{
	size_t gpu_free = 0;
	size_t gpu_total = 0;
	cudaSetDevice(gpu_index);
	cudaMemGetInfo(&gpu_free, &gpu_total);
	return gpu_free;
}

// 计算所有的p
void cal_p(double* atom_slice, int slice, double* absorb_n, int* sumSeries_n, double* proj_coff_mat, dim3 blockSize_P,
            cufftComplex* cuda_result1, int p_height, int p_width,
           double* atom_slice_i, double* absorb_n_i, double* proj_coff_mat_i,
           cufftComplex* cuda_result2, double* cuda_s_2, int* sumSeries_n_i,
           cufftComplex* p_PitchedPtr, size_t pitch2, double* cuda_kx, double* cuda_ky,
           double parameter, int *allAtoms, int *allAtomsi)
{
	//定义计算离子性和原子性矩阵所需要的参数矩阵
	float* coffMat;
	float* coffAbsob;
	float* coffMat_cuda;
	float* coffAbsob_cuda;

	//cuda_result_absorb保存考虑吸收
	cufftComplex* cuda_result_absorb;
	//一层中仅有离子性或原子性元素的叠加
	// slice从1开始
	int atom_nums = atom_slice[slice - 1];

	size_t ava = get_memory(max_gpu_index);
	printf("%d-------%d-------%zd MB\n",slice, atom_nums, ava/1024/1024);
	*allAtoms += atom_nums;
	int p_vol = p_width * p_height;
	int x, y;

	int p_batchs = maxThreads/p_vol;
	// printf("p_batchs: %d\n", p_batchs);
	if (atom_nums > 0)
	{
		int ci = atom_nums / p_batchs + 1;
		if (absorb_n == NULL)
		{
			int len = 14 * atom_nums;
			coffMat = (float *)malloc(sizeof(float) * len);
			//一次获取到当前传播层计算原子性和离子性矩阵的所有参数
			for (int m = 0; m < atom_nums; m++)
			{
				for (int n = 0; n < 14; n++)
					coffMat[m * 14 + n] = proj_coff_mat[*sumSeries_n * 14 + m * 14 + n];
			}
			cudaMalloc((void**)&coffMat_cuda, sizeof(float) * len);
			cudaMemcpy(coffMat_cuda, coffMat, sizeof(float) * len, cudaMemcpyHostToDevice);
			set_Zero << < blockSize_P, threadsPerBlock >> >(cuda_result1, p_height, p_width);
			//数据已从内存复制到显存中，故释放在内存中开辟的内存。
			if (coffMat!=NULL)
			{
				free(coffMat);
			}
			
			for(int i = 0;i < ci; i++)
			{
				int start_atom = i * p_batchs;
				int end_atom = (i+1) * p_batchs;
				if(end_atom > atom_nums)
				{
					end_atom = atom_nums;
				}
				int batch_atom = end_atom - start_atom;
				// printf("start_atom:%d end_atom:%d batch_atom:%d\n",start_atom, end_atom, batch_atom);
				if(batch_atom > 0)
				{
					setThreadGrid(batch_atom, p_vol, &x, &y);
					dim3 grid(x, y);
					cal_atomORion_fun << < grid, threadsPerBlock >> >(cuda_result1, coffMat_cuda+start_atom*14, cuda_s_2, batch_atom
																			, cuda_kx, cuda_ky, p_height, p_width, p_vol);
					cudaDeviceSynchronize();
				}
			}
		}
		else
		{
			cudaMalloc((void**)&cuda_result_absorb, sizeof(cufftComplex) * p_height * p_width);
			//计算该层系数所需要的空间长度
			int len1 = 14 * atom_nums;
			int len2 = 10 * atom_nums;
			coffMat = (float *)malloc(sizeof(float) * len1);
			coffAbsob = (float *)malloc(sizeof(float) * len2);

			//一次获取到当前传播层计算原子性和离子性矩阵的所有参数
			for (int m = 0; m < atom_nums; m++)
			{
				for (int n = 0; n < 14; n++)
					coffMat[m * 14 + n] = proj_coff_mat[*sumSeries_n * 14 + m * 14 + n];
				for (int n = 0; n < 10; n++)
					coffAbsob[m * 10 + n] = absorb_n[*sumSeries_n * 10 + m * 10 + n];
			}
			cudaMalloc((void**)&coffMat_cuda, sizeof(float) * len1);
			cudaMalloc((void**)&coffAbsob_cuda, sizeof(float) * len2);
			cudaMemcpy(coffMat_cuda, coffMat, sizeof(float) * len1, cudaMemcpyHostToDevice);
			cudaMemcpy(coffAbsob_cuda, coffAbsob, sizeof(float) * len2, cudaMemcpyHostToDevice);
			set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result1, p_height, p_width);
			set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result_absorb, p_height, p_width);
			//数据已从内存复制到显存中，故释放在内存中开辟的内存。
			if (coffMat!=NULL)
			{
				free(coffMat);
			}
			if (coffAbsob!=NULL)
			{
				free(coffAbsob);
			}
			for(int i = 0;i < ci; i++)
			{
				int start_atom = i * p_batchs;
				int end_atom = (i+1) * p_batchs;
				if(end_atom > atom_nums)
				{
					end_atom = atom_nums;
				}
				int batch_atom = end_atom - start_atom;
				if(batch_atom>0)
				{
					setThreadGrid(batch_atom, p_vol, &x, &y);
					dim3 grid(x, y);
					// 进行第silce层透射时，cal_atomORion_fun进行某个元素的 occupy_rate * ∑i(1-5)Ai*exp(-Bi * s^2).* exp(-(kx*X + ky*Y)*i)
					// 计算不吸收的一部分
					cal_atomORion_fun << <grid, threadsPerBlock >> >(cuda_result1, coffMat_cuda + start_atom * 14, cuda_s_2, batch_atom
																				, cuda_kx, cuda_ky, p_height, p_width, p_vol);
																				cudaDeviceSynchronize();
																				//cal_absorb_fun 计算 occupy_rate * ∑i(1-5)j * Ci * exp(-Di * s^2).* exp(-(kx*X + ky*Y)*i)
					cal_absorb_fun << <grid, threadsPerBlock >> >(cuda_result_absorb, coffMat_cuda + start_atom * 14, coffAbsob_cuda + start_atom * 10
						, cuda_s_2, batch_atom, cuda_kx, cuda_ky, p_height, p_width, p_vol);
						cudaDeviceSynchronize();
				}
			}		
		}
	}
	//考虑吸收，将非吸收项和吸收项相加
	if ((absorb_n != NULL) && (atom_nums != 0))
	{
		add_fun_absob << <blockSize_P, threadsPerBlock >> >(cuda_result1, cuda_result_absorb, p_height, p_width);
		cudaDeviceSynchronize();
		//考虑吸收，才释放相应的显存
		if (cuda_result_absorb!=NULL)
		{
			cudaFree(cuda_result_absorb);
		}
		
	}
	//释放当前层在显存开辟的空间
				
	if(coffMat_cuda!=NULL)
	{
		cudaFree(coffMat_cuda);
	}
	if (coffAbsob_cuda!=NULL)
	{
		cudaFree(coffAbsob_cuda);
	}
	//累加一层后的偏移值
	*sumSeries_n += atom_nums;
	int atomi_nums = atom_slice_i[slice - 1];
	
	allAtomsi += atomi_nums;
	//一层中所有既有原子性又有离子性元素的叠加
	if (atomi_nums > 0)
	{
		int ci = atomi_nums / p_batchs + 1;
		if (absorb_n_i == NULL)
		{
			//从proj_coff_mat_i中获取到计算所需要的所有系数
			//计算该层系数所需要的空间长度
			int len = 26 * atomi_nums;
			coffMat = (float *)malloc(sizeof(float) * len);
			//一次获取到当前传播层计算原子性和离子性矩阵的所有参数
			for (int m = 0; m < atomi_nums; m++)
			{
				for (int n = 0; n < 26; n++)
					coffMat[m * 26 + n] = proj_coff_mat_i[*sumSeries_n_i * 26 + m * 26 + n];
			}
			cudaMalloc((void**)&coffMat_cuda, sizeof(float) * len);
			cudaMemcpy(coffMat_cuda, coffMat, sizeof(float) * len, cudaMemcpyHostToDevice);
			set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result2, p_height, p_width);
			if (coffMat!=NULL)
			{
				free(coffMat);
			}

			for(int i = 0;i < ci; i++)
			{
				int start_atom = i * p_batchs;
				int end_atom = (i+1) * p_batchs;
				if(end_atom > atomi_nums)
				{
					end_atom = atomi_nums;
				}
				int batch_atom = end_atom - start_atom;
				if(batch_atom>0)
				{
					setThreadGrid(batch_atom, p_vol, &x, &y);
					dim3 grid(x, y);
					cal_atomAndion_fun << <grid, threadsPerBlock >> >(cuda_result2, coffMat_cuda + start_atom * 26, cuda_s_2, batch_atom
																		, cuda_kx, cuda_ky, p_height, p_width, p_vol);
																		cudaDeviceSynchronize();
				}
			}
		}
		else
		{
			//从proj_coff_mat和absorp中获取到进行计算所需要的所有系数
			cudaMalloc((void**)&cuda_result_absorb, sizeof(cufftComplex) * p_height * p_width);

			//计算该层系数所需要的空间长度
			int len1 = 26 * atomi_nums;
			int len2 = 10 * atomi_nums;

			coffMat = (float *)malloc(sizeof(float) * len1);
			coffAbsob = (float *)malloc(sizeof(float) * len2);

			//一次获取到当前传播层计算原子性和离子性矩阵的所有参数
			for (int m = 0; m < atomi_nums; m++)
			{
				for (int n = 0; n < 26; n++)
					coffMat[m * 26 + n] = proj_coff_mat_i[*sumSeries_n_i * 26 + m * 26 + n];
				for (int n = 0; n < 10; n++)
					coffAbsob[m * 10 + n] = absorb_n_i[*sumSeries_n_i * 10 + m * 10 + n];
			}

			cudaMalloc((void**)&coffMat_cuda, sizeof(float) * len1);
			cudaMalloc((void**)&coffAbsob_cuda, sizeof(float) * len2);
			cudaMemcpy(coffMat_cuda, coffMat, sizeof(float) * len1, cudaMemcpyHostToDevice);
			cudaMemcpy(coffAbsob_cuda, coffAbsob, sizeof(float) * len2, cudaMemcpyHostToDevice);
			set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result2, p_height, p_width);
			set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result_absorb, p_height, p_width);
			//数据已从内存复制到显存中，故释放在内存中开辟的内存。
			if (coffMat!=NULL)
			{
				free(coffMat);
			}
			if(coffAbsob!=NULL)
			{
				free(coffAbsob);
			}

			for(int i = 0;i < ci; i++)
			{
				int start_atom = i * p_batchs;
				int end_atom = (i+1) * p_batchs;
				if(end_atom > atomi_nums)
				{
					end_atom = atomi_nums;
				}
				int batch_atom = end_atom - start_atom;
				if(batch_atom>0)
				{
					setThreadGrid(batch_atom, p_vol, &x, &y);
					dim3 grid(x, y);
				
					cal_atomAndion_fun << <grid, threadsPerBlock >> >(
							cuda_result2, coffMat_cuda + start_atom * 26, cuda_s_2, batch_atom
						, cuda_kx, cuda_ky
						, p_height, p_width, p_vol);
						cudaDeviceSynchronize();
					cal_absorb_fun2 << <grid, threadsPerBlock >> >(
						cuda_result_absorb, coffMat_cuda + start_atom * 26, coffAbsob_cuda + start_atom * 10, cuda_s_2, batch_atom
						, cuda_kx, cuda_ky
						, p_height, p_width, p_vol);	
						cudaDeviceSynchronize();
				}
			}
		}
	}
	
	//如果存在吸收
	//将既有离子又有原子性的元素与它的吸收相加
	if ((absorb_n_i != NULL) && (atomi_nums != 0))
	{
		add_fun_absob << <blockSize_P, threadsPerBlock >> >(cuda_result2, cuda_result_absorb, p_height,
		                                                      p_width);
		//考虑吸收才释放相应的显存
		if(cuda_result_absorb!=NULL)
		{
			cudaFree(cuda_result_absorb);
		}
		
	}
	//释放当前层在显存开辟的空间
	if (coffMat_cuda!=NULL)
	{
		cudaFree(coffMat_cuda);
	}
	if(coffAbsob_cuda!=NULL)
	{
		cudaFree(coffAbsob_cuda);
	}
	//累计当前层后的偏移量
	*sumSeries_n_i += atomi_nums;

	//将一层所有仅有离子或原子性跟既有离子又有原子性的元素相加
	if (atom_nums == 0)
		set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result1, p_height, p_width);
	if (atomi_nums == 0)
		set_Zero << <blockSize_P, threadsPerBlock >> >(cuda_result2, p_height, p_width);
	if (atom_nums != 0 || atomi_nums != 0)
	{
		add_fun_absob << <blockSize_P, threadsPerBlock >> >(cuda_result1, cuda_result2, p_height, p_width);
		//将计算好的一层存入三维矩阵P中
		cufftComplex* rowHead2;
		rowHead2 = (cufftComplex *)((char *)p_PitchedPtr + pitch2 * (slice - 1));
		copy_result_to_p << <blockSize_P, threadsPerBlock >> >(cuda_result1, rowHead2, p_height, p_width,
																parameter);
	}
}


void initKernel_P(double* atom_slice, int slices, double* absorb_n, double* proj_coff_mat,
				  int p_height, int p_width,
                  double* atom_slice_i, double* absorb_n_i,
                  double* proj_coff_mat_i, double* cuda_s_2,
                  cufftComplex* p_PitchedPtr, size_t pitch2, double* cuda_kx,
                  double* cuda_ky, double parameter, double* cuda_aperTrue, double sigma, dim3 blockSize_P
)
{

	// cuda_result1保存一层中所有仅存在离子性或原子性的元素累加结果
	// cuda_resul2保存一层中所有既存在离子性又存在原子性的元素累加结果
	cufftComplex *cuda_result1, *cuda_result2;
	cudaMalloc((void**)&cuda_result1, sizeof(cufftComplex) * p_height * p_width);
	cudaMalloc((void**)&cuda_result2, sizeof(cufftComplex) * p_height * p_width);
	// 定义cal_atomORion_fun、plan、accumulate、set_Zero所需要的线程数目
	// dim3 threadsPerBlock_P(16, 16);
	// dim3 blockSize_P((p_height + threadsPerBlock.y - 1) / threadsPerBlock.y,
	//                  (p_width + threadsPerBlock.x - 1) / threadsPerBlock.x);
	// 记录每一层传播完成后proj_coff_mat和proj_coff_mat_i的下标所需要的移动的偏移量
	int sumSeries_n = 0;
	int sumSeries_n_i = 0;
	//程序开始计算Slices次透射传播的所有P
	int allAtoms=0, allAtomsi=0;
	for (int slice = 1; slice <= slices; slice++)
	{
		cal_p(atom_slice, slice, absorb_n, &sumSeries_n, proj_coff_mat, blockSize_P,
		      cuda_result1, p_height, p_width,
		      atom_slice_i, absorb_n_i, proj_coff_mat_i,
		      cuda_result2,  cuda_s_2, &sumSeries_n_i,
		      p_PitchedPtr, pitch2, cuda_kx, cuda_ky,
			  parameter, &allAtoms, &allAtomsi);
	}
	// printf_2d(p_PitchedPtr, p_width, p_height, slices, pitch2);
	printf("atoms: %d, atomsi: %d\n", allAtoms, allAtomsi);
	if(cuda_result1!=NULL)
	{
		cudaFree(cuda_result1);
	}
	if(cuda_result2!=NULL)
	{
		cudaFree(cuda_result2);
	}
	//定义P进行ifftShift所需要的线程数目
	// dim3 p_threadsPerBlock(32, 32);
	int p_blockX, p_blockY;
	//定义P与AperTrue点乘、scaler、p_exp_sigma所需要的线程数目
	setThreadGrid(slices, p_height * p_width, &p_blockX, &p_blockY);
	dim3 p1_dimGrid(p_blockX, p_blockY);
	p_aperTrue_pointMul << <p1_dimGrid, threadsPerBlock >> >(p_PitchedPtr, pitch2, cuda_aperTrue, p_height,
															   p_width, slices, p_height * p_width);
	setThreadGrid(slices, (p_height * p_width) / 2, &p_blockX, &p_blockY);
	dim3 p_dimGrid(p_blockX, p_blockY);
	ifftShift << <p_dimGrid, threadsPerBlock >> >(p_PitchedPtr, pitch2, p_height, p_width, slices,
	                                                p_height * p_width);
	cufftHandle p;
	cufftComplex* rowHead;
	
	for (int i = 0; i < slices; i++)
	{
		rowHead = (cufftComplex *)((char *)p_PitchedPtr + i * pitch2);
		cufftPlan2d(&p, p_height, p_width, CUFFT_C2C);
		cufftExecC2C(p, rowHead, rowHead, CUFFT_INVERSE);
		cufftDestroy(p);
	}
	p_exp_sigma << <p1_dimGrid, threadsPerBlock >> >(p_PitchedPtr, pitch2, p_height, p_width, slices,
	                                                   p_height * p_width, sigma / 1000);
}


size_t men_check(int gpu_count)
{
	size_t gpu_free = 0;
	size_t gpu_total = 0;
	for (int i = 0; i < gpu_count; i++)
	{
		cudaSetDevice(i);
		size_t free = 0;
		size_t total = 0;
		cudaMemGetInfo(&free, &total);
		gpu_free += free;
		gpu_total += total;
	}
	if (gpu_free > gpu_total)
	{
		printf("Out Of Memory\n");
		return 0;
	}
	return gpu_free;
}


__global__ void accumulate_result(double* result_sum, double* result, int w_num, int aper_num, int depth,
                                  double* gfsf)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	if (id < aper_num * depth)
	{
		double temp = 0;
		for (int j = 0; j < w_num; j++)
		{
			temp += result_sum[id + j * aper_num * depth] * gfsf[j];
		}
		result[id] = temp;
	}
}





__global__ void accumulate_result2(double* result_sum, double* result, int w_num, int aper_num, int depth, int len,
                                  double* gfsf)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	if (id < aper_num * depth * len)
	{
		double temp = 0;
		for (int j = 0; j < w_num; j++)
		{
			temp += result_sum[id + j * aper_num * depth * len] * gfsf[j];
		}
		result[id] = temp;
	}
}


void myProjectfun(double* my_real_w, double* my_img_w, double* aper,
 				  double* atom_slice, double* absorb_n,double* atom_slice_i,
                  double* absorb_n_i, double* proj_coff_mat, 
				  double* proj_coff_mat_i, double* kx, double* ky,double* s_2,
                  double* my_real_const, double* my_img_const,
				  double* gfsf, int height, int width,
                  int depth, int p_height, int p_width, int slices, int step, int beginrow,
                  int begincol, int width_red, double sigma, int w_num, 
				  int aper_num, double parameter,double* aper_true,
				  int* mid_layer, int layer_len, double *aper2,
				  double* series_n_corr, double* series_n_i_corr,
				  double* ele_n_corr, double* ele_n_i_corr,
				  double* corr_info_matrix, int cim_len,
				  double* return_result,double* mid_layer_mat, double* potentialx, double* potentialy)
{
	/*
   变量说明：
		  my_real_w:初始的w矩阵的实部矩阵,二维矩阵，维度为：height * width * w_num
		  my_imag_w:初始的w矩阵的虚部矩阵,二维矩阵，维度为：height * width * w_num
		  aper:维度为：height * width * aper_num
		  atom_slice:保存每层透射的仅有离子性或原子性元素的个数--> atom_slice[n1,n2,...,nk],维度为：Slices
		  atom_slice_i:保存每层透射的既有离子性又有原子性元素的个数--> atom_slice_i[n1,n2,...,nk],维度为：Slices
		  proj_coff_mat:保存每层透射仅有离子性或原子性元素计算相关系数
		  proj_coff_mat_i:保存每层透射既有离子性又有原子性元素计算相关系数
		  absorp_n:保存仅有离子性或原子性元素考虑吸收性的相关系数
		  absorp_n_i:保存仅有离子性或原子性元素考虑吸收性的相关系数
		  kx和ky:二维矩阵，维度为p_height * p_width
		  s_2:维度为p_height * p_width
		  my_real_const和my_img_const:分别是真空传播矩阵的实部矩阵和虚部矩阵，二维矩阵，维度：height * width
		  gfsf:w的比例系数，维度为：1*w_num
		  height：矩阵的行数
		  width:矩阵的列数
		  depth:一次要计算的点数
		  Slices:透射层数
		  step:在大P中的行列移动步长
		  beginRow和beginCol:大P中起始坐标的值
		  width_red:大P中横向移动的次数
		  sigma:用于计算P = exp(i*sigma*P)
		  w_num:不同的初始w矩阵的个数
		  aper_num：aper矩阵的个数
		  return_result：保存返回结果，维度为：w_num * aper_num * depth
		  mid_layer_mat：保存中间的结果，维度为：aper_num * depth * w_num *layer_len
*/

	printf("----------------CUDA START---------------\n");
	int gpu_count;
	cudaGetDeviceCount(&gpu_count);
	size_t gpu_free = 0;
	// 多GPU改进,只用一个GPU

	size_t max_gpu_free = 0;
	size_t max_gpu_total = 0;
	cudaDeviceProp prop;
	for(int i=0;i<gpu_count;i++)
    {
		cudaSetDevice(i);
   	 	cudaMemGetInfo( &avail, &total);
		if (avail > max_gpu_free){
			max_gpu_free = avail;
			max_gpu_total = total;
			max_gpu_index = i;
		}
	}
	cudaSetDevice(max_gpu_index);
	printf("gpu count:%d, select gpu: %d\n", gpu_count, max_gpu_index);
	cudaGetDeviceProperties(&prop,max_gpu_index);
	printf("name:%s\n",prop.name);
	printf("multiProcessorCount:%d\n",prop.multiProcessorCount);
	printf("maxThreadsPerBlock:%d\n",prop.maxThreadsPerBlock);
	printf("maxThreadsPerMultiProcessor:%d\n", prop.maxThreadsPerMultiProcessor);
	
	printf( "Max thread dimensions:  (%d, %d, %d)\n",
	prop.maxThreadsDim[0], prop.maxThreadsDim[1],
	prop.maxThreadsDim[2] );
	printf( "Max grid dimensions:  (%d, %d, %d)\n",
	prop.maxGridSize[0], prop.maxGridSize[1],
	prop.maxGridSize[2] );
	
	maxThreads = 8 * prop.maxThreadsDim[1]  * prop.maxGridSize[1];
	// printf("maxThreads: %zu\n", maxThreads);

	// printf("%d\n",prop.multiProcessorCount*prop.maxThreadsPerMultiProcessor);
	int maxThread = prop.maxThreadsPerBlock;
	// 每块使用最大线程数
	threadsPerBlock.x = sqrt(maxThread);
	threadsPerBlock.y = sqrt(maxThread);
	printf("max_gpu_free,max_gpu_total MB is %zu, %zu\n", max_gpu_free / 1024 / 1024, max_gpu_total / 1024 / 1024);
	printf("w_num:%d, slices:%d, aper_num:%d, layer_num:%d, depth:%d\n",w_num, slices, aper_num, layer_len, depth);
	printf("p_height:%d, p_width:%d, height:%d, width:%d\n", p_height, p_width, height, width);
	printf("mid_layer:");
	for(int i=0;i<layer_len;i++){
		printf("%d\t", mid_layer[i]);
	}
	printf("\n");
	//创建cufft执行计划变量p
	cufftHandle p;
	size_t pitch1, pitch2;
	size_t pitch3, pitch4;

	//用于在显存中保存Kx和Ky矩阵
	double *cuda_kx, *cuda_ky;
	//在显存中申请空间
	cudaMalloc((void**)&cuda_kx, sizeof(double) * p_height * p_width);
	cudaMalloc((void**)&cuda_ky, sizeof(double) * p_height * p_width);
	//将内存中的kx,ky复制到显存中
	cudaMemcpy(cuda_kx, kx, sizeof(double) * p_height * p_width, cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_ky, ky, sizeof(double) * p_height * p_width, cudaMemcpyHostToDevice);

	//w_PitchedPtr是所有W组成的三维矩阵在显存中的首地址指针
	//p_PitchedPtr是所有P组成的三维矩阵在显存中的首地址指针
	cufftComplex* w_PitchedPtr;
	cufftComplex* p_PitchedPtr;
	cudaMallocPitch((void**)&p_PitchedPtr, &pitch2, p_width * p_height * sizeof(cufftComplex), slices);
	long p_size = p_width * p_height * sizeof(cufftComplex) * slices;
	printf("p_size: %ld MB\n", p_size/1024/1024);
	//cuda_s_2在显存中用于保存s^2矩阵
	double *cuda_s_2;
	//为cuda_s_2和cuda_aperTrue变量申请显存
	cudaMalloc((void**)&cuda_s_2, sizeof(double) * p_height * p_width);
	cudaMemcpy(cuda_s_2, s_2, sizeof(double) * p_height * p_width, cudaMemcpyHostToDevice);

	//cuda_aperTrue在显存中保存aperTrue矩阵
	double* cuda_aperTrue;
	cudaMalloc((void**)&cuda_aperTrue, sizeof(double) * p_height * p_width);
	//将内存中cuda_aperTrue复制到显存中
	cudaMemcpy(cuda_aperTrue, aper_true, sizeof(double) * p_height * p_width, cudaMemcpyHostToDevice);
	dim3 blockSize_P((p_height + threadsPerBlock.y - 1) / threadsPerBlock.y,
	                 (p_width + threadsPerBlock.x - 1) / threadsPerBlock.x);
	//定义W这个三维矩阵中每一个二维的首指针
	cufftComplex* rowHead1;
	//定义大P这个三维矩阵中每一个二维的首指针
	cufftComplex* rowHead2;
	for (int i = 0; i < slices ; i++){
		rowHead2 = (cufftComplex *)((char *)p_PitchedPtr + pitch2 * i);
		set_Zero << < blockSize_P, threadsPerBlock >> >(rowHead2, p_height, p_width);
	}
	printf("potential start\n");
	start = clock();
	// 所有层的P
	initKernel_P(atom_slice, slices, absorb_n, proj_coff_mat,
	             p_height, p_width, atom_slice_i, absorb_n_i, proj_coff_mat_i,
	             cuda_s_2,p_PitchedPtr, pitch2, cuda_kx,
	             cuda_ky, parameter, cuda_aperTrue, sigma, blockSize_P);

	stop = clock();
	duration=(double)(stop-start)/CLK_TCK;
	printf("potential finished\ntime=%.2lf s\n", duration);
	
	// cufftComplex* potential_c = (cufftComplex*)malloc(p_size);
	// cudaMemcpy2D(potential_c, p_width * p_height * sizeof(cufftComplex), p_PitchedPtr, pitch2,
	//              p_width * p_height * sizeof(cufftComplex), slices, cudaMemcpyDeviceToHost);
	// for(int i=0;i<p_width * p_height * slices;i++){
	// 	potentialx[i] = potential_c[i].x;
	// 	potentialy[i] = potential_c[i].y;
	// }
	// printf("X-----%.6f, Y-------%.6f\n", potentialx[0], potentialy[0]);



	if (cuda_s_2 != NULL)
	{
		cudaFree(cuda_s_2);
	}
	if (cuda_kx != NULL)
	{
		cudaFree(cuda_kx);
	}
	if (cuda_ky != NULL)
	{	
		cudaFree(cuda_ky);
	}
	if (cuda_aperTrue != NULL)
	{
		cudaFree(cuda_aperTrue);
	}	
	cudaDeviceSynchronize();
	
	// 输出p正确

	// 传播变量开辟空间
	// 给w赋初值所需要的空间
	double *cuda_real, *cuda_img;
	cudaMalloc((void**)&cuda_real, sizeof(double) * height * width * w_num);
	cudaMalloc((void**)&cuda_img, sizeof(double) * height * width * w_num);
	//把W*w_num的值复制到显存中
	cudaMemcpy(cuda_real, my_real_w, sizeof(double) * height * width * w_num, cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_img, my_img_w, sizeof(double) * height * width * w_num, cudaMemcpyHostToDevice);

	//constMat为真空传播矩阵
	cufftComplex* constMat = (cufftComplex*)malloc(width * height * sizeof(cufftComplex));
	//t_constMat是显存中的首地址指针，constMat是主机端内存的首地址指针
	cufftComplex* cuda_constMat;
	cudaMallocPitch((void**)&cuda_constMat, &pitch4, width * sizeof(cufftComplex), height);
	//将外部传来的真空传播矩阵，赋值到显存中
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			constMat[i * width + j].x = my_real_const[i * width + j];
			constMat[i * width + j].y = my_img_const[i * width + j];
		}
	}
	cudaMemcpy(cuda_constMat, constMat, sizeof(cufftComplex) * height * width, cudaMemcpyHostToDevice);
	if(constMat != NULL)
	{
		free(constMat);
	}

	//给aper分配显存空间并把值从内存中复制到显存中	
	double* aper_cuda;
	cudaMalloc((void**)&aper_cuda, sizeof(double) * height * width * aper_num);
	cudaMemcpy(aper_cuda, aper, sizeof(double) * height * width * aper_num, cudaMemcpyHostToDevice);

	double* result_sum_dev;
	cudaMalloc((void**)&result_sum_dev, sizeof(double) * aper_num * depth * w_num);
	double* result_layer_dev;
	double* aper2_cuda;

	if (layer_len > 0)
	{
		cudaMalloc((void**)&result_layer_dev, sizeof(double) * aper_num * depth * w_num * layer_len);
		cudaMalloc((void**)&aper2_cuda, sizeof(double) * height * width * aper_num);
		cudaMemcpy(aper2_cuda, aper2, sizeof(double) * height * width * aper_num, cudaMemcpyHostToDevice);
	}

	double* aper_result;
	// dim3 threadsPerBlock(32, 32);
	int VOL = height * width;
	int W_threadsize = VOL / 2;
	int blockX, blockY;

	gpu_free  = get_memory(max_gpu_index);
	int w_batchs = gpu_free / VOL / (sizeof(cufftComplex) * w_num +  sizeof(double) * aper_num) >> 1;

	printf("GPU Free %zu MB\n", gpu_free / 1024 / 1024);
	if (gpu_free == 0)
	{
		return;
	}
	int n_cufft[2] = { width, height};
	// 如果超过显卡硬件内存，depth需要分批次进行传播
	// int w_batchs = gpu_free / VOL / sizeof(cufftComplex) / 6; //每次计算多少个点 
	int batch_n = depth / w_batchs; // 要计算多少次

	for (int w_series = 0; w_series < w_num; w_series++)
	{
		//假设这里要用w_series不同的w给W赋初值，这里需要得出每个w的首地址
		int size = height * width * w_series;
		double* cuda_real_head = cuda_real + size;
		double* cuda_img_head = cuda_img + size;
		// 将depth个点分批次进行运算,将w_series最后结果保存在第几个位置
		int rsd_num = w_series * aper_num * depth; //0 24 48
		
		for (int b = 0; b <= batch_n; b++)
		{
			// 20201112 保存的第几层
			int layer_i = 0;
			int batch_start = b * w_batchs; //0 6 9 8 2
			int batch_end = (b + 1) * w_batchs; //3 9 12 10
			if (batch_end > depth)
				batch_end = depth; // 8 8
			// 每次计算的点数
			int batch_size = batch_end - batch_start; //2 -1 0
			start = clock();
			if (batch_size > 0)
			{
				// 在这里开辟空间
				if (b == 0)
				{
					cudaDeviceSynchronize();
					cudaMallocPitch((void**)&w_PitchedPtr, &pitch1, VOL * sizeof(cufftComplex), batch_size);
					cudaDeviceSynchronize();
				}
				else if (b == batch_n)
				{
					cudaDeviceSynchronize();
					if (w_PitchedPtr != NULL)
					{
						cudaFree(w_PitchedPtr);
					}
					cudaDeviceSynchronize();
					cudaMallocPitch((void**)&w_PitchedPtr, &pitch1, VOL * sizeof(cufftComplex), batch_size);
					cudaDeviceSynchronize();
				}
				//定义进行w与p点乘、w与真空传播矩阵点乘、scaler、ifft2、fft2所需要的线程数目
				setThreadGrid(batch_size, VOL, &blockX, &blockY);
				dim3 blockSize(blockX, blockY);

				//给程序一开始的W赋值
				cudaDeviceSynchronize();
				initKernel_W << <blockSize, threadsPerBlock >> >(w_PitchedPtr, pitch1, height, width, batch_size, VOL,
				                                                 cuda_real_head,
				                                                 cuda_img_head);
				cudaDeviceSynchronize();
				// 计算每一层的传播直到w4
				for (int slice = 0; slice < slices; slice++)
				{
					// 计算每一层透射p的首地址
					rowHead2 = (cufftComplex *)((char *)p_PitchedPtr + pitch2 * slice);
					// 所有的点同时算slice层
					cudaDeviceSynchronize();
					w_p_pointMul << <blockSize, threadsPerBlock >> >(w_PitchedPtr, rowHead2, pitch1, height, width,
																	 batch_size, height * width, p_width, 
																	 beginrow, begincol, width_red, step, batch_start);
					cudaDeviceSynchronize();
					rowHead1 = (cufftComplex *)((char *)w_PitchedPtr);
					cufftPlanMany(&p, 2, n_cufft, NULL, 1, VOL, NULL, 1, VOL, CUFFT_C2C, batch_size);
					cufftExecC2C(p, rowHead1, rowHead1, CUFFT_FORWARD);
					cudaDeviceSynchronize();
					cufftDestroy(p);
					// 把data数据给w_PitchedPtr
					//将w1，进行fft2操作后的结果记为w2
					//将上一步的w2，与真空传播矩阵(这里设为constM)进行点乘操作后的结果记为w3
					w_constM_pointMul << <blockSize, threadsPerBlock >> >(
						w_PitchedPtr, pitch1, cuda_constMat, height, width,
						batch_size, height * width);
					cudaDeviceSynchronize();
					//将上一步的w3进行ifft2操作得到w4

					// 20201112 modify by ypj
					// 保存中间矩阵
					cudaDeviceSynchronize();
					if (layer_i < layer_len && (slice + 1) == mid_layer[layer_i])
					{	
						// printf("cur_mid_layer: %d\n", mid_layer[layer_i]);
						cudaMallocPitch((void**)&aper_result, &pitch3, sizeof(double) * VOL, batch_size);
						for (int i = 0; i < aper_num; i++)
						{
							//  获取到每个aper的首地址
							double* aper2_head = aper2_cuda + i * VOL; //第几个aper
							absw_and_aperpm<< <blockSize, threadsPerBlock >> >(w_PitchedPtr, aper_result, aper2_head, pitch1, pitch3, height, width, batch_size, VOL);
							// w_final矩阵与aper矩阵进行点乘保存到result_PitchedPtr中												
							cudaDeviceSynchronize();
							int cur_start = i * depth + batch_start + rsd_num + layer_i * aper_num * depth * w_num;
							double* result_layer_dev_head = result_layer_dev + cur_start;
							// 矩阵元素求和
							int h = log(VOL) / log(2); // 16
							if (pow(2, h) != VOL)
							{
								h++;
							}
							int thread_size = pow(2, h) / 2; // 2**15
							while (thread_size > 0)
							{
								int t1 = 1024;
								int t2 = thread_size / t1; // 2 ** 5 == 32
								if (thread_size / 1024 == 0)
								{
									t1 = thread_size;  // 2 ** 15
									t2 = 1; // 1
								}
								dim3 bthread(1, t1); //	(1, 2**15)
								dim3 g(batch_size, t2); // (2999, 1)
								// aper_result与aper点乘以后的矩阵求和结果存放在result_sum_dev_head中
								sum_w_kernel << < g, bthread >> >(aper_result, pitch3, VOL, thread_size
															, result_layer_dev_head, batch_size);	
								cudaDeviceSynchronize();
								thread_size /= 2; // 2的15次方除16次
							}
						} // aper循环结束
						cudaDeviceSynchronize();
						if (aper_result != NULL)
						{
							cudaFree(aper_result);
						}
						cudaDeviceSynchronize();
						layer_i += 1;
					}
					cudaDeviceSynchronize();
					rowHead1 = (cufftComplex *)((char *)w_PitchedPtr);
					cufftPlanMany(&p, 2, n_cufft, NULL, 1, VOL, NULL, 1, VOL, CUFFT_C2C, batch_size);
					cufftExecC2C(p, rowHead1, rowHead1, CUFFT_INVERSE);
					cudaDeviceSynchronize();
					cufftDestroy(p);
					cudaDeviceSynchronize();
					//因为cufft包本身设计的原因，上一步的得到w4矩阵所得到的数值的大小是真实数值的height*width??????????????г???height*width????,?????????????????
					scaler << <blockSize, threadsPerBlock >> >(w_PitchedPtr, pitch1, height, width, batch_size,
					                                           height * width);
					cudaDeviceSynchronize();

				} // 总层数计算完成
				//开始进行与不同的aper矩阵进行点乘操作，然后将点乘后的结果进行矩阵元素的累加
				//定义W进行ifftShift所需要的线程数目
				int w_blockX, w_blockY;
				setThreadGrid(batch_size, W_threadsize, &w_blockX, &w_blockY);
				dim3 w_dimGrid(w_blockX, w_blockY);

				// 进行fft2和fftshift
				cudaDeviceSynchronize();
				rowHead1 = (cufftComplex *)((char *)w_PitchedPtr);
				cufftPlanMany(&p, 2, n_cufft, NULL, 1, VOL, NULL, 1, VOL, CUFFT_C2C, batch_size);
				cufftExecC2C(p, rowHead1, rowHead1, CUFFT_FORWARD);
				cudaDeviceSynchronize();
				cufftDestroy(p);

				cudaDeviceSynchronize();
				ifftShift << <w_dimGrid, threadsPerBlock >> >(w_PitchedPtr, pitch1, height, width, batch_size, VOL);
				cudaDeviceSynchronize();
				// VOL * batch_size  ==> aper_num * depth * w_num
				cudaMallocPitch((void**)&aper_result, &pitch3, VOL * sizeof(double), batch_size);
				for (int i = 0; i < aper_num; i++)
				{
					//  获取到每个aper的首地址
					double* aper_head = aper_cuda + i * VOL; //第一个aper
					absw_and_aperpm<< <blockSize, threadsPerBlock >> >(w_PitchedPtr, aper_result, aper_head, pitch1, pitch3, height, width, batch_size, VOL);
					// w_final矩阵与aper矩阵进行点乘保存到result_PitchedPtr中				
					cudaDeviceSynchronize();
					// 矩阵元素求和
					int h = log(VOL) / log(2); // 16
					if (pow(2, h) != VOL)
					{
						h++;
					}
					int thread_size = pow(2, h) / 2;
					// aper_num * depth * w_num
					double* result_sum_dev_head = result_sum_dev + i * depth  + rsd_num + batch_start;
					while (thread_size > 0)
					{
						int t1 = 1024;
						int t2 = thread_size / t1;
						if (thread_size / 1024 == 0)
						{
							t1 = thread_size;
							t2 = 1;
						}
						dim3 bthread(1, t1); //	(1,1024)
						dim3 g(batch_size, t2); // (3,32)
						//	aper_result与aper点乘以后的矩阵求和结果存放在result_sum_dev_head中
						sum_w_kernel << < g, bthread >> >(aper_result, pitch3, VOL, thread_size
						                            , result_sum_dev_head, batch_size);
						cudaDeviceSynchronize();
						thread_size /= 2; // 2的15次方除16次
					}
					//printf_dd(result_sum_dev, aper_num * depth * w_num, i * depth  + rsd_num + batch_start);
				} // aper循环结束
				cudaDeviceSynchronize();
				if (aper_result != NULL)
				{
					cudaFree(aper_result);
				}
				cudaDeviceSynchronize();
			}
			stop = clock();
			duration=(double)(stop-start)/CLK_TCK;
			printf("w: (%d/%d), Pixel complete: (%d/%d), time per pixel: %.2lf sec\n", w_series+1, w_num, batch_start, depth, duration/batch_size);
		} // 分块循环结束
		cudaDeviceSynchronize();
		if (w_PitchedPtr != NULL)
		{
			cudaFree(w_PitchedPtr);
		}
		cudaDeviceSynchronize();
	} // w_num循环结束

	rowHead1 = NULL;
	rowHead2 = NULL;
	cufftDestroy(p);
	
	if (cuda_real!=NULL)
	{
		cudaFree(cuda_real);
	}
	if(cuda_img!=NULL)
	{
		cudaFree(cuda_img);
	}
	if(cuda_constMat!=NULL)
	{
		cudaFree(cuda_constMat);
	}
	if(p_PitchedPtr!=NULL)
	{
		cudaFree(p_PitchedPtr);
	}
	if(aper_cuda!=NULL)
	{
		cudaFree(aper_cuda);
	}
	if(aper2_cuda!=NULL)
	{
		cudaFree(aper2_cuda);
	}

	setThreadGrid(aper_num, depth, &blockX, &blockY);
	dim3 blockSize(blockX, blockY);
	double* cuda_temp2;
	cudaMalloc((void **)&cuda_temp2, sizeof(double) * aper_num * depth);
	double* cuda_gfsf;
	cudaMalloc((void**)&cuda_gfsf, sizeof(double) * w_num);
	cudaMemcpy(cuda_gfsf, gfsf, sizeof(double) * w_num, cudaMemcpyHostToDevice);
	accumulate_result << <blockSize, threadsPerBlock >>
		>(result_sum_dev, cuda_temp2, w_num, aper_num, depth, cuda_gfsf);
	double* result_sum = (double*)malloc(sizeof(double) * aper_num * depth);
	cudaMemcpy(result_sum, cuda_temp2, sizeof(double) * aper_num * depth, cudaMemcpyDeviceToHost);
	for (int i = 0; i < aper_num * depth; i++)
	{
		return_result[i] = result_sum[i];
	}
	// printf("\n");
	if (cuda_temp2!=NULL)
	{
		cudaFree(cuda_temp2);
	}
	if (result_sum!=NULL)
	{
		free(result_sum);
	}
	if(result_sum_dev!=NULL)
	{
		cudaFree(result_sum_dev);
	}
	// 20201114 create by ypj 
	if (layer_len > 0){
		double* cuda_temp3;
		cudaMalloc((void **)&cuda_temp3, sizeof(double) * aper_num * depth * layer_len);
		setThreadGrid(aper_num * layer_len, depth , &blockX, &blockY);
		dim3 blockSize(blockX, blockY);
		accumulate_result2 << <blockSize, threadsPerBlock >>
			>(result_layer_dev , cuda_temp3, 
			w_num, aper_num, depth, layer_len, cuda_gfsf);
		double* result_sum2 = (double*)malloc(sizeof(double) * aper_num * depth * layer_len);
		cudaMemcpy(result_sum2, cuda_temp3, sizeof(double) * aper_num * depth * layer_len, cudaMemcpyDeviceToHost);
		for (int i = 0; i < aper_num * depth * layer_len; i++)
		{
			mid_layer_mat[i] = result_sum2[i];
		}
		if (result_sum2!=NULL)
		{
			free(result_sum2);
		}
		if (cuda_temp3!=NULL)
		{
			cudaFree(cuda_temp3);
		}
		if (result_layer_dev!=NULL)
		{
			cudaFree(result_layer_dev);
		}
	}
	if(cuda_gfsf!=NULL)
	{
		cudaFree(cuda_gfsf);
	}
	printf("----------------CUDA Finished---------------\n");
}
