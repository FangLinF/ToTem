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
		P[row * Width + col].x = result[row * Width + col].x * parameter;
		P[row * Width + col].y = result[row * Width + col].y * parameter;
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
	// size_t gpu_free = 0;
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
	// cufftHandle p;
	size_t pitch2;

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
	// cufftComplex* w_PitchedPtr;
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
	// cufftComplex* rowHead1;
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
	
	cufftComplex* potential_c = (cufftComplex*)malloc(p_size);
	cudaMemcpy2D(potential_c, p_width * p_height * sizeof(cufftComplex), p_PitchedPtr, pitch2,
	             p_width * p_height * sizeof(cufftComplex), slices, cudaMemcpyDeviceToHost);
	for(int i=0;i<p_width * p_height * slices;i++){
		potentialx[i] = potential_c[i].x;
		potentialy[i] = potential_c[i].y;
	}
	// printf("X-----%.6f, Y-------%.6f\n", potentialx[0], potentialy[0]);

	if (potential_c != NULL){
		free(potential_c);
	}
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
	
	printf("----------------CUDA Finished---------------\n");
}
