import multiprocessing
import pycuda.driver as cuda
import numpy as np
from pycuda.compiler import SourceModule
from pycuda import gpuarray, tools, autoinit, curandom
import scipy.io as scio
import skcuda.fft as cu_fft
from matplotlib import pyplot as plt
import sys
import os
import time
from multiprocessing import Process, Queue, Manager, Pool
import multiprocessing as mp
import tifffile as tiff

block = (32, 32, 1)
batches = 1
file_name = 0
N = 320


class ComplexMatrix:
    def __init__(self, a, b):
        self.nums, self.row, self.col = a.shape
        self.a = a
        self.b = b

    def add(self):
        c = np.zeros_like(self.a, dtype=np.complex128)
        for i in range(self.nums):
            for j in range(self.row):
                c[i][j] = list(map(lambda x, y: np.complex(
                    x, y), self.a[i][j], self.b[i][j]))
        return c

    def sub(self):
        c = np.zeros_like(self.a, dtype=np.complex128)
        for i in range(self.nums):
            for j in range(self.row):
                c[i][j] = list(map(lambda x, y: np.complex(
                    x, -y), self.a[i][j], self.b[i][j]))
        return c


class stem:
    def __init__(self, para_part1, para_part2, myrandrange, U, V, mywave, nums):
        self.para_part1 = para_part1
        self.pp2 = para_part2.astype(np.float64)
        self.mrr_list = myrandrange.astype(np.float64)
        self.u = U.astype(np.float64)
        self.v = V.astype(np.float64)
        self.mywave = mywave.astype(np.complex128)
        self.NUMS = nums
        global N
        N = self.mywave.shape[0]

    def setThreadGrid(self, nums, VOL):
        blockThread = block[0] * block[1]
        blockSum = (nums * VOL + blockThread - 1) // blockThread
        if (blockSum < blockThread):
            blockX = 1
            blockY = blockSum
        else:
            blockX = (blockSum + blockThread - 1) // blockThread
            blockY = (blockSum + blockX - 1) // blockX
        return blockX, blockY

    def WholeTCC2D_newTEM_forpar(self, pp1, para, batches):
        len_p = len(para[0])
        u_gpu = gpuarray.to_gpu(self.u)
        v_gpu = gpuarray.to_gpu(self.v)
        a = pp1['lambda'] * self.u
        b = pp1['lambda'] * self.v
        oumigau = np.repeat(a[np.newaxis, :, :], batches, axis=0)
        oumigav = np.repeat(b[np.newaxis, :, :], batches, axis=0)
        x_mat, y_mat = self.setThreadGrid(batches, N * N)
        oumiga = ComplexMatrix(oumigau, oumigav)
        oumiga_a = oumiga.add()
        oumiga_s = oumiga.sub()
        oa = gpuarray.to_gpu(oumiga_a)
        os = gpuarray.to_gpu(oumiga_s)
        # 提取变量
        oa2 = oa.__mul__(oa)
        oa3 = oa2.__mul__(oa)
        oa4 = oa2.__mul__(oa2)
        oa5 = oa3.__mul__(oa2)
        os2 = os.__mul__(os)
        os3 = os2.__mul__(os)
        os4 = os2.__mul__(os2)
        os5 = os3.__mul__(os2)
        os6 = os5.__mul__(os)

        kai = gpuarray.zeros((batches, N, N), dtype=np.float64)

        self.getKai(
            kai, os3, oa2.__mul__(os), os4, os.__mul__(
                oa3), oa2.__mul__(os2), os5,
            oa4.__mul__(os), oa3.__mul__(os2), os6, oa4.__mul__(os2),
            oa3.__mul__(os3), oa5.__mul__(os),
            cuda.In(para), np.int32(len_p), np.int32(N), np.int32(N),
            np.int32(N * N), np.int32(batches),
            block=block, grid=self.setThreadGrid(batches, N * N)
        )

        W = gpuarray.zeros_like(kai)

        self.getW(
            W, kai, os2, oa.__mul__(os), cuda.In(para),
            np.int32(len_p), np.int32(N), np.int32(N),
            np.int32(N * N), np.int32(batches),
            block=block, grid=(x_mat, y_mat)
        )

        gradient_kai_u = gpuarray.zeros_like(kai)
        self.gradient_u(
            gradient_kai_u, np.float64(pp1['lambda']), cuda.In(para),
            np.int32(len_p), u_gpu, v_gpu,
            np.int32(N), np.int32(N), np.int32(
                N * N), np.int32(batches),
            block=block, grid=(x_mat, y_mat)
        )

        # gradient_kai_u = gradient_kai_u.transpose((0, 2, 1))
        gradient_kai_v = gpuarray.zeros_like(kai)
        self.gradient_v(
            gradient_kai_v,
            np.float64(pp1['lambda']), cuda.In(para), np.int32(len_p),
            u_gpu, v_gpu, np.int32(N), np.int32(N),
            np.int32(N * N), np.int32(batches),
            block=block, grid=(x_mat, y_mat)
        )
        # gradient_kai_v = gradient_kai_v.transpose((0, 2, 1))

        E_s_coh = gpuarray.zeros_like(kai)
        self.getESCoh(
            E_s_coh, np.float64(pp1['alafa']), np.float64(pp1['lambda']),
            gradient_kai_u, gradient_kai_v,
            np.int32(N * N), np.int32(batches), np.float64(np.pi),
            block=block, grid=(x_mat, y_mat)
        )
        # E_s_coh = E_s_coh.transpose((0, 2, 1))
        E_s_coh_cpu = E_s_coh.transpose((0, 2, 1)).get().copy()
        E_s_coh = gpuarray.to_gpu(E_s_coh_cpu)
        t = gpuarray.zeros_like(kai, dtype=np.complex128)
        self.getT(t, W, np.float64(pp1['lambda']),
                  np.int32(N * N), np.int32(batches), np.float64(np.pi),
                  block=block, grid=(x_mat, y_mat)
                  )

        mtotal = np.int32(pp1['mtotal'])
        tlittle_all = gpuarray.zeros(
            (mtotal * batches, N, N), dtype=np.complex128)

        self.getTlittleAll(
            tlittle_all, E_s_coh, t, mtotal, np.float64(np.pi),
            np.float64(pp1['delta_yita']), u_gpu, v_gpu, np.float64(
                pp1['lambda']),
            np.int32(N), np.int32(N), np.int32(
                N * N), np.int32(batches),
            block=block, grid=(x_mat, y_mat)
        )
        return tlittle_all

    def producer(self, q, bat, d=0):
        cuda.init()
        device = cuda.Device(d)
        ctx = device.make_context()

        mod = SourceModule("""
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cufft.h>


__global__ void createNewPP2(double* result, double* pp2, double* mrr,double* ai,
int width, int batches, double pi)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) +
        (threadIdx.y * blockDim.x + threadIdx.x);
    if (id < batches)
    {
        result[id*width+0] = pp2[0];
        result[id*width+1] = pp2[1]+(ai[id*21+0]-0.5)*2*mrr[0];
        result[id*width+2] = pp2[2]+(ai[id*21+1]-0.5)*2*mrr[1];
        result[id*width+3] = pp2[3]+(ai[id*21+2]-0.5)*2*mrr[2];
        result[id*width+4] = pp2[4]+(ai[id*21+3]-0.5)*2*mrr[3];
        result[id*width+5] = pp2[5]+(ai[id*21+4]-0.5)*2*mrr[4];
        result[id*width+6] = pp2[6]+(ai[id*21+5]-0.5)*2*mrr[5];
        result[id*width+7] = pp2[7]+(ai[id*21+6]-0.5)*2*mrr[6];
        result[id*width+8] = pp2[8]+(ai[id*21+7]-0.5)*2*mrr[7];
        result[id*width+9] = pp2[9]+(ai[id*21+8]-0.5)*2*mrr[8];

        result[id*width+10] = pp2[10]+(ai[id*21+9]-0.5)*2*mrr[9];
        result[id*width+11] = pp2[11]+(ai[id*21+10]-0.5)*2*mrr[10];
        result[id*width+12] = pp2[12]+(ai[id*21+11]-0.5)*2*mrr[11];

        result[id*width+13] = pp2[13]+(ai[id*21+12]-0.5)*2*mrr[12];
        result[id*width+14] = pp2[14]+(ai[id*21+13]-0.5)*2*mrr[13];
        result[id*width+15] = pp2[15]+(ai[id*21+14]-0.5)*2*mrr[14];
        result[id*width+16] = pp2[16]+(ai[id*21+15]-0.5)*2*mrr[15];

        result[id*width+17] = pp2[17]+(ai[id*21+16]-0.5)*2*mrr[16];
        result[id*width+18] = pp2[18]+(ai[id*21+17]-0.5)*2*mrr[17];

        result[id*width+19] = pp2[19]+(ai[id*21+18]-0.5)*2*mrr[18];
        result[id*width+20] = pp2[20]+(ai[id*21+19]-0.5)*2*mrr[19];
        result[id*width+21] = pp2[21]+(ai[id*21+20]-0.5)*2*mrr[20];
        
        // S5 - phiD5
        result[id*width+22] = pp2[22];
        result[id*width+23] = pp2[23];
        result[id*width+24] = pp2[24];
        result[id*width+25] = pp2[25];

        result[id*width+26] = result[id*width+2]*cos(result[id*width+3]/180*pi);
        result[id*width+27] = result[id*width+2]*sin(result[id*width+3]/180*pi);
        result[id*width+28] = result[id*width+4]*cos(result[id*width+5]/180*pi);
        result[id*width+29] = result[id*width+4]*sin(result[id*width+5]/180*pi);
        result[id*width+30] = result[id*width+6]*cos(result[id*width+7]/180*pi);
        result[id*width+31] = result[id*width+6]*sin(result[id*width+7]/180*pi);
        result[id*width+32] = result[id*width+9]*cos(result[id*width+10]/180*pi);
        result[id*width+33] = result[id*width+9]*sin(result[id*width+10]/180*pi);
        result[id*width+34] = result[id*width+11]*cos(result[id*width+12]/180*pi);
        result[id*width+35] = result[id*width+11]*sin(result[id*width+12]/180*pi);
        result[id*width+36] = result[id*width+8];
        result[id*width+37] = result[id*width+13]*cos(result[id*width+14]/180*pi);
        result[id*width+38] = result[id*width+13]*sin(result[id*width+14]/180*pi);
        result[id*width+39] = result[id*width+15]*cos(result[id*width+16]/180*pi);
        result[id*width+40] = result[id*width+15]*sin(result[id*width+16]/180*pi);
        result[id*width+41] = result[id*width+17]*cos(result[id*width+18]/180*pi);
        result[id*width+42] = result[id*width+17]*sin(result[id*width+18]/180*pi);
        result[id*width+43] = result[id*width+19]*cos(result[id*width+20]/180*pi);
        result[id*width+44] = result[id*width+19]*sin(result[id*width+20]/180*pi);
        result[id*width+45] = result[id*width+22]*cos(result[id*width+23]/180*pi);
        result[id*width+46] = result[id*width+22]*sin(result[id*width+23]/180*pi);
        result[id*width+47] = result[id*width+24]*cos(result[id*width+25]/180*pi);
        result[id*width+48] = result[id*width+24]*sin(result[id*width+25]/180*pi);
    }
}

__global__ void fun1(double* result, double* uv, double* para, int width, int height, int plen, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int cur_id = id / VOL;
	int row = (id % VOL) / width;
	int col = (id % VOL) % width;
	if (cur_id < batches)
	{
        double pa = para[cur_id * plen];
        if(row < height && col < width)
        {
            result[row * width + col] = pa * uv[row * width + col];
        }
	}
}
__global__ void getKai(double* result,
    cufftDoubleComplex* a, cufftDoubleComplex* b, cufftDoubleComplex* c, cufftDoubleComplex* d, cufftDoubleComplex* e,
    cufftDoubleComplex* f, cufftDoubleComplex* g, cufftDoubleComplex* h, cufftDoubleComplex* i, cufftDoubleComplex* j,
    cufftDoubleComplex* k, cufftDoubleComplex* l,
    double* para, int plen, int width, int height, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int cur_id = id / VOL;
	int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if (cur_id < batches)
    {
        double t1 = para[cur_id * plen + 28] / 3;
        double t2 = para[cur_id * plen + 29] / 3;
        double t3 = para[cur_id * plen + 30];
        double t4 = para[cur_id * plen + 31];
        double t5 = para[cur_id * plen + 32] / 4;
        double t6 = para[cur_id * plen + 33] / 4;
        double t7 = para[cur_id * plen + 34];
        double t8 = para[cur_id * plen + 35];
        double t9 = para[cur_id * plen + 8];
        double t10 = para[cur_id * plen + 37] / 5;
        double t11 = para[cur_id * plen + 38] / 5;
        double t12 = para[cur_id * plen + 39];
        double t13 = para[cur_id * plen + 40];
        double t14 = para[cur_id * plen + 41];
        double t15 = para[cur_id * plen + 42];
        double t16 = para[cur_id * plen + 43] / 6;
        double t17 = para[cur_id * plen + 44] / 6;
        double t18 = para[cur_id * plen + 45];
        double t19 = para[cur_id * plen + 46];
        double t20 = para[cur_id * plen + 21];
        double t21 = para[cur_id * plen + 47];
        double t22 = para[cur_id * plen + 48];
        if (row < height && col < width)
        {
            int now_id = cur_id * VOL + row * width + col;
            double a1 = t1 *  a[now_id].x - t2 *  a[now_id].y;

            double a2 = t3 *  b[now_id].x - t4 *  b[now_id].y;
            double a3 = t5 *  c[now_id].x - t6 *  c[now_id].y;
            double a4 = t7 *  d[now_id].x - t8 *  d[now_id].y;
            double a5 = t9 *  e[now_id].x / 4;
            double a6 = t10 * f[now_id].x - t11 * f[now_id].y;
            double a7 = t12 * g[now_id].x - t13 * g[now_id].y;
            double a8 = t14 * h[now_id].x - t15 * h[now_id].y;
            double a9 = t16 * i[now_id].x - t17 * i[now_id].y;
            double a10 = t18* j[now_id].x - t19 * j[now_id].y;
            double a11 = t20* k[now_id].x / 6;
            double a12 = t21* l[now_id].x - t22 * l[now_id].y;
            result[now_id] = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12;
        }
    }
}

__global__ void getW(double* result, double* kai, cufftDoubleComplex* a,
                     cufftDoubleComplex* b, double* para,int plen, int width,
                     int height, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int cur_id = id / VOL;
	int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if (cur_id < batches)
    {
        double t1 = para[cur_id * plen + 26];
        double t2 = para[cur_id * plen + 27];
        double t3 = para[cur_id * plen + 1];
        if (row < height && col < width)
        {
            int now_id = cur_id * VOL + row * width + col;
            double a1 = t1 * a[now_id].x - t2 * a[now_id].y;
            double a2 = t3 * b[now_id].x;
            result[now_id] = a1 / 2 + a2 / 2  + kai[now_id];
        }
    }
}

__device__ void tempgradient_u1(double* result, double la, double* para,
                            int plen, double* uu, double* vv,
                            int width, int height, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int cur_id = id / VOL;
	int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if(cur_id < batches)
    {
        double A1x = para[cur_id * plen + 26];
        double A1y = para[cur_id * plen + 27];
        double A2x = para[cur_id * plen + 28];
        double A2y = para[cur_id * plen + 29];
        double B2x = para[cur_id * plen + 30];
        double B2y = para[cur_id * plen + 31];
        double A3x = para[cur_id * plen + 32];
        double A3y = para[cur_id * plen + 33];
        double S3x = para[cur_id * plen + 34];
        double S3y = para[cur_id * plen + 35];
        double A4x = para[cur_id * plen + 37];
        double A4y = para[cur_id * plen + 38];
        double C3  = para[cur_id * plen + 8 ];
        if(row < height && col < width)
        {
            double u = uu[row * width + col];
            double v = vv[row * width + col];

            result[cur_id * VOL + row * width + col]=la*(A1x*u+A1y*v)
                +A2x*pow(la, 2)/3*(3*pow(u, 2)-3*pow(v, 2))+A2y*pow(la, 2)/3*6*u*v
                +B2x*(3*pow(u, 2)+pow(v, 2))*pow(la, 2) - B2y*pow(la, 2)*2*u*v
                +1.0/4*pow(la, 2)*la * (A3x*(4*pow(u, 2)*u-12*u*pow(v, 2))+A3y*(12*pow(u, 2)*v-4*pow(v, 2)*v))
                +pow(la, 2)*la* (2*u*(S3x*(pow(u, 2)-pow(v, 2))-S3y*2*u*v) +(pow(u, 2)+pow(v, 2))*(2*S3x*u-2*S3y*v))
                +C3*pow(la, 2)*la*(pow(u, 2)+pow(v, 2))*u
                +1.0/5*pow(la, 4) * (A4x*(5*pow(u, 4)+5*pow(v, 4)-30*pow(u, 2)*pow(v, 2))+A4y*(20*pow(u, 3)*v-20*u*pow(v, 3)));
        }
    }
}

__device__ void tempgradient_u2(double* result, double la, double* para,
    int plen, double* uu, double* vv,
    int width, int height, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
	int cur_id = id / VOL;
	int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if(cur_id < batches)
    {
        double C1  = para[cur_id * plen + 1];

        double D4x = para[cur_id * plen + 39];
        double D4y = para[cur_id * plen + 40];
        double B4x = para[cur_id * plen + 41];
        double B4y = para[cur_id * plen + 42];
        double A5x = para[cur_id * plen + 43];
        double A5y = para[cur_id * plen + 44];
        double S5x = para[cur_id * plen + 45];
        double S5y = para[cur_id * plen + 46];
        double C5  = para[cur_id * plen + 21];
        double D5x = para[cur_id * plen + 47];
        double D5y = para[cur_id * plen + 48];

        if(row < height && col < width)
        {
            double u = uu[row * width + col];
            double v = vv[row * width + col];

            result[cur_id * VOL + row * width + col]=result[cur_id * VOL + row * width + col]
                +pow(la, 4) * ( (2*u) * (D4x*(pow(u, 3)-3*u*pow(v, 2))-D4y*(3*pow(u, 2)*v-pow(v, 3))) + (pow(u, 2)+pow(v, 2))*(D4x*(3*pow(u, 2)-3*pow(v, 2))-D4y*6*u*v) )
                +pow(la, 4) * (  2*(pow(u, 2)+pow(v, 2))*2*u *(B4x*u-B4y*v) + pow((u*u+v*v),  2)*B4x)
                +1.0/6*pow(la, 5) * (A5x*(4*pow(u, 3)-14*v)*(pow(u, 2)-pow(v, 2))+A5x*(pow(u, 4)+pow(v, 4)-14*u*v)*2*u - A5y*2*v*(3*pow(u, 4)+3*pow(v, 4)-10*pow(u, 2)*pow(v, 2))-A5y*2*v*u*(12*pow(u, 3)-20*u*pow(v, 2)) )
                +pow(la, 5) * (2*(pow(u, 2)+pow(v, 2))*2*u*(S5x*(pow(u, 2)-pow(v, 2))-2*u*v*S5y) + pow((u*u+v*v),  2)*(2*S5x*u-2*v*S5y) )
                +1.0/6*pow(la, 5) * C5 * 3*pow((u*u+v*v),  2)*2*u
                +pow(la, 5) * (2*u*(D5x*(pow(u, 4)+pow(v, 4)-6*pow(u, 2)*pow(v, 2))-D5y*4*u*v*(pow(u, 2)-pow(v, 2))) + (pow(u, 2)+pow(v, 2))*(D5x*(4*pow(u, 3)-12*u*pow(v, 2))-D5y*(12*pow(u, 2)*v-4*pow(v, 3))) )
                + C1 * la * u;
        }
    }
}



__global__ void gradient_u(double* result, double la, double* para,
    int plen, double* uu, double* vv,
    int width, int height, int VOL, int batches)
{
    tempgradient_u1(result,  la,  para, plen, uu, vv, width,  height,  VOL,  batches);
    tempgradient_u2(result,  la,  para, plen, uu, vv, width,  height,  VOL,  batches);
}


__device__ void tempgradient_v1(double* result, double la, double* para,
    int plen, double* uu, double* vv,
    int width, int height,int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
    );
    int cur_id = id / VOL;
    int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if(cur_id < batches)
    {
        double A1x = para[cur_id * plen + 26];
        double A1y = para[cur_id * plen + 27];
        //double A2x = para[cur_id * plen + 28];
        //double A2y = para[cur_id * plen + 29];
        double B2x = para[cur_id * plen + 30];
        double B2y = para[cur_id * plen + 31];
        double A3x = para[cur_id * plen + 32];
        double A3y = para[cur_id * plen + 33];
        double S3x = para[cur_id * plen + 34];
        double S3y = para[cur_id * plen + 35];
        double C3  = para[cur_id * plen + 8 ];
        double A4x = para[cur_id * plen + 37];
        double A4y = para[cur_id * plen + 38];
        if(row < height && col < width)
        {
            double u = uu[row * width + col];
            double v = vv[row * width + col];

            result[cur_id * VOL + row * width + col] = la * (A1y * u - A1x * v)
            - A3x * pow(la, 2) / 3 * 6 * u * v
            + A3y * pow(la, 2) / 3 * (-3 * pow(v, 2) + 3 * pow(u, 2))
            + B2x * pow(la, 2) * 2 * u * v
            - B2y * pow(la, 2) * (pow(u, 2) + 3 * pow(v, 2))
            + 1.0 / 4 * pow(la, 3) * (A3x * (4 * pow(v, 3) - 12 * pow(u, 2) * v)
            + A3y * (4 * pow(u, 3) - 12 * u * pow(v, 2)))
            + pow(la, 3) * (2 * v * (S3x * (pow(u, 2) - pow(v, 2)) - S3y * 2 * u * v) + (pow(u, 2) + pow(v, 2)) * (-2 * S3x * v - 2 * S3y * u))
            + C3 * pow(la, 3) * (pow(u, 2) + pow(v, 2)) * v
            + 1.0 / 5 * pow(la, 4) * (A4x * (20 * u * pow(v, 3) - 20 * pow(u, 3) * v) + A4y * (5 * pow(v, 4) + 5 * pow(u, 4) - 30 * pow(u, 2) * pow(v, 2)));
        }
    }
}
__device__ void tempgradient_v2(double* result, double la, double* para,
    int plen, double* uu, double* vv,
    int width, int height,int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
    );
    int cur_id = id / VOL;
    int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if(cur_id < batches)
    {
        double C1  = para[cur_id * plen + 1];
        double D4x = para[cur_id * plen + 39];
        double D4y = para[cur_id * plen + 40];
        double B4x = para[cur_id * plen + 41];
        double B4y = para[cur_id * plen + 42];
        double A5x = para[cur_id * plen + 43];
        double A5y = para[cur_id * plen + 44];
        double S5x = para[cur_id * plen + 45];
        double S5y = para[cur_id * plen + 46];
        double C5  = para[cur_id * plen + 21];
        double D5x = para[cur_id * plen + 47];
        double D5y = para[cur_id * plen + 48];

        if(row < height && col < width)
        {
            double u = uu[row * width + col];
            double v = vv[row * width + col];
            result[cur_id * VOL + row * width + col] = result[cur_id * VOL + row * width + col]
            + pow(la, 4) * ((2 * v) * (D4x * (pow(u, 3) - 3 * u * pow(v, 2)) - D4y * (3 * pow(u, 2) * v - pow(v, 3))) + (pow(u, 2) + pow(v, 2)) * (-D4x * 6 * u * v - D4y * (3 * pow(u, 2) - 3 * pow(v, 2))))
            + pow(la, 4) * (2 * (pow(u, 2) + pow(v, 2)) * 2 * v *(B4x * u - B4y * v) - pow((pow(u, 2) + pow(v, 2)), 2) * B4y)
            + 1.0 / 6 * pow(la, 5) * (A5x * (4 * pow(v, 3) - 14 * u) * (pow(u, 2) - pow(v, 2))
            + A5x * (pow(u, 4) + pow(v, 4) - 14 * u * v) * (-2 * v) - A5y * 2 * u * (3 * pow(u, 4) + 3 * pow(v, 4) - 10 * pow(u, 2) * pow(v, 2)) - A5y * 2 * v * u * (12 * pow(v, 3) - 20 * pow(u, 2) * v))
            + pow(la, 5) * (2 * (pow(u, 2) + pow(v, 2)) * 2 * v * (S5x * (pow(u, 2) - pow(v, 2)) -2 * u * v * S5y) + pow((pow(u, 2) + pow(v, 2)), 2) * (-2 * S5x * v - 2 * u * S5y))
            + 1.0 / 6 * pow(la, 5) * C5 * 3 * pow((pow(u, 2) + pow(v, 2)), 2) * 2 * v
            + pow(la, 5) * (2 * v * (D5x * (pow(u, 4) + pow(v, 4) - 6 * pow(u, 2) * pow(v, 2)) - D5y * 4 * u * v * (pow(u, 2) - pow(v, 2))) + (pow(u, 2) + pow(v, 2)) * (D5x * (4 * pow(v, 3) - 12 * pow(u, 2) * v) - D5y * (4 * pow(u, 3) - 12 * u * pow(v, 2))))
            + C1 * la * v;
        }
    }
}

__global__ void gradient_v(double* result, double la, double* para,
    int plen, double* uu, double* vv,
    int width, int height,int VOL, int batches)
{
    tempgradient_v1( result,  la,  para, plen, uu,  vv, width,  height, VOL,  batches);
    tempgradient_v2( result,  la,  para, plen, uu,  vv, width,  height, VOL,  batches);
}

__global__ void getESCoh(double* result, double al, double la, double* gradient_kai_u, double* gradient_kai_v,
                        int VOL, int batches, double pi)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y)
    + (threadIdx.y * blockDim.x + threadIdx.x);
    if(id < batches * VOL )
    {
        double u = pow(gradient_kai_u[id], 2);
        double v = pow(gradient_kai_v[id], 2);
        result[id] = expf(-1 * al * al * pi * pi / la / la * (u + v));
    }
}

__global__ void getT(cufftDoubleComplex* result, double* w, double la,
                     int VOL, int batches, double pi)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y)
            + (threadIdx.y * blockDim.x + threadIdx.x);
    if (id < batches * VOL)
    {
        //double real = 0;
        double img =  2 * pi * w[id] / la;
        result[id].x = expf(0) * cos(-1 * img);
        result[id].y = expf(0) * sin(-1 * img);
    }
}

__global__ void  getTlittleAll(cufftDoubleComplex* result, double* E_s_coh, cufftDoubleComplex* t,
                            int mt, double pi, double de, double* u, double* v,
                            double la, int height, int width, int VOL, int batches)
{
    int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y)
            + (threadIdx.y * blockDim.x + threadIdx.x);
    int cur_id = id / VOL;
    int row = (id % VOL) / width;
    int col = (id % VOL) % width;
    if (cur_id < batches)
    {
        if (row<height && col<width)
        {
            int now_id = cur_id * VOL + row * width + col;
            for(int m=0; m < mt; m++)
            {
                double t1 = -1 * pi * (m+1-(mt+1)*1.0 / 2) * de
                            * (pow(u[row * width + col], 2) + pow(v[row * width + col], 2)) * la;
                double a = expf(0) * cos(t1);
                double b = expf(0) * sin(t1);
                double c = t[now_id].x;
                double d = t[now_id].y;
                double e = E_s_coh[now_id];
                result[m * VOL * batches + now_id].x = (a*c-b*d) * e;
                result[m * VOL * batches + now_id].y = (b*c+a*d) * e;
            }
        }
    }
}



__global__ void myIfftShift(cufftDoubleComplex* P_PitchedPtr, int pitch, int P_Height, int P_Width, int P_Slices, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
	int d_id = id / (VOL / 2);
	int b_id = id % (VOL / 2) / (VOL / 4);
	int e_id = id % (VOL / 2) % (VOL / 4);

	int row = e_id / (P_Width / 2);
	int col = e_id % (P_Width / 2);

	int dest_row, dest_col;

	if (d_id < P_Slices)
	{
        //printf("%d\\n", d_id);
		if (b_id == 0)
		{
			dest_row = row + (P_Height / 2);
			dest_col = col + (P_Width / 2);
		}
		else if (b_id == 1)
		{
			col = col + (P_Width / 2);
			dest_row = row + (P_Height / 2);
			dest_col = col - (P_Width / 2);
		}
		double e1_real = P_PitchedPtr[d_id * pitch+row * P_Width + col].x;
		double e1_img = P_PitchedPtr[d_id * pitch+row * P_Width + col].y;

		double e2_real = P_PitchedPtr[d_id * pitch+dest_row * P_Width + dest_col].x;
		double e2_img = P_PitchedPtr[d_id * pitch+dest_row * P_Width + dest_col].y;

		P_PitchedPtr[d_id * pitch+row * P_Width + col].x = e2_real;
		P_PitchedPtr[d_id * pitch+row * P_Width + col].y = e2_img;

		P_PitchedPtr[d_id * pitch+dest_row * P_Width + dest_col].x = e1_real;
		P_PitchedPtr[d_id * pitch+dest_row * P_Width + dest_col].y = e1_img;
	}
}


__global__ void accumulate_result(double* result, double* gfsf, double* result1,
                                 int gfsf_num, int batches, int height, int width, int VOL)
{
	int id = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x
	);
    int cur_id = id / VOL; // 0-batches
    int row = (id % VOL) / width;
    int col = (id % VOL) % width;
	if (row < height && col < width)
	{
	    int now_id = cur_id * VOL + row * width + col;
		double temp = 0;
		for (int j = 0; j < gfsf_num; j++)
		{
			temp += result1[cur_id * gfsf_num * VOL + j * VOL + row * width + col] * gfsf[j];
		}
		result[now_id] = temp;
	}
}
""")
        # ctx.push()
        self.createNewPP2 = mod.get_function("createNewPP2")
        self.getKai = mod.get_function("getKai")
        self.getW = mod.get_function("getW")
        self.gradient_u = mod.get_function("gradient_u")
        self.gradient_v = mod.get_function("gradient_v")
        self.getESCoh = mod.get_function('getESCoh')
        self.getT = mod.get_function('getT')
        self.getTlittleAll = mod.get_function('getTlittleAll')
        self.myIfftShift = mod.get_function('myIfftShift')
        self.accumulate_result = mod.get_function('accumulate_result')

        # 每批次计算batches个
        global batches
        batches = bat
        count = (self.NUMS + batches - 1) // batches
        lenPP2 = len(self.pp2)

        mtotal = np.int32(self.para_part1['mtotal'])
        plan = cu_fft.Plan((N, N), np.complex128,
                           np.complex128, batches * mtotal)
        for _ in range(count):
            newPp2 = np.zeros((batches, lenPP2)).astype(np.float64)
            # newPp2 = np.repeat(self.pp2, batches, 1).astype(np.float64)
            # 随机数
            ai = curandom.rand((batches, 21)).astype(np.float64)
            self.createNewPP2(
                cuda.InOut(newPp2), cuda.In(self.pp2), cuda.In(self.mrr_list), ai,
                np.int32(lenPP2), np.int32(batches), np.float64(np.pi),
                block=block, grid=self.setThreadGrid(batches, 1)
            )

            mytcc = self.WholeTCC2D_newTEM_forpar(self.para_part1, newPp2, batches)

            # myinten = np.zeros((batches, N, N))

            myware = np.repeat(
                self.mywave[np.newaxis, :, :], batches * mtotal, axis=0)
            myware_gpu = gpuarray.to_gpu(myware)
            result1 = myware_gpu.__mul__(mytcc)
            self.myIfftShift(
                result1, np.int32(N * N), np.int32(N), np.int32(N),
                np.int32(batches), np.int32(N * N),
                block=block, grid=self.setThreadGrid(batches, N * N // 2)
            )
            result2 = gpuarray.zeros_like(result1, np.complex128)

            cu_fft.ifft(result1, result2, plan, True)
            result3 = result2.__abs__()
            result4 = result3.__mul__(result3)
            gfsf = np.array([self.para_part1['gfsf']])
            result = gpuarray.zeros((batches, N, N), dtype=np.float64)

            self.accumulate_result(
                result, cuda.In(gfsf), result4,
                np.int32(len(gfsf)), np.int32(batches),
                np.int32(N), np.int32(N), np.int32(N * N),
                block=block, grid=self.setThreadGrid(batches, N * N)
            )
            myinten = result.get()
            q.put((myinten, newPp2, batches))
        q.put((0, 0, 0))
        ctx.pop()


def consumer(q, commpara, file_name, file_path, position):
    # global file_name
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    while True:
        myinten, newPp2, batches = q.get()
        if batches > 0:
            for b in range(batches):
                # myresul = myinten[b, 123:151, 144:188]
                myresul = myinten[b, position[0]:position[1], position[2]:position[3]]
                # myresul = myinten[b]
                myresul = (myresul - myresul.min()) / \
                          (myresul.max() - myresul.min())
                savimage = np.repeat(myresul[:, :, np.newaxis], 3, axis=2)
                plt.imsave(os.path.join(file_path, str(
                    file_name) + '.jpg'), savimage)
                para = newPp2[b]
                savpara = np.array([para[1], para[2], -para[3], para[4], -para[5], para[6], -para[7],
                                    para[8] / 1000, para[9] / 1000, -
                                    para[10], para[11] / 1000, -para[12],
                                    para[13] / 1000, -para[14],
                                    para[15] / 1000, -para[16],
                                    para[17] / 1000, -para[18],
                                    para[21] / 1000000,
                                    para[19] / 1000000, -para[20],
                                    para[24] / 1000000, -para[25],
                                    para[22] / 1000000, -para[23]
                                    ])
                savpara = np.hstack((commpara, savpara)).astype(np.float)
                savpara.tofile(os.path.join(
                    file_path, str(file_name) + ".para"))
                sys.stdout.write(
                    "\r" + ">> {0} image and para is writing".format(file_name))
                sys.stdout.flush()
                file_name += 1

        else:
            break


if __name__ == '__main__':
    matlab = scio.loadmat("pydata.mat")

    othervar = matlab['othervar'][0]
    process_num = int(othervar[0])
    batches = int(othervar[1])
    allnums = int(othervar[2])
    top, left = othervar[3], othervar[4]
    width, height = othervar[5], othervar[6]
    # startX, startY = map(int, str(matlab['othervar'][3][0]).split())
    # width, height = map(int, str(matlab['othervar'][4][0]).split())
    file_path = str(matlab['savePathName'][0])

    position = [left, left + width, top, top + height]

    # allnums = 50
    # batches = 1
    # process_num = 1
    # position = [0, 320, 0, 320]

    cuda.init()
    # 选择显卡id
    dev = cuda.Device(0)
    ctx = dev.make_context()
    cc = multiprocessing.cpu_count()
    if process_num > cc // 2:
        process_num = cc // 2 - 1
    while True:
        if allnums % process_num == 0:
            break
        process_num -= 1
    if process_num < 1:
        process_num = 1
    nums = allnums // process_num

    print("进程数量:", process_num)
    print("每个进程生成的数量:", nums)
    print("每个进程每批次生成的数量:", batches)
    print("保存路径:", file_path)

    para_part1 = list(matlab['para_part1'][0][0])
    ppl = []
    # gfsf = []
    for p in para_part1:
        # if len(p[0]) > 1:
        #     gfsf = p[0]
        #     continue
        ppl.append(p[0][0])
    ppl = np.array(ppl).astype(np.float64)
    pp1_name = ['sampling', 'mtotal', 'alafa',
                'lambda', 'gmax', 'gfsf',
                'delta_yita', 'tilt', 'phitilt',
                'tiltx', 'tilty']
    pp1 = dict(zip(pp1_name, ppl))
    # pp1['gfsf'] = gfsf

    para_part2 = matlab['para_part2'][0][0]
    pp2 = []
    for p in para_part2:
        pp2.append(p[0][0])
    pp2 = np.array(pp2).astype(np.float64)
    myrandrange = matlab['myrandrange'][0][0]
    mar = []
    for p in myrandrange:
        mar.append(p[0][0])
    mar = np.array(mar).astype(np.float64)
    U = matlab['U']
    V = matlab['V']
    commpara = matlab['commpara'][0]
    mywave = matlab['mywave']

    q_list = []
    for i in range(process_num):
        q = Manager().Queue()
        q_list.append(q)
    start_time = time.time()
    pool = Pool(process_num)
    for i in range(process_num):
        s = stem(pp1, pp2, mar, U, V, mywave, nums)
        pool.apply_async(s.producer, args=(q_list[i], batches,))
    poolc = Pool(process_num)
    for i in range(process_num):
        poolc.apply_async(consumer, args=(
            q_list[i], commpara, i * nums, file_path, position,))

    pool.close()
    poolc.close()
    pool.join()
    poolc.join()

    '''
    q = Queue()
    s = stem(pp1, pp2, mar, U, V, A1_phi,
            A2_phi, B2_phi, A3_phi, S3_phi, B4_phi, mywave, nums)
    p = Process(target=s.producer, args=(q, batches, 0))
    c = Process(target=consumer, args=(q, commpara, 0, file_path,))
    p.start()
    c.start()
    p.join()
    c.terminate()
    '''
    ctx.pop()
    end_time = time.time()
    print("\n%d images finished" % (allnums))
    print("cost time:", end_time - start_time)
