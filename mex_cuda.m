if exist('myresul.mat','file')
    delete('myresul.mat');
end
if exist('allresul.mat','file')
    delete('allresul.mat');
end
if exist('paraslicing.txt','file')
     delete('paraslicing.txt');
end


system("nvcc -c cuda_source_lobato.cu");
system("nvcc -c cuda_source_peng_and_corr.cu");
system("nvcc -c cuda_source_lobato_potential.cu");
system("nvcc -c cuda_source_peng_and_corr_potential.cu");
cuda_name = cuda_v;
if ispc
    CUDA_PATH = getenv('CUDA_PATH');
    cuda_cmd = sprintf('%s"%s%s" %s%s%s"',"mex  waveCuda.cpp cuda_source_lobato.obj -L",CUDA_PATH,"\lib\x64", "-lcudart -lcufft -output waveCudalobato_",cuda_name,".mexw64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s"%s%s" %s%s%s"',"mex  waveCuda.cpp cuda_source_peng_and_corr.obj -L",CUDA_PATH,"\lib\x64", "-lcudart -lcufft -output waveCudaPengAndCorr_",cuda_name,".mexw64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s"%s%s" %s%s%s"',"mex  waveCuda.cpp cuda_source_lobato_potential.obj -L",CUDA_PATH,"\lib\x64", "-lcudart -lcufft -output waveCudalobatoP_",cuda_name,".mexw64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s"%s%s" %s%s%s"',"mex  waveCuda.cpp cuda_source_peng_and_corr_potential.obj -L",CUDA_PATH,"\lib\x64", "-lcudart -lcufft -output waveCudaPengAndCorrP_",cuda_name,".mexw64");
    eval(cuda_cmd);
    %编译命令
    mcc -m ToTEM_submit.m -a *.mexw64;
end
if isunix
    cuda_cmd = sprintf('%s%s%s',"mex waveCuda.cpp cuda_source_lobato.o -lcudart -lcufft -L/usr/local/cuda/lib64 -output waveCudalobato_",cuda_name,".mexa64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s%s%s',"mex waveCuda.cpp cuda_source_peng_and_corr.o -lcudart -lcufft -L/usr/local/cuda/lib64 -output waveCudaPengAndCorr_",cuda_name,".mexa64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s%s%s',"mex waveCuda.cpp cuda_source_lobato_potential.o -lcudart -lcufft -L/usr/local/cuda/lib64 -output waveCudalobatoP_",cuda_name,".mexa64");
    eval(cuda_cmd);
    cuda_cmd = sprintf('%s%s%s',"mex waveCuda.cpp cuda_source_peng_and_corr_potential.o -lcudart -lcufft -L/usr/local/cuda/lib64 -output waveCudaPengAndCorrP_",cuda_name,".mexa64");
    eval(cuda_cmd);
    %编译命令
    mcc -m ToTEM_submit.m -a *.mexa64;
end
