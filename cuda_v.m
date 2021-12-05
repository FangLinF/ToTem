function cuda_name=cuda_v
    [~, result] = system("nvcc -V");
    cuda_version = extractAfter(result,"release ");
    cuda_list = strsplit(cuda_version,'.');
    cuda_name = cuda_list{1};
end
