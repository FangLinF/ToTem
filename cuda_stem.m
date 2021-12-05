
slicenum=max(length(series_n), length(series_n_i));
if isempty(series_n)
    series_n=zeros(1,slicenum);
else
    if length(series_n)<slicenum
        temp=zeros(1,slicenum);
        temp(1:length(series_n))=series_n;
        series_n=temp;
    end
        
end
if isempty(series_n_i)
    series_n_i=zeros(1,slicenum);
else
    if length(series_n_i)<slicenum
        temp=zeros(1,slicenum);
        temp(1:length(series_n_i))=series_n_i;
        series_n_i=temp;
    end
end
my_real_w = permute(real(handles.probe),[2,1,3]);
my_real_w = reshape(my_real_w,256,[]);
my_img_w = permute(imag(handles.probe),[2,1,3]);
my_img_w = reshape(my_img_w,256,[]);
[Height,Width,w_num] = size(handles.probe);
[aper_height,aper_width,aper_num] = size(APER);
Depth = handles.width_red * handles.hight_red;
aperTrue = APERTURE';
aper = permute(APER,[2,1,3]);
aper = reshape(aper,256,[]);


APER2=zeros(aper_height, aper_width, aper_num);
for i=1:aper_num;
    APER2(:,:,i) = ifftshift(APER(:,:,i));
end
aper2 = permute(APER2,[2,1,3]);
aper2 = reshape(aper2,256,[]);

atom_slice = series_n;
atom_slice_i = series_n_i;
proj_coff_mat = ele_n';
proj_coff_mat_i = ele_n_i';
Kx = gx_green';
Ky = gy_green';
S_2 = s2_green';
my_real_const = real(psf_fft)';
my_img_const = imag(psf_fft)';
gfsf = handles.gfsf;
P_Height = handles.green_Ncol;
P_Width = handles.green_Nrow;
beginRow = handles.probe_ingreenNrow;
beginCol = handles.probe_ingreenNcol;
width_red = handles.width_red;
parameter = PARAMETER;
Slices = length(series_n);
Step = handles.probestep;
absorp_flag = 0;
absorp_i_flag = 0;
if ~isempty(absorp_n) % ���absorp_n
    absorp_flag = 1;
end
if ~isempty(absorp_n_i) % ���absorp_n_i
    absorp_i_flag = 1;
end
mid_ceng = uint32(handles.mid_slice_num);
ceng_len = length(mid_ceng(:));
absorp_n2 = absorp_n';
absorp_n_i2 = absorp_n_i';

% corr
seriesn_corr = series_n_corr;
seriesn_i_corr = series_n_i_corr;
elen_corr = ele_n_corr';
elen_i_corr = ele_n_i_corr';
corrinfo_matrix = 0;
cimlen = 0;


cuda_ver = cuda_v;
%peng's scattering factor
if paraflag == 'p'
    waveCuda = ['waveCudaPengAndCorr_',cuda_ver];
end
if paraflag == 'n'
    corrinfo_matrix = permute(corr_info_matrix,[2,1,3]);
    [a, b, cimlen] = size(corr_info_matrix);
    waveCuda = ['waveCudaPengAndCorr_',cuda_ver];
end
% lobato's scattering factor
if paraflag == 'l'
   waveCuda = ['waveCudalobato_',cuda_ver];  
end


run_cmd = sprintf("%s(%s)",waveCuda,"my_real_w,my_img_w,aper,atom_slice,absorp_n2,atom_slice_i,absorp_n_i2,proj_coff_mat,proj_coff_mat_i,Kx,Ky,S_2,my_real_const,my_img_const,gfsf,Height,Width,Depth,P_Height,P_Width,Slices,Step,beginRow,beginCol,width_red,sigma,w_num,aper_num,parameter,aperTrue,absorp_flag,absorp_i_flag,mid_ceng,ceng_len,aper2, seriesn_corr, seriesn_i_corr, elen_corr, elen_i_corr, corrinfo_matrix, cimlen");
% [myresul,mid_ceng_mat] = waveCuda(my_real_w,my_img_w,aper,atom_slice,absorp_n',atom_slice_i, ...
%                absorp_n_i',proj_coff_mat,proj_coff_mat_i,Kx,Ky,S_2,my_real_const,...
%                my_img_const,gfsf,Height,Width,Depth,P_Height,P_Width,Slices,Step,...
%                beginRow,beginCol,width_red,sigma,w_num,aper_num,parameter,aperTrue,....
%                absorp_flag,absorp_i_flag,mid_ceng,ceng_len,aper2);
           %mid_ceng_mat aper,ceng,matrix size
% 如果是算结果，potentialx 和 potentialy是空值
% 如果是算势场，myresul 和 mid_ceng_mat 是空值
[myresul, mid_ceng_mat, potentialx, potentialy] = eval(run_cmd);
eval(sprintf("%s %s","clear",waveCuda));



