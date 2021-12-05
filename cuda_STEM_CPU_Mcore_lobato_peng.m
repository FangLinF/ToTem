%cuda的核心计算，多cpu的算法
function [aamyresul, aamidresul]=cuda_STEM_CPU_Mcore_lobato_peng(PARAMETER, sigma, ...%包含了哪些原子，每层原子多少个；参数parameter算potential的;相互作用系数sigama，算透过函数需要的
    ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
    ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
    series_n, series_n_i, ...  %原子排列次序
    gx_green, gy_green, s2, APERTURE,...   %绿色区域计算potential时候，gx和gy是用来算原子位置相关的倒格矢，s2是算peng potential的,APERTURE是防止2/3光阑外发生wrap，详见笔记
    probe, gfsf, psf_fft, ...%probe的形式，还有几个矩阵的比例分配，传播矩阵
    green_Ncol, green_Nrow, ...%绿色的尺寸，绿色倒空间的sx和sy，
    probe_ingreenNrow, probe_ingreenNcol, ...
    width_red, hight_red, probestep, ...  %选取的相对位置的左上角x和y坐标，以及右下角x和y坐标
    APER,paraflag, midslice);  %最后增加的光阑函数；采用lobato还是peng的参数来计算

%计算每层的势场
potential=GetPotential4AllSlice_multicore_lotabo_peng(green_Ncol, green_Nrow,... 
    ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
    ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
    series_n, series_n_i, ...  %原子排列次序
    s2, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。

%传播
[probesx, probesy]=size(probe);


for i=0:width_red-1;
     tempresul(i+1).matrix = zeros(hight_red,length(APER(1,1,:)));  %清零
     midtempresul(i+1).matrix= zeros(hight_red,length(APER(1,1,:)),length(midslice));  %清零
end
% parfor i=0:width_red-1;
%     i
%     for j=0:hight_red-1;
%         for nn=1:length(gfsf);
%             mywave=probe(:,:,nn);  %读入波函数
%             for k=1:length(potential(1,1,:));
%                 mywave =  mywave.*potential(j*probestep+probe_ingreenNrow:j*probestep+probe_ingreenNrow+probesx-1, ...
%                    i*probestep+probe_ingreenNcol:i*probestep+probe_ingreenNcol+probesx-1,k);
%                    mywave=fft2(mywave);   %注意，节省掉一个fftshift
%                    mywave=mywave.*psf_fft;
%                    mywave=ifft2(mywave);
%             end
%             mywave=fftshift(fft2(mywave));
%             for kk=1:length(APER(1,1,:))    
%                 tempresul(i+1).matrix(j+1,kk)=tempresul(i+1).matrix(j+1,kk)+gfsf(nn).*sum(sum(abs(mywave.*APER(:,:,kk)).^2));%低通、带通，或高通函数
%                 %myresul(j+1,i+1,kk)=myresul(j+1,i+1,kk)+gfsf(nn).*sum(sum(abs(mywave.*APER(:,:,kk)).^2));  %低通、带通，或高通函数
%             end
%         end
%     end
% end

for kk=1:length(APER(1,1,:))    
       APER(:,:,kk)=ifftshift(APER(:,:,kk));  %把光阑的象限重新分配一下
end 

mypar=parpool;
parfor i=0:width_red-1;
    i
    for j=0:hight_red-1;
        for nn=1:length(gfsf);
            mywave=probe(:,:,nn);  %读入波函数
            for k=1:length(potential(1,1,:));
                mywave =  mywave.*potential(j*probestep+probe_ingreenNrow:j*probestep+probe_ingreenNrow+probesx-1, ...
                   i*probestep+probe_ingreenNcol:i*probestep+probe_ingreenNcol+probesx-1,k);
                   mywave=fft2(mywave);   %注意，节省掉一个fftshift
                   mywave=mywave.*psf_fft;
                   
                   if ~isempty(find(midslice==k, 1))
                       oo=find(midslice==k, 1);
                       for kk=1:length(APER(1,1,:))    
                          midtempresul(i+1).matrix(j+1,kk,oo)=midtempresul(i+1).matrix(j+1,kk,oo)+gfsf(nn).*sum(sum(abs(mywave.*APER(:,:,kk)).^2));%低通、带通，或高通函数
                       end 
                   end
                   mywave=ifft2(mywave);
            end
            mywave=fft2(mywave);
            for kk=1:length(APER(1,1,:))    
                 tempresul(i+1).matrix(j+1,kk)=tempresul(i+1).matrix(j+1,kk)+gfsf(nn).*sum(sum(abs(mywave.*APER(:,:,kk)).^2));%低通、带通，或高通函数
            end 
        end
    end
end
delete(mypar)

myresul=zeros(width_red, hight_red, length(APER(1,1,:)));
midresul=zeros(width_red, hight_red, length(APER(1,1,:)),length(midslice));
for i=0:width_red-1;
    for j=1:length(midslice)
        midresul(i+1,:,:,:) = midtempresul(i+1).matrix(:,:,:);
    end
    myresul(i+1,:,:) = tempresul(i+1).matrix(:,:);
end

clear tempresul
clear midtempresul
aamidresul = zeros(hight_red*width_red*length(APER(1,1,:)),length(midslice));
aamyresul = zeros(hight_red*width_red,length(APER(1,1,:)));
aamidresul(:) = midresul(:);
aamyresul(:) = myresul(:);
clear tempresul
% for kk=1:length(APER(1,1,:)) 
%     figure;imshow(myresul(:,:,kk)/handles.width_red/handles.height_red,[]);colorbar
% end
return