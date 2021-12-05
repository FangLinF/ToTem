function multislice_unwrapphase(mywave, potential, psf_fft, sampling, lambda, eachthick)
[e,f]=find(fftshift(abs(ifft2(psf_fft))) == max(max(fftshift(abs(ifft2(psf_fft))))))
%找到psf_fft的中心位置。
psf_realspace = fftshift(ifft2(psf_fft));
sizexy = 10;  %规定一下psf的尺寸
smallpsf = psf_realspace(e-sizexy:e+sizexy, f-sizexy:f+sizexy);  %取一小部分进行卷积 
%上述的smallpsf只能作为参考，因为wrap太多次了，很难解的
[sx,sy]=size(psf_fft); 
psf=zeros(sx, sy);  
%用这个psf来乘以吧。首先做一个和上述时域psf相同尺寸的矩阵，再用时域公式区计算时域的复数矩阵，看看能不能对得上。
[xx,yy]=meshgrid([1:sx]-e, [1:sy]-f);
xx=xx*sampling;
yy=yy*sampling;
%做出x和y的矩阵,R=(x,y),以下根据6.90式子写出psf在时域的表示
psf=(1/sqrt(-1)/lambda/eachthick)*exp(sqrt(-1)*pi/lambda/eachthick*(xx.^2+yy.^2));
%这么算出来的psf，多半和psf_realspace不一样。因为在主干程序里面有这句psf_fft(find(p2>1/(16*handles.sampling*handles.sampling)))=0;
%所以比较的时候，考虑把上面那句不要了，然后进行比较
%其实在图像域直接写psf，不需要考虑滤波。因为在傅里叶做法里面，不加这个光阑限制的话，左边和右边，上边和下边的，会出现混叠效果。但是如果是图像域直接写卷积的话，不会有这样的干扰，外围填0了的。
%psf对应出来后，用相位的显示来写！
%--------------------得到psfA和psfP
%这里省略了一句，需要对后才能写的吧。要保证和psfsmall的复数值一样，而相位存在反转


mywaveA = abs(mywave);
mywaveP = angle(mywave);
mywave = mywaveA.*exp(sqrt(-1)*mywaveP);  %在计算势场的时候，exp(sqrt(-1)*V)是散射势；
%靠近原子中心，电子云大，相位落后。所以这里定义，波函数中的相位部分的正的数值，代表相位落后
for k=1:length(potential(1,1,:));
        mywave =  mywave.*potential(:,:,k);
        
        potentialA = abs(potential(:,:,k));
        potentialP = unwrap( angle(potential(:,:,k)) );  
        %把波函数的相位正确的解开。其实对应了GetPotential4AllSlice_multicore_lotabo_peng_corr
        %246行的V部分。注意，V只能是正值的。
        mywaveA = getmultiple_AP(mywaveA, mywaveP, potentialA, potentialP);   %这个函数你来写
        %把mywaveA和mywaveP来作为乘法的输入信息，表示一个复数波函数的振幅和相位。波函数和势场点乘
        
        [mywaveA, mywaveP] = getconv_AP(mywaveA, mywaveP, psfA, psfP);  %这里除了复数的乘法，还有复数的加法问题要解决
%上式是在图像域里面处理的，所以不存在傅里叶变换了。
%         mywave=fft2(mywave);   %注意，节省掉一个fftshift
%         mywave=mywave.*psf_fft;
%         mywave=ifft2(mywave);
end