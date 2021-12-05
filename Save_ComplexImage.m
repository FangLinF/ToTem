function Save_ComplexImage(myI,filename)
[Ncol,Nrow]=size(myI);
%保存复数图像的方法。
temp=myI.';  %存储复数矩阵的方法  需要图像做转置
tempall(1:2:2*Nrow-1,:)=real(temp);
tempall(2:2:2*Nrow,:)=imag(temp);
fid=fopen(filename,'w');
fwrite(fid,tempall,'float');
fclose(fid);