function [all_nuclear,all_nuc_ion]  = getpotentialpara(x); 
[PengPara, PengIon]=ABSFpara;  %提出所有的参数

%x包含所有原子的信息，第1到3列是原子坐标；第4列是原子的质子数；第5列是价态；
%第6列是离子性；第7列是debye；第8列是占据率
          % %需要认真得到x，王晨做

all_nuclear=[]; all_nuc_ion=[];  %结果先赋空
all_ion=[];

%case 1，只有原子性或只有离子性
b=find(x(:,6)==0);  %找出只有原子性的原子
if ~isempty(b)
    all_nuclear(:,1) = x(b,1);  %获得质子数
    all_nuclear(:,2:4) = x(b,2:4); %获得坐标
    all_nuclear(:,5) = x(b,7); %DB
    all_nuclear(:,6) = x(b,8).*(1-x(b,6)); %占据率*离子性
    all_nuclear(:,7:16) =  [PengPara(all_nuclear(:,1), 1), PengPara(all_nuclear(:,1), 6), PengPara(all_nuclear(:,1), 2), ...
                        PengPara(all_nuclear(:,1), 7), PengPara(all_nuclear(:,1), 3),...
                        PengPara(all_nuclear(:,1), 8), PengPara(all_nuclear(:,1), 4), PengPara(all_nuclear(:,1), 9), ...
                        PengPara(all_nuclear(:,1), 5), PengPara(all_nuclear(:,1), 10)]; %每个原子的散射参数
end

b=find(x(:,6)==1);  %找出只有离子性的原子
len=length(b);  %个数有几个离子性的
if ~isempty(b)
    all_ion(:,1) = x(b,1);  %获得质子数
    all_ion(:,2) = x(b,5);  %获得价态
    all_ion(:,3:5) = x(b,2:4); %获得坐标
    all_ion(:,6) = x(b,7); %DB
    all_ion(:,7) = x(b,8).*x(b,6); %占据率*离子性
    
    [tempjiatai, aa, xuhao] =unique(all_ion(:,1:2),'rows','stable');
    clear aa
    %把质子数和价态对应的参数，都赋值到新的参数里面
    for i=1:length(tempjiatai(:,1))
        k = find(tempjiatai(i,1)==PengIon(:,1) & tempjiatai(i,2)==PengIon(:,2));
        tempionpara(i, :)  = [PengIon(k,3), PengIon(k,8), PengIon(k,4), PengIon(k,9), PengIon(k,5), ...
                   PengIon(k,10), PengIon(k,6), PengIon(k,11), PengIon(k,7), PengIon(k,12)];
    end
    all_ion(:,8:17) = tempionpara(xuhao, :);   %把所有原子的离子性都赋值。从tempionpara中提取
end

%这里if的判断，得到了单纯离子性或原子性的规定矩阵。详见 程序说明.doc文件
if ~isempty(all_ion)   
    all_nuclear(end+1:end+len,:)=all_ion(:,2:end);
    clear all_ion
end

% case 2
%找出既有原子性又有离子性的原子。---------------------------------
b=find(x(:,6)>0 & x(:,6)<1);  %找出只有离子性的原子
if ~isempty(b)
    all_nuc_ion(:,1) = x(b,1);  %获得质子数
    all_nuc_ion(:,2) = x(b,5);  %获得价态
    all_nuc_ion(:,3:5) = x(b,2:4); %获得坐标
    all_nuc_ion(:,6) = x(b,7); %DB
    all_nuc_ion(:,7) = x(b,8); %占据率
    all_nuc_ion(:,8) = x(b,6); %离子性
    
    [tempjiatai, aa, xuhao] =unique(all_nuc_ion(:,1:2),'rows','stable');
    clear aa
    %把质子数和价态对应的参数，都赋值到新的参数里面
    for i=1:length(tempjiatai(:,1))
        k = find(tempjiatai(i,1)==PengIon(:,1) & tempjiatai(i,2)==PengIon(:,2));
        tempionpara(i, 1:10)  = [PengIon(k,3), PengIon(k,8), PengIon(k,4), PengIon(k,9), PengIon(k,5), ...
                   PengIon(k,10), PengIon(k,6), PengIon(k,11), PengIon(k,7), PengIon(k,12)];
    end
    all_nuc_ion(:,9:18) = tempionpara(xuhao, :);   %把所有原子的离子性都赋值。从tempionpara中提取
    
    all_nuc_ion(:,19) = 1- x(b,6); %原子性
    all_nuc_ion(:,20:29) = [PengPara(all_nuc_ion(:,1), 1), PengPara(all_nuc_ion(:,1), 6), PengPara(all_nuc_ion(:,1), 2), ...
                        PengPara(all_nuc_ion(:,1), 7), PengPara(all_nuc_ion(:,1), 3),...
                        PengPara(all_nuc_ion(:,1), 8), PengPara(all_nuc_ion(:,1), 4), PengPara(all_nuc_ion(:,1), 9), ...
                        PengPara(all_nuc_ion(:,1), 5), PengPara(all_nuc_ion(:,1), 10)]; %每个原子的散射参数
    all_nuc_ion(:,2)=[];
end
return;