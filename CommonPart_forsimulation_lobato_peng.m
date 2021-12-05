function [all_nuclear, all_nuc_ion, absorp_n, absorp_n_i, corr_peng_nuc, corr_peng_nuc_ion]=CommonPart_forsimulation_lotabo_peng(hObject, eventdata, handles, flag);
%首先确定是否考虑吸收，如果考虑吸收，增加一个吸收矩阵。
%基础矩阵，考虑原子性和（或）离子性；对应一个（或空矩阵）吸收矩阵
%是否把B做衰减，如果是，所有的吸收和原子性和离子性，都在b上修改参数。
%step 0
%根据所选的模拟区域，计算修正原子的坐标位置;需要结合show的运行结果
load allresul.mat   %现在是从存储的allresul读取x_new结果
x=allresul.x_new;
% x=handles.x_new;

%x(:,2:4)=handles.atompos;   
x=x(find(x(:,2)>=handles.tl_green(1) & x(:,2)<=handles.rd_green(1) ...
    & x(:,3)>=handles.tl_green(2) & x(:,3)<=handles.rd_green(2)),:);%把位于绿色框外的原子都删掉
x(:,2:3)=x(:,2:3)-repmat(handles.tl_green, length(x(:,1)),1);   %需要转化到绿色框里面的有效位置。
%tl_green和rd_green变量都是一样的，对于STRM,HRTEM和CBED


%step 1
%计算弹性势的矩阵，分弹性势1和弹性势2。
              %弹性势1只考虑原子性或离子性
              %弹性势2同时考虑原子势和离子性
[all_nuclear,all_nuc_ion]  = getpotentialpara_lotabo_peng(x, flag); 
%修改这里
%如果是correct，这里加上修正的函数-------------------
corr_peng_nuc=[]; corr_peng_nuc_ion=[];
if flag == 'n'
    [corr_peng_nuc, corr_peng_nuc_ion] = getpotentialpara_corr_peng(x);  %获得关于peng老师的修正
    if ~get(handles.radiobutton14, 'value') & ~get(handles.radiobutton4, 'value')  %如果是不使用DW衰减和AP方法。
        if ~isempty(corr_peng_nuc)
            corr_peng_nuc(:,5) = 0;  %第五列debye waller因子不影响衰减的势场,于是把dw因子视为0，这样衰减加到势场上不影响结果
        end
        if ~isempty(corr_peng_nuc_ion)
            corr_peng_nuc_ion(:,5) = 0;  %第五列debye waller因子不影响衰减的势场
        end
    end
end
%--------------------------------------------

%计算原子数。
len1=0; len2=0;
if ~isempty(all_nuclear)
    len1=length(all_nuclear(:,1));
end
if ~isempty(all_nuc_ion)
    len2=length(all_nuc_ion(:,1));
end
    
%step 2
%如果考虑吸收absorptive potential
absorp_n=[];  absorp_n_i=[];  % 设置初值
element_n=[]; element_n_i=[];
if get(handles.radiobutton4, 'value')   %计算吸收，与输入电压以及DebyeWaller因子有关
    if ~isempty(all_nuclear)
        absorp_n=zeros(len1, 10);
        %给出所有原子的质子数
        element_n=unique([all_nuclear(:,1).';all_nuclear(:,5).'].','row');      
    end
    if ~isempty(all_nuc_ion)
        absorp_n_i=zeros(len2, 10);
        %给出所有原子的质子数
        element_n_i=unique([all_nuc_ion(:,1).';all_nuc_ion(:,5).'].','row');
    end
    %给出离子性和原子性总共包含的所有元素有哪些
    allelement=unique([element_n;element_n_i],'row');
    for i=1:length(allelement(:,1))  %对每种原子，求他们的吸收势场  20210223这一句是不一样的
          myresul(i,:)=mynewftds_ab_lotabo_peng(allelement(i,1),allelement(i,2),str2num(get(handles.edit_Vol, 'string')),flag);  
          %4个变量，第一个是原子质子数；第二个是B因子，第三个是电压单位是千伏，第四个是计算方法，p和l分别代表peng和lotabo
    end
    
    if ~isempty(element_n)
        absorp_n=zeros(length(all_nuclear(:,1)),10);
        for i=1:length(allelement(:,1));  %每种原子的吸收
            e=find(all_nuclear(:,1)==allelement(i, 1) & all_nuclear(:,5)==allelement(i, 2));   %找到质子数和B一致的原子的序号
            absorp_n(e,:)=repmat(myresul(i,:),length(e),1);  %把算出来的ab系数，都拷贝到对应的原子序号上
        end
    end
    if ~isempty(element_n_i)
        absorp_n_i=zeros(length(all_nuc_ion(:,1)),10);
        for i=1:length(allelement(:,1));  %每种原子的吸收
            e=find(all_nuc_ion(:,1)==allelement(i, 1) & all_nuc_ion(:,5)==allelement(i, 2));   %找到质子数和B一致的原子的序号
            absorp_n_i(e,:)=repmat(myresul(i,:),length(e),1);  %把算出来的ab系数，都拷贝到对应的原子序号上
        end
    end
    
end

%step 3
if flag == 'p'|| flag=='n'  %20210223  %只有彭的系数才做如下操作，或者peng的修正才做b+B
   %如果考虑debye waller的衰减包罗，所有的b相关的系数都需要加上B的常数
   %20201230修改
    %做相应操作。b_re和b_im都加B
    if get(handles.radiobutton14, 'value') || get(handles.radiobutton4, 'value')  %如果是DW选择，就需要加入DW的影响。
        if ~isempty(all_nuclear)
            all_nuclear(:,8) = all_nuclear(:,8)+all_nuclear(:,5);
            all_nuclear(:,10) = all_nuclear(:,10)+all_nuclear(:,5);
            all_nuclear(:,12) = all_nuclear(:,12)+all_nuclear(:,5);
            all_nuclear(:,14) = all_nuclear(:,14)+all_nuclear(:,5);
            all_nuclear(:,16) = all_nuclear(:,16)+all_nuclear(:,5);
        end

        if ~isempty(all_nuc_ion)
            all_nuc_ion(:,9) = all_nuc_ion(:,9)+all_nuc_ion(:,5);
            all_nuc_ion(:,11) = all_nuc_ion(:,11)+all_nuc_ion(:,5);
            all_nuc_ion(:,13) = all_nuc_ion(:,13)+all_nuc_ion(:,5);
            all_nuc_ion(:,15) = all_nuc_ion(:,15)+all_nuc_ion(:,5);
            all_nuc_ion(:,17) = all_nuc_ion(:,17)+all_nuc_ion(:,5);
        
            all_nuc_ion(:,20) = all_nuc_ion(:,20)+all_nuc_ion(:,5);
            all_nuc_ion(:,22) = all_nuc_ion(:,22)+all_nuc_ion(:,5);
            all_nuc_ion(:,24) = all_nuc_ion(:,24)+all_nuc_ion(:,5);
            all_nuc_ion(:,26) = all_nuc_ion(:,26)+all_nuc_ion(:,5);
            all_nuc_ion(:,28) = all_nuc_ion(:,28)+all_nuc_ion(:,5);
        
        end
    end
end
if flag == 'l'  %20210223  %如果时lotabo的系数不做任何操作
    if get(handles.radiobutton14, 'value') || get(handles.radiobutton4, 'value')  %如果是DW选择，就需要加入DW的影响。
         disp('Lotabo parameters cannot add debye waller in this step, the potential will be add exp(-Bs^2)， not here')
    end
end
    
if get(handles.radiobutton4, 'value')   %如果考虑AP近似的曲线  
    %这个方法是无论peng老师还是lotabo系数拟合出来的高斯线，都需要做的b的系数的加0.5B   
    if ~isempty(absorp_n)
        absorp_n(:,2)=absorp_n(:,2) + 0.5*all_nuclear(:,5);   %注意，这里有0.5的系数
        absorp_n(:,4)=absorp_n(:,4) + 0.5*all_nuclear(:,5);
        absorp_n(:,6)=absorp_n(:,6) + 0.5*all_nuclear(:,5);
        absorp_n(:,8)=absorp_n(:,8) + 0.5*all_nuclear(:,5);
        absorp_n(:,10)=absorp_n(:,10) + 0.5*all_nuclear(:,5);
    end
    
    if ~isempty(absorp_n_i)
        absorp_n_i(:,2)=absorp_n_i(:,2) + 0.5*all_nuc_ion(:,5);
        absorp_n_i(:,4)=absorp_n_i(:,4) + 0.5*all_nuc_ion(:,5);
        absorp_n_i(:,6)=absorp_n_i(:,6) + 0.5*all_nuc_ion(:,5);
        absorp_n_i(:,8)=absorp_n_i(:,8) + 0.5*all_nuc_ion(:,5);
        absorp_n_i(:,10)=absorp_n_i(:,10) + 0.5*all_nuc_ion(:,5);
    end
end

return