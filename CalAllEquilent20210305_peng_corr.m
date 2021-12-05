function [corr_ele_n, corr_ele_n_i, corr_info, series_n_corr, series_n_i_corr] = ...    %多增加一个参数，表示出射波还要走多长的真空，到达共同特定的底部
            CalAllEquilent20210305_peng_corr(corr_nuclear_copy, corr_peng_nuc_ion_copy, eachthick, slicethick, multiceng);   
discon=0.015;  %如果与上下层的距离在discon以内，就分2层，如果是大于这个值，就是1层内
        %   这里增加势场的修正的原子位置分析。
%   直接根据层的厚度区间以及原子的坐标，做划分；
%   划分时，顶多分2层
%   corr_nuclear_copy, corr_peng_nuc_ion_copy的输入的种类，原子坐标，参数，
%   slicethick为0 0.2 0.4 0.6 0.8，如果样品是从0-1.0的厚度的话
%20201102修改成新的方式，根据论文里面的公式，把等价原子位置找到，并且给出a和b的系数，以及特定原子的矩阵。
%系数详见cuda编程以及matlab编程    第六点

%这里使用的是和lobato一样的处理方法，所以忽略这个原子了
% % % %DBmode是计算Debye的方式，只针对lobato的参数。这是因为lobato的参数在不选择DW衰减和absorptive的时候，即，只用frozen时候，B都等于0了
% % % %DBmode如果是0,就是frozen方式，如果是1，就是DB衰减方法



%输出待计算的原子的信息矩阵，注意ele_n是doc里面扣掉红色的那些；series表示每层包含多少个原子，
corr_ele_n=[]; corr_ele_n_i=[]; 
series_n_corr=[]; series_n_i_corr=[];%每层包含了多少的原子
%测试
% load randpos8
% corr_nuclear_copy(:,2:4)=corr_nuclear_copy(:,2:4)+randpos8*0.7;

if ~isempty(corr_nuclear_copy)
    [a, b, c] = unique(corr_nuclear_copy(:,7:26),'rows');  %返回的a和b，a是不重复的系数；c是归类后的序号
    
    corr_info = a ;  %把元素的信息存这里
    corr_nuclear_copy(:,7)=c;
    corr_nuclear_copy(:,8:end)=[];   %把第八列之后的都删掉
    
    for i=1:length(slicethick)
        jj = find(corr_nuclear_copy(:,4) >= slicethick(i) & corr_nuclear_copy(:,4) < slicethick(i)+eachthick);  %找到厚度满足要求的原子序号
        templen = length(jj);  %记录有多少个原子位于这一层内
        corr_ele_n(end+1:end+length(jj),1:11) = [corr_nuclear_copy(jj,:).';... %7个元素
            1.0.*ones(1,templen);...%计算的权重
            i.*ones(1,templen); ...%位于哪个原子层数
            corr_nuclear_copy(jj,4).'.*ones(1,templen) - slicethick(i);...%距离上表面的距离
            corr_nuclear_copy(jj,4).'.*ones(1,templen) - slicethick(i)-eachthick].';  %距离下表面的距离
    end
    if eachthick > 2*discon & multiceng>=1 %如果要求分层；且，如果层厚大于2倍的discon，就接着将原子能量进行分层；否则有可能会出现2次能量增大问题
        %如果原子距离上表面小于discon，且，不是位于第1层，就分一半的能量
        jj =  find( corr_ele_n(:,end-1) < discon & corr_ele_n(:,end-2) >1 );
        %如果原子距离下表面大于-discon，且，不是位于最后层，就分一半的能量
        jjjj =  find( corr_ele_n(:,end) > -discon & corr_ele_n(:,end-2) <length(slicethick) );
        
        if ~isempty(jj)
            corr_ele_n(jj,8) = 0.5; %分一半能量             
            corr_ele_n(end+1:end+length(jj),:) = corr_ele_n(jj,:); %做拷贝
            corr_ele_n(end:-1:end-length(jj)+1,end-2) = corr_ele_n(end:-1:end-length(jj)+1,end-2)-1;  %层号减1
        end
        
        if ~isempty(jjjj)
            corr_ele_n(jjjj,8) = 0.5; %分一半能量
            corr_ele_n(end+1:end+length(jjjj),:) = corr_ele_n(jjjj,:); %做拷贝
            corr_ele_n(end:-1:end-length(jjjj)+1,end-2) = corr_ele_n(end:-1:end-length(jjjj)+1,end-2)+1;  %层号减1
        end
    end
    corr_ele_n(:,end-1:end) = [];  %去掉距离上下表面的信息 %就剩下9个元素
    [a,c] = sort(corr_ele_n(:,end)); %根据最后的层信息，重新排序；
    corr_ele_n = corr_ele_n(c,:);
    
    for i=1:length(slicethick)
        series_n_corr(i) = length( find(corr_ele_n(:,end)==i));  %找到每层包含多少个原子
    end
    
        %这里去掉多余的信息
    corr_ele_n(:,6) = corr_ele_n(:,6).*corr_ele_n(:,8);
    corr_ele_n(:,[9,8,4,1])=[];
end


%这部分没法测试是否正确，关于离子性和原子性都有
if ~isempty(corr_peng_nuc_ion_copy)
    [a, b, c] = unique(corr_peng_nuc_ion_copy(:,8:27),'rows');  %返回的a和b，a是不重复的系数；c是归类后的序号   
    corr_info = a ;  %把元素的信息存这里  
    kindnum = length(corr_info(:,1)); %总共有多少类
    corr_peng_nuc_ion_copy(:,8)=c;
    
    [a, b, c] = unique(corr_peng_nuc_ion_copy(:,29:48),'rows');  %返回的a和b，a是不重复的系数；c是归类后的序号   
    %corr_info_n_i = a ;  %把元素的信息存这里 
    corr_info(end+1:end+length(a(:,1)),:) = a+kindnum;  %每个序号数+原本离子性的序号
    corr_peng_nuc_ion_copy(:,29)=c+kindnum;  %找到对应的离子序号
    
    corr_peng_nuc_ion_copy(:,30:end)=[];   %把第16列之后的都删掉
    corr_peng_nuc_ion_copy(:,9:27)=[];   %把第9列之后的都删掉  只留下原来第28列是离子性比例；和29列赋值的离子序号
    
    for i=1:length(slicethick)
        jj = find(corr_peng_nuc_ion_copy(:,4) >= slicethick(i) & corr_peng_nuc_ion_copy(:,4) < slicethick(i)+eachthick);  %找到厚度满足要求的原子序号
        templen = length(jj);  %记录有多少个原子位于这一层内
        corr_ele_n_i(end+1:end+length(jj),1:14) = [corr_peng_nuc_ion_copy(jj,:).';... %10个元素
            1.0.*ones(1,templen);...%计算的权重
            i.*ones(1,templen); ...%位于哪个原子层数
            corr_peng_nuc_ion_copy(jj,4).'.*ones(1,templen) - slicethick(i);...%距离上表面的距离
            corr_peng_nuc_ion_copy(jj,4).'.*ones(1,templen) - slicethick(i)-eachthick].';  %距离下表面的距离
    end
    if eachthick > 2*discon  & multiceng>=1 %如果要求分层，且，如果层厚大于2倍的discon，就接着将原子能量进行分层；否则有可能会出现2次能量增大问题
        %如果原子距离上表面小于discon，且，不是位于第1层，就分一半的能量
        jj =  find( corr_peng_nuc_ion_copy(:,end-1) < discon & corr_peng_nuc_ion_copy(:,end-2) >1 );
        %如果原子距离下表面大于-discon，且，不是位于最后层，就分一半的能量
        jjjj =  find( corr_peng_nuc_ion_copy(:,end) > -discon & corr_peng_nuc_ion_copy(:,end-2) <length(slicethick) );
        
        if ~isempty(jj)
            corr_ele_n_i(jj,11) = 0.5; %分一半能量             
            corr_ele_n_i(end+1:end+length(jj),:) = corr_ele_n_i(jj,:); %做拷贝
            corr_ele_n_i(end:-1:end-length(jj)+1,end-2) = corr_ele_n_i(end:-1:end-length(jj)+1,end-2)-1;  %层号减1
        end
        
        if ~isempty(jjjj)
            corr_ele_n_i(jjjj,11) = 0.5; %分一半能量
            corr_ele_n_i(end+1:end+length(jjjj),:) = corr_ele_n_i(jjjj,:); %做拷贝
            corr_ele_n_i(end:-1:end-length(jjjj)+1,end-2) = corr_ele_n_i(end:-1:end-length(jjjj)+1,end-2)+1;  %层号减1
        end
        
    end
    %corr_ele_n_i(:,end-1:end) = [];  %去掉距离上下表面的信息 %就剩下12个元素
    [a,c] = sort(corr_ele_n_i(:,end)); %根据最后的层信息，重新排序；
    corr_ele_n_i = corr_ele_n_i(c,:);
    
      
    for i=1:length(slicethick)
        series_n_i_corr(i) = length( find(corr_ele_n_i(:,end)==i)); 
    end
    
    %这里去掉多余的信息
    corr_ele_n_i(:,6) = corr_ele_n_i(:,6).*corr_ele_n_i(:,11);  %占据率
    corr_ele_n_i(:,[12,11,4,1])=[];  %上下层的距离，实际高度，质子数
end

   

pp=1;




        
function mypercent = project_poten_contribute(B, thick1, thick2);  %高斯函数的A和B，以及高斯函数积分的范围
% if thick2==inf   %从某个高度到无穷大的积分。积分针对的函数。详见DOC
%     mypercent=0.5* erfc(pi*thick1./sqrt(B));
% else
% %   ---------- thick1
% % 
% %       o atom sympble
% %  ----------- thick1+eachthick   给出原子到上顶的距离；到下底的距离就是zheight_2
%     mypercent =0.5* erf(pi*thick1./sqrt(B)) + 0.5* erf(pi*thick2./sqrt(B));  %将本层内的强度进行叠加
% end
%  如果离上表面和下表面足够远，这个数为1；另外，0.5乘上去，是因为erf函数中自变量如果很大，是等于1；也就是，默认偶函数两倍放大
% 而本函数，对上下表面的积分的数值并不一样。

    mypercent=abs(0.5* erf(2*pi*thick2(:)./sqrt(B(:))) - 0.5* erf(2*pi*thick1(:)./sqrt(B(:)))) ;


