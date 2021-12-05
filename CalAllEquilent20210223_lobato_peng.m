function [ele_n, ele_n_i, newabsorp_n, newabsorp_n_i,series_n, series_n_i] = ...    %多增加一个参数，表示出射波还要走多长的真空，到达共同特定的底部
            CalAllEquilent20210223_lotabo_peng(all_n, all_n_i, absorp_n, absorp_n_i, eachthick, slicethick, multiceng, paraflag, DBmode);  
%   直接根据层的厚度区间以及原子的坐标，做划分；
%   划分时，multiceng表示上下影响多少层，一般是2*M+1层
%   all_n的参数，
%   slicethick为0 0.2 0.4 0.6 0.8，如果样品是从0-1.0的厚度的话
%20201102修改成新的方式，根据论文里面的公式，把等价原子位置找到，并且给出a和b的系数，以及特定原子的矩阵。
%系数详见cuda编程以及matlab编程    第六点

%DBmode是计算Debye的方式，只针对lobato的参数。这是因为lobato的参数在不选择DW衰减和absorptive的时候，即，只用frozen时候，B都等于0了
%DBmode如果是0,就是frozen方式，如果是1，就是DB衰减方法



%输出待计算的原子的信息矩阵，注意ele_n是doc里面扣掉红色的那些；series表示每层包含多少个原子，
series_n=[]; series_n_i=[]; %每层包含了多少的原子
ele_n=[]; ele_n_i=[]; 
newabsorp_n=[]; newabsorp_n_i=[];
%测试
% load randpos8
% all_n(:,2:4)=all_n(:,2:4)+randpos8*0.7;

if ~isempty(all_n)
    for i=1:length(slicethick)
        jj = find(all_n(:,4) >= slicethick(i)-multiceng*eachthick & all_n(:,4) < slicethick(i)+multiceng*eachthick+eachthick);  %找到厚度满足要求的原子序号
        series_n(i) = length(jj);  %记录有多少个原子位于这一层内
        ele_n(end+1:end+length(jj),1:16) = all_n(jj,:);
        
         %根据论文10.b式子
        z1(1:length(jj),1) = ele_n(end-series_n(i)+1:end, 4) -(slicethick(i)+eachthick);
        z2(1:length(jj),1) = ele_n(end-series_n(i)+1:end, 4) - slicethick(i);
        if i==1
            z2(:)=inf;
        end
        if i==length(slicethick)
            z1(:)=-inf;  
        end
        z2(find(abs(z2)>=multiceng*eachthick)) = inf;
        z1(find(abs(z1)>=multiceng*eachthick)) = -inf;
        ele_n(end-series_n(i)+1:end, 17) = project_poten_contribute(ele_n(end-series_n(i)+1:end,8),  z1, z2);
        ele_n(end-series_n(i)+1:end, 18) = project_poten_contribute(ele_n(end-series_n(i)+1:end,10),  z1, z2);
        ele_n(end-series_n(i)+1:end, 19) = project_poten_contribute(ele_n(end-series_n(i)+1:end,12),  z1, z2);
        ele_n(end-series_n(i)+1:end, 20) = project_poten_contribute(ele_n(end-series_n(i)+1:end,14),  z1, z2);
        ele_n(end-series_n(i)+1:end, 21) = project_poten_contribute(ele_n(end-series_n(i)+1:end,16),  z1, z2);
        
        
        if ~isempty(absorp_n)  %检查是否要计算吸收
            newabsorp_n(end+1:end+length(jj),1:10) = absorp_n(jj,:);
            newabsorp_n(end-series_n(i)+1:end, 11) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,2),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 12) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,4),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 13) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,6),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 14) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,8),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 15) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,10),  z1, z2);
        end
        z1=[];  %需要清零，否则数据可能会串
        z2=[];
    end
    %乘之后，按照要求输出
    ele_n(:,7)=ele_n(:,7).*ele_n(:,17);
    ele_n(:,9)=ele_n(:,9).*ele_n(:,18);
    ele_n(:,11)=ele_n(:,11).*ele_n(:,19);
    ele_n(:,13)=ele_n(:,13).*ele_n(:,20);
    ele_n(:,15)=ele_n(:,15).*ele_n(:,21);
    ele_n(:,17:21) = [];

    if paraflag == 'l'  
        temp=ele_n(:,5);
    end
    ele_n(:,5)=[]; %删掉DB
    ele_n(:,4)=[]; %删掉Z
    ele_n(:,1)=[]; %删掉原子序数
    if paraflag == 'l' 
        ele_n(:,end+1)=temp*DBmode;
    end
    
    %这里之上，输出结果为：
    % ele_n中，分别是：
    % 原子序数，坐标1-3，DW因子，原子性*占据率，ab序数的序号，每个高斯的比例
    %之后调整为本程序的输出：
    % ele_n必须是13列
    %坐标1-2，占据率*离子3，ab系数   
 
    if ~isempty(absorp_n)  %检查是否要计算吸收
        %乘之后，按照要求输出
       newabsorp_n(:,1)=newabsorp_n(:,1).*newabsorp_n(:,11);
       newabsorp_n(:,3)=newabsorp_n(:,3).*newabsorp_n(:,12);
       newabsorp_n(:,5)=newabsorp_n(:,5).*newabsorp_n(:,13);
       newabsorp_n(:,7)=newabsorp_n(:,7).*newabsorp_n(:,14);
       newabsorp_n(:,9)=newabsorp_n(:,9).*newabsorp_n(:,15);

       newabsorp_n(:,11:15)=[];
    end
end



% xx=[1 2 3; 1 3 2; 1 2 3; 1 0 0];
% [a,b,c]=unique(xx,'rows')
% 
% a =
% 
%      1     0     0
%      1     2     3
%      1     3     2
% 
% 
% b =
% 
%      4
%      1
%      2
% 
% 
% c =
% 
%      2
%      3
%      2
%      1


if ~isempty(all_n_i)
    for i=1:length(slicethick)
        jj = find(all_n_i(:,4) >= slicethick(i)-multiceng*eachthick & all_n_i(:,4) < slicethick(i)+multiceng*eachthick+eachthick);  %找到厚度满足要求的原子序号
        series_n_i(i) = length(jj);  %记录有多少个原子位于这一层内
        ele_n_i(end+1:end+length(jj),1:28) = all_n_i(jj,:);
        
         %根据论文10.b式子
        z1(1:length(jj),1) = ele_n_i(end-series_n_i(i)+1:end, 4) -(slicethick(i)+eachthick);
        z2(1:length(jj),1) = ele_n_i(end-series_n_i(i)+1:end, 4) - slicethick(i);
        if i==1
            z2(:)=inf;
        end
        if i==length(slicethick)
            z1(:)=-inf;  
        end
        z2(find(abs(z2)>=multiceng*eachthick)) = inf;
        z1(find(abs(z1)>=multiceng*eachthick)) = -inf;
        
        ele_n_i(end-series_n_i(i)+1:end, 29) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,9),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 30) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,11),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 31) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,13),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 32) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,15),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 33) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,17),  z1, z2);

        ele_n_i(end-series_n_i(i)+1:end, 34) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,20),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 35) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,22),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 36) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,24),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 37) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,26),  z1, z2);
        ele_n_i(end-series_n_i(i)+1:end, 38) = project_poten_contribute(ele_n_i(end-series_n_i(i)+1:end,28),  z1, z2);

        
        if ~isempty(absorp_n_i)  %检查是否要计算吸收
            newabsorp_n_i(end+1:end+length(jj),1:10) = absorp_n_i(jj,:);
            newabsorp_n_i(end-series_n_i(i)+1:end, 11) = project_poten_contribute(newabsorp_n_i(end-series_n_i(i)+1:end,2),  z1, z2);
            newabsorp_n_i(end-series_n_i(i)+1:end, 12) = project_poten_contribute(newabsorp_n_i(end-series_n_i(i)+1:end,4),  z1, z2);
            newabsorp_n_i(end-series_n_i(i)+1:end, 13) = project_poten_contribute(newabsorp_n_i(end-series_n_i(i)+1:end,6),  z1, z2);
            newabsorp_n_i(end-series_n_i(i)+1:end, 14) = project_poten_contribute(newabsorp_n_i(end-series_n_i(i)+1:end,8),  z1, z2);
            newabsorp_n_i(end-series_n_i(i)+1:end, 15) = project_poten_contribute(newabsorp_n_i(end-series_n_i(i)+1:end,10),  z1, z2);
        end
        z1=[];  %需要清零，否则数据可能会串
        z2=[];
        
    end
    %乘之后，按照要求输出
       ele_n_i(:,8)=ele_n_i(:,8).*ele_n_i(:,29);
       ele_n_i(:,10)=ele_n_i(:,10).*ele_n_i(:,30);
       ele_n_i(:,12)=ele_n_i(:,12).*ele_n_i(:,31);
       ele_n_i(:,14)=ele_n_i(:,14).*ele_n_i(:,32);
       ele_n_i(:,16)=ele_n_i(:,16).*ele_n_i(:,33);
       
       ele_n_i(:,19)=ele_n_i(:,19).*ele_n_i(:,34);
       ele_n_i(:,21)=ele_n_i(:,21).*ele_n_i(:,35);
       ele_n_i(:,23)=ele_n_i(:,23).*ele_n_i(:,36);
       ele_n_i(:,25)=ele_n_i(:,25).*ele_n_i(:,37);
       ele_n_i(:,27)=ele_n_i(:,27).*ele_n_i(:,38);

       if paraflag == 'l'
          temp=ele_n_i(:,5);
       end
       ele_n_i(:,29:end)=[];
       ele_n_i(:,5)=[]; %删掉DB
       ele_n_i(:,4)=[]; %删掉Z
       ele_n_i(:,1)=[]; %删掉原子序数
       
       if paraflag == 'l'
          ele_n_i(:,end+1)=temp*DBmode;
       end
       
          

    if ~isempty(newabsorp_n_i)  %检查是否要计算吸收       
        %乘之后，按照要求输出
       newabsorp_n_i(:,1)=newabsorp_n_i(:,1).*newabsorp_n_i(:,11);
       newabsorp_n_i(:,3)=newabsorp_n_i(:,3).*newabsorp_n_i(:,12);
       newabsorp_n_i(:,5)=newabsorp_n_i(:,5).*newabsorp_n_i(:,13);
       newabsorp_n_i(:,7)=newabsorp_n_i(:,7).*newabsorp_n_i(:,14);
       newabsorp_n_i(:,9)=newabsorp_n_i(:,9).*newabsorp_n_i(:,15);

       newabsorp_n_i(:,11:15)=[];
    end
    
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


