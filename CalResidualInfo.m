function [resiinfo, series_resi] = ...
             CalResidualInfo(all_n, all_n_i, eachthick, slicethick, multiceng);  
         
         series_n=[]; series_n_i=[]; %每层包含了多少的原子
ele_n=[]; ele_n_i=[]; 
newabsorp_n=[]; newabsorp_n_i=[];


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

    ele_n(:,5)=[]; %删掉DB
    ele_n(:,4)=[]; %删掉Z
    ele_n(:,1)=[]; %删掉原子序数
    
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