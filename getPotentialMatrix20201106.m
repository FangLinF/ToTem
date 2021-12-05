function [ele_n_potentialMatix, ele_n_i_potentialMatix, absorp_n_potentialMatix, absorp_n_i_potentialMatix] = getPotentialMatrix20201106( s2_green, PARAMETER, ele_n_para4matrix, ele_n_i_para4matrix, newabsorp_n_para4matrix, newabsorp_n_i_para4matrix );
%计算势场矩阵，返回的是4维矩阵，N个原子*5*sx*sy
ele_n_potentialMatix=[]; 
ele_n_i_potentialMatix=[]; 
absorp_n_potentialMatix=[];
absorp_n_i_potentialMatix=[];

if ~isempty(ele_n_para4matrix) 
    ele_n_potentialMatix = getPotentialComm( s2_green, PARAMETER, ele_n_para4matrix );
end
if ~isempty(ele_n_i_para4matrix) 
    ele_n_i_potentialMatix = getPotentialComm( s2_green, PARAMETER, ele_n_i_para4matrix );
end
if ~isempty(newabsorp_n_para4matrix) 
    absorp_n_potentialMatix = getPotentialComm( s2_green, PARAMETER, newabsorp_n_para4matrix );
end
if ~isempty(newabsorp_n_i_para4matrix) 
    absorp_n_i_potentialMatix = getPotentialComm( s2_green, PARAMETER, newabsorp_n_i_para4matrix );
end

return;

%计算五个高斯矩阵，并且加入了parameter的参数
function potentialMatrix = getPotentialComm(s2, PARAMETER, ab_coff);
[sx, sy]=size(s2);
potentialMatrix = zeros(sx, sy, 5, length(ab_coff(:,1)));
for i=1:length(ab_coff(:,1))
     potentialMatrix(1:sx, 1:sy, 1, i) = PARAMETER.* ab_coff(i,1).*exp(-s2.*ab_coff(i,2)); 
     potentialMatrix(1:sx, 1:sy, 2, i) = PARAMETER.* ab_coff(i,3).*exp(-s2.*ab_coff(i,4)); 
     potentialMatrix(1:sx, 1:sy, 3, i) = PARAMETER.* ab_coff(i,5).*exp(-s2.*ab_coff(i,6)); 
     potentialMatrix(1:sx, 1:sy, 4, i) = PARAMETER.* ab_coff(i,7).*exp(-s2.*ab_coff(i,8)); 
     potentialMatrix(1:sx, 1:sy, 5, i) = PARAMETER.* ab_coff(i,9).*exp(-s2.*ab_coff(i,10)); 
end
 
