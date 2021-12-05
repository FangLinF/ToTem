function resul_pos = randshiftpos(ini_pos, eachthick, slicethick, mid_slice_num);
rng('shuffle');
resul_pos =  ini_pos;
randpos=randn(length(ini_pos(:,1)),3).*repmat(sqrt(ini_pos(:,5)/8/pi/pi),1,3);   %随机坐标的变化量
resul_pos(:,2:4)=randpos+ini_pos(:,2:4);  %变成新的随机坐标，之后使用allnuclear_copy进行计算  %考虑是3个方向的振动，所以三个方向分量需要扣掉sqrt(3)的部分 20210103
%resul_pos(:,2:3)=randpos(:,1:2)+ini_pos(:,2:3);
while length( find(resul_pos(:,4)<slicethick(1) | resul_pos(:,4)>=slicethick(end)+eachthick) ) >0  
    %if any atom shift to the boundaries of crystal, its position will be
    %random again
    i = find((resul_pos(:,4)<slicethick(1) | resul_pos(:,4)>=slicethick(end)+eachthick));
    randpos=randn(length(i),3).*repmat(sqrt(ini_pos(i,5)/8/pi/pi),1,3);
    resul_pos(i,2:4)=randpos+ini_pos(i,2:4);
    %resul_pos(i,2:3)=randpos(:,1:2)+ini_pos(i,2:3);
end

if ~isempty(mid_slice_num)
    %atoms must be within their special slices
    for j = 1 : length(mid_slice_num)
        %确保位于哪层内的原子，振动后依然位于哪层
        
        %找到哪些原子在振动前是属于midslice层内的
        k = find( ini_pos(:,4) < slicethick(mid_slice_num(j)) + eachthick ); % find atoms that within the mid_slice_num thickness
        %for example, 10 12 14 15 atoms should be within mid_slice_num
        
        %找到哪些原子在振动后是属于midslice层内的
        kk = find( resul_pos(:,4) < slicethick(mid_slice_num(j)) + eachthick );
        %for example after shift, 10 12 15 17 atoms are within mid_slice_num
        
        while ~isequaltwo( k, kk)  %振动前后应该一致
            [kandkk, ik, ikk]= intersect(k,kk);   %k和kk的交集
            % for example, 10, 12, 15 will be the intersection
            
            somekk = setdiff(kk,kandkk);   %找到那些振动出层的原子的序号,找到补集
            % for example, 17 means that the 17th atoms should be shifted
            somek  = setdiff(k, kandkk);   %找到那些振动出层的原子的序号,找到补集
            % for example, 14 means that the 14th atoms should be shifted
            
            randpos=randn(length([somek; somekk]),3).*repmat(sqrt(ini_pos([somek; somekk],5)/8/pi/pi),1,3);
            resul_pos([somek; somekk],2:4) = randpos+ini_pos([somek; somekk],2:4);
            %resul_pos([somek; somekk],2:3) = randpos(:,1:2)+ini_pos([somek; somekk],2:3);
            kk = find( resul_pos(:,4) < slicethick(mid_slice_num(j)) + eachthick );
        end
    end
end
%output the number of atoms before and after vibration
%考察这些层内，分别有多少个原子，振动后的
num = unique([mid_slice_num,length(slicethick)]);
for i = 1: length(num)
    %哪些原子位于第kk层；
    disp(strcat('After vibration, in selected region, within the ', num2str(num(i)), 'th slice, there are '));
    
    kk = find( resul_pos(:,4) >= slicethick(1) & resul_pos(:,4) < slicethick(num(i))+eachthick );

    allelement = unique( resul_pos(kk,1) );  %找到有多少种原子
    for j= 1 : length(allelement)
        elementnum = length ( find( resul_pos(kk,1) == allelement(j)));
        disp(strcat(' Element ', num2str(allelement(j)), ':', num2str(elementnum)));
    end
    
    
    
    disp(strcat('Before vibration, within the ', num2str(num(i)), 'th slice, there are '));
    
    kk = find( ini_pos(:,4) >= slicethick(1) & ini_pos(:,4) < slicethick(num(i))+eachthick );

    allelement = unique( resul_pos(kk,1) );  %找到有多少种原子
    for j= 1 : length(allelement)
        elementnum = length ( find( ini_pos(kk,1) == allelement(j)));
        disp(strcat(' Element ', num2str(allelement(j)), ':', num2str(elementnum)));
    end
end
return;

function flag = isequaltwo( k, kk)  %判断两个矩阵是否完全一样 find whether two vectors are the same; 
%same flag=1; different flag=0
if length(k)~=length(kk)  %两个矩阵不一样
    flag = 0;
    return;
else
    if sum( k == kk) == length(k)  %两个矩阵一样
        flag = 1;
    else
        flag = 0;  %两个矩阵不一样
    end
end
    
    