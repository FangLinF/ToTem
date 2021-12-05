

function Drawsupercell_figure(x_new, atomsize, x_mysavenewz, top, bottom)  %所有原子的坐标， 原子的放大尺寸, 存的z轴高度，top和bottom的值
%绘制原子结构
%需要包含比较多的信息，但是现在这个pdb信息较少,第一列是元素质子数，二到四是坐标，第五列占据率，第六列假装是价态
atompos = x_new(:,2:4);
atompos(:,4) = x_mysavenewz+top; 
%画图显示一下结构
figure; 
hold on;
view([0 0 -1])

maxx=max(atompos(:,1));
maxy=max(atompos(:,2));
maxz=max(atompos(:,3))+bottom;
mm=max([maxx, maxy, maxz]); 

hold on; line([0 0], [0 0], [0, maxz],'color','b')
line([0 0], [0 maxy], [0, 0],'color','g')
line([0 maxx], [0 0], [0, 0],'color','r')
line([0 0], [0 0], [0, -mm/10],'color','b','linestyle','--')
line([0 0], [0 -mm/10], [0, 0],'color','g','linestyle','--')
line([0 -mm/10], [0 0], [0, 0],'color','r','linestyle','--')

item_ele=unique(x_new(:,1));  %考察有几种原子
for i=1:length(item_ele)  %元素的种类
    r=1-mod(item_ele(i),5)*0.25;    %根据原子序数，设置该原子特有的颜色，r,g,b为三基色的值
    g=1-mod(floor(item_ele(i)/5),5)*0.25;
    b=1-mod(floor(item_ele(i)/25),5)*0.25;
    %居中画图
    plot3(atompos(find(x_new(:,1)==item_ele(i)),1),...
        atompos(find(x_new(:,1)==item_ele(i)),2), ...
        atompos(find(x_new(:,1)==item_ele(i)),3),'o',...  %绘点，形状为‘o’
         'MarkerEdgeColor',[0 0 0],...    %设置边缘颜色为黑色
         'MarkerSize',atomsize.*(20-log2(2)),...  %设置原子的大小,随着显示cell数目的增加大小会变小
         'MarkerFaceColor',[r g b]);    %设置原子的颜色
end
axis equal
xlabel('x'); ylabel('y'); zlabel('z');