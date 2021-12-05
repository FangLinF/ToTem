function varargout = ToTEM_submit_v1(varargin)
% TOTEM_SUBMIT_V1 MATLAB code for ToTEM_submit_v1.fig

% Last Modified by GUIDE v2.5 27-May-2021 11:47:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ToTEM_submit_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @ToTEM_submit_v1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ToTEM_submit_v1 is made visible.
function ToTEM_submit_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ToTEM_submit_v1 (see VARARGIN)

% Choose default command line output for ToTEM_submit_v1
handles.output = hObject;

handles.exepath=cd; %记住exe文件所在的文件夹
handles.execd=cd;  %记录exe文件所在的目录
handles.tempcd=cd;  %用来记录上一个打开的目录位置，方便文件读写
handles.saveresult = cd; % remember the current dir to save the result

set(handles.text21, 'visible', 'off');   % aper的mrad和1/nm单位
set(handles.text125, 'visible', 'on');
set(handles.text24, 'visible', 'on');   % convergence & aperture的mrad和1/nm单位
set(handles.text134, 'visible', 'off');
    
axis(handles.axes1); axis off;axis(handles.axes2); axis off;


% %x和y方向，需要左右扩展4埃
handles.outside_ext=4;   
%如果是计算每类原子的势场，不是每个原子的计算
handles.TOTALELEment=10;
%考虑振动的次数
handles.vibration=30;
%显示原子时，原子的尺寸
handles.dis_atomsize = 1;

%一次能够计算多少个probe的传播。
handles.mynode=70;

%STEM激活
set(handles.radiobutton11,'value',1)
% Update handles structure

handles.conv_source=0; 
handles.conv_sampling=1; %把卷积的点源尺寸，以及放大率抽样；0和1的初值不进行任何操作

set(handles.GPURB,'value',1);  %默认GPU计算


%把批量生图的一些信息存下来，设定初值
%文件路径
handles.batchs.PathName = cd;
handles.batchs.wavesx=320;  handles.batchs.wavesy=handles.batchs.wavesx;
handles.batchs.showorrun = 'R';
handles.batchs.totalnumber =500;
handles.batchs.processnum = 2;
handles.batchs.GPUBatch = 5;
handles.batchs.top_left = [101,51];
handles.batchs.width_heigh = [50,70];

guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = ToTEM_submit_v1_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
function Untitled_4_Callback(hObject, eventdata, handles)
function Untitled_10_Callback(hObject, eventdata, handles)
function Untitled_1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %进入上一次打开的目录
[FileName,PathName]=uigetfile({'*.pdb','Program Database File (*.pdb)'});%打开文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

cd(handles.execd);  %回到旧的文件夹
handles.tempcd=PathName;  %记录上一个打开的文件夹

%读取pdb或者cif的具体内容
if sum(FileName(end-2:end)=='pdb')==3  %表示是pdb格式
    xxx=readpdb(strcat(PathName, FileName));
else
    disp('Can not read this file');
    return;  
end
handles.x=[];
handles.x(:,1:4)=xxx(:,1:4);  %element，coordinates
handles.x(:,5) = 0;  %价态
handles.x(:,6) = 0;  %离子性
handles.x(:,7) = xxx(:,6); %dw %will input from interface
handles.x(:,8) = xxx(:,5); %occupy

%input DW factors
    z=unique(handles.x(:,1));
    for i=1:length(z)
        prompt={sprintf('DW factor for the %.0f Atom ',z(i))};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'0'};
        DW=inputdlg(prompt,dlg_title,num_lines,defaultans);
        handles.x(find(handles.x(:,1)==z(i)),7)=str2num(cell2mat(DW));
    end
    temp = handles.x;
%     xlswrite(strcat(PathName, FileName(1:end-4)),temp);
%     disp(strcat('Parameters of atoms are saved in the file', strcat(PathName, FileName(1:end-4)),'.xlsx'));

%handles.x is a variable including all information of atoms, the atoms'
%coordinators is in the range of [0, MaxValue]. Specially, handles.x will
%not be changed in following codes.------begin 202012290930
atompos=handles.x(:,2:4);  % 单位变埃~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % 原子高度，最小值变为0，且增加一层真空层厚度，方便分层
                                                            % 并且保证不会位于第0层。因为这里最少是位于第二层，下两行取了floor的行数。
                                                            % 比如厚度是3,0-3-6-9；加了1.5*eachslice后，值为4.5。分在第一层正中间

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；
                %只有handles.x中的坐标需要改变一下；只有在设置层厚时候；
handles.x(:,2:4)=atompos;  %重新给原子的坐标，保证都是大于0的，且有边界
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %把旋转前的结构记录下来；handles.x_new记录的是旋转后
Drawsupercell(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290935
guidata(hObject, handles);


function Drawsupercell(hObject, eventdata, handles)
%绘制原子结构 draw atoms of crystal in axis1 and axis2 from handles.x_new
%需要包含比较多的信息，但是现在这个pdb信息较少,第一列是元素质子数，二到四是坐标，第五列占据率，第六列假装是价态
atompos=handles.x_new(:,2:4);
axis(handles.axes1); hold off
cla(handles.axes1);   %清除原有图像
cla(handles.axes2);   %清除原有图像

%画图显示一下结构
axes(handles.axes1);   
hold on;
view([0 0 -1])

maxx=max(atompos(:,1));
maxy=max(atompos(:,2));
maxz=max(atompos(:,3));
mm=max([maxx, maxy, maxz]);

hold on; line([0 0], [0 0], [0, mm],'color','b')
line([0 0], [0 mm], [0, 0],'color','g')
line([0 mm], [0 0], [0, 0],'color','r')
line([0 0], [0 0], [0, -mm/10],'color','b','linestyle','--')
line([0 0], [0 -mm/10], [0, 0],'color','g','linestyle','--')
line([0 -mm/10], [0 0], [0, 0],'color','r','linestyle','--')

item_ele=unique(handles.x(:,1));  %考察有几种原子
for i=1:length(item_ele)  %元素的种类 there are many elements
    r=1-mod(item_ele(i),5)*0.25;    %根据原子序数，设置该原子特有的颜色，r,g,b为三基色的值
    g=1-mod(floor(item_ele(i)/5),5)*0.25;
    b=1-mod(floor(item_ele(i)/25),5)*0.25;
    %居中画图
    axes(handles.axes1);  %画原子 draw atoms
    plot3(atompos(find(handles.x(:,1)==item_ele(i)),1),...
        atompos(handles.x(:,1)==item_ele(i),2), ...
        atompos(find(handles.x(:,1)==item_ele(i)),3),'o',...  %绘点，形状为‘o’
         'MarkerEdgeColor',[0 0 0],...    %设置边缘颜色为黑色
         'MarkerSize',handles.dis_atomsize.*(20-log2(2)),...  %设置原子的大小,随着显示cell数目的增加大小会变小
         'MarkerFaceColor',[r g b]);    %设置原子的颜色
     axes(handles.axes2);  hold on;%画图显示一下图例 draw atoms' label
     plot(0.1*i,0.7,'o',...     %显示原子
     'MarkerEdgeColor',[0 0 0],...
     'MarkerSize',15,...
     'MarkerFaceColor',[r g b]);
     text(0.1*i-0.01,0.3,num2str(item_ele(i)));  %在原子下方显示对应的原子质子数
end
axes(handles.axes1);  %画原子
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
axes(handles.axes2);
box off;axis off;   %不显示边框和坐标轴
xlim([0 1]);ylim([0 1]);  %设置x轴与y轴的取值范围
guidata(hObject, handles);
   


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %进入上一次打开的目录
[FileName,PathName]=uigetfile({'*.xls','Excel File (*.xls)';'*.xlsx','Excel File (*.xlsx)'});%打开文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

cd(handles.execd);  %回到旧的文件夹
handles.tempcd=PathName;  %记录上一个打开的文件夹

%读取pdb或者cif的具体内容
if sum(FileName(end-2:end)=='xls')==3 || sum(FileName(end-3:end)=='xlsx')==4 %表示是pdb格式
    handles.x=xlsread(strcat(PathName, FileName));
else
    disp('Cannot read this file!');
    return;
end

atompos=handles.x(:,2:4);  % 单位变埃~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % 原子高度，最小值变为0，且增加一层真空层厚度，方便分层
                                                            % 并且保证不会位于第0层。因为这里最少是位于第二层，下两行取了floor的行数。
                                                            % 比如厚度是3,0-3-6-9；加了1.5*eachslice后，值为4.5。分在第一层正中间

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；
                %只有handles.x中的坐标需要改变一下；只有在设置层厚时候；
handles.x(:,2:4)=atompos;  %重新给原子的坐标，保证都是大于0的，且有边界
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %把旋转前的结构记录下来；handles.x_new记录的是旋转后
Drawsupercell(hObject, eventdata, handles)
guidata(hObject, handles);


function radiobutton2_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton2,'value')==1
%     set(handles.radiobutton4, 'value', 0);
% end
guidata(hObject, handles);




function radiobutton4_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton4,'value')==1
%     set(handles.radiobutton2, 'value', 0);
% end
guidata(hObject, handles);



% --- Executes on button press in checkbox2.
function checkbox1_Callback(hObject, eventdata, handles)
if get(handles.popupmenu5,'value') == 3
    set(handles.popupmenu5,'value',2);  %如果需要分层计算，但是又旋转了lobato的系数；需要强行换为peng的
end
guidata(hObject, handles);


function edit4_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
function edit18_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Vol_Callback(hObject, eventdata, handles)
function edit_Vol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Con_Callback(hObject, eventdata, handles)
function edit_Con_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ape_Callback(hObject, eventdata, handles)
function edit_Ape_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Spr_Callback(hObject, eventdata, handles)
function edit_Spr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tilt_1_Callback(hObject, eventdata, handles)
function edit_Tilt_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tilt_2_Callback(hObject, eventdata, handles)
function edit_Tilt_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Focus_Callback(hObject, eventdata, handles)
function edit_Focus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A1_1_Callback(hObject, eventdata, handles)
function edit_A1_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A1_2_Callback(hObject, eventdata, handles)
function edit_A1_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A2_1_Callback(hObject, eventdata, handles)
function edit_A2_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A2_2_Callback(hObject, eventdata, handles)
function edit_A2_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B2_1_Callback(hObject, eventdata, handles)
function edit_B2_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B2_2_Callback(hObject, eventdata, handles)
function edit_B2_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Cs_Callback(hObject, eventdata, handles)
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
%给出scherzer条件
vol=str2num(get(handles.edit_Vol,'string'));
lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %计算波长，单位nm
str=num2str(-1.2*sqrt( 1000*str2num(get(handles.edit_Cs, 'string'))*lambda));
disp(strcat('Scherzer focus is:',str, 'nm'));

% --- Executes during object creation, after setting all properties.
function edit_Cs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A3_1_Callback(hObject, eventdata, handles)
function edit_A3_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A3_2_Callback(hObject, eventdata, handles)
function edit_A3_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S3_1_Callback(hObject, eventdata, handles)
function edit_S3_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S3_2_Callback(hObject, eventdata, handles)
function edit_S3_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A4_1_Callback(hObject, eventdata, handles)
function edit_A4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A4_2_Callback(hObject, eventdata, handles)
function edit_A4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D4_1_Callback(hObject, eventdata, handles)
function edit_D4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D4_2_Callback(hObject, eventdata, handles)
function edit_D4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B4_1_Callback(hObject, eventdata, handles)
function edit_B4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B4_2_Callback(hObject, eventdata, handles)
function edit_B4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_C5_Callback(hObject, eventdata, handles)
function edit_C5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A5_1_Callback(hObject, eventdata, handles)
function edit_A5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A5_2_Callback(hObject, eventdata, handles)
function edit_A5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D5_1_Callback(hObject, eventdata, handles)
function edit_D5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D5_2_Callback(hObject, eventdata, handles)
function edit_D5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S5_1_Callback(hObject, eventdata, handles)
function edit_S5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S5_2_Callback(hObject, eventdata, handles)
function edit_S5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit75_Callback(hObject, eventdata, handles)
function edit75_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit77_Callback(hObject, eventdata, handles)
function edit77_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------------------------------
function resul=Gaussian_focal(mtotal,rms,lambda,gmax);   %计算高斯分布的值
M=(mtotal-1)/2;
if M==0
    resul.delta_yita=0;  %离焦偏离是多少
    resul.gfsf=1;      %高斯取值是多少
else
    resul.delta_yita=energyspread(lambda, mtotal, rms, gmax);  %认为单张图由多张图叠加，确定离焦间距的具体数值
    resul.gfsf=gaussian_focal(mtotal,rms,resul.delta_yita);   %确定最优离焦条件下，高斯离散的百分比比例 
end


function f_delta=gaussian_focal(mtotal,delta,delta_yita);%the calculation of a Gaussian focal spread function
M=(mtotal-1)/2;
if M==0
    f_delta=1;
else
   f_delta=zeros(1,2*M+1);
   a=zeros(1,2*M+1);
   a=(-M:M)*delta_yita;
   f_delta=delta_yita*exp(-a.*a/(2*delta*delta))/(sqrt(2*pi)*delta);
end


function delta_yita=energyspread(lambda, mtotal, rms, gmax)   %计算最优的deltayita的值，即离焦偏离是多少
M=(mtotal-1)/2;
%yita ----- yita is the focal_step for rms of gaussian function
%pre  ----- the precision for the yita
pre=0.1; %unit nm
fieldmin=0.5;fieldmax=2.0;  %----- search yita from rms-field,unit is nm
uu=0:0.01:gmax;
fixvalue=exp(-0.5*(pi*rms*lambda).^2*uu.^4);
tempdelta=fieldmin.*rms:pre:fieldmax.*rms;
for i=1:length(tempdelta)
    gfsf=gaussian_focal(mtotal,rms,tempdelta(i));
    f_delta_value(i)=max( abs( fixvalue-real( ( gfsf*exp( -sqrt( -1 ) *pi*lambda*( ( -M:M )' )*tempdelta(i)*uu.^2) ) ) ) );   %相减之后的差值,注意这里在论文里面有录入错误，程序是正确的
end
delta_yita=tempdelta(find(f_delta_value==min(f_delta_value)));
pp=1;



function tccpara=readtccfromhandles_newSTEM(hObject, handles, lambda)  %本次是根据TEM的参数来读取的
tccpara.lambda=lambda;
%if get(handles.checkbox_polar,'Value')==1  %必定是极坐标系
   %读取说明：%参考Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;但是圈出来的系数在UM  1998 72 PP109-119以及UM 1996 64 249-264中，B2 S3 D4 B4
%S5 R5(D5)，并没有系数，因此，读取数据后，应该要分别乘以 3  4  5  5  6  6

tccpara.focus=str2num(get(handles.edit_Focus,'String'));
   
   %A1的数值不变；A1的角度是电镜显示的数值的1/2 
tccpara.A1=str2num(get(handles.edit_A1_1,'String')); tccpara.phiA1=-str2num(get(handles.edit_A1_2,'String'));
   
   %A2的数值不变；A2的角度是电镜显示的数值的1/3
   tccpara.A2=str2num(get(handles.edit_A2_1,'String')); tccpara.phiA2=-str2num(get(handles.edit_A2_2,'String'));
  
   %B2的数值不变；B2的角度与电镜显示的取负数
   tccpara.B2=str2num(get(handles.edit_B2_1,'String')); %这个参数引入时候和rew相差了3倍。例如rew里面为2000，这里必须是6000；因为程序中的系数多乘了1/3
   tccpara.phiB2=-str2num(get(handles.edit_B2_2,'String'));  
  
   %Cs的单位从um到nm
   tccpara.Cs=str2num(get(handles.edit_Cs,'String'))*(10^3); %UNIT from um to nm
   
   %A3的单位从um到nm；A3的角度与电镜显示的除以4
   tccpara.A3=str2num(get(handles.edit_A3_1,'String'))*(10^3); %UNIT from um to nm
   tccpara.phiA3=-str2num(get(handles.edit_A3_2,'String'));
   
   %S3的单位从um到nm；S3的角度与电镜显示的取负数除以2
   tccpara.S3=str2num(get(handles.edit_S3_1,'String'))*(10^3); 
   tccpara.phiS3=-str2num(get(handles.edit_S3_2,'String'));
   
   %A4的单位从um到nm；A4的角度与电镜显示的除以5
   tccpara.A4=str2num(get(handles.edit_A4_1,'String'))*(10^3);
   tccpara.phiA4=-str2num(get(handles.edit_A4_2,'String'));
   
   %D4的单位从um到nm；D4的角度与电镜显示的取负数除3
   tccpara.D4=str2num(get(handles.edit_D4_1,'String'))*(10^3); 
   tccpara.phiD4=-str2num(get(handles.edit_D4_2,'String'));
   
   %B4的单位从um到nm；B4的角度与电镜显示的一样   %需要再研究的角度关系！！！
   tccpara.B4=str2num(get(handles.edit_B4_1,'String'))*(10^3); 
   tccpara.phiB4=-str2num(get(handles.edit_B4_2,'String'));
   
    %A5的单位从mm到nm；A4的角度与电镜显示的除以6
   tccpara.A5=str2num(get(handles.edit_A5_1,'String'))*(10^6); 
   tccpara.phiA5=-str2num(get(handles.edit_A5_2,'String'));
   
   %A5的单位从mm到nm；
   tccpara.C5=str2num(get(handles.edit_C5,'String'))*(10^6);
   
   %S5的单位从mm到nm； %这个数值在hrtem里面还没有出现，因此角度关系暂时没有管
   tccpara.S5=str2num(get(handles.edit_S5_1,'String'))*(10^6); tccpara.phiS5=-str2num(get(handles.edit_S5_2,'String'));
  
   %D5的单位从mm到nm；%这个数值在hrtem里面还没有出现，因此角度关系暂时没有管
   tccpara.D5=str2num(get(handles.edit_D5_1,'String'))*(10^6); tccpara.phiD5=-str2num(get(handles.edit_D5_2,'String'));

   
   
   tccpara.A1x=tccpara.A1*cos(tccpara.phiA1/180*pi);   tccpara.A1y=tccpara.A1*sin(tccpara.phiA1/180*pi);
   %C1=focus;  focus是在输入的defocs附近有偏离的
   tccpara.A2x=tccpara.A2*cos(tccpara.phiA2/180*pi);   tccpara.A2y=tccpara.A2*sin(tccpara.phiA2/180*pi);
   tccpara.B2x=tccpara.B2*cos(tccpara.phiB2/180*pi);     tccpara.B2y=tccpara.B2*sin(tccpara.phiB2/180*pi);
   tccpara.A3x=tccpara.A3*cos(tccpara.phiA3/180*pi);   tccpara.A3y=tccpara.A3*sin(tccpara.phiA3/180*pi);
   tccpara.S3x=tccpara.S3*cos(tccpara.phiS3/180*pi);   tccpara.S3y=tccpara.S3*sin(tccpara.phiS3/180*pi);
   tccpara.C3=tccpara.Cs;
   tccpara.A4x=tccpara.A4*cos(tccpara.phiA4/180*pi);   tccpara.A4y=tccpara.A4*sin(tccpara.phiA4/180*pi);
   tccpara.D4x=tccpara.D4*cos(tccpara.phiD4/180*pi);   tccpara.D4y=tccpara.D4*sin(tccpara.phiD4/180*pi);
   tccpara.B4x=tccpara.B4*cos(tccpara.phiB4/180*pi);     tccpara.B4y=tccpara.B4*sin(tccpara.phiB4/180*pi);
   tccpara.A5x=tccpara.A5*cos(tccpara.phiA5/180*pi);   tccpara.A5y=tccpara.A5*sin(tccpara.phiA5/180*pi);
   tccpara.S5x=tccpara.S5*cos(tccpara.phiS5/180*pi);   tccpara.S5y=tccpara.S5*sin(tccpara.phiS5/180*pi);
   tccpara.C5=tccpara.C5;
   tccpara.D5x=tccpara.D5*cos(tccpara.phiD5/180*pi);   tccpara.D5y=tccpara.D5*sin(tccpara.phiD5/180*pi);

pp=1;
   
    

function x=myaperture(wave,sx,sy,kvector1,kvector2,shiftx,shifty,flag);%g_w is the radius  NOT dia
%sx、sy为两个方向倒易空间单位波矢，例如** 埃-1；kvector1为光阑的波矢长度低频，kvector2为波矢高频，例如10 埃-1；shiftx为平移的波矢，例如从中心处移动一点点,flag表示操作方式，取kvector1到kvector2之间的信息，还是不取这之间的信息。
%例如低频滤波：myaperture(newdd,1,1,0,kvector,0,0,1);
%例如取带通：myaperture(newdd,1,1,kvector1,kvector,0,0,1);
%例如取高频：myaperture(newdd,1,1,0,kvector,0,0,0);
tx=shiftx/sx;
ty=shifty/sy;
 
[m,n]=size(wave);
uu=-round((n-1)/2):round(n/2)-1;
vv=-round((m-1)/2):round(m/2)-1;
[uu,vv]=meshgrid((uu-ty).*sy,(vv-tx).*sx);
uuvv=sqrt(uu.*uu+vv.*vv);
clear uu
clear vv
if flag==1
   wave(find(uuvv<kvector1 | uuvv>kvector2))=0;    %  滤波
else
   wave(find(uuvv>=kvector1 & uuvv<=kvector2))=0;    %  滤波
end
x=wave;
return;

function flag = getparaflag(handles)
%得到计算参数的选择，分别是p，l和n三种
num = get(handles.popupmenu5,'value');
if num == 1
    flag = 'p';
elseif num == 2;
    flag = 'n';
elseif num ==3;
    flag = 'l';
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
load allresul;
eachthick = allresul.eachthick;
handles.eachthick= eachthick;
slicethick = allresul.slicethick;

%增加中间层的描述
display(strcat('Thickness of eachslice is ', num2str(eachthick),' angstrom, there are totally ', num2str(length(slicethick)), ' slices'));

handles.mid_slice_num = sort( str2num(get(handles.edit105,'string')));  % the mid slice's output
handles.mid_slice_num(find(handles.mid_slice_num>=length(slicethick))) = [];
display(strcat( 'Output the wave at the middle slices:', num2str(handles.mid_slice_num)))

paraflag = getparaflag(handles); 
%set(handles.pushbutton3,'enable','off');
%if paraflag == 'p'|| paraflag == 'l'  %彭老师公式 20210223
    [all_nuclear, all_nuc_ion, absorp_n, absorp_n_i, corr_peng_nuc, corr_peng_nuc_ion]=CommonPart_forsimulation_lobato_peng(hObject, eventdata, handles,paraflag);  %step0-step3移动到这里
%elseif paraflag == 'n' %彭老师公式加修正
     %如果是peng的修正，会得到corr相关的两个矩阵。计算思路，是把这两个矩阵计算势场，之后作为原势场的初值
%end
    
%外层循环是所有原子的坐标随机变化30次
%内层循环：是否考虑分层效应，如果是，则需要多搭几层，再计算每层的贡献。

%all_nuclear 第一列是元素种类，第2-4列是坐标，第5列是B因子，第6列是占据率
%step 4
%提取一些常数，与电压、抽样率无关
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
handles.lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %计算波长，单位nm
handles.vol=vol*1000;

%step 5
%修改2020.04.29
if get(handles.radiobutton11,'value') ||  get(handles.radiobutton10,'value') ||  get(handles.radiobutton13,'value')  %stem像 or idpc or CBED像 
   probesx=str2num(get(handles.edit4,'String'));
   [para_part1, U, V]=CommonPara_TEM(handles, probesx, probesx);  %画CTF和PhasePlate用固定的大小即可  
   para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码
    %参数分为两个部分，lambda需要都有
   residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %计算phase plate visualing the residual aberrations，kai'
    %-------------------------------------
   mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,2);
   
   myap=ones(probesx);  %构造光阑函数，kmax是外面定的
   myap=myaperture(myap,1/(para_part1.sampling.*probesx),1/(para_part1.sampling.*probesx),0,para_part1.gmax*0.01/para_part1.lambda,0,0,1); %如果10mrad，读入就是1；再乘以0.01，换到到rad，再除以nm单位的波长
   num=length(mytcc(1,1,:));  %考察有多少个probe需要扫过图像
   for i=1:num  %光阑在中间，mytcc的重心移边上了。
        guiyihua=sum(sum(abs(ifft2(ifftshift(myap.*mytcc(:,:,i)))).^2));  %20201226 add one aper to normalize 
        guiyihua=sqrt(guiyihua);  %20201231需要求开根号，才能让probe满足3.59和3.68式
        handles.probe(:,:,i) =  fftshift(ifft2(ifftshift(myap.*mytcc(:,:,i))))./guiyihua;
   end
   handles.gfsf=para_part1.gfsf;
end
%__________修改2020.04.29 主要删除原有的STEM后缀的那些程序
if get(handles.radiobutton9,'value')  %HRTEM代码  翻译代码
    sy=handles.green_Nrow;sx=handles.green_Ncol;   %绿色框格的尺寸
    [para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %画CTF和PhasePlate用固定的大小即可
    U=U-para_part1.tiltx;  %add 20210119 the incident wave is tilted but not the CTF is tilted
    V=V-para_part1.tilty;  %the tilted CTF is shown for CTF display. but in simulation, only a tilted incident beam is required.
 
    para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码
    %参数分为两个部分，lambda需要都有
    residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %计算phase plate visualing the residual aberrations
    %-------------------------------------
    mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,1);

    myap=ones(size(U));
    myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);  
    
    handles.gfsf=para_part1.gfsf;
end


%step 6  成像
%如果考虑原子实际的vibration，需要生成很多的坐标
if get(handles.radiobutton2, 'value')   %if vibration is required
    allvib = str2num(get(handles.edit99,'string'));
else
    allvib = 1;
end

for vib = 1 : allvib;  %vibration
    %随机产生原子的坐标偏离,需要与每个原子的Debye有关；之后叠加到原子的坐标上
    all_nuclear_copy=all_nuclear;
    all_nuc_ion_copy=all_nuc_ion;
    if paraflag=='n'
        corr_nuclear_copy = corr_peng_nuc;
        corr_peng_nuc_ion_copy = corr_peng_nuc_ion;
    end
        
    if ~isempty(all_nuclear) & get(handles.radiobutton2, 'value') % if vibration is required
        % make random positions
        % satisfy the condition that the atoms before special slices will
        % not be shifted to its slices
        % input atoms' parameters, thickness of each slice, the thickness
        % of the top, the slice number of middle slice.
        all_nuclear_copy = randshiftpos(all_nuclear_copy, eachthick, slicethick, handles.mid_slice_num); 
        if paraflag == 'n'
            if ~isempty(corr_peng_nuc)  %如果需要修正势场，把修正部分的参数的原子位置更新一下
                corr_nuclear_copy = corr_peng_nuc;
                corr_nuclear_copy(:,2:4) = all_nuclear_copy(:,2:4);
            end
        end
    end
    if ~isempty(all_nuc_ion) & get(handles.radiobutton2, 'value')
        all_nuc_ion_copy = randshiftpos(all_nuc_ion_copy, eachthick, slicethick, handles.mid_slice_num);
        if paraflag == 'n'
            if ~isempty(corr_peng_nuc_ion) %如果需要修正势场，把修正部分的参数的原子位置更新一下
                corr_peng_nuc_ion_copy = corr_peng_nuc_ion;
                corr_peng_nuc_ion_copy(:,2:4) = all_nuc_ion_copy(:,2:4);
            end
        end
    end
    
    %需要或者不需要计算高度带来的差别,并且还需要清除没有落在范围内的原子。    ~~~~~~~~~~~~~~~~~~~~~~~
    %improved multi-slice method
    flag=get(handles.checkbox1, 'value');  %是否需要考虑原子划分在不同的层
    if flag==1
        multiceng=str2num(get(handles.edit98, 'string'));
        flag=multiceng;
        if paraflag == 'l'
            disp('Lobato parameters cannot be sliced into multiple slices in this program')
            flag=0;
        end
    end
    %增加两种参数的选择，peng和lobato的
%     [ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i] = ...
%              CalAllEquilent20201106(all_nuclear_copy, all_nuc_ion_copy, absorp_n, absorp_n_i, eachthick, slicethick, flag);  
    DBmode=double(get(handles.radiobutton4,'value')|get(handles.radiobutton14,'value'));  %任何一个数是1，就是1
    [ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i] = ...
             CalAllEquilent20210223_lobato_peng(all_nuclear_copy, all_nuc_ion_copy, absorp_n, absorp_n_i, eachthick, slicethick, flag, paraflag, DBmode);  
    
    ele_n_corr=[]; ele_n_i_corr=[]; corr_info=[]; series_n_corr=[]; series_n_i_corr=[];
    if paraflag == 'n'
        [ele_n_corr, ele_n_i_corr, corr_info, series_n_corr, series_n_i_corr] =  ...
            CalAllEquilent20210305_peng_corr(corr_nuclear_copy, corr_peng_nuc_ion_copy, eachthick, slicethick, flag);   
    end   
    
    %得到第vib张图像  ~~~~~~~~~~~~~~~~~~~~
    if vib==1
         if get(handles.radiobutton11,'value') || get(handles.radiobutton13,'value')  %stem像
             [myresul, midresul] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton10,'value')  %cbed像, 
             myresul=CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton9,'value')  %hrtem像
             myresul=HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
             
         end
     else
         vib
         if get(handles.radiobutton11,'value') || get(handles.radiobutton13,'value')  %stem像
             [myresultemp, midresultemp] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
             myresul=myresul+myresultemp;
             midresul = midresul+midresultemp; 
         end
         if get(handles.radiobutton10,'value')  %cbed像
             myresul=myresul+CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i,ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton9,'value')  %hrtem像
             myresul=myresul+HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc, ele_n_corr, ele_n_i_corr, corr_info,series_n_corr, series_n_i_corr);
         end
      end
end

disp(strcat('Results will be saved in the folder:', handles.saveresult));
nname = get(handles.edit106,'string');

allstemnum = 1; %所有stem图像的个数
steminfo = [];
if get(handles.radiobutton11,'value')  %stem像 需要重新reshape矩阵大小
    aper_dect = str2num(get(handles.edit79,'string'));%制造接收光阑 

   %中间结果 & 最后结果
   if isempty(midresul)
       midresul=myresul(:);
   else
      midresul(:,end+1)=myresul(:);
   end
   handles.mid_slice_num = [handles.mid_slice_num, length(slicethick)];
   midimg_result_temp = zeros( handles.width_red , handles.hight_red, length(myresul(1,:)));
   for i=1:length(handles.mid_slice_num)
       midimg_result_temp(:) = midresul(:,i)/probesx/probesx/allvib;
              
       for j=1:length(myresul(1,:))
            figure;imshow(midimg_result_temp(:,:,j).', 'XData',...
           str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
           'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
          'DisplayRange',[]);axis on;xlabel('Angstrom');
      
            distitle = strcat('STEM images (thickness:', ...
                num2str(handles.mid_slice_num(i)*eachthick) , ...
                '; Dectector: ', ...
                num2str(aper_dect(j,1)), ...
                'to ', num2str(aper_dect(j,2))   );colorbar
            title( distitle )
            
            discontent = strcat('Mid result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), '_Float Data in',...
                strcat(handles.saveresult, '\',nname ,'_', 'STEM_', num2str(handles.mid_slice_num(i)*eachthick), ...
                'slice_Dect_',num2str(round(aper_dect(j,1))),'_',num2str(round(aper_dect(j,2))),'.dat'))
            disp(discontent);
            
            fname = strcat(handles.saveresult, '\',nname ,'_', 'STEM_', num2str(handles.mid_slice_num(i)*eachthick), ...
                'slice_Dect_',num2str(round(aper_dect(j,1))),'_',num2str(round(aper_dect(j,2))),'.dat')
            fid = fopen(fname,'w');
            fwrite(fid, midimg_result_temp(:,:,j), 'float'); 
            fclose(fid);
            
            steminfo(allstemnum).fname = fname;
            steminfo(allstemnum).discontent = discontent;
            steminfo(allstemnum).title = distitle;
            allstemnum = allstemnum+1;
       end
   end
   
end

if get(handles.radiobutton9,'value')  %hrtem像 需要重新reshape矩阵大小
    for i=1:length(myresul(1,1,:))
      figure;imshow(myresul(:,:,i)./allvib, 'XData',...
           str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
           'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
          'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar
      
      disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
          '_Float Data in', strcat(handles.saveresult, '\', nname ,'_', 'HRTEM','.dat')));
      fid = fopen(strcat(handles.saveresult, '\', nname ,'_', 'HRTEM','.dat'),'w');
      fwrite(fid, myresul(:,:,i).'/allvib, 'float'); 
      fclose(fid);
    end
end
if get(handles.radiobutton10,'value')  %CBED像 需要重新reshape矩阵大小
    for i=1:length(myresul(1,1,:))
      figure;imshow(myresul(:,:,i)./allvib,[])
      
      disp(strcat('Result _W', num2str(handles.CBEDprobesx),' * H', num2str(handles.CBEDprobesx), ...
          '_Float Data in', strcat(handles.saveresult, '\', nname ,'_', 'CBED','.dat')));
      fid = fopen(strcat(handles.saveresult, '\', nname ,'_', 'CBED','.dat'),'w');
      fwrite(fid, myresul(:,:,i).'/allvib, 'float'); 
      fclose(fid);
    end
end

% step=2;
% sampling=0.1;
% [x,y] = meshgrid( ([1:newy]-mycen(newy))/(newy*sampling*step/handles.conv_sampling) ,([1:newx]-mycen(newx))/(newx*sampling*step/handles.conv_sampling));
%         ss = 0.7/(sqrt(-log(0.5))*2);
%Kepeng
if get(handles.radiobutton13,'value')  %idpc像
    aper_dect = str2num(get(handles.edit79,'string'));%制造接收光阑 

    [kx,ky]=meshgrid((-handles.width_red/2:handles.width_red/2-1)./(handles.width_red.*para_part1.sampling*str2num(get(handles.edit77,'string'))),(-handles.hight_red/2:handles.hight_red/2-1)/(handles.hight_red.*para_part1.sampling*str2num(get(handles.edit77,'string'))));
    k2=kx.^2+ky.^2;  %数值范围不一样，要检查
    ss = 0.1/(sqrt(-log(0.5))*2);

   %中间结果 & 最后结果
     if isempty(midresul)
         midresul=myresul(:);
     else
         midresul(:,end+1)=myresul(:);
     end
     handles.mid_slice_num = [handles.mid_slice_num, length(slicethick)];
     midimg_result_temp = zeros( handles.width_red,handles.hight_red, length(myresul(1,:)));
     for i=1:length(handles.mid_slice_num)
          midimg_result_temp(:) = midresul(:,i);
          
          ta=zeros(handles.width_red,handles.hight_red);    
          for j=1:(length(myresul(1,:))/4)
              ta(:)=midimg_result_temp(:,:,4*(j-1)+1)/probesx/probesx/allvib;
              aa(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+2)/probesx/probesx/allvib;
              bb(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+3)/probesx/probesx/allvib;
              cc(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+4)/probesx/probesx/allvib;
              dd(:,:,j)=ta';
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+1)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));aa(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+2)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));bb(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+3)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));cc(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+4)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));dd(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
       
        
%                 px=(fftshift(fft2(bb(:,:,j)-dd(:,:,j))).*kx)*(2*pi*sqrt(-1)); 
%                 px(handles.hight_red/2+1,handles.width_red/2+1)=mean(mean(bb(:,:,j)-dd(:,:,j)));
%                 py=(fftshift(fft2(cc(:,:,j)-aa(:,:,j))).*ky)*(2*pi*sqrt(-1)); 
%                 py(handles.hight_red/2+1,handles.width_red/2+1)=mean(mean(cc(:,:,j)-aa(:,:,j)));
%                 px=real( ifft2(ifftshift(px)) );
%                 py=real( ifft2(ifftshift(py)) );
                px = bb(:,:,j)-dd(:,:,j);
                py = cc(:,:,j)-aa(:,:,j);
                
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'aa_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, aa.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'bb_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, bb.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'cc_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, cc.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'dd_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, dd.', 'float'); 
                fclose(fid);
                
                
                
                figure;subplot(2,2,1);imshow(px, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('x vector');subplot(2,2,2);imshow(py, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('y vector');
                [sx,sy]=size(px);newp=zeros(sx,sy,3);newp(:,:,1)=px;newp(:,:,2)=py;
                subplot(2,2,3);image(newp*255);axis on;xlabel('Angstrom');
                subplot(2,2,4);imshow(-(px).^2+(py).^2, 'XData',...
                       str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                       'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                       'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('x^2+y^2');    
                title( strcat('DDPC (thick:', ...
                      num2str(handles.mid_slice_num(i)*eachthick) , ...
                      '; Dectector: ', ...
                      num2str(aper_dect(j,1)), ...
                     'to ', num2str(aper_dect(j,2))   ))
                disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red),...
                    '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'DDPC_', num2str(handles.mid_slice_num(i)),'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat')));
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'DDPCx_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, px.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'DDPCy_', num2str(handles.mid_slice_num(i)),'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, py.', 'float'); 
                fclose(fid);
                 
                 
                p=(fftshift(fft2(px)).*kx+fftshift(fft2(py)).*ky)./(2*pi*sqrt(-1)*k2); 
                p(handles.hight_red/2+1,handles.width_red/2+1)=(mean(mean(px))+mean(mean(py)));
                p=real( ifft2(ifftshift(p)) );
                
                figure;imshow(p, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar;
                title( strcat('IDPC (thick:', ...
                      num2str(handles.mid_slice_num(i)*eachthick) , ...
                      '; Dectector: ', ...
                      num2str(aper_dect(j,1)), ...
                     'to ', num2str(aper_dect(j,2))   ))
                disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
                     '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'IDPC_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat')));
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'IDPC_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, p.', 'float'); 
                fclose(fid);
            
%                 figure;imshow((aa(:,:,j)+cc(:,:,j)+bb(:,:,j)+dd(:,:,j)), 'XData',...
%                        str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
%                        'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
%                        'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar;
%               
%                 title( strcat('Traditional STEM (thick:', ...
%                       num2str(handles.mid_slice_num(i)*eachthick) , ...
%                       '; Dectector: ', ...
%                       num2str(aper_dect(j,1)), ...
%                      'to ', num2str(aper_dect(j,2))   ))
%                  disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
%                       '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'IDPC_STEM', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_', num2str(round(aper_dect(j,1))), ...
%                       '_', num2str(round(aper_dect(j,2))),'.dat')));
%                  fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'IDPC_STEM_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
%                         '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
%                  fwrite(fid, aa(:,:,j)+cc(:,:,j)+bb(:,:,j)+dd(:,:,j), 'float'); 
%                  fclose(fid);
          end
   end

end
    
 
%Kepeng
% handles.steminfo = steminfo; %save all file information to do convolution
% guidata(hObject, handles);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%KepengWu1124   begin
function APER =makeAPER(handles);
aper_dect = str2num(get(handles.edit79,'string'));%制造接收光阑
[px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
    (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
p2_probe=(px.^2+py.^2);
aper_range = aper_dect*0.001/(handles.lambda*10); %1/A 单位
num=length(aper_dect(:,1));
APER=zeros(handles.probesx, handles.probesx,num);

[theta, rho] =  cart2pol(px,py);  % range of theta is from -pi to pi
for i=1:num  %构建接收光阑
    temp=zeros(handles.probesx, handles.probesx);
    temp(find(p2_probe>=aper_range(i,1).^2 & p2_probe<=aper_range(i,2).^2 )) = 1;

    if aper_dect(i,3) >180  %把0-360的书写范围转变到-180 到180
        aper_dect(i,3) = aper_dect(i,3)-180; 
    end
    if aper_dect(i,3) >= -180 & aper_dect(i,3) < -90
        dectvalue = aper_dect(i,3);
    elseif aper_dect(i,3) >= -90 & aper_dect(i,3) < 0
        dectvalue = aper_dect(i,3)-90;
    elseif aper_dect(i,3) >= 0 & aper_dect(i,3) < 90
        dectvalue = aper_dect(i,3)-180;
    elseif aper_dect(i,3) >= 90 & aper_dect(i,3) < 180
        dectvalue = aper_dect(i,3)-270;
    end
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= dectvalue*pi/180 & theta<(dectvalue+90)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+1) = tempaper;
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= (dectvalue+90)*pi/180 & theta<(dectvalue+180)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+2) = tempaper;
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= (dectvalue+180)*pi/180 & theta<(dectvalue+270)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+3) = tempaper;
    
    APER(:,:,4*(i-1)+4) = 1- APER(:,:,4*(i-1)+1) -APER(:,:,4*(i-1)+2) - APER(:,:,4*(i-1)+3);

    APER(:,:,4*(i-1)+1) = temp.* APER(:,:,4*(i-1)+1);  %加上带通了
    APER(:,:,4*(i-1)+2) = temp.* APER(:,:,4*(i-1)+2);
    APER(:,:,4*(i-1)+3) = temp.* APER(:,:,4*(i-1)+3);
    APER(:,:,4*(i-1)+4) = temp.* APER(:,:,4*(i-1)+4);
    % APER(:,:,4*(i-1)+1)=temp;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  KepengWu1124 end
function myresul=HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc,...
    ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr)
%需要计算green的尺寸，以及所有原子的相对位置，计算的区域，以及拼凑回原图时候的坐标
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom单位
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %单位与书一样，/kV*A

%APERTURE  % 光阑矩阵，保证传播不会产生wrap效应
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
%重大改进，20210225，否则高频的晶格信息丢失了
clear tx ty

paraflag = getparaflag(handles);
if get(handles.GPURB,'value')  %使用GPU来计算
    %GPU必须是方形
    %%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布
    gncol = handles.green_Ncol; gnrow = handles.green_Nrow;
    green=max(handles.green_Nrow, handles.green_Ncol);
    handles.green_Ncol=green; handles.green_Nrow=green;
    [gx_green,gy_green]=meshgrid((-green/2:(green/2-1))./(green*handles.sampling), ...
        (-green/2:(green/2-1))./(green*handles.sampling)); %REPRO单位是1/A)^2,不是(1/nm)^2
    sx_green=gx_green/2;  % peng 的s参数
    sy_green=gy_green/2;
    s2_green=sx_green.^2+sy_green.^2;
    APERTURE=ones(green, green);
    
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%赋初值
        for i = 1:length(corr_info(:,1))  %把每种的修正势场都算出来，种类很少，所以算与元素个数相当的矩阵就好
            r2 = s2_green;
            r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
            r2(find(r2<corr_info(i,5).^2))=0;
            r2(find(r2>=corr_info(i,6).^2))=1;
            corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
        end  %gpu,输入这个矩阵更方便一些
    end
    
    handles.probe=[];
    APER=[];
    psf_fft=[];
    handles.probe_ingreenNrow=[]; handles.probe_ingreenNcol=[]; handles.probestep=[];
    tic
    cuda_stem_potential;
    toc
    
    handles.green_Ncol=gncol; handles.green_Nrow=gnrow;
    potential_temp = zeros(green,green, max(length(series_n), length(series_n_i)) );
    potential_temp(:) = potentialx(:) + sqrt(-1)*potentialy(:);
    
    for i=1:length(potential_temp(1,1,:)); %需要转置
        potential_temp(:,:,i) = potential_temp(:,:,i).';
    end
    potential = potential_temp(1:handles.green_Ncol,1:handles.green_Nrow, :);
    clear potential_temp;
end
%CPU可以是矩形
    %%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
        (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO单位是1/A)^2,不是(1/nm)^2
    sx_green=gx_green/2;  % peng 的s参数
    sy_green=gy_green/2;
    s2_green=sx_green.^2+sy_green.^2;
    

if get(handles.CPURB,'value')  %使用CPU来计算 
    
%计算每层的势场
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
%     ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
%     series_n, series_n_i, ...  %原子排列次序   
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
potential=GetPotential4AllSlice_multicore_lotabo_peng_corr(handles.green_Ncol, handles.green_Nrow,... 
    ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收  %本函数调用修正peng的参数
    ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
    series_n, series_n_i, ...  %原子排列次序   
    ele_n_corr, ele_n_i_corr, corr_info, ... %如果修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr,...
    s2_green, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。
end


%成像光阑
myap=ones(size(gx_green));
myap=myaperture(myap,1/(handles.sampling.*handles.green_Ncol),1/(handles.sampling.*handles.green_Nrow),0,str2num(get(handles.edit_Ape,'string'))*0.1,0,0,1);  %埃的单位

%像的初始化
myinten=zeros(size(gx_green));

%add 20210119
Vol=str2num(get(handles.edit_Vol,'String')); %电压
lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %计算波长，单位nm
tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);
phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % 弧度单位 例如 10mrad  
tiltx=tilt/lambda*cos(phitilt/180*pi); %换算到倒易格式的单位 1/nm
tilty=tilt/lambda*sin(phitilt/180*pi);
tiltxnum=tiltx*handles.sampling*0.1*handles.green_Ncol;  %nm和sampling 埃的换算
tiltynum=tilty*handles.sampling*0.1*handles.green_Ncol;

%波函数传播  %20210119 入射倾斜的波函数----------------

%   波函数propogation;与zheight有关
p2=((gx_green+tiltx*0.1).^2+(gy_green+tilty*0.1).^2);  %传播子
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2);   %都在单位A计算。
psf_fft(find(p2>1/(16*handles.sampling*handles.sampling)))=0;  %加入半光阑
psf_fft=ifftshift(psf_fft);

myfftwave = zeros(size(gx_green)); 
if tilt > 0
     myfftwave(mycen(handles.green_Ncol)-round(tiltynum), mycen(handles.green_Nrow)-round(tiltxnum)) = 1;
else
     myfftwave(mycen(handles.green_Ncol), mycen(handles.green_Nrow)) = 1;  %赋值为平面波
end
   %this wave shift is the same as the U*U+V*V before about line 885 U=U-para_part1.tiltx;   
mywave=ifft2(ifftshift(myfftwave))*handles.green_Nrow*handles.green_Ncol;  %读入波函数
%----------end 20210119

for k=1:length(potential(1,1,:));
    if rem(k,10)==0
        pp=1;
    end
        mywave =  mywave.*potential(:,:,k);
        mywave=fft2(mywave);   %注意，节省掉一个fftshift
        mywave=mywave.*psf_fft;
        mywave=ifft2(mywave);
end
%存波函数-----add for Yao fenfa
[sx, sy]= size(gx_green);
disp(strcat('Wave saved in HRTEM_exitwave.dat (complex 8), H',num2str(sx), '*W',num2str(sy)))
a=zeros(2*sx,sy);
a(1:2:end,:)=real(mywave);
a(2:2:end,:)=imag(mywave);
fid = fopen('HRTEM_exitwave.dat','w');
fwrite(fid, a, 'float')
fclose(fid);

mywave=fftshift(fft2(mywave));

    
for nn=1:length(handles.gfsf);
    myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn).*myap))).^2;
end
myresul=myinten(handles.HRTEM_ingreenNcol: handles.HRTEM_ingreenNcol+handles.hight_red-1, ...
                handles.HRTEM_ingreenNrow: handles.HRTEM_ingreenNrow+handles.width_red-1);
return;


function myresul=CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
%需要计算green的尺寸，以及所有原子的相对位置，计算的区域，以及拼凑回原图时候的坐标
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom单位
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %单位与书一样，/kV*A

%APERTURE  % 光阑矩阵，保证传播不会产生wrap效应
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
clear tx ty

%%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
    (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO单位是1/A)^2,不是(1/nm)^2
sx_green=gx_green/2;  % peng 的s参数
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;



%   波函数propogation;与zheight有关
[px, py]=meshgrid( (-handles.CBEDprobesx/2:(handles.CBEDprobesx/2-1))./(handles.CBEDprobesx*handles.sampling), ...
    (-handles.CBEDprobesx/2:(handles.CBEDprobesx/2-1))./(handles.CBEDprobesx*handles.sampling));
p2_probe=(px.^2+py.^2);  %传播子
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %都在单位A计算。
%psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %加入半光阑
psf_fft(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %加入2/3的光阑
psf_fft=ifftshift(psf_fft);

paraflag = getparaflag(handles);
if get(handles.GPURB,'value')  %使用GPU来计算
    %GPU必须是方形
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%赋初值
        for i = 1:length(corr_info(:,1))  %把每种的修正势场都算出来，种类很少，所以算与元素个数相当的矩阵就好
            r2 = s2_green;
            r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
            r2(find(r2<corr_info(i,5).^2))=0;
            r2(find(r2>=corr_info(i,6).^2))=1;
            corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
        end  %gpu,输入这个矩阵更方便一些
    end
    
    APERTURE=ones(handles.green_Nrow, handles.green_Ncol);
    
    APER=[];
    handles.probe_ingreenNrow=[]; handles.probe_ingreenNcol=[]; handles.probestep=[];
    tic
    cuda_stem_potential;
    toc

    potential_temp = zeros(handles.green_Nrow,handles.green_Ncol, max(length(series_n), length(series_n_i)) );
    potential_temp(:) = potentialx(:) + sqrt(-1)*potentialy(:);
    
    for i=1:length(potential_temp(1,1,:)); %需要转置
        potential(:,:,i) = potential_temp(:,:,i).';
    end
    clear potential_temp
end
%CPU可以是矩形
    %%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布


if get(handles.CPURB,'value')  %使用CPU来计算 
    
%计算每层的势场
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
%     ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
%     series_n, series_n_i, ...  %原子排列次序   
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
potential=GetPotential4AllSlice_multicore_lobato_peng_corr(handles.green_Ncol, handles.green_Nrow,... 
    ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收  %本函数调用修正peng的参数
    ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
    series_n, series_n_i, ...  %原子排列次序   
    ele_n_corr, ele_n_i_corr, corr_info, ... %如果修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr,...
    s2_green, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。
end

%计算每层的势场 %使用多线程计算
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
%     ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
%     series_n, series_n_i, ...  %原子排列次序
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。

% %计算每层的势场
% potential=GetPotential4AllSlice(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
%     ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
%     series_n, series_n_i, ...  %原子排列次序
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE);   %为了STEM计算不出错，这里只带入到HRTEM和CBED。

%成像光阑
myap=ones(handles.CBEDprobesx);
myap(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %16是加入半光阑 %9是加入2/3的光阑

%传播
myresul=zeros(handles.CBEDprobesx);
[probesx, probesy]=size(handles.probe);
for nn=1:length(handles.gfsf);
    mywave=handles.probe(:,:,nn);  %读入波函数
    for k=1:length(potential(1,1,:));
        mywave =  mywave.*potential(handles.CBED_ingreenNrow:handles.CBED_ingreenNrow+handles.CBEDprobesx-1, ...
                   handles.CBED_ingreenNcol:handles.CBED_ingreenNcol+handles.CBEDprobesx-1,k);
        mywave=fft2(mywave);   %注意，节省掉一个fftshift
        mywave=mywave.*psf_fft;
        mywave=ifft2(mywave);
    end
    mywave=fftshift(fft2(mywave));
    myresul=myresul+handles.gfsf(nn)*(abs(mywave).^2);
end
%myresul=myresul.*myap;  %加上光阑。
return;

 

        
function [myresul,mid_ceng_mat] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
    series_n_corr, series_n_i_corr);
%需要计算green的尺寸，以及所有原子的相对位置，计算的区域，以及拼凑回原图时候的坐标
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
%#############从这里开始，gpu和cpu都只对方形的形状进行处理，否则会出错
    if handles.green_Nrow<handles.green_Ncol handles.green_Nrow=handles.green_Ncol;else handles.green_Ncol=handles.green_Nrow; end   %20201107
%#################################################################

PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom单位
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %单位与书一样，/kV*A
% %测试：与kirkland书一样图像，图5-2
% vol=10:10:1000;
% lambda=1.0e+9*h.*c./sqrt(e.*vol*1000.*(2*me*c*c+e.*vol*1000));
% sigma=2.*pi./(lambda.*10.*vol).*(me*c*c+e.*vol*1000)./(2*me*c*c+e.*vol*1000);

%APERTURE  % 光阑矩阵，保证传播不会产生wrap效应
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
clear tx ty

%%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
    (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO单位是1/A)^2,不是(1/nm)^2
sx_green=gx_green/2;  % peng 的s参数
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;

% 
%   probe_propogation;与zheight有关 old codes20210119
% [px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
%     (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
% p2_probe=(px.^2+py.^2);  %传播子
% psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %都在单位A计算。
% psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %加入半光阑
% psf_fft=ifftshift(psf_fft);
%add 20210119
Vol=str2num(get(handles.edit_Vol,'String')); %电压
lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %计算波长，单位nm
tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);
phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % 弧度单位 例如 10mrad  
tiltx=tilt/lambda*cos(phitilt/180*pi); %换算到倒易格式的单位 1/nm
tilty=tilt/lambda*sin(phitilt/180*pi);

%20210119 入射倾斜的波函数在多层法传播时的注意点----------------
%   波函数propogation;与zheight有关
[px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
    (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
p2_probe=((px+tiltx*0.1).^2+(py+tilty*0.1).^2);  %传播子
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %都在单位A计算。
%psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %加入半光阑
psf_fft(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %加入2/3的光阑
psf_fft=ifftshift(psf_fft);

%end 20210119



aper_dect = str2num(get(handles.edit79,'string'));%制造接收光阑
aper_range = aper_dect*0.001/(handles.lambda*10); %1/A 单位
%计算有几个接收光阑

%Kepeng
if get(handles.radiobutton11,'value')  %stem像
    num=length(aper_dect(:,1));
    APER=zeros(handles.probesx, handles.probesx,num);
    for i=1:num  %构建接收光阑
        temp=zeros(handles.probesx, handles.probesx);
        temp(find((px.^2+py.^2)>=min(aper_range(i,:)).^2 & (px.^2+py.^2)<=max(aper_range(i,:)).^2 )) = 1;
        APER(:,:,i)=temp;
    end
end

if get(handles.radiobutton13,'value')
    APER=makeAPER(handles);
end

% if 1  %3d stem
%     [dx,dy] = meshgrid( [1:handles.probesx]-mycen(handles.probesx) ,[1:handles.probesx]-mycen(handles.probesx) );
%     APER(:,:,end+1) = dx.*sum(APER,3); APER(:,:,end+1) = dy.*sum(APER(:,:,1:end-1),3);
% end
%Kepeng

% 
% myresul1 = cuda_STEM_core_complex(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%        handles.width_red, handles.hight_red, handles.probestep,APER); 

paraflag = getparaflag(handles); %which parameter will be used for the scatting factor
if get(handles.CPURB, 'value') 
%     [myresul,mid_ceng_mat] = cuda_STEM_CPU_Mcore_lobato_peng(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%         handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%         handles.width_red, handles.hight_red, handles.probestep,APER, paraflag, handles.mid_slice_num);
    [myresul,mid_ceng_mat] = cuda_STEM_CPU_Mcore_lobato_peng_corr(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, ...
        ele_n_corr, ele_n_i_corr, corr_info, ... %增加修正peng的系数，带入的修正势场的量
        series_n_corr, series_n_i_corr,...
        gx_green,gy_green, s2_green,APERTURE, ... ...
        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
        handles.width_red, handles.hight_red, handles.probestep,APER, paraflag, handles.mid_slice_num);
    
    
    % myresul1 = cuda_STEM_core_complex(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%        handles.width_red, handles.hight_red, handles.probestep,APER); 
end
if get(handles.GPURB, 'value')  %如果使用GPU计算
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%赋初值
%         for i = 1:length(corr_info(:,1))  %把每种的修正势场都算出来，种类很少，所以算与元素个数相当的矩阵就好
%             r2 = s2_green;
%             r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
%             r2(find(r2<corr_info(i,5).^2))=0;
%             r2(find(r2>=corr_info(i,6).^2))=1;
%             corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
%         end  %gpu,输入这个矩阵更方便一些
        for i = 1:length(corr_info(:,1))  %把每种的修正势场都算出来，种类很少，所以算与元素个数相当的矩阵就好
            s2 = s2_green;
            g2 = 4*s2;
            corr_info_matrix(:,:,i) =(corr_info(i,11)*(2+corr_info(i,12)*g2)./(1+corr_info(i,12).*g2).^2 + ...
                                  corr_info(i,13)*(2+corr_info(i,14)*g2)./(1+corr_info(i,14).*g2).^2 + ...
                                  corr_info(i,15)*(2+corr_info(i,16)*g2)./(1+corr_info(i,16).*g2).^2 + ...
                                  corr_info(i,17)*(2+corr_info(i,18)*g2)./(1+corr_info(i,18).*g2).^2 + ...
                                  corr_info(i,19)*(2+corr_info(i,20)*g2)./(1+corr_info(i,20).*g2).^2 ) -...
                                 (corr_info(i,1).*exp(-s2.*corr_info(i,2)) ...
                                 + corr_info(i,3).*exp(-s2.*corr_info(i,4)) ...
                                 + corr_info(i,5).*exp(-s2.*corr_info(i,6)) ...
                                 + corr_info(i,7).*exp(-s2.*corr_info(i,8)) ...
                                 + corr_info(i,9).*exp(-s2.*corr_info(i,10)));                           
           end  %gpu,输入这个矩阵更方便一些
    end
    %cuda计算，需要输入的变量是 corr_info_matrix; ele_n_corr, ele_n_i_corr;series_n_corr, series_n_i_corr
    tic
    cuda_stem;
    toc
end
 %save myresul myresul
 %stop;




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
if get(handles.radiobutton11, 'value') | get(handles.radiobutton13, 'value')   %如果是STEM和IDPC的计算，显示相应的区域
    handles.sampling=str2num(get(handles.edit75,'string'));  %埃/pixel   %这个是传播矩阵和投影势函数矩阵的抽样率，不是probe的扫描抽样率
    if ~isempty ( strfind(get(handles.edit77,'string'), '.'))
        msgbox('''Scan Step'' Must be an integer');
        return;
    end
    handles.probestepsampling=str2num(get(handles.edit77,'string')) * handles.sampling;  %表示扫描时候的分辨率
    handles.probesx=str2num(get(handles.edit4,'string'));  %读取probe的尺寸
    if rem(handles.probesx,2)==1 handles.probesx=handles.probesx+1;end
    
    %设定扫描区域topleft-rightdown的坐标
        %红色框的宽度和高度分别为多少像素
    width_red=str2num(get(handles.edit20,'string')); hight_red=str2num(get(handles.edit19,'string'));
    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    rd_red=tl_red+[width_red, hight_red]*handles.probestepsampling;   %单位是anstrong
    %左上角和右下角的坐标（x，y）值。单位为anstrong。笔记中的红色框

    %放在开始的时候赋初值outside_ext=4;  %再向外延拓4埃，这样保证投影势场做傅里叶变换时候，周期延拓的效果不会影响probe扫描区域。
                %笔记中的绿色框，单位也是埃
    tl_green=tl_red-handles.sampling*round(handles.outside_ext./handles.sampling)-handles.probesx/2*handles.sampling;  %第二个减法的操作，是为了保证在扫描时，所选的区域左上角坐标是整数个像素
    rd_green=rd_red+handles.sampling*round(handles.outside_ext./handles.sampling)+handles.probesx/2*handles.sampling;  %第二个加法的操作，是为了保证在扫描时，所选的区域左上角坐标是整数个像素

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %绿色区域的图像尺寸，单位是像素
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
      %给出probe位于绿色框内的起始坐标点位置;搜索的尺寸大小见width_red和hight_red
    probe_in_green=[round(handles.outside_ext./handles.sampling)+1, round(handles.outside_ext./handles.sampling)+1];  

    %需要存一下绿色框的像素数，抽样率；probe的尺寸；红色框的宽和高。~~~~~~~~~~~~~
    load allresul
    atompos= allresul.x_new(:,2:4);  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；

    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %行和列
    handles.probe_ingreenNrow=probe_in_green(1);  
    handles.probe_ingreenNcol=probe_in_green(2);   %在绿色框中选取的坐标位置
    handles.width_red=width_red; 
    handles.hight_red=hight_red;  %扫描区的像素个数
    handles.probestep=str2num(get(handles.edit77,'string')); %扫描时的步长(单位是倍数）
    handles.tl_green=tl_green;   %把绿色框的左上角坐标（埃单位）保存下来，需要之后换算所有的原子位置，并做转换
    handles.rd_green=tl_green+greensize*handles.sampling;%绿色框的右下角坐标，需要考虑哪些原子在这个范围外的，就不带入stem成像的后续计算。

    maxz=max(atompos(:,3));
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         %画一个probe的形状，这个probe的中心位于红色区域的左上角坐标
    tl_topleftprobe=tl_green+(probe_in_green-1)*handles.sampling; rd_topleftprobe=tl_green+(probe_in_green-1)*handles.sampling+handles.probesx*handles.sampling;
    hold on; %plot(tl_topleftprobe(1), tl_topleftprobe(2), 'ko');plot(rd_topleftprobe(1), rd_topleftprobe(2), 'ko')
         line([tl_topleftprobe(1),rd_topleftprobe(1)],[tl_topleftprobe(2),tl_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([tl_topleftprobe(1),rd_topleftprobe(1)],[rd_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([tl_topleftprobe(1),tl_topleftprobe(1)],[tl_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([rd_topleftprobe(1),rd_topleftprobe(1)],[tl_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
    view([0,0,-1])
    axis equal
end

if get(handles.radiobutton10, 'value')  %如果是CBED图
    handles.sampling=str2num(get(handles.edit75,'string'));  %埃/pixel   %这个是传播矩阵和投影势函数矩阵的抽样率，不是probe的扫描抽样率
    handles.CBEDprobesx=str2num(get(handles.edit4,'string'));  %读取probe尺寸
    if rem(handles.CBEDprobesx,2)==1 handles.CBEDprobesx=handles.CBEDprobesx+1;end


    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    %红色的是作图区域，这里是probe的尺寸
    rd_red=tl_red+handles.CBEDprobesx*handles.sampling; 

    tl_green=tl_red-handles.sampling*round(handles.outside_ext./handles.sampling);  %埃单位  %绿色的左上角和右下角坐标。图像需要扩展一定尺寸，否则会有wrap
    rd_green=rd_red+handles.sampling*round(handles.outside_ext./handles.sampling);  

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %绿色区域的图像尺寸，单位是像素
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
    %给出imaging位于绿色框内的起始坐标点位置;之后从大图中裁剪最终成像结果
    imaging_in_green=[round(handles.outside_ext./handles.sampling)+1, round(handles.outside_ext./handles.sampling)+1];

    load allresul
    x=allresul.x_new;
    atompos=x(:,2:4);  % 单位变埃~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    handles.atompos= atompos;  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；

    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %行和列
    handles.CBED_ingreenNrow=imaging_in_green(1); %注意已经+1了，所以从这点开始直接取位置 
    handles.CBED_ingreenNcol=imaging_in_green(2);   %在绿色框中选取的坐标位置
    handles.tl_green=tl_green;   %把绿色框的左上角坐标（埃单位）保存下来，需要之后换算所有的原子位置，并做转换
    handles.rd_green=tl_green+greensize*handles.sampling;%绿色框的右下角坐标，需要考虑哪些原子在这个范围外的，就不带入stem成像的后续计算。

    maxz=max(atompos(:,3));
    load allresul
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
    view([0,0,-1])
    axis equal
end

if get(handles.radiobutton9, 'value')  %如果是HRTEM的计算，显示相应的区域
    handles.sampling=str2num(get(handles.edit75,'string'));  %埃/pixel   %这个是传播矩阵和投影势函数矩阵的抽样率，不是probe的扫描抽样率
    handles.HRTEMedgesx=str2num(get(handles.edit4,'string'));  %读取外围有多大区域需要裁掉

    %设定成像区域topleft-rightdown的坐标
    %红色框的宽度和高度分别为多少像素
    width_red=str2num(get(handles.edit20,'string')); hight_red=str2num(get(handles.edit19,'string'));
    
    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    rd_red=tl_red+[width_red, hight_red]*handles.sampling;   %单位是anstrong  %这里与STEM不同
    %左上角和右下角的坐标（x，y）值。单位为anstrong。

    %绿色框，单位也是埃；是实际成像时的左上角和右下角坐标，之后要扣掉剪裁区；并且这个是边界绝对的尺寸，不需要除以2；也就是说，扣掉256的话，成像时候选择的大小至少要512
    tl_green=tl_red-handles.HRTEMedgesx*handles.sampling;  %埃单位
    rd_green=rd_red+handles.HRTEMedgesx*handles.sampling;  

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %绿色区域的图像尺寸，单位是像素
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
    %给出imaging位于绿色框内的起始坐标点位置;之后从大图中裁剪最终成像结果
    imaging_in_green=[handles.HRTEMedgesx+1, handles.HRTEMedgesx+1];  

    %需要存一下绿色框的像素数，抽样率；probe的尺寸；红色框的宽和高。~~~~~~~~~~~~~


    x=handles.x;
    atompos=x(:,2:4);  % 单位变埃~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %需要换算出，所有的原子位置都放在绿色的框里面的起始位置是什么

    handles.atompos= atompos;  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；

%     handles.eachthick=eachthick;  %记录厚度
    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %行和列
    handles.HRTEM_ingreenNrow=imaging_in_green(1); %注意已经+1了，所以从这点开始直接取位置 
    handles.HRTEM_ingreenNcol=imaging_in_green(2);   %在绿色框中选取的坐标位置
    handles.tl_green=tl_green;   %把绿色框的左上角坐标（埃单位）保存下来，需要之后换算所有的原子位置，并做转换
    handles.rd_green=tl_green+greensize*handles.sampling;%绿色框的右下角坐标，需要考虑哪些原子在这个范围外的，就不带入stem成像的后续计算。
    handles.width_red=width_red;
    handles.hight_red=hight_red;
    
    maxz=max(atompos(:,3));
    load allresul
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    %axes(handles.axes1);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
    view([0,0,-1])
    axis equal
end
guidata(hObject, handles);





function edit78_Callback(hObject, eventdata, handles)
function edit78_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton4_CreateFcn(hObject, eventdata, handles)



function edit79_Callback(hObject, eventdata, handles)
%求波长，求mrad，以确定最小的sampling rate应该是多少
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %计算波长，单位nm

aper_dect = str2num(get(handles.edit79,'string'));  %detector 的光阑尺寸
kmax=max(aper_dect(:))*10^(-3);

kk=kmax./(lambda*10);  %最高的倒空间的频率是多少，单位是1/埃

%由于正空间的最小间隔单位为 sampling，所以最高的倒空间频率为1/sampling；
%另外，由于s的最大有效值为6 埃分之一值；划到g的值的话，为2*s=12 埃分之一；
%但是计算波传播的时候，会加上1/2的光阑（或2/3的光阑）以保证不会发生 wrap 的情况
%所以，推荐 1/sampling=kk 要小于 6
%由于6之后的数据是不准的，至少用现在方法，但是为了保证高频信息还有，就需要samping rate足够小。

%计算最高的空间rad值为多少。因为到6是有效的，所以乘以波长后，就是最高的mrad值
disp(strcat('Max aperture of detector is ', num2str(6*lambda*10 *1000),' mrad'));
disp('Because 6 1/A is maximum reciprocal-lattice-vector to calculate the projected potential')
%计算中，为了防止wrap，计算过程倒易空间都加了半宽的光阑。
%为了保证到了较高空间还能有倒空间的分量，0.5/sampling 是最高的倒空间频率， 0.5*(0.5*1/samping)*lambda=alpha mrad
% %所以 sampling=lambda*0.25/alpha ,注意化到埃和rad的单位;其中一个0.5是-0.5/sampling 到0.5/sampling的最高频率
% disp(strcat('Imaging sampling rate should be smaller than ', num2str(0.5*0.5*lambda*10/kmax),' A/pixel'));
%所以 sampling=lambda*0.25/alpha ,注意化到埃和rad的单位;其中一个0.5是-0.5/sampling
%到0.5/sampling的最高频率;现在取的是2/3的光阑
disp(strcat('Imaging sampling rate should be smaller than ', num2str(0.5*2/3*lambda*10/kmax),' A/pixel'));
pp=1;



% --- Executes during object creation, after setting all properties.
function edit79_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
if get(handles.radiobutton10,'value')==1
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
    set(handles.radiobutton13, 'value', 0);

end
if get(handles.radiobutton10,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM相关关闭
    set(handles.text106, 'visible', 'off');  %与STEM相关
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'off');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture的mrad和1/nm单位
    set(handles.text134, 'visible', 'off');
end
guidata(hObject, handles);


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
if get(handles.radiobutton9,'value')==1
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
    set(handles.radiobutton13, 'value', 0);

end

if get(handles.radiobutton9,'value')==1
    set(handles.text115, 'visible', 'on');  %HRTEM相关
    
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
    
   
    set(handles.text106, 'visible', 'off');  %与STEM相关
    set(handles.text108, 'visible','off');
    set(handles.edit77, 'visible','off');
    
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'off');
    set(handles.text134, 'visible', 'on');   % convergence & aperture的mrad和1/nm单位
    set(handles.text24, 'visible', 'off');

else
    
    set(handles.text115, 'visible', 'off');  %HRTEM相关

end

guidata(hObject, handles);



% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
%20200429修改，整个代码全改了
probesx=str2num(get(handles.edit4,'String'));
   [para_part1, U, V]=CommonPara_TEM(handles, probesx, probesx);  %画CTF和PhasePlate用固定的大小即可
   para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码
    %参数分为两个部分，lambda需要都有
   residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %计算phase plate visualing the residual aberrations
    %-------------------------------------
   mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,2);
   
   myap=ones(probesx);  %构造光阑函数，kmax是外面定的
   myap=myaperture(myap,1/(para_part1.sampling.*probesx),1/(para_part1.sampling.*probesx),0,para_part1.gmax*0.01/para_part1.lambda,0,0,1); %如果10mrad，读入就是1；再乘以0.01，换到到rad，再除以nm单位的波长
   num=length(mytcc(1,1,:));  %考察有多少个probe需要扫过图像
   sum_probe=zeros(probesx);
   
   handles.gfsf=para_part1.gfsf;
   for i=1:num
        guiyihua=sum(sum(abs(ifft2(ifftshift(myap.*mytcc(:,:,i)))).^2));  %20201226 add one aper to normalize 
        guiyihua=sqrt(guiyihua);  %20201231需要求开根号，才能让probe满足3.59和3.68式
        sum_probe = sum_probe + handles.gfsf(i)*myap.*mytcc(:,:,i)./guiyihua;
   end
   

  [sx,sy]=size(myap);
  temp = abs(fftshift(ifft2(ifftshift(sum_probe)))).^2;

figure;imshow(abs(fftshift(ifft2(ifftshift(sum_probe)))).^2,'XData',(-probesx/2:probesx/2-1).*para_part1.sampling,'YData',(-probesx/2:probesx/2-1).*para_part1.sampling,'DisplayRange',[]);axis on;xlabel('nm');
title('Intensity of Probe Wave')
figure;plot((-probesx/2:probesx/2-1).*para_part1.sampling, temp(sx/2+1,:) ,'r','linewidth', 2);
grid on
title('Intensity profile of Probe Wave')
figure;imshow(-angle(sum_probe.*myap),'XData',U(1, 1:end),'YData',V(1:end, 1),'DisplayRange',[]);axis on;xlabel('1/nm');
title('Phase of Probe Wave')


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
CTFflag=0;
CTF_Callback(hObject, eventdata, handles, CTFflag)



function CTF_Callback(hObject, eventdata, handles, CTFflag)
sx=256;sy=256;
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %画CTF和PhasePlate用固定的大小即可
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码
%参数分为两个部分，lambda需要都有
residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %计算phase plate visualing the residual aberrations
%-------------------------------------
%mytcc=WholeTCC2D(tccpara,defocus, U,V);
mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,1);
myresul=zeros(sx,sy);
for i=1:length(para_part1.gfsf);
     myresul=myresul+para_part1.gfsf(i).*mytcc(:,:,i);
end

figure;hold off;
if CTFflag==0;
    myresul=real(myresul);
elseif CTFflag==1
    myresul=imag(myresul);
elseif CTFflag==2
    myresul=abs(myresul);
end
myap=ones(size(U));
myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);
disimage=myresul.*myap;

imshow(disimage,'XData',U(1, 1:end).*para_part2.lambda*1000,'YData',V(1:end, 1).*para_part2.lambda*1000,'DisplayRange',[]);axis on;xlabel('mrad');

colorbar;
if CTFflag==0;
    title('Real CTF');
elseif CTFflag==1
    title('Imaginary CTF');
elseif CTFflag==2
    title('Damping CTF');
end
return;


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton11,'value')==1 
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton13, 'value', 0);
end
if get(handles.radiobutton11,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM相关关闭
    
    set(handles.text21, 'visible', 'off');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture的mrad和1/nm单位
    set(handles.text134, 'visible', 'off');
    
    set(handles.text106, 'visible', 'on');  %STEM相关
    set(handles.text108, 'visible','on');
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit77, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
else
    set(handles.text106, 'visible', 'off');
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'off');
    set(handles.text24, 'visible', 'off');   % convergence & aperture的mrad和1/nm单位
    set(handles.text134, 'visible', 'on');
    
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton11


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
sx=256; sy=256;
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %画CTF和PhasePlate用固定的大小即可
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码
%参数分为两个部分，lambda需要都有

residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %计算phase plate visualing the residual aberrations
residual_phase=angle(exp(sqrt(-1)*2*pi*residual_aberr/para_part1.lambda));


figure; hold off; imshow(residual_phase,'XData',U(1, 1:end).*para_part2.lambda*1000,'YData',V(1:end, 1).*para_part2.lambda*1000,'DisplayRange',[]);axis on;xlabel('mrad');title('ctg(kai)');
axis equal
    
prompt = {'Radiu 1'; 'Radiu 2'};
dlg_title = '2D Phase Plate';
num_lines = 1;
def = {'16';'28'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer) 
    return; 
end
rad1=str2num(answer{1}); rad2=str2num(answer{2}); %画两个环，把残留像散的图像圈出来

title(strcat('2D Phase Plate ',num2str(rad1), '&',num2str(rad2)));
colorbar
hold on; grid on  %这里把mrad的圆画上
x=linspace(-rad1,rad1,sx*5);
y=sqrt(rad1.*rad1-x.*x);
plot(x,y,'r-','LineWidth',3)
plot(x,-y,'r-','LineWidth',3)
%rad2=0.5*rad1;  %以1/2倍aper的数值来画图，只是换算到mrad单位
x=linspace(-rad2,rad2,sx*20);
y=sqrt(rad2.*rad2-x.*x);
plot(x,y,'b-','LineWidth',3)
plot(x,-y,'b-','LineWidth',3)
%title(strcat(num2str(rad1),' and  ', num2str(rad2), 'mrad phase plate due to the residual aberrations'))

guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
CTFflag=1;
CTF_Callback(hObject, eventdata, handles, CTFflag);
pp=1;


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
CTFflag=2;
CTF_Callback(hObject, eventdata, handles, CTFflag);



function [tccpara, U, V]=CommonPara_TEM(handles, sx, sy);
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
tccpara.sampling=str2num(get(handles.edit75,'String'))*0.1; %  抽样率 %UNIT NM/PIXEL  %这里与EW_RECONSTRUCT不一样的单位
tccpara.mtotal=(get(handles.popupmenu4,'Value')-1)*2+1; %  总共有多少个离焦偏离，比如这个值为5，就必须有1 2 3 4 5个mytcc矩阵，其中第三张图无离焦偏离

tccpara.alafa=str2num(get(handles.edit_Con,'String'))*(10^(-3)); % 半角宽度
Vol=str2num(get(handles.edit_Vol,'String')); %电压
tccpara.lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %计算波长，单位nm
strcat( 'wavelength of incident electron (nm): ', num2str(tccpara.lambda) )

tcccpara.yita=str2num(get(handles.edit_Spr,'String')); 
tccpara.gmax=str2num(get(handles.edit_Ape,'String'))*0.1;   %换算到A的单位
 
resul=Gaussian_focal(tccpara.mtotal,tcccpara.yita,tccpara.lambda,tccpara.gmax); %为离焦图像所占总图像的百分之几
tccpara.gfsf=resul.gfsf;
tccpara.delta_yita=resul.delta_yita;

%束倾斜多少角度
tccpara.tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);tccpara.phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % 弧度单位 例如 10mrad  
tccpara.tiltx=tccpara.tilt/tccpara.lambda*cos(tccpara.phitilt/180*pi); %换算到倒易格式的单位 1/nm
tccpara.tilty=tccpara.tilt/tccpara.lambda*sin(tccpara.phitilt/180*pi);


v=-round((sx+1)/2)+1:sx-(round((sx+1)/2));
u=-round((sy+1)/2)+1:sy-(round((sy+1)/2));
u=u.*1/(tccpara.sampling.*sy);
v=v.*1/(tccpara.sampling.*sx);
[U,V]=meshgrid(u,v);
U=U+tccpara.tiltx;
V=V+tccpara.tilty;


function tccpara=readtccfromhandles_newTEM(hObject, handles, lambda)  %本次是根据TEM的参数来读取的
tccpara.lambda=lambda;
%if get(handles.checkbox_polar,'Value')==1  %必定是极坐标系
   %读取说明：%参考Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;但是圈出来的系数在UM  1998 72 PP109-119以及UM 1996 64 249-264中，B2 S3 D4 B4
%S5 R5(D5)，并没有系数，因此，读取数据后，应该要分别乘以 3  4  5  5  6  6

   tccpara.focus=str2num(get(handles.edit_Focus,'String'));
   
   %A1的数值不变；A1的角度是电镜显示的数值的1/2 
   tccpara.A1=str2num(get(handles.edit_A1_1,'String')); tccpara.phiA1=-str2num(get(handles.edit_A1_2,'String'));
   
   %A2的数值不变；A2的角度是电镜显示的数值的1/3
   tccpara.A2=str2num(get(handles.edit_A2_1,'String')); tccpara.phiA2=-str2num(get(handles.edit_A2_2,'String'));
  
   %B2的数值不变；B2的角度与电镜显示的取负数
   tccpara.B2=str2num(get(handles.edit_B2_1,'String')); %这个参数引入时候和rew相差了3倍。例如rew里面为2000，这里必须是6000；因为程序中的系数多乘了1/3
   tccpara.phiB2=-str2num(get(handles.edit_B2_2,'String'));  
  
   %Cs的单位从um到nm
   tccpara.Cs=str2num(get(handles.edit_Cs,'String'))*(10^3); %UNIT from um to nm
   
   %A3的单位从um到nm；A3的角度与电镜显示的除以4
   tccpara.A3=str2num(get(handles.edit_A3_1,'String'))*(10^3); %UNIT from um to nm
   tccpara.phiA3=-str2num(get(handles.edit_A3_2,'String'));
   
   %S3的单位从um到nm；S3的角度与电镜显示的取负数除以2
   tccpara.S3=str2num(get(handles.edit_S3_1,'String'))*(10^3); 
   tccpara.phiS3=-str2num(get(handles.edit_S3_2,'String'));
   
   %A4的单位从um到nm；A4的角度与电镜显示的除以5
   tccpara.A4=str2num(get(handles.edit_A4_1,'String'))*(10^3);
   tccpara.phiA4=-str2num(get(handles.edit_A4_2,'String'));
   
   %D4的单位从um到nm；D4的角度与电镜显示的取负数除3
   tccpara.D4=str2num(get(handles.edit_D4_1,'String'))*(10^3); 
   tccpara.phiD4=-str2num(get(handles.edit_D4_2,'String'));
   
   %B4的单位从um到nm；B4的角度与电镜显示的一样   %需要再研究的角度关系！！！
   tccpara.B4=str2num(get(handles.edit_B4_1,'String'))*(10^3); 
   tccpara.phiB4=-str2num(get(handles.edit_B4_2,'String'));
   
    %A5的单位从mm到nm；A4的角度与电镜显示的除以6
   tccpara.A5=str2num(get(handles.edit_A5_1,'String'))*(10^6); 
   tccpara.phiA5=-str2num(get(handles.edit_A5_2,'String'));
   
   %A5的单位从mm到nm；
   tccpara.C5=str2num(get(handles.edit_C5,'String'))*(10^6);
   
   %S5的单位从mm到nm； %这个数值在hrtem里面还没有出现，因此角度关系暂时没有管
   tccpara.S5=str2num(get(handles.edit_S5_1,'String'))*(10^6); tccpara.phiS5=-str2num(get(handles.edit_S5_2,'String'));
  
   %D5的单位从mm到nm；%这个数值在hrtem里面还没有出现，因此角度关系暂时没有管
   tccpara.D5=str2num(get(handles.edit_D5_1,'String'))*(10^6); tccpara.phiD5=-str2num(get(handles.edit_D5_2,'String'));

   
   
   tccpara.A1x=tccpara.A1*cos(tccpara.phiA1/180*pi);   tccpara.A1y=tccpara.A1*sin(tccpara.phiA1/180*pi);
   %C1=focus;  focus是在输入的defocs附近有偏离的
   tccpara.A2x=tccpara.A2*cos(tccpara.phiA2/180*pi);   tccpara.A2y=tccpara.A2*sin(tccpara.phiA2/180*pi);
   tccpara.B2x=tccpara.B2*cos(tccpara.phiB2/180*pi);     tccpara.B2y=tccpara.B2*sin(tccpara.phiB2/180*pi);
   tccpara.A3x=tccpara.A3*cos(tccpara.phiA3/180*pi);   tccpara.A3y=tccpara.A3*sin(tccpara.phiA3/180*pi);
   tccpara.S3x=tccpara.S3*cos(tccpara.phiS3/180*pi);   tccpara.S3y=tccpara.S3*sin(tccpara.phiS3/180*pi);
   tccpara.C3=tccpara.Cs;
   tccpara.A4x=tccpara.A4*cos(tccpara.phiA4/180*pi);   tccpara.A4y=tccpara.A4*sin(tccpara.phiA4/180*pi);
   tccpara.D4x=tccpara.D4*cos(tccpara.phiD4/180*pi);   tccpara.D4y=tccpara.D4*sin(tccpara.phiD4/180*pi);
   tccpara.B4x=tccpara.B4*cos(tccpara.phiB4/180*pi);     tccpara.B4y=tccpara.B4*sin(tccpara.phiB4/180*pi);
   tccpara.A5x=tccpara.A5*cos(tccpara.phiA5/180*pi);   tccpara.A5y=tccpara.A5*sin(tccpara.phiA5/180*pi);
   tccpara.S5x=tccpara.S5*cos(tccpara.phiS5/180*pi);   tccpara.S5y=tccpara.S5*sin(tccpara.phiS5/180*pi);
   tccpara.C5=tccpara.C5;
   tccpara.D5x=tccpara.D5*cos(tccpara.phiD5/180*pi);   tccpara.D5y=tccpara.D5*sin(tccpara.phiD5/180*pi);

pp=1;

function W=WholeTCC_residual_newTEM(para,u,v);   %按照公式，把TCC的所有系数都考虑进去
%配合 formular公式，该公式把所有的高阶像散形式都表现出来，且求得了对u和v方向的梯度

oumigau=u.*para.lambda;  %
oumigav=v.*para.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para.A1x;   A1y=para.A1y;
%C1=focus;  focus是在输入的defocs附近有偏离的
A2x=para.A2x;   A2y=para.A2y;
B2x=para.B2x;   B2y=para.B2y;
A3x=para.A3x;   A3y=para.A3y;
S3x=para.S3x;   S3y=para.S3y;
C3=para.Cs;
A4x=para.A4x;   A4y=para.A4y;
D4x=para.D4x;   D4y=para.D4y;
B4x=para.B4x;   B4y=para.B4y;
A5x=para.A5x;   A5y=para.A5y;
S5x=para.S5x;   S5y=para.S5y;
C5=para.C5;
D5x=para.D5x;   D5y=para.D5y;

%参考Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;公式2.4式子，oumiga=u*lambda+i*v*lambda;  结合Ultramicroscopy 64
%1996, pp.249-2264的equation 11
W= real( ((A2x/3 + (A2y*i)/3).*(oumigau - i*oumigav).^3) ) ...
    +real (((B2x + (B2y*i)).*(oumigau + i*oumigav).^2.*(oumigau - i*oumigav)) ) ...
    +real ( ((A3x/4 + (A3y*i)/4).*(oumigau - i*oumigav).^4) ) ...
    +real ( ((oumigau - i*oumigav).*(oumigau + i*oumigav).^3.*(S3x + (S3y*i)))  ) ...
    +  real(C3*(oumigau + i*oumigav).^2.*(oumigau - i*oumigav).^2)/4 ...
    + real( ((A4x/5 + (A4y*i)/5).*(oumigau - i*oumigav).^5) ) ...
    + real( ((D4x + (D4y*i)).*(oumigau + i*oumigav).^4.*(oumigau - i*oumigav)) ) ...
    + real( ((B4x + (B4y*i)).*((oumigau + i*oumigav).^3).*((oumigau - i*oumigav).^2)))  ...
    + real( ((A5x/6 + (A5y*i)/6).*(oumigau - i*oumigav).^6) ) ...
    + real( (oumigau + i*oumigav).^4.*(oumigau - i*oumigav).^2.*(S5x + (S5y*i)) ) ...
    + real(C5*(oumigau + i*oumigav).^3.*(oumigau - i*oumigav).^3)/6 ...
    + real( ((D5x + (D5y*i)).*(oumigau + i*oumigav).^5.*(oumigau - i*oumigav)) )  ;

    
    pp=1;


function [tlittle_all, damping]=WholeTCC2D_newTEM(para1, para2, u,v, flag);   %按照公式，把TCC的所有系数都考虑进去
%配合 formular公式，该公式把所有的高阶像散形式都表现出来，且求得了对u和v方向的梯度

oumigau=u.*para1.lambda;  %
oumigav=v.*para1.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para2.A1x;   A1y=para2.A1y;
C1=para2.focus;
A2x=para2.A2x;   A2y=para2.A2y;
B2x=para2.B2x;   B2y=para2.B2y;
A3x=para2.A3x;   A3y=para2.A3y;
S3x=para2.S3x;   S3y=para2.S3y;
C3=para2.Cs;
A4x=para2.A4x;   A4y=para2.A4y;
D4x=para2.D4x;   D4y=para2.D4y;
B4x=para2.B4x;   B4y=para2.B4y;
A5x=para2.A5x;   A5y=para2.A5y;
S5x=para2.S5x;   S5y=para2.S5y;
C5=para2.C5;
D5x=para2.D5x;   D5y=para2.D5y;

%参考Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;公式2.4式子，oumiga=u*lambda+i*v*lambda;
kai=WholeTCC_residual_newTEM(para2,u,v);
W=real(((A1x + (A1y*i)).*(oumigau - i*oumigav).^2))/2 ...   %'2-fold astigmatism A1'
    + real(C1.*(oumigau + i*oumigav).*(oumigau - i*oumigav))/2 + kai;

flag=2;
if flag==2  %与原程序相比，B等都有变化
%new
tempgradient_kai_u=para1.lambda.*(A1x*u+A1y*v) ...%    +C1.*para.lambda.*u ...
    +A2x*para1.lambda.^2/3*(3*u.*u-3*v.*v)+A2y*para1.lambda.^2/3*6*u.*v ...
    +B2x.*(3*u.^2+v.^2).*para1.lambda.^2 - B2y.*para1.lambda.^2.*2.*u.*v ...
    +1/4*para1.lambda.^3 * (A3x*(4*u.^3-12*u.*v.*v)+A3y*(12*u.*u.*v-4*v.^3)) ...
    +para1.lambda.^3* (2*u.*(S3x*(u.*u-v.*v)-S3y*2*u.*v) +(u.*u+v.*v).*(2*S3x*u-2*S3y*v)) ...
    +C3*para1.lambda^3*(u.^2+v.^2).*u ...
    +1/5*para1.lambda^4 * (A4x*(5*u.^4+5*v.^4-30.*u.^2.*v.^2)+A4y*(20*u.^3.*v-20*u.*v.^3)) ...
    +para1.lambda^4 * ( (2*u) .* (D4x*(u.^3-3*u.*v.^2)-D4y*(3*u.^2.*v-v.^3)) + (u.^2+v.^2).*(D4x*(3.*u.^2-3*v.^2)-D4y*6*u.*v) ) ...
    +para1.lambda^4 * (  2*(u.^2+v.^2).*2.*u .*(B4x*u-B4y*v) + (u.^2+v.^2).^2*B4x) ...
    +1/6*para1.lambda^5 * (A5x*(4.*u.^3-14*v).*(u.^2-v.^2)+A5x*(u.^4+v.^4-14*u.*v)*2.*u   - A5y*2*v.*(3*u.^4+3*v.^4-10*u.^2.*v.^2)-A5y*2.*v.*u.*(12*u.^3-20*u.*v.^2) ) ...
    +para1.lambda^5 * (2*(u.^2+v.^2).*2.*u.*(S5x*(u.^2-v.^2)-2*u.*v*S5y) + (u.^2+v.^2).^2.*(2*S5x*u-2*v*S5y) ) ...
    +1/6*para1.lambda^5 * C5 * 3*(u.^2+v.^2).^2.*2.*u ...
    +para1.lambda^5 * (2*u.*(D5x*(u.^4+v.^4-6*u.^2.*v.^2)-D5y*4*u.*v.*(u.^2-v.^2)) + (u.^2+v.^2).*(D5x.*(4.*u.^3-12*u.*v.^2)-D5y*(12*u.^2.*v-4*v.^3)) ) ;
    
    
tempgradient_kai_v=para1.lambda.*(A1y*u-A1x*v) ...%    +C1.*para.lambda.*v ...
    -A3x.*para1.lambda.^2/3*6*u.*v + A3y*para1.lambda.^2/3*(-3*v.^2+3*u.^2) ...
    + B2x.*para1.lambda.*para1.lambda.*2*u.*v - B2y.*para1.lambda.*para1.lambda.*(u.^2+3*v.^2) ...
    +1/4*para1.lambda.^3 * (A3x*(4*v.^3-12*u.*u.*v)+A3y*(4*u.^3-12*u.*v.*v)) ...
    +para1.lambda.^3* (2*v.*(S3x*(u.*u-v.*v)-S3y*2*u.*v) +(u.*u+v.*v).*(-2*S3x*v-2*S3y*u) ) ...
    +C3*para1.lambda.^3.*(u.^2+v.^2).*v ...
    +1/5*para1.lambda^4 * (A4x*(20*u.*v.^3-20.*u.^3.*v)+A4y*(5*v.^4+5*u.^4-30*u.^2.*v.^2)) ...
    +para1.lambda^4 * ( (2*v) .* (D4x*(u.^3-3*u.*v.^2)-D4y*(3*u.^2.*v-v.^3)) + (u.^2+v.^2).*(-D4x*6*u.*v-D4y*(3*u.^2-3*v.^2))) ...
    +para1.lambda^4 * (  2*(u.^2+v.^2).*2.*v .*(B4x*u-B4y*v) - (u.^2+v.^2).^2*B4y) ...
    +1/6*para1.lambda^5 * (A5x*(4.*v.^3-14*u).*(u.^2-v.^2)+A5x*(u.^4+v.^4-14*u.*v).*(-2.*v)   - A5y*2*u.*(3*u.^4+3*v.^4-10*u.^2.*v.^2)-A5y*2.*v.*u.*(12*v.^3-20*u.^2.*v) ) ...
    +para1.lambda^5 * (2*(u.^2+v.^2).*2.*v.*(S5x*(u.^2-v.^2)-2*u.*v*S5y) + (u.^2+v.^2).^2.*(-2*S5x*v-2*u*S5y) ) ...
    +1/6*para1.lambda^5 * C5 * 3*(u.^2+v.^2).^2.*2.*v ...
    +para1.lambda^5 * (2*v.*(D5x*(u.^4+v.^4-6*u.^2.*v.^2)-D5y*4*u.*v.*(u.^2-v.^2)) + (u.^2+v.^2).*(D5x.*(4.*v.^3-12*u.^2.*v)-D5y*(4*u.^3-12*u.*v.^2)) ) ;

end
if flag==1
    tempgradient_kai_u=0;
    tempgradient_kai_v=0;
end

for m=1:para1.mtotal
%    C1=focus+(m-(para.mtotal+1)/2)*para.delta_yita;
    
    gradient_kai_u=tempgradient_kai_u+C1.*para1.lambda.*u;
    gradient_kai_v=tempgradient_kai_v+C1.*para1.lambda.*v;

%gradient_W=gradient_W_u + i* gradient_W_v;
%参考um 1996 vol 64:109-135，公式5
E_s_coh=exp(-para1.alafa.*para1.alafa.*pi*pi/para1.lambda/para1.lambda*(gradient_kai_u.^2+gradient_kai_v.^2));
kai=2*pi*W/para1.lambda;
%参考Phil. Trans. R. Soc. A 2009, vol 367:3755-3771;公式2.2下方
%E_s_coh=exp(-para.alpha.*para.alpha/4/para.lambda/para.lambda.*gradient_W.*conj(gradient_W).*(2*pi*2*pi/para.lambda/para.lambda));  %空间相干性 这里para.alpha为半角;gradient_kai实际上是书上的gradient_W，需要再乘以w*pi的
%本表式根据Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;公式2.1写出；注意gradient_W_u_v实际上是对W求导的。
%与参考um 1996 vol 64:109-135相比，该文献中的 kai=W/lambda，W为phil的文献变量，
%所以参考um 1996 vol 64:109-135，公式5,与本表示结果一致
E_f_coh=exp(-i*pi*(m-(para1.mtotal+1)/2)*para1.delta_yita*(u.^2+v.^2).*para1.lambda);  %色差   这里的u和v，是倒格矢量
%参考um 1996 vol 64:109-135，公式35,与本表式来源
t=exp(-i*kai);
% if m==1
% figure; imshow(real(t),[]);title('real t')
% end
% figure; imshow(real(E_f_coh),[]);title('E f coh real')
% figure; imshow(E_s_coh,[]);title('E s coh')
%参考Phil. Trans. R. Soc. A 2009, vol 367:3755-3771;公式2.2下方
tlittle_all(:,:,m)=E_s_coh.*E_f_coh.*t;
damping(:,:,m)=E_s_coh.*E_f_coh;
%滤波
[sx,sy]=size(u);
if min(sx,sy)>2  %如果不是一维数据，就做带通滤波
%   tlittle_all(:,:,m)=myaperture(tlittle_all(:,:,m),1/(para.sampling.*sx),1/(para.sampling.*sy),0,para.gmax,0,0,1);
else 
    
   tlittle_all(1,round(para1.gmax./para1.sampling):end,m)=0;  %一维函数的带通，高频处都为0，此时para.sampling是表示倒空间的间隔频率是0.01nm_1的值，固定下来了
end
end

   
function W=WholeTCC2D_newTEMphase(para1, para2, u,v);   %按照公式，把TCC的所有系数都考虑进去
%配合 formular公式，该公式把所有的高阶像散形式都表现出来，且求得了对u和v方向的梯度

oumigau=u.*para1.lambda;  %
oumigav=v.*para1.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para2.A1x;   A1y=para2.A1y;
C1=para2.focus;
A2x=para2.A2x;   A2y=para2.A2y;
B2x=para2.B2x;   B2y=para2.B2y;
A3x=para2.A3x;   A3y=para2.A3y;
S3x=para2.S3x;   S3y=para2.S3y;
C3=para2.Cs;
A4x=para2.A4x;   A4y=para2.A4y;
D4x=para2.D4x;   D4y=para2.D4y;
B4x=para2.B4x;   B4y=para2.B4y;
A5x=para2.A5x;   A5y=para2.A5y;
S5x=para2.S5x;   S5y=para2.S5y;
C5=para2.C5;
D5x=para2.D5x;   D5y=para2.D5y;

%参考Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;公式2.4式子，oumiga=u*lambda+i*v*lambda;
kai=WholeTCC_residual_newTEM(para2,u,v);
W=real(((A1x + (A1y*i)).*(oumigau - i*oumigav).^2))/2 ...   %'2-fold astigmatism A1'
    + real(C1.*(oumigau + i*oumigav).*(oumigau - i*oumigav))/2 + kai;

% 
% function edit83_Callback(hObject, eventdata, handles)
% function edit83_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% 
% 
% function edit84_Callback(hObject, eventdata, handles)
% function edit84_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end






function edit92_Callback(hObject, eventdata, handles)
function edit92_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)  
% in dialog screen, atoms in initial crystal of handles.x will be drawed
% and the view point is along the at&el direction 
at = str2num(get(handles.edit91,'string'));
el = str2num(get(handles.edit92,'string'));
step1forshowcrystal(at, el, hObject, eventdata, handles);

function step1forshowcrystal(at, el, hObject, eventdata, handles) 
handles.x_new(:,2:4) = handles.x(:,2:4);       % read coordinates from initial crystal     
Drawsupercell(hObject, eventdata, handles)
axes(handles.axes1)
hold on;
view(at, el);  %view along the at&el direction


%the coordinates of atoms in projection will be edited, which will be
%sliced by multi-slice method
atompos=handles.x(:,2:4);  % read atoms' coordinators
% rotate atoms via RR matrix
at=-at;
el=90+el;

RR=[cosd(at) -sind(at) 0; cosd(el)*sind(at)  cosd(el)*cosd(at)  -sind(el);
sind(at)*sind(el) sind(el)*cosd(at) cosd(el)];
newatom=(RR*atompos.').';

%these coordinates will be sliced in the following codes 
atompos(:,3) = newatom(:,3)-min(newatom(:,3));  % the beginning height is 0 
atompos(:,1) = newatom(:,1)-min(newatom(:,1))+str2num(get(handles.edit100,'string'));
atompos(:,2) = newatom(:,2)-min(newatom(:,2))+str2num(get(handles.edit101,'string'));
handles.x_new(:,2:4)=atompos;
handles.x_mysavenewz = atompos(:,3);  %remember z to slicing crystal better 把z轴的值记录住，这样在划分层后，可以随时显示正确的原子高度；但如果分层不如意，下次再分层就调用这个高度

Drawsupercell_figure(handles.x_new, handles.dis_atomsize, handles.x_mysavenewz, 0, 0)
title('Rotate crystal and view from top(Red-X, Green-Y, Blue-Z)')
view([0 0 -1])  
guidata(hObject, handles);

% old codes about rotation of crystal are in v5 version.



% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
%according to the screen show to project crystal
[at,el] = view();
set(handles.edit91,'string',num2str(at));
set(handles.edit92,'string',num2str(el));
step1forshowcrystal(at, el, hObject, eventdata, handles);



function edit91_Callback(hObject, eventdata, handles)
function edit91_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nx=mycen(x);
if rem(x,2)==0;
    nx=x/2+1;
else
    nx=(x+1)/2;
end



function edit94_Callback(hObject, eventdata, handles)
function edit94_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit99_Callback(hObject, eventdata, handles)
function edit99_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit98_Callback(hObject, eventdata, handles)
function edit98_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)

function Untitled_24_Callback(hObject, eventdata, handles)
prompt='Atom Size';
dlg_title='Atomsize( constant>0 ):';
num_lines = 1;%输入对话框的行数;
default_val={ num2str(handles.dis_atomsize) };%默认的值;
answer=inputdlg(prompt,dlg_title,num_lines,default_val);
handles.dis_atomsize=str2num(answer{1});

guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
Drawsupercell(hObject, eventdata, handles)

function edit100_Callback(hObject, eventdata, handles)
function edit100_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit101_Callback(hObject, eventdata, handles)
function edit101_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_26_Callback(hObject, eventdata, handles)  %显示分层的图像
cd(handles.tempcd)  %进入上一次打开的目录
[FileName,PathName]=uigetfile({'*.mat','Slicing Crystal File (*.mat)'});%打开文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    cd(handles.execd)
    return;
end
load(strcat(PathName, FileName));
cd(handles.execd)
%load allresul.mat
%画出三维的slicing的平面；
Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);  %
hold on;   

maxx = max(allresul.x_new(:,2));
maxy = max(allresul.x_new(:,3));
[yy,xx] = meshgrid(linspace(0, maxy, max(10, round( maxy/ (20* handles.dis_atomsize)))), linspace(0, maxx, max(10, round( maxx/ (20* handles.dis_atomsize)))));
for i=0:length( allresul.slicethick )-1
    mesh(xx,yy,i*allresul.eachthick.*ones(size(xx)));
end

%%再多画一层底面，黑色的，作为出射面；
surf(xx,yy,(i+1)*allresul.eachthick.*ones(size(xx)), zeros(size(xx)));

disp('Project on the top slice')
title('Sices of 3D show');


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ypj
if exist('paraslicing.txt','file') == 0
    x=[10 0 0 0.2 1];
    save paraslicing.txt -ascii x
end
load paraslicing.txt
myvar.x_new = handles.x_new;
myvar.slicenum = paraslicing(1);  %这几个按钮都不要了。
myvar.top = paraslicing(2);
myvar.bottom = paraslicing(3);
myvar.x_mysavenewz = handles.x_mysavenewz;
disp(strcat('From top to bottom, distance is :', num2str(max(handles.x_mysavenewz)-min(handles.x_mysavenewz))))
%Kepeng
myvar.atomsize = paraslicing(4);
myvar.Projected = paraslicing(5);
% myvar.Pick = paraslicing(6);
%Kepeng
slicingdisplay(myvar);
guidata(hObject, handles);


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
if get(handles.radiobutton13,'value')==1 
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
end
if get(handles.radiobutton13,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM相关关闭
    
    set(handles.text21, 'visible', 'off');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture的mrad和1/nm单位
    set(handles.text134, 'visible', 'off');
    
    set(handles.text106, 'visible', 'on');  %STEM相关
    set(handles.text108, 'visible','on');
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit77, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
        
else
    set(handles.text106, 'visible', 'off');
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper的mrad和1/nm单位
    set(handles.text125, 'visible', 'off');
    set(handles.text24, 'visible', 'off');   % convergence & aperture的mrad和1/nm单位
    set(handles.text134, 'visible', 'on');
    
end
guidata(hObject, handles);



function edit105_Callback(hObject, eventdata, handles)
function edit105_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radiobutton14_Callback(hObject, eventdata, handles)



function edit106_Callback(hObject, eventdata, handles)
disp(strcat('Results will be saved in the folder:', handles.saveresult));
disp(strcat('Name as:', get(handles.edit106,'string'),'_*'));
guidata(hObject, handles);

function edit106_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
selpath = uigetdir();  % select one dir to save results
handles.saveresult = selpath;
disp(strcat('Results will be saved in the folder:', handles.saveresult));
disp(strcat('Name as:', get(handles.edit106,'string'),'_*'));
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_27_Callback(hObject, eventdata, handles)

prompt = {'Source Size (Angstrom)'; 'Sampling (Integer number)'};
dlg_title = 'Converlution for STEM images';
num_lines = 1;
def = {num2str(handles.conv_source);num2str(handles.conv_sampling)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer) 
    return; 
end
handles.conv_source=str2num(answer{1}); 
handles.conv_sampling=str2num(answer{2}); %把卷积的点源尺寸，以及放大率抽样


cd(handles.saveresult)  %进入上一次打开的目录
[FileName,PathName]=uigetfile({'*.dat','DATA File (*.dat)';'*.*','All File (*.*)'});%打开文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    cd(handles.execd)
    return;
end

cd(handles.execd);  %回到旧的文件夹
sx = str2num(get(handles.edit19, 'string'));
sy = str2num(get(handles.edit20, 'string'));

fid = fopen(strcat(PathName, FileName), 'r');
tempimage = fread(fid, [sy, sx], 'float');
fclose(fid);
tempimage= tempimage.';

sampling = str2num(get(handles.edit75, 'string'));% image sampling rate
step = str2num(get(handles.edit77, 'string'));  %step

if handles.conv_sampling>0
    [oldx, oldy]=size(tempimage);
    newx = round(oldx*handles.conv_sampling);
    newy = round(oldy*handles.conv_sampling);
    newimage = zeros(newx, newy);
    tempfft = fftshift(fft2(tempimage))/oldx/oldy;
    newimage(mycen(newx)-mycen(oldx)+1:mycen(newx)-mycen(oldx)+oldx, mycen(newy)-mycen(oldy)+1:mycen(newy)-mycen(oldy)+oldy) = tempfft;

    if handles.conv_source>0  %the input image sampling rate is sampling*step; the zoomed image sampling rate is sampling*step/handles.conv_sampling
        [x,y] = meshgrid( ([1:newy]-mycen(newy))/(newy*sampling*step/handles.conv_sampling) ,([1:newx]-mycen(newx))/(newx*sampling*step/handles.conv_sampling));
        ss = handles.conv_source/(sqrt(-log(0.5))*2);
        newimage = newimage.*exp(-((pi*ss)^2*(x.^2+y.^2)));
    else
        disp('Input a value (source size) larger than 0')
        return;
    end

    newimageresul = real( ifft2(ifftshift(newimage))*newx*newy);
else
    disp('Input a value (converlution sampling) larger than 0')
    return;
end

figure;imshow(newimageresul, 'XData',...
           str2num(get(handles.edit17,'string')) + [1:newy-1].*sampling*step/handles.conv_sampling, ...
           'YData', str2num(get(handles.edit18,'string')) + [0:newx-1].*sampling*step/handles.conv_sampling,...
          'DisplayRange',[]);axis on;
      
tempname = FileName; tempname(strfind(tempname,'_')) = ' ';
title(strcat('File name:', tempname))

fid=fopen(strcat(PathName, FileName(1:end-4),'_source_',num2str( handles.conv_source), '_zoom',num2str(handles.conv_sampling),'.dat'), 'w');
fwrite(fid, newimageresul.', 'float');
fclose(fid);
figure;surf(newimageresul)
shading interp
view([0 -1 0])
guidata(hObject, handles);
% 
% NyO = round(handles.Oversampling*Ny);%Ny,Nx为图象尺寸
% NxO = round(handles.Oversampling*Nx);
%     imgft = zeros(NyO,NxO);
%     NxOMid = floor(handles.Oversampling*Nx/2)+1;
%     NyOMid = floor(handles.Oversampling*Ny/2)+1;
%     NxMid = floor(Nx/2)+1;
%     NyMid = floor(Ny/2)+1;
%     imgft(NyOMid-NyMid+[1:Ny],NxOMid-NxMid+[1:Nx]) = fftshift(fft2(img));
%     % Apply effective source size for STEM images:
%     if ss > 0%ss与sorcesize有关,ss=sourcesize     ss = ss/(sqrt(-log(0.5))*2);
%         [qx,qy] = meshgrid((-NxOMid+[1:NxO])/(Nx*dx),(-NyOMid+[1:NyO])/(Ny*dy));
%         img = ifft2(ifftshift(imgft.*exp(-((pi*ss)^2*(qx.^2+qy.^2)))));        
%     else
%         img = ifft2(ifftshift(imgft));
%     end
%     if realFlag
%        img = real(img); 
%     end
%     scale = 1/(Nx*Ny);
%     Nx = round(handles.Oversampling*Nx);
%     Ny = round(handles.Oversampling*Ny);
%     dx = dx/handles.Oversampling;
%     dy = dy/handles.Oversampling;
%     scale = (scale*Nx*Ny);
% else
%     if ss > 0
%         NxMid = floor(Nx/2)+1;
%         NyMid = floor(Ny/2)+1;
%         [qx,qy] = meshgrid((-NxMid+[1:Nx])/(Nx*dx),(-NyMid+[1:Ny])/(Ny*dy));
%         img = real(ifft2(fft2(img).*ifftshift(exp(-((pi*ss)^2*(qx.^2+qy.^2))))));        
%     end    
% end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
if get(handles.popupmenu5,'value')==3  %如果采用lobato的系数，就不能分层
    set(handles.checkbox1, 'value', 0);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CPURB.
function CPURB_Callback(hObject, eventdata, handles)
% hObject    handle to CPURB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.CPURB, 'value')
    set(handles.GPURB,'value',0);
else
    set(handles.GPURB,'value',1);
end
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of CPURB


% --- Executes on button press in GPURB.
function GPURB_Callback(hObject, eventdata, handles)
% hObject    handle to GPURB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.GPURB, 'value')
    set(handles.CPURB,'value',0);
else
    set(handles.CPURB,'value',1);
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of GPURB


% --- Executes during object creation, after setting all properties.
function pushbutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Untitled_28_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %进入上一次打开的目录
[FileName,PathName]=uigetfile({'*.cif','CIF File (*.cif)'});%打开文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

lenxyz = [0,  0,  0];
prompt={sprintf(strcat('Cube''s width heigh thickness',num2str(lenxyz)))};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0';'0';'0'};

lenxyz=inputdlg(prompt,dlg_title,num_lines,defaultans);
lenxyz = str2num( cell2mat(lenxyz) );       

        
    z=unique(handles.x(:,1));
    for i=1:length(z)
        prompt={sprintf('DW factor for the %.0f Atom ',z(i))};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'0'};
        DW=inputdlg(prompt,dlg_title,num_lines,defaultans);
        handles.x(find(handles.x(:,1)==z(i)),7)=str2num(cell2mat(DW));
    end
    temp = handles.x;
    xlswrite(strcat(PathName, FileName(1:end-4)),temp);
    disp(strcat('Parameters of atoms are saved in the file', strcat(PathName, FileName(1:end-4)),'.xlsx'));


cd(handles.execd);  %回到旧的文件夹
handles.tempcd=PathName;  %记录上一个打开的文件夹

%读取pdb或者cif的具体内容
if sum(FileName(end-2:end)=='xls')==3 || sum(FileName(end-3:end)=='xlsx')==4 %表示是pdb格式
    handles.x=xlsread(strcat(PathName, FileName));
else
    disp('Cannot read this file!');
    return;
end

atompos=handles.x(:,2:4);  % 单位变埃~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % 原子高度，最小值变为0，且增加一层真空层厚度，方便分层
                                                            % 并且保证不会位于第0层。因为这里最少是位于第二层，下两行取了floor的行数。
                                                            % 比如厚度是3,0-3-6-9；加了1.5*eachslice后，值为4.5。分在第一层正中间

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %重新修改一下原子的坐标，其他信息需要查看handles.x里面读取的结果；
                %只有handles.x中的坐标需要改变一下；只有在设置层厚时候；
handles.x(:,2:4)=atompos;  %重新给原子的坐标，保证都是大于0的，且有边界
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %把旋转前的结构记录下来；handles.x_new记录的是旋转后
Drawsupercell(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_29_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_30_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_31_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_32_Callback(hObject, eventdata, handles)
%导入波函数，需要尺寸信息。
[FileName,handles.batchs.PathName]=uigetfile({'*.dat','Data Only Format (*.dat)'});%打开后缀为txt的文件
if length(FileName)==1 & FileName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    mycddir(handles.execd)
    return;
end
%选择波函数后，波函数的图像要显示
    prompt={'Enter wave size: width (Cubic)';...
       'Input flag: Show wave(S) and Run simulation(R)'; 'Images'' number of ';'Process number';...
         'Batch size for one process';...
         'TOP-LEFT coordinate (x-y)';...  %需要一行内输入两个值
         'Intercept size (width & hight)' };   %需要一行内输入两个值
    name='Input Image Info.';
    numlines=1;
    defaultanswer={num2str(handles.batchs.wavesx),...
        handles.batchs.showorrun,...
        num2str(handles.batchs.totalnumber),...
        num2str(handles.batchs.processnum),...
        num2str(handles.batchs.GPUBatch),...
        num2str(handles.batchs.top_left),...
        num2str(handles.batchs.width_heigh)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    options.Resize='on';
    options.WindowStyle='normal';

    handles.batchs.wavesx=str2num(answer{1});handles.batchs.wavesy=handles.batchs.wavesx;
    handles.batchs.showorrun = answer{2};
    handles.batchs.totalnumber = str2num(answer{3});
    handles.batchs.processnum = str2num(answer{4});
    handles.batchs.GPUBatch = str2num(answer{5});
    handles.batchs.top_left = str2num(answer{6});
    handles.batchs.width_heigh = str2num(answer{7});
    
   

fid=fopen(strcat(handles.batchs.PathName,FileName),'r');  %读取波函数。读石墨烯  %1204把sampling改为0.1，不用0.03，计算太慢，图像大小是320
a=fread(fid,[handles.batchs.wavesx*2,handles.batchs.wavesy],'float');%
fclose(fid);
disimage=a(1:2:end,:)+sqrt(-1)*a(2:2:end,:);   %
mywave=disimage.';

if handles.batchs.showorrun == 'S'  %只是先看一下所选的区域是否合适。
    %把波函数显示出来，并且画出截取的区域位置
    figure;imshow(angle(disimage),[]); title('phase of wave');   %把波函数显示出来，并且要截取的区域，左上角用红圈，右下角用黄星表示
    hold on; plot(handles.batchs.top_left(2),handles.batchs.top_left(1),'ro');   
    plot(handles.batchs.top_left(2)+handles.batchs.width_heigh(1),handles.batchs.top_left(1)+handles.batchs.width_heigh(2),'y*');
    %显示对角线
    line([handles.batchs.top_left(2), handles.batchs.top_left(1)], [handles.batchs.top_left(2)+handles.batchs.width_heigh(1), handles.batchs.top_left(1)+handles.batchs.width_heigh(2)],'color','r')
    %显示四条边
    axis on
    return;
end

%存储图像的文件夹名称
handles.batchs.savePathName=uigetdir();


%部分计算和读数据开始
sx=handles.batchs.wavesx;
sy=handles.batchs.wavesy;
sxy=handles.batchs.wavesx

%提取一些常数，与电压、抽样率无关
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
handles.lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %计算波长，单位nm
handles.vol=vol*1000;
handles.sampling=str2num(get(handles.edit75,'string'));  %埃/pixel  

%%绿色的尺寸，绿色倒格矢。gx用于计算原子位置的，s2使用peng的参数计算原子势场分布
[gx_green,gy_green]=meshgrid((-sxy/2:(sxy/2-1))./(sxy*handles.sampling), ...
    (-sxy/2:(sxy/2-1))./(sxy*handles.sampling)); %REPRO单位是1/A)^2,不是(1/nm)^2
sx_green=gx_green/2;  % peng 的s参数
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;

mywave=fftshift(fft2(mywave));   %构造波函数，如果从外界读取也可以;到频率空间


%读取像散的参数。
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %画CTF和PhasePlate用固定的大小即可
U=U-para_part1.tiltx;  %add 20210119 the incident wave is tilted but not the CTF is tilted
V=V-para_part1.tilty;  %the tilted CTF is shown for CTF display. but in simulation, only a tilted incident beam is required.
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %只用一个函数来把框格中的数据都读取了，这样保证其他地方需要读取数据时候，用统一的代码

myap=ones(size(U));
myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);  
handles.gfsf=para_part1.gfsf;

%为每个数值设置范围
%读取A1-等其他数值的幅度%以及角度范围
[txtName,txtPathName]=uigetfile({'*.txt','TXT Only Format (*.txt)'});%打开后缀为txt的文件
if length(txtName)==1 & txtName(1)==0;  %如果选择cancal没打开文件的话，程序终止
    mycddir(handles.execd)
    return;
else
    %读取所有参数的数值
%newpara_part2.phiA1=newpara_part2.phiA1+(rand-0.5)*2*myrandrange.phiA1;
%%代码变成，所有的角度，也是在一定范围内，加减一个数值。如果加减180，就是正负180绕一圈
allvalue = load(strcat(txtPathName,txtName));
myrandrange.focus=allvalue(1); %4.8;%7.5; %units: nm will be chanced with both positive and nagative directions
myrandrange.A1=allvalue(2);%0.96;%6; %units:nm
myrandrange.phiA1=allvalue(3);%0.96;%6; %units:nm
myrandrange.A2=allvalue(4);%10;%49.268;%30; %units:nm
myrandrange.phiA2=allvalue(5);%0.96;%6; %units:nm
myrandrange.B2=allvalue(6);%10;%15;%16.4;%15; %units:nm ****  有加*的都是0709修改的范围
myrandrange.phiB2=allvalue(7);%0.96;%6; %units:nm
myrandrange.Cs=allvalue(8);%1000; %3000;%2220;%3000; %unit nm  ****
myrandrange.A3=allvalue(9);%700; %50;%2220;%50; %units:nm  ****
myrandrange.phiA3=allvalue(10);%0.96;%6; %units:nm
myrandrange.S3=allvalue(11);%600; %100;%558;%100; %units:nm ****
myrandrange.phiS3=allvalue(12);%allvalue(3);%0.96;%6; %units:nm
myrandrange.A4=allvalue(13);%50000; %5000;%94800;%5000; %units:nm ****
myrandrange.phiA4=allvalue(14);%0.96;%6; %units:nm
myrandrange.D4=allvalue(15);10000; %20000;%19000;%20000; %units:nm ****
myrandrange.phiD4=allvalue(16);%0.96;%6; %units:nm
myrandrange.B4=allvalue(17);%40000; %5000;%19000;%5000; %units:nm ****
myrandrange.phiB4=allvalue(18);%0.96;%6; %units:nm
myrandrange.A5=allvalue(19);%1000000; %5000000;%3800000;%5000000; %units:nm ****
myrandrange.phiA5=allvalue(20);%0.96;%6; %units:nm
myrandrange.C5=allvalue(21);%1000000; %2000000;%3800000;%2000000; %units:nm ****  %其实也可以设置百分比。

% myrandrange.focus=4;%4.8;%7.5; %units: nm will be chanced with both positive and nagative directions
% myrandrange.A1=2.5;%0.96;%6; %units:nm
% myrandrange.A2=10;%49.268;%30; %units:nm
% myrandrange.B2=10;%15;%16.4;%15; %units:nm ****  有加*的都是0709修改的范围
% myrandrange.Cs=1000; %3000;%2220;%3000; %unit nm  ****
% myrandrange.A3=700; %50;%2220;%50; %units:nm  ****
% myrandrange.S3=600; %100;%558;%100; %units:nm ****
% myrandrange.A4=50000; %5000;%94800;%5000; %units:nm ****
% myrandrange.D4=10000; %20000;%19000;%20000; %units:nm ****
% myrandrange.B4=40000; %5000;%19000;%5000; %units:nm ****
% myrandrange.A5=1000000; %5000000;%3800000;%5000000; %units:nm ****
% myrandrange.C5=1000000; %2000000;%3800000;%2000000; %units:nm ****  %其实也可以设置百分比。

% for k=1:91
%      A1_phi(k)=-184+4*k;  
% end
% A2_phi=A1_phi;
% B2_phi=A1_phi;
% 
% for k=1:37
%      A3_phi(k)=-190+10*k;
% end
% S3_phi=A3_phi;
% A4_phi=A3_phi;
% B4_phi=A3_phi;
% D4_phi=A3_phi;
% A5_phi=A3_phi;

end
%rand_value=0*rand_value;

yita=str2num(get(handles.edit_Spr,'String')); 


%随机产生参数的变化
beginnumold=0;
tic

    
myap=ones(size(gx_green));
myap=myaperture(myap,1/(handles.sampling.*sxy),1/(handles.sampling.*sxy),0,0.1*str2num(get(handles.edit_Ape,'string')),0,0,1);
mywave=mywave.*myap;

commpara=[para_part1.sampling, vol, ... %抽样率，kv单位的电压
        para_part1.gmax, yita, para_part1.mtotal, para_part1.alafa*1000, ...  %最高频率，focus spread，计算的高斯函数的个数， beam convergence
        para_part1.tilt*1000, -para_part1.phitilt];
    
rand('seed',sum(100*clock));

othervar = [handles.batchs.processnum, handles.batchs.GPUBatch,...
    handles.batchs.totalnumber, handles.batchs.top_left, handles.batchs.width_heigh];
savePathName = handles.batchs.savePathName;
save('pydata.mat', 'para_part1', 'para_part2', 'myrandrange', ...
    'U', 'V', 'mywave', 'commpara', 'othervar', 'savePathName');
dos("python cuda_stem64.py");

% neixunhuan=6;
% for j=1:100
%     beginnum=beginnumold+(j-1)*neixunhuan;
%     %allimage = zeros(28*44,neixunhuan);
%     %allpara = zeros(33, neixunhuan);
%     
%   %  rand_value=(rand(1000,21)-0.5)*2;  %控制产生的个数，与下2行一起改
%     for i=1:neixunhuan; %length(rand_value(:,1))
%         newpara_part2=para_part2;
%         %newpara_part2.focus=round(newpara_part2.focus+(rand-0.5)*2*myrandrange.focus,1);  %加入随机数，并且保留小数点后两位
%         newpara_part2.focus=newpara_part2.focus+(rand-0.5)*2*myrandrange.focus; %add at 20200820，实验像
%         %newpara_part2.A1=round(newpara_part2.A1+(rand-0.5)*2*myrandrange.A1,1);  %测试像
%         newpara_part2.A1=newpara_part2.A1+(rand-0.5)*2*myrandrange.A1; %add at 20200820，实验像
%         newpara_part2.phiA1=A1_phi(round(rand*(90))+1); %实验像
%         %newpara_part2.phiA1=round(rand*180); %add at 20200820，测试像
%         %newpara_part2.A2=round(newpara_part2.A2+(rand-0.5)*2*myrandrange.A2,0);  %测试像
%         newpara_part2.A2=newpara_part2.A2+(rand-0.5)*2*myrandrange.A2; %add at 20200820，实验像
%         newpara_part2.phiA2=A2_phi(round(rand*(90))+1); %实验像
%         %newpara_part2.phiA2=round(rand*180); %add at 20200820
%         %newpara_part2.B2=round(newpara_part2.B2+(rand-0.5)*2*myrandrange.B2,0);  %加入随机数，并且保留小数点后1位
%         newpara_part2.B2=newpara_part2.B2+(rand-0.5)*2*myrandrange.B2; %add at 20200820,实验像
%         newpara_part2.phiB2=B2_phi(round(rand*(90))+1); %实验像
%         %newpara_part2.phiB2=round(rand*180); %add at 20200820，测试像
%         %newpara_part2.Cs=round(newpara_part2.Cs+(rand-0.5)*2*myrandrange.Cs,-2);  %加入随机数，并且保留小数点后1位
%         newpara_part2.Cs=newpara_part2.Cs+(rand-0.5)*2*myrandrange.Cs; %add at 20200820，实验像
%         %newpara_part2.A3=round(newpara_part2.A3+(rand-0.5)*2*myrandrange.A3,-2);  %加入随机数，并且保留小数点后1位
%         newpara_part2.A3=newpara_part2.A3+(rand-0.5)*2*myrandrange.A3; %add at 20200820，实验像
%         newpara_part2.phiA3=A3_phi(round(rand*(36))+1); %实验像
%         %newpara_part2.phiA3=round(rand*180); %add at 20200820，测试像
%         %newpara_part2.S3=round(newpara_part2.S3+(rand-0.5)*2*myrandrange.S3,-2);  %加入随机数，并且保留小数点后1位
%         newpara_part2.S3=newpara_part2.S3+(rand-0.5)*2*myrandrange.S3; %add at 20200820，实验像
%         newpara_part2.phiS3=S3_phi(round(rand*(36))+1); %实验像
%         %newpara_part2.phiS3=round(rand*180); %add at 20200820，测试像
% %不变的参数，就不加随机        %newpara_part2.A4=round(newpara_part2.A4+(rand-0.5)*2*myrandrange.A4,-3);  %加入随机数，并且保留小数点后1位
%         %newpara_part2.phiA4=A4_phi(round(rand*(36))+1);
%         %newpara_part2.D4=round(newpara_part2.D4+(rand-0.5)*2*myrandrange.D4,-3);  %加入随机数，并且保留小数点后1位
%         %newpara_part2.phiD4=D4_phi(round(rand*(36))+1);
%         %newpara_part2.B4=round(newpara_part2.B4+(rand-0.5)*2*myrandrange.B4,-3);  %加入随机数，并且保留小数点前4位
%         newpara_part2.B4=newpara_part2.B4+(rand-0.5)*2*myrandrange.B4; %add at 20200820，实验像
%         newpara_part2.phiB4=B4_phi(round(rand*(36))+1); %实验像
%         %newpara_part2.phiB4=round(rand*180); %测试像
%        %newpara_part2.A5=round(newpara_part2.A5+(rand-0.5)*2*myrandrange.A5,-5);  %加入随机数，并且保留小数点前4位
%       %newpara_part2.phiA5=A5_phi(round(rand*(36))+1);
%         %newpara_part2.C5=round(newpara_part2.C5+(rand-0.5)*2*myrandrange.C5,-5);  %加入随机数，并且保留小数点前4位
% 
%     
%       newpara_part2.A1x=newpara_part2.A1*cos(newpara_part2.phiA1/180*pi);   newpara_part2.A1y=newpara_part2.A1*sin(newpara_part2.phiA1/180*pi);
%     %C1=focus;  focus是在输入的defocs附近有偏离的
%       newpara_part2.A2x=newpara_part2.A2*cos(newpara_part2.phiA2/180*pi);   newpara_part2.A2y=newpara_part2.A2*sin(newpara_part2.phiA2/180*pi);
%       newpara_part2.B2x=newpara_part2.B2*cos(newpara_part2.phiB2/180*pi);   newpara_part2.B2y=newpara_part2.B2*sin(newpara_part2.phiB2/180*pi);
%         newpara_part2.A3x=newpara_part2.A3*cos(newpara_part2.phiA3/180*pi);   newpara_part2.A3y=newpara_part2.A3*sin(newpara_part2.phiA3/180*pi);
%        newpara_part2.S3x=newpara_part2.S3*cos(newpara_part2.phiS3/180*pi);   newpara_part2.S3y=newpara_part2.S3*sin(newpara_part2.phiS3/180*pi);
%        newpara_part2.C3=newpara_part2.Cs;
%       newpara_part2.A4x=newpara_part2.A4*cos(newpara_part2.phiA4/180*pi);   newpara_part2.A4y=newpara_part2.A4*sin(newpara_part2.phiA4/180*pi);
%       newpara_part2.D4x=newpara_part2.D4*cos(newpara_part2.phiD4/180*pi);   newpara_part2.D4y=newpara_part2.D4*sin(newpara_part2.phiD4/180*pi);
%       newpara_part2.B4x=newpara_part2.B4*cos(newpara_part2.phiB4/180*pi);   newpara_part2.B4y=newpara_part2.B4*sin(newpara_part2.phiB4/180*pi);
%       newpara_part2.A5x=newpara_part2.A5*cos(newpara_part2.phiA5/180*pi);   newpara_part2.A5y=newpara_part2.A5*sin(newpara_part2.phiA5/180*pi);
%       newpara_part2.S5x=newpara_part2.S5*cos(newpara_part2.phiS5/180*pi);   newpara_part2.S5y=newpara_part2.S5*sin(newpara_part2.phiS5/180*pi);
%        newpara_part2.C5=newpara_part2.C5;
%        newpara_part2.D5x=newpara_part2.D5*cos(newpara_part2.phiD5/180*pi);   newpara_part2.D5y=newpara_part2.D5*sin(newpara_part2.phiD5/180*pi);
%    
%        % para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);%如果测试面上的数据，就用这一条
%        % newpara_part2=para_part2; %如果测试面上的数据，就用这一条
%        
%       %  residual_aberr=WholeTCC_residual_newTEM(newpara_part2, U,V);  %计算phase plate visualing the residual aberrations
%         mytcc=WholeTCC2D_newTEM_forpar(para_part1, newpara_part2, U,V);
%     
% 
%        %像的初始化
%        myinten=zeros(sxy,sxy);
%        %波函数传播
% %        for nn=1:para_part1.mtotal;
% %           myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn).*myap))).^2;
% %        end
%        for nn=1:para_part1.mtotal;
%           myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn)))).^2;  %在myap后面直接把myap乘上去了
%        end
%        
%        savimage=zeros(28,44,3);
%        myresul=myinten(124:151,145:188); 
%        myresul=(myresul-min(myresul(:)))./(max(myresul(:))-min(myresul(:)));
%        savimage(:,:,1)=myresul;
%        savimage(:,:,2)=myresul;
%        savimage(:,:,3)=myresul;
%        %20201208新增
%        imwrite(savimage, strcat(num2str(beginnum+i),'.jpg'));
%        
%        %allimage(:,i)= myresul(:);
%    
%  
% 
%     
% %     %加入调制sin
% %     %myresul = sin(myresul*255./10);  %%%%%%%%%%%%%%%%0619
% %     myresul = sin(myresul*25);  %%%%%%%%%%%%%%%%0619
% %     myresul = (myresul-min(myresul(:)))./(max(myresul(:))-min(myresul(:))); %%%%%%%%%%%%%%%%0619
% %    
% %     savimage=repmat(myresul,1,1,3);
% %     %savimage(:,:,1)=myresul;savimage(:,:,2)=myresul;savimage(:,:,3)=myresul;
% %     imwrite(savimage,strcat('C:\Users\EM_Lab\Desktop\cos_0820_2i\',num2str(beginnum+i),'.jpg')); 
%     
%     myallpara=[commpara, ...
%         newpara_part2.focus, ...
%         newpara_part2.A1, -newpara_part2.phiA1, ...
%         newpara_part2.A2, -newpara_part2.phiA2, ...
%         newpara_part2.B2, -newpara_part2.phiB2, ...
%         newpara_part2.Cs/1000, ...
%         newpara_part2.A3/1000, -newpara_part2.phiA3, ...
%         newpara_part2.S3/1000, -newpara_part2.phiS3, ...
%         newpara_part2.A4/1000, -newpara_part2.phiA4, ...
%         newpara_part2.D4/1000, -newpara_part2.phiD4, ...        
%         newpara_part2.B4/1000, -newpara_part2.phiB4, ...
%         newpara_part2.C5/1000000, ...
%         newpara_part2.A5/1000000, -newpara_part2.phiA5, ...
%         newpara_part2.D5/1000000, -newpara_part2.phiD5, ...        
%         newpara_part2.S5/1000000, -newpara_part2.phiS5];  %通过tccread命令的参数，角度需要取负
%         
%         %allpara(:,i) = myallpara;
%        
%         fid=fopen(strcat(num2str(beginnum+i),'.para'),'w');
%         fwrite(fid, myallpara, 'float');
%         fclose(fid);
%     end
%         %fid=fopen(strcat('H:\goku\1204\batch3\',num2str(j),'img.dat'),'w');
%         %fwrite(fid, allimage, 'float');
%         %fclose(fid);    
%         %fid=fopen(strcat('H:\goku\1204\batch3\',num2str(j),'.para'),'w');
%         %fwrite(fid, allpara, 'float');
%         %fclose(fid);
% end


%figure;imshow(myinten,[]);  %显示图像
toc
