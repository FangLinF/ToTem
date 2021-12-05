function varargout = slicingdisplay(varargin)
%这个面板只供显示分层的情况，主要只导出顶部和底部填充的高度，以及每层分层的厚度是多少
%this dialog provide the top and bottom extends, and the thickness of each
%slice
% And it give the user the detail about atoms in every slices

% SLICINGDISPLAY MATLAB code for slicingdisplay.fig
%      SLICINGDISPLAY, by itself, creates a new SLICINGDISPLAY or raises the existing
%      singleton*.
%
%      H = SLICINGDISPLAY returns the handle to a new SLICINGDISPLAY or the handle to
%      the existing singleton*.
%
%      SLICINGDISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLICINGDISPLAY.M with the given input arguments.
%
%      SLICINGDISPLAY('Property','Value',...) creates a new SLICINGDISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before slicingdisplay_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to slicingdisplay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help slicingdisplay


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @slicingdisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @slicingdisplay_OutputFcn, ...
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


% --- Executes just before slicingdisplay is made visible.
function slicingdisplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to slicingdisplay (see VARARGIN)

% Choose default command line output for slicingdisplay
handles.output = hObject;

handles.extendis=5;

handles.x_new=varargin{1}.x_new;  %display image
handles.top=varargin{1}.top;  %上表面
handles.bottom=varargin{1}.bottom;  %下表面
handles.slicenum=varargin{1}.slicenum;  %层数
handles.x_mysavenewz = varargin{1}.x_mysavenewz;

%Kepeng
handles.atomsize = varargin{1}.atomsize;%原子显示尺寸
handles.Projected = varargin{1}.Projected;
%Kepeng

set(handles.edit3, 'string', num2str(handles.slicenum));
set(handles.edit4, 'string', num2str(handles.top));
set(handles.edit5, 'string', num2str(handles.bottom));
set(handles.edit7, 'string', num2str(max(handles.x_mysavenewz)-min(handles.x_mysavenewz)));
set(handles.edit8, 'string', strcat('1:',num2str(handles.slicenum)) );
%Kepeng
set(handles.edit6, 'string', num2str(handles.atomsize));
set(handles.edit1, 'string', num2str(handles.Projected));

set(handles.edit10, 'string', num2str( (max(handles.x_mysavenewz)-min(handles.x_mysavenewz)+handles.top+handles.bottom)/handles.slicenum));  %thickness of each slice

uipushtool4_ClickedCallback(hObject, eventdata, handles)

guidata(hObject, handles);

% UIWAIT makes slicingdisplay wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function edit3_KeyPressFcn(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = slicingdisplay_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
varargin{1}.top = str2num(get(handles.edit4,'string'));
varargin{1}.bottom = str2num(get(handles.edit5,'string'));
varargin{1}.slicenum = str2num(get(handles.edit3,'string'));



function edit3_Callback(hObject, eventdata, handles)
handles.top = str2num(get(handles.edit4,'string'));
handles.bottom = str2num(get(handles.edit5,'string'));
handles.slicenum = str2num(get(handles.edit3,'string'));
set(handles.edit10, 'string', num2str( (max(handles.x_mysavenewz)-min(handles.x_mysavenewz)+handles.top+handles.bottom)/handles.slicenum));  %thickness of each slice
if str2num(get(handles.edit3, 'string'))<0;  %预留空间是整数
    errdlg(' Input a Positive Integer Number ')
end
slider1_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
handles.top = str2num(get(handles.edit4,'string'));
handles.bottom = str2num(get(handles.edit5,'string'));
handles.slicenum = str2num(get(handles.edit3,'string'));
set(handles.edit10, 'string', num2str( (max(handles.x_mysavenewz)-min(handles.x_mysavenewz)+handles.top+handles.bottom)/handles.slicenum));  %thickness of each slice

if str2num(get(handles.edit4, 'string'))<0;  %预留空间是整数
    errdlg(' Input a Positive Number ')
end


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
handles.top = str2num(get(handles.edit4,'string'));
handles.bottom = str2num(get(handles.edit5,'string'));
handles.slicenum = str2num(get(handles.edit3,'string'));
set(handles.edit10, 'string', num2str( (max(handles.x_mysavenewz)-min(handles.x_mysavenewz)+handles.top+handles.bottom)/handles.slicenum));  %thickness of each slice
if str2num(get(handles.edit5, 'string'))<0;  %预留空间是整数
    errdlg(' Input a Positive Number ')
end

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function Updateshow(hObject, eventdata, handles)





% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

numfrom = round(get(handles.slider1, 'value')*(str2num(get(handles.edit3,'string')) - str2num(get(handles.edit1,'string')))) +1;  %显示第几个slice
set(handles.edit2,'string', strcat( num2str(numfrom),'-',num2str(numfrom+str2num(get(handles.edit1,'string'))-1)));
uipushtool4_ClickedCallback(hObject, eventdata, handles)
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit6, 'string'))<0
    errdlg(' Input a Number [0,1]')
end


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
%画出三维的slicing的平面；显示每层的原子排布
maxx = max(handles.x_new(:,2));
maxy = max(handles.x_new(:,3));
maxz = max(handles.x_new(:,4));

top = str2num(get(handles.edit4,'string'));
bottom = str2num(get(handles.edit5,'string'));

handles.x_new(:,4) = handles.x_mysavenewz + top; %顶部加上一些空间，由用户定，这样可以保证划分原子时候，即使在某一层，也可以位于不同高度
eachthick=(max(handles.x_new(:,4)) + bottom)./str2num(get(handles.edit3,'string'));   % 单位变为埃
handles.slicethick = (0:str2num(get(handles.edit3,'string'))) * eachthick;

%画出三维的slicing的平面；
% Drawsupercell_figure(handles.x_new, str2num(get(handles.edit6, 'string')), handles.x_mysavenewz, top, bottom)
% view([2 1 0])
Drawsupercell_figure(handles.x_new, str2num(get(handles.edit6, 'string')), handles.x_mysavenewz, top, bottom); 
hold on;   
dis_atomsize = str2num(get(handles.edit6, 'string'));
%20201227修改网格的长宽，之前颠倒了。
[yy,xx] = meshgrid(linspace(0, maxy, max(10, round( maxy/ (20* dis_atomsize)))), linspace(0, maxx, max(10, round( maxx/ (20* dis_atomsize)))));
for i=0:str2num(get(handles.edit3,'string'))-1
    surf(xx,yy,i*eachthick.*ones(size(xx)));
end

%%再多画一层底面，黑色的，作为出射面；
surf(xx,yy,(i+1)*eachthick.*ones(size(xx)), zeros(size(xx)));
view([2 1 0])
hidden off
quiver3(-max(handles.x_new(:,2))/10,-max(handles.x_new(:,3))/10, 0, 0,0,max(handles.x_new(:,4))/2,  1, 'linewidth', 3)  %在前三个值坐标处，显示001方向的矢量
text(-max(handles.x_new(:,2))/10,-max(handles.x_new(:,3))/10, 0,'Beam')
disp('Project from the top slice')
title('3D Slicing')
allresul.x_new = handles.x_new;
allresul.eachthick = eachthick;
allresul.slicethick = handles.slicethick;
disp(strcat( 'Thickness of each slice is:', num2str(eachthick)))
save allresul allresul



% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject, eventdata, handles)  
%show atoms of specified slices in handels.axis1
top = str2num(get(handles.edit4,'string'));
bottom = str2num(get(handles.edit5,'string'));

drawatoms = handles.x_new(:,1:4) ; %元素种类和第几个需要画出来
drawatoms(:,4) = handles.x_mysavenewz + top; %顶部加上一些空间，由用户定，这样可以保证划分原子时候，即使在某一层，也可以位于不同高度
eachthick=(max(drawatoms(:,4)) + bottom)./str2num(get(handles.edit3,'string'));   % 单位变为埃

numfrom = round(get(handles.slider1, 'value')*(str2num(get(handles.edit3,'string')) - str2num(get(handles.edit1,'string')))) +1;  %显示第几个slice
numto = numfrom+str2num(get(handles.edit1,'string'))-1;

i = find(drawatoms(:,4)>=(numfrom-1).*eachthick & drawatoms(:,4)<numto.*eachthick);
drawatoms=drawatoms(i,:); %只保留符合要求的原子来画图
%显示每个面上的层内的平面原子；


cla(handles.axes1);   %清除原有图像
% hold off;
cla(handles.axes2);   %清除原有图像
% hold off;
%画图显示一下结构

% xrag=get(gca,'Xlim');  %记录当前图像显示时候，x和y的坐标范围
% yrag=get(gca,'Ylim');
   
item_ele=unique(drawatoms(:,1));  %考察有几种原子
for i=1:length(item_ele)  %元素的种类
    r=1-mod(item_ele(i),5)*0.25;    %根据原子序数，设置该原子特有的颜色，r,g,b为三基色的值
    g=1-mod(floor(item_ele(i)/5),5)*0.25;
    b=1-mod(floor(item_ele(i)/25),5)*0.25;
    %居中画图
    axes(handles.axes1);  %画原子
    plot(drawatoms(find(drawatoms(:,1)==item_ele(i)),2),...    
        - drawatoms(find(drawatoms(:,1)==item_ele(i)),3), 'o',...  %绘点，形状为‘o’  %20201230原子位置的坐标，必须y轴是负的，这样是从上往下看的角度
         'MarkerEdgeColor',[0 0 0],...    %设置边缘颜色为黑色
         'MarkerSize',str2num(get(handles.edit6, 'string')).*(20-log2(2)),...  %设置原子的大小,随着显示cell数目的增加大小会变小
         'MarkerFaceColor',[r g b]);    %设置原子的颜色   
    hold on;axis equal
     axes(handles.axes2); %画图显示一下图例
     plot(0.1*i,0.7,'o',...     %显示原子
     'MarkerEdgeColor',[0 0 0],...
     'MarkerSize',15,...
     'MarkerFaceColor',[r g b]);
     text(0.1*i-0.01,0.3,num2str(item_ele(i)));  %在原子下方显示对应的原子质子数
      hold on;
end
axis(handles.axes1);  %画原子
set(handles.axes1, 'XLim',[min(handles.x_new(:,2)),max(handles.x_new(:,2))],'YLim',[-max(handles.x_new(:,3)),-min(handles.x_new(:,3))+0.0001])
guidata(hObject, handles);

% 
% % --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
x=[];   
x(1) = str2num(get(handles.edit3,'string'));  %这几个按钮都不要了。
x(2) = str2num(get(handles.edit4,'string'));
x(3) = str2num(get(handles.edit5,'string'));
x(4) = str2num(get(handles.edit6,'string'));
x(5) = 1;
save paraslicing.txt -ascii x

delete(hObject);




function edit8_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
top = str2num(get(handles.edit4,'string'));
bottom = str2num(get(handles.edit5,'string'));
%每层的厚度
%所有层的原子被投影在上表面，那么这个上表面的位置是哪里
handles.x_new(:,4) = handles.x_mysavenewz + top; %顶部加上一些空间，由用户定，这样可以保证划分原子时候，即使在某一层，也可以位于不同高度
eachthick=max((max(handles.x_new(:,4)) + bottom)./str2num(get(handles.edit3,'string')),0.0001);   % 单位变为埃 20201227 修改
handles.slicethick = (0:str2num(get(handles.edit3,'string'))) * eachthick;  %把每层的厚度都标记出来，注意，底面的高度放在最后一个变量

%如果次序是 1 2 5 6，则堆叠成新的1 2 3 4层
num = str2num( get(handles.edit8, 'string') );  %得到需要拿来接着计算的序号是多少

%重新打造一个结构，它们只取这个原子的特定层拿来做新的层传播
tempatoms = []; 
for i = 1:length(num)
     jj = find(handles.x_new(:,4) >= (num(i)-1)*eachthick & handles.x_new(:,4) < num(i)*eachthick);  %找到厚度满足要求的原子序号
     if isempty(tempatoms)  %如果待记录保留原子的数据还是空的 
          tempatoms = handles.x_new(jj, :);
          tempatoms(:,4) = tempatoms(:,4) - (num(i)-1)*eachthick + (i-1)*eachthick; 
     else
          temp = handles.x_new(jj, :);
          temp(:,4) = temp(:,4) - (num(i)-1)*eachthick + (i-1)*eachthick; 
          tempatoms(end+1 :end+length(jj), :) = temp;
     end
end
handles.slicethick = (0:length(num)-1) * eachthick;  %新的层厚

%画出三维的slicing的平面；
if isempty(tempatoms) 
    msgbox('no atom');
    return;
end
% Drawsupercell_figure(tempatoms, str2num(get(handles.edit6, 'string')), tempatoms(:,4), 0, 0);  %
% view([2 1 0])
Drawsupercell_figure(tempatoms, str2num(get(handles.edit6, 'string')), tempatoms(:,4), 0, 0); 
hold on;   
dis_atomsize = str2num(get(handles.edit6, 'string'));
maxx = max(tempatoms(:,2));
maxy = max(tempatoms(:,3));  
[yy,xx] = meshgrid(linspace(0, maxy, max(10, round( maxy/ (20* dis_atomsize)))), linspace(0, maxx, max(10, round( maxx/ (20* dis_atomsize))))); %20201227修改
for i=0:length(num)-1
    surf(xx,yy,i*eachthick.*ones(size(xx)));
end

%%再多画一层底面，黑色的，作为出射面；
surf(xx,yy,(i+1)*eachthick.*ones(size(xx)), zeros(size(xx)));
view([2 1 0])
hidden off
quiver3(-max(tempatoms(:,2))/10,-max(tempatoms(:,3))/10, 0, 0,0,max(tempatoms(:,4))/2,  1, 'linewidth', 3)  %在前三个值坐标处，显示001方向的矢量
text(-max(tempatoms(:,2))/10,-max(tempatoms(:,3))/10, 0,'Beam')
disp('Project from the top slice')
title('Select some slices and 3D show');

allresul.x_new = tempatoms;
allresul.eachthick = eachthick;
allresul.slicethick = handles.slicethick;
disp(strcat( 'Thickness of each slice is:', num2str(eachthick)))
save allresul allresul



function edit9_Callback(hObject, eventdata, handles)
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
