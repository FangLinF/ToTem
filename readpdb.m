function resul=readpdb(filename);
%读取PDB文件中原子的位置等信息
%元素种类、坐标信息、占据、DW信息读取出来
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%判断filename最后三个字母为PDB
if (sum(filename(end-2:end)=='pdb')==3  || sum(filename(end-2:end)=='PDB')==3)
    fid=fopen(filename,'r');
    str=fread(fid,'char');
    fclose(fid);
%     k=1; %计算有多少个ATOM
%     for i=1:length(str)-4;   %找到ATOM的字符
%         if sum(char(str(i:i+3))==['A';'T';'O';'M'])==4;
%             ATOMcharPOS(k)=i;   %记录ATOM的字符的位置
%             k=k+1;
%         end   
%     end
    ATOMcharPOS=strfind(str','ATOM');  %%%%%%上一段直接简化的写法
    %ENDPos=strfind(str','TER');   %结束位置的位置
    
    %ATOM到ATOM之间必定是等间距的。
    InfoLenth=ATOMcharPOS(end)-ATOMcharPOS(end-1);   %存储的信息长度必定是这个数值，拿最后一个原子的信息长度为准
    tempstr=str(ATOMcharPOS(end-1)+4:ATOMcharPOS(end)-1)';    %拿最后这个原子信息来考察哪些字段开始可能是原子种类、原子坐标信息 
    
    %寻找到元素名称起始的偏离量
    tempNumpos=find(tempstr>='0'&tempstr<='9');  %这些字符是数字
    tempCappos=[find(tempstr>='a'&tempstr<='z'), find(tempstr>='A'&tempstr<='Z')];  %这些字符是字母
    tempCappos=sort(tempCappos);
    tempspacepos=findstr(tempstr,' '); %这些字符是空格
    tempspotpos=findstr(tempstr,'.');  %这些字符是 .
    
    %前面的结构如图~~~~~~~~~~~~~~~~~~~~
    %空格区 数字区 空格区 符号区 空格区 其他信息  数值信息
    %这一段要找到符号区前面的字符~~~~~~~~~~~~~~
    
    %找到数字后面跟着的元素名称起始和中止位置
    elenameshift=max(tempspacepos( find(tempspacepos<tempCappos(1)) ))   ;  %空格位置上，在第一个数字位置后、在第一个字母前的位置;另外，考虑上ATOM的4个字符
    elenameend=min(tempspacepos(find(tempspacepos>tempCappos(1))))   ;   %找到字母中，比字母开始靠后的
    %找到数据区的起始和中止位置
    valueshift=max(tempspacepos(find(tempspacepos<tempspotpos(1))))  ;  %找到第一个小数点数的前面的空格
    if tempCappos(end) < valueshift   %如果字母都在此段数值之前，表示valueshift之后全部都是数字，则数值边界设定为数据段结尾
        valueend=InfoLenth-4;
    else
        valueend=max(tempspacepos(find(tempspacepos<tempCappos(end)))) -1;   %最后一个字母前，一定是5个数据的截止处
    end
  
    %这个信息可以不要！只要有ATOM就一定是包括了原子的信息
%     %把最后一个TER的信息也记录进去
%     for i=length(str)-200:length(str)-2;   %找到ATOM的字符
%           if sum(char(str(i:i+2))==['T';'E';'R'])==3;
%              ATOMcharPOS(k)=i;   %记录结束字符的位置，所以，最后一个字节信息内，是不存在原子信息的
%              k=k+1;
%           end
%     end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~以上得到的ATOMcharPOS给出了“每行元素信息的行开始位置以及行结束位置”
    
    %给出三个坐标、占据率、振动因子的信息，与ATOM偏离的关系。
    
    
%     %判断至少有10个连续的字符是数字和.，这样才能算是真实的原子  %最后一个ATOM默认包含了原子信息
%     for kk=1:length(ATOMcharPOS)-1;
%         if sum(str(ATOMcharPOS(kk):ATOMcharPOS(kk+1)-1)>='0' & str(ATOMcharPOS(kk):ATOMcharPOS(kk+1)-1)<='9' ) <10   %如果少于十个数字，就一定不是原子信息
%             ATOMcharPOS(kk)=[];
%         end
%     end
%     for kk=1:length(ATOMcharPOS)-1;   %将一段内的元素种类、坐标信息、占据、DW信息读取出来
%         [ele(kk),x(kk),y(kk),z(kk),ocp(kk),dw(kk)]=getatominfo(str(ATOMcharPOS(kk):ATOMcharPOS(kk+1)-1));
%     end
%     resul=[ele;x;y;z;ocp;dw]';
        %输入每段信息，判断是否是原子
    for i=1:length(ATOMcharPOS)
        smallstr=str(ATOMcharPOS(i)+4:ATOMcharPOS(i)+InfoLenth-1)';
        if length(findstr(smallstr,'.'))>=5;   %%  如果是原子，这段信息里面包含三个点.理论上是应该有5个小数点的
            [ele(i),x(i),y(i),z(i),ocp(i),dw(i)]=GetAtomInfo(smallstr, elenameshift,elenameend,valueshift,valueend);
        end
    end
    resul=[ele;x;y;z;ocp;dw]';
    return;
else  %不是PDB文件则退出
    printf('Not PDB File');
    return;
end

pp=1;


function [ele,x,y,z,ocp,dw]=getatominfo(str);

numpos=find(str=='.');
spacepos=find(str==' ');
if length(numpos)==5   %有5个含小数点的数值
   st=max(spacepos(find(min(numpos)-spacepos>0)));
   en=min(spacepos(find(spacepos-max(numpos)>0)));
   temp=str2num(char(str(33:67)'));
   x=temp(1);y=temp(2);z=temp(3);ocp=temp(4);dw=temp(5);
   
   %求解元素的序号，从最后几个数值中读取吧，前面的字符可能包含了H这个比较尴尬的名字
   temp=str(end-20:end);
   ele=getatomelement(temp);
end

function [ele,x,y,z,ocp,dw]=GetAtomInfo(str, elenameshift,elenameend,valueshift,valueend);
% temp=str2num(char(str(valueshift:valueend)));
charfirst=char(str(valueshift:valueend));
spacepos=max(find(charfirst=='-'));
charfirst(spacepos+1:end+1)=charfirst(spacepos:end);
charfirst(spacepos)=(' ');
temp=str2num(charfirst);
x=temp(1);y=temp(2);z=temp(3);ocp=temp(4);dw=temp(5);
   
   %求解元素的序号，从最后几个数值中读取吧，前面的字符可能包含了H这个比较尴尬的名字
ele=getatomelement(str(elenameshift:elenameend));
pp=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%这里帮忙补齐所有元素的名称
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ele=getatomelement(str);   %从字符串中找到108个元素的序号名称，输出量为原子序号
strr=char(str()); %转成字符串格式
%在strr字符串中寻找特定的元素名称
strr(find(isspace(strr)))=[];%去空格
allelementA(1,:)=' H';
allelementA(2,:)='He';
allelementA(3,:)='Li';
allelementA(4,:)='Be';
allelementA(5,:)=' B';
allelementA(6,:)=' C';
allelementA(7,:)=' N';
allelementA(8,:)=' O';
allelementA(9,:)=' F';
allelementA(10,:)='Ne';
allelementA(11,:)='Na';
allelementA(12,:)='Mg';
allelementA(13,:)='Al';
allelementA(14,:)='Si';
allelementA(15,:)=' P';
allelementA(16,:)=' S';
allelementA(17,:)='Cl';
allelementA(18,:)='Ar';
allelementA(19,:)=' K';
allelementA(20,:)='Ca';
allelementA(21,:)='Sc';
allelementA(22,:)='Ti';
allelementA(23,:)=' V';
allelementA(24,:)='Cr';
allelementA(25,:)='Mn';
allelementA(26,:)='Fe';
allelementA(27,:)='Co';
allelementA(28,:)='Ni';
allelementA(29,:)='Cu';
allelementA(30,:)='Zn';
allelementA(31,:)='Ga';
allelementA(32,:)='Ge';
allelementA(33,:)='As';
allelementA(34,:)='Se';
allelementA(35,:)='Br';
allelementA(36,:)='Kr';
allelementA(37,:)='Rb';
allelementA(38,:)='Sr';
allelementA(39,:)=' Y';
allelementA(40,:)='Zr';
allelementA(41,:)='Nb';
allelementA(42,:)='Mo';
allelementA(43,:)='Tc';
allelementA(44,:)='Ru';
allelementA(45,:)='Rh';
allelementA(46,:)='Pd';
allelementA(47,:)='Ag';
allelementA(48,:)='Cd';
allelementA(49,:)='In';
allelementA(50,:)='Sn';
allelementA(51,:)='Sb';
allelementA(52,:)='Te';
allelementA(53,:)=' I';
allelementA(54,:)='Xe';
allelementA(55,:)='Cs';
allelementA(56,:)='Ba';
allelementA(57,:)='La';
allelementA(58,:)='Ce';
allelementA(59,:)='Pr';
allelementA(60,:)='Nd';
allelementA(61,:)='Pm';
allelementA(62,:)='Sm';
allelementA(63,:)='Eu';
allelementA(64,:)='Gd';
allelementA(65,:)='Tb';
allelementA(66,:)='Dy';
allelementA(67,:)='Ho';
allelementA(68,:)='Er';
allelementA(69,:)='Tm';
allelementA(70,:)='Yb';
allelementA(71,:)='Lu';
allelementA(72,:)='Hf';
allelementA(73,:)='Ta';
allelementA(74,:)=' W';
allelementA(75,:)='Re';
allelementA(76,:)='Os';
allelementA(77,:)='Ir';
allelementA(78,:)='Pt';
allelementA(79,:)='Au';
allelementA(80,:)='Hg';
allelementA(81,:)='Tl';
allelementA(82,:)='Pb';
allelementA(83,:)='Bi';
allelementA(84,:)='Po';
allelementA(85,:)='At';
allelementA(86,:)='Rn';
allelementA(87,:)='Fr';
allelementA(88,:)='Ra';
allelementA(89,:)='Ac';
allelementA(90,:)='Th';
allelementA(91,:)='Pa';
allelementA(92,:)=' U';
allelementA(93,:)='Nq';
allelementA(94,:)='Pu';
allelementA(95,:)='Am';
allelementA(96,:)='Cm';
allelementA(97,:)='Bk';
allelementA(98,:)='Cf';
allelementA(99,:)='Es';
allelementA(100,:)='Fm';
allelementA(101,:)='Md';
allelementA(102,:)='No';
allelementA(103,:)='Lr';
allelementA(104,:)='Rf';
allelementA(105,:)='Db';
allelementA(106,:)='Sg';
allelementA(107,:)='Bh';
allelementA(108,:)='Hs';
allelementA(109,:)='Mt';
allelementA(110,:)='Ds';
allelementA(111,:)='Rg';

allelementB(1,:)=' H';
allelementB(2,:)='HE';
allelementB(3,:)='LI';
allelementB(4,:)='BE';
allelementB(5,:)=' B';
allelementB(6,:)=' C';
allelementB(7,:)=' N';
allelementB(8,:)=' O';
allelementB(9,:)=' F';
allelementB(10,:)='NE';
allelementB(11,:)='NA';
allelementB(12,:)='MG';
allelementB(13,:)='AL';
allelementB(14,:)='SI';
allelementB(15,:)=' P';
allelementB(16,:)=' S';
allelementB(17,:)='CL';
allelementB(18,:)='AR';
allelementB(19,:)=' K';
allelementB(20,:)='CA';
allelementB(21,:)='SC';
allelementB(22,:)='TI';
allelementB(23,:)=' V';
allelementB(24,:)='CR';
allelementB(25,:)='MN';
allelementB(26,:)='FE';
allelementB(27,:)='CO';
allelementB(28,:)='NI';
allelementB(29,:)='CU';
allelementB(30,:)='ZN';
allelementB(31,:)='GA';
allelementB(32,:)='GE';
allelementB(33,:)='AS';
allelementB(34,:)='SE';
allelementB(35,:)='BR';
allelementB(36,:)='KR';
allelementB(37,:)='RB';
allelementB(38,:)='SR';
allelementB(39,:)=' Y';
allelementB(40,:)='ZR';
allelementB(41,:)='NB';
allelementB(42,:)='MO';
allelementB(43,:)='TC';
allelementB(44,:)='RU';
allelementB(45,:)='RH';
allelementB(46,:)='PD';
allelementB(47,:)='AG';
allelementB(48,:)='CD';
allelementB(49,:)='IN';
allelementB(50,:)='SN';
allelementB(51,:)='SB';
allelementB(52,:)='TE';
allelementB(53,:)=' I';
allelementB(54,:)='XE';
allelementB(55,:)='CS';
allelementB(56,:)='BA';
allelementB(57,:)='LA';
allelementB(58,:)='CE';
allelementB(59,:)='PR';
allelementB(60,:)='ND';
allelementB(61,:)='PM';
allelementB(62,:)='SM';
allelementB(63,:)='EU';
allelementB(64,:)='GD';
allelementB(65,:)='TB';
allelementB(66,:)='DY';
allelementB(67,:)='HO';
allelementB(68,:)='ER';
allelementB(69,:)='TM';
allelementB(70,:)='YB';
allelementB(71,:)='LU';
allelementB(72,:)='HF';
allelementB(73,:)='TA';
allelementB(74,:)=' W';
allelementB(75,:)='RE';
allelementB(76,:)='OS';
allelementB(77,:)='IR';
allelementB(78,:)='PT';
allelementB(79,:)='AU';
allelementB(80,:)='HG';
allelementB(81,:)='TL';
allelementB(82,:)='PB';
allelementB(83,:)='BI';
allelementB(84,:)='PO';
allelementB(85,:)='AT';
allelementB(86,:)='RN';
allelementB(87,:)='FR';
allelementB(88,:)='RA';
allelementB(89,:)='AC';
allelementB(90,:)='TH';
allelementB(91,:)='PA';
allelementB(92,:)=' U';
allelementB(93,:)='NQ';
allelementB(94,:)='PU';
allelementB(95,:)='AM';
allelementB(96,:)='CM';
allelementB(97,:)='BK';
allelementB(98,:)='CF';
allelementB(99,:)='ES';
allelementB(100,:)='FM';
allelementB(101,:)='MD';
allelementB(102,:)='NO';
allelementB(103,:)='LR';
allelementB(104,:)='RF';
allelementB(105,:)='DB';
allelementB(106,:)='SG';
allelementB(107,:)='BH';
allelementB(108,:)='HS';
allelementB(109,:)='MT';
allelementB(110,:)='DS';
allelementB(111,:)='RG';


allelementC(1,:)=' H';
allelementC(2,:)='hE';
allelementC(3,:)='lI';
allelementC(4,:)='bE';
allelementC(5,:)=' B';
allelementC(6,:)=' C';
allelementC(7,:)=' N';
allelementC(8,:)=' O';
allelementC(9,:)=' F';
allelementC(10,:)='nE';
allelementC(11,:)='nA';
allelementC(12,:)='mG';
allelementC(13,:)='aL';
allelementC(14,:)='sI';
allelementC(15,:)=' P';
allelementC(16,:)=' S';
allelementC(17,:)='cL';
allelementC(18,:)='aR';
allelementC(19,:)=' K';
allelementC(20,:)='cA';
allelementC(21,:)='sC';
allelementC(22,:)='tI';
allelementC(23,:)=' V';
allelementC(24,:)='cR';
allelementC(25,:)='mN';
allelementC(26,:)='fE';
allelementC(27,:)='cO';
allelementC(28,:)='nI';
allelementC(29,:)='cU';
allelementC(30,:)='zN';
allelementC(31,:)='gA';
allelementC(32,:)='gE';
allelementC(33,:)='aS';
allelementC(34,:)='sE';
allelementC(35,:)='bR';
allelementC(36,:)='kR';
allelementC(37,:)='rB';
allelementC(38,:)='sR';
allelementC(39,:)=' Y';
allelementC(40,:)='zR';
allelementC(41,:)='nB';
allelementC(42,:)='mO';
allelementC(43,:)='tC';
allelementC(44,:)='rU';
allelementC(45,:)='rH';
allelementC(46,:)='pD';
allelementC(47,:)='aG';
allelementC(48,:)='cD';
allelementC(49,:)='iN';
allelementC(50,:)='sN';
allelementC(51,:)='sB';
allelementC(52,:)='tE';
allelementC(53,:)=' I';
allelementC(54,:)='xE';
allelementC(55,:)='cS';
allelementC(56,:)='bA';
allelementC(57,:)='lA';
allelementC(58,:)='cE';
allelementC(59,:)='pR';
allelementC(60,:)='nD';
allelementC(61,:)='pM';
allelementC(62,:)='sM';
allelementC(63,:)='eU';
allelementC(64,:)='gD';
allelementC(65,:)='tB';
allelementC(66,:)='dY';
allelementC(67,:)='hO';
allelementC(68,:)='eR';
allelementC(69,:)='tM';
allelementC(70,:)='yB';
allelementC(71,:)='lU';
allelementC(72,:)='hF';
allelementC(73,:)='tA';
allelementC(74,:)=' W';
allelementC(75,:)='rE';
allelementC(76,:)='oS';
allelementC(77,:)='iR';
allelementC(78,:)='pT';
allelementC(79,:)='aU';
allelementC(80,:)='hG';
allelementC(81,:)='tL';
allelementC(82,:)='pB';
allelementC(83,:)='bI';
allelementC(84,:)='pO';
allelementC(85,:)='aT';
allelementC(86,:)='rN';
allelementC(87,:)='fR';
allelementC(88,:)='rA';
allelementC(89,:)='aC';
allelementC(90,:)='tH';
allelementC(91,:)='pA';
allelementC(92,:)=' U';
allelementC(93,:)='nQ';
allelementC(94,:)='pU';
allelementC(95,:)='aM';
allelementC(96,:)='cM';
allelementC(97,:)='bK';
allelementC(98,:)='cF';
allelementC(99,:)='eS';
allelementC(100,:)='fM';
allelementC(101,:)='mD';
allelementC(102,:)='nO';
allelementC(103,:)='lR';
allelementC(104,:)='rF';
allelementC(105,:)='dB';
allelementC(106,:)='sG';
allelementC(107,:)='bH';
allelementC(108,:)='hS';
allelementC(109,:)='mT';
allelementC(110,:)='dS';
allelementC(111,:)='rG';


allelementD(1,:)=' h';
allelementD(2,:)='he';
allelementD(3,:)='li';
allelementD(4,:)='be';
allelementD(5,:)=' b';
allelementD(6,:)=' c';
allelementD(7,:)=' n';
allelementD(8,:)=' o';
allelementD(9,:)=' f';
allelementD(10,:)='ne';
allelementD(11,:)='na';
allelementD(12,:)='mg';
allelementD(13,:)='al';
allelementD(14,:)='si';
allelementD(15,:)=' p';
allelementD(16,:)=' s';
allelementD(17,:)='cl';
allelementD(18,:)='ar';
allelementD(19,:)=' k';
allelementD(20,:)='ca';
allelementD(21,:)='sc';
allelementD(22,:)='ti';
allelementD(23,:)=' v';
allelementD(24,:)='cr';
allelementD(25,:)='mn';
allelementD(26,:)='fe';
allelementD(27,:)='co';
allelementD(28,:)='ni';
allelementD(29,:)='cu';
allelementD(30,:)='zn';
allelementD(31,:)='ga';
allelementD(32,:)='ge';
allelementD(33,:)='as';
allelementD(34,:)='se';
allelementD(35,:)='br';
allelementD(36,:)='kr';
allelementD(37,:)='rb';
allelementD(38,:)='sr';
allelementD(39,:)=' y';
allelementD(40,:)='zr';
allelementD(41,:)='nb';
allelementD(42,:)='mo';
allelementD(43,:)='tc';
allelementD(44,:)='ru';
allelementD(45,:)='rh';
allelementD(46,:)='pd';
allelementD(47,:)='ag';
allelementD(48,:)='cd';
allelementD(49,:)='in';
allelementD(50,:)='sn';
allelementD(51,:)='sb';
allelementD(52,:)='te';
allelementD(53,:)=' i';
allelementD(54,:)='xe';
allelementD(55,:)='cs';
allelementD(56,:)='ba';
allelementD(57,:)='la';
allelementD(58,:)='ce';
allelementD(59,:)='pr';
allelementD(60,:)='nd';
allelementD(61,:)='pm';
allelementD(62,:)='sm';
allelementD(63,:)='eu';
allelementD(64,:)='gd';
allelementD(65,:)='tb';
allelementD(66,:)='dy';
allelementD(67,:)='ho';
allelementD(68,:)='er';
allelementD(69,:)='tm';
allelementD(70,:)='yb';
allelementD(71,:)='lu';
allelementD(72,:)='hf';
allelementD(73,:)='ta';
allelementD(74,:)=' w';
allelementD(75,:)='re';
allelementD(76,:)='os';
allelementD(77,:)='ir';
allelementD(78,:)='pt';
allelementD(79,:)='au';
allelementD(80,:)='hg';
allelementD(81,:)='tl';
allelementD(82,:)='pb';
allelementD(83,:)='bi';
allelementD(84,:)='po';
allelementD(85,:)='at';
allelementD(86,:)='rn';
allelementD(87,:)='fr';
allelementD(88,:)='ra';
allelementD(89,:)='ac';
allelementD(90,:)='th';
allelementD(91,:)='pa';
allelementD(92,:)=' u';
allelementD(93,:)='nq';
allelementD(94,:)='pu';
allelementD(95,:)='am';
allelementD(96,:)='cm';
allelementD(97,:)='bk';
allelementD(98,:)='cf';
allelementD(99,:)='es';
allelementD(100,:)='fm';
allelementD(101,:)='md';
allelementD(102,:)='no';
allelementD(103,:)='lr';
allelementD(104,:)='rf';
allelementD(105,:)='db';
allelementD(106,:)='sg';
allelementD(107,:)='bh';
allelementD(108,:)='hs';
allelementD(109,:)='mt';
allelementD(110,:)='ds';
allelementD(111,:)='rg';


for i=1:106;   %注意要3个字符，第一个标记这个格式是字符形式；H与前面的空格两位构成元素名称，因为有些元素名字是两个字母
%   if ~isempty(findstr(strr,strcat([' ',allelementA(i,:)])))  ||  ~isempty(findstr(strr,strcat([' ',allelementB(i,:)])))  ||  ~isempty(findstr(strr,strcat([' ',allelementC(i,:)])))   ||  ~isempty(findstr(strr,strcat([' ',allelementD(i,:)])))
%    if ~isempty(findstr(strr,allelementA(i,:)))  ||  ~isempty(findstr(strr,allelementB(i,:)))  ||  ~isempty(findstr(strr,allelementC(i,:)))   ||  ~isempty(findstr(strr,allelementD(i,:)))
%         ele=i; return;
% %    end
%     strr(find(isspace(strr)))=[];
    s1= allelementA(i,:);
    s1(find(isspace(s1)))=[];%去空格
    s2= allelementB(i,:);
    s2(find(isspace(s2)))=[];
    s3= allelementA(i,:);
    s3(find(isspace(s3)))=[];
    s4= allelementA(i,:);
    s4(find(isspace(s4)))=[];
    if strcmp(strr,s1) || strcmp(strr,s2) || strcmp(strr,s3) || strcmp(strr,s4)
        ele=i;return
    end
end


%如果无法辨认，输出
'Can NOT regonize ATOM ELEMENT'
ele=0;
return;

