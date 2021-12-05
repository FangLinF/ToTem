function myresul=mynewftds_ab(z,B,V);  %V的单位是kV
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
%h=4.13566743e-15; h_ba=h/2/pi;  %%%
e=1.602102e-19;
me=9.1093807e-31; % me=0.5110*10^6;
%The electron rest mass (symbol: m e) is the mass of a stationary electron, also known as the invariant mass of the electron. It is one of the fundamental constants of physics.It has a value of about 9.109 × 10 ?31 kilograms or about 5.486 × 10 ?4 atomic mass units, equivalent to an energy of about 8.187 × 10 ?14 joules or about 0.5110 MeV.
c=2.9979251e+8;
bili=sqrt(e*V*1000.*(1000*e.*V+2*me*c*c))./(me*c*c+1000*e.*V);
ele_v=bili*c;  %电子运动速度
para=2*h/(me*ele_v);

%a=cputime;
ab_re=ABSFpara;
ab=ab_re(z,:);

myfe6=myfe(ab,6);
alpha=sqrt(((0.023933754*z)/myfe6-36));  %由peng1996的公式16求的


[sx,sy]=meshgrid(-30:0.1:30, -30:0.1:30);
ss=sx+sqrt(-1)*sy;
s=0:0.1:30;
f=zeros(1,301);
for i=1:301
% %     f(i)=sum(fe(z,B,abs(s/2+s1(i))).*fe(z,B,abs(s/2-s1(i))).*(1-exp(-2*B*(s1(i).^2-s.^2/4))));  
     f(i)=0.1*0.1*sum(sum( myfe(ab,abs(s(i)/2+ss),alpha, z) .* myfe(ab,abs(s(i)/2-ss),alpha, z) .*(1-exp(-2*B*(ss.*conj(ss)-s(i).^2/4))).*exp(-B*s(i).^2/2)))  ;
end
%clear s s1 temp

f=f.*para*10^10;  %需要化简到埃的单位。
centf=f(1);  %0.01是积分里面ds的值
f=f./centf;
%b=cputime-a
% ft = fittype(@(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,x) a1*exp(-b1*x.^2)+ a2*exp(-b2*x.^2)+a3*exp(-b3*x.^2)+a4*exp(-b4*x.^2)+a5*exp(-b5*x.^2));
ft = fittype(@(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,x) a1*exp(-b1*x.^2)+ a2*exp(-b2*x.^2)+a3*exp(-b3*x.^2)+a4*exp(-b4*x.^2)+a5*exp(-b5*x.^2));
% 设置初始值
opts=fitoptions(ft);
%      F = FITOPTIONS('METHOD',VALUE)
%      NearestInterpolant     - Nearest neighbor interpolation
%      LinearInterpolant      - Linear interpolation
%      PchipInterpolant       - Piecewise Cubic Hermite interpolation.
%      CubicSplineInterpolant - Cubic Spline interpolation
%      BiharmonicInterpolant  - Biharmonic Surface Interpolation
%      SmoothingSpline        - Smoothing Spline 
%      LowessFit              - Lowess Fit
%      LinearLeastSquares     - Linear Least Squares
%      NonlinearLeastSquares  - Nonlinear Least Squares
opts.StartPoint=zeros(1,10);
opts.Lower=[-10,0,-10,0,0,0,0,0,-10,0];
% opts.Lower=[my_start(z,:)];
opts.Upper=[0,inf,10,inf,10,inf,10,inf,10,inf];
% startp=[my_start(z,:)];
% fitresult = lsqcurvefit(ft,startp,s,f);
f( find(isinf(f)) ) =0; f( find(isnan(f)) ) =0; 
fitresult =fit(s',f',ft,opts);
% ,'StartPoint',[my_start(z,:)]);
%%%%%%figure;plot(fitresult,s,f);  这条有用
% fitresult = fit(s',ftds',ft,opts);
% figure;plot(fitresult,s,ftds);
x=coeffvalues(fitresult);
myresul=zeros(1,10);
myresul(1)=x(1).*centf;  %a1
myresul(3)=x(2).*centf;
myresul(5)=x(3).*centf;
myresul(7)=x(4).*centf;
myresul(9)=x(5).*centf;  %a5
myresul(2)=x(6);         %b1
myresul(4)=x(7);
myresul(6)=x(8);
myresul(8)=x(9);
myresul(10)=x(10);       %b5

% figure;plot(s,f.*centf)
% hold on;
% ss=0:0.01:10;
% yy=myresul(1).*exp(-myresul(2)*ss.*ss)+myresul(3).*exp(-myresul(4)*ss.*ss)+myresul(5).*exp(-myresul(6)*ss.*ss)+...
% myresul(7).*exp(-myresul(8)*ss.*ss)+myresul(9).*exp(-myresul(10)*ss.*ss);
% hold on; plot(ss, yy,'g')
pp=1;


function fele=myfe(ab,s,alpha,z)
% e=find(abs(s)<=6);
% fele=zeros(size(s));
% fele(e)=ab(1)*exp(-ab(6)*s(e).^2)+ ab(2)*exp(-ab(7)*s(e).^2)+ab(3)*exp(-ab(8)*s(e).^2)+ab(4)*exp(-ab(9)*s(e).^2)+ab(5)*exp(-ab(10)*s(e).^2);


fele=zeros(size(s));
fele=ab(1)*exp(-ab(6)*s.^2)+ ab(2)*exp(-ab(7)*s.^2)+ab(3)*exp(-ab(8)*s.^2)+ab(4)*exp(-ab(9)*s.^2)+ab(5)*exp(-ab(10)*s.^2);
e=find(abs(s)>6);
if ~isempty(e)
    fele(e)=0.023933754*(z./(s(e).^2+alpha.^2));
end
% if s<=6.0
% % z为原子序数（1-98）
% % ab_re=ABSFpara;
% % ab=ab_re(z,:);
% 
% %ab是输入的10个系数
% % fe= ab(1)*exp(-(ab(6)+B/2)*s.^2)+ ab(2)*exp(-(ab(7)+B/2)*s.^2)+ab(3)*exp(-(ab(8)+B/2)*s.^2)+ab(4)*exp(-(ab(9)+B/2)*s.^2)+ab(5)*exp(-(ab(10)+B/2)*s.^2);
%     fele= ab(1)*exp(-ab(6)*s.^2)+ ab(2)*exp(-ab(7)*s.^2)+ab(3)*exp(-ab(8)*s.^2)+ab(4)*exp(-ab(9)*s.^2)+ab(5)*exp(-ab(10)*s.^2);
% else
%     fele=0.023933754*(z./(s.^2+alpha.^2));
% end