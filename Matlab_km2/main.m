clc;clear;
tspans=[0,4*pi];%区间端点;
h=0.1;%区间等分数目;
N=5000;%方程的个数;
% w=linspace(-2,2,N);
% W=(2./(1+w.^2))';%表示频率分布
w=linspace(-2,2,N);
W=zeros(N,1);
gama=1;
for i=1:N
    t=rand;
  W(i)=gama*tan(pi*(i-1/2)/t-1/2);
end
K=4;%表示coupling parameter
figure(K*1000+1);
Y0=rand(N,1);%初始值;
[x,Y]=Kuramoto_ODE45('interaction',tspans,h,W,K,Y0,N);
plot(x,Y);
colorbar;
xlabel('t');
ylabel('\theta');
title('The numerical solution of Kuramoto model');
%%%%%%%%%%%%%%%%%量化Kuramoto model 的耦合程度R(序参数).
figure(100*K+2)
theta=0:0.01:2*pi;
x0=cos(theta);
y0=sin(theta);
plot(x0,y0,'r-','LineWidth',1);
hold on
Theta=Y(:,end);
rend=sqrt((sum(cos(Theta))/N).^2+(sum(sin(Theta))/N).^2); %
scatter(cos(Theta),sin(Theta),'filled','cdata',w);%散点图。
title(['order parameter r=',num2str(rend)]);
colorbar
%%%%%%%%%%%%%%%%%参序数r与时间t的关系%%%%%%%%%%%%%
r=sqrt((sum(cos(Y))/N).^2+(sum(sin(Y))/N).^2);
figure(100*K+4);
plot(x,r);
colorbar;
xlabel('t');
ylabel('r');
title('参序数r与时间t的关系');
%%%%%%%%%%%%%%%当频率分布函数为g(w)=r/(pi*(r^2+w^2))时，r=(-1+sqrt(1+k^2))/k;
% figure(100*K+3)
% k=0:0.01:20;
% rr=(-1+sqrt(1+k.^2))./k;
% plot(k,rr,'r-');
% xlabel('coupling parameter k');
% ylabel('order parameter r');

% figure(2)
% clc;clear;
% W=[1,2]';po=[1,1]';c=1;
% dpt=@(P,W,c)W+c*sum(sin(meshgrid(P)-meshgrid(P)'),2);
% [t,p]=ode45(@(t,p)dpt(p,W,c),[0,pi],po);
% plot(t,p)

