function [x,Y]=Kuramoto_ODE451(interaction1,tspans,h,W,C,a,Y0,N)
%EuLer格式,
%求解方程y'=vectorfun(t,y);其中t \in[a,b];Y0为初始值；n为自变量的离散个数；Y为求解结果
%a1表示变量t前面的系数，A表示Y前面的系数。N表示向量分量的个数；
x=tspans(1):h:tspans(2);
n=length(x);
Y=zeros(N,n);%存放数值的解
x(1)=tspans(1);
Y(:,1)=Y0;
for i=1:n-1
K1=feval(interaction1,x(i),W,C,a,Y(:,i));
K2=feval(interaction1,x(i)+h/2,W,C,a,Y(:,i)+h/2.*K1);
K3=feval(interaction1,x(i)+h/2,W,C,a,Y(:,i)+h/2.*K2);%%Rugekutta四阶格式
K4=feval(interaction1,x(i)+h,W,C,a,Y(:,i)+h.*K3);
Y(:,i+1)=Y(:,i)+(h/6).*(K1+2.*K2+2.*K3+K4);
end





