function F=interaction1(x,W,k,a,Y)
n=length(Y);
A=zeros(n,n);
for i=1:n
    for j=1:n
    A(i,j)=(sin(Y(j)-Y(i))).*exp(-a.*(abs(i-j)));
    end
end
F=sum(A,2).*k+x*0+W;



    