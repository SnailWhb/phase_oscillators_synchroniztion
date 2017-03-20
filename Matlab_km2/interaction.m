function F=interaction(x,W,K,Y)
n=length(Y);
A=zeros(n,n);
for i=1:n
    for j=1:n
    A(i,j)=sin(Y(j)-Y(i));
    end
end
F=sum(A,2).*K+x*0+W;



    