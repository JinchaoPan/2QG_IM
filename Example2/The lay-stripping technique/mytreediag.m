function X = mytreediag(b,a,c,f)
%b主对角线，a下对角线，c上对角线；列向量
%要求abs（b）>abs(a)+abs(c),严格对角占优
n=length(b);
a=a(:);a=[0;a];
beta=zeros(n-1,1);
y=zeros(n,1);
X=zeros(n,1);
beta(1)=c(1)/b(1);
for i=2:1:n-1
    beta(i)=c(i)/(b(i)-a(i)*beta(i-1));
end
y(1)=f(1)/b(1);
for i=2:1:n
    y(i)=(f(i)-a(i)*y(i-1))/(b(i)-a(i)*beta(i-1));
end
X(n)=y(n);
for i=n-1:-1:1
    X(i)=y(i)-beta(i)*X(i+1);
end
end