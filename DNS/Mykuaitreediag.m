function X = Mykuaitreediag(Bm,Bd,Bu,A,C,F)
% Bm为B中主对角元素矩阵排列Bm:n x m; n为每个块三角矩阵的维度，m为整个块三角矩阵维度
% Bu为B中上对角元素矩阵排列Bu:(n-1) x (m); 砍去最后一行
% Bd为B中下对角元素矩阵排列Bd:(n-1) x (m); 砍去第一一行
% A为块三对角矩阵中下块对角线矩阵A:(n)x(m-1);为对角矩阵; 砍去第一一列
% C为块三对角矩阵中上块对角线矩阵C:(n)x(m-1);为对角矩阵; 砍去最后一列
% F为右端F:n x m
% B主对角线，A下对角线，C上对角线
%要求abs（b）>abs(a)+abs(c),严格对角占优
[n,m]=size(Bm);
Beta=zeros(n,n,m-1);
y=zeros(n,m);
X=zeros(n,m);
A=[zeros(n,1) A];
% beta(1)=c(1)/b(1);
% Beta(:,:,1)=(diag(Bm(:,1))+diag(Bd(:,1),-1)+diag(Bu(:,1),1))\diag(C(:,1));
for i=1:1:n  % 用追赶法，所用时间差不多
    CCi=zeros(n,1);CCi(i)=C(i,1);
    Beta(:,i,1)=mytreediag(Bm(:,1),Bd(:,1),Bu(:,1),CCi);
end

for i=2:1:m-1
%     Beta(:,i)=mytreediag(Bm(:,i)-A(:,i).*Beta(:,i-1),Bd(:,i),Bu(:,i),C(:,i));
    Beta(:,:,i)=(diag(Bm(:,i))+diag(Bd(:,i),-1)+diag(Bu(:,i),1)-A(:,i).*Beta(:,:,i-1))\diag(C(:,i));
end
% y(1)=f(1)/b(1);
y(:,1)=mytreediag(Bm(:,1),Bd(:,1),Bu(:,1),F(:,1));
for i=2:1:m
%     y(i)=(f(i)-a(i)*y(i-1))/(b(i)-a(i)*beta(i-1));
%     y(:,i)=mytreediag(Bm(:,i)-A(:,i).*Beta(:,i-1),Bd(:,i),Bu(:,i),F(:,i)-A(:,i).*y(:,i-1));
    y(:,i)=(diag(Bm(:,i))+diag(Bd(:,i),-1)+diag(Bu(:,i),1)-A(:,i).*Beta(:,:,i-1))\(F(:,i)-A(:,i).*y(:,i-1));

end
X(:,m)=y(:,m);
for i=m-1:-1:1
    X(:,i)=y(:,i)-Beta(:,:,i)*X(:,i+1);
end
end