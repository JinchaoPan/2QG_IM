%%
Ro=Ro;
Atau_hx=tau/(hx)/hx/2/Re*sign(A); Atau_hy=tau/(hy)/hy/2/Re*sign(A);
%%
q1NUM=zeros((n-1)*(m-1),k+1);
q2NUM=zeros((n-1)*(m-1),k+1);
% q1_00=q1(Xio,Yjo,tk(1));   q1_00=q1_00(:); q1NUM(:,1)=q1_00;
% q2_00=q2(Xio,Yjo,tk(1));   q2_00=q2_00(:); q2NUM(:,1)=q2_00;
q1_00=q1_00T;   q1_00=q1_00(:); q1NUM(:,1)=q1_00;
q2_00=q2_00T;   q2_00=q2_00(:); q2NUM(:,1)=q2_00;
%%
SXQ_1k=zeros(2,(m+1),k+1);
SYQ_1k=zeros((n+1),2,k+1);
SXQ_2k=zeros(2,(m+1),k+1);
SYQ_2k=zeros((n+1),2,k+1);
FF1_0k=zeros((n-1)*(m-1),k+1);
FF2_0k=zeros((n-1)*(m-1),k+1);
tic
%%
for kk=1:1:k
q1_0=q1NUM(:,kk);
q2_0=q2NUM(:,kk);

ta=tk(kk); ta1=tk(kk+1);
f1_0=f1(Xio,Yjo,ta+tau/2)';  f1_0=f1_0(:);   % x排列  
% f2_0=f2(Xio,Yjo,ta+tau/2)'-Sigma*DD_psi2(Xio,Yjo,ta+tau/2)';  f2_0=f2_0(:);   % x排列  
% f2_0=f2(Xio,Yjo,ta+tau/2)';  f2_0=f2_0(:)-Sigma*DDpsi2k(:,kk);   % x排列  
f2_0=f2(Xio,Yjo,ta+tau/2)';  f2_0=f2_0(:);   % x排列  

Sx_q1=[q1(Xi(1,:),Yj(1,:),ta);q1(Xi(end,:),Yj(end,:),ta)];
Sy_q1=[q1(Xi(:,1),Yj(:,1),ta) q1(Xi(:,end),Yj(:,end),ta)];
Sx_q11=[q1(Xi(1,:),Yj(1,:),ta1);q1(Xi(end,:),Yj(end,:),ta1)];
Sy_q11=[q1(Xi(:,1),Yj(:,1),ta1) q1(Xi(:,end),Yj(:,end),ta1)];

Sx_q2=[q2(Xi(1,:),Yj(1,:),ta);q2(Xi(end,:),Yj(end,:),ta)]; 
Sy_q2=[q2(Xi(:,1),Yj(:,1),ta) q2(Xi(:,end),Yj(:,end),ta)];
Sx_q22=[q2(Xi(1,:),Yj(1,:),ta1);q2(Xi(end,:),Yj(end,:),ta1)];
Sy_q22=[q2(Xi(:,1),Yj(:,1),ta1) q2(Xi(:,end),Yj(:,end),ta1)];

SXQ_1k(:,:,kk)=Sx_q1;  SYQ_1k(:,:,kk)=Sy_q1;
SXQ_2k(:,:,kk)=Sx_q2;  SYQ_2k(:,:,kk)=Sy_q2;
SXQ_1k(:,:,kk+1)=Sx_q11;  SYQ_1k(:,:,kk+1)=Sy_q11;
SXQ_2k(:,:,kk+1)=Sx_q22;  SYQ_2k(:,:,kk+1)=Sy_q22;
FF1_0k(:,kk)=f1_0;     FF2_0k(:,kk)=f2_0;
%% psi偏微分
dxpsi1=reshape(dxpsi0_1k(:,kk),[n-1,m-1]); 
dxpsi2=reshape(dxpsi0_2k(:,kk),[n-1,m-1]);
dypsi1=reshape(dypsi0_1k(:,kk),[n-1,m-1]); 
dypsi2=reshape(dypsi0_2k(:,kk),[n-1,m-1]);

dx0psi1=dxpsi0m_1k(:,1,kk);
dxmpsi1=dxpsi0m_1k(:,end,kk);

dx0psi2=dxpsi0m_2k(:,1,kk);
dxmpsi2=dxpsi0m_2k(:,end,kk);
%% psi在x=0,x=m处偏微分
dxpsi1=dxpsi1*(-tau/(4*hy)); dypsi1=dypsi1*(tau/(4*hx)); 
dx0psi1=dx0psi1*(-tau/(4*hy));  dxmpsi1=dxmpsi1*(-tau/(4*hy));

dxpsi2=dxpsi2*(-tau/(4*hy)); dypsi2=dypsi2*(tau/(4*hx)); 
dx0psi2=dx0psi2*(-tau/(4*hy));  dxmpsi2=dxmpsi2*(-tau/(4*hy));

% 考虑变号
dxpsi1=((-1)^(Jk))*dxpsi1;    dypsi1=((-1)^(Jk))*dypsi1; 
dx0psi1=((-1)^(Jk))*dx0psi1;  dxmpsi1=((-1)^(Jk))*dxmpsi1;

dxpsi2=((-1)^(Jk))*dxpsi2;    dypsi2=((-1)^(Jk))*dypsi2; 
dx0psi2=((-1)^(Jk))*dx0psi2;  dxmpsi2=((-1)^(Jk))*dxmpsi2;
%% 右端qq**边界值
Ddxpsi1=[dx0psi1 dxpsi1 dxmpsi1]; % n-1 x m+1
Ddxpsi2=[dx0psi2 dxpsi2 dxmpsi2]; % n-1 x m+1
q1_0=reshape(q1_0,[n-1,m-1]);  
q2_0=reshape(q2_0,[n-1,m-1]); 
q1_0=[Sx_q1(1,:);[Sy_q1(2:end-1,1) q1_0 Sy_q1(2:end-1,end)];Sx_q1(end,:)]; %(n+1) x (m+1)
q2_0=[Sx_q2(1,:);[Sy_q2(2:end-1,1) q2_0 Sy_q2(2:end-1,end)];Sx_q2(end,:)]; %(n+1) x (m+1)

qq10=zeros(n-1,m+1);
qq20=zeros(n-1,m+1);
for i=1:1:m+1
    for j=1:1:n-1
        qq10(j,i)=(-Ddxpsi1(j,i)+Atau_hy)*q1_0((j+1)-1,i)+(1-2*Atau_hy)*q1_0((j+1),i)+(Ddxpsi1(j,i)+Atau_hy)*q1_0((j+1)+1,i);     
        qq20(j,i)=(-Ddxpsi2(j,i)+Atau_hy)*q2_0((j+1)-1,i)+(1-2*Atau_hy)*q2_0((j+1),i)+(Ddxpsi2(j,i)+Atau_hy)*q2_0((j+1)+1,i);
    end
end

BQQL1=zeros(m-1,n-1);
BQQL2=zeros(m-1,n-1);
for i=1:1:m-1
    for j=1:1:n-1
        BQQL1(i,j)=(-dypsi1(j,i)+Atau_hx)*qq10(j,(i+1)-1)+(1-2*Atau_hx)*qq10(j,(i+1))+(dypsi1(j,i)+Atau_hx)*qq10(j,(i+1)+1);
        BQQL2(i,j)=(-dypsi2(j,i)+Atau_hx)*qq20(j,(i+1)-1)+(1-2*Atau_hx)*qq20(j,(i+1))+(dypsi2(j,i)+Atau_hx)*qq20(j,(i+1)+1);
    end
end

%% q*初值与边界
QQx01=zeros(n-1,1);  QQxm1=zeros(n-1,1);
Sx0Q1=Sy_q11(:,1);   SxmQ1=Sy_q11(:,end);                                  % x=0,x=m边界处q^(k+1)值 (n+1) x 1

QQx02=zeros(n-1,1);  QQxm2=zeros(n-1,1);
Sx0Q2=Sy_q22(:,1);   SxmQ2=Sy_q22(:,end);                                  % x=0,x=m边界处q^(k+1)值 (n+1) x 1
for j=1:1:n-1
    QQx01(j)=(dx0psi1(j)-Atau_hy)*Sx0Q1((j+1)-1)+(1+2*Atau_hy)*Sx0Q1((j+1))-(dx0psi1(j)+Atau_hy)*Sx0Q1((j+1)+1);
    QQxm1(j)=(dxmpsi1(j)-Atau_hy)*SxmQ1((j+1)-1)+(1+2*Atau_hy)*SxmQ1((j+1))-(dxmpsi1(j)+Atau_hy)*SxmQ1((j+1)+1);

    QQx02(j)=(dx0psi2(j)-Atau_hy)*Sx0Q2((j+1)-1)+(1+2*Atau_hy)*Sx0Q2((j+1))-(dx0psi2(j)+Atau_hy)*Sx0Q2((j+1)+1);
    QQxm2(j)=(dxmpsi2(j)-Atau_hy)*SxmQ2((j+1)-1)+(1+2*Atau_hy)*SxmQ2((j+1))-(dxmpsi2(j)+Atau_hy)*SxmQ2((j+1)+1);
end

%% 边界构造
% BQQS=[dypsi1(:,1).*QQx01 zeros(n-1,m-3) -dypsi1(:,end).*QQxm1]';
% f1_0=reshape(f1_0,[m-1,n-1]);
% BQQ=BQQL-BQQS+tau*f1_0;

BQQS1=[(dypsi1(:,1)-Atau_hx).*QQx01 (-dypsi1(:,end)-Atau_hx).*QQxm1]';
BQQL1([1 end],:)=BQQL1([1 end],:)-BQQS1;

f1_0=f1_0+F1*DDpsi1d2k(:,kk)/Re;
f1_0=reshape(f1_0,[m-1,n-1]);
BQQ1=BQQL1+tau*f1_0;

BQQS2=[(dypsi2(:,1)-Atau_hx).*QQx02 (-dypsi2(:,end)-Atau_hx).*QQxm2]';
BQQL2([1 end],:)=BQQL2([1 end],:)-BQQS2;

f2_0=f2_0-Sigma*DDpsi2k(:,kk)-F2*DDpsi1d2k(:,kk)/Re-DD_RS/Re;
f2_0=reshape(f2_0,[m-1,n-1]);
BQQ2=BQQL2+tau*f2_0;
%% q*初值+边界构造
Qstart1=zeros(m-1,n-1);
Qstart2=zeros(m-1,n-1);

for i=1:1:n-1
    Qstart1(:,i)=mytreediag(ones(m-1,1)+2*Atau_hx,dypsi1(i,2:end)-Atau_hx,-dypsi1(i,1:end-1)-Atau_hx,BQQ1(:,i)); 
    Qstart2(:,i)=mytreediag(ones(m-1,1)+2*Atau_hx,dypsi2(i,2:end)-Atau_hx,-dypsi2(i,1:end-1)-Atau_hx,BQQ2(:,i));
end
Qstart1=Qstart1'; %(n-1) x (m-1)
Qstart2=Qstart2'; %(n-1) x (m-1)
%% Sx_q11
Sy0Q1=Sx_q11(1,2:end-1); SynQ1=Sx_q11(end,2:end-1);
BFQ1=[(dxpsi1(1,:)-Atau_hy).*Sy0Q1; (-dxpsi1(end,:)-Atau_hy).*SynQ1];
Qstart1([1 end],:)=Qstart1([1 end],:)-BFQ1;

Sy0Q2=Sx_q22(1,2:end-1); SynQ2=Sx_q22(end,2:end-1);
BFQ2=[(dxpsi2(1,:)-Atau_hy).*Sy0Q2; (-dxpsi2(end,:)-Atau_hy).*SynQ2];
Qstart2([1 end],:)=Qstart2([1 end],:)-BFQ2;

q1_num=zeros(n-1,m-1);
q2_num=zeros(n-1,m-1);
for i=1:1:m-1
    q1_num(:,i)=mytreediag(ones(n-1,1)+2*Atau_hy,dxpsi1(2:end,i)-Atau_hy,-dxpsi1(1:end-1,i)-Atau_hy,Qstart1(:,i)); 
    q2_num(:,i)=mytreediag(ones(n-1,1)+2*Atau_hy,dxpsi2(2:end,i)-Atau_hy,-dxpsi2(1:end-1,i)-Atau_hy,Qstart2(:,i)); 
end
q1_num=q1_num(:);
q2_num=q2_num(:);
% 滤波
% q1_num=MyQG_Filter_Linear(q1_num,Sx_q1(:,2:end-1),Sy_q1(2:end-1,:),alpha1,hx,hy,FG);
% q2_num=MyQG_Filter_Linear(q2_num,Sx_q2(:,2:end-1),Sy_q2(2:end-1,:),alpha1,hx,hy,FG);
%
q1NUM(:,kk+1)=q1_num;
q2NUM(:,kk+1)=q2_num;
end
%% t=T/2
% cb1=[-4.61357302560226,-1.60237378196137];
% cb2=[-4.78528749156675	-0.680360329543651];
tk_ms=k/2+1;
ta1=tk(tk_ms);
q1_k1=q1(Xi,Yj,tk(tk_ms));  
q2_k1=q2(Xi,Yj,tk(tk_ms)); 
tk_ms2=k+1;
q1_k2=q1(Xi,Yj,tk(tk_ms2));  
q2_k2=q2(Xi,Yj,tk(tk_ms2)); 
cbq1=[min(min(q1_k1(:)),min(q1_k2(:))),max(max(q1_k1(:)),max(q1_k2(:)))];
cbq2=[min(min(q2_k1(:)),min(q2_k2(:))),max(max(q2_k1(:)),max(q2_k2(:)))];

q1_num=q1NUM(:,tk_ms);     q1_num1=reshape(q1_num,[n-1,m-1]);
q2_num=q2NUM(:,tk_ms);     q2_num1=reshape(q2_num,[n-1,m-1]);
q1_num1=[SYQ_1k(:,1,tk_ms) [SXQ_1k(1,2:end-1,tk_ms);q1_num1;SXQ_1k(end,2:end-1,tk_ms)] SYQ_1k(:,end,tk_ms)];
q2_num1=[SYQ_2k(:,1,tk_ms) [SXQ_2k(1,2:end-1,tk_ms);q2_num1;SXQ_2k(end,2:end-1,tk_ms)] SYQ_2k(:,end,tk_ms)];
abs_error1=abs(q1_num1-q1_k1);
abs_error2=abs(q2_num1-q2_k1);


figure(1)
subplot(2,2,1)
surf(Xi,Yj,q1_k1,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Exact solution of q1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q1_num1,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Numerical result of q1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error1,'EdgeColor','none')
colorbar
title(['Absolution error of q1: t=' num2str(ta1)])
view(0,90)
% q1_01=reshape(q1_00,[n-1,m-1]);
% subplot(2,2,4)
% surf(Xio,Yjo,q1_01,'EdgeColor','none')
% colorbar
% title(['Initial value: q1: t=' num2str(ta1)])
% view(0,90)

figure(2)
subplot(2,2,1)
surf(Xi,Yj,q2_k1,'EdgeColor','none')
colorbar
caxis(cbq2)
% caxis(cb)
title(['Exact solution of q2: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q2_num1,'EdgeColor','none')
colorbar
% caxis(cb)
title(['Numerical result of q2: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error2,'EdgeColor','none')
colorbar
title(['Absolution error of q2: t=' num2str(ta1)])
view(0,90)
% q2_01=reshape(q2_00,[n-1,m-1]);
% subplot(2,2,4)
% surf(Xio,Yjo,q2_01,'EdgeColor','none')
% colorbar
% % caxis(cb)
% title(['Initial value: q2: t=' num2str(ta1)])
% view(0,90)

%% t=T
tk_ms=k+1;
ta1=tk(tk_ms);
q1_k1=q1(Xi,Yj,tk(tk_ms));  
q2_k1=q2(Xi,Yj,tk(tk_ms)); 
q1_num=q1NUM(:,tk_ms);     q1_num1=reshape(q1_num,[n-1,m-1]);
q2_num=q2NUM(:,tk_ms);     q2_num1=reshape(q2_num,[n-1,m-1]);
q1_num1=[SYQ_1k(:,1,tk_ms) [SXQ_1k(1,2:end-1,tk_ms);q1_num1;SXQ_1k(end,2:end-1,tk_ms)] SYQ_1k(:,end,tk_ms)];
q2_num1=[SYQ_2k(:,1,tk_ms) [SXQ_2k(1,2:end-1,tk_ms);q2_num1;SXQ_2k(end,2:end-1,tk_ms)] SYQ_2k(:,end,tk_ms)];
abs_error1=abs(q1_num1-q1_k1);
abs_error2=abs(q2_num1-q2_k1);


figure(3)
subplot(2,2,1)
surf(Xi,Yj,q1_k1,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Exact solution of q1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q1_num1,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Numerical result of q1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error1,'EdgeColor','none')
colorbar
title(['Absolution error of q1: t=' num2str(ta1)])
view(0,90)
% q1_01=reshape(q1_00,[n-1,m-1]);
% subplot(2,2,4)
% surf(Xio,Yjo,q1_01,'EdgeColor','none')
% colorbar
% title(['Initial value: q1: t=' num2str(ta1)])
% view(0,90)
%%
figure(4)
subplot(2,2,1)
surf(Xi,Yj,q2_k1,'EdgeColor','none')
colorbar
caxis(cbq2)
title(['Exact solution of q2: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q2_num1,'EdgeColor','none')
colorbar
caxis(cbq2)
title(['Numerical result of q2: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error2,'EdgeColor','none')
colorbar
title(['Absolution error of q2: t=' num2str(ta1)])
view(0,90)
% q2_01=reshape(q2_00,[n-1,m-1]);
% subplot(2,2,4)
% surf(Xio,Yjo,q2_01,'EdgeColor','none')
% colorbar
% % caxis(cb)
% title(['Initial value: q2: t=' num2str(ta1)])
% view(0,90)
toc
