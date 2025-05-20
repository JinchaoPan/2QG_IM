tic
Ro=Ro;% 考虑A，sigma不为0 Ro不为1
%% 网格
Lx=1; Ly=1; T=4*3;
m=50*1;n=50*1;k=400*1;
hx=Lx/m;hy=Ly/n;tau=T/k;
xi=0:hx:Lx; yj=0:hy:Ly; tk=0:tau:T;
[Xi,Yj]=meshgrid(xi,yj);
Xio=Xi(2:end-1,2:end-1);    Yjo=Yj(2:end-1,2:end-1);

DCy=diag(-ones(n+1,1))+diag(ones(n-1,1),2); DCy=DCy(1:end-2,:)/(2*hy);             % (n-1)*(n+1) x (n+1)*(m-1) = (n-1)*(m-1)
DCx=diag(-ones(m+1,1))+diag(ones(m-1,1),-2); DCx=DCx(:,1:end-2)/(2*hx);            % (n-1)*(m+1) x (m+1)*(m-1) = (n-1)*(m-1)

delta=0;
psi001=@(x,y,t)(psi1(x,y,t)+delta*((x.*y).*(x-Lx).*(y-Ly)));
psi002=@(x,y,t)(psi2(x,y,t)+delta*((x.*y).*(x-Lx).*(y-Ly)));

Atau_hx=tau/(hx)/hx/2/Re*sign(A); Atau_hy=tau/(hy)/hy/2/Re*sign(A);
%%
q1NUM=zeros((n-1)*(m-1),k);
q2NUM=zeros((n-1)*(m-1),k);
psi_NUM1=zeros((n-1)*(m-1),k);
psi_NUM2=zeros((n-1)*(m-1),k);
psi0_1=zeros(n+1,m+1);
psi0_2=zeros(n+1,m+1);

ta=tk(1);
q1_00=q1(Xio,Yjo,ta);   q1_00=q1_00(:);
q2_00=q2(Xio,Yjo,ta);   q2_00=q2_00(:);

q1NUM(:,1)=q1_00;
q2NUM(:,1)=q2_00;
%% 边界矩阵
bfive=Ro*1/(hx^2); cfive=Ro*1/(hy^2);  afive=-2*(bfive+cfive);
SX_PSI1=zeros(2,m+1,k+1); SY_PSI1=zeros(n+1,2,k+1);
SX_PSI2=zeros(2,m+1,k+1); SY_PSI2=zeros(n+1,2,k+1);
BBF=zeros(2*(n-1)*(m-1),k+1);
for kkf=1:1:k+1
    ta1=tk(kkf);
    Sx_psi1=[psi1(Xi(1,:),Yj(1,:),ta1);psi1(Xi(end,:),Yj(end,:),ta1)];
    Sy_psi1=[psi1(Xi(:,1),Yj(:,1),ta1) psi1(Xi(:,end),Yj(:,end),ta1)];

    Sx_psi2=[psi2(Xi(1,:),Yj(1,:),ta1);psi2(Xi(end,:),Yj(end,:),ta1)];
    Sy_psi2=[psi2(Xi(:,1),Yj(:,1),ta1) psi2(Xi(:,end),Yj(:,end),ta1)];

    SX_PSI1(:,:,kkf)=Sx_psi1; SY_PSI1(:,:,kkf)=Sy_psi1;
    SX_PSI2(:,:,kkf)=Sx_psi2; SY_PSI2(:,:,kkf)=Sy_psi2;

    Bxf1=[Sx_psi1(1,2:end-1);zeros((n-1)-2,m-1);Sx_psi1(end,2:end-1)];
    Bxf1=cfive*Bxf1(:);

    Bxf2=[Sx_psi2(1,2:end-1);zeros((n-1)-2,m-1);Sx_psi2(end,2:end-1)];
    Bxf2=cfive*Bxf2+Rs(Xio,Yjo,ta1); Bxf2=Bxf2(:);

    Byf1=bfive*[Sy_psi1(2:end-1,1);zeros((n-1)*(m-3),1);Sy_psi1(2:end-1,end)];
    Byf2=bfive*[Sy_psi2(2:end-1,1);zeros((n-1)*(m-3),1);Sy_psi2(2:end-1,end)];

    BBF(:,kkf)=[Bxf1;Bxf2]+[Byf1;Byf2]+betaK*[Yjo(:);Yjo(:)];
end

% Bf=[q1NUM(:,kkb);q2NUM(:,kkb)]-BBF(:,kkb);

%% 循环追赶法
Am1=afive*ones((n-1),(m-1))-F1;   Au1=cfive*ones((n-2),(m-1)); Ad1=cfive*ones((n-2),(m-1));
Auu1=bfive*ones((n-1),(m-1)-1);   Add1=bfive*ones((n-1),(m-1)-1);

Am2=afive*ones((n-1),(m-1))-F2;   Au2=cfive*ones((n-2),(m-1)); Ad2=cfive*ones((n-2),(m-1));
Auu2=bfive*ones((n-1),(m-1)-1);   Add2=bfive*ones((n-1),(m-1)-1);

FF1=F1;FF2=F2;
diagC=cfive*ones((n-2),(m-1)); diagC=[diagC;zeros(1,(m-1))]; diagC=diagC(:);

%% 追赶法7对角矩阵
% tzhui=tic;
% 五对角矩阵
% BETAA=(diag(Am1(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
%       +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1)))\diag(F1*ones((n-1)*(m-1),1));
Temp1=(diag(Am1(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
    +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1)));
BETAA=Temp1\diag(F1*ones((n-1)*(m-1),1));
Temp2=Temp1; Temp2(1:(n-1)*(m-1)+1:end)=afive-F2;

% YY1=BETAA*(Bf(1:(n-1)*(m-1),:)/F1);
%
% Temp2=Temp1; Temp2(1:(n-1)*(m-1)+1:end)=afive-F2;
%
% psi_NUM2=(Temp2-F2*BETAA)\(Bf((n-1)*(m-1)+1:end,:)-F2*YY1);
%
% psi_NUM1=YY1-BETAA*psi_NUM2;
% % toc(tzhui)
%
% %%
% psi_NUM=[psi_NUM1;psi_NUM2];

%%
SXQ_1k=zeros(2,(m+1),k+1);
SYQ_1k=zeros((n+1),2,k+1);
SXQ_2k=zeros(2,(m+1),k+1);
SYQ_2k=zeros((n+1),2,k+1);
FF1_0k=zeros((n-1)*(m-1),k+1);
FF2_0k=zeros((n-1)*(m-1),k+1);

YG1=zeros(k,1);
YG2=zeros(k,1);
for kk=1:1:k
    Bf=[q1NUM(:,kk);q2NUM(:,kk)]-BBF(:,kk);
    YY1=BETAA*(Bf(1:(n-1)*(m-1))/F1);


    psi_NUM2k=(Temp2-F2*BETAA)\(Bf((n-1)*(m-1)+1:end)-F2*YY1);
    psi_NUM1k=YY1-BETAA*psi_NUM2k;

    psi_NUM1(:,kk)=psi_NUM1k;
    psi_NUM2(:,kk)=psi_NUM2k;


    %%
    psi0_1(2:end-1,2:end-1)=reshape(psi_NUM1k,[n-1,m-1]);
    psi0_2(2:end-1,2:end-1)=reshape(psi_NUM2k,[n-1,m-1]);

    Sx_psi1=SX_PSI1(:,:,kk);
    Sy_psi1=SY_PSI1(:,:,kk);
    Sx_psi2=SX_PSI2(:,:,kk);
    Sy_psi2=SY_PSI2(:,:,kk);

    psi0_1([1 n+1],:)=Sx_psi1;  psi0_1(:,[1 m+1])=Sy_psi1;
    psi0_2([1 n+1],:)=Sx_psi2;  psi0_2(:,[1 m+1])=Sy_psi2;

    dxpsi1=psi0_1(2:end-1,:)*DCx;
    dypsi1=DCy*psi0_1(:,2:end-1);

    dxpsi2=psi0_2(2:end-1,:)*DCx;
    dypsi2=DCy*psi0_2(:,2:end-1);

    %%
    ta1=tk(kk+1);
    q1_0=q1NUM(:,kk); q2_0=q2NUM(:,kk);
    f1_0=f1(Xio,Yjo,ta1);  f1_0=f1_0(:); f2_0=f2(Xio,Yjo,ta1);  f2_0=f2_0(:);

    Sx_q1=[q1(Xi(1,:),Yj(1,:),ta1);q1(Xi(end,:),Yj(end,:),ta1)];
    Sy_q1=[q1(Xi(:,1),Yj(:,1),ta1) q1(Xi(:,end),Yj(:,end),ta1)];

    Sx_q2=[q2(Xi(1,:),Yj(1,:),ta1);q2(Xi(end,:),Yj(end,:),ta1)];
    Sy_q2=[q2(Xi(:,1),Yj(:,1),ta1) q2(Xi(:,end),Yj(:,end),ta1)];

    SXQ_1k(:,:,kk+1)=Sx_q1;  SYQ_1k(:,:,kk+1)=Sy_q1;
    SXQ_2k(:,:,kk+1)=Sx_q2;  SYQ_2k(:,:,kk+1)=Sy_q2;
    FF1_0k(:,kk+1)=f1_0;     FF2_0k(:,kk+1)=f2_0;

    %%
    DDpsi1d2k=q1NUM(:,kk)-q2NUM(:,kk)+(F1+F2)*(psi_NUM1(:,kk)-psi_NUM2(:,kk));
    DDpsi1d2k=DDpsi1d2k/Ro;

    DDpsi2k=q2NUM(:,kk)+F2*(psi_NUM2(:,kk)-psi_NUM1(:,kk))-Yjo(:);
    DDpsi2k=DDpsi2k/Ro;
%     DDpsi2k=DD_psi2(Xio,Yjo,tk(kk));  DDpsi2k= DDpsi2k(:);
    f1_0=f1_0+F1*DDpsi1d2k/Re*sign(A);
    f2_0=f2_0-Sigma*DDpsi2k-F2*DDpsi1d2k/Re*sign(A);
    %% psi偏微分
    % dxpsi1=reshape(dxpsi0_1k(:,kk),[n-1,m-1]);
    % dxpsi2=reshape(dxpsi0_2k(:,kk),[n-1,m-1]);
    % dypsi1=reshape(dypsi0_1k(:,kk),[n-1,m-1]);
    % dypsi2=reshape(dypsi0_2k(:,kk),[n-1,m-1]);

    % 考虑变号
    dxpsi1=((-1)^(Jk))*dxpsi1;    dypsi1=((-1)^(Jk))*dypsi1;
    % dx0psi1=((-1)^(Jk))*dx0psi1;  dxmpsi1=((-1)^(Jk))*dxmpsi1;

    dxpsi2=((-1)^(Jk))*dxpsi2;    dypsi2=((-1)^(Jk))*dypsi2;
    % dx0psi2=((-1)^(Jk))*dx0psi2;  dxmpsi2=((-1)^(Jk))*dxmpsi2;


    dxpsi1=dxpsi1*(-tau/(2*hy));   dypsi1=dypsi1*(tau/(2*hx));
    dxpsi2=dxpsi2*(-tau/(2*hy));   dypsi2=dypsi2*(tau/(2*hx));
    %% 已知量=边界+右端
    BS=[dxpsi1(1,:)-Atau_hy;-dxpsi1(end,:)-Atau_hy].*Sx_q1(:,2:end-1);
    BS=[BS(1,:);zeros((n-1)-2,m-1);BS(end,:)]; BS=BS(:);

    BS0=[dxpsi2(1,:)-Atau_hy;-dxpsi2(end,:)-Atau_hy].*Sx_q2(:,2:end-1);
    BS0=[BS0(1,:);zeros((n-1)-2,m-1);BS0(end,:)]; BS0=BS0(:);

    diagr1=dypsi1;
    BS1=[(diagr1(:,1)-Atau_hx).*Sy_q1(2:end-1,1);zeros((n-1)*(m-3),1);(-diagr1(:,end)-Atau_hx).*Sy_q1(2:end-1,end)];
    diagr10=dypsi2;
    BS10=[(diagr10(:,1)-Atau_hx).*Sy_q2(2:end-1,1);zeros((n-1)*(m-3),1);(-diagr10(:,end)-Atau_hx).*Sy_q2(2:end-1,end)];

    B=q1_0-BS-BS1+tau*f1_0;

    B0=q2_0-BS0-BS10+tau*f2_0;

    %% 追赶法求解
    Gm=ones((n-1),(m-1))+2*(Atau_hx+Atau_hy); Gd=dxpsi1(2:end,:)-Atau_hy;Gu=-dxpsi1(1:end-1,:)-Atau_hy;
    Gdd=dypsi1(:,2:end)-Atau_hx; Guu=-dypsi1(:,1:end-1)-Atau_hx;
    q1_num=Mykuaitreediag(Gm,Gd,Gu,Gdd,Guu,reshape(B,[n-1,m-1]));
    q1_num=q1_num(:);

    Gm0=ones((n-1),(m-1))+2*(Atau_hx+Atau_hy); Gd0=dxpsi2(2:end,:)-Atau_hy;Gu0=-dxpsi2(1:end-1,:)-Atau_hy;
    Gdd0=dypsi2(:,2:end)-Atau_hx; Guu0=-dypsi2(:,1:end-1)-Atau_hx;
    q2_num=Mykuaitreediag(Gm0,Gd0,Gu0,Gdd0,Guu0,reshape(B0,[n-1,m-1]));
    q2_num=q2_num(:);

    % YG1(kk)=max(sum(abs(G),2));
    % YG2(kk)=max(sum(abs(G0),2));

    %%
    q1NUM(:,kk+1)=q1_num;
    q2NUM(:,kk+1)=q2_num;
end

%%
Bf=[q1NUM(:,k+1);q2NUM(:,k+1)]-BBF(:,k+1);
YY1=BETAA*(Bf(1:(n-1)*(m-1))/F1);

Temp2=Temp1; Temp2(1:(n-1)*(m-1)+1:end)=afive-F2;

psi_NUM2k=(Temp2-F2*BETAA)\(Bf((n-1)*(m-1)+1:end)-F2*YY1);
psi_NUM1k=YY1-BETAA*psi_NUM2k;

psi_NUM1(:,k+1)=psi_NUM1k;
psi_NUM2(:,k+1)=psi_NUM2k;
%%
tk_ms=k+1;
% ta1=tk(tk_ms);
q1_k=q1(Xio,Yjo,tk(tk_ms));  q1_k=q1_k(:);
q2_k=q2(Xio,Yjo,tk(tk_ms));  q2_k=q2_k(:);

abs_error1=abs(q1_num-q1_k);
abs_error2=abs(q2_num-q2_k);

q1_num1=reshape(q1_num,[n-1,m-1]);
q1_k1=reshape(q1_k,[n-1,m-1]);
abs_error1=reshape(abs_error1,[n-1,m-1]);
q1_01=reshape(q1_00,[n-1,m-1]);

q2_num1=reshape(q2_num,[n-1,m-1]);
q2_k1=reshape(q2_k,[n-1,m-1]);
abs_error2=reshape(abs_error2,[n-1,m-1]);
q2_01=reshape(q2_00,[n-1,m-1]);

aaeq1=max(abs_error1(:));
aaeq2=max(abs_error2(:));

figure(1)
subplot(2,2,1)
surf(Xio,Yjo,q1_k1,'EdgeColor','none')
colorbar
title('True: q1')
view(0,90)
subplot(2,2,2)
surf(Xio,Yjo,q1_num1,'EdgeColor','none')
colorbar
title('Num: q1')
view(0,90)
subplot(2,2,4)
surf(Xio,Yjo,abs_error1,'EdgeColor','none')
colorbar
title('Absolution error of q1')
view(0,90)
subplot(2,2,3)
surf(Xio,Yjo,q1_01,'EdgeColor','none')
colorbar
title('Initial value: q1')
view(0,90)

%%
figure(2)
subplot(2,2,1)
surf(Xio,Yjo,q2_k1,'EdgeColor','none')
colorbar
title('True: q2')
view(0,90)
subplot(2,2,2)
surf(Xio,Yjo,q2_num1,'EdgeColor','none')
colorbar
title('Num: q2')
view(0,90)
subplot(2,2,4)
surf(Xio,Yjo,abs_error2,'EdgeColor','none')
colorbar
title('Absolution error of q2')
view(0,90)
subplot(2,2,3)
surf(Xio,Yjo,q2_01,'EdgeColor','none')
colorbar
title('Initial value: q2')
view(0,90)

%%
psi_NUM=[psi_NUM1;psi_NUM2];
psi_num=psi_NUM(:,tk_ms);
ta1=tk(tk_ms);
psi1_0=psi1(Xio,Yjo,ta1);  psi1_0=psi1_0(:);  cb1=[min(psi1_0(:))-1e-10,max(psi1_0(:))];
psi2_0=psi2(Xio,Yjo,ta1);  psi2_0=psi2_0(:);  cb2=[min(psi2_0(:))-1e-10,max(psi2_0(:))];
psi_0=[psi1_0;psi2_0];

error_absf=abs(psi_num-psi_0);

psi_01=reshape(psi_0(1:(n-1)*(m-1)),[n-1,m-1]);
psi_02=reshape(psi_0((n-1)*(m-1)+1:end),[n-1,m-1]);

psi_num1=reshape(psi_num(1:(n-1)*(m-1)),[n-1,m-1]);
psi_num2=reshape(psi_num((n-1)*(m-1)+1:end),[n-1,m-1]);

error_absf1=reshape(error_absf(1:(n-1)*(m-1)),[n-1,m-1]);
error_absf2=reshape(error_absf((n-1)*(m-1)+1:end),[n-1,m-1]);

psi1_001=psi001(Xi,Yj,ta1);
psi2_001=psi002(Xi,Yj,ta1);


figure(3)
subplot(2,2,1)
surf(Xio,Yjo,psi_01,'EdgeColor','none')
colorbar
caxis(cb1)
title('True: \psi_1')
view(0,90)
subplot(2,2,2)
surf(Xio,Yjo,psi_num1,'EdgeColor','none')
colorbar
caxis(cb1)
title('Num: \psi_1')
view(0,90)
subplot(2,2,3)
surf(Xio,Yjo,error_absf1,'EdgeColor','none')
colorbar
title('Absolution error of \psi_1')
view(0,90)
subplot(2,2,4)
surf(Xi,Yj,psi1_001,'EdgeColor','none')
colorbar
title('Initial value of \psi_1')
view(0,90)

figure(4)
subplot(2,2,1)
surf(Xio,Yjo,psi_02,'EdgeColor','none')
colorbar
title('True: \psi_2')
view(0,90)
caxis(cb2)
subplot(2,2,2)
surf(Xio,Yjo,psi_num2,'EdgeColor','none')
colorbar
title('Num: \psi_2')
view(0,90)
caxis(cb2)
subplot(2,2,3)
surf(Xio,Yjo,error_absf2,'EdgeColor','none')
colorbar
title('Absolution error of \psi_2')
view(0,90)
subplot(2,2,4)
surf(Xi,Yj,psi2_001,'EdgeColor','none')
colorbar
title('Initial value of \psi_2')
view(0,90)
% clear G G0
toc