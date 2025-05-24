%% 已知量处理
tic
%% 边界矩阵
% RS=Rs(Xio,Yjo,ta1);
% bfive=1/(hx^2); cfive=1/(hy^2);  afive=-2*(bfive+cfive);
% SX_PSI1=zeros(2,m+1,k+1); SY_PSI1=zeros(n+1,2,k+1);
% SX_PSI2=zeros(2,m+1,k+1); SY_PSI2=zeros(n+1,2,k+1);
% BBF=zeros(2*(n-1)*(m-1),k+1);
% for kkf=1:1:k+1
%     ta1=tk(kkf);
%     Sx_psi1=[psi1(Xi(1,:),Yj(1,:),ta1);psi1(Xi(end,:),Yj(end,:),ta1)];
%     Sy_psi1=[psi1(Xi(:,1),Yj(:,1),ta1) psi1(Xi(:,end),Yj(:,end),ta1)];
% 
%     Sx_psi2=[psi2(Xi(1,:),Yj(1,:),ta1);psi2(Xi(end,:),Yj(end,:),ta1)];
%     Sy_psi2=[psi2(Xi(:,1),Yj(:,1),ta1) psi2(Xi(:,end),Yj(:,end),ta1)];
% 
%     SX_PSI1(:,:,kkf)=Sx_psi1; SY_PSI1(:,:,kkf)=Sy_psi1;
%     SX_PSI2(:,:,kkf)=Sx_psi2; SY_PSI2(:,:,kkf)=Sy_psi2;
% 
%     Bxf1=[Sx_psi1(1,2:end-1);zeros((n-1)-2,m-1);Sx_psi1(end,2:end-1)];
%     Bxf1=cfive*Bxf1(:);
% 
%     Bxf2=[Sx_psi2(1,2:end-1);zeros((n-1)-2,m-1);Sx_psi2(end,2:end-1)];
%     Bxf2=cfive*Bxf2+Rs(Xio,Yjo,ta1); Bxf2=Bxf2(:);
%     
%     Byf1=bfive*[Sy_psi1(2:end-1,1);zeros((n-1)*(m-3),1);Sy_psi1(2:end-1,end)];
%     Byf2=bfive*[Sy_psi2(2:end-1,1);zeros((n-1)*(m-3),1);Sy_psi2(2:end-1,end)];
% 
%     BBF(:,kkf)=[Bxf1;Bxf2]+[Byf1;Byf2]+betaK*[Yjo(:);Yjo(:)];
% end

% Bf=[q1_3;q2_3]-([Bxf1;Bxf2]+[Byf1;Byf2]+betaK*[Yjo(:);Yjo(:)]);

Bf=[q1NUM;q2NUM]-BBF;


%% 循环追赶法
% Am1=afive*ones((n-1),(m-1))-F1;   Au1=cfive*ones((n-2),(m-1)); Ad1=cfive*ones((n-2),(m-1));
% Auu1=bfive*ones((n-1),(m-1)-1);   Add1=bfive*ones((n-1),(m-1)-1);         
% 
% Am2=afive*ones((n-1),(m-1))-F2;   Au2=cfive*ones((n-2),(m-1)); Ad2=cfive*ones((n-2),(m-1));
% Auu2=bfive*ones((n-1),(m-1)-1);   Add2=bfive*ones((n-1),(m-1)-1); 
% 
% FF1=F1;FF2=F2; 
% diagC=cfive*ones((n-2),(m-1)); diagC=[diagC;zeros(1,(m-1))]; diagC=diagC(:);

%% 追赶法7对角矩阵
% tzhui=tic;
% 五对角矩阵
% BETAA=(diag(Am1(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
%       +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1)))\diag(F1*ones((n-1)*(m-1),1));
% Temp1=(diag(Am1(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
%       +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1)));
% BETAA=Temp1\diag(F1*ones((n-1)*(m-1),1));
% Temp2=Temp1; Temp2(1:(n-1)*(m-1)+1:end)=afive-F2;

YY1=BETAA*(Bf(1:(n-1)*(m-1),:)/F1);
psi_NUM2=(Temp2-F2*BETAA)\(Bf((n-1)*(m-1)+1:end,:)-F2*YY1);

psi_NUM1=YY1-BETAA*psi_NUM2;
% toc(tzhui)

%% 迭代法
% Bf1=Bf(1:(n-1)*(m-1),:);   Bf2=Bf((n-1)*(m-1)+1:end,:);
% Bf1=Bf1-TKK*F1*psi0_2k;
% % Temp1=(diag(Am1(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
% %       +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1))); 
% % Temp1=TKK*Temp1;
% psi_NUM1=Temp1\Bf1;
% 
% % Temp2=(diag(Am2(:))+diag(diagC(1:end-1,1),1)+diag(diagC(1:end-1,1),-1)...
% %       +diag(Auu1(:),n-1)+diag(Add1(:),-(n-1))); 
% % Temp2=TKK*Temp2;
% Bf2=Bf2-TKK*F2*psi_NUM1;
% psi_NUM2=Temp2\Bf2;
%%
if any(isnan(Bf(:)))
    psi_NUM1=NaN*ones((n-1)*(m-1),k+1);
    psi_NUM2=NaN*ones((n-1)*(m-1),k+1);
end
psi_NUM=[psi_NUM1;psi_NUM2];
psi0_1k=psi_NUM(1:(n-1)*(m-1),:); psi0_2k=psi_NUM((n-1)*(m-1)+1:end,:);

%% 下一步psi偏微分
dxpsi0_1k=zeros((n-1)*(m-1),k+1);
dxpsi0_2k=zeros((n-1)*(m-1),k+1);
dypsi0_1k=zeros((n-1)*(m-1),k+1);
dypsi0_2k=zeros((n-1)*(m-1),k+1);

dxpsi0m_1k=zeros((n-1),2,k);
dxpsi0m_2k=zeros((n-1),2,k);

psi0_1=zeros(n+1,m+1);psi0_2=zeros(n+1,m+1);
psi0_11=zeros(n+1,m+1);psi0_22=zeros(n+1,m+1);
DCy=diag(-ones(n+1,1))+diag(ones(n-1,1),2); DCy=DCy(1:end-2,:)/(2*hy);             % (n-1)*(n+1) x (n+1)*(m-1) = (n-1)*(m-1)
DCx=diag(-ones(m+1,1))+diag(ones(m-1,1),-2); DCx=DCx(:,1:end-2)/(2*hx);            % (n-1)*(m+1) x (m+1)*(m-1) = (n-1)*(m-1)

%
DDpsi2k=zeros((n-1)*(m-1),k);
DDpsi2kk=zeros(m-1,n-1);

DDpsi1d2k=zeros((n-1)*(m-1),k);
DDpsi1d2kk=zeros(m-1,n-1);
for kk=1:1:k
    psi0_1(2:end-1,2:end-1)=reshape(psi0_1k(:,kk),[n-1,m-1]);
    psi0_2(2:end-1,2:end-1)=reshape(psi0_2k(:,kk),[n-1,m-1]);
    psi0_11(2:end-1,2:end-1)=reshape(psi0_1k(:,kk+1),[n-1,m-1]);
    psi0_22(2:end-1,2:end-1)=reshape(psi0_2k(:,kk+1),[n-1,m-1]);

    Sx_psi1=SX_PSI1(:,:,kk);
    Sy_psi1=SY_PSI1(:,:,kk);
    Sx_psi2=SX_PSI2(:,:,kk);
    Sy_psi2=SY_PSI2(:,:,kk);

    Sx_psi11=SX_PSI1(:,:,kk+1);
    Sy_psi11=SY_PSI1(:,:,kk+1);
    Sx_psi22=SX_PSI2(:,:,kk+1);
    Sy_psi22=SY_PSI2(:,:,kk+1);

    psi0_1([1 n+1],:)=Sx_psi1;  psi0_1(:,[1 m+1])=Sy_psi1;
    psi0_2([1 n+1],:)=Sx_psi2;  psi0_2(:,[1 m+1])=Sy_psi2;

    psi0_11([1 n+1],:)=Sx_psi11;  psi0_11(:,[1 m+1])=Sy_psi11;
    psi0_22([1 n+1],:)=Sx_psi22;  psi0_22(:,[1 m+1])=Sy_psi22;

    dxpsi1=psi0_1(2:end-1,:)*DCx;
    dypsi1=DCy*psi0_1(:,2:end-1);

    dxpsi2=psi0_2(2:end-1,:)*DCx;
    dypsi2=DCy*psi0_2(:,2:end-1);

    dxpsi11=psi0_11(2:end-1,:)*DCx;
    dypsi11=DCy*psi0_11(:,2:end-1);

    dxpsi22=psi0_22(2:end-1,:)*DCx;
    dypsi22=DCy*psi0_22(:,2:end-1);

    dxpsi0_1k(:,kk)=0.5*(dxpsi1(:)+dxpsi11(:)); dypsi0_1k(:,kk)=0.5*(dypsi1(:)+dypsi11(:));
    dxpsi0_2k(:,kk)=0.5*(dxpsi2(:)+dxpsi22(:)); dypsi0_2k(:,kk)=0.5*(dypsi2(:)+dypsi22(:));
    %% 边界psi值
    dxpsix01=(psi0_1(2:end-1,2)-psi0_1(2:end-1,1))/hx;
    dxpsix011=(psi0_11(2:end-1,2)-psi0_11(2:end-1,1))/hx;
    dxpsi0m_1k(:,1,kk)=0.5*(dxpsix01+dxpsix011);

    dxpsixm1=(psi0_1(2:end-1,end)-psi0_1(2:end-1,end-1))/hx;
    dxpsixm11=(psi0_11(2:end-1,end)-psi0_11(2:end-1,end-1))/hx;
    dxpsi0m_1k(:,end,kk)=0.5*(dxpsixm1+dxpsixm11);

    dxpsix02=(psi0_2(2:end-1,2)-psi0_2(2:end-1,1))/hx;
    dxpsix022=(psi0_22(2:end-1,2)-psi0_22(2:end-1,1))/hx;
    dxpsi0m_2k(:,1,kk)=0.5*(dxpsix02+dxpsix022);

    dxpsixm2=(psi0_2(2:end-1,end)-psi0_2(2:end-1,end-1))/hx;
    dxpsixm22=(psi0_22(2:end-1,end)-psi0_22(2:end-1,end-1))/hx;
    dxpsi0m_2k(:,end,kk)=0.5*(dxpsixm2+dxpsixm22);

    %% 加拉普拉斯项
%     DDpsi2kk=DD_psi2(Xio,Yjo,tk(kk)+tau/2)';    % x排列
%     DDpsi2k(:,kk)=DDpsi2kk(:);
    % 精确算
%     ta=tk(kk);
%     DDpsi2kk=q2(Xio,Yjo,ta+tau/2)'+F2*(psi2(Xio,Yjo,ta+tau/2)'-psi1(Xio,Yjo,ta+tau/2)')-betaK*Yjo'-Rs(Xio,Yjo,ta+tau/2)';    % x排列
%     DDpsi2k(:,kk)=DDpsi2kk(:);
%     DDpsi1d2kk=DD_psi1(Xio,Yjo,ta+tau/2)'-DD_psi2(Xio,Yjo,ta+tau/2)';
%     DDpsi1d2k(:,kk)=DDpsi1d2kk(:);

    % 数值算 Laplcian psi2
    psi0_2S=0.5*(psi0_2+psi0_22);
    for i=1:1:m-1
        for j=1:1:n-1
            DDpsi2kk(i,j)=(psi0_2S((j+1),(i+1)-1)+psi0_2S((j+1),(i+1)+1))*Dxfive+...
                (psi0_2S((j+1)-1,(i+1))+psi0_2S((j+1)+1,(i+1)))*Dyfive+...
                psi0_2S((j+1),(i+1))*Dxyfive;
        end
    end
    DDpsi2k(:,kk)=DDpsi2kk(:);

    % 数值算 Laplcian (psi1-psi2)
    psi0_1d2S=0.5*(psi0_1+psi0_11)-psi0_2S;
    for i=1:1:m-1
        for j=1:1:n-1
            DDpsi1d2kk(i,j)=(psi0_1d2S((j+1),(i+1)-1)+psi0_1d2S((j+1),(i+1)+1))*Dxfive+...
                (psi0_1d2S((j+1)-1,(i+1))+psi0_1d2S((j+1)+1,(i+1)))*Dyfive+...
                psi0_1d2S((j+1),(i+1))*Dxyfive;
        end
    end
    DDpsi1d2k(:,kk)=DDpsi1d2kk(:);
end
%% 数值计算需要添加
DDpsi2k=DDpsi2k/hx/hy; %数值计算需要添加
DDpsi1d2k=DDpsi1d2k/hx/hy;   %数值计算需要添加

% 等式计算拉普拉斯项
% DDpsi1d2k=q1NUM-q2NUM+(F1+F2)*(psi0_1k-psi0_2k);
% DDpsi1d2k=0.5*(DDpsi1d2k(:,1:end-1)+DDpsi1d2k(:,2:end))/Ro;
% 
% DDpsi2k=q2NUM+F2*(psi0_2k-psi0_1k-Yjo(:));
% DDpsi2k=0.5*(DDpsi2k(:,1:end-1)+DDpsi2k(:,2:end))/Ro;

%% t=T/2
cbp1=[-0.694407288240248	1.51245125508571];
cbp2=[0.456067242804967	1.58236089335150];
tk_ms=k/2+1;
% tk_ms=k;
psi_num=psi_NUM(:,tk_ms);
ta1=tk(tk_ms);
psi_01=psi1(Xi,Yj,ta1); 
psi_02=psi2(Xi,Yj,ta1); 

psi_num1=reshape(psi_num(1:(n-1)*(m-1)),[n-1,m-1]);
psi_num2=reshape(psi_num((n-1)*(m-1)+1:end),[n-1,m-1]);
psi_num1=[SY_PSI1(:,1,tk_ms) [SX_PSI1(1,2:end-1,tk_ms);psi_num1;SX_PSI1(end,2:end-1,tk_ms)] SY_PSI1(:,end,tk_ms)];
psi_num2=[SY_PSI2(:,1,tk_ms) [SX_PSI2(1,2:end-1,tk_ms);psi_num2;SX_PSI2(end,2:end-1,tk_ms)] SY_PSI2(:,end,tk_ms)];

error_absf1=abs(psi_num1-psi_01);
error_absf2=abs(psi_num2-psi_02);

figure(5)
subplot(2,2,1)
surf(Xi,Yj,psi_01,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Exact solution of \psi_1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,psi_num1,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Numerical result of \psi_1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,error_absf1,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_1: t=' num2str(ta1)])
view(0,90)

figure(6)
subplot(2,2,1)
surf(Xi,Yj,psi_02,'EdgeColor','none')
colorbar
title(['Exact solution of \psi_2: t=' num2str(ta1)])
view(0,90)
caxis(cbp2)
subplot(2,2,2)
surf(Xi,Yj,psi_num2,'EdgeColor','none')
colorbar
title(['Numerical result of \psi_2: t=' num2str(ta1)])
view(0,90)
caxis(cbp2)
subplot(2,2,3)
surf(Xi,Yj,error_absf2,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_2: t=' num2str(ta1)])
view(0,90)
% subplot(2,2,4)
% surf(Xi,Yj,psi2_001,'EdgeColor','none')
% colorbar
% title('Initial value of \psi_2')
% view(0,90)


%% t=T
tk_ms=k+1;
% tk_ms=k/2+1;
psi_num=psi_NUM(:,tk_ms);
ta1=tk(tk_ms);
psi_01=psi1(Xi,Yj,ta1); 
psi_02=psi2(Xi,Yj,ta1); 

psi_num1=reshape(psi_num(1:(n-1)*(m-1)),[n-1,m-1]);
psi_num2=reshape(psi_num((n-1)*(m-1)+1:end),[n-1,m-1]);
psi_num1=[SY_PSI1(:,1,tk_ms) [SX_PSI1(1,2:end-1,tk_ms);psi_num1;SX_PSI1(end,2:end-1,tk_ms)] SY_PSI1(:,end,tk_ms)];
psi_num2=[SY_PSI2(:,1,tk_ms) [SX_PSI2(1,2:end-1,tk_ms);psi_num2;SX_PSI2(end,2:end-1,tk_ms)] SY_PSI2(:,end,tk_ms)];

error_absf1=abs(psi_num1-psi_01);
error_absf2=abs(psi_num2-psi_02);

figure(7)
subplot(2,2,1)
surf(Xi,Yj,psi_01,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Exact solution of \psi_1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,psi_num1,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Numerical result of \psi_1: t=' num2str(ta1)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,error_absf1,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_1: t=' num2str(ta1)])
view(0,90)
% subplot(2,2,4)
% surf(Xi,Yj,psi1_001,'EdgeColor','none')
% colorbar
% % caxis(cbp)
% title('Initial guass of \psi_1')
% view(0,90)

figure(8)
subplot(2,2,1)
surf(Xi,Yj,psi_02,'EdgeColor','none')
colorbar
title(['Exact solution of \psi_2: t=' num2str(ta1)])
view(0,90)
caxis(cbp2)
subplot(2,2,2)
surf(Xi,Yj,psi_num2,'EdgeColor','none')
colorbar
title(['Numerical result of \psi_2: t=' num2str(ta1)])
view(0,90)
caxis(cbp2)
subplot(2,2,3)
surf(Xi,Yj,error_absf2,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_2: t=' num2str(ta1)])
view(0,90)
% subplot(2,2,4)
% surf(Xi,Yj,psi2_001,'EdgeColor','none')
% colorbar
% title('Initial value of \psi_2')
% view(0,90)
toc
% clear psi_NUM q1NUM q2NUM  Bf

%% 全局误差
% Q1EX=zeros((n-1)*(m-1),k+1);
% Q2EX=zeros((n-1)*(m-1),k+1);
% PSI_EX1=zeros((n-1)*(m-1),k+1);
% PSI_EX2=zeros((n-1)*(m-1),k+1);
% 
% for kk=1:1:k+1
%     TempQ=q1(Xio,Yjo,tk(kk));
%     Q1EX(:,kk)=TempQ(:);
%     TempQ=q2(Xio,Yjo,tk(kk));
%     Q2EX(:,kk)=TempQ(:);
% 
%     TempP=psi1(Xio,Yjo,tk(kk));
%     PSI_EX1(:,kk)=TempP(:);
%     TempP=psi2(Xio,Yjo,tk(kk));
%     PSI_EX2(:,kk)=TempP(:);
% end

error_absf1=abs(PSI_EX1-psi_NUM1); 
error_absf2=abs(PSI_EX2-psi_NUM2);
abs_error1=abs(Q1EX-q1NUM);         
abs_error2=abs(Q2EX-q2NUM);

psi_num=psi_NUM(:);