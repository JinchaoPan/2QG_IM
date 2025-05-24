%% Set initial guesses (\psi^0) 
delta=20;
% psi001=@(x,y,t)(psi1(x,y,t)+delta*sin((x.*y).*(x-Lx).*(y-Ly)));
% psi002=@(x,y,t)(psi2(x,y,t)+delta*sin((x.*y).*(x-Lx).*(y-Ly)));

psi001=@(x,y,t)(psi1(x,y,t)+delta*sin((x.*y).*(x-Lx).*(y-Ly))*cos(t));
psi002=@(x,y,t)(psi2(x,y,t)+delta*((x.*y).*(x-Lx).*(y-Ly))*sin(t));

dxpsi0_1k=zeros((n-1)*(m-1),k);
dxpsi0_2k=zeros((n-1)*(m-1),k);
dypsi0_1k=zeros((n-1)*(m-1),k);
dypsi0_2k=zeros((n-1)*(m-1),k);

dxpsi0m_1k=zeros((n-1),2,k);
dxpsi0m_2k=zeros((n-1),2,k);

% psi0_1k=psi_NUM(1:(n-1)*(m-1),:); psi0_2k=psi_NUM((n-1)*(m-1)+1:end,:);

psi0_1=zeros(n+1,m+1);    psi0_2=zeros(n+1,m+1);
psi0_11=zeros(n+1,m+1);   psi0_22=zeros(n+1,m+1);
DCy=diag(-ones(n+1,1))+diag(ones(n-1,1),2); DCy=DCy(1:end-2,:)/(2*hy);             % (n-1)*(n+1) x (n+1)*(m-1) = (n-1)*(m-1)
DCx=diag(-ones(m+1,1))+diag(ones(m-1,1),-2); DCx=DCx(:,1:end-2)/(2*hx);            % (n-1)*(m+1) x (m+1)*(m-1) = (n-1)*(m-1)

%
Dxfive=hy/(hx); Dyfive=hx/(hy);  Dxyfive=-2*(Dxfive+Dyfive);
DDpsi2k=zeros((n-1)*(m-1),k);
DDpsi2kk=zeros(m-1,n-1);

DDpsi1d2k=zeros((n-1)*(m-1),k);
DDpsi1d2kk=zeros(m-1,n-1);

psi0_1k=zeros((n-1)*(m-1),k+1);
psi0_2k=zeros((n-1)*(m-1),k+1);
temp1=psi001(Xio,Yjo,tk(1)); psi0_1k(:,1)=temp1(:);
temp2=psi002(Xio,Yjo,tk(1)); psi0_2k(:,1)=temp2(:);
%% 求半点偏导数 
for kk=1:1:k
    ta=tk(kk); ta1=tk(kk+1);
    psi0_1(2:end-1,2:end-1)=psi001(Xio,Yjo,ta);
    psi0_2(2:end-1,2:end-1)=psi002(Xio,Yjo,ta);
    psi0_11(2:end-1,2:end-1)=psi001(Xio,Yjo,ta1);
    psi0_22(2:end-1,2:end-1)=psi002(Xio,Yjo,ta1);
    temp1=psi0_11(2:end-1,2:end-1); psi0_1k(:,kk+1)=temp1(:);
    temp2=psi0_22(2:end-1,2:end-1); psi0_2k(:,kk+1)=temp2(:);

    Sx_psi1=[psi1(Xi(1,:),Yj(1,:),ta);psi1(Xi(end,:),Yj(end,:),ta)];
    Sy_psi1=[psi1(Xi(:,1),Yj(:,1),ta) psi1(Xi(:,end),Yj(:,end),ta)];
    Sx_psi2=[psi2(Xi(1,:),Yj(1,:),ta);psi2(Xi(end,:),Yj(end,:),ta)];
    Sy_psi2=[psi2(Xi(:,1),Yj(:,1),ta) psi2(Xi(:,end),Yj(:,end),ta)];

    Sx_psi11=[psi1(Xi(1,:),Yj(1,:),ta1);psi1(Xi(end,:),Yj(end,:),ta1)];
    Sy_psi11=[psi1(Xi(:,1),Yj(:,1),ta1) psi1(Xi(:,end),Yj(:,end),ta1)];
    Sx_psi22=[psi2(Xi(1,:),Yj(1,:),ta1);psi2(Xi(end,:),Yj(end,:),ta1)];
    Sy_psi22=[psi2(Xi(:,1),Yj(:,1),ta1) psi2(Xi(:,end),Yj(:,end),ta1)];

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

    %% 边界psi偏导
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
    
    %% 加拉普拉斯项： x排列
    % 精确算
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
DDpsi2k=DDpsi2k/hx/hy;   %数值计算需要添加
DDpsi1d2k=DDpsi1d2k/hx/hy;   %数值计算需要添加
%%
DD_RS=Rs_D(Xio,Yjo,tk(kk)+tau/2)';
DD_RS=DD_RS(:);