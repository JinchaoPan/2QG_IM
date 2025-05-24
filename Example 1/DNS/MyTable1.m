%% Output some data about (EMQ,EMP) of DNS in Table 1
Q1EX=zeros((n-1)*(m-1),k+1); 
Q2EX=zeros((n-1)*(m-1),k+1);
PSI_EX1=zeros((n-1)*(m-1),k+1);
PSI_EX2=zeros((n-1)*(m-1),k+1);

for kk=1:1:k+1
    TempQ=q1(Xio,Yjo,tk(kk));
    Q1EX(:,kk)=TempQ(:);
    TempQ=q2(Xio,Yjo,tk(kk));
    Q2EX(:,kk)=TempQ(:);

    TempP=psi1(Xio,Yjo,tk(kk));
    PSI_EX1(:,kk)=TempP(:);
    TempP=psi2(Xio,Yjo,tk(kk));
    PSI_EX2(:,kk)=TempP(:);
end

%% 相对最大误差
AErrPSI1=abs(PSI_EX1-psi_NUM1); EMP1=max(AErrPSI1(:))./max(abs(PSI_EX1(:)));
AErrPSI2=abs(PSI_EX2-psi_NUM2); EMP2=max(AErrPSI2(:))./max(abs(PSI_EX2(:)));
AErrQ1=abs(Q1EX-q1NUM);         EMQ1=max(AErrQ1(:))./max(abs(Q1EX(:)));
AErrQ2=abs(Q2EX-q2NUM);         EMQ2=max(AErrQ2(:))./max(abs(Q2EX(:)));

%%
Table1=[EMQ1 EMQ2 EMP1 EMP2 ]