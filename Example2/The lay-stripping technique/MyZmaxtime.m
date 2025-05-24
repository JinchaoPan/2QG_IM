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
AErrPSI1=abs(PSI_EX1-psi_NUM1); AErrPSI1=max(AErrPSI1(:));
AErrPSI2=abs(PSI_EX2-psi_NUM2); AErrPSI2=max(AErrPSI2(:));
AErrQ1=abs(Q1EX-q1NUM);         AErrQ1=max(AErrQ1(:));
AErrQ2=abs(Q2EX-q2NUM);         AErrQ2=max(AErrQ2(:));

% AErrPSI1=abs(PSI_EX1-psi_NUM1); AErrPSI1=max(AErrPSI1(:,end));
% AErrPSI2=abs(PSI_EX2-psi_NUM2); AErrPSI2=max(AErrPSI2(:,end));
% AErrQ1=abs(Q1EX-q1NUM);         AErrQ1=max(AErrQ1(:,end));
% AErrQ2=abs(Q2EX-q2NUM);         AErrQ2=max(AErrQ2(:,end));