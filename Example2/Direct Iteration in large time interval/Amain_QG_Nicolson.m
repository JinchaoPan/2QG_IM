% We show the direct implementations of our iteration algorithm 
% in large time interval (0,T) cannot obtain convergent iteration sequences.

% You can change time interval T in file "MyQG_Mesh_Time".

clear   
clc
tstart=tic;
Myfunction5
MyQG_Mesh_Time
MyQG_Mesh_Nicolson
%% Iteration
MyQG_Nicolson
MyQG_FiveLapNicolson
%%
kqg=2; XK=10;
ERROR_absp1=zeros(1,XK);
ERROR_absp2=zeros(1,XK);
ERROR_absp1(:,1)=max(error_absf1(:));
ERROR_absp2(:,1)=max(error_absf2(:));

ERROR_absQ1=zeros(1,XK);
ERROR_absQ2=zeros(1,XK);
ERROR_absQ1(:,1)=max(abs_error1(:));
ERROR_absQ2(:,1)=max(abs_error2(:));

PSI_NUMTK=zeros(2*(n-1)*(m-1)*(k+1),XK);
DiffPsi=zeros(2*(n-1)*(m-1)*(k+1),XK);
max_errork=1;
while kqg<=XK && max_errork>=1e-5
    MyQG_Nicolson1
    MyQG_FiveLapNicolson1
    ERROR_absp1(:,kqg)=max(error_absf1(:));
    ERROR_absp2(:,kqg)=max(error_absf2(:));

    ERROR_absQ1(:,kqg)=max(abs_error1(:));
    ERROR_absQ2(:,kqg)=max(abs_error2(:));

    PSI_NUMTK(:,kqg)=psi_num;
    DiffPsi(:,kqg)=abs(PSI_NUMTK(:,kqg)-PSI_NUMTK(:,kqg-1));
    max_errork=max(DiffPsi(:,kqg));
    kqg=kqg+1;
end
if T==12
    MyFigure6
end
MyTable2
toc(tstart)

