% Applying the lay-stripping technique
% If you notice any difference with the Figure in the manuscript, 
% please maximize the Figure window to observe more details.

% The total time interval T and time sub-interval TL=T/L  
% can be seen in file "MyQG_Mesh_Time".
clear   
clc
tstart=tic;
Myfunction5
MyQG_Mesh_Time
ERROR_ABSP1=zeros(L,XK);ERROR_ABSP2=zeros(L,XK);
ERROR_ABSQ1=zeros(L,XK);ERROR_ABSQ2=zeros(L,XK);
for Ti=0:TL:(TT-TL)       
    MyQG_Mesh_Nicolson
    %%
    MyQG_Nicolson
    MyQG_FiveLapNicolson
    %%
    kqg=2; 
    ERROR_absp1=zeros(1,XK); 
    ERROR_absp2=zeros(1,XK);
    ERROR_absp1(:,1)=max(error_absf1(:));
    ERROR_absp2(:,1)=max(error_absf2(:));
    
    ERROR_absQ1=zeros(1,XK);
    ERROR_absQ2=zeros(1,XK);
    ERROR_absQ1(:,1)=max(abs_error1(:));
    ERROR_absQ2(:,1)=max(abs_error2(:));

%     PSI_NUMTK=zeros(2*(n-1)*(m-1),XK);
%     DiffPsi=zeros(2*(n-1)*(m-1),XK);
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
    ERROR_ABSP1(Li,:)=ERROR_absp1;
    ERROR_ABSP2(Li,:)=ERROR_absp2;
    ERROR_ABSQ1(Li,:)=ERROR_absQ1;
    ERROR_ABSQ2(Li,:)=ERROR_absQ2;
    %% 下一段初值
    q1_00T=reshape(q1NUM(:,end),[n-1,m-1]);
    q2_00T=reshape(q2NUM(:,end),[n-1,m-1]);
   
    Li=Li+1;
end
MyQG_QPhfigure
MyFigure9

toc(tstart)

