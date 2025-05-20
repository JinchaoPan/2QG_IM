clear   % 考虑sigma不为0，A不为0情况，加入拉普拉斯项(差分)；雅可比项变号+考虑Ro不为1
clc
tstart=tic;
Myfunction5
% Myfunction7
MyQG_Mesh_Time

for TTk=0:DT:(TT-DT)       % TTk 为每一段时间的初始时刻
    MyQG_Mesh_Nicolson
    %%
    MyQG_Nicolson
    % MyQG_Filter_FVM0
    MyQG_FiveLapNicolson
    %%
    kqg=2; XK=10*1;
    ERROR_absp1=zeros(1,XK);
    ERROR_absp2=zeros(1,XK);
    ERROR_absp1(:,1)=max(error_absf1(:));
    ERROR_absp2(:,1)=max(error_absf2(:));

    ERROR_absQ1=zeros(1,XK);
    ERROR_absQ2=zeros(1,XK);
    ERROR_absQ1(:,1)=max(abs_error1(:));
    ERROR_absQ2(:,1)=max(abs_error2(:));

    PSI_NUMTK=zeros(2*(n-1)*(m-1),XK);
    DiffPsi=zeros(2*(n-1)*(m-1),XK);
    max_errork=1;
    while kqg<=XK && max_errork>=1e-6
%     while kqg<=XK
        MyQG_Nicolson1
        %     MyQG_Filter_FVM0
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
    %% 下一段初值
    q1_00T=reshape(q1NUM(:,end),[n-1,m-1]);
    q2_00T=reshape(q2NUM(:,end),[n-1,m-1]);
end
toc(tstart)

