%% q:t=T/2
figure(1)
subplot(2,2,1)
surf(Xi,Yj,q1_k1h,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Exact solution of q1: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q1_num1h,'EdgeColor','none')
colorbar
caxis(cbq1)
title(['Numerical result of q1: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error1h,'EdgeColor','none')
colorbar
title(['Absolution error of q1: t=' num2str(ta1h)])
view(0,90)


figure(2)
subplot(2,2,1)
surf(Xi,Yj,q2_k1h,'EdgeColor','none')
colorbar
caxis(cbq2)
% caxis(cb)
title(['Exact solution of q2: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,q2_num1h,'EdgeColor','none')
colorbar
% caxis(cb)
title(['Numerical result of q2: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,abs_error2h,'EdgeColor','none')
colorbar
title(['Absolution error of q2: t=' num2str(ta1h)])
view(0,90)

%% \psi:t=T/2
figure(5)
subplot(2,2,1)
surf(Xi,Yj,psi_01h,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Exact solution of \psi_1: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,2)
surf(Xi,Yj,psi_num1h,'EdgeColor','none')
colorbar
caxis(cbp1)
title(['Numerical result of \psi_1: t=' num2str(ta1h)])
view(0,90)
subplot(2,2,3)
surf(Xi,Yj,error_absf1h,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_1: t=' num2str(ta1h)])
view(0,90)

figure(6)
subplot(2,2,1)
surf(Xi,Yj,psi_02h,'EdgeColor','none')
colorbar
title(['Exact solution of \psi_2: t=' num2str(ta1h)])
view(0,90)
caxis(cbp2)
subplot(2,2,2)
surf(Xi,Yj,psi_num2h,'EdgeColor','none')
colorbar
title(['Numerical result of \psi_2: t=' num2str(ta1h)])
view(0,90)
caxis(cbp2)
subplot(2,2,3)
surf(Xi,Yj,error_absf2h,'EdgeColor','none')
colorbar
title(['Absolution error of \psi_2: t=' num2str(ta1h)])
view(0,90)