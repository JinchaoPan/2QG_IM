%% Time Mesh
TT=12;       % Total time interval
TL=4;        % Time Sub-interval
L=TT/TL;     % Number of Time Sub-intervals

tau=1/60;    % Time step size \tau
DTk=TL/tau;  % Number of time steps in Time Sub-intervals

%% Spatial Mesh
Lx=1; Ly=1*1; 
m=50*1;n=50*1;
hx=Lx/m;hy=Ly/n;
xi=0:hx:Lx; yj=(0:hy:Ly); 
[Xi,Yj]=meshgrid(xi,yj);    
Xio=Xi(2:end-1,2:end-1);    Yjo=Yj(2:end-1,2:end-1);

%% Initial value of q
q1_00T=q1(Xio,Yjo,0);
q2_00T=q2(Xio,Yjo,0);
%%
Lih=ceil(L/2);Li=1;
odev=mod(L,2);
XK=18;