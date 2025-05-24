%% Time Mesh
% if T>=10,the iterative process will diverge since our theoretical result (Theorem 2.3).
% if T<=8, the iterative process will converge.
T=12;         % Total time interval   
tau=1/60;    % Time step size \tau
k=T/tau;     % Number of time steps
tk=0:tau:T;

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