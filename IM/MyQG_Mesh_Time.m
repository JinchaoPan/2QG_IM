%% 网格 给定psi初值
Lx=1; Ly=1*1; 
m=50*1*1;n=50*1*1;
hx=Lx/m;hy=Ly/n;
xi=0:hx:Lx; yj=(0:hy:Ly); 
[Xi,Yj]=meshgrid(xi,yj);    
Xio=Xi(2:end-1,2:end-1);    Yjo=Yj(2:end-1,2:end-1);
%% 
TT=12;       % Total time
DT=12;       % Sub-time interval
kT=TT/DT;    % [0,TT] is divided into kT sub-time intervals where kT is Integer.

DTk=400;     % 子时间分段次数，时间步长： DT/DTK %400
%% 初始q
q1_00T=q1(Xio,Yjo,0);
q2_00T=q2(Xio,Yjo,0);


