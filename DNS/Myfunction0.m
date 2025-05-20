bpk=1; BCC=0;
% betaK=1;
betaK=1;
Rs=@(x,y,t)(0*x);
Rs_x=@(x,y,t)(0*x); Rs_y=@(x,y,t)(0*x); 
Rs_D=@(x,y,t)(0*x);

%%
Jk=1;        % Jk为偶数时为原方程，Jk为奇数时为雅可比项和dq/dt在同一边
Sigma=0.005*1;
Ro=1.5e-0; Re=450; A=Ro/Re;
% Ro=1.5; Re=450; A=0;
%%
% F1=1;F2=2;
Fr=0.1*1e-0;         deltaF=2/3;
F1=Fr/deltaF;   F2=Fr/(1-deltaF);
%% 配置psi
pa=2;
psi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*y.*exp(x)+cos(y).*(2-x.^2)))*bpk+BCC;
psi2=@(x,y,t)(exp(sin(t)).*(cos(pa*t).*(x-0.5).*exp(y)+cos(x).*(2-(y-0.5).^2)))*bpk+BCC;
Dx_psi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*y.*exp(x)-2*x.*cos(y)))*bpk; 
Dy_psi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*exp(x)-sin(y).*(2-x.^2)))*bpk;
Dx_psi2=@(x,y,t)(exp(sin(t)).*(cos(pa*t).*exp(y)-sin(x).*(2-(y-0.5).^2)))*bpk; 
Dy_psi2=@(x,y,t)(exp(sin(t)).*(cos(pa*t).*(x-0.5).*exp(y)-2*(y-0.5).*cos(x)))*bpk;

DD_psi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*y.*exp(x)-2*cos(y)-cos(y).*(2-x.^2)))*bpk;  
DD_psi2=@(x,y,t)(exp(sin(t)).*(-cos(x).*(2-(y-0.5).^2)+cos(pa*t).*(x-0.5).*exp(y)-2*cos(x)))*bpk;

%% 
q1=@(x,y,t)(Ro*DD_psi1(x,y,t)-F1.*(psi1(x,y,t)-psi2(x,y,t))+betaK.*y);
q2=@(x,y,t)(Ro*DD_psi2(x,y,t)-F2.*(psi2(x,y,t)-psi1(x,y,t))+betaK.*y+Rs(x,y,t));

%% 配置q的导数
DT1_q=@(x,y,t)(y.*exp(x).*cos(pa*t).*pa);
DT2_q=@(x,y,t)((x-0.5).*exp(y).*(-sin(pa*t).*pa));
Dt_q1=@(x,y,t)((Ro*DD_psi1(x,y,t)-F1.*(psi1(x,y,t)-psi2(x,y,t))).*cos(t)+exp(sin(t)).*...
               (Ro*DT1_q(x,y,t)-F1.*(DT1_q(x,y,t)-DT2_q(x,y,t))))*bpk; % 与psi关于t导数有关
Dt_q2=@(x,y,t)((Ro*DD_psi2(x,y,t)-F2.*(psi2(x,y,t)-psi1(x,y,t))).*cos(t)+exp(sin(t)).*...
               (Ro*DT2_q(x,y,t)-F2.*(DT2_q(x,y,t)-DT1_q(x,y,t))))*bpk;

Dx_DDpsi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*y.*exp(x)-cos(y).*(-2*x)));
Dy_DDpsi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*exp(x)+2*sin(y)+sin(y).*(2-x.^2)));

Dx_DDpsi2=@(x,y,t)(exp(sin(t)).*(sin(x).*(2-(y-0.5).^2)+cos(pa*t).*exp(y)+2*sin(x)));
Dy_DDpsi2=@(x,y,t)(exp(sin(t)).*(-cos(x).*(-2*(y-0.5))+cos(pa*t).*(x-0.5).*exp(y)));

Dx_q1=@(x,y,t)(Ro*Dx_DDpsi1(x,y,t)-F1*(Dx_psi1(x,y,t)-Dx_psi2(x,y,t)))*bpk; 
Dy_q1=@(x,y,t)(Ro*Dy_DDpsi1(x,y,t)-F1*(Dy_psi1(x,y,t)-Dy_psi2(x,y,t)))*bpk+betaK; 

Dx_q2=@(x,y,t)(Ro*Dx_DDpsi2(x,y,t)-F2*(Dx_psi2(x,y,t)-Dx_psi1(x,y,t)))*bpk+Rs_x(x,y,t); 
Dy_q2=@(x,y,t)(Ro*Dy_DDpsi2(x,y,t)-F2*(Dy_psi2(x,y,t)-Dy_psi1(x,y,t)))*bpk+betaK+Rs_y(x,y,t); 

%%
f1=@(x,y,t)(Dt_q1(x,y,t)-((-1)^(Jk))*(Dy_psi1(x,y,t).*Dx_q1(x,y,t)-Dx_psi1(x,y,t).*Dy_q1(x,y,t)));
f2=@(x,y,t)(Dt_q2(x,y,t)-((-1)^(Jk))*(Dy_psi2(x,y,t).*Dx_q2(x,y,t)-Dx_psi2(x,y,t).*Dy_q2(x,y,t)))+Sigma*DD_psi2(x,y,t);
%%
DD_DDpsi1=@(x,y,t)(exp(sin(t)).*(sin(pa*t).*y.*exp(x)+(6-x.^2).*cos(y)));
DD_DDpsi2=@(x,y,t)(exp(sin(t)).*(cos(pa*t).*(x-0.5).*exp(y)+(6-(y-0.5).^2).*cos(x)));

f1=@(x,y,t)(f1(x,y,t)-A*DD_DDpsi1(x,y,t));
f2=@(x,y,t)(f2(x,y,t)-A*DD_DDpsi2(x,y,t));