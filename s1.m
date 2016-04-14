% construct A

% error stuff
error=1;
u_2_old=0;
% number of sections
N=25;
K=1; %3.156
Te=20;

u_0=450;
r_start=1;
r_end=2;

while error>exp(-10)
    
h=(r_end-r_start)/N;

% from r0 to rN
r_raw=(1:h:2);

% from r1 to rN
r=r_raw(2:end);

% from r1 to rN-1 this is for diag plus 1
r_p1=r(1:end-1);

% from r2 to rN-1 this is from diag minus 1
r_m1=r(2:end);
r_m1(end)=0;

v=factor2(r,h);
v_p1=factor3(r_p1,h);
v_m1=factor1(r_m1,h);

D=diag(v);
D_p1=diag(v_p1,1);
D_m1=diag(v_m1,-1);

A=D+D_p1+D_m1;
A(N,N)=-(2*r(end)/h^2+2*K*r(end)/h+K);
A(N,N-1)=2*r(end)/h^2;

% construct b
b=zeros(N,1);
b(1)=-(r(1)/h^2-1/(2*h))*450;

%-KTe(2h/h+1)
b(N)=-K*Te*(2*r(end)/h+1);

u=A\b;
error=abs(u(end)-u_2_old);
u_2_old=u(end);
N=N*2;
end

figure
plot(r,u)
title('Plot of Temperature with respect to Radius, K=1')
xlabel('radius, cm')
ylabel('temperature, Celsius')


