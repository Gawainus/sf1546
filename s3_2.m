% construct A

% error stuff
error=1;
u_2_old=0;
% number of sections
N=200;
K=2.18034%2.192something
Te=20;

u_0=458;
r_start=1;
r_end=2;

% polynomial definition
P =  @(r)(-r.^3 + 3*r.^2 - 2*r)./(-(1/3)*r.^3 + (3/2)*r.^2 - 2*r + 1)

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
v_p1=factor3_s3(r_p1,h);
v_m1=factor1_s3(r_m1,h);

D=diag(v);
D_p1=diag(v_p1,1);
D_m1=diag(v_m1,-1);

A=D+D_p1+D_m1;
A(N,N)=-(2*r(end)/h^2+2*K*r(end)/h+K);
A(N,N-1)=2*r(end)/h^2;

% construct b
b=zeros(N,1);
b(1)=-(r(1)/h^2-(1+P(r(1)))/(2*h))*u_0;

%-KTe(2h/h+1)
b(N)=-K*Te*(1 + P(r(end)) + (2*r(end))/h);

u=A\b;
error=abs(u(end)-u_2_old);
u_2_old=u(end);
N=N*2;
end

display(u(end));


