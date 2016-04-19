% construct A

% number of sections
N=3200;
Te=20;

u_0=450;
r_start=1;
r_end=2;
    
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

D=diag([v 0]);
D_p1=diag([v_p1 0],1);
D_m1=diag([v_m1 0],-1);

A=D+D_p1+D_m1;

A(N,N-1)=0;
A(N,N)=1;

A(N+1,N-1)=2*r(end)/h^2;
A(N+1,N+1)=Te*(2*r(end)/h)-200*r(end)/h-100;


% construct b
b=zeros(N+1,1);
b(1)=-(r(1)/h^2-1/(2*h))*u_0;

b(N)=100;
% 200r/h^2
b(N+1)=200*r(end)/h^2;


u=A\b;


display(u);


