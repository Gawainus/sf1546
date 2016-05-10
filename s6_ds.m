u_0=100;
r_start=1;
r_end=2;

alpha = 100
deltabetaepsilon = (5.67 * 10^-8)*0.3
k = 20


N=25
u = ones(N*2,1)
u_old = zeros(N*2,1)

while abs(u(end)-u_old(end)) > 0.5
u_old = u
% number of sections
N=N*2;

h=(r_end-r_start)/N

startguess_u = 100:200/(N-1):300

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
A(N,N)=-((2*r(end))/h^2 + alpha + (2*h*r(end)*alpha)/(k*h^2));
A(N,N-1)=(2*r(end))/h^2;

corr = ones(N,1)



u = startguess_u'

itercheck = 1
dominance_radiation = -1

Te = 200

while dominance_radiation <= 0
  Te = Te + 100
  supplement = alpha*Te - deltabetaepsilon*u(end)^4 + deltabetaepsilon*Te^4 + (2*h*r(end)*alpha*Te)/(k*h^2) - (2*h*r(end)*deltabetaepsilon*u(end)^4)/(k*h^2) + (2*h*r(end)*deltabetaepsilon*Te^4)/(k*h^2);
  supplement_jaco = -4*deltabetaepsilon*u(end)^3 - (8*h*r(end)*deltabetaepsilon*u(end)^3)/(k*h^2);
  corr = 1
 while norm(corr, Inf) > 1e-9
    F = A*u; F(1) = F(1) + (r(1)/h^2 - 1/(2*h))*u_0; F(end) = F(end) + supplement;
    J = A; J(end) = J(end) + supplement_jaco;
    corr = J\F;
    u = u-corr;
  end
  dominance_radiation = abs(deltabetaepsilon*(u(end)^4 - Te^4)) - abs(alpha*(u(end) - Te));
end
end
figure
plot(r,u)
title('Plot of Temperature with respect to Radius, Te=1200')
xlabel('radius, cm')
ylabel('temperature, Celsius')



