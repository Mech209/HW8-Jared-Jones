%HW8 Crank Nicolson Method to solve Heat Diffusion Equation
clear 
clc
close
dx = .5;
dt = .1;
F = 0;
D = .1;
T = 0:dt:10;
x = 0:dx:pi();
lambda = D*(dt/dx^2);
%initial condition
k = 1;
u = sin(k*x);
%boundary conditions
g0 = 0;
gL = 0;
%tridiagonal solver
for j = 2:length(x)-1
    for n = 2:length(x)-1
   u(j) = lambda/2*(u(j-1)) + (1+lambda)*u(j) -lambda/2*(u(j+1));
    end
end
%RHS
a = zeros(length(x));
a(1) = lambda*g0 + (1-lambda)*u(1) + lambda/2*(u(2));
for j = 2:length(x)-1
 a(j) = lambda/2*(u(j-1)) + (1-lambda)*u(j) + u(j+1);
end
a(end) = lambda/2*u(j) + (1-lambda)*(u(j+1)) + lambda*gL;
a = a + dt*F;