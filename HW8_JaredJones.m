%HW8 Crank Nicolson Method to solve Heat Diffusion Equation
clear;clc;close

F = 0;
N = 8;
L = pi();
dx = L/(N+1);
T = 10;
dt = T/(N+1);
x = 0:dx:L;
t = 0:dt:T;

%initial condition
k = 1;

%boundary conditions
g0 = 0;
gL = 0;
 
D = 0.1;


nx=N+2;   % or nx=[(10-0)/dx]+1 with L=10 cm
nt=N+2;  % to compute 10 time steps k=[1,11]
% ---Constant Coefficients of the tridiagonal system
b = D/(2*dx^2);   % Super diagonal: coefficients of u(i+1)
c = b;               % Subdiagonal: coefficients of u(i-1)
a = 1/dt+b+c;        % Main Diagonal: coefficients of u(i)
% Boundary conditions and Initial Conditions
Uo(1)=0; Uo(nx)=0;
Un =sin(k*x); 
% Store results for future use
UUU(1,:)=Uo;
% Loop over time
for k =2:nt
    for ii =1:nx-2
      if ii==1
   d(ii)= c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2)+c*Un(1);
    elseif ii == nx-2
   d(ii)=c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2)+b*Un(nx);
     else
   d(ii)=c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2);
      end
    end

% note that d is row vector
% Transform a, b, c constants in column vectors:
bb=b*ones(nx-3,1);
cc=bb;
aa=a*ones(nx-2,1);
% Use column vectors to construct tridiagonal matrices
AA=diag(aa)+ diag(-bb,1)+ diag(-cc,-1);
% Find the solution for interior nodes i=2,3,4,5
UU=AA\d';
% Build the whole solution as row vector
Un=[Un(1),UU',Un(nx)];
% Store results for future use
UUU(k,:)=Un;
end


%exact solution u(x,t) 
u_exact = exp(-D*(k^2)*t')*sin(k*x);
u_exact1 = exp(-D*(k^2)*t(T/5))*sin(k*x); %t = T/5
u_exact2 = exp(-D*(k^2)*t(T/2))*sin(k*x); %t = T/2
u_exact3 = exp(-D*(k^2)*t(T))*sin(k*x); %t = T


%error
err = 1/N*(abs((UUU - u_exact)/u_exact));

%plot

plot(x, UUU(:,[T/5 T/2 T]), x, u_exact(:,[T/5 T/2 T]))

