%Code Method Gauss_Seidel METH
% Abraham Blanco  1223970
clear all; clc;
tic
%% Parameters
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;
% Define the number of points on the interior (this does not include the
% exterior boundary points)
M=input('Value of X Intenal Nodes=');
N=input('Value of Y Internal Nodes=');
M1=M+2;
N1=N+2;
% this generates the x and y values that will be used to calculate 
xvalues = linspace(0,2*pi,M+2);
yvalues = linspace(0,2*pi,N+2);
%%
%U matrix (guess)
U = ones(M+2,N+2);
%solving for right hand side with F equation
for i=1:length(xvalues);
    for j=1:length(yvalues);
F(i,j) = cos ( (0.5*pi)* (2*((xvalues(i)-ax) / (bx - ax))+1 )).*sin( pi*((yvalues(j)-ay) / (by -ay)));
    end
end

%% Boundary Conditions for "top" and "bottom"

% Bottom boundary values
phi_ab = ((yvalues - ay).^2 ) .* sin( pi *(yvalues - ay) / (2*(by-ay)) ) ; 

% Top boundary values
psy_ab = cos (pi*(yvalues-ay)).*cosh(by-yvalues);

% place these known values in the solution grid
U(1,:) = phi_ab;
U(end,:) = psy_ab; 

DX = 2*pi/(M+1);
A = 1/DX.^2
DY = 2*pi/(N+1);
B = 1/DY.^2
DEN = -2*(A+B)

% normalize elements
A = A/DEN;
B = B/DEN;
F = F/DEN;
DEN = 1;
error=10;
error_iterations=0
% check for diagonal dominance of elements
abs(DEN) >= abs(2*A+2*B)
while error>10^-10;
    W=U;
for P = 1:1000;

for j = 2:M+1;
    
    % Left boundary
    U(j,1) = DEN*(  F(j,1) - (2*B)*U(j,2) - A*U(j-1,1) - A*U(j+1,1) );
    % Right Boundary
    U(j,end) = DEN*(  F(j,end) - (2*B)*U(j,end-1) - A*U(j-1,end) - A*U(j+1,end) );
end
%% Main Sweep of Gauss-Siedel

for k = 2:N+1;
    for j = 2:M+1;
        U(j,k) = DEN*(  F(j,k) - B*U(j,k-1) - B*U(j,k+1)- A*U(j-1,k) - A*U(j+1,k) );
    end
end
end
error=abs(max(max(((W-U)./W))));
error_iterations=error_iterations+1
end
toc
figure
subplot(1,2,1),surf(U)
subplot(1,2,2),contour(U)