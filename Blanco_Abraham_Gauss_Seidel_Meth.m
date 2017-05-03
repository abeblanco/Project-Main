%Code Method Gauss_Seidel METHOD
% Abraham Blanco  1223970
clear all; clc;
%% Parameters
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;
% Define the number of points on the interior (this does not include the exterior boundary points)
M=input('Value of X Intenal Nodes=');
N=input('Value of Y Internal Nodes=');
tic %time begins here after inputs are set
M1=M+2;
N1=N+2;
% this generates the x and y values that will be used to calculate the F matrix
xvalues = linspace(0,2*pi,M+2);
yvalues = linspace(0,2*pi,N+2);
%%
%U matrix ( initial guess)
U = ones(M+2,N+2);
%H1= ones(M+2,N+2);
%solving for right hand side (F equation)
for i=1:length(xvalues);
    for j=1:length(yvalues);
F(i,j) = cos ( (0.5*pi)* (2*((xvalues(i)-ax) / (bx - ax))+1 )).*sin( pi*((yvalues(j)-ay) / (by -ay)));
%H(i,j) = 0; %last simulation with F=0
    end
end

%% Boundary Conditions for "Left" and "Right" side of Matrix

% Left boundary values (Dirchelet Condition)
phi_ab = ((yvalues - ay).^2 ) .* sin( pi *(yvalues - ay) / (2*(by-ay)) ) ; 

% Right boundary values (Dirchelet Condition)
psy_ab = cos (pi*(yvalues-ay)).*cosh(by-yvalues);

% place these known values in the solution grid
%solveing for Left and Right Boundaries
U(1,:) = phi_ab;
U(end,:) = psy_ab; 
%H1(1,:) = phi_ab;
%H1(end,:) = psy_ab;
%%
DX = 2*pi/(M+1);
A = 1/DX.^2;
DY = 2*pi/(N+1);
B = 1/DY.^2;
R = -2*(A+B);

% normalize elements
A = A/R;
B = B/R;
F = F/R;
%H = H/R;

R = 1;
error=10;
error_iterations=0;
save(METH.mat);
%%
load(METH.mat);
% check for diagonal dominance of elements
abs(R) >= abs(2*A+2*B);
while error>10^-10;
   W1=U;
   %W2=H1;
for j = 2:M+1;
    
    % Left boundary
    U(j,1) = (  F(j,1) - (2*B)*U(j,2) - A*U(j-1,1) - A*U(j+1,1) );
    %H1(j,1) = (  H(j,1) - (2*B)*H1(j,2) - A*H1(j-1,1) - A*H1(j+1,1) );
    % Right Boundary
    U(j,end) = (  F(j,end) - (2*B)*U(j,end-1) - A*U(j-1,end) - A*U(j+1,end) );
    %H1(j,end) = (  H(j,end) - (2*B)*H1(j,end-1) - A*H1(j-1,end) - A*H1(j+1,end) );
end
%% Main Sweep of Gauss-Siedel

for k = 2:N+1;
    for j = 2:M+1;
        U(j,k) = (  F(j,k) - B*U(j,k-1) - B*U(j,k+1)- A*U(j-1,k) - A*U(j+1,k) );
       % H1(j,k) = (  H(j,k) - B*H1(j,k-1) - B*H1(j,k+1)- A*H1(j-1,k) - A*H1(j+1,k) );
    end
end
error=abs(max(max(((W1-U)./W1))));
%error=abs(max(max(((W2-H1)./W2))));
error_iterations=error_iterations+1;
end
toc
error_iterations
figure
subplot(1,2,1),surf(U), xlabel('Y axis'), ylabel('X axis'), zlabel('Z axis'), title('F=cosx*siny')
subplot(1,2,2),contour(U), xlabel('Y axis'), ylabel('X axis'), title('F=cosx*siny')
%figure
%subplot(1,2,1),surf(H1), xlabel('Y axis'), ylabel('X axis'), zlabel('Z axis'), title('F=0')
%subplot(1,2,2),contour(H1), xlabel('Y axis'), ylabel('X axis'), title('F=0')
