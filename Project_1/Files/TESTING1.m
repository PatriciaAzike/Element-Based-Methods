clc
clear all
close all

f=@(x) cos(pi/2 * x);
n=64;
x_npoints=linspace(-1,1,n+1); 
%[x_npoints,wgl] = legendre_gauss_lobatto(n+1);
k = 50; % equispaced points
xg = linspace(-1,1,k); 

fx_npoints = f(x_npoints); 
df=@(x) -pi/2*sin(pi/2*x);

[L_matrix,dL_matrix]=langrange_function(xg,x_npoints);


f_N =fx_npoints*L_matrix;


 

%figure(1)
plot(xg, f(xg), 'b--', 'LineWidth', 2); hold on;
plot(xg, f_N, 'r-', 'LineWidth', 2);
legend('Analytic', 'Numerical')
xlabel('x')
ylabel('f(x)')
%plot(f_N)







