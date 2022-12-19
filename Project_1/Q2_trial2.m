clc
clear all
close all

f=@(x)(cos(pi/2.*(x.^2)));
n=63;
x_npoints=linspace(-1,1,n+1);
k = 50; % equispaced points
xg = linspace(-1,1,k); 
f_N = zeros(size(xg)); 
fx_npoints = f(x_npoints); 
 for i = 0:n      
     dL = deriv(xg, x_npoints,i); 
     f_N = f_N + fx_npoints(i+1)*dL; 
 end 
 
 
df=@(x)(-pi/2*sin(pi/2.*(x.^2)));
plot(xg, df(xg), 'b-', 'LineWidth', 2); hold on;
plot(xg, f_N, 'r-', 'LineWidth', 2);
%plot(f_N)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = lagrange_basis(x, x_npoints, i)
    n = length(x_npoints)-1; 
    xi = x_npoints(i+1); 
    L = 1; 
    for j = 0:n 
       if j~=i  
            L = L.*(x-x_npoints(j+1))/(xi-x_npoints(j+1)); 
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DL = deriv(x, x_npoints, i)
    n = length(x_npoints)-1; 
    xi = x_npoints(i+1); 
    DL = 0; 
    for k = 0:n 
        xk = x_npoints(k+1);
       if k~=i 
           L = lagrange_basis(x, x_npoints, i);
            DL = DL + L./(xi-xk); 
       end
    end
end