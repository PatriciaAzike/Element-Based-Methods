clc
clear all
close all

f=@(x) cos(pi/2 * x);
n=64;
x_npoints=linspace(-1,1,n+1)'; 
k = 50; % equispaced points
xg = linspace(-1,1,k)'; 
%f_N = zeros(size(xg)); 
fx_npoints = f(x_npoints); 

D = DDL(xg,x_npoints);
f_N = D*f(x_npoints);

 
df=@(x) -pi/2*sin(pi/2*x);
%figure(1)
plot(xg, df(xg), 'b-', 'LineWidth', 2); hold on;
plot(xg, f_N, 'r-', 'LineWidth', 2);
%plot(f_N)


function D = DDL(x,x_npoints)
%Q=50;
    
    Q = length(x);
    n=length(x_npoints);
    D = zeros(Q,n);
%Q=n;
    for i = 1:n
        xi=x_npoints(i);
        for j=1:n
            xj=x_npoints(j);
            prod=ones(Q,1);
            if j~=i
                for k=1:n
                    xk=x_npoints(k);
                    if ((k~=i) && (k~=j))
                        prod=prod.*(x-xk)./(xi-xk);
                    end
                end
                D(:,i) = D(:,i) + prod./(xi-xj);
            end
        end
    end
end



