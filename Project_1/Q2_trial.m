% clc; clear all; close all;
% 
% n = 64;
% x_npoints = linspace(-1,1,n+1);
% 
% f=@(x)(cos(pi/2.*(x.^2)));
% 
% df=@(x)(-pi/2*sin(pi/2.*(x.^2)));
% 
% fx_npoints = f(x_npoints); 
% k = 50; % 50 evenly spaced points
% xg = linspace(-1,1,k); 
% df_N = zeros(size(xg)); 
% 
% for i = 1:n      
%     L = derivative(xg, x_npoints, i); 
%     df_N = df_N + fx_npoints(i)*L; 
% end 
% %d = numerics(x_npoints,n);
% plot(xg, df(xg), 'b-', 'LineWidth', 2); hold on;
% plot(xg, df_N, 'r-', 'LineWidth', 2); 
% ylim([-1500,1000])
% xlim([-1,1])
% function error_num=numerics(x_npoints,n)
% f=@(x)(cos(pi/2.*(x.^2)));
% 
% df=@(x)(-pi/2*sin(pi/2.*(x.^2)));
% 
% fx_npoints = f(x_npoints); 
% k = 50; % 50 evenly spaced points
% xg = linspace(-1,1,k); 
% df_N = zeros(size(xg)); 
% 
% for i = 0:n      
%     L = derivative(xg, x_npoints, i); 
%     df_N = df_N + fx_npoints(i)*L; 
% end 
% error_num=sum(abs(df_N-df(xg)))/sum(abs(df(xg)));
% end

% for i = 0:n      
%     L = derivative(xg, x_npoints, i); 
%     %f_N = f_N + fx_npoints(i+1)*L; 
% end 
% plot(xg, f(xg), 'b-', 'LineWidth', 2); hold on;
% plot(xg, f_N, 'r-', 'LineWidth', 2); 
% plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)
% legend("Analytical", 'Numerical', 'Nodes')
%error_num=sum(abs(df_N-f(xg)))/sum(abs(df(xg)));
%end

% % plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)




% function dL = derivative(x, x_npoints,i)
%     n = length(x_npoints);
%     tmp1 = 0;
%     for k = 1:n
%         if i==k
%             continue;
%         else
%             tmp2 = 1;
%             for j = 1:n
%                 if j~=i && j~=k
%                     tmp2 = tmp2.* (x - x_npoints(j))/(x_npoints(i)-x_npoints(j));
%                 end
%             end
%         end
%         tmp1 = tmp1 + tmp2*(1/(x_npoints(i) - x_npoints(k)));
%     end
%     dL = tmp1;
% 
% 
% end


clc
close all
clear all
format long



% %%%%%%%%%%%%% A1: Equally spaced points %%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for i=1:64 % order of polynomial
      xgl = linspace(-1,1,i+1); %% equi-spaced points
      e_norm(i)=numerics(xgl,i);
end
semilogy(1:1:64,e_norm,'LineWidth',2);
axis tight



%%%%%%%%%%%%%%%% A2: Lobatto points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

for i=1:64 % order of polynomial
      [xgl,wgl] = legendre_gauss_lobatto(i+1);
      e_norm(i)=numerics(xgl,i);
end
semilogy(1:1:64,e_norm,'LineWidth',2);
axis tight

%%%%%%%%%%%%%%%% A3: Legendre points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

for i=1:64 % order of polynomial
      [xgl,wgl] = legendre_gauss(i+1);
      e_norm(i)=numerics(xgl,i);
end
semilogy(1:1:64,e_norm,'LineWidth',2);
axis tight




%%%%%%%%%%%%%%%%%%%Defining the function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error_num=numerics(x_npoints,n)
f=@(x)(cos(pi/2.*(x)));

fx_npoints = f(x_npoints); 


k = 50; % 50 evenly spaced points
xg = linspace(-1,1,k); 
f_N = zeros(size(xg)); 

for i = 0:n      
    L = lagrange_basis(xg, x_npoints, i); 
    f_N = f_N + fx_npoints(i+1)*L; 
end 
% plot(xg, f(xg), 'b-', 'LineWidth', 2); hold on;
% plot(xg, f_N, 'r-', 'LineWidth', 2); 
% plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)
% legend("Analytical", 'Numerical', 'Nodes')
error_num=sum(abs(f_N-f(xg)))/sum(abs(f(xg)));
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xgl,wgl] = legendre_gauss_lobatto(P)

p=P-1; %Order of the Polynomials
ph=floor( (p+1)/2 );

for i=1:ph
   x=cos( (2*i-1)*pi/(2*p+1) );
   for k=1:20
      [L0,L0_1,L0_2]=legendre_poly(p,x); %Compute Nth order Derivatives of Legendre Polys
      
      %Get new Newton Iteration
      dx=-(1-x^2)*L0_1/(-2*x*L0_1 + (1-x^2)*L0_2);
      x=x+dx;
      if (abs(dx) < 1.0e-20) 
         break
      end
   end
   xgl(p+2-i)=x;
   wgl(p+2-i)=2/(p*(p+1)*L0^2);
end

%Check for Zero Root
if (p+1 ~= 2*ph)
   x=0;
   [L0,L0_1,L0_2]=legendre_poly(p,x);
   xgl(ph+1)=x;
   wgl(ph+1)=2/(p*(p+1)*L0^2);
end
   
%Find remainder of roots via symmetry
for i=1:ph
   xgl(i)=-xgl(p+2-i);
   wgl(i)=+wgl(p+2-i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xgl,wgl] = legendre_gauss(P)

p=P-1; %Order of Polynomials (P is the number of points: p+1)
ph=floor( (p+1)/2 );

for i=1:ph
   x=cos( (2*i-1)*pi/(2*p+1) );
   for k=1:20
      [L0,L0_1,L0_2]=legendre_poly(p+1,x); %Compute the N+1 Legendre Polys
      
      %Get new Newton Iteration
      dx=-L0/L0_1;
      x=x+dx;
      if (abs(dx) < 1.0e-20) 
         break
      end
   end
   xgl(p+2-i)=x;
   wgl(p+2-i)=2/( (1-x^2)*L0_1^2 );
end

%Check for Zero Root
if (p+1 ~= 2*ph)
   x=0;
   [L0,L0_1,L0_2]=legendre_poly(p+1,x);
   xgl(ph+1)=x;
   wgl(ph+1)=2/( (1-x^2)*L0_1^2 );
end
   
%Find remainder of roots via symmetry
for i=1:ph
   xgl(i)=-xgl(p+2-i);
   wgl(i)=+wgl(p+2-i);
end
end
