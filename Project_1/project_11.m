clc
close all
clear all
format long

N=64;
x_npoints=linspace(-1,1,N+1);
[ error_L1,error_L2,error_Linf]=numerics(x_npoints);
% for i=1:N % order of polynomial
%       xgl = linspace(-1,1,i+1); %% equi-spaced points
%       [error_L1_EQ(i),error_L2_EQ(i), error_Linf_EQ(i)]=numerics(xgl,i);
%       
%       [xgl,wgl] = legendre_gauss_lobatto(i+1); %% LGL points
%       [error_L1_LGL(i),error_L2_LGL(i), error_Linf_LGL(i)]=numerics(xgl,i);
%       
%       [xgl,wgl] = legendre_gauss(i+1); %% LG points
%       [error_L1_LG(i),error_L2_LG(i), error_Linf_LG(i)]=numerics(xgl,i);
%       
% end

% 
% %%Plotting L1 norm
% figure(1)
% semilogy(1:1:N,error_L1_EQ,'r', 'LineWidth',2); hold on;
% semilogy(1:1:N,error_L1_LGL,'b', 'LineWidth',2);hold on;
% semilogy(1:1:N,error_L1_LG,'m', 'LineWidth',2);
% xlabel('N')
% ylabel('L_1 Norm')
% legend('Equi-spaced', 'LGL', 'LG', location='N')
% title('L_1 Norms of Function Interpolation using various Points')
% axis tight
% 
% 
% %%Plotting L2 norm
% figure(2)
% semilogy(1:1:N,error_L2_EQ,'r', 'LineWidth',2); hold on;
% semilogy(1:1:N,error_L2_LGL,'b', 'LineWidth',2);hold on;
% semilogy(1:1:N,error_L2_LG,'m', 'LineWidth',2);
% xlabel('N')
% ylabel('L_2 Norm')
% legend('Equi-spaced', 'LGL', 'LG', location='N')
% title('L_2 Norms of Function Interpolation using various Points')
% axis tight
% 
% 
% %%Plotting L-inf norm
% figure(3)
% semilogy(1:1:N,error_Linf_EQ,'r', 'LineWidth',2); hold on;
% semilogy(1:1:N,error_Linf_LGL,'b', 'LineWidth',2);hold on;
% semilogy(1:1:N,error_Linf_LG,'m', 'LineWidth',2);
% xlabel('N')
% ylabel('L_{\infty} Norm')
% legend('Equi-spaced', 'LGL', 'LG', location='N')
% title('L_{\infty} Norms of Function Interpolation using various Points')
% axis tight
% 

%%%%%%%%%%%%%%%%%%%%%%Defining the error function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ error_L1,error_L2,error_Linf]=numerics(x_npoints)
    f=@(x)(cos(pi/2.*(x)));
    P=f(x_npoints);
    k = 50; % 50 evenly spaced points
    xg = linspace(-1,1,k); 
    P=f(x_npoints)
    L_matrix=langrange_basis(xg,x_npoints);
    
    f_N =L_matrix*f(x_npoints)';
 

    error_L1=sum(abs(f_N-f(xg)))/sum(abs(f(xg)));
    error_L2=sqrt(sum((f_N-f(xg)).^2))/sqrt(sum((f(xg)).^2));
    error_Linf=max(abs(f_N-f(xg)))/max(abs(f(xg)));
end

%%%%%%%%%%%%%%%%%%%%%%%%Lagrange Basis Function%%%%%%%%%%%%%%%%%%%%%
function L_matrix=langrange_basis(x,x_npoints)
L_matrix=zeros(length(x),length(x_npoints));
N=length(x_npoints)-1;
for i = 0:N
    xi = x_npoints(i+1); 
    L = 1;
    for j = 0:N 
       if j~=i  
            L = L.*(x-x_npoints(j+1))/(xi-x_npoints(j+1)); 
       end
    end
    L_matrix(:,i+1)=L;
       
end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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