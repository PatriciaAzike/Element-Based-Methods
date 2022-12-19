clc
clear all
close all
% 
% 
% 
% % %%%%%%%%%%%%% A1: Equally spaced points %%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% for i=1:64 % order of polynomial
%       xgl = linspace(-1,1,i+1); %% equi-spaced points
%       [error_L1(i),error_L2(i), error_Linf(i)]=numerics(xgl,i);
%      
% end
% semilogy(1:1:64,error_L1,'r', 'LineWidth',2); hold on;
% semilogy(1:1:64,error_L2,'b', 'LineWidth',2);hold on;
% semilogy(1:1:64,error_Linf,'g', 'LineWidth',2);
% xlabel('N')
% ylabel('Error Norm')
% legend('L_1', 'L_2', 'L-inf', location='North')
% title('Error Norms of Derivative Interpolation using Equi-spaced Points')
% axis tight
% 
% 
% %%%%%%%%%%%%%%% A2: Lobatto points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% 
% for i=1:64 % order of polynomial
%       [xgl,wgl] = legendre_gauss_lobatto(i+1);
%       [error_L1(i),error_L2(i), error_Linf(i)]=numerics(xgl,i);
% end
% semilogy(1:1:64,error_L1,'r', 'LineWidth',2); hold on;
% semilogy(1:1:64,error_L2,'b', 'LineWidth',2);hold on;
% semilogy(1:1:64,error_Linf,'g', 'LineWidth',2);
% xlabel('N')
% ylabel('Error Norm')
% legend('L_1', 'L_2', 'L-inf', location='NE')
% title('Error Norms of Derivative Interpolation using LGL Points')
% axis tight
% 
% 
% %%%%%%%%%%%%%%%% A3: Legendre points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% 
% for i=1:64 % order of polynomial
%       [xgl,wgl] = legendre_gauss(i+1);
%       [error_L1(i),error_L2(i), error_Linf(i)]=numerics(xgl,i);
% end
% semilogy(1:1:64,error_L1,'r', 'LineWidth',2); hold on;
% semilogy(1:1:64,error_L2,'b', 'LineWidth',2);hold on;
% semilogy(1:1:64,error_Linf,'g', 'LineWidth',2);
% xlabel('N')
% ylabel('Error Norm')
% legend('L_1', 'L_2', 'L-inf', location='NE')
% title('Error Norms of Derivative Interpolation using LG Points')
% axis tight
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [error_L1,error_L2, error_Linf]=numerics(x_npoints,n)
%     f=@(x)(cos(pi/2.*(x)));
%     df=@(x)(-pi/2*sin(pi/2.*x));
%     fx_npoints = f(x_npoints); 
%     k = 50; % 50 evenly spaced points
%     xg = linspace(-1,1,k); 
%     f_N = zeros(size(xg)); 
% 
%     for i = 0:n      
%          dL = DLB(xg, x_npoints,i); 
%          f_N = f_N + dL * fx_npoints(i+1); 
%     end 
% % plot(xg, f(xg), 'b-', 'LineWidth', 2); hold on;
% % plot(xg, f_N, 'r-', 'LineWidth', 2); 
% % plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)
% % legend("Analytical", 'Numerical', 'Nodes')
%     error_L1=sum(abs(f_N-df(xg)))/sum(abs(df(xg)));
%     error_L2=sqrt(sum((f_N-df(xg)).^2))/sqrt(sum((df(xg)).^2));
%     error_Linf=max(abs(f_N-df(xg)))/max(abs(df(xg)));
% end
% 
% 
% function dL=DLB(x,x_npoints,i)
%     dL=0;
%     n=length(x_npoints)-1;
%     xi=x_npoints(i+1);
%     for j=0:n
%         xj=x_npoints(j+1);
%         prod=1;
%         if j~=i
%             for k=0:n
%                 xk=x_npoints(k+1);
%                 if (k~=i && k~=j)
%                     prod=prod.*(x-xk)/(xi-xk);
%                 end
%              end
%              dL=dL+prod/(xi-xj);
%          end
%      end
%     
% end
% 
% 
% 

figure(1)
i = 70;
xgl = linspace(-1,1,i+1); %% equi-spaced points
[I,L]=numerics(xgl,i);

function [f_N,L]=numerics(x_npoints,n)
    f=@(x)(cos(pi/2.*(x)));
    fx_npoints = f(x_npoints); 

    k = 50; % 50 evenly spaced points
    xg = linspace(-1,1,k); 
    f_N = zeros(size(xg)); 
    for i = 0:n      
        L = lagrange_basis(xg, x_npoints, i); % Lagrange basis
        f_N = f_N + fx_npoints(i+1)*L; 
    end 

    plot(xg, f(xg), 'b-', 'LineWidth', 4); hold on;
    plot(xg, f_N, 'r-', 'LineWidth', 2); 
    %plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)
    legend("Analytical", 'Numerical')
    xlabel('x')
    ylabel('f(x)')
    axis tight
end

%% %%%%%%%%%%%%%%%%%%%%%Lagrange Basis Function%%%%%%%%%%%%%%%%%%%%%
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