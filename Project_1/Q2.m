clc
clear all
close all



N=64;
for i=1:N % order of polynomial
      xgl = linspace(-1,1,i+1); %% equi-spaced points
      [error_L1_EQ(i),error_L2_EQ(i), error_Linf_EQ(i)]=numerics(xgl,i);
      
      [xgl,wgl] = legendre_gauss_lobatto(i+1); %% LGL points
      [error_L1_LGL(i),error_L2_LGL(i), error_Linf_LGL(i)]=numerics(xgl,i);
      
      [xgl,wgl] = legendre_gauss(i+1); %% LG points
      [error_L1_LG(i),error_L2_LG(i), error_Linf_LG(i)]=numerics(xgl,i);
  
end

%% Plotting L1 norm
figure(1)
semilogy(1:1:N,error_L1_LG,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_L1_EQ,'b', 'LineWidth',2);hold on;
semilogy(1:1:N,error_Linf_LG,'m', 'LineWidth',2);
xlabel('N')
ylabel('L_1 Norm')
legend('Equi-spaced', 'LGL', 'LG', location='N')
title('L_1 Norms of Derivative using various Points')
%axis tight


%% Plotting L2 norm
% figure(2)
% semilogy(1:1:N,error_L2_EQ,'r', 'LineWidth',2); hold on;
% semilogy(1:1:N,error_L2_LGL,'b', 'LineWidth',2);hold on;
% semilogy(1:1:N,error_L2_LG,'m', 'LineWidth',2);
% xlabel('N')
% ylabel('L_2 Norm')
% legend('Equi-spaced', 'LGL', 'LG', location='N')
% title('L_2 Norms of Derivative using various Points')
%axis tight


%% Plotting L-inf norm
% figure(3)
% semilogy(1:1:N,error_Linf_EQ,'r', 'LineWidth',2); hold on;
% semilogy(1:1:N,error_Linf_LGL,'b', 'LineWidth',2);hold on;
% semilogy(1:1:N,error_Linf_LG,'m', 'LineWidth',2);
% xlabel('N')
% ylabel('L_{\infty} Norm')
% legend('Equi-spaced', 'LGL', 'LG', location='N')
% title('L_{\infty} Norms of Derivative using various Points')
%axis tight



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%Definition of the error function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [error_L1,error_L2, error_Linf]=numerics(x_npoints,n)
    f=@(x)(cos(pi/2.*(x)));
    df=@(x)(-pi/2*sin(pi/2.*x));
    fx_npoints = f(x_npoints); 
    k = 101; % 50 evenly spaced points
    xg = linspace(-1,1,k); 
    f_N = zeros(size(xg)); 

    for i = 0:n      
         dL = DLB(xg, x_npoints,i); 
         f_N = f_N + dL * fx_npoints(i+1); 
    end 
% plot(xg, f(xg), 'b-', 'LineWidth', 2); hold on;
% plot(xg, f_N, 'r-', 'LineWidth', 2); 
% plot(x_npoints, fx_npoints, 'ko', 'LineWidth', 3)
% legend("Analytical", 'Numerical', 'Nodes')
    error_L1=sum(abs(f_N-df(xg)))/sum(abs(df(xg)));
    error_L2=sqrt(sum((f_N-df(xg)).^2))/sqrt(sum((df(xg)).^2));
    error_Linf=max(abs(f_N-df(xg)))/max(abs(df(xg)));
end


function dL=DLB(x,x_npoints,i)
    dL=0;
    n=length(x_npoints)-1;
    xi=x_npoints(i+1);
    for j=0:n
        xj=x_npoints(j+1);
        prod=1;
        if j~=i
            for k=0:n
                xk=x_npoints(k+1);
                if (k~=i && k~=j)
                    prod=prod.*(x-xk)/(xi-xk);
                end
             end
             dL=dL+prod/(xi-xj);
         end
     end
    
end



