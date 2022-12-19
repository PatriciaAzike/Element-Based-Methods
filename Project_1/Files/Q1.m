clc
close all
clear all
format long

N=64;

error_L1_EQ=zeros(64,1);
error_L2_EQ=zeros(64,1);
error_Linf_EQ=zeros(64,1);


error_L1_LGL=zeros(64,1);
error_L2_LGL=zeros(64,1);
error_Linf_LGL=zeros(64,1);

error_L1_LG=zeros(64,1);
error_L2_LG=zeros(64,1);
error_Linf_LG=zeros(64,1);




for i=1:N % order of polynomial
      xgl = linspace(-1,1,i+1); %% equi-spaced points
      [error_L1_EQ(i),error_L2_EQ(i), error_Linf_EQ(i)]=numerics(xgl);
      
      [xgl,wgl] = legendre_gauss_lobatto(i+1); %% LGL points
      [error_L1_LGL(i),error_L2_LGL(i), error_Linf_LGL(i)]=numerics(xgl);
      
      [xgl,wgl] = legendre_gauss(i+1); %% LG points
      [error_L1_LG(i),error_L2_LG(i), error_Linf_LG(i)]=numerics(xgl);
      
end


%%Plotting L1 norm
figure(1)
semilogy(1:1:N,error_L1_EQ,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_L1_LGL,'b', 'LineWidth',2);hold on;
semilogy(1:1:N,error_L1_LG,'m', 'LineWidth',2);
xlabel('N')
ylabel('L_1 Norm')
legend('Equi-spaced', 'LGL', 'LG', location='N')
title('L_1 Norms of Function Interpolation using various Points')
%axis tight


%%Plotting L2 norm
figure(2)
semilogy(1:1:N,error_L2_EQ,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_L2_LGL,'b', 'LineWidth',2);hold on;
semilogy(1:1:N,error_L2_LG,'m', 'LineWidth',2);
xlabel('N')
ylabel('L_2 Norm')
legend('Equi-spaced', 'LGL', 'LG', location='N')
title('L_2 Norms of Function Interpolation using various Points')
%axis tight


%%Plotting L-inf norm
figure(3)
semilogy(1:1:N,error_Linf_EQ,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_Linf_LGL,'b', 'LineWidth',2);hold on;
semilogy(1:1:N,error_Linf_LG,'m', 'LineWidth',2);
xlabel('N')
ylabel('L_{\infty} Norm')
legend('Equi-spaced', 'LGL', 'LG', location='N')
title('L_{\infty} Norms of Function Interpolation using various Points')
%axis tight


%%%%%%%%%%%%%%%%%%%%%%Defining the error function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ error_L1,error_L2,error_Linf]=numerics(x_npoints)
    f=@(x)(cos(pi/2.*(x)));
 
    k = 50; % 50 evenly spaced points
    xg = linspace(-1,1,k); 
    [L_matrix,dL_matrix]=langrange_function(xg,x_npoints);
    
    f_N =f(x_npoints)*L_matrix;
 

    error_L1=sum(abs(f_N-f(xg)))/sum(abs(f(xg)));
    error_L2=sqrt(sum((f_N-f(xg)).^2))/sqrt(sum((f(xg)).^2));
    error_Linf=max(abs(f_N-f(xg)))/max(abs(f(xg)));
end



