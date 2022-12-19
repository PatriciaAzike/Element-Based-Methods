close all
clear all
format long




error_L1_LGL=zeros(64,1);
error_L2_LGL=zeros(64,1);

error_L1_LG=zeros(64,1);
error_L2_LG=zeros(64,1);


N=64;
for i=1:N % order of polynomial
     
      
      [xgl,wgl] = legendre_gauss_lobatto(i+1); %% LGL points
      [error_L1_LGL(i),error_L2_LGL(i)]=numerics2(xgl,xgl,wgl,i+1);
      
      [xgl,wgl] = legendre_gauss(i+1); %% LG points
      [error_L1_LG(i),error_L2_LG(i)]=numerics2(xgl,xgl,wgl,i+1);
  
end


figure(1)
semilogy(1:1:N,error_L1_LGL,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_L2_LGL,'b--', 'LineWidth',2);hold on;
xlabel('N')
ylabel('Error Norm')
legend('L_1', 'L_2', location='NE')
title('Error Norms of Integration using LGL Points')


figure(2)
semilogy(1:1:N,error_L1_LG,'r', 'LineWidth',2); hold on;
semilogy(1:1:N,error_L1_LG,'b--', 'LineWidth',2);hold on;
xlabel('N')
ylabel('Error Norm')
legend('L_1', 'L_2', location='NE')
title('Error Norms of Integration using LG Points')




function [error_L1,error_L2]=numerics2(xgl,xnq,wnq,Q)

    f=@(x)(cos(pi/2.*(x)));
    
    Int_exact = 4/pi;
    
    [L_matrix,dL_matrix]= langrange_function(xnq, xgl);
    
    fx = f(xgl);
    
    f_sum =fx*L_matrix;
    
    Int = 0;
    
    for k = 1:Q
        
        Int = Int + wnq(k)*f_sum(k);
    end
 
    error_L1= abs(Int_exact - Int)/abs(Int_exact);
    

    error_L2=sqrt((Int_exact - Int).^2)/sqrt(Int_exact).^2;

end
    

