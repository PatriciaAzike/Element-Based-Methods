%clc
close all
clear all
format long

l1_norm_LGL = zeros(64,1);
l2_norm_LGL = zeros(64,1);
l1_norm_LG = zeros(64,1);
l2_norm_LG = zeros(64,1);

N = (1:64)';

for ipt = 1:2
    
    for nop = 1:64
        
        noq = nop;
        ngl = nop+1;
        nq = noq+1;
        
        if(ipt == 1)
            [xgl,wgl] = legendre_gauss_lobatto(ngl);


            [xnq,wnq] = legendre_gauss_lobatto(nq);

            [error_L1,error_L2]=numerics2(xgl,xnq,wnq,ngl,nq);

            l1_norm_LGL(nop) = error_L1;
            l2_norm_LGL(nop) = error_L2;
        else
            [xgl,wgl] = legendre_gauss(ngl);


            [xnq,wnq] = legendre_gauss(nq);

            [error_L1,error_L2]=numerics2(xgl,xnq,wnq,ngl,nq);

            l1_norm_LG(nop) = error_L1;
            l2_norm_LG(nop) = error_L2;
        end
    end
    
end

figure(1)
semilogy(N,l1_norm_LGL,'r', 'LineWidth',2); hold on;
semilogy(N,l2_norm_LGL,'b--', 'LineWidth',2);hold on;
xlabel('N')
ylabel('Error Norm')
legend('L_1', 'L_2', location='NE')
title('Error Norms of Interpolation using LGL Points')
%axis tight

figure(2)
semilogy(N,l1_norm_LG,'r', 'LineWidth',2); hold on;
semilogy(N,l2_norm_LG,'b--', 'LineWidth',2);hold on;
%semilogy(1:1:64,error_Linf,'g', 'LineWidth',2);
xlabel('N')
ylabel('Error Norm')
legend('L_1', 'L_2', location='NE')
title('Error Norms of Interpolation using LG Points')
%axis tight


% %%%%%%%%%%%%%%%% A3: Legendre points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% 
% for i=1:64 % order of polynomial
%       [xgl,wgl] = legendre_gauss(i+1);
%       [error_L1(i),error_L2(i), error_Linf(i)]=numerics(xgl,wgl,i);
% end
% semilogy(1:1:64,error_L1,'r', 'LineWidth',2); hold on;
% semilogy(1:1:64,error_L2,'b', 'LineWidth',2);hold on;
% semilogy(1:1:64,error_Linf,'g', 'LineWidth',2);
% xlabel('N')
% ylabel('Error Norm')
% legend('L_1', 'L_2', 'L-inf', location='NE')
% title('Error Norms of Interpolation using LG Points')
% axis tight

function [error_L1,error_L2]=numerics2(xgl,xnq,wnq,P,Q)

    f=@(x)(cos(pi/2.*(x)));
    
    Int_exact = 4/pi;
    
    [L,dL] = lagrange_basis2(xnq, xgl, P,Q);
    
    fx = f(xgl)';
    
    f_sum = L'*fx;
    
    Int = 0;
    
    for k = 1:Q
        
        Int = Int + wnq(k)*f_sum(k);
    end
 
    error_L1= abs(Int_exact - Int)/abs(Int_exact);
    

    error_L2=sqrt((Int_exact - Int).^2)/sqrt(Int_exact).^2;

end
    

%%%%%%%%%%%%%%%%%%%Defining the function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
function [L,dL] = lagrange_basis2(x, x_npoints, P,Q)
    
    L = zeros(P,Q);
    dL = zeros(P,Q);
    
    for l = 1:Q

        for i = 1:P
            Li = 1;
            dLi = 0;
            for j = 1:P
               if j~=i  
                   
                    Li = Li.*(x(l)-x_npoints(j))/(x_npoints(i)-x_npoints(j));
                    
                    for k=1:P
                        %xk=x_npoints(k);
                        if (k~=i && k~=j)
                            prod=prod.*(x(l)-x_npoints(k))/(x_npoints(i)-x_npoints(k));
                        end
                    end
                     
                     dLi=dLi+prod/(x_npoints(i)-x_npoints(j));

               end

            end
            L(i,l) = Li;
            dL(i,l) = dLi;
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