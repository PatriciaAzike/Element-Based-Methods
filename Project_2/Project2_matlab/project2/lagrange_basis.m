% This code constructs the lagrange basis polynomial as well as the
% derivative matrices


function [psi,dpsi] = lagrange_basis(ngl,nq,xgl,xnq)
    
    psi=zeros(ngl,nq); % lagrange basis matrix
    dpsi=zeros(ngl,nq); % derivative matrix
    N=ngl-1;
    for i = 0:N
        xi = xgl(i+1); 
        L = 1;

        for j = 0:N 
           xj=xgl(j+1);
           prod=ones(1,nq);

           if j~=i  
                L = L.*(xnq-xj)./(xi-xj);
                for k=0:N
                    xk=xgl(k+1);
                    if (k~=i && k~=j)
                        prod=prod.*(xnq-xk)./(xi-xk);
                    end

                end
                dpsi(i+1,:) = dpsi(i+1,:)+prod./(xi-xj);


           end
                psi(i+1,:)=L;  

        end
    
    
    
    end 

end

