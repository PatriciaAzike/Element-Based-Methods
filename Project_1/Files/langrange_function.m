%%% Computes Lagrange basis matrix, as well as the derivative matrix



function [L_matrix,dL_matrix]=langrange_function(x,x_npoints)
L_matrix=zeros(length(x_npoints),length(x));
dL_matrix=zeros(length(x_npoints),length(x));
N=length(x_npoints)-1;
for i = 0:N
    xi = x_npoints(i+1); 
    L = 1;
  
    for j = 0:N 
       xj=x_npoints(j+1);
       prod=ones(1,length(x));
      
       if j~=i  
            L = L.*(x-x_npoints(j+1))./(xi-x_npoints(j+1));
            for k=0:N
                xk=x_npoints(k+1);
                if (k~=i && k~=j)
                    prod=prod.*(x-xk)./(xi-xk);
                end
                
            end
            dL_matrix(i+1,:) = dL_matrix(i+1,:)+prod./(xi-xj);
           
                  
       end
       
      
       L_matrix(i+1,:)=L;  
               
    end
    
    
    
end 

end
