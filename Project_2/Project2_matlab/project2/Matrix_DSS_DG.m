% This code create the global mass and differentiation matrices for the DG
% method using the Direct Stiffness (DSS) operator. Note that periodicty is
% not considered here because the Flux matrix accounts for it



function [Mmatrix, Dmatrix] = Matrix_DSS_DG(Me,De,intma,ngl,Ne,npoin)
    Mmatrix = zeros(npoin,npoin); % global mass matrix for DG
    Dmatrix = zeros(npoin,npoin); % global differentiation matrix for DG
    Dhat = De.'; % Dhat is simply the transpose of D(e) from the CG method
    
     for e = 1:Ne % go over all elements
         for j = 1:ngl % go over all columns of M(e)
             J = intma(j,e); % point to the global gridpoint number
             
             for i = 1:ngl % go over all rows of M(e)
                 I = intma(i,e); %point to the global gridpoint number
                 
                 Mmatrix(I,J) = Mmatrix(I,J) + Me(i,j);
                 Dmatrix(I,J) = Dmatrix(I,J) + Dhat(i,j);
                
             end
 
         end
            
        
    end
    
end




