% This code create the global mass and differentiation matrices for the CG
% method using the Direct Stiffness (DSS) operator. We take periodicity into
% account because of the boundary conditions



function [Mmatrix, Dmatrix, Lmatrix] = Matrix_DSS(Me,De,Le,intma,ngl,Ne,npoin)
    Mmatrix = zeros(npoin,npoin); % global mass matrix 
    Dmatrix = zeros(npoin,npoin); % global differentiation matrix 
    Lmatrix = zeros(npoin,npoin); % global Laplacian matrix 
    Dhat = De.';
    
     for e = 1:Ne % go over all elements
         for j = 1:ngl % go over all columns of M(e), D(e) and L(e)
             J = intma(j,e); % point to global gridpoint number
             
             
             for i = 1:ngl % go over all rows of M(e), D(e) and L(e)
                 I = intma(i,e); % point to global gridpoint number
                 
                 Mmatrix(I,J) = Mmatrix(I,J) + Me(i,j); 
                 Dmatrix(I,J) = Dmatrix(I,J) + Dhat(i,j);
                 Lmatrix(I,J) = Lmatrix(I,J) + Le(i,j);
             end
             
         end
            
        
    end
    
end


