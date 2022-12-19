% This code create the global mass and differentiation matrices for the CG
% method using the Direct Stiffness (DSS) operator. We take periodicity into
% account because of the boundary conditions



function [Mmatrix, Dmatrix] = Matrix_DSS(Me,De,intma,periodicity,ngl,Ne,npoin)
    Mmatrix = zeros(npoin,npoin); % global mass matrix for CG
    Dmatrix = zeros(npoin,npoin); % global differentiation matrix for CG
     for e = 1:Ne % go over all elements
         for j = 1:ngl % go over all columns of M(e)
             Z = intma(j,e); % point to global gridpoint number
             J = periodicity(Z); %periodicity
             
             for i = 1:ngl % go over all rows of M(e)
                 Y = intma(i,e); % point to global gridpoint number
                 I = periodicity(Y);
                 
                 Mmatrix(I,J) = Mmatrix(I,J) + Me(i,j);
                 Mmatrix(end,end) = 1;
                 Dmatrix(I,J) = Dmatrix(I,J) + De(i,j);
             end
             
         end
            
        
    end
    
end




