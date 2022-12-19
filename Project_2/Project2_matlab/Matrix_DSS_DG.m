function [Mmatrix, Dmatrix] = Matrix_DSS_DG(Me,De,intma,ngl,Ne,npoin)
    Mmatrix = zeros(npoin,npoin); 
    Dmatrix = zeros(npoin,npoin);
    Dhat = De.';
    
     for e = 1:Ne
         for j = 1:ngl
             J = intma(j,e);
             
             for i = 1:ngl
                 I = intma(i,e);
                 Mmatrix(I,J) = Mmatrix(I,J) + Me(i,j);
                 Dmatrix(I,J) = Dmatrix(I,J) + Dhat(i,j);
                
             end
 
         end
            
        
    end
    
end




