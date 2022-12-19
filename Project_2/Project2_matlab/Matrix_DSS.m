function [Mmatrix, Dmatrix] = Matrix_DSS(Me,De,intma,periodicity,ngl,Ne,npoin)
    Mmatrix = zeros(npoin,npoin); %for CG
    Dmatrix = zeros(npoin,npoin);
     for e = 1:Ne
         for j = 1:ngl
             Z = intma(j,e);
             J = periodicity(Z);
             
             for i = 1:ngl
                 Y = intma(i,e);
                 I = periodicity(Y);
                 Mmatrix(I,J) = Mmatrix(I,J) + Me(i,j);
                 Mmatrix(end,end) = 1;
                 Dmatrix(I,J) = Dmatrix(I,J) + De(i,j);
             end
             
         end
            
        
    end
    
end




