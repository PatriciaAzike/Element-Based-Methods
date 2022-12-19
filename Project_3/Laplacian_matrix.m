% This code constructs the 1D element Laplacian matrix


function Le = Laplacian_matrix(intma,coord,Ne,ngl,nq,wnq,dpsi)
    
     for e = 1:Ne % go over all elements
        Le =  zeros(ngl,ngl);
        xe =  zeros(ngl,1); 
        for i = 1:ngl
            I = intma(i,e);
            xe(i) = coord(I); % coordinates of the element
        end
        
        dxe = xe(ngl) - xe(1);  % element size
       
        for k = 1:nq % go over quadrature points
            for j = 1:ngl % loop over columns
                for i = 1:ngl % loop over rows
                    Le(i,j) = Le(i,j) + (2/dxe)*(wnq(k) * dpsi(i,k)*dpsi(j,k));
                end
            end
            
        end
    end


end

