% This code constructs the 1D element mass matrix


function Me = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,psi)
    
    for e = 1:Ne % go over all elements
        Me =  zeros(ngl,ngl);
        xe =  zeros(ngl,1); 
        for i = 1:ngl
            I = intma(i,e);
            xe(i) = coord(I); % coordinates of the element
        end
        
        dxe = xe(ngl) - xe(1);  % element size
       
        
        for k = 1:nq % go over quadrature points
            for j = 1:ngl % loop over columns
                for i = 1:ngl % loop over all rows
                    
                    Me(i,j) = Me(i,j)+ 0.5*dxe*wnq(k) * psi(i,k)*psi(j,k);
                end
            end
            
        end
    end
end
