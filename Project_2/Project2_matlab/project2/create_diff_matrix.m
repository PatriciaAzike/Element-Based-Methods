% This code constructs the 1D element differentiation matrix for both CG
% and DG methods



function De = create_diff_matrix(Ne,ngl,nq,wnq,psi,dpsi)
    
    for e = 1:Ne % go over all elements
        De =  zeros(ngl,ngl);
        for k = 1:nq % go over quadrature points
            for j = 1:ngl % loop over columns
                for i = 1:ngl % loop over rows
                    De(i,j) = De(i,j) + wnq(k) * psi(i,k)*dpsi(j,k);
                end
            end
            
        end
    end


end