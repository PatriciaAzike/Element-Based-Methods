function De = create_diff_matrix(Ne,ngl,nq,wnq,psi,dpsi)
    
    for e = 1:Ne
        De =  zeros(ngl,ngl);
        for k = 1:nq
            for j = 1:ngl
                for i = 1:ngl
                    De(i,j) = De(i,j) + wnq(k) * psi(i,k)*dpsi(j,k);
                end
            end
            
        end
    end


end