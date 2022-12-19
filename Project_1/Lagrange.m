x_points = linspace(-1,1,50);
f_true = cos(pi/2.*x_k);

N = size(x_k,2);
L = ones(N,N);
for i





function f = Lagrange_basis(x,xpoints)
    N = length(xpoints);
    
    for l = 1:N
        x_l = xpoints(l);
        for i = 1:N
            x_i = x(i);
            L(i,l) = 1;
            for j = 1:N
                x_j = x(j);
                if (j~=i)
                    L(i,l) = (x_l-x_j)/(x_i-x_j);
                end
            end   
        end    
        
        

    end
end    
