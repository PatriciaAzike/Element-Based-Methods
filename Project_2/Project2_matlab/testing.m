clear all; 
close all;
clc

tic

%Input Data
%---------------------------------------------------------------------%
u = 2; % veolcity
time_final=1; %final time 

Ne=4; %Number of Elements
N=4;    %Interpolation Order
integration_points=2; %=1 for LGL and =2 for LG - LGL is default, but feel free to experiment
integration_type=2; %=1 is inexact and =2 is exact - used to compute how many quadrature points we need
space_method_type='cg'; %CG or DG - assumes that the code can switch between CG and DG
flux_type=2; %1=centered flux and 2=upwind - default is upwind, but feel free to experiment

Courant_max=0.25; %dt controlled by Courant_max
%---------------------------------------------------------------------%

%Store Constants
ngl=N + 1;
npoin_cg=N*Ne + 1; %number of points for CG 
npoin_dg=ngl*Ne; %number of points for DG

%Compute Interpolation Points
[xgl,wgl]=legendre_gauss_lobatto(ngl); %LGL interpolation points is default

%Compute Integration Points
if (integration_points == 1) %LGL integration points
    integration_text='LGL';
    if (integration_type == 1) %inexact
        noq=N;
    elseif (integration_type == 2) %exact
        noq=N+1;
    end
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss_lobatto(nq);
elseif (integration_points == 2) %LG integration points
    integration_text='LG';
    noq=N;
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss(nq);
end


% ---> Compute Lagrange Polynomial and derivatives
% Use the lagrange_basis routine you created in Project 1
% which returns values of basis functions psi at quadrature points ksi_k
% and values of derivative of basis functions dpsi at quadrature points
% both psi and dpsi are arrays of size ngl x nq, where ngl = N+1 is the 
% number of nodal points and nq is the number of quadrature points
% Suggested format:
[psi,dpsi] = lagrange_basis(ngl,nq,xgl,xnq);


%Create Grid - different arrays for CG and DG
[coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg]=create_grid(ngl,Ne,npoin_cg,npoin_dg,xgl);

if strcmp(space_method_type,'cg')
    npoin=npoin_cg;
    coord=coord_cg;
    intma=intma_cg;
    periodicity=periodicity_cg;
elseif strcmp(space_method_type,'dg')
    npoin=npoin_dg;
    coord=coord_dg;
    intma=intma_dg;
    periodicity=periodicity_dg;
end

%Choose time-step such that it does not violate Courant_max, but also
%divides evenly into the final time
dx=coord(2)-coord(1);
dt=Courant_max*dx/u;
ntime=round(time_final/dt);
dt=time_final/ntime;
Courant=u*dt/dx;
disp(['Courant = ',num2str(Courant),' dt = ',num2str(dt),' ntime = ',num2str(ntime),' time_final = ',num2str(time_final)])


% ---> Create Local/Element Mass and Differentiation Matrices
%
% Create function which computes element mass matrix Me (ngl x ngl size) 
Me = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,psi);
%
% Create function which computes element differentiation matrix De (ngl x ngl size) 
De = create_diff_matrix(Ne,ngl,nq,wnq,psi,dpsi);

% Form Global Mass and Differentiation Matrices (both npoin x npoin size, or N_p x N_p)

[Mmatrix, Dmatrix] = Matrix_DSS(Me,De,intma,periodicity,ngl,Ne,npoin_cg); % working on cg for now

% ---> Create right-hand-size matrix Rmatrix (npoin x npoin size)
% Derive Rmatrix from the matrix form of governing equation 
% using mass, differentiation and flux matrices
% You need to modify the equation so you have
% dq/dt = Rmatrix*q
%
Rmatrix = -inv(Mmatrix)*Dmatrix*u;

%Compute Initial Condition
time=0;
qe = initial_condition(coord,npoin);

%Initialize the solution with the initial solution qe
q=qe;

%Time Integration - time-loop is hidden inside the time-integration routine
[q,time] = time_integration(q,Rmatrix,periodicity,time,ntime,dt);


% ---> Compute Norm
% l2_norm= ...


%Plot Solution
    h=figure;
    figure(h);
    plot_handle=plot(coord,q,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(coord,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    title_text=[space_method_type,': Ne = ' num2str(Ne) ', N = ' num2str(N) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);






function [psi,dpsi] = lagrange_basis(ngl,nq,xgl,xnq)
    
    psi=zeros(ngl,nq); % lagrange basis matrix
    dpsi=zeros(ngl,nq); % derivative matrix
    N=ngl-1;
    for i = 0:N
        xi = xgl(i+1); 
        L = 1;

        for j = 0:N 
           xj=xgl(j+1);
           prod=ones(1,nq);

           if j~=i  
                L = L.*(xnq-xj)./(xi-xj);
                for k=0:N
                    xk=xgl(k+1);
                    if (k~=i && k~=j)
                        prod=prod.*(xnq-xk)./(xi-xk);
                    end

                end
                dpsi(i+1,:) = dpsi(i+1,:)+prod./(xi-xj);


           end
                psi(i+1,:)=L;  

        end
    
    
    
    end 

end





function Me = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,psi)
    %dxe = coord(2) - coord(1);
    %dxe
    xe =  zeros(ngl,1); 
    Me =  zeros(ngl,ngl,Ne);
    
    
    for e = 1:Ne
        
        for i = 1:ngl
            I = intma(i,e);
            xe(i) = coord(I); % coordinates of the element
        end
        
        dxe = xe(ngl) - xe(1);  % element size
       
        
        for k = 1:nq
            for j = 1:ngl
                for i = 1:ngl
                    
                    Me(i,j,e) = Me(i,j,e)+ 0.5*dxe*wnq(k) * psi(i,k)*psi(j,k);
                end
            end
            
        end
    end
end


function De = create_diff_matrix(Ne,ngl,nq,wnq,psi,dpsi)
    De =  zeros(ngl,ngl,Ne);
    for e = 1:Ne
        for k = 1:nq
            for j = 1:ngl
                for i = 1:ngl
                    De(i,j,e) = De(i,j,e) + wnq(k) * psi(i,k)*dpsi(j,k);
                end
            end
            
        end
    end


end


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
                 Mmatrix(I,J) = Me(i,j,e);
                 Mmatrix(end,end) = 1;
                 Dmatrix(Y,Z) = Dmatrix(Y,Z) + De(i,j,e);
             end
         end
            
        
    end
    
end




