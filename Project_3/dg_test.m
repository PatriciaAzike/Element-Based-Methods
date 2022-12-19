%%%% Driver code for DG method: takes in the number of elements, Ne, the
%%%% order of the polynomial and the integration type, and returns the L_2
%%%% norm

clc
clear all
close all
format short



Ne=4; %Number of Elements
N=1;    %Interpolation Order
integration_points=1; %=1 for LGL and =2 for LG - LGL is default, but feel free to experiment
integration_type=2; %=1 is inexact and =2 is exact - used to compute how many quadrature points we need
space_method_type='dg'; %CG or DG - assumes that the code can switch between CG and DG
flux_type=1; %1=centered flux 

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


% Create function which computes element mass matrix Me (ngl x ngl size) 
Me = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,psi);

% Create function which computes element differentiation matrix De (ngl x ngl size) 
De = create_diff_matrix(Ne,ngl,nq,wnq,psi,dpsi);
%
% Create function which computes element Laplacian matrix De (ngl x ngl size) 
Le = Laplacian_matrix(intma,coord,Ne,ngl,nq,wnq,dpsi);
% Form Global Mass and Differentiation Matrices (both npoin x npoin size, or N_p x N_p)
[Mmatrix,Dmatrix, Lmatrix] = Matrix_DSS(Me,De,Le,intma,ngl,Ne,npoin); 


%Apply BCs

if (flux_type == 1)
        Fmatrix = Fmatrix_centered_flux(intma,Ne,npoin,ngl);
elseif (flux_type == 2)
        Fmatrix = Fmatrix_upwind_flux(intma,Ne,npoin,ngl);
end
% ---> Create right-hand-size matrix Rmatrix (npoin x npoin size)

F_q = Fmatrix;
F_q(1,end) = 0; F_q(end,1) = 0;
F_q(end,end) = F_q(end,end) + 0.5;

F_Q = Fmatrix;
F_Q(1,1) = F_Q(1,1) - 0.5;
F_Q(1,end) = 0; F_Q(end,1) = 0;
F_Q(end,end) = F_Q(end,end);

B_Q = zeros(npoin,1);
B_Q(end) = - 0.5*pi;

D_DG_q = F_q - Dmatrix;
D_DG_Q = F_Q - Dmatrix;

L_DG = D_DG_Q * (Mmatrix\D_DG_q);
f = f_function(coord,npoin);

% Compute the numerical solution
q = (L_DG + Mmatrix)\(Mmatrix * f - B_Q);



%Compute exact solution
qe = q_exact(coord,npoin);


% ---> Compute Norm
l2_norm = sqrt(sum((q - qe).^2)/sum((qe).^2));

%Plot Solution
    h=figure;
    figure(h);
    plot_handle=plot(coord,q,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(coord,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x)','FontSize',18);
    legend('Numerical', 'Exact',location='NW')

    title_text=[space_method_type,': Ne = ' num2str(Ne), ', N = ' num2str(N), ', Q = ' num2str(noq), ', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    axis tight











