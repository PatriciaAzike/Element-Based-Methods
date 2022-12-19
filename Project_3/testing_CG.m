% This code reads information from the driver code set and
% executes it for the set variables. It also produces the convergence plots
% for both exact and inexact integration

warning off

clear all; 
close all;
format shortE
clc

tic

space_method_type = 'cg';
for i=1:2 %integration_type: 1- inexact, 2- exact
    if i==1
        N=[1 2 4 8 16];
        NeM=[4 8 16 32 64; 4 8 16 32 64; 2 4 8 12 16; 1 2 4 6 8; 1 2 3 4 8]; %Matrix containing Ne
        ng_points=zeros(length(N),size(NeM,2));
        L2_inexact=zeros(length(N),size(NeM,2));
        for j =1:length(N)
            for k=1:size(NeM,2)
                ng_points(j,k)=N(j)*NeM(j,k) + 1;
                [L2_inexact(j,k)]=driver_code(i,N(j),  NeM(j,k), space_method_type);
            end     
        end
         
    else
        N=[1 2 4 8 16];
        NeM=[4 8 16 32 64; 4 8 16 32 64; 2 4 8 12 16; 1 2 4 6 8; 1 2 3 4 8];
        L2_exact=zeros(length(N),size(NeM,2));
        for j =1:length(N)
            for k=1:size(NeM,2)
                [L2_exact(j,k)]=driver_code(i,N(j),  NeM(j,k), space_method_type);
            end   
        end
    end
end



figure(1)% Plotting L_2 norm of inexact and exact integration
slope_inex = zeros(size(ng_points,1)-1,1); % slope of the convergence
Label={'N=1','N=2','N=4','N=8','N=16', 'DG exact'};
color = {'r-*', 'g-*', 'b-*', 'c-*', 'm-*', 'k-o'};

for i=1:5
    
    loglog(ng_points(i,:),L2_inexact(i,:),color{i}, 'LineWidth',2), hold on
    if i<5
        line = polyfit(log10(ng_points(i,:)), log10(L2_inexact(i,:)),1);
        slope_inex(i) = line(1);
    else
        line_N16_in = polyfit(log10(ng_points(i,1:2)), log10(L2_inexact(i,1:2)),1); %for N=16
    
        
    end  
    
end
slope_inexact = [slope_inex;line_N16_in(1)]; 

% % comment these lines to see only the plot of the inexact integration
for i=1:5
    loglog(ng_points(i,:),L2_exact(i,:),color{end}, 'LineWidth',2), hold on
end

xlabel('N_p')
ylabel('L_2 Norm')
legend(Label, location='SW')
title('CG Inexact and Exact Integration')


figure(2) % Plotting L_2 norm of exact integration
slope_ex = zeros(size(ng_points,1)-1,1); % slope of the convergence

for i=1:5
    loglog(ng_points(i,:),L2_exact(i,:),'*-', 'LineWidth',2), hold on 
    if i<5
        line = polyfit(log10(ng_points(i,:)), log10(L2_exact(i,:)),1);
        slope_ex(i) = line(1);
    else
        line_N16_ex = polyfit(log10(ng_points(i,1:2)), log10(L2_exact(i,1:2)),1); %for N=16
    end
end

slope_exact = [slope_ex;line_N16_ex(1)]; 
xlabel('N_p')
ylabel('L_2 Norm')
legend(Label, location='NE')
title('CG Exact Integration')




