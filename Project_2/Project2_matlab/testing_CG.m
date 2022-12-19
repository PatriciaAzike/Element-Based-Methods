clear all; 
close all;
format long
clc

tic


for i=1:2 %integration_type: 1- inexact, 2- exact
    if i==1
        N=[1 4 8 16];
        NeM=[16 32 64; 4 8 16; 2 4 8;1 2 4]; %Matrix containing Ne
        ng_points=zeros(length(N),size(NeM,2));
        L2_inexact=zeros(length(N),size(NeM,2));
        for j =1:length(N)
            for k=1:size(NeM,2)
                ng_points(j,k)=N(j)*NeM(j,k) + 1;
                [L2_inexact(j,k)]=driver_CG(i,N(j),  NeM(j,k));
            end     
        end
         
    else
        N=[1 4 8 16];
        NeM=[16 32 64; 4 8 16; 2 4 8;1 2 4];
        L2_exact=zeros(length(N),size(NeM,2));
        for j =1:length(N)
            for k=1:size(NeM,2)
                [L2_exact(j,k)]=driver_CG(i,N(j),  NeM(j,k));
            end   
        end
    end
end



figure(1)% Plotting L_2 norm of inexact integration
Label={'N=1','N=4','N=8','N=16'};
for i=1:4
    semilogy(ng_points(1,:),L2_inexact(i,:),'*-', 'LineWidth',2), hold on
    
end

xlabel('N_p')
ylabel('L_2 Norm')
legend(Label, location='SW')
title('CG Inexact Integration')

figure(2) % Plotting L_2 norm of exact integration
for i=1:4
    semilogy(ng_points(1,:),L2_exact(i,:),'*-', 'LineWidth',2), hold on  
end
xlabel('N_p')
ylabel('L_2 Norm')
legend(Label, location='SW')
title('CG Exact Integration')




