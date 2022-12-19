%This function computes the exact Solution.
%---------------------------------------------------------------------%
function [qe] = q_exact(coord,npoin)

%Initialize
qe=zeros(npoin,1);

%Generate Grid Points
for i=1:npoin
  x=coord(i);  
  
  qe(i)=sin(pi*x); 
  
end %i      


      
