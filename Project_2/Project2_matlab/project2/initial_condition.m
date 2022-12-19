%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe] = initial_condition(coord,npoin)

%Initialize
qe=zeros(npoin,1);

%Generate Grid Points
for i=1:npoin
  x=coord(i);  
  
  qe(i)=exp( -64.0*x^2 ); %IC used in Fig 5.9
  
end %i      


      
