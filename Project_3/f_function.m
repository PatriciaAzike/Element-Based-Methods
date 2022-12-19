% This code creates the right hand function f(x) in the Helmoltz equation

function f = f_function(coord,npoin)

%Initialize
f=zeros(npoin,1);

%Generate Grid Points
for i=1:npoin
  x=coord(i);  
  
  f(i)=(1-pi^2)*sin(pi*x); 
  
end %i      
