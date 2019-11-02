clear;
%
%dt=0.01; 
%t=0:dt:300; 
%
%p1 = 0;
%p2 = 2.5;
%
%p1s = zeros(1, length(t));
%p2s = zeros(1, length(t));
%
%i1 = i2 = 0;
%alpha1 = 3; alpha2 = 2.5;
%beta = gama = 4;

%for i = 1:length(t)
%  
%  dpdt = (alpha1 / (1 + (p2 / (1 + i2) ) ^ beta )) - p1;
%  p1 = dpdt * dt + p1;
%  p1s(i) = p1;
%  
%  dpdt = (alpha2 / (1 + (p1 / (1 + i1) ) ^ gama)) - p2;
%  p2 = dpdt * dt + p2;
%  p2s(i) = p2;  
%
%endfor

%figure;
%hold on;
%grid on;
%
%%plot(t p1s, 'b,P1', t p2s, 'r,P2');
%xlabel('t');
%ylabel('Concentration');