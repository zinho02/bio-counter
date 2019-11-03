clear;

dt=0.01; 
t=0:dt:200; 

%i = 1 ate 5000 

p1 = 0;
p2 = 2.5;
ii1 = 0;
ii2 = 0;

% maximal expression rates
km1 = 1;
km2 = 3;

% inibicao entre os iis
kq1 = 1;
kq2 = 1;

%dissociation constant
kp1 = 1;
kp2 = 2;

% teste
beta2 =3;
gama2 = 4;


p1s = zeros(1, length(t));
p2s = zeros(1, length(t));

% pulso representa uma explosao na concentracao de mRNAs que possuem dois operators
% o primeiro precisa de p1 AND (NOT p2) e o segundo precisa de p2 AND (NOT p1)
pulso = zeros(1, length(t));

i1s = zeros(1, length(t));
i2s = zeros(1, length(t));

alpha1 = 10; alpha2 = 5;
beta = 3;
gama = 4;

p3 = 0;
p3s = zeros(1, length(t));

for i = 2:length(t)
   
   if (mod(i, 2000) >800 && mod(i,2000)<1400)
     pulso(i) = pulso(i-1)+0.05;
   endif
   
   if (mod(i, 2000) >=1400 && mod(i,2000)<1700)
     pulso(i) = pulso(i-1)-0.1;
   endif
  
  div = 1 + p1/kp1 + p2/kp2 + (p1*p2)/(kp1*kp2);
  
  
%  verde, depende do p1 (azul) 
%  dpdt = alpha0 + km1 * pulso(i) * ((1 + p1 / (kp1)) / div) - kd1*ii1;
  dpdt = (    pulso(i) * ((1 + p1 / (kp1)) / div)     /   (1 + (ii2 / (1 + p1s(i) ) ) ^ beta2)   ) - ii1;
  ii1 = dpdt * dt + ii1;
  i1s(i) = ii1; 
  
  
%  amarelo deve ser o primeiro
% usei collins pois precisava que o crescimento de amarelo reprimisse o crescimento do verde.
%  dpdt = alpha0 + km2 * pulso(i) * ((1 + (p2) / (kp2 )) / div) - kd2*ii2;

% (pulso controla a execucao da  AND entre p1 e not p2) -> representa o maximal expression rate de collins           collins  para verde reprimir amarelo
  dpdt = (   pulso(i) * ((1 + p1 / (kp1)) / div)       /              (1 + (ii1 / (1 + p2s(i) ) ) ^ gama2)     ) - ii2;
  ii2 = dpdt * dt + ii2;
  i2s(i) = ii2; 
 
  dpdt = (alpha1 / (1 + (p2 / (1 + i2s(i) ) ) ^ beta )) - p1;
  p1 = dpdt * dt + p1;
  p1s(i) = p1;
  
  dpdt = (alpha2 / (1 + (p1 / (1 + i1s(i)) ) ^ gama)) - p2;
  p2 = dpdt * dt + p2;
  p2s(i) = p2;

endfor

figure;
hold on;
grid on;
%plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , pulso , 'k;pulso;');
plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , i2s, 'y;i2;', t ,i1s,'g;i1;', t , pulso , 'k;pulso;');
xlabel('t');
ylabel('Concentration');