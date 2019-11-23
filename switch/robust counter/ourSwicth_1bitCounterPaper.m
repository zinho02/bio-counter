clear;

dt=0.01; 
t=0:dt:80; 

%i = 1 ate 20000 

p1 = 0;
p2 = 2.5;
ii1 = 0;
ii2 = 0;

% maximal expression rates
km1 = 4;
km2 = 4;

%dissociation constant
kp1 = 1;
kp2 = 2;

% parametros do collins para os iis
beta2 =4;
gama2 = 4;


p1s = zeros(1, length(t));
p2s = zeros(1, length(t));

pulso = zeros(1, length(t));

i1s = zeros(1, length(t));
i2s = zeros(1, length(t));

alpha1 = 10; alpha2 = 5;
beta = 3;
gama = 4;

p3 = 0;
p3s = zeros(1, length(t));

soma = 0.1
for i = 2:length(t)
   
   if (mod(i, 2000) >800 && mod(i,2000)<1400)
     pulso(i) = pulso(i-1)+soma;
   endif
   
   if (mod(i, 2000) >=1400 && mod(i,2000)<1700)
     pulso(i) = pulso(i-1)-soma*2;
   endif
  
  div = 1 + p1/kp1 + p2/kp2 + (p1*p2)/(kp1*kp2);
  
%  i1
  dpdt = (  km1  *pulso(i) * ((1 + p1 / (kp1)) / div)     /   (1 + (ii2 / (1 + p1s(i) ) ) ^ beta2)   ) - ii1;
  ii1 = dpdt * dt + ii1;
  i1s(i) = ii1; 
  
%  i2  
  dpdt = (  km2  *pulso(i) * ((1 + p2 / (kp2)) / div)       /              (1 + (ii1 / (1 + p2s(i) ) ) ^ gama2)     ) - ii2;
  ii2 = dpdt * dt + ii2;
  i2s(i) = ii2; 
 
% p1
  dpdt = (alpha1 / (1 + (p2 / (1 + i2s(i) ) ) ^ beta )) - p1;
  p1 = dpdt * dt + p1;
  p1s(i) = p1;

% p2  
  dpdt = (alpha2 / (1 + (p1 / (1 + i1s(i)) ) ^ gama)) - p2;
  p2 = dpdt * dt + p2;
  p2s(i) = p2;

endfor

%figure;
%hold on;
%plot(t, p1s, 'b;R1;', t ,p2s, 'r;R2;', t , i2s, 'm;i2;', t ,i1s,'g;i1;', t , pulso , 'k;Pulse;');
%xlabel('Time');
%ylabel('Concentration');

figure;
hold on;
plot(t, p1s, ':k;R1;', t , p2s, '-k;R2;', t , pulso , '--k;Pulse;');
xlabel('Time');
ylabel('Concentration');

figure;
hold on;
plot(t,i1s, 'k;i1;', t , i2s, ':k;i2;', t , pulso , '--k;Pulse;');
xlabel('Time');
ylabel('Concentration');