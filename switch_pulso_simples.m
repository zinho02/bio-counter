clear;

dt=0.01; 
t=0:dt:80; 

%i = 1 ate 5000 

p1 = 0;
p2 = 2.5;
ii1 = 0;
ii2 = 0;

%dissociation constant
kp1 = 1;
kp2 = 2;


%leak from repression
alpha0 = 0.01;

%degradation rate
kd1 = 0.3;
kd2 = 0.3;

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

for i = 2:length(t)
   
   if (mod(i, 1000) >800 || mod(i,1000)<200  && i>700)
     pulso(i) = pulso(i-1)+0.05;
   endif
   
   if ((mod(i, 1000) >=200 && mod(i,1000)<400) && i>1000 )
     pulso(i) = pulso(i-1)-0.1;
   endif
  
  div = 1 + p1/kp1 + p2/kp2 + (p1*p2)/(kp1*kp2);
  
  dpdt = alpha0 + pulso(i) * ((1 + p1 / kp1) / div) - kd1*ii1;
  ii1 = dpdt * dt + ii1;
  i1s(i) = ii1; 
  
  dpdt = alpha0 + pulso(i) * ((1 + p2 / kp2) / div) - kd2*ii2;
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