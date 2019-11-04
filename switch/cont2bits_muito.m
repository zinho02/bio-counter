clear;

dt=0.01; 
t=0:dt:600; 

%i = 1 ate 60000

p1 = 0;
p2 = 2.5;
ii1 = 0;
ii2 = 0;

%
alpha0 = 0.1;

% maximal expression rates
km1 = 2;
km2 = 2;

%dissociation constant
kp1 = 1;
kp2 = 2;

% parametros do collins para os iis
beta2 =3;
gama2 = 4;

%mass action constant
k1 = 2;

%degradation rate
kd1 = 0.8;

p1s = zeros(1, length(t));
p2s = zeros(1, length(t));

% pulso representa uma explosao na concentracao de mRNAs que possuem dois operators
% o primeiro precisa de p1 AND (NOT p2) e o segundo precisa de p2 AND (NOT p1)
pulso = zeros(1, length(t));

i1s = zeros(1, length(t));
i2s = zeros(1, length(t));


%posso alterar de diversas estranhas maneiras os alphas
alpha1 = 10; alpha2 = 5;
beta = 3;
gama = 4;

%-------------2bits---------------
% pulso para o segundo switch
pulso2 = zeros(1, length(t));
pulso2aux = 0;

sep1 = 0;
sep2 = 2.5;
seii1 = 0;
seii2 = 0;

% maximal expression rates
sekm1 = 1;
sekm2 = 1;

%dissociation constant
sekp1 = 1;
sekp2 = 2;

% parametros do collins para os iis
sebeta2 =3;
segama2 = 4;

sep1s = zeros(1, length(t));
sep2s = zeros(1, length(t));

sei1s = zeros(1, length(t));
sei2s = zeros(1, length(t));

sealpha1 = 10; sealpha2 = 5;
sebeta = 3;
segama = 4;

for i = 2:length(t)
   
   if (mod(i, 4000) >800 && mod(i,4000)<1400)
     pulso(i) = pulso(i-1)+0.05;
   endif
   
   if (mod(i, 4000) >=1400 && mod(i,4000)<1700)
     pulso(i) = pulso(i-1)-0.1;
   endif
  
  
div = 1 + p1/kp1 + p2/kp2 + (p1*p2)/(kp1*kp2);
  
  
%  verde, depende do p1 (azul) 
%  dpdt = alpha0 + km1 * pulso(i) * ((1 + p1 / (kp1)) / div) - kd1*ii1;
  dpdt = (  km1  *pulso(i) * ((1 + p1 / (kp1)) / div)     /   (1 + (ii2 / (1 + p1s(i) ) ) ^ beta2)   ) - ii1;
  ii1 = dpdt * dt + ii1;
  i1s(i) = ii1; 
  
  
%  amarelo deve ser o primeiro
% usei collins pois precisava que o crescimento de amarelo reprimisse o crescimento do verde.
%  dpdt = alpha0 + km2 * pulso(i) * ((1 + (p2) / (kp2 )) / div) - kd2*ii2;

% (pulso controla a execucao da  AND entre p1 e not p2) -> representa o maximal expression rate de collins           collins  para verde reprimir amarelo
  dpdt = (  km2  *pulso(i) * ((1 + p2 / (kp2)) / div)       /              (1 + (ii1 / (1 + p2s(i) ) ) ^ gama2)     ) - ii2;
  ii2 = dpdt * dt + ii2;
  i2s(i) = ii2; 
 
  dpdt = (alpha1 / (1 + (p2 / (1 + i2s(i) ) ) ^ beta )) - p1;
  p1 = dpdt * dt + p1;
  p1s(i) = p1;
  
  dpdt = (alpha2 / (1 + (p1 / (1 + i1s(i)) ) ^ gama)) - p2;
  p2 = dpdt * dt + p2;
  p2s(i) = p2;
  
  
%%  ----------2bits-----------
%%  mRNA que vai funcionar como pulso2
  dpdt =  k1*i1s(i) -kd1*pulso2aux;
  pulso2aux = dpdt * dt + pulso2aux;
  pulso2(i) = pulso2aux;
  
  sediv = 1 + sep1/sekp1 + sep2/sekp2 + (sep1*sep2)/(sekp1*sekp2);
  
  dpdt = (  sekm1  *pulso2(i) * ((1 + sep1 / (sekp1)) / sediv)     /   (1 + (seii2 / (1 + sep1s(i) ) ) ^ sebeta2)   ) - seii1;
  seii1 = dpdt * dt + seii1;
  sei1s(i) = seii1; 

  dpdt = (  sekm2  *pulso2(i) * ((1 + sep2 / (sekp2)) / sediv)       /              (1 + (seii1 / (1 + sep2s(i) ) ) ^ segama2)     ) - seii2;
  seii2 = dpdt * dt + seii2;
  sei2s(i) = seii2; 

%  dpdt = ( 1.01* pulso2(i) * ((1 + sep1 / (sekp1)) / sediv)     /   (1 + (seii2 / (1 + sep1s(i) ) ) ^ sebeta2)   ) - seii1;
%  seii1 = dpdt * dt + seii1;
%  sei1s(i) = seii1; 
%    
%  dpdt = ( pulso2(i) * ((1 + sep1 / (sekp1)) / sediv)       /              (1 + (seii1 / (1 + sep2s(i) ) ) ^ segama2)     ) - seii2;
%  seii2 = dpdt * dt + seii2;
%  sei2s(i) = seii2; 
 
  dpdt = (sealpha1 / (1 + (sep2 / (1 + sei2s(i) ) ) ^ sebeta )) - sep1;
  sep1 = dpdt * dt + sep1;
  sep1s(i) = sep1;
  
  dpdt = (sealpha2 / (1 + (sep1 / (1 + sei1s(i)) ) ^ segama)) - sep2;
  sep2 = dpdt * dt + sep2;
  sep2s(i) = sep2;


endfor

figure;
hold on;
grid on;
%plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , pulso , 'k;pulso;');
plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , i2s, 'y;i2;', t ,i1s,'g;i1;', t , pulso , 'k;pulso;');
%
figure;
hold on;
grid on;
%plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , pulso , 'k;pulso;');
plot(t, sep1s, 'b;seP1;', t ,sep2s, 'r;seP2;', t , sei2s, 'y;sei2;', t ,sei1s,'g;sei1;', t , pulso , 'k;pulso;', t, pulso2, 'm;pulso2;');
xlabel('t');
ylabel('Concentration');

figure;
hold on;
grid on;
%plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , pulso , 'k;pulso;');
plot(t, p1s, 'b;P1;', t, sep1s, 'r;seP1;', t , pulso , 'k;pulso;');
xlabel('t');
ylabel('Concentration');