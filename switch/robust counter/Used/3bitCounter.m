clear;

dt=0.01; 
t=0:dt:250; 

%i = 1 ate 26000

p1 = 0;
p2 = 2.5;
ii1 = 0;
ii2 = 0;

%
alpha0 = 0.1;

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
beta = 4;
gama = 4;


%2 and 3 bits
%constants
k2 = 2;
k3 = 1.5;

%degradation rate
kd2 = 0.8;
kd3 = 0.9;
%end 2 and 3 bits

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
sebeta2 =4;
segama2 = 4;

sep1s = zeros(1, length(t));
sep2s = zeros(1, length(t));

sei1s = zeros(1, length(t));
sei2s = zeros(1, length(t));

sealpha1 = 10; sealpha2 = 5;
sebeta =4;
segama = 4;

%-------------3bits---------------
% pulso para o terceiroswitch
pulso3 = zeros(1, length(t));
pulso3aux = 0;

thp1 = 0;
thp2 = 2.5;
thii1 = 0;
thii2 = 0;

% maximal expression rates
thkm1 = 1;
thkm2 = 1;

%dissociation constant
thkp1 = 1;
thkp2 = 2;

% parametros do collins para os iis
thbeta2 =4;
thgama2 = 4;

thp1s = zeros(1, length(t));
thp2s = zeros(1, length(t));

thi1s = zeros(1, length(t));
thi2s = zeros(1, length(t));

thalpha1 = 10; thalpha2 = 5;
thbeta = 4;
thgama = 4;


soma = 0.1;
for i = 2:length(t)
   
   if (mod(i, 3000) >800 && mod(i,3000)<1400)
     pulso(i) = pulso(i-1)+soma;
   endif
   
   if (mod(i, 3000) >=1400 && mod(i,3000)<1700)
     pulso(i) = pulso(i-1)-soma*2;
   endif
  
  
div = 1 + p1/kp1 + p2/kp2 + (p1*p2)/(kp1*kp2);
  
  dpdt = (  km1  *pulso(i) * ((1 + p1 / (kp1)) / div)     /   (1 + (ii2 / (1 + p1s(i) ) ) ^ beta2)   ) - ii1;
  ii1 = dpdt * dt + ii1;
  i1s(i) = ii1; 
  
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
  dpdt =  k2*i1s(i) -kd2*pulso2aux;
  pulso2aux = dpdt * dt + pulso2aux;
  pulso2(i) = pulso2aux;
  
  sediv = 1 + sep1/sekp1 + sep2/sekp2 + (sep1*sep2)/(sekp1*sekp2);
  
  dpdt = (  sekm1  *pulso2(i) * ((1 + sep1 / (sekp1)) / sediv)     /   (1 + (seii2 / (1 + sep1s(i) ) ) ^ sebeta2)   ) - seii1;
  seii1 = dpdt * dt + seii1;
  sei1s(i) = seii1; 

  dpdt = (  sekm2  *pulso2(i) * ((1 + sep2 / (sekp2)) / sediv)       /              (1 + (seii1 / (1 + sep2s(i) ) ) ^ segama2)     ) - seii2;
  seii2 = dpdt * dt + seii2;
  sei2s(i) = seii2; 
 
  dpdt = (sealpha1 / (1 + (sep2 / (1 + sei2s(i) ) ) ^ sebeta )) - sep1;
  sep1 = dpdt * dt + sep1;
  sep1s(i) = sep1;
  
  dpdt = (sealpha2 / (1 + (sep1 / (1 + sei1s(i)) ) ^ segama)) - sep2;
  sep2 = dpdt * dt + sep2;
  sep2s(i) = sep2;

  
 %%  ----------3bits-----------
  dpdt =  k3*sei1s(i) -kd3*pulso3aux;
  pulso3aux = dpdt * dt + pulso3aux;
  pulso3(i) = pulso3aux;
  
  thdiv = 1 + thp1/thkp1 + thp2/thkp2 + (thp1*thp2)/(thkp1*thkp2);

  dpdt = (  thkm1  *pulso3(i) * ((1 + thp1 / (thkp1)) / thdiv)     /   (1 + (thii2 / (1 + thp1s(i) ) ) ^ thbeta2)   ) - thii1;
  thii1 = dpdt * dt + thii1;
  thi1s(i) = thii1; 

  dpdt = (  thkm2  *pulso3(i) * ((1 + thp2 / (thkp2)) / thdiv)       /              (1 + (thii1 / (1 + thp2s(i) ) ) ^ thgama2)     ) - thii2;
  thii2 = dpdt * dt + thii2;
  thi2s(i) = thii2;  
 
  dpdt = (thalpha1 / (1 + (thp2 / (1 + thi2s(i) ) ) ^ thbeta )) - thp1;
  thp1 = dpdt * dt + thp1;
  thp1s(i) = thp1;
  
  dpdt = (thalpha2 / (1 + (thp1 / (1 + thi1s(i)) ) ^ thgama)) - thp2;
  thp2 = dpdt * dt + thp2;
  thp2s(i) = thp2;
  
endfor

figure;
hold on;
plot(t, p1s, ':k;R1;', t, sep1s, '--k;R3;', t, thp1s, 'k;R5;', t , pulso , '-.k;Pulse;');
xlabel('Time');
ylabel('Concentration');