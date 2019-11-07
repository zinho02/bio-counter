clear;

dt=0.01; 
t=0:dt:300; 

% parametros mRNAs e proteinas
alfa0 = 0.03;
alfa = 298.2;
beta = 0.2;
n = 2;

% valores inicias mRNAs e proteinas e promoters em cada estado
m1 = 0;
m2 = 0;
m3 = 0;
p1 = 0;
p2 = 5;
p3 = 15;
con1 = 0; con2 = 0; con3 = 0; con4 = 0; con5 = 0; con6 = 0;

% concentracao de mRNAs e proteinas ao passar do tempo
m1s = zeros(1, length(t));
m2s = zeros(1, length(t));
m3s = zeros(1, length(t));
p1s = zeros(1, length(t));
p2s = zeros(1, length(t));
p3s = zeros(1, length(t));

% states
s = zeros(1, length(t));

% concentracao de promoters  ao passar do tempo
cons1 = zeros(1, length(t)); cons2 = zeros(1, length(t));cons3 = zeros(1, length(t));
cons4 = zeros(1, length(t));cons5 = zeros(1, length(t));cons6 = zeros(1, length(t)); 


%parametros da producao
k1 = 10;
k2 = 10;
k3 = 10;

%fazer um kq para cada e testar os valores
kq = 5;

%par deg
kd1 = 0.2;
kd2 = 0.2;
kd3 = 0.2;
kd4 = 0.2;
kd5 = 0.2;
kd6 = 0.2;

% maximal expressions rates
alpha1 = alpha3 = alpha5 = 1;
alpha2 = alpha4 = alpha6 = 3.77;


ajuste = 10;

for i = 1:length(t)
  dm1dt = alfa0 .+ (alfa ./ (1 .+ p3.^n)) .- m1;
  m1 = dm1dt * dt + m1;
  m1s(i) = m1;
  
  dm2dt = alfa0 .+ (alfa ./ (1 .+ p1.^n)) .- m2;
  m2 = dm2dt * dt + m2;
  m2s(i) = m2;
  
  dm3dt = alfa0 .+ (alfa ./ (1 .+ p2.^n)) .- m3;
  m3 = dm3dt * dt + m3;
  m3s(i) = m3;
  
  dp1dt = (beta .* m1 .- beta .* p1);
  p1 = dp1dt * dt + p1;
  p1s(i) = p1;
  
  dp2dt = beta .* m2 .- beta .* p2;
  p2 = dp2dt * dt + p2;
  p2s(i) = p2;
  
  dp3dt = (beta .* m3 .- beta .* p3);
  p3 = dp3dt * dt + p3;
  p3s(i) = p3;
  
  if p1 > 25 && p2 > 25 && p3 > 25
    s(i) = 8*ajuste;
  endif
  
  if p1 > 25 && p2 > 25 && p3 < 25
    s(i) = 2.5*ajuste;
  endif
  
  if p1 > 25 && p2 < 25 && p3 > 25
    s(i) = 3.5*ajuste;
  endif
  
  if p1 > 25 && p2 < 25 && p3 < 25
    s(i) = 3*ajuste;
  endif
  
  if p1 < 25 && p2 > 25 && p3 > 25
    s(i) = 4.5*ajuste;
  endif
  
  if p1 < 25 && p2 > 25 && p3 < 25
    s(i) = 2*ajuste;
  endif
  
  if p1 < 25 && p2 < 25 && p3 > 25
    s(i) = 4*ajuste;
  endif
  
  if p1 < 25 && p2 < 25 && p3 < 25
    s(i) = 7*ajuste;
  endif
  
div = (1 + (p3 / k3) + (p2 / k2) + (p1 / k1) + (p1 * p2 * p3 / (k1 * k2 * k3 * kq)) + (p1 * p2 / (k1 * k2 * kq)) + (p1 * p3 / (k1 * k3 * kq)) + (p2 * p3 / (k2 * k3 * kq))); 

%  ------2------
  dconsdt = alpha1 * (p2 / k2) / div; 
  dconsdt -= con1*kd1;
  con1 = dt * dconsdt + con1;
  cons1(i) = con1;

%  ------2.5----
  dconsdt = alpha2 * (p1 * p2 / (k1 * k2 * kq)) / div;
  dconsdt -= con2*kd2;
  con2 = dt * dconsdt + con2;
  cons2(i) = con2;
  
%  ------3----
  dconsdt = alpha3 * (p1 / k1) / div;
  dconsdt -= con3*kd3;
  con3 = dt * dconsdt + con3;
  cons3(i) = con3;

%  ------3.5----
  dconsdt = alpha4 * (p1 * p3 / (k1 * k3 * kq)) / div;
  dconsdt -= con4*kd4;
  con4 = dt * dconsdt + con4;
  cons4(i) = con4;

%  ------4----
  dconsdt = alpha5 * (p3 / k3) / div;
  dconsdt -= con5*kd5;
  con5 = dt * dconsdt + con5;
  cons5(i) = con5;

%  ------4.5----
  dconsdt = alpha6 * (p3 * p2 / (k3 * k2 * kq)) / div;
  dconsdt -= con6*kd6;
  con6 = dt * dconsdt + con6;
  cons6(i) = con6;

endfor

figure;
hold on;
grid on;
plot(t,s, 'm;States;', t, p1s, 'r;Pi;',  t, p2s, 'b;P2;',  t,  p3s,'g;P3;');
%plot(t, s, 'm;States;', t, cons1, 'r;P2;',  t, cons2, 'b;P1P2;',  t,  cons3, 'g;P1;',  t, cons4, 'y;P1P3;', t, cons5, 'c;P3;',  t, cons6, 'k;P2P3;');
xlabel('t');
ylabel('Concentration');