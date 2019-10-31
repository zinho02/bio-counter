clear;

dt=0.01; 
t=0:dt:800; 

alfa0 = 0.03;
alfa = 298.2;
beta = 0.2;
n = 2;

m1 = 0;
m2 = 0;
m3 = 0;
p1 = 0;
p2 = 5;
p3 = 15;

m1s = zeros(1, length(t));
m2s = zeros(1, length(t));
m3s = zeros(1, length(t));
p1s = zeros(1, length(t));
p2s = zeros(1, length(t));
p3s = zeros(1, length(t));

s = zeros(1, length(t));

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
    s(i) = 800;
  endif
  
  if p1 > 25 && p2 > 25 && p3 < 25
    s(i) = 250;
  endif
  
  if p1 > 25 && p2 < 25 && p3 > 25
    s(i) = 350;
  endif
  
  if p1 > 25 && p2 < 25 && p3 < 25
    s(i) = 300;
  endif
  
  if p1 < 25 && p2 > 25 && p3 > 25
    s(i) = 450;
  endif
  
  if p1 < 25 && p2 > 25 && p3 < 25
    s(i) = 200;
  endif
  
  if p1 < 25 && p2 < 25 && p3 > 25
    s(i) = 400;
  endif
  
  if p1 < 25 && p2 < 25 && p3 < 25
    s(i) = 700;
  endif
endfor

figure;
hold on;
grid on;

plot(t, s, 'm;S;');
xlabel('t');
ylabel('Concentration');