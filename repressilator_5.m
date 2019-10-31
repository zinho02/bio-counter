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
m4 = 0;
m5 = 0;
p1 = 0;
p2 = 5;
p3 = 15;
p4 = 0;
p5 = 0;

m1s = zeros(1, length(t));
m2s = zeros(1, length(t));
m3s = zeros(1, length(t));
m4s = zeros(1, length(t));
m5s = zeros(1, length(t));
p1s = zeros(1, length(t));
p2s = zeros(1, length(t));
p3s = zeros(1, length(t));
p4s = zeros(1, length(t));
p5s = zeros(1, length(t));

s = zeros(1, length(t));

for i = 1:length(t)
  dm1dt = alfa0 .+ (alfa ./ (1 .+ p5.^n)) .- m1;
  m1 = dm1dt * dt + m1;
  m1s(i) = m1;
  
  dm2dt = alfa0 .+ (alfa ./ (1 .+ p1.^n)) .- m2;
  m2 = dm2dt * dt + m2;
  m2s(i) = m2;
  
  dm3dt = alfa0 .+ (alfa ./ (1 .+ p2.^n)) .- m3;
  m3 = dm3dt * dt + m3;
  m3s(i) = m3;
  
  dm4dt = alfa0 .+ (alfa ./ (1 .+ p3.^n)) .- m4;
  m4 = dm4dt * dt + m4;
  m4s(i) = m4;
  
  dm5dt = alfa0 .+ (alfa ./ (1 .+ p4.^n)) .- m5;
  m5 = dm5dt * dt + m5;
  m5s(i) = m5;
  
  dp1dt = (beta .* m1 .- beta .* p1);
  p1 = dp1dt * dt + p1;
  p1s(i) = p1;
  
  dp2dt = beta .* m2 .- beta .* p2;
  p2 = dp2dt * dt + p2;
  p2s(i) = p2;
  
  dp3dt = (beta .* m3 .- beta .* p3);
  p3 = dp3dt * dt + p3;
  p3s(i) = p3;
  
  dp4dt = (beta .* m4 .- beta .* p4);
  p4 = dp4dt * dt + p4;
  p4s(i) = p4;
  
  dp5dt = (beta .* m5 .- beta .* p5);
  p5 = dp5dt * dt + p5;
  p5s(i) = p5;
  
endfor

figure;
hold on;
grid on;

plot(t, p1s, 'm;p1;', t, p2s, 'g;p2;', t, p3s, 'b;p3;', t, p4s, 'k;p4;', t, p5s, 'r;p5;');
xlabel('t');
ylabel('Concentration');