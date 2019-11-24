clear;

dt=0.05; 
t=0:dt:400; 

alfa0 = 0.03;
alfa = 298.2;
beta = 0.2;
n = 2;

ms = zeros(1,5);
ps = zeros(1,5); ps(1) = 5; ps(3) = 15; ps(4)  = 30;

mcons = zeros(5, length(t));
pcons = zeros(5, length(t));

for i = 1:length(t)
  for j = 1:5
%    mRNA concentration
    dmdt = alfa0 + (alfa ./ (1 .+ ps(1 + mod((j+3), 5)) .^n)) - ms(j);
    ms(j) = dmdt * dt + ms(j);
    mcons(j,i) = ms(j);
  endfor
  
  for j = 1:5
%    protein concentration
    dpdt = (beta .* ms(j) .- beta .* ps(j));
    ps(j) = dpdt * dt + ps(j);
    pcons(j,i) = ps(j);
  endfor
  
endfor

figure;
hold on;
grid on;
plot(t, pcons(1,:), 'g;P1;',t, pcons(2,:), 'b;P2;',t, pcons(3,:), 'm;P3;',t, pcons(4,:), 'k;P4;',t, pcons(5,:), 'c;P5;');
plot(t,pcons(1,:));
xlabel('t');
ylabel('Concentration');