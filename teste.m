clear;

dt=0.01; 
t=0:dt:100; 

p1 = 0;
p2 = 2.5;

p1s = zeros(1, length(t));
p2s = zeros(1, length(t));

i1 = i2 = 0;
alpha1 = 5; alpha2 = 2.5;
beta = 2;
gama = 4;

for i = 1:length(t)
  
  dpdt = (alpha1 / (1 + (p2 / (1 + i2) ) ^ beta )) - p1;
  p1 = dpdt * dt + p1;
  p1s(i) = p1;
  
  dpdt = (alpha2 / (1 + (p1 / (1 + i1) ) ^ gama)) - p2;
  p2 = dpdt * dt + p2;
  p2s(i) = p2;  
  
  if (i==idivide(length(t), 9))
    i2 = 5;
    disp("oi");
   endif
  
  if (i==2*idivide(length(t), 9))
    i2 = 0;
    disp("oi");
   endif
   
   if (i==3*idivide(length(t), 9))
    i1 = 5;
    disp("oi");
   endif
   
   if (i==4*idivide(length(t), 9))
    i1 = 0;
    disp("oi");
   endif
   
   if (i==5*idivide(length(t), 9))
    i2 = 5;
    disp("oi");
   endif
  
  if (i==6*idivide(length(t), 9))
    i2 = 0;
    disp("oi");
   endif
   
   if (i==7*idivide(length(t), 9))
    i1 = 5;
    disp("oi");
   endif
   
   if (i==8*idivide(length(t), 9))
    i1 = 0;
    disp("oi");
   endif

   
endfor

figure;
hold on;
grid on;
plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;');
xlabel('t');
ylabel('Concentration');