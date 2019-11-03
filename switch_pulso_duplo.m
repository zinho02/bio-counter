clear;

dt=0.01; 
t=0:dt:50; 

%i = 1 ate 5000 

p1 = 0;
p2 = 2.5;

p1s = zeros(1, length(t));
p2s = zeros(1, length(t));

i1s = zeros(1, length(t));
i2s = zeros(1, length(t));

alpha1 = 10; alpha2 = 5;
beta = 3;
gama = 4;

p3 = 0;
p3s = zeros(1, length(t));

for i = 1:length(t)
     
   if (i> 800 && i <1200)
    i2s(i) = i2s(i-1)+0.05;
   endif
   
   if (i>= 1200 && i <1400)
    i2s(i) = i2s(i-1)-0.1;
   endif
   
   if (i>1800 && i <2200)
    i1s(i) = i1s(i-1)+0.05;
   endif
   
   if (i>= 2200 && i <2400)
    i1s(i) = i1s(i-1)-0.1;
   endif
   
   if (i>2800 && i <3200)
    i2s(i) = i2s(i-1)+0.05;
   endif
   
   if (i>= 3200 && i <3400)
    i2s(i) = i2s(i-1)-0.1;
   endif
   
   if (i>3800 && i <4200)
    i1s(i) = i1s(i-1)+0.05;
   endif
   
   if (i>= 4200 && i <4400)
    i1s(i) = i1s(i-1)-0.1;
   endif
   
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
plot(t, p1s, 'b;P1;', t ,p2s, 'r;P2;', t , i2s, 'y;i2;', t ,i1s,'g;i1;');
xlabel('t');
ylabel('Concentration');