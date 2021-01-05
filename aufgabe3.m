close all;

%Variablen
L       = 1;
T       = 1;   %Betrachtungszeit
n       = 100;
K       = 1;
a       = 2;
tMax    = 100;
xMax    = 100;
%%%%%%%%%%%%%
v       = zeros(xMax, tMax);
t       = linspace(0, T, tMax);
x       = linspace(0, L, xMax);
deltaX  = 1;
deltaT  = T/tMax;
d       = K*deltaT/deltaX^2;
An      = zeros(n, 1);
u       = zeros(xMax, tMax);
funAE    = @(t) a*exp(-t);
funE    = @(n, t) exp(-(n^2)*(pi^2)*t);
funSin  = @(n, x) sin(n*pi*x);

%Randwerte
for i = 1:tMax
    v(1, i) = 1;
    v(xMax, i) = funAE(t(i));
end

%Numerische Lösung (explizit)
for k = 2:tMax-1
    for j = 2:xMax-1
        v(j,k)     = d*(v(j+1, k-1) - 2*v(j, k-1) + v(j-1, k-1)) + v(j, k-1) - funAE(t(k))*deltaT;
    end
end

%Bausteine berechnen
for i = 1:n
    An(i, 1)      = (-2*(a*(-1)^(i+1)+1))/(i*pi);
end

%Exakte Lösung
for k = 2:tMax
    for j = 1:xMax
        temp  = 0;
        for i = 1:n
            temp   = temp + An(i, 1)*funE(n, t(k))*funSin(n, x(j));
        end
        u(j,k)     = (funAE(t(k))-1)*x(j) + 1 + temp;
    end
end

%Plot i)
figure('Name', 'exakt und explizite numerische Lösung (i)','NumberTitle','off')
plot(x, u(1:xMax,1),'b')

hold on 
plot(x, v(1:xMax,1),'b+')

plot(x, u(1:xMax,  ceil(tMax/6)),'r')
plot(x, v(1:xMax,  ceil(tMax/6)),'ro')

plot(x, u(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, v(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, u(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, v(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, u(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, v(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, u(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, v(1:xMax, 5* ceil(tMax/6)),'ys')
hold off
