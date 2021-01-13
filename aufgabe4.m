close all;
clear;

%veränderbare Variablen
T       = 1;        %obere Grenze für t (Betrachtungszeit)
n       = 100;      %Anzahl der Summanden pro Zeitschritt
eps     = 1;        %Epsilon (aus Aufgabenstellung)
a       = 0;        %variabel (aus Aufgabenstellung)
alpha   = 1;        %Konstante aus Separationsansatz (X * T = -alpha^2)
xMax    = 100;      %Anzahl Wegschritte
%abhängige Variablen (bzw. vorgegeben aus Aufgabenstellung)
L       = 1;                                %obere Intervallgrenze für x
x       = linspace(0, L, xMax);             %Unterteilung Weg
deltaX  = L/xMax;                           %Abstand zwischen 2 Weg-Schritten
deltaT  = (0.5*deltaX^2)/eps;               %Abstand zwischen 2 Zeit-Schritten
tMax    = T/deltaT;                         %Anzahl Zeitschritte
t       = linspace(0, T, tMax);             %Unterteilung Zeit
u       = zeros(xMax, tMax);                %exakte Lösung
v       = zeros(xMax, tMax);                %numerische Lösung
%%%%%%%%%%%%%%%%%
%Funktionen für exakte Lösung
funSin  = @(x) sin((3*pi*x)/2); 
funE    = @(x) exp(a*x);
funEx   = @(x, t) exp(a*x - alpha^2*t);

%Randwerte
for i = 1:xMax
   u(i, 1) = funE(x(i))*funSin(x(i));
end
for i = 1:tMax
   u(1, i) = 0;
end

%Numerische Löung (explizit)
for k = 2:tMax-1
%     if (k > 2)  
%         u(xMax, k-1) = 2*a*deltaX*u(xMax-1,k-2)+u(xMax-2,k-2);
%     end
    for j = 2:xMax-1
        u(j, k) = (eps*(u(j+1, k-1) - 2*u(j, k-1) - u(j-1, k-1))/deltaX - 2*a*(u(j, k-1) - u(j-1, k-1)))*(deltaT/deltaX) + u(j, k-1);     
    end
    u(xMax, k) = 2*a*deltaX*u(xMax-1,k-1)+u(xMax-2,k-1);
end

%Exakte Lösung
for k = 1:tMax
    for j = 1:xMax
        v(j, k) = funEx(x(j), t(k))*funSin(x(j));
    end
end

%Plot
figure('Name', 'exakte Lösung und explizite numerische Lösung','NumberTitle','off')
plot(x, v(1:xMax,1),'b')
xlabel('Betrachtetes Objekt','FontAngle','italic');
ylabel('Temperatur','FontAngle','italic');
hold on 
plot(x, u(1:xMax,1),'b+')

plot(x, v(1:xMax,  ceil(tMax/6)),'r')
plot(x, u(1:xMax,  ceil(tMax/6)),'ro')

plot(x, v(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, u(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, v(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, u(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, v(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, u(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, v(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, u(1:xMax, 5* ceil(tMax/6)),'ys')
hold off