close all;

%Variablen
L       = 1;
T       = 0.005;   %Betrachtungszeit (<= 0.005)
n       = 100;
K       = 1;
tMax    = 100;
xMax    = 100;
%%%%%%%%%
t       = linspace(0, T, tMax);
x       = linspace(0, L, xMax);
deltaX  = L/xMax;
deltaT  = T/tMax;
% deltaT  = (1/(2*K))*deltaX^2;
d       = K*deltaT/deltaX^2;
%Bausteine
B_nI    = zeros(n, 1);
B_nII   = zeros(n, 1);
B_nIII  = zeros(n, 1);
%r-Vektor
rI      = zeros(xMax, 1);
%numerisch explizit
vI      = zeros(xMax, tMax);
vII     = zeros(xMax, tMax);
vIII    = zeros(xMax, tMax);
%exakt
uI      = zeros(xMax, tMax);
uII     = zeros(xMax, tMax);
uIII    = zeros(xMax, tMax);
%implizit
wI      = zeros(xMax, tMax);
wII     = zeros(xMax, tMax);
wIII    = zeros(xMax, tMax);
%LGS-Matrix (tridiagonal)
A       = zeros(xMax, tMax);
%Funktionen für exakte Lösung
funE    = @(n, t) exp(-((n*pi/L)^2)*K*t);
funSin  = @(n, x) sin((n*pi*x)/L); 
%Anfangsbedingungen
funI    = @(x) 3*sin(pi*x) + 5*sin(4*pi*x);
funII   = @(x) 1;
funIII  = @(x) x;

%Anfangswerte
for j = 1:xMax
    uI(j, 1)    = funI(x(j));
    uII(j, 1)   = funII(x(j));
    uIII(j, 1)  = funIII(x(j));
    vI(j, 1)    = funI(x(j));
    vII(j, 1)   = funII(x(j));
    vIII(j, 1)  = funIII(x(j));
end

%Bausteine berechnen
for i = 1:n
    B_nI(i, 1)      = (2/L) * integral(@(x) funI(x).*funSin(i, x),0,L);
    B_nII(i, 1)     = (2/L) * integral(@(x) funII(x).*funSin(i, x),0,L);
    B_nIII(i, 1)    = (2/L) * integral(@(x) funIII(x).*funSin(i, x),0,L);
end

%Exakte Lösung
for j = 2:xMax-1
    for k = 2:tMax-1
        tempI   = 0;
        tempII  = 0;
        tempIII = 0;
        for i = 1:n
            tempI       = tempI + (B_nI(i, 1)*funE(i, t(k))*funSin(i, x(j)));
            tempII      = tempII + (B_nII(i, 1)*funE(i, t(k))*funSin(i, x(j)));
            tempIII     = tempIII + (B_nIII(i, 1)*funE(i, t(k))*funSin(i, x(j)));
        end
        uI(j,k)     = tempI;
        uII(j,k)    = tempII;
        uIII(j,k)   = tempIII;
    end
end

%Numerische Lösung (explizit)
for k = 2:tMax-1
    for j = 2:xMax-1
        vI(j,k)     = d*(vI(j+1, k-1) - 2*vI(j, k-1) + vI(j-1, k-1)) + vI(j, k-1);
        vII(j,k)    = d*(vII(j+1, k-1) - 2*vII(j, k-1) + vII(j-1, k-1)) + vII(j, k-1);
        vIII(j,k)   = d*(vIII(j+1, k-1) - 2*vIII(j, k-1) + vIII(j-1, k-1)) + vIII(j, k-1);
    end
end

%Tridiagonalmatrix A
vec1    = zeros(xMax-2, 1);
vec2    = zeros(xMax-3, 1);

for j = 1:xMax
    vec1(j, 1) = 1+2*d;
end
for j = 1:xMax-1
    vec2(j, 1) = -d;
end
A1 = diag(vec1);
A2 = diag(vec2, 1);
A3 = diag(vec2, -1);
A = A1 + A2 + A3;

RI = zeros(xMax, 1);
RII = zeros(xMax, 1);
RIII = zeros(xMax, 1);
for i = 1:xMax-1
  RI(i, 1) = funI(x(i+1));
  RII(i, 1) = funII(x(i+1));
  RIII(i, 1) = funIII(x(i+1));
end

for i = 1:xMax-1
    xVecI = A\RI(:,i); % erster x-Wert, wird neuer r-Vektor
    RI = [RI, xVecI];
    xVecII = A\RII(:,i); % erster x-Wert, wird neuer r-Vektor
    RII = [RII, xVecII];
    xVecIII = A\RIII(:,i); % erster x-Wert, wird neuer r-Vektor
    RIII = [RIII, xVecIII];
end

for i = 2:xMax-1
    wI(i,:) = RI(i-1,:);
    wII(i,:) = RII(i-1,:);
    wIII(i,:) = RIII(i-1,:);
end

%%%Plots
%Plot i)
figure('Name', 'exakt und explizite numerische Lösung (i)','NumberTitle','off')
plot(x, uI(1:xMax,1),'b')

hold on 
plot(x, vI(1:xMax,1),'b+')

plot(x, uI(1:xMax,  ceil(tMax/6)),'r')
plot(x, vI(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uI(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, vI(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uI(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, vI(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uI(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, vI(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uI(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, vI(1:xMax, 5* ceil(tMax/6)),'ys')
hold off

%Plot ii)
figure('Name', 'exakt und explizite numerische Lösung (ii)','NumberTitle','off')
plot(x, uII(1:xMax,1),'b')

hold on 
plot(x, vII(1:xMax,1),'b+')

plot(x, uII(1:xMax,  ceil(tMax/6)),'r')
plot(x, vII(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uII(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, vII(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uII(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, vII(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uII(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, vII(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uII(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, vII(1:xMax, 5* ceil(tMax/6)),'ys')
hold off

%Plot ii)
figure('Name', 'exakt und explizite numerische Lösung (iii)','NumberTitle','off')
plot(x, uIII(1:xMax,1),'b')

hold on 
plot(x, vIII(1:xMax,1),'b+')

plot(x, uIII(1:xMax,  ceil(tMax/6)),'r')
plot(x, vIII(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uIII(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, vIII(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uIII(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, vIII(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uIII(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, vIII(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uIII(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, vIII(1:xMax, 5* ceil(tMax/6)),'ys')
hold off

%Plot i)
figure('Name', 'exakt und implizite numerische Lösung (i)','NumberTitle','off')
plot(x, uI(1:xMax,1),'b')

hold on 
plot(x, vI(1:xMax,1),'b+')

plot(x, uI(1:xMax,  ceil(tMax/6)),'r')
plot(x, wI(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uI(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, wI(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uI(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, wI(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uI(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, wI(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uI(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, wI(1:xMax, 5* ceil(tMax/6)),'ys')
hold off

%Plot ii)
figure('Name', 'exakt und implizite numerische Lösung (ii)','NumberTitle','off')
plot(x, uII(1:xMax,1),'b')

hold on 
plot(x, vII(1:xMax,1),'b+')

plot(x, uII(1:xMax,  ceil(tMax/6)),'r')
plot(x, wII(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uII(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, wII(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uII(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, wII(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uII(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, wII(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uII(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, wII(1:xMax, 5* ceil(tMax/6)),'ys')
hold off

%Plot iii)
figure('Name', 'exakt und implizite numerische Lösung (iii)','NumberTitle','off')
plot(x, uIII(1:xMax,1),'b')

hold on 
plot(x, vIII(1:xMax,1),'b+')

plot(x, uIII(1:xMax,  ceil(tMax/6)),'r')
plot(x, wIII(1:xMax,  ceil(tMax/6)),'ro')

plot(x, uIII(1:xMax, 2* ceil(tMax/6)),'m')
plot(x, wIII(1:xMax, 2* ceil(tMax/6)),'m*')

plot(x, uIII(1:xMax, 3* ceil(tMax/6)),'k')
plot(x, wIII(1:xMax, 3* ceil(tMax/6)),'k.')

plot(x, uIII(1:xMax, 4* ceil(tMax/6)),'g')
plot(x, wIII(1:xMax, 4* ceil(tMax/6)),'gd')

plot(x, uIII(1:xMax, 5* ceil(tMax/6)),'y')
plot(x, wIII(1:xMax, 5* ceil(tMax/6)),'ys')
hold off
