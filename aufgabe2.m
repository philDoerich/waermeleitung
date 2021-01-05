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
d       = K*deltaT/deltaX^2;
%Bausteine
A_nI    = zeros(n, 1);
A_nII   = zeros(n, 1);
%r-Vektor
rI      = zeros(xMax, 1);
%numerisch explizit
vI      = zeros(xMax, tMax);
vII     = zeros(xMax, tMax);
%exakt
uI      = zeros(xMax, tMax);
uII     = zeros(xMax, tMax);
%implizit
wI      = zeros(xMax, tMax);
wII     = zeros(xMax, tMax);
%LGS-Matrix (tridiagonal)
A       = zeros(xMax, tMax);
%Funktionen für exakte Lösung
funE    = @(n, t) exp(-((n*pi/L)^2)*K*t);
funCos  = @(n, x) cos((n*pi*x)/L); 
%Anfangsbedingungen
funI    = @(x) 9+3*cos(pi*x) + 5*cos(4*pi*x);
funII   = @(x) x;


%Randwerte
for j = 1:xMax
    uI(j, 1)    = funI(x(j));
    uII(j, 1)   = funII(x(j));
    vI(j, 1)    = funI(x(j));
    vII(j, 1)   = funII(x(j));   
end

%Bausteine berechnen nach Neumann
for i = 1:n
    A_nI(i, 1)      = (2/L) * integral(@(x) funI(x).*funCos(i, x),0,L);
    A_nII(i, 1)     = (2/L) * integral(@(x) funII(x).*funCos(i, x),0,L);
end

%Exakte Lösung mit der Formel von Neumann
for j = 1:xMax
    for k = 2:tMax-1
        tempI   = 0;
        tempII  = 0;
        for i = 1:n
            tempI       = tempI + (A_nI(i, 1)*funE(i, t(k))*funCos(i, x(j)));
            tempII      = tempII + (A_nII(i, 1)*funE(i, t(k))*funCos(i, x(j)));
        end
        uI(j,k)     = tempI + 0.5*(2/L) * integral(@(x) funI(x),0,L);
        uII(j,k)    = tempII + 0.5*(2/L) * integral(@(x) funII(x),0,L);
    end
end


%Numerische Lösung (explizit)
for k = 2:tMax-1
    vI(1,k)     = d*(vI(2, k-1) - 2*vI(1, k-1) + vI(2, k-1)) + vI(1, k-1);
    vII(1,k)    = d*(vII(2, k-1) - 2*vII(1, k-1) + vII(2, k-1)) + vII(1, k-1);
    
    vI(xMax,k)     = d*(vI(xMax-1, k-1) - 2*vI(xMax, k-1) + vI(xMax-1, k-1)) + vI(xMax, k-1);
    vII(xMax,k)    = d*(vII(xMax-1, k-1) - 2*vII(xMax, k-1) + vII(xMax-1, k-1)) + vII(xMax, k-1);
    for j = 2:xMax-1
        vI(j,k)     = d*(vI(j+1, k-1) - 2*vI(j, k-1) + vI(j-1, k-1)) + vI(j, k-1);
        vII(j,k)    = d*(vII(j+1, k-1) - 2*vII(j, k-1) + vII(j-1, k-1)) + vII(j, k-1);
    end
end


%Numerische Lösung (implizit)
%Tridiagonalmatrix A
vec1    = zeros(xMax-2, 1);
vec2    = zeros(xMax-3, 1);

for j = 1:xMax+2
    vec1(j, 1) = 1+2*d;
end
for j = 1:xMax+1
    vec2(j, 1) = -d;
end
A1 = diag(vec1);
A2 = diag(vec2, 1);
A3 = diag(vec2, -1);
A = A1 + A2 + A3;

RI = zeros(xMax+2, 1);
RII = zeros(xMax+2, 1);

for i = 2:xMax+1
  RI(i, 1) = funI(x(i-1));
  RII(i, 1) = funII(x(i-1));
end
RI(1, 1) = funI(x(2));
RII(1, 1) = funII(x(2));
RI(xMax+2, 1) = funI(x(xMax-1));
RII(xMax+2, 1) = funII(x(xMax-1));

for i = 1:xMax+1
    xVecI = A\RI(:,i); % erster x-Wert, wird neuer r-Vektor
%     xVecI(1, 1) = xVecI(3, 1);
%     xVecI(xMax+2, 1) = xVecI(xMax, 1);
    RI = [RI, xVecI];
    xVecII = A\RII(:,i); % erster x-Wert, wird neuer r-Vektor
%     xVecII(1, 1) = xVecII(3, 1);
%     xVecII(xMax+2, 1) = xVecII(xMax, 1);
    RII = [RII, xVecII];
end

for i = 2:xMax+1
    wI(i-1,:) = RI(i,2:xMax+1);
    wII(i-1,:) = RII(i,2:xMax+1);
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

