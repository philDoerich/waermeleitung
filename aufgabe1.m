%Variablen
L       = 1;
T       = 6;   %Betrachtungszeit
n       = 100;
K       = 1;
tMax    = 100;
xMax    = 100;
%%%%%%%%%
t       = linspace(0, T, tMax);
x       = linspace(0, L, xMax);
deltaX  = L/xMax;
deltaT  = (1/(2*K))*deltaX^2;
d       = K*deltaT/deltaX^2;
B_nI    = zeros(n, 1);
B_nII   = zeros(n, 1);
B_nIII  = zeros(n, 1);
rI      = zeros(xMax, 1);
vI      = zeros(xMax, tMax);
vII     = zeros(xMax, tMax);
vIII    = zeros(xMax, tMax);
uI      = zeros(xMax, tMax);
uII     = zeros(xMax, tMax);
uIII    = zeros(xMax, tMax);
wI      = zeros(xMax, tMax);
wII     = zeros(xMax, tMax);
wIII    = zeros(xMax, tMax);
A       = zeros(xMax, tMax);
funE    = @(n, t) exp(-((n*pi/L)^2)*K*t);
funSin  = @(n, x) sin((n*pi*x)/L); 
funI    = @(x) 3*sin(pi*x) + 5*sin(4*pi*x);
funII   = @(x) 1;
funIII  = @(x) x;

%Bausteine berechnen
for i = 1:n
    B_nI(i, 1)      = (2/L) .* integral(@(x) funI(x).*funSin(i, x),0,L);
    B_nII(i, 1)     = (2/L) .* integral(@(x) funII(x).*funSin(i, x),0,L);
    B_nIII(i, 1)    = (2/L) .* integral(@(x) funIII(x).*funSin(i, x),0,L);
end

%Exakte Lösung
for j = 1:xMax
    for k = 1:tMax
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

%Anfangswerte
for j = 1:xMax
    vI(j, 1)    = funI(x(j));
    vII(j, 1)   = funII(x(j));
    vIII(j, 1)  = funIII(x(j));
end

%Numerische Lösung (explizit)
for j = 2:xMax-1
    for k = 2:tMax-1
        vI(j,k)     = d*(vI(j+1, k-1) - 2*vI(j, k-1) + vI(j-1, k-1)) + vI(j, k-1);
        vII(j,k)    = d*(vII(j+1, k-1) - 2*vII(j, k-1) + vII(j-1, k-1)) + vII(j, k-1);
        vIII(j,k)   = d*(vIII(j+1, k-1) - 2*vIII(j, k-1) + vIII(j-1, k-1)) + vIII(j, k-1);
    end
end

%Numerische Lösung (implizit)
%for j = 2:xMax-1
    %for k = 2:tMax-1
    %    wI(j,k)     = -d*wI(j-1, k+1) + (1 + 2*d)*wI(j, k+1) - d*wI(j+1, k+1);
   %     wII(j,k)    = -d*wII(j-1, k+1) + (1 + 2*d)*wII(j, k+1) - d*wII(j+1, k+1);
  %      wIII(j,k)   = -d*wIII(j-1, k+1) + (1 + 2*d)*wIII(j, k+1) - d*wIII(j+1, k+1);
 %   end
%end



vec1    = zeros(xMax, 1);
vec2    = zeros(xMax-1, 1);

for j = 1:xMax
    vec1(j, 1) = 1+2*d;
end
for j = 1:xMax-1
    vec2(j, 1) = -d;
end

A1 = diag(vec1);
A2 = diag(vec2, 1);
A3 = diag(vec2, -1);

%A = A1 + A2 + A3;



figure('Name', 'Plot i)','NumberTitle','off')
plot(x, uI(1:xMax, 1));

figure('Name', 'Plot ii)','NumberTitle','off')
plot(x, uII(1:xMax, 1));

figure('Name', 'Plot iii)','NumberTitle','off')
plot(x, uIII(1:xMax, 1));

figure('Name', 'Plot numerical i)','NumberTitle','off')
plot(x, vI);

figure('Name', 'Plot numerical ii)','NumberTitle','off')
plot(x, vII);

figure('Name', 'Plot numerical iii)','NumberTitle','off')
plot(x, vIII);