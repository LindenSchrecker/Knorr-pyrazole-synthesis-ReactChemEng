%DM 
k1f = 2.23589;
k1r = 0;
k2f = 1753.989;
k2r = 0.25614;
k3f = 13860.73;
k3r = 0;
k4 = 0.63061;
k5 = 0;
k6 = 0.0399;
k7 = 11.87965;
k8 = 0;
k9 = 3.69715;

tspan = linspace(0, 21);                                 % Create Constant ‘tspan’
ratio=0.66:0.02:1.5;                                     % Vector of ratio values
ya2 = zeros(numel(tspan), numel(ratio));                 % Preallocate
for n = 1:numel(ratio)
    ratio1 = ratio(n);                                   % ratio1 is the single value version of ratio which is [1]/[2]
    perB = (100*ratio1)/(1+ratio1);
    A0 = (100-perB)/100*0.4;
    B0 = perB/100*0.4;
[t,ya] = ode45(@(t,y) odefcn(t,y,k1f,k1r,k2f,k2r,k3f,k3r,k4,k5,k6,k7,k8,k9), tspan, [0 A0 B0 0 0 0 0]);
ya2(:,n) = ya(:,1);                                      % this selects pyrazole concentrations to plot (1:P,2:A,3:B,4:M1,5:M2,6:D,7:H) 
end

dataconcme = readtable('LS-P01-066-conc-ramps.csv');
datadiffme = readtable('LS-P01-066-different-excess-all.csv');
datamultime = readtable('LS-P01-066-multi-variate-ramps.csv');



%DE parameters
k1f = 6.22879;
k1r = 0;
k2f = 81.44262;
k2r = 0.01174;
k3f = 551.6668;
k3r = 0;
k4 = 0.74228;
k5 = 0;
k6 = 0.44781;
k7 = 3.26931;
k8 = 0;
k9 = 4.97573;

tspan = linspace(0, 21);                                 % Create Constant ‘tspan’
ratio=0.66:0.02:1.5;                                     % Vector of ratio values
yb2 = zeros(numel(tspan), numel(ratio));                 % Preallocate
for n = 1:numel(ratio)
    ratio1 = ratio(n);                                   % ratio1 is the single value version of ratio which is [1]/[2]
    perB = (100*ratio1)/(1+ratio1);
    A0 = (100-perB)/100*0.4;
    B0 = perB/100*0.4;
[t,yb] = ode45(@(t,y) odefcn(t,y,k1f,k1r,k2f,k2r,k3f,k3r,k4,k5,k6,k7,k8,k9), tspan, [0 A0 B0 0 0 0 0]);
yb2(:,n) = yb(:,1);                                      % this selects pyrazole concentrations to plot (1:P,2:A,3:B,4:M1,5:M2,6:D,7:H) 
end

dataconcet = readtable('LS-P01-072-conc-ramps.csv');
datadiffet = readtable('LS-P01-069-different-excess-all.csv');


figure;
hold on
surf(t,ratio,ya2.','FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor', 'none')
plot3( dataconcme.tau, dataconcme.ratio, dataconcme.pyr, 'ro', 'MarkerFaceColor','#FFFFFF', 'MarkerSize',5)
plot3( datadiffme.tau, datadiffme.ratio, datadiffme.pyr, 'ro', 'MarkerFaceColor','#FFFFFF', 'MarkerSize',5)
plot3( datamultime.tau, datamultime.ratio, datamultime.pyr, 'ro', 'MarkerFaceColor','#FFFFFF', 'MarkerSize',5)

surf(t,ratio,yb2.','FaceColor','b','FaceAlpha',0.5, 'EdgeColor', 'none')
plot3( dataconcet.tau, dataconcet.ratio, dataconcet.pyr, 'bo', 'MarkerFaceColor','#FFFFFF', 'MarkerSize',5)
plot3( datadiffet.tau, datadiffet.ratio, datadiffet.pyr, 'bo', 'MarkerFaceColor','#FFFFFF', 'MarkerSize',5)

hold off
grid on

set(gca,'FontSize',20)
xlabel('tau / min', 'FontSize',24)
ylabel('[1]/[2]', 'FontSize',24)
zlabel('[pyrazole 4] / M', 'FontSize',24)

hold off

legend('Surface R = Me','Data R = Me','','','Surface R = Et','Data R = Et', 'Location','southeast')

function dydt = odefcn(t,y,k1f,k1r,k2f,k2r,k3f,k3r,k4,k5,k6,k7,k8,k9)
dydt = zeros(7,1);                                       %1:P,2:A,3:B,4:M1,5:M2,6:D,7:H
RXN15 = k1f*y(2)*y(3)-k1r*y(4)*y(7);
RXN25 = k2f*y(4)-k2r*y(5);
RXN35 = k3f*y(4)*y(2)-k3r*y(6)*y(7);
RXN45 = k4*y(5)*y(3);
RXN55 = k5*y(5);
RXN65 = k6*y(5)*y(1);
RXN75 = k7*y(6)*y(3);
RXN85 = k8*y(6);
RXN95 = k9*y(6)*y(1);

dydt(1) = RXN45+RXN55+RXN65+RXN75+RXN85+RXN95;
dydt(2) = -RXN15-RXN35+RXN75+RXN85+RXN95;
dydt(3) = -RXN15;
dydt(4) = RXN15-RXN25-RXN35;
dydt(5) = RXN25-RXN45-RXN55-RXN65;
dydt(6) = RXN35-RXN75-RXN85-RXN95;
dydt(7) = RXN15+RXN35+RXN45+RXN55+RXN65;
end