%kinetic rate constants fitted in Berkeley Madonna
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

%read in reaction data
dataconc = readtable('LS-P01-066-conc-ramps.csv');
datadiff = readtable('LS-P01-066-different-excess-all.csv');
datamulti = readtable('LS-P01-066-multi-variate-ramps.csv');
dataSS = readtable('LS-P01-066-SS-points.csv');

%solve microkinetic model ODEs across a grid for response surface construction
tspan = linspace(0, 21, 10000);                          % Create constant tspan 
ratio = 0.66:0.02:1.5;                                   % Vector of ratio values
ya2 = zeros(numel(tspan), numel(ratio));                 % Preallocate
for n = 1:numel(ratio)
    ratio1 = ratio(n);                                   % ratio1 is the single value version of ratio which is [1]/[2]
    perB = (100*ratio1)/(1+ratio1);
    A0 = (100-perB)/100*0.4;
    B0 = perB/100*0.4;
[t,ya] = ode45(@(t,y) odefcn(t,y,k1f,k1r,k2f,k2r,k3f,k3r,k4,k5,k6,k7,k8,k9), tspan, [0 A0 B0 0 0 0 0]);
ya2(:,n) = ya(:,1);                                      % this selects pyrazole concentrations to plot (1:P,2:A,3:B,4:M1,5:M2,6:D,7:H) 
end

%plotting
figure
surf(t,ratio,ya2.','FaceAlpha',0.7)                      %.' transposes ya2 to a vector
hold on
grid on
plot3( datadiff.tau, datadiff.ratio, datadiff.pyr, 'ro', 'MarkerFaceColor','#D9FFFF')
plot3( dataconc.tau, dataconc.ratio, dataconc.pyr, 'bo', 'MarkerFaceColor','#D9FFFF')
plot3( datamulti.tau, datamulti.ratio, datamulti.pyr, 'mo', 'MarkerFaceColor','#D9FFFF')
plot3( dataSS.tau, dataSS.ratio, dataSS.pyr, 'k*', 'MarkerSize',20)

set(gca,'FontSize',20)
xlabel('tau / min', 'FontSize',24)
ylabel('[1a]/[2]', 'FontSize',24)
zlabel('[pyrazole 4a] / M', 'FontSize',24)
shading('interp')
hold off

legend('','Residence time ramps','Reactant stoichiometry ramps','Multivariate ramps','Steady state points', 'Location','southeast')

%Microkinetic model function
function dydt = odefcn(t,y,k1f,k1r,k2f,k2r,k3f,k3r,k4,k5,k6,k7,k8,k9)
dydt = zeros(7,1);                                        %1:P,2:A,3:B,4:M1,5:M2,6:D,7:H
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