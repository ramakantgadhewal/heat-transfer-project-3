Mean = table2array(readtable('Mean.dat'));

n = numel(Mean(:,1));

x = linspace(Mean(1,1), Mean(n,1), n);

figure(3);
plot(x, Mean(:,2),'LineWidth',1);
title({'Mean Temperature along a Heated Section of Pipe', 'Radial Nodes = 201, Heated Section = 300D, Re = 2000, Pr = 0.7'}, 'FontSize', 16);
xlabel('Heated Section in Pipe Diameters, x', 'FontSize', 16);
xticks(0:20:300);
ylabel('Dimensionless Mean Temperature, \theta_m', 'FontSize', 16);
legend('\theta_m(x)', 'FontSize', 16);
legend('location','east');
grid on;
grid minor;
set(gca, 'MinorGridLineStyle','-','FontSize', 16);