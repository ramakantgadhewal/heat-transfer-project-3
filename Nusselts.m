Nusselt = table2array(readtable('Nusselt.dat'));

n = numel(Nusselt(:,1));

x = linspace(Nusselt(1,1), Nusselt(n,1), n);

figure(1);
plot(x,Nusselt(:,2), 'LineWidth',1);
title({'Nusselt Number along a Heated Section of Pipe', 'Radial Nodes = 201, Heated Section = 300D, Re = 2000, Pr = 0.7'}, 'FontSize', 16);
xticks(0:20:300);
xlabel('Heated Section in Pipe Diameters, x', 'FontSize', 16);
yticks(0:10:60);
ylabel('Nusselt Number, Nu', 'FontSize', 16);
legend('Nu(x)', 'FontSize', 16);
legend('location','east');
grid on;
grid minor;
set(gca, 'MinorGridLineStyle','-','FontSize', 16);



