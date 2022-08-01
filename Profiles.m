Profile = table2array(readtable('Profile.dat'));

rows = numel(Profile(:,1));
columns = numel(Profile(1,:));
n = 2 * rows - 1;

TP = zeros(n, columns);

temp = flip(Profile);

for i=1:rows
    TP(i,:) = temp(i,:);
end

j = 1;
for i=rows:n
    TP(i,:) = Profile(j,:);
    j = j + 1; 
end

figure(2);
r = linspace(-1, 1, n);
plot(TP, r,'LineWidth',1);
title({'Temperature Profiles along a Heated Section of Pipe', 'Radial Nodes = 201, Heated Section = 300D, Re = 2000, Pr = 0.7'}, 'FontSize', 16);
xlabel('Dimensionless Temperature, \theta', 'FontSize', 16);
ylabel('Dimensionless Radius, r', 'FontSize', 16);
legend('60D','120D','180D','240D','300D', 'FontSize', 16);
legend('location','northwest');
grid on;
grid minor;
set(gca, 'MinorGridLineStyle','-','FontSize', 16);