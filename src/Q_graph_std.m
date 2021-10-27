% CODICE NON RILEVANTE --> GESTIONE GRAFICI NELLA APPENDICE 
%
% codice generazione grafici finali | modello standard 
% i grafici sono nelle appendici
%

close all
clc

Qdata = importdata('Heatfluxgrab.txt');

Qdata(:,1) = Qdata(:,1)./(abs(Qdata(1,1)) + abs(Qdata(end,1))); 
Qdata(:,2) = 1.635e+6*Qdata(:,2);

figure(1)
plot(Qdata(:,1),Qdata(:,2), 'k', 'LineWidth', 2);
hold on
plot(x/abs(x(end)-x(1)),Qvec,'b','LineWidth',2)
title('heatstd01.m - Q^{\prime}');
grid on
grid minor
xlabel('X - asse ugello %')
ylabel('Q^{\prime} [W/m^2]')
legend('file NASA', 'heatstd01.m');
saveas(1, 'Qprime_std','jpg')

TdataWall = importdata('WallTempgrab.txt');

TdataWall(:,1) = TdataWall(:,1)./( abs(TdataWall(1,1)) + abs(TdataWall(end,1)));
TdataWall(:,2) = (TdataWall(:,2) - 491)/1.8;

figure(2)
plot(TdataWall(:,1), TdataWall(:,2), 'k', 'LineWidth', 2);
hold on
plot(x/abs(x(end)-x(1)), Tmatr(:,1)-273,'b','LineWidth',2)
yticks(-200:50:1500)
title('heatstd01.m - T_{wh}')
xlabel('X - asse ugello %')
ylabel('T [°C]')
grid on
grid minor
legend('file NASA', 'heatstd01.m');
saveas(2, 'TWall_std','jpg')

TdataH2 = importdata('CoolantTemp.txt');

TdataH2(:,1) = TdataH2(:,1)./( abs(TdataH2(1,1)) + abs(TdataH2(end,1)));
TdataH2(:,2) = (TdataH2(:,2) - 491)/1.8;

figure(3)
plot(TdataH2(:,1), TdataH2(:,2)+272, 'k', 'LineWidth', 2);
hold on
plot(x/abs(x(end)-x(1)),Tmatr(:,3),'b','LineWidth',2)
title('heatstd01.m - H_{2}');
legend('file NASA', 'heatstd01.m');
grid on
grid minor
yticks(20:20:300)
xlabel('X - asse ugello %')
ylabel('T [K]')
saveas(3, 'TH2_std','jpg')