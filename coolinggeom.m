function [coolingarea] = coolinggeom(option,x,Area)
%
% programma calcolo andamento area tubi di raffreddamento da file
% TUBEgeometry_300pt_cm.txt
%
% il programma relaziona le ascisse con le aree della geometria dei tubi e dell'ugello RL10
%
% il problema sorge dal fatto che la discretizzazione della geometria dell'ugello é
% diversa da quella del file con l'andamento dell'area dei tubi
%
% option:
%   1) geometria reale ugello RL10 => file NASA
%   2) geometria lineare           => file NASA
%
% !!! non si é scelto di utilizzare i dati provenienti dal metodo di RAO
% per la discritizzazione dell'ugello in quanto non esiste relazione
% diametro ugello vs diametro tubi !!!
%

    % import dei dati da file NASA -> geometria aree tubi di raffreddamento
    coolingdata = importdata('TUBEgeometry300pt_cm.txt'); % dati (X  area)
    if option == 1
        nozzledata =  importdata('RvsX_299pt_cm.txt'); % dati (X   R)
        coolingdata = [linspace(nozzledata(1,1),coolingdata(1,1),100)', coolingdata(1,2)*ones(100,1); coolingdata];
        % determino intervallo in cui comparo le ascisse dei tubi e dell'ugello
        interval   = 0.4; % valore massimo di errore che si trova tra l'ordinata di un punto ed il suo successivo
        tubenozzle = zeros(size(nozzledata));

        % ciclo di comparazione ascisse
        for j = 1:length(nozzledata)
            for i = 1:length(coolingdata)
                if (nozzledata(j,1)+interval >= coolingdata(i,1) && nozzledata(j,1)-interval <= coolingdata(i,1))
                    % se le condizioni sono soddisfatte scambio nella posizione
                    % j-esima il valore alla posizione i-esima
                    tubenozzle(j,1) = nozzledata(j,1);
                    tubenozzle(j,2) = coolingdata(i,2);
                    ii = j;
                end
            end
        end

        % siccome si ha a disposizione 299 punti per i tubi puó essere che l'array non sará pieno
        % le posizioni vuote con 0 
        % sostituzione con l'ultimo valore dell'area dei tubi -> errore dovuto ai dati disponibili sui tubi
        tubenozzle(ii:end,1) = nozzledata(ii:end,1);
        tubenozzle(ii:end,2) = tubenozzle(ii,2);

        % PLOT
        figure(21)
        plot(nozzledata(:,1)*1e-2,(nozzledata(:,2)*1e-2).^2*pi,'Linewidth',2)
        hold on
        plot(tubenozzle(:,1)*1e-2,tubenozzle(:,2),'k','LineWidth',2);
        title('ANDAMENTO AREE RL10 REALI');
        grid on
        grid minor
        xlabel('X [m]')
        ylabel('Area')
        legend('RL10 nozzle   [m^{2}]','H2 supercritic [cm^{2}]','Location','best');
    
        tubenozzle(:,2) = tubenozzle(:,2).*1e-4;
        
        coolingarea = tubenozzle;

    elseif option == 2

        interval = 0.4; % valore massimo di errore che si trova tra l'ordinata di un punto ed il suo successivo
        tube     = zeros(length(x),2);
        
        x = x.*1e+2;
       
        coolingdata = [linspace(x(1,1),coolingdata(1,1),100)', coolingdata(1,2)*ones(100,1); coolingdata];
        
        for j = 1:length(x)
            for i = 1:length(coolingdata)
                if (x(j)+interval >= coolingdata(i,1) && x(j)-interval <= coolingdata(i,1))
                    tube(j,1) = x(j);
                    tube(j,2) = coolingdata(i,2);
                    ii = j;
                end
            end
        end
        
        tube(ii:end,1) = x(ii:end);
        tube(ii:end,2) = tube(ii,2);
        
        x = x.*1e-2;
        
        % PLOT
        figure(21)
        plot(x(:),Area,'Linewidth',2)
        hold on
        plot(tube(:,1)*1e-2,tube(:,2),'k','LineWidth',2);
        grid on
        grid minor
        xlabel('X [m]')
        ylabel('Area')
        legend('RL10 linear nozzle [m^{2}]','H2 supercritic         [cm^{2}]','Location','best');
        
        tube(:,2) = tube(:,2) .* 1e-4;
        
        coolingarea = tube;
    
    else

        error('incorrect input')

    end
end