function [Mvect,Area,Pressure,T,x,Astar,throat_position] = geomfunc(Pcc,Tcc,gamma,option)
%
% funzione calcolo geometria ugello RL10
%
% opzioni: 
% 1) geometria => da file NASA
% 2) approccio lineare con dati noti a priori => file NASA
%
% calcolo: 
%   M
%   A
%   P
%   T
%   x -> discretizzazione lunghezza assiale ugello
%   Astar
%   throat_position -> posizione gola
%
% tutti i valori calcolati sono relazionati a posizione ed area
%

    if option == 1
        
        % TUTTI I DATI ASSEGNATI SONO PRESI DA FILE NASA RL10
        
        geomdata = importdata('RvsX_299pt_cm.txt'); % dati (X   R)
        geomdata = geomdata/1e+2; % trasformazione in metri
        
        [Rstar,throat_position] = min(geomdata(:,2)); % trova gola e la posizione nel vettore
        
        x = geomdata(:,1);
        % divisione dell'analisi in 2 sezioni = conv+div
        A_vec_conv = geomdata(1:throat_position,2); 
        A_vec_div  = geomdata(throat_position+1:end,2);
        % calcolo effettivo area -> aggiornamento vettori
        A_vec_conv = A_vec_conv.^2 *pi; % area convergente
        A_vec_div  = A_vec_div.^2 * pi; % area divergente
        Astar      = Rstar^2 * pi; 
    
        figure(24)
        plot(geomdata(:,1),geomdata(:,2),'LineWidth',2)
        hold on
        grid on
        grid minor
        title('GEOMETRIA ESATTA')
        xlabel('X [m]')
        ylabel('R [m]')
        
    elseif option == 2

        Astar = 0.0119;
        Lcc   = 0.0596;     % injector -> parte cilindrica
        epsilon_ratio = 61; % Ae/Astar
        
        %studio convergente
        L_conv = 0.308-Lcc; % convergente ugello
        Rstar  = sqrt(Astar/pi);
        Din    = 0.273;
        Rin    = Din/2;
        
        %studio divergente
        Ae    = Astar * epsilon_ratio;
        Re    = sqrt(Ae/pi);
        L_div = 1.182;
                
        % scrittura in parti del vettore che descrive la geometria
        x = [linspace(0,Lcc,20)';linspace(Lcc,L_conv+Lcc,100)';linspace(L_conv+Lcc,L_conv+Lcc+L_div,200)'];
        R = [Rin*ones(20,1);flip(linspace(Rstar,Rin,1e+2))';linspace(Rstar,Re,200)'];
        x = x - Lcc-L_conv;
        A = R.^2 * pi;
        
        [~,throat_position] = min(A);
        
        A_vec_conv = A(1:throat_position);
        A_vec_div  = A(throat_position+1:end);
        
        figure(24)
        plot(x,R,'LineWidth',2);
        hold on
        title('GEOMETRIA LINEARE');
        grid on
        grid minor
        title('parete ugello isoentropico -> lineare con dati file NASA')
        xlabel('X [m]')
        ylabel('R [m]')
        
    else
        error('incorrect input');
    end 

    area_vs_mach = @(M,A) A/Astar - 1./M .*(2/(gamma+1)*(1+(gamma-1)/2.*M.^2))^((gamma+1)/(2*(gamma-1)));
    
    Mach_vector_conv = zeros(length(A_vec_conv),1);

    %studio A_vec_conv
    %partenza da: 
    a = 0.001; 
    b = 1;
    for j = 1:length(A_vec_conv)
        Mach_vector_conv(j) = bisez(a,b,area_vs_mach,A_vec_conv(j));
    end

    Mach_vector_div = zeros(length(A_vec_div),1);

    %studio A_vec_div
    %partenza da:
    a = 1;
    b = 10;
    for j = 1:length(A_vec_div)
        Mach_vector_div(j) = bisez(a,b,area_vs_mach,A_vec_div(j));
    end
    
    % immagazzinamento di tutti i dati in vettori
    Mvect    = [Mach_vector_conv; 1; Mach_vector_div];
    T        = Tcc ./ (1 + (gamma-1)/2 * Mvect.^2);                    % temperature relative
    Pressure = Pcc ./ (1 + (gamma-1)/2 * Mvect.^2).^(gamma/(gamma-1)); % pressioni relative   -> !!!Pa!!!
    Area     = [A_vec_conv;Astar;A_vec_div];
    x        = [x;x(end)];
    
    % PLOT
    figure(25)
    plot(Mvect, Area./Astar, 'LineWidth',2);
    hold on
    grid on
    grid minor
    xlabel('Mach number');
    ylabel('A/A*');
    title('Area(Mach)');
    figure(26)
    plot(Mvect,Pressure/Pcc, 'LineWidth',2);
    hold on
    grid on
    grid minor
    xlabel('Mach number');
    ylabel('P/P_{0}');  
    title('P(M)');
    figure(27)
    plot(x,Pressure/Pcc,'LineWidth',2);
    grid on
    grid minor
    xlabel('X [m]');
    hold on
    plot(x,T/Tcc,'LineWidth',2);
    legend('Pressure','Temperature');
    
end