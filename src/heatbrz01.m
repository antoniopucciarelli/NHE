%% REFERENCES  
%   MECHANICS & THERMODYNAMICS OF PROPULSION, HILL & PETERSON (+++) [EQUATIONS]
%   LIQUID PROPELLANTS ROCKETS, ALTMAN (+)                          [REFERENCE TEMPERATURES]
%   FILE NASA -> ARTICLE - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19970010379.pdf - , p. 131-132
%
% MODELLO SCAMBIO TERMICO UGELLO RL10 TRAMITE MODELLO SEMIEMPIRICO ACCURATO BARTZ

%% CALCOLO GEOMETRIA E DATI
% selezione tipo di analisi tramite option:
%   1) ugello RL10    => geom file NASA
%   2) ugello lineare => geom agli estremi dal file NASA
%

close all
% option
option = 1;

% importazione directories
addpath(..\data\P&Wdata)

% proprietá camera di combustione
Pc    = 3.278e+6;
Tc    = 3.222e+3;
gamma = 1.2105;

% chiamata geomfunc per generazione geometria -> (P,T,M) in relazione ad (x,A)
[Mvect,Area,Pressure,T,xvec,Astar,throat_position] = geomfunc(Pc,Tc,gamma,option);      

% chiamata coolinggeom per generazione geometria tubi concorde con andamento
% dei punti del nozzle
[coolingarea] = coolinggeom(option,xvec,Area);

% modello matematico scambio termico H2 supercritico -> (11.35)
hl = @(G,D,mu_b,cp_b,cp,k_b) 0.023 * G*cp* (G*D/mu_b)^(-0.2) * (mu_b*cp_b/k_b)^(-0.67);                          % (1)
% modello semiempirico scambio termico H2+O2 in ugello -> (11.38)
hg    = @(sigma, A, D, mu0, Pr, p0, cstar, cp0, rc, Astar) (0.026/D^0.2 * ...
    mu0^0.2*cp0/Pr^0.6 * (p0/cstar)^0.8 * (D/rc)^0.1) * (Astar/A)^0.9 * sigma;                                   % (2)
% sigma => variabile per il calcolo di hg
sigma = @(gamma, M, Twh, T0) 1/((0.5 * Twh/T0 * (1+(gamma-1)/2 * M^2) + 0.5)^(-0.4) * (1+(gamma-1)/2* M^2)^1.2); % (3)
rc    = 0.4; % letteratura -> raggio di curvatura gola [ref. Mechanics & thermodynamics of propulsion, p. 530]

% modello trasporto calore 
Q = @(T0, TL, hg, hl, kw, dL) (T0-TL)/(1/hg + dL/kw + 1/hl);

% definizione proprietá wall -> acciaio
kw_steel = 46.6;  % conducibilità termica [W/m^2]
dL       = 0.001; % spessore parete       [m] --> immagini brazing/tubi NASA

% definizione portata
mH2 =  3.5; % [kg/s]
mO2 = 17.5; % [kg/s]

len    = length(xvec);                % dummy variable
cstar  = Astar*Pressure(1)/(mO2+mH2); % utilizzato in calcolo HG
Dbartz = sqrt(Astar/pi)*2;            % modello BARTZ -> D* nella (2) | (11.38)

% definizione vettori raccolta risultati
Qvec_brz  = zeros(len-1,1); % vettore flusso termico
Tmatr_brz = zeros(len-1,3); % matrice temperature [wall,wall H2,H2]

% matrici di riscontro per eventuali incongrunenze/errori/studio andamento proprietá
Hvec_brz  = zeros(len-1,2); % matrice riscontro dati
H1vec_brz = zeros(len-1,1); % matrice riscontro dati

%
% scelta pressioni idrogeno in scambiatore
% determinare option_p per aprire file con relativi dati
% tabelle dati disponibili 
% dati ricavati da REFPROP
% opzione option_p:
%   1) 50bar
%   2) 55bar
%   3) 60bar
%   4) 65bar
%   5) 70bar
%
option_p = 3;

dataFILE = ["hydrogen_data50bar.txt"; ...
            "hydrogen_data55bar.txt"; ...
            "hydrogen_data60bar.txt"; ...
            "hydrogen_data65bar.txt"; ...
            "hydrogen_data70bar.txt"];
P_option = dataFILE(option_p);

% caricamento dati proprietá H2 supercritico T(30K->300K) 
% terza colonna dati -> CpH2 fase supercritica
H2data = importdata(P_option);

TL  = 240; % uscita H2 -> T=240K                      -> file NASA
Twh = 600; % temperatura parete camera di combustione -> file NASA

% waitbar di comoditá per andamento stato codice
w = waitbar(0,'0 %');  

% temperatura di film che viene adottata in scambio di calore convettivo
TrefH2 = (TL + 60)/2; % 60K é la temperatura di ingresso del liquido -> file NASA

% calcolo proprietá temperatura di riferimento LH2 -> da usare poi in (1) | (11.35)
[cp0H2,k0H2,visc0H2] = FINDH2data(TrefH2,H2data.data);
% calcolo proprietá miscela assolute -> [.] in (2) | (11.38)
% !!!imponendo trasformazione adiabatica, i termini 0 (assoluti) non variano
% e di conseguenza anche il valore in [.] della (2) | (11.38)
[cp0,visc_misc0,k_misc0,pr0,gamma0] = CEAcase(mH2,mO2,T(1),Pressure(1)*1e-5);

% ciclo studio ugello da entrata ad uscita 
for i = 1:len-1
     
    % calcolo dell'unico valore necessario per descrivere sigma (2) | (11.38)
    [~,~,~,~,gamma] = CEAcase(mH2,mO2,T(i),Pressure(i)*1e-5);
    
    % ciclo ricerca valore Cp idrogeno in jacket da tabella 
    [cpH2,k_H2,visc_H2] = FINDH2data(TL,H2data.data);
    
    % definizione geometria tubo alla posizione i-esima
    DH2 = 2*sqrt(coolingarea(i,2)/pi);
    GH2 = mH2/(180*coolingarea(i,2));
    
    % chiamata funzioni scambiatore di calore
    % hl -> determina il coeff di scambio termico del jacket
    HL = hl(GH2,DH2,visc0H2,cp0H2,cpH2,k0H2); 
    
    % calcolo simga in (11.38)
    SIGMA = sigma(gamma, Mvect(i), Twh, T(1));
    % calcolo hl -> (11.38)
    HG = hg(SIGMA, Area(i), Dbartz, visc_misc0, pr0, Pressure(1), cstar, cp0, rc, Astar);
    
    % salvataggio dati di riscontro errori
    Hvec_brz(i,:)  = [HL, HG];                    % consente di controllare i 2 coeff che descrivono lo scambio termico -> main values
    H1vec_brz(i,:) = (1/HG + dL/kw_steel + 1/HL); % denominatore (11.30) 
    
    % calcolo scambio termico
    Qp = Q(T(1), TL, HG, HL, kw_steel, dL);
    
    % calcolo area di scambio termico, posizione i-esima
    dA = abs((xvec(i+1)-xvec(i)))* 2*sqrt(Area(i)/pi) * pi; 
    
    % calcolo temperature 
    Twh = T(1) - Qp/HG;            % temperatura a parete
    Twc = Twh  - Qp*dL/kw_steel;   % temperatura a contatto con liquido
    TL  = TL   - Qp*dA/(mH2*cpH2); % temperatura liquido
    % le temperatura del liquido viene espressa come variazione di
    % temperatura dovuta al flusso di calore Qp attraverso la parete
    % laterale dell'ugello
    
    % salvataggio dati temperature    [K]
    Tmatr_brz(i,:) = [Twh,Twc,TL];
    % salvataggio dati flusso termico [W/m^2]
    Qvec_brz(i)    = Qp;
    
    % aggiornamento waitbar
    waitbar(i/(len-1),w,sprintf('%d %%', floor(i/(len-1) * 100)));
   
end
delete(w);
cd DATA
% salvataggio dati risultati
save('xvec.mat','xvec')
save('Q_brz.mat','Qvec_brz');
save('Tmatr_brz.mat','Tmatr_brz')
save('Hvec_brz.mat','Hvec_brz')
cd ..
%% PLOT & SAVING
% !!!nel caso in cui si bloccasse il programma!!!
%    il programma é diviso in sezioni -> é possibile fare il plot dei dati calcolati ed immagazzinati fino al blocco del programma!!!
%
% ogni plot é seguito dal salvataggio del grafico in formato .jpg
% i file verranno salvati in directories /FIGURES_**
%
% !!! ogni run del programma sovrascrive le .jpg
% 
a = length(Qvec_brz);
x = xvec(1:a);
switch(P_option) 
    case dataFILE(1)   
        cd FIGURES_50
    case dataFILE(2)   
        cd FIGURES_55
    case dataFILE(3)   
        cd FIGURES_60
    case dataFILE(4)   
        cd FIGURES_65
    case dataFILE(5)   
        cd FIGURES_70
    otherwise
        cd FIGURES
end
figure(31)
plot(x/abs(x(end)-x(1)),Qvec_brz,'LineWidth',2)
hold on
title('Q^{\prime}')
xlabel('X')
ylabel('Q^{\prime} [W/m^2]')
grid on
grid minor
yticks(0:2e+6:5e+7)
saveas(31,'Q_brz','jpg')
figure(32)
plot(x/abs(x(end)-x(1)),Tmatr_brz(:,1)-273,'LineWidth',2)
hold on
title('T_{wh} BARTZ')
yticks(-200:50:1500)
xlabel('X')
ylabel('T [°C]')
grid on
grid minor
saveas(32,'Tchamber_brz','jpg')
figure(33)
plot(x/abs(x(end)-x(1)),Tmatr_brz(:,2),'LineWidth',2)
hold on
title('T_{wc} BARTZ')
xlabel('X')
ylabel('T [K]')
grid on 
grid minor
yticks(50:25:900)
saveas(33,'Twall_brz','jpg')
figure(34)
plot(x/abs(x(end)-x(1)),Tmatr_brz(:,3),'LineWidth',2)
hold on
title('T_{H_{2}}')
xlabel('X')
ylabel('T [K]')
axis('auto')
yticks(50:5:300)
grid on
grid minor
saveas(34,'TH2_brz','jpg')
figure(35)
plot(H2data.data(:,1),H2data.data(:,3),'LineWidth',2)
hold on
xlabel('T [K]')
ylabel('Cp [kJ/kg K]')
title('Cp_{H_{2}} supercritic')
yticks(0:0.5:20)
xticks(0:25:300)
grid on
grid minor
saveas(35,'CPH2_brz','jpg')
cd ..