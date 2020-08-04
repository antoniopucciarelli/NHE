%% REFERENCES  
%   MECHANICS & THERMODYNAMICS OF PROPULSION, HILL & PETERSON (+++) [EQUATIONS]
%   LIQUID PROPELLANTS ROCKETS, ALTMAN (+)                          [REFERENCE TEMPERATURES]
%   FILE NASA -> ARTICLE - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19970010379.pdf - , p. 131-132
%
% MODELLO SCAMBIO TERMICO UGELLO RL10 TRAMITE MODELLO SEMIEMPIRICO SEMPLIFICATO BARTZ

%% CALCOLO GEOMETRIA E DATI
% selezione tipo di analisi tramite option:
%   1) ugello RL10    => geom file NASA
%   2) ugello lineare => geom agli estremi dal file NASA
%

close all
% option
option = 1;

% proprietá camera di combustione => da file riugello.m
Pc    = 3.278e+6;
Tc    = 3.222e+3;
gamma = 1.2105;

% chiamata geomfunc per generazione geometria -> (P,T,M) in relazione ad (x,A)
[Mvect,Area,Pressure,T,xvec,Astar,throat_position] = geomfunc(Pc,Tc,gamma,option);

% chiamata coolinggeom per generazione geometria tubi concorde con andamento
% dei punti del nozzle
[coolingarea] = coolinggeom(option,xvec,Area);

% modello matematico scambio termico H2 supercritico -> (11.35)
hl = @(G,D,mu_b,cp_b,cp,k_b) 0.023 * G*cp* (G*D/mu_b)^(-0.2) * (mu_b*cp_b/k_b)^(-0.67); % (1)

% modello matematico meno accurato bartz --> utilizzato in questo programma -> (11.36)
hg = @(G,D,cp,mu_b,Pr) 0.023* G*cp *(G*D/mu_b)^(-0.2) * Pr^(-0.67); % (2)
% modello semiempirico parte supersonica -> (11.37) 
% il modello (11.37), non viene utilizzato da subito ma dopo 30 iterazioni
% dopo la gola
hgssc = @(a,k,L,G,mu,cp,Pr) k/L * a * (G*L/mu)^0.8 * Pr^0.33; % (3)

% modello trasporto calore -> (11.32)
% il modello trascura la parte di irragiamento
Q = @(T0, Tl, hg, hl, kw, dL) (T0-Tl)/(1/hg + dL/kw + 1/hl);

% definizione proprietá wall -> acciaio
kw_steel = 46.6;  % conducibilità termica [W/m^2]
dL       = 0.001; % spessore parete       [m] --> immagini brazing/tubi NASA

% definizione portata
mH2 =  3.5; % [kg/s]
mO2 = 17.5; % [kg/s]

len = length(xvec); % variabile dummy
% definizione vettori raccolta risultati
Qvec   = zeros(len-1,1); % vettore flusso termico
Tmatr  = zeros(len-1,3); % matrice temperature [wall,wall-H2,H2]
Cp     = zeros(len-1,1); % vettore Cp miscela
K_f    = zeros(len-1,1); % vettore k miscela
Visc_f = zeros(len-1,1); % vettore visc miscela
Pran_f = zeros(len-1,1); % vettore pr miscela

% matrici di riscontro dati -> da analizzare per controllare lo scambio
% termico
Hvec  = zeros(len-1,2); % matr riscontro dati
H1vec = zeros(len-1,1); % matr riscontro dati

%
% scelta pressioni idrogeno in scambiatore
% determinare pressure_p per aprire file con relativi dati
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

% chimata programma per descrive le proprietá di riferimento di LH2 nello
% scambio termico -> media Tingresso + Tuscita
TrefH2 = (TL + 60)/2; % 60K é la temperatura di ingresso del liquido -> file NASA
[cp0H2,k0H2,visc0H2] = FINDH2data(TrefH2,H2data.data);

% waitbar di comoditá per andamento stato codice
w = waitbar(0,'0 %');

% r => letteratura cap11 -> valore di riferimento sperimentale
% descrizione temperatura parete camera vs temp totale di stagnazione
r = 0.91; % recovery factor  

% ciclo studio ugello da entrata ad uscita 
for i=1:len-1
    
    % calcolo con NASACEA le proprietá del fluido durante espansione
    % adiabatica nella sezione i-esima
    [cp,visc,k_misc,pr] = CEAcase(mH2,mO2,T(i),Pressure(i)*1e-5);
    
    % il valore di cp miscela viene immagazzinato in un vettore per un
    % possibile riscontro
    Cp(i) = cp;  
    
    % REFERENCE TEMPERATURE  
    % dalla letteratura -> parte supersonica e prima della gola descrivibile 
    % con formula (1)
    % parte iniettori - gola con (2)
    %
    % scelte: si é utilizzato il modello (1) 30 iterazioni prima della gola
    if i>(throat_position-30)
        % --------(1)----------
        Tr = r*(T(1)-T(i)) + T(i);
        Tf = Twh + 0.23*(T(i)-Twh) + 0.19*(Tr-Twh);
    else
        % --------(2)----------
        Tf = (T(i)+Twh)/2;
    end
    % REFERENCE TEMPERATURE PROPERTIES      
    % calcolo tramite NASACEA delle proprietá della miscela rispetto Pressione
    % i-esima e Temperatura di film/reference
    % gli output sono solo quelli necessari per il ciclo 
    [~,visc_misc_f,~,pr_f] = CEAcase(mH2,mO2,Tf,Pressure(i)*1e-5);
    
    % salvataggio valori in matrice di studio -> eventuali problemi
    Visc_f(i) = visc_misc_f; 
    Pran_f(i) = pr_f;
    
    % ciclo ricerca valore Cp idrogeno in jacket da tabella 
    [cpH2] = FINDH2data(TL,H2data.data);
    
    % calcolo portate e diam relativi alla miscela propellente 
    % (3) varia con posizione i-esima
    Dmiscela = 2*sqrt(Area(i)/pi); % (3)
    Gmiscela = (mO2+mH2)/Area(i);

    % definizione geometria tubo alla posizione i-esima
    % (4) varia con posizione i-esima
    DH2 = 2*sqrt(coolingarea(i,2)/pi); % (4)
    GH2 = mH2/(180*coolingarea(i,2));
    
    % chiamata funzioni scambiatore di calore
    % hl -> determina il coeff di scambio termico del jacket
    HL = hl(GH2,DH2,visc0H2,cp0H2,cpH2,k0H2); 
    
    % studio coefficiente di scambio convettivo in base alla sezione di
    % riferimento
    % per flusso subsionico e per 30 iterazioni dopo gola si é scelta (2)
    % per flusso supersonico (3)
    if i<throat_position+30
       % hg -> determina il coeff di scambio termico della miscela    
       HG = hg(Gmiscela,Dmiscela,cp,visc_misc_f,pr_f); % (2)
    else 
       % (3) funziona con distanza dalla gola -> espressa con L 
       L = abs(xvec(i) - xvec(throat_position));
       % hg -> determina il coeff di scambio termico della miscela    
       HG = hgssc(0.028,k_misc,L,Gmiscela,visc,cp,pr);
    end
    
    % calcolo scambio termico
    Qp = Q(T(1), TL, HG, HL, kw_steel, dL);
    
    % calcolo area di scambio termico, posizione i-esima
    % !!! si utilizza abs(.) in quanto abbiamo che xvec é un vettore (>0 .and. <0)
    dA = abs((xvec(i+1)-xvec(i)))* 2*sqrt(Area(i)/pi) * pi; 
    
    % calcolo temperature 
    Twh = T(1) - Qp/HG;            % temperatura a parete
    Twc = Twh  - Qp*dL/kw_steel;   % temperatura a contatto con liquido
    TL  = TL   - Qp*dA/(mH2*cpH2); % temperatura liquido modificata a TL
    % le temperatura del liquido viene espressa come variazione di
    % temperatura dovuta al flusso di calore Qp attraverso la parete
    % laterale dell'ugello
    
    % salvataggio matrici di riscontro
    Hvec(i,:)  = [HL, HG];                    % consente di controllare i 2 coeff che descrivono lo scambio termico -> main values
    H1vec(i,:) = (1/HG + dL/kw_steel + 1/HL); % denominatore (11.30) 
    
    % salvataggio valori temperature [camera, contatto parete fluido, fluido supercritico]
    Tmatr(i,:) = [Twh,Twc,TL];
    % salvataggio valori flusso termico [W/m^2]
    Qvec(i)    = Qp;
    
    % aggiornamento waitbar
    waitbar(i/(len-1),w,sprintf('%d %%', floor(i/(len-1) * 100)));

end
delete(w);
cd DATA
% salvataggio dati
save('xvec.mat','xvec')
save('Q_std.mat','Qvec');
save('Tmatr_std.mat','Tmatr')
save('Hvec_std.mat','Hvec')
save('CpMIX_std.mat','Cp')
save('ViscMIX_std.mat','Visc_f')
save('PrMIX_std.mat','Pran_f')
cd ..
%% PLOT & SAVING
% !!!nel caso in cui si bloccasse il tutto!!!
%    il programma é diviso in sezioni -> é possibile fare il plot dei dati calcolati ed immagazzinati fino al blocco del programma!!!
%
% ogni plot é seguito dal salvataggio del grafico in formato .jpg
% i file verranno salvati in directories /FIGURES_**
%
% !!! ogni run del programma sovrascrive le .jpg
% 
a = length(Qvec);
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
plot(x/abs(x(end)-x(1)),Qvec,'LineWidth',2);
hold on
title('Q^{\prime}')
xlabel('X')
ylabel('Q^{\prime} [W/m^2]')
grid on
grid minor
yticks(0:2e+6:5e+7)
saveas(31,'Q_std','jpg')
figure(32)
plot(x/abs(x(end)-x(1)),Tmatr(:,1)-273,'LineWidth',2)
hold on
title('T_{wh}')
yticks(-200:50:1500)
xlabel('X')
ylabel('T [°C]')
grid on
grid minor
saveas(32,'Tchamber_std','jpg')
figure(33)
plot(x/abs(x(end)-x(1)),Tmatr(:,2),'LineWidth',2)
hold on
title('T_{wc}')
xlabel('X')
ylabel('T [K]')
grid on 
grid minor
yticks(100:25:900)
saveas(33,'Twall_std','jpg')
figure(34)
plot(x/abs(x(end)-x(1)),Tmatr(:,3),'LineWidth',2)
hold on
title('T_{H_{2}}')
xlabel('X')
ylabel('T [K]')
axis('auto')
yticks(50:5:250)
grid on
grid minor
saveas(34,'TH2_std','jpg')
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
saveas(35,'CPH2_std','jpg')
figure(36)
plot(x/abs(x(end)-x(1)),Cp/1e+3,'LineWidth',2)
hold on
grid on
grid minor
title('Cp propellente OF5')
xlabel('X')
ylabel('Cp [kJ/kg K]')
saveas(36,'CPMIX_std','jpg')
figure(37)
plot(x/abs(x(end)-x(1)),Visc_f,'LineWidth',2)
hold on
grid on
grid minor
title('\mu REFERENCE TEMPERATURE');
saveas(37,'muMIX_std','jpg')
figure(38)
plot(x/abs(x(end)-x(1)),Pran_f,'LineWidth',2)
hold on
grid on
grid minor
title('Prandtl REFERENCE TEMPERATURE')
saveas(38,'Pr_std','jpg')
cd ..