function [cp, visc_misc, k_misc, pr, gamma] = CEAcase(mH2, mO2, T, Pressure)
%
% funzione chiamante NASACEA -> problema tp 
% input  -> portata H2, portata O2, T, Pressione
% output -> Cp, vis, k, Pr, gamma
% 
    
    % analisi del fluido, mi sposto nella directory NASA_CEA
    cd NASA_CEA/

    % scrittura dati relativi miscela in DATAfile
    DATAfile = fopen('DATAfile.txt','w');
    fprintf(DATAfile,'%d\n',mH2);
    fprintf(DATAfile,'%d\n',mO2);
    fprintf(DATAfile,'%d\n',T);
    fprintf(DATAfile,'%d\n',Pressure);
    fclose(DATAfile);
    
    % chiamata programma generazione file RL10nozzlecase.inp 
    % -> descrizione in file make_case.f95
    system('./make_case');
    
    % chiamata programma che esegue automaticamente file RL10nozzlecase.inp
    % -> descrizione in file FCEA.f
    system('./FCEA');
    
    % chiamata programma che trova il relativo valore di Cp miscela da
    % RL10nozzlecase.plt -> file con i dati necessari allo studio
    %   1 cp    [kJ/kg K]
    %   2 k     [mW/cm K] MILLIWATTS/(CM)(K)
    %   3 visc  [mP] MILLIPOISE
    %   4 Pr    [-]
    %   5 gamma [-]
    
    % caricamento dati da output caso CEA
    MISCdata = importdata('RL10nozzlecase.plt');
   
    % CONVERSIONE

    % conversione Cp [J/kg K]
    cp = MISCdata.data(1);
    cp = cp * 1e+3;
 
    % conversione k [mW/cm K] < k > [W/m K]
    k_misc = MISCdata.data(2);
    k_misc = k_misc * 1e-1;
 
    % conversione visc [mP]< visc >[Pa s] 
    visc_misc = MISCdata.data(3);
    visc_misc = visc_misc * 1e-4;
    
    % Pr
    pr = MISCdata.data(4);
    
    % gamma
    gamma = MISCdata.data(5);
    
    % main directory
    cd .. 
    
end