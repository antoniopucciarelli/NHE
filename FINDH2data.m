function [cpH2,k_H2,visc_H2] = FINDH2data(TL,data)
%
% ciclo ricerca valore Cp idrogeno in jacket da tabella scelta tramite
% option_p
%
% dati tabella:
% 1) T
% 2) P
% 3) Cp
% 4) K
% 5) mu
% 
% i dati sono stati ricavati da REFPROP
%

    k    = 1;
    resp = 1;
    toll = 1;
    
    while resp == 1
        a = abs(data(k,1)-TL);
        if(a<=toll)
            cpH2    = data(k,3);
            k_H2    = data(k,4);
            visc_H2 = data(k,5);
            resp    = 0;
        end
        k = k+1;
    end
    
    cpH2 = cpH2*1e+3;
    
end