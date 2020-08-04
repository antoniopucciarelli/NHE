function [x] = bisez(a, b, fun, A)
% metodo di bisezione applicato alla funzione relazione M ad Area
% input:
% estermi [a,b]
% funzione area_vs_mach
% A = valore Area della sezione di studio
%

it   = -1;
toll = 1e-4;
err  = 1;
nmax = 100;

% controllo sulla validita' degli estremi
if (fun (a,A) * fun (b,A) > 0)
    error ('function error \n geometry error');
end

while (it <= nmax-1 && err > toll) 
    it = it+1;
    x  = (b+a)/2;    %stima dello zero
    fc = fun(x,A);   
    
    if (fc == 0)
       err = 0;
    else
       err = abs(fc); 
    end    

    % scelta del nuovo estremo per l'eventuale ciclo successivo       
    if (fc*fun(a,A) > 0) 
          a = x; 
    else 
          b = x; 
    end      
end

end