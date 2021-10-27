! programma creazione caso CEA 
!
! 1) il programma legge i dati scritti su un file .txt -> DATAfile.txt generato da MATLAB 
!
! 2) i dati di pressione e temperatura nella posizione i-esima dell'ugello vengono scritti su
!       un file .inp -> RL10nozzlecase.inp
!
! 3) il file .inp detta al CEA di generare un caso di studio
!       problema temperature pressione
!
! 4) l'output del CEA verrá contenuto in RL10nozzlecase.out e RL10nozzlecase.plt          
!       out -> Cp K mu Pr GAMMA    
!  
! 5) il file .plt verrá utilizzato per lo studio dei coefficienti di scambio termico
!
program make_case
    implicit none

    real :: T, P, mH2, mO2 
    
    open(unit=1, file='DATAfile.txt',       status='old', action='read')
    open(unit=2, file='RL10nozzlecase.inp', status='replace')

    read(1,*) mH2
    read(1,*) mO2
    read(1,*) T
    read(1,*) P

    write(2,*) 'problem'   
    write(2,*) '    tp    t,k=',T,'    p,bar=', P,','
    write(2,*) 'react'  
    write(2,*) '    fuel=H2 wt=', mH2,'   t,k=', T  
    write(2,*) '    oxid=O2 wt=', mO2,'   t,k=', T  
    write(2,*) 'output  transport' 
    write(2,*) '    plot cp cond vis pran gamma'  
    write(2,*) 'end'

    close(1)
    close(2)

end program make_case