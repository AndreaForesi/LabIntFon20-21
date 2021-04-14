import argparse
import itertools
import numpy as np
import logging
from scipy.signal import argrelmin,find_peaks
import matplotlib.pyplot as plt
from glob import glob

"""Uso questo script per fare cambiamenti a quello principale e tenere fuznioni
    che al momento non utilizzo"""

RecLen=1030

###     Library     ###
"""-signal"""
def GestioneoffIterativa(channel):
    sigma=11
    l=0
    r=50
    while ( (sigma>10) & (r<RecLen-50) ):

        a=channel[l:r+1]
        offset=np.mean(a)
        sigma=np.std(a)
        l=l+50
        r=r+50

    return offset,sigma
def GestioneoffPrimiCamp(channel):
    """-Si esegue media e devstd dell'offest considerando i primi punti"""
    Dati    =   channel[0:201]
    off     =   np.mean(Dati)
    Sigma   =   np.std(Dati)
    return off,Sigma
def GestioneoffDueSoglie(channel):
    """-Gestione offset"""
    #impostare una soglia per gestire l'offset
    SOGLIASEGNALEDown=14650
    SOGLIASEGNALEUp=14775
    RumoreChannel=np.ma.masked_where((channel>SOGLIASEGNALEUp)|(channel<SOGLIASEGNALEDown),channel)    #ricerco il rumore
    offset=RumoreChannel.mean()    #calcolo
    Doffset=RumoreChannel.std()  #verrà usato per definire quando potremo considerare il segnale
    return offset, Doffset
def TimeStampMeas(Signal,min_pos,sigma):
    """-Time stamp measurement"""

    ###Media negative slope
    """
    negativeslope=np.array([])
    min_val=Signal[min_pos]
    i=min_pos

    while ((Signal[i]<min_val*0.2) & (i>1)):
        if Signal[i]>0.8*min_val:
            negativeslope=np.append(negativeslope,i)
        i-=1

    time_stamp = np.mean(negativeslope)
    """
    ###Primo valore oltre n*sigma
    i=min_pos
    while (Signal[i]<-5*sigma):
        i-=1
    time_stamp = i
    """
    ###media tra alcuni valori attorno a 5 sigma
    time = np.arange(len(Signal))
    i=min_pos
    MAX=True

    while (Signal[i]<-3*sigma):
        if (MAX==True):
            if (Signal[i]<-7*sigma):
                max = i
            else:
                MAX=False
                max=i-1
        min=i
        i-=1

    if (i==min_pos):
        return min_pos

    time_stamp = (min+max)/2
    """

    return time_stamp
def ChargeMeas(Signal):
    """-Charge measurement"""
    Charge=np.abs(np.trapz(Signal))
    #Charge=np.abs(np.min(Signal))
    return Charge
def findstart(Sign_arr,min_pos,sigma,limit):
    """- find the first data position of the negative slope of impulse"""
    if (min_pos-30>=limit):
        i=min_pos-30
    else:
        i=limit
    """
    i=min_pos
    while(Sign_arr[i]<0.2*Sign_arr[min_pos]):
        if i==limit:
            break
        i-=1
    """
    return i
def findstop(Sign_arr,min_pos,sigma,limit):
    """- find the last data position of the negative slope of impulse"""
    if (min_pos+350<=limit):
        i=min_pos+350
    else:
        i=limit
    """
    i=min_pos
    if i==limit: return i
    while(Sign_arr[i]<0.2*Sign_arr[min_pos]):
        if i>=(limit):
            break
        i=i+1
    """
    return i
def EstraiSegnale(channel,min_pos,Sigma,inf,sup):
    ##time array
    time=np.arange(len(channel))
    ##Seleziona la fine e l'inizio del segnale trmite findstart e findstop (minore di 1 sigma)
    start_time=findstart(channel,min_pos,Sigma,inf)
    stop_time=findstop(channel,min_pos,Sigma,sup)
    ##Estrazione segnale start -->to--> stop
    Signal = channel[start_time:(stop_time+1)]
    Signal_time = time[start_time:(stop_time+1)]
    return Signal,Signal_time
"""-Threshold for discrimination of signal"""
def ThrsDet(PMTDataFile):
    import re
    logging.basicConfig(level=logging.INFO)
    Re_Num=re.compile('\d+')
    Re_Record=re.compile('\ARecord')
    Re_EOF=re.compile('\A\Z')
    """-Restituisce una soglia proporzionale al rumore e ai minimi letti dei segnali provenienti da un PMT"""
    minRapp=np.array([])
    ContErr=evento=0
    EOF=False
    data = np.array([])
    with open(PMTDataFile) as file:
        line=file.readline()
        while(not EOF):
            line=file.readline()
            if Re_Num.match(line):
                #incrementa la lista di dati dell'evento ricercato
                data=np.append(data,np.array([int(line)]))
            if Re_Record.match(line):
                if (evento%500==0):
                    logging.info(f'Lettura evento {evento}')
                try:
                    Sign, Sigma = Gestioneoff(data)
                    minRapp=np.append(minRapp,np.array([np.abs(np.min(Sign))/Sigma]))
                except:
                    ContErr+=1
                    logging.warning(f'warning. Tot={ContErr}')
                evento+=1
                data=np.array([])
            if Re_EOF.match(line):
                logging.info(f'Lettura evento {evento}')
                Sign, Sigma = Gestioneoff(data)
                try:
                    minRapp=np.append(minRapp,np.array([np.abs(np.min(Sign))/Sigma]))
                except:
                    ContErr+=1
                data=np.array([])
                EOF=True

    print(f'Errori: {ContErr}, eventi: {evento}, ContErr/eventi: {ContErr/evento*100:.2f}%')
    return minRapp

"""-Coinicidence"""
###Attenzione alle wave segnate tra cui fa le coincidenze
def DoubleRecognize(channel1,channel2):
    """-Trigger per coincidenza tripla basata su minimi dei segnali"""
    time=np.arange(0,len(channel1))
    off1, Sigma1=GestioneoffIterativa(channel1)
    off2, Sigma2=GestioneoffIterativa(channel2)

    """-impostazione threshold ricerca minimi"""
    #le scelte della motliplicazione di offset sono state fatte off programma
    Thrs1=150*Sigma1
    Thrs2=50*Sigma2

    """-Ricerca del segnale impulsivo"""
    #minimirelativi
    min_pos1,minval1 = find_peaks(-channel1, height=Thrs1)
    min_pos2,minval2 = find_peaks(-channel2, height=Thrs2)

    Double=False
    M1=M2=0
    for m1 in min_pos1:
        inf,sup = CreaInt(m1,13)
        for m2 in min_pos2:
            if (m2>=inf) & (m2<=sup):
                if not Double:
                    if M1<m1: M1=m1
                    if M2<m2: M2=m2
                    Double = True

    if Double:
        inf,sup = CreaInt(M1,13)
        TS1,C1=TimeStampMeas(channel1,M1,Sigma1,inf,sup)
        TS2,C2=TimeStampMeas(channel2,M2,Sigma2,inf,sup)
        return TS1,C1,TS2,C2
    else:
        return None
def TripleRecognize(channel1, channel2,channel3):
    logging.debug('TripleRecognize')
    """-Trigger per coincidenza tripla basata su minimi dei segnali"""
    time=np.arange(0,len(channel1))
    channel1, Sigma1=Gestioneoff(channel1)
    channel2, Sigma2=Gestioneoff(channel2)
    channel3, Sigma3=Gestioneoff(channel3)

    """-impostazione threshold ricerca minimi"""
    #le scelte della motliplicazione di offset sono state fatte off programma
    Thrs1=15*Sigma1
    Thrs2=15*Sigma2
    Thrs3=15*Sigma3

    """-Ricerca del segnale impulsivo"""
    #minimirelativi
    min_pos1,minval1 = find_peaks(-channel1, height=Thrs1)
    min_pos2,minval2 = find_peaks(-channel2, height=Thrs2)
    min_pos3,minval3 = find_peaks(-channel3, height=Thrs3)

    """-Trigger tripla"""
    Tripla=False
    M1=M2=M3=0
    for m1 in min_pos1:
        inf,sup = CreaInt(m1,13)
        for m2 in min_pos2:
            if (m2>=inf) & (m2<=sup):
                for m3 in min_pos3:
                    if (m3>=inf) & (m3<=sup):
                        Tripla = True
                        if channel1[M1]>channel1[m1]: M1=m1
                        if channel2[M2]>channel2[m2]: M2=m2
                        if channel3[M3]>channel3[m3]: M3=m3

    if Tripla:
        inf,sup = CreaInt(M1,13)
        TS1,C1=TimeStampMeas(channel1,M1,Sigma1,inf,sup)
        TS2,C2=TimeStampMeas(channel2,M2,Sigma2,inf,sup)
        TS3,C3=TimeStampMeas(channel3,M3,Sigma3,inf,sup)
        return TS1,C1,TS2,C2,TS3,C3
    else:
        return None
def QuadRecognize(channel1, channel2,channel3,channel4):
    """-Trigger per coincidenza tripla basata su minimi dei segnali"""
    time=np.arange(0,len(channel1))
    channel1, Sigma1=Gestioneoff(channel1)
    channel2, Sigma2=Gestioneoff(channel2)
    channel3, Sigma3=Gestioneoff(channel3)
    channel4, Sigma4=Gestioneoff(channel4)

    """-impostazione threshold per ricerca minimi, filtraggio minimi"""
    #le scelte della motliplicazione di offset sono state fatte off programma
    Thrs1=15*Sigma1
    Thrs2=15*Sigma2
    Thrs3=15*Sigma3
    Thrs4=15*Sigma4

    """-Ricerca del segnale impulsivo"""
    #minimirelativi
    min_pos1,minval1 = find_peaks(-channel1, height=Thrs1)
    min_pos2,minval2 = find_peaks(-channel2, height=Thrs2)
    min_pos3,minval3 = find_peaks(-channel3, height=Thrs3)
    min_pos4,minval4 = find_peaks(-channel4, height=Thrs4)

    """-Trigger quadrupla"""
    Quadrupla=False
    M1=M2=M3=M4=0
    for m1 in min_pos1:
        inf,sup = CreaInt(m1,13)
        for m2 in min_pos2:
            if (m2>=inf) & (m2<=sup):
                for m3 in min_pos3:
                    if (m3>=inf) & (m3<=sup):
                        for m4 in min_pos4:
                            if (m4>=inf) & (m4<=sup):
                                Quadrupla=True
                                if channel1[M1]>channel1[m1]: M1=m1
                                if channel2[M2]>channel2[m2]: M2=m2
                                if channel3[M3]>channel3[m3]: M3=m3
                                if channel4[M4]>channel4[m4]: M4=m4

    if Quadrupla:
        inf,sup = CreaInt(M1,13)
        TS1,C1=TimeStampMeas(channel1,M1,Sigma1,inf,sup)
        TS2,C2=TimeStampMeas(channel2,M2,Sigma2,inf,sup)
        TS3,C3=TimeStampMeas(channel3,M3,Sigma3,inf,sup)
        TS4,C4=TimeStampMeas(channel4,M4,Sigma4,inf,sup)
        return TS1,C1,TS2,C2,TS3,C3,TS4,C4
    return None

"""Funzioni Read"""
def ReadChannelWave(OutFilePath,dir,wave):

    """-funzione Read wave file"""

    Err = Data = C = Ev = TS = np.array([])
    Tot = Tbuff = 0

    logging.debug('apertura file wave.txt')
    with open(wave) as file:
        CurrLine=0
        StartEv=7
        CurrEv=0
        conterr=0
        print('\n\t\tStart!!\n')
        for line in file:
            if (CurrLine>=StartEv) & (CurrLine<StartEv+RecLen):
                Data = np.append(Data,int(line))
            if re_TimeLine.match(line):
                #lettura triggerstamp
                TCurr = int(re_num.search(line)[0])
                Dtbuff = (TCurr - Tbuff)
                Tbuff = TCurr
                if Dtbuff>0:    #su un manuale è spiegato come misurarlo precisamente
                    Tot+=Dtbuff
            CurrLine+=1
            if CurrLine==StartEv+RecLen:
                try:
                    #signal recognize
                    Data,Sigma=GestioneoffPrimiCamp(Data)
                    min_pos=np.argmin(Data)
                    TSbuff = TimeStampMeas(Data,min_pos,Sigma)
                    DataSign,time = EstraiSegnale(Data,min_pos,Sigma,1,RecLen-1)
                    Cbuff = ChargeMeas(DataSign)
                    TS = np.append(TS,TSbuff)
                    C=np.append(C,int(Cbuff))
                    Ev=np.append(Ev,int(CurrEv))
                except:
                    conterr+=1
                    Err=np.append(Err,(CurrEv))
                    logging.warning(f'problemi con l evento: {CurrEv}')
                Data=np.array([])
                CurrEv+=1
                StartEv=CurrLine+7
                if (CurrEv%2500==0):
                    logging.info(CurrEv)

    logging.debug(f'creazione file {OutFilePath}')
    Res=np.column_stack((Ev,C,TS))
    #Res=np.column_stack((Ev,C))
    np.savetxt(OutFilePath,Res)

    print('\n\t\tFinish!!\n')

    if len(Err) != 0:
        print(f'sono stati trovati {conterr} errori su un massimo di {CurrEv} eventi')
        print(f'rumori/segnali={conterr/CurrEv}')
        #np.savetxt(f'Problemi{wave}',Err)

    print(f'Time {wave}: ', Tot*8, '\tns\n\n')

    return
def Read2ChannelWave(TCOFilePath,dir,app):
    """-funzione read3Channel"""
    logging.debug('apertura file wave.txt')
    #richiede la modifica della directory ed eventale appendice al nome
    w1 = f'{dir}/wave0.txt'
    w2 = f'{dir}/wave1.txt'
    TCOutFile=TCOFilePath

    logging.info(f'creazione file {TCOutFile}')
    with open(TCOutFile,'w') as tcof:
        tcof.write('#\tEvento\t\tTSPM1\t\tcPMT1\t\tTSPM2\t\tcPMT2\n')
    tcof=open(TCOutFile,'a')

    Err=Data1=Data2=np.array([])

    with open(w1) as file1, open(w2) as file2:
        CurrLine=0
        StartEv=7
        CurrEv=0
        conterr=0
        print('Start!!')
        for line1,line2 in itertools.zip_longest(file1,file2):
            if (CurrLine>=StartEv) & (CurrLine<StartEv+RecLen):
                if (line1==None) | (line2==None):
                    break
                Data1=np.append(Data1,np.array([int(line1)]))
                Data2=np.append(Data2,np.array([int(line2)]))
            CurrLine+=1
            if CurrLine==StartEv+RecLen:
                #logging.debug(f'Fine evento, {line1}')
                try:
                    TS1,C1,TS2,C2=DoubleRecognize(Data1, Data2)
                    tcof.write(f'\t{CurrEv}\t{TS1}\t{C1}\t{TS2}\t{C2}\n')
                except:
                    conterr+=1
                    Err=np.append(Err,np.array([CurrEv]))
                    logging.warning(f'problemi con la doppia nell evento: {CurrEv}')
                Data1=Data2=np.array([])
                CurrEv+=1
                StartEv=CurrLine+7
                if (CurrEv%1000==0):
                    logging.info(CurrEv)

    #output
    logging.info('Chiusura file..')
    tcof.close()
    print('Finish!!')
    print(f'sono stati trovati {conterr} errori su un massimo di {CurrEv} eventi')
    print(f'rumori/segnali={conterr/CurrEv}')

    np.savetxt(f'{dir}/ProblemiTripla{app}.txt',Err)

    return TCOutFile
def Read3ChannelWave(TCOFilePath,dir,app):
    """-funzione read3Channel"""
    logging.debug('apertura file wave.txt')
    #richiede la modifica della directory ed eventale appendice al nome
    w1 = f'{dir}/wave0.txt'
    w2 = f'{dir}/wave1.txt'
    w3 = f'{dir}/wave2.txt'
    TCOutFile=TCOFilePath

    logging.info(f'creazione file {TCOutFile}')
    with open(TCOutFile,'w') as tcof:
        tcof.write('#\tEvento\t\tTSPM1\t\tcPMT1\t\tTSPM2\t\tcPMT2\t\tTSPM3\t\tcPMT3\n')
    tcof=open(TCOutFile,'a')

    Err=Data1=Data2=Data3=np.array([])

    with open(w1) as file1, open(w2) as file2, open(w3) as file3:
        CurrLine=0
        StartEv=7
        CurrEv=0
        conterr=0
        print('Start!!')
        for line1,line2,line3, in itertools.zip_longest(file1,file2,file3):
            if (CurrLine>=StartEv) & (CurrLine<StartEv+1030):
                if (line1==None) | (line2==None) | (line3==None):
                    break
                Data1=np.append(Data1,np.array([int(line1)]))
                Data2=np.append(Data2,np.array([int(line2)]))
                Data3=np.append(Data3,np.array([int(line3)]))
            CurrLine+=1
            if CurrLine==StartEv+1030:
                #logging.debug(f'Fine evento, {line1}')
                try:
                    TS1,C1,TS2,C2,TS3,C3=TripleRecognize(Data1, Data2, Data3)
                    tcof.write(f'\t{CurrEv}\t{TS1}\t{C1}\t{TS2}\t{C2}\t{TS3}\t{C3}\n')
                except:
                    conterr+=1
                    Err=np.append(Err,np.array([CurrEv]))
                    #logging.warning(f'problemi con la tripla nell evento: {CurrEv}')
                Data1=Data2=Data3=np.array([])
                CurrEv+=1
                StartEv=CurrLine+7
                if (CurrEv%1000==0):
                    logging.info(CurrEv)

    #output
    logging.info('Chiusura file..')
    tcof.close()
    print('Finish!!')
    print(f'sono stati trovati {conterr} errori su un massimo di {CurrEv} eventi')
    print(f'rumori/segnali={conterr/CurrEv}')

    np.savetxt(f'{dir}/ProblemiTripla{app}.txt',Err)

    return TCOutFile
def Read4ChannelWave(TCOFilePath,dir,app):
    """-Legge gli eventi letti in 4 PMT e restituisce i Time stamp e charge
        degli eventi nel quale èstata riconosciuta una coincidenza tra i 4 segnali"""

    #richide modifiche
    logging.debug('apertura file wave.txt')
    w1 = f'{dir}/wave0.txt'
    w2 = f'{dir}/wave1.txt'
    w3 = f'{dir}/wave2.txt'
    w4 = f'{dir}/wave3.txt'
    TCOutFile = TCOFilePath


    logging.info(f'creazione file {TCOutFile}')
    with open(TCOutFile,'w') as tcof:
        tcof.write('#\tEvento\t\tTSPM1\t\tcPMT1\t\tTSPM2\t\tcPMT2\t\tTSPM3\t\tcPMT3\t\tTSPMT4\t\tcPMT4\n')
    tcof=open(TCOutFile,'a')

    Err=Data1=Data2=Data3=Data4=np.array([])
    with open(w1) as file1, open(w2) as file2, open(w3) as file3, open(w4) as file4:
        CurrLine=0
        StartEv=7
        CurrEv=0
        conterr=0
        print('Start!!')
        for line1,line2,line3,line4 in itertools.zip_longest(file1,file2,file3,file4):
            if (CurrLine>=StartEv) & (CurrLine<StartEv+1030):
                if (line1==None) | (line2==None) | (line3==None) | (line4==None):
                    break
                #incrementa la lista di dati dell'evento ricercato
                Data1=np.append(Data1,np.array([int(line1)]))
                Data2=np.append(Data2,np.array([int(line2)]))
                Data3=np.append(Data3,np.array([int(line3)]))
                Data4=np.append(Data4,np.array([int(line4)]))
            CurrLine+=1
            if CurrLine==StartEv+1030:
                StartEv=CurrLine+7
                try:
                    TS1,C1,TS2,C2,TS3,C3,TS4,C4=QuadRecognize(Data1, Data2, Data3, Data4)
                    tcof.write(f'\t{CurrEv}\t{TS1}\t{C1}\t{TS2}\t{C2}\t{TS3}\t{C3}\t{TS4}\t{C4}\n')
                except:
                    conterr+=1
                    Err=np.append(Err,np.array([CurrEv]))
                    #logging.warning(f'problemi con la quadrupla nell evento: {CurrEv}')
                if (CurrEv%1000==0):
                    logging.info(CurrEv)
                Data1=Data2=Data3=Data4=np.array([])
                CurrEv+=1

    #output
    print('Finish!!')
    logging.info('Chiusura file..')
    tcof.close()
    print(f'sono stati trovati {conterr} errori su un massimo di {CurrEv} eventi')
    print(f'rumori/segnali={conterr/CurrEv}')


    np.savetxt(f'{dir}/ProblemiQuadrupla{app}.txt',Err)

    return TCOutFile

"""-Dt measure function"""
###Le differenze sono fatte tra unità temporali dell'adc
def DtMeasure(TCFile,ResultFile, NumCh):
    """-Lettura file result"""
    logging.info(f'Lettura {TCFile}')
    if NumCh==2:
        ev,t0,c0,t1,c1 = np.loadtxt(TCFile,unpack=True)
        Dt10=t1-t0
        Intestazione='#\tEvNum\tDtCh(1-0)\n'
    if NumCh==3:
        ev,t0,c0,t1,c1,t2,c2 = np.loadtxt(TCFile,unpack=True)
        Dt10=t1-t0
        Dt20,Dt21=t2-t0,t2-t1
        Intestazione='#\tEvNum\t\tDt(1-0)\t\tDt(2-0)\t\tDt(2-1)\n'
    if NumCh==4:
        ev,t0,c0,t1,c1,t2,c2,t3,c3 = np.loadtxt(TCFile,unpack=True)
        Dt10 = t1-t0
        Dt20 , Dt21 = t2-t0 , t2-t1
        Dt30 , Dt31 , Dt32 = t3-t0 , t3-t1 , t3-t2
        Intestazione='#\tEvNum\t\tDt(1-0)\t\tDt(2-0)\t\tDt(2-1)\t\tDt(3-0)\t\tDt(3-1)\t\tDt(3-2)\n'

    """-scrittura su file"""
    if NumCh==2:
        logging.info(f'Scrittura risultati su {ResultFile}')
        with open(ResultFile,'w')as dati:
            dati.write(Intestazione)
            for i in range(len(t1)):
                dati.write(f'\t{ev[i]}\t{Dt10[i]}\n')
    if NumCh==3:
        logging.info(f'Scrittura risultati su {ResultFile}')
        with open(ResultFile,'w')as dati:
            dati.write(Intestazione)
            for i in range(len(t1)):
                dati.write(f'\t{ev[i]}\t\t{Dt10[i]}\t\t{Dt20[i]}\t\t{Dt21[i]}\n')
    if NumCh==4:
        logging.info(f'Scrittura risultati su {ResultFile}')
        with open(ResultFile,'w')as dati:
            dati.write(Intestazione)
            for i in range(len(t1)):
                dati.write(f'\t{ev[i]}\t\t{Dt10[i]}\t\t{Dt20[i]}\t\t{Dt21[i]}\t\t{Dt30[i]}\t\t{Dt31[i]}\t\t{Dt32[i]}\n')

    print(f'\nI Delta t tra i vari canali per ogni evento rivelato sono nel file {ResultFile}\n')

    return ResultFile
def Dt2ChMeasure(TSFile1, TSFile2, ResultFile):
    """-Lettura file result"""
    logging.debug(f'Lettura {TSFile1} , {TSFile2}\n')

    ev1,c1,t1 = np.loadtxt(TSFile1,unpack=True)
    ev2,c2,t2 = np.loadtxt(TSFile2,unpack=True)

    #conto Dt
    Dt = (t1-t2)*4

    """-scrittura su file"""
    Res = np.column_stack((ev1,Dt))
    np.savetxt(ResultFile,Res)
    logging.info(f'Scrittura risultati su {ResultFile}')

    return

def TEven(FileName):
    """-Misura il tempo impiegato per una presa dati CONTINUA"""

    import re

    ##Per ricrcare nel file
    re_TimeLine = re.compile('Trigger')
    re_num = re.compile(r'\b\d+\b')
    ##Inizializzaizoni
    TCurr = Tbuff = Tot = 0
    yy=[]
    i = 0
    ##Lettura
    with open(FileName) as File:
        for line in File:
            if re_TimeLine.match(line):
                i+=1
                #lettura triggerstamp
                TCurr = int(re_num.search(line)[0])
                Dtbuff = (TCurr - Tbuff)
                Tbuff = TCurr
                #if Dtbuff>0:
                Tot+=Dtbuff
                if i%1000 == 0:
                    logging.info(i)
                yy.append(TCurr)

    print(Tot*8)

    yy=np.array(yy)
    xx=np.arange(len(yy))
    plt.plot(xx,yy,'.')
    plt.show()


    return Tot

def TriImp(FilePathPMT1,FilePathPMT2, dir, orario):

    Ev1,c1,Ts1 = np.loadtxt(FilePathPMT1,unpack=True)
    Ev2,c2,Ts2 = np.loadtxt(FilePathPMT2,unpack=True)

    ##RettaZero
    #a1 = 207.7
    #a2 = 69.772
    ##Retta
    a1 = 197.6
    b1 = 8485
    a2 = 70
    b2 = -156


    Ev = somma = differenza = np.array([])

    for i in range(len(Ev1)):
        c1[i] = (c1[i]-b1)/a1
        c2[i] = (c2[i]-b2)/a2
        somma=np.append(somma,(c1[i]+c2[i]))
        differenza=np.append(differenza,(c1[i]-c2[i]))
        Ev=np.append(Ev,Ev1[i])

    somma=np.column_stack((Ev,somma))
    differenza=np.column_stack((Ev,differenza))

    OutFilePathSum = f'{dir}/Somma_3Imp_Retta_{orario}.txt'
    OutFilePathDiff = f'{dir}/Diff_3Imp_Retta_{orario}.txt'

    np.savetxt( OutFilePathSum , somma )
    np.savetxt( OutFilePathDiff , differenza )

from scipy.optimize import curve_fit

def ConvSig(x,C,s,d,start):
    Con = -C*(np.exp(-s*(x-start))-np.exp(-d*(x-start)))
    f = np.heaviside(x-start,0)*Con
    return f
def findstart4Fit(Sign_arr,min_pos,meno,limit):
    """- find the first data position of the negative slope of impulse"""
    if (min_pos-meno>=limit):
        i=min_pos-meno
    else:
        i=limit
    """
    i=min_pos
    while(Sign_arr[i]<0.001*Sign_arr[min_pos]):
        if i==limit:
            break
        i-=1
    """
    return i
def fitsignal(Signal,min_pos):

    time=np.linspace(0,RecLen-1,RecLen)
    #print(f'minimo sig\t{Signal[min_pos]:.6}')
    #print(f'min pos\t\t{min_pos}')
    #print(f'start\t\t{start}')

    w=np.ones(len(Signal))

    chi2r=10000
    meno=10
    while ((chi2r>5000)&(meno<60)):
        start=findstart4Fit(Signal,min_pos,meno,1)
        initC=[Signal[min_pos],0,0,start]
        parsC,covmC = curve_fit(ConvSig,time,Signal,initC,w,absolute_sigma=False)
        C,s,d,start=parsC[0],parsC[1],parsC[2],parsC[3]
        #print(f'cost\t\t\t{C:.6}')
        #print(f'salita\t\t\t{s:.4}')
        #print(f'discesa\t\t\t{d:.4}')
        #print(f'start fit\t\t{start:.5}')
        chi2 = ((1/w**2*(Signal-ConvSig(time,parsC[0],parsC[1],parsC[2],parsC[3]))**2)).sum()
        ndof=len(Signal)-len(initC)
        #print(chi2, ndof)
        #print(f'chi2:\t\t{chi2/ndof:.7}')
        chi2r=chi2/ndof
        meno+=5

    yy = ConvSig(time,C,s,d,start)
    ##Creare il segnale con cui farci poi l'integrale
    plt.plot(time,yy,label='FitCurve',color='orange')
    #plt.legend()
    #plt.show()

    return yy

###     PROCESS     ###
def Positroni():
    Data = np.array([])

    wave='1a_misura_me_20210305/PMT2_cesio_1729.txt'
    #Eventoricerca=0
    TEven(wave)

    """
    with open(wave) as file:
        CurrLine=CurrEv=0
        StartEv=7
        for line in file:
            if (CurrLine>=StartEv) & (CurrLine<StartEv+RecLen):
                Data = np.append(Data,int(line))
            CurrLine+=1
            if (CurrLine==StartEv+RecLen):
                if (CurrEv==Eventoricerca):
                    chi2r=10000
                    meno=10
                    #print(CurrEv)
                    #Try:
                    #signal recognize
                    off,Sigma=GestioneoffIterativa(Data)
                    Data-=off
                    plt.plot(np.linspace(0,len(Data),len(Data)),Data,label='Data')
                    min_pos=np.argmin(Data)
                    yy=fitsignal(Data,min_pos,meno)
                    #except:
                    #logging.warning(f'Error {CurrEv}')

                CurrEv+=1
                Data=np.array([])
                StartEv=CurrLine+7

                if (CurrEv==Eventoricerca+1):
                    break
                    print(CurrEv)
    """

if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-help', help='print this message?')

    parser.add_argument('-info',action='store_true', help='print info messages')
    parser.add_argument('-debug',action='store_true', help='print debug messages')


    args = parser.parse_args()

    if args.info:
        logging.basicConfig(level=logging.INFO)
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    Positroni()
