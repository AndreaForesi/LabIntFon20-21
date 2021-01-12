import argparse
import itertools
import numpy as np
import logging
from scipy.signal import argrelmin,find_peaks

###     Library     ###
"""-signal"""
def Gestioneoff(channel):
    """-Gestione offset"""
    #impostare una soglia per gestire l'offset
    SOGLIASEGNALEDown=14650
    SOGLIASEGNALEUp=14775
    RumoreChannel=np.ma.masked_where((channel>SOGLIASEGNALEUp)|(channel<SOGLIASEGNALEDown),channel)    #ricerco il rumore
    offset=RumoreChannel.mean()    #calcolo
    Doffset=RumoreChannel.std()  #verrà usato per definire quando potremo considerare il segnale
    channel=channel-offset    #sposto il segnale

    return channel, Doffset
def TimeStampMeas(channel,min_pos,Sigma,inf,sup):

    """-Time stamp and charge measurement"""
    ####Sia Charge che Time-Stamp sono calcolati con le scelte di start e stop time
    ##time array
    time=np.arange(len(channel))

    ###inizio e fine del segnale da considerare
    ##Seleziona la fine e l'inizio del segnale trmite findstart e findstop (minore di 1 sigma)
    start_time=findstart(channel,min_pos,Sigma,inf)
    stop_time=findstop(channel,min_pos,Sigma,sup)
    #start_time=
    #stop_time=

    ##Estrazione del segnale
    """-NEGATIVE SLOPE WITH 10-90% selection                    #commenta i metodi da non usare
    negative_slope_time=np.ma.masked_where((time>(min_pos))|(time<start_time), time)   #|(channel<min_pos*0.1)|(channel>min_pos*0.9)
    negative_slope_value=np.ma.masked_where((time>(min_pos))|(time<start_time), channel)     #|(channel<min_pos*0.1)|(channel>min_pos*0.9)
    """

    #Signal Selection
    Signal,Signal_time = estrai_segnale(channel,time,(start_time-1),(min_pos+3))


    ###CALCOLO TIME-STAMP
    """-Media di qualcosa"""
    #time_stamp=negative_slope_time.mean()
    """-Baricentro del segnale"""
    time_stamp=np.trapz(Signal_time*Signal)/np.trapz(Signal)

    """-Charge"""
    Charge=np.trapz(Signal)

    """-Ritorna Time-Stamp e Charge dell'impulso"""
    return time_stamp, Charge
def findstart(Sign_arr,min_pos,sigma,limit):
    """- find the first data position of the negative slope of impulse"""
    i=min_pos
    while(Sign_arr[i]<(-sigma)):
        if i==limit:
            break
        i=i-1
    return i
def findstop(Sign_arr,min_pos,sigma,limit):
    """- find the last data position of the negative slope of impulse"""
    i=min_pos
    while(Sign_arr[i]<(-sigma)):
        if i==(limit-1):
            break
        i=i+1
    return i
def estrai_segnale(channel,time,start_time,stop_time):
    Signal = channel[start_time:(stop_time+1)]
    Signal_time = time[start_time:(stop_time+1)]
    return Signal,Signal_time


"""-Threshold for discrimination of signal"""
def ThrsDet(PMTDataFile):
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
def DoubleRecognize(channel1,channel2):
    """-Trigger per coincidenza tripla basata su minimi dei segnali"""
    time=np.arange(0,len(channel1))
    channel1, Sigma1=Gestioneoff(channel1)
    channel2, Sigma2=Gestioneoff(channel2)

    """-impostazione threshold ricerca minimi"""
    #le scelte della motliplicazione di offset sono state fatte off programma
    Thrs1=15*Sigma1
    Thrs2=15*Sigma2

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
            if (CurrLine>=StartEv) & (CurrLine<StartEv+1030):
                if (line1==None) | (line2==None):
                    break
                Data1=np.append(Data1,np.array([int(line1)]))
                Data2=np.append(Data2,np.array([int(line2)]))
            CurrLine+=1
            if CurrLine==StartEv+1030:
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

"""-Altro"""
def CreaInt(Pos,Dt):
    Sup = Pos+Dt
    if Sup>1030:
        Sup=1030
    inf = Pos-Dt
    if inf<0:
        inf=0
    return inf,Sup


"""-Measurement"""
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
"""-Velocity calculate"""
###Inserire i ritardi giusti e alcune misure
def VMuMeasure(DTMeasureFile,VResultFilePath,dir,app):
    """-calcolo velocità particelle"""
    A=np.loadtxt(DTMeasureFile)
    dt21=A[:,1]*4
    dt31=A[:,2]*4
    dt32=A[:,3]*4

    """-Ritardi nel sistema di acquisizione"""
    R12=-2.74
    R13=10.0
    R23=12.8


    ##conti per la velocità
    L=280
    vl=15.28        #calcolata usando l'oscilloscopio
    #vl=            #calcolata con l'adc
    PosX=-(dt21+R12)*vl/2    #cateto
    h=175.25
    d=np.sqrt((PosX**2)+(h**2))
    #TOF
    TOF1=(dt31+PosX/vl+R13)
    TOF2=(dt32-PosX/vl+R23)
    TOF=(TOF1+TOF2)/2
    TOF=(dt31+dt32+R13+R23)/2
    #v
    v1=d/TOF1
    v2=d/TOF2
    v=d/TOF

    cosT=PosX/d

    ##scrivo su txt per farci l' istogramma con root
    np.savetxt(f'{dir}/TOF{app}.txt',TOF)
    np.savetxt(f'{dir}/PosX{app}.txt',PosX)
    np.savetxt(f'{dir}/cosT{app}.txt',cosT)
    np.savetxt(f'{dir}/cosT2{app}.txt',cosT**2)
    #np.savetxt(f'{dir}/V1{app}.txt',v1)
    #np.savetxt(f'{dir}/V2{app}.txt',v2)
    #np.savetxt(f'{dir}/V{app}.txt',v)
    """
    """
    ev=A[:,0]
    ev=np.array(ev,dtype=int)
    Resu=np.column_stack((ev,v))
    #np.savetxt(VResultFilePath,Resu)


###     PROCESS     ###
def TimeOfFlight(TC,Dt,V):
    """-Time of flight experience"""
    ###     COSE DA CAMBIARE    ###
    ##Inserire la coincidenza tra 3 o 4 canali
    ch=3               #indicare quanti caali si vogliono mettere in coincidenza
    ##Inserire la cartella dal quale prendere i dati
    dir='2020_12_15_doppie'              #str
    #opzionale, unappedice col quale riconoscere i dati
    app='TriplaFore'              #str

    ###Resto delprocesso###
    """-File Stamp"""
    TCFile=f'{dir}/TsC{app}.txt'
    if TC == True:
        if ch==2:
            #exitfile=dir/TCResult3PMT(app).txt
            TCFile=Read2ChannelWave(TCFile,dir,app)
        if ch==3:
            #exitfile=dir/TCResult3PMT(app).txt
            TCFile=Read3ChannelWave(TCFile,dir,app)
        if ch==4:
            #exitfile=dir/TCResult4PMT(app).txt
            TCFile=Read4ChannelWave(TCFile,dir,app)
    """-Dt File"""
    DtFile=f'{dir}/Dt{app}.txt'
    if Dt == True:
        DtFile=DtMeasure(TCFile,DtFile,ch)

    """-V measure"""
    if V == True:
        VResultFilePath=f'{dir}/EvV{app}.txt'
        VMuMeasure(DtFile,VResultFilePath,dir,app)


if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-help', help='print this message?')

    parser.add_argument('-info',action='store_true', help='print info messages')
    parser.add_argument('-debug',action='store_true', help='print debug messages')

    parser.add_argument('-TC',default='False',action='store_true', help='Abilitate Time Charge measure of signal')
    parser.add_argument('-Dt',default='False',action='store_true', help='Abilitate the measure of Delta time about every channel')
    parser.add_argument('-V',default='False',action='store_true', help='Abilitate the Velocity measure')

    args = parser.parse_args()

    if args.info:
        logging.basicConfig(level=logging.INFO)
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    TimeOfFlight(args.TC,args.Dt,args.V)
