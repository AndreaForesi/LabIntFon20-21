import argparse
import numpy as np
import logging
from itertools import zip_longest
from glob import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import re


re_TimeLine = re.compile('Trigger')
re_num = re.compile(r'\b\d+\b')


###     Library     ###
"""Gestione offset"""
def GestioneoffIterativa(channel):
    sigma=11
    l=0
    r=50
    while ( (sigma>10) & (r<RecLen-50) ):
        a=channel[l:r]
        offset=np.mean(a)
        sigma=np.std(a)
        l=l+50
        r=r+50
    return offset,sigma
def GestioneoffPrimiCamp(channel):
    """-Si esegue media e devstd dell'offest considerando i primi punti"""
    Dati    = channel[0:201]
    off     = np.mean(Dati)
    Sigma   = np.std(Dati)
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

    """-signal"""
def EstraiSegnale(channel,min_pos):
    """-Dai dati channel ricava una porzione precisa di segnale"""
    #time array
    start_time=findstart4Fit(channel,min_pos,30)
    stop_time=findstop(channel,min_pos)
    Signal = np.array(channel[start_time:(stop_time)])
    Signal_time = np.arange(start_time,stop_time)
    return Signal,Signal_time
def findstart(Sign_arr,min_pos,sigma):
    """- find the first data position of the negative slope of impulse"""
    """
    if (min_pos-30>=1):
        i=min_pos-30
    else:
        i=1
    """
    i=min_pos
    while(Sign_arr[i]<-sigma):
        if i==1:
            break
        i=i-1
    return i
def findstop(Sign_arr,min_pos):
    """- find the last data position of the negative slope of impulse"""
    if (min_pos+350<=RecLen-1):
        i=min_pos+350
    else:
        i=RecLen
    """
    i=min_pos
    while(Sign_arr[i]<0.2*Sign_arr[min_pos]):
        if i==(limit-1):
            break
        i=i+1
    """
    return i

"""-Fit"""
def retta(x,m,q):
    return m*x + q
def ConvSig(x,C,s,d,start):
    """-Curve of the fit"""
    Con = -C*(np.exp(-s*(x-start))-np.exp(-d*(x-start)))
    f = np.heaviside(x-start,0)*Con
    return f
def findstart4Fit(Sign_arr,min_pos,meno):
    """- find the Start for the fit"""
    if (min_pos-meno>=1):
        i=min_pos-meno
    else:
        i=1
    return i
def fitsignal(Signal,min_pos):
    """-Fit of the data event"""
    time=np.linspace(0,RecLen-1,RecLen)

    meno=10
    chi2r=10000
    w=np.ones(len(Signal))
    while ((chi2r>5000)&(meno<60)):
        start=findstart4Fit(Signal,min_pos,meno)
        initC=[Signal[min_pos],0,0,start]
        #print(meno)
        parsC,covmC = curve_fit(ConvSig,time,Signal,initC,w,absolute_sigma=False)
        C,a,b,startfit=parsC[0],parsC[1],parsC[2],parsC[3]
        chi2 = ((1/w**2*(Signal-ConvSig(time,parsC[0],parsC[1],parsC[2],parsC[3]))**2)).sum()
        ndof=len(Signal)-len(initC)
        chi2r=chi2/ndof
        #print(chi2r)
        meno+=5


    yy = ConvSig(time,C,a,b,startfit)
    #plt.plot(time,yy);plt.show()

    if ((chi2r>5000)|(meno==60)):
        logging.warning('No fit')
        return


    #Charge = -C*(np.exp(b*start)/b*(np.exp(-b*yy[int(start)+1])-np.exp(-b*yy[int(start)+380]))-np.exp(a*start)/a*(np.exp(-a*yy[int(start)+1])-np.exp(-a*yy[int(start)+380])))

    time=np.linspace(0,RecLen-1,2081)
    Charge=np.trapz(yy[int(startfit):int(startfit)+1442])
    time_stamp=startfit

    return Charge , time_stamp

"""-Measurement"""
def ChargeMeas(Signal):
    """-Charge measurement"""
    Charge=np.abs(np.trapz(Signal))
    #Charge=np.abs(np.min(Signal))
    return Charge
def TimeStampMeas(Signal,min_pos,sigma):
    """-Time stamp measurement"""

    """
    ###Media negative slope
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

    ###FIT RETTA
    init=[0,0]
    y=Signal[start:stop]
    time=np.arange(start,stop,1)
    w=np.ones(len(y))
    pars,covm = curve_fit(retta,time,y,init,w,absolute_sigma=False)
    m,q=pars[0],pars[1]
    chi2 = ((1/w**2*(y-retta(time,pars[0],pars[1]))**2)).sum()
    ndof=len(y)-len(init)
    chi2r=chi2/ndof
    """
    """
    start=findstart(Signal,min_pos,5*sigma)
    stop=start+1
    x1,y1=start,Signal[start]
    x2,y2=start+1,Signal[start+1]
    m=(y1-y2)/(x1-x2)
    q=y1-m*x1

    time_stamp = (-5*sigma-q)/m

    #print(time_stamp)

    #xx=np.linspace(start-20,stop+20,10000)
    #plt.plot(xx,retta(xx,m,q),label='retta')
    #plt.plot(time_stamp,retta(time_stamp,m,q),'.')
    #plt.legend()
    #plt.show()
    """

    return time_stamp
def Timedata(WaveFileName):
    """-Misura il tempo impiegato per una presa dati CONTINUA"""
    logging.info(f'\tInizio file \n\t{WaveFileName} \n ')
    ##Inizializzaizoni
    TCurr = Tbuff = Tot = 0
    ##Lettura
    with open(WaveFileName) as File:
        for line in File:
            if re_TimeLine.match(line):
                #lettura triggerstamp
                TCurr = int(re_num.search(line)[0])
                Dtbuff = (TCurr - Tbuff)
                Tbuff = TCurr
                if Dtbuff>0:
                    Tot+=Dtbuff

    print(f'Time {WaveFileName}: ', Tot*8, '\tns\n\n')

    return
def Dt2ChMeasure(TSFile1, TSFile2, ResultFile):
    """-Lettura file result"""
    logging.debug(f'Lettura {TSFile1} , {TSFile2}\n')

    t1 = np.loadtxt(TSFile1,unpack=True)
    t2 = np.loadtxt(TSFile2,unpack=True)

    Dt=np.array([])
    EV=np.array([])

    #conto Dt
    #for i in range(len(ev1)):
    #    if ((c1[i]!=0)&(c2[i]!=0)):
    #        Dt = np.append(Dt, (t1[i]-t2[i])*4)
    #        EV = np.append(EV, ev1[i])

    Dt=t1-t2

    """-scrittura su file"""
    #Res = np.column_stack((EV,Dt))
    np.savetxt(ResultFile,Dt)
    #logging.info(f'Scrittura Ev,Dt(1-2)[ns],C1,C2 su {ResultFile}')

    return
def Dt3ChMeasure(TSFile1, TSFile2, TSFile3, ResultFile):
    """-Lettura file result"""
    logging.debug(f'Lettura {TSFile1} , {TSFile2}, {TSFile3} \n')

    ev1,c1,t1 = np.loadtxt(TSFile1,unpack=True)
    ev2,c2,t2 = np.loadtxt(TSFile2,unpack=True)
    ev3,c3,t3 = np.loadtxt(TSFile3,unpack=True)

    #conto Dt
    Dt12 = (t1-t2)*4
    Dt13 = (t1-t3)*4
    Dt23 = (t2-t3)*4

    """-scrittura su file"""
    Res = np.column_stack((ev1,Dt12,Dt13,Dt23,c1,c2,c3))
    np.savetxt(ResultFile,Res)
    logging.info(f'Scrittura Ev,Dt(1-2)[ns],Dt(1-3)[ns],Dt(2-3)[ns],C1,C2,C3 su\n\t\t{ResultFile}')

    return


def Correzione(FilePath,lenFile):

    Ev,C,Ts=np.loadtxt(FilePath,unpack=True)
    #7839
    NewEv=np.zeros(lenFile)
    NewC=np.zeros(lenFile)
    NewTs=np.zeros(lenFile)

    for i in range(len(Ev)):

        NewC[int(Ev[i])]=C[i]

        print(NewC[int(Ev[i])])

        NewTs[int(Ev[i])]=Ts[i]

        NewEv[int(Ev[i])]=Ev[i]

    print(NewC)

    Res=np.column_stack((NewEv,NewC,NewTs))

    np.savetxt(f'New{FilePath}',Res)

"""Histogram"""
def HistCheck(wave):
    """-funzione Read wave file"""

    Data = Hist = Ev = np.array([])

    logging.debug('apertura file wave.txt')
    with open(wave) as file:
        CurrLine=0
        StartEv=7
        CurrEv=0
        print('\n\t\tStart!!\n')
        for line in file:
            if (CurrLine>=StartEv) & (CurrLine<StartEv+1030):
                Data = np.append(Data,int(line))
            CurrLine+=1
            if CurrLine==StartEv+RecLen:
                off,Sigma=GestioneoffIterativa(Data)
                #Data-=off

                Hist=np.append(Hist,Sigma)
                #Ev=np.append(Ev,int(CurrEv))

                Data=np.array([])
                CurrEv+=1
                StartEv=CurrLine+7

                if (CurrEv%1000==0):
                    logging.info(CurrEv)
                #if (CurrEv==1): break

    OutFilePath=f'Hist{wave}.txt'
    logging.debug(f'creazione file {OutFilePath}')
    #Res=np.column_stack((Ev,Hist))
    np.savetxt(OutFilePath,Hist)

    print('\n\t\tFinish!!\n')

    return
def SogliaCheck(DataHist):
    Ev,Data,Ts=np.loadtxt(DataHist,unpack=True)
    PartEvent=np.array([])
    Soglia=200000
    for i in range(len(Data)):
        if Data[i]>Soglia:
            PartEvent=np.append(PartEvent,Ev[i])

    np.savetxt('PartEvent.txt',PartEvent)

"""Funzioni Read"""
def ReadChannelWave(OutFileAppend,wave):

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
            CurrLine+=1

            if CurrLine==StartEv+RecLen:
                try:
                    ###signal recognize
                    off,Sigma=GestioneoffIterativa(Data)
                    Data-=off
                    min_pos=np.argmin(Data)
                    TSbuff = TimeStampMeas(Data,min_pos,Sigma)
                    #DataSign,time = EstraiSegnale(Data,min_pos)
                    #Cbuff=ChargeMeas(DataSign)
                    TS = np.append(TS,TSbuff)
                    #C=np.append(C,int(Cbuff))
                    Ev=np.append(Ev,int(CurrEv))
                except:
                    conterr+=1
                    Err=np.append(Err,(CurrEv))
                    logging.warning(f'problemi con l evento: {CurrEv}')

                Data=np.array([])
                CurrEv+=1
                StartEv=CurrLine+7

                if (CurrEv%1000==0):
                    logging.info(CurrEv)

    OutFilePath=f'{OutFileAppend}{wave}'
    #Res=np.column_stack((Ev,C,TS))
    #Res=np.column_stack((Ev,C))
    Res=np.column_stack((Ev,TS))
    logging.debug(f'creazione file {OutFilePath}')
    np.savetxt(OutFilePath,Res)

    print('\n\t\tFinish!!\n')

    if len(Err) != 0:
        print(f'sono stati trovati {conterr} errori su un massimo di {CurrEv} eventi')
        print(f'rumori/segnali={conterr/CurrEv}')
        np.savetxt(f'Problemi{wave}',Err)

    return


###     PROCESS     ###
RecLen=1030
def Positroni():
    """-Positron Annihilation experience"""

    File=glob('PMT*_doppie_1725.txt')
    for wave in File:

        OutFileAppend='Out_Ts5Sigma_'
        ReadChannelWave(OutFileAppend,wave)

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
