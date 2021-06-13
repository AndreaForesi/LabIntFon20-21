import numpy as np
import argparse
import logging
import matplotlib.pyplot as plt

"""-Efficienza PMT e telescopi"""
tau1=tau2=tau3=50e-9
Doppie=1000
def EffPMT(dir,ConteggiFilePath,Plot,Tabella,PMT):
    """-Ricerca punti lavoro PMT"""


    logging.debug(f'Apertura file {ConteggiFilePath}...')
    VAlim,time,C1,C2,C3,Triple=np.loadtxt(f'{dir}/{ConteggiFilePath}',unpack=True,skiprows=1)


    logging.debug('Misure...')
    ###Frequenze singole
    f1=C1/time;f2=C2/time;f3=C3/time
    Df1=np.sqrt(C1)/time;Df2=np.sqrt(C2)/time;Df3=np.sqrt(C3)/time
    ###F doppie
    fDo=Doppie/time
    DfDo=np.sqrt(Doppie)/time
    ###Frequenze casuali triple
    fcas=2*tau1*fDo*f3
    Dfcas=(3/100 + Df3/f3)*fcas
    Pcas=fcas*tau3*100
    ###efficenza
    eff=(Triple-fcas*time)/Doppie
    Deff=np.sqrt((eff*(1-eff))/Doppie)
    ###F triple
    fTr=Triple/time
    DfTr=Deff*Doppie/time/100
    #DPcas=Dfcas*100*tau2
    eff*=100
    Deff*=100


    logging.debug('Result Stack...')
    Res=np.column_stack((VAlim,f1,f2,f3,fDo,fTr,eff,fcas))
    ResultFilePath=f'{dir}/Eff_{ConteggiFilePath}'
    np.savetxt(ResultFilePath,Res)

    if (PMT==1)|(PMT==4):
        fx,Dfx=f1,Df1
    if (PMT==2)|(PMT==5):
        fx,Dfx=f2,Df2
    if (PMT==3)|(PMT==6):
        fx,Dfx=f3,Df3


    if Tabella:
        TabFilePath = f'{dir}/Tab_{ConteggiFilePath}'
        file = open(TabFilePath,'w')
        for i in range(len(VAlim)):
            file.writelines(f'\t{VAlim[i]:.1f} & {eff[i]:.3f}\,$\pm$\,{Deff[i]:.3f} & {fDo[i]:.3f}\,$\pm$\,{DfDo[i]:.3f} & {fTr[i]:.3f}\,$\pm$\,{DfTr[i]:.3f}  & {fx[i]:.3f}\,$\pm$\,{Dfx[i]:.3f}   \\\ \n\t\\hline\n')
        file.close()


    if Plot:
        logging.debug('Plotting data..')

        fig, (graph1,graph2)=plt.subplots(2)
        fig.suptitle(f'PMT-{PMT} Vs V Alimentazione')
        graph1.set(xlabel='', ylabel='Efficienza [%]')
        graph2.set(xlabel='Alimentazione [V]', ylabel='frequenza [Hz]')
        graph1.errorbar(VAlim,eff,Deff,linestyle='',color='b',marker='.',label='Soglia = -35mV')
        graph2.errorbar(VAlim,fx,Dfx,fmt='.',label='Soglia = -35mV')
        graph1.grid();graph1.legend();graph2.grid();graph2.legend()

        plt.show()
def MultiPlotEff(dir,FilePath1,FilePath2,PMT):
    """-Ricerca punti lavoro PMT"""


    logging.info(f'Apertura file {FilePath1}...')
    VAlim1,time,C1,C2,C3,Triple=np.loadtxt(f'{dir}/{FilePath1}',unpack=True,skiprows=1)

    logging.debug('Misure 1...')
    ###Frequenze singole
    if (PMT==1)|(PMT==4):
        f1=C1/time
        Df1=np.sqrt(C1)/time
    if (PMT==2)|(PMT==5):
        f1=C2/time
        Df1=np.sqrt(C2)/time
    if (PMT==3)|(PMT==6):
        f1=C3/time
        Df1=np.sqrt(C3)/time

    ###efficenza
    eff1=Triple/Doppie
    Deff1=np.sqrt((eff1*(1-eff1))/Doppie)
    eff1*=100
    Deff1*=100

    logging.info(f'Apertura file {FilePath2}...')
    VAlim2,time,C1,C2,C3,Triple=np.loadtxt(f'{dir}/{FilePath2}',unpack=True,skiprows=1)

    logging.debug('Misure 2...')
    ###Frequenze singole
    if (PMT==1)|(PMT==4):
        f2=C1/time
        Df2=np.sqrt(C1)/time
    if (PMT==2)|(PMT==5):
        f2=C2/time
        Df2=np.sqrt(C2)/time
    if (PMT==3)|(PMT==6):
        f2=C3/time
        Df2=np.sqrt(C3)/time

    ###efficenza
    eff2=Triple/Doppie
    Deff2=np.sqrt((eff2*(1-eff2))/Doppie)
    eff2*=100
    Deff2*=100


    logging.debug('Plotting data..')
    fig, (graph1,graph2)=plt.subplots(2)
    if (PMT==1)|(PMT==4):
        graph1.errorbar(VAlim1,eff1,Deff1,linestyle='',color='b',marker='.',label='Soglia = -25mV')
        graph1.errorbar(VAlim2,eff2,Deff2,linestyle='',color='r',marker='.',label='Soglia = -30mV')
    if (PMT==2)|(PMT==5):
        graph1.errorbar(VAlim1,eff1,Deff1,linestyle='',color='b',marker='.',label='Soglia = -35mV')
        graph1.errorbar(VAlim2,eff2,Deff2,linestyle='',color='r',marker='.',label='Soglia = -40mV')
    if (PMT==3)|(PMT==6):
        graph1.errorbar(VAlim1,eff1,Deff1,linestyle='',color='b',marker='.',label='Soglia = -35mV')
        graph1.errorbar(VAlim2,eff2,Deff2,linestyle='',color='r',marker='.',label='Soglia = -40mV')

    fig.suptitle(f'PMT-{PMT} Vs V Alimentazione')
    graph1.set(xlabel='', ylabel='Efficienza [%]')
    graph2.set(xlabel='Alimentazione [V]', ylabel='frequenza [Hz]')
    if (PMT==1)|(PMT==4):
        graph2.errorbar(VAlim1,f1,Df1,color='b',fmt='.',label='Soglia = -25mV')
        graph2.errorbar(VAlim2,f2,Df2,color='r',fmt='.',label='Soglia = -30mV')
    if (PMT==2)|(PMT==5):
        graph2.errorbar(VAlim1,f1,Df1,color='b',fmt='.',label='Soglia = -35mV')
        graph2.errorbar(VAlim2,f2,Df2,color='r',fmt='.',label='Soglia = -40mV')
    if (PMT==3)|(PMT==6):
        graph2.errorbar(VAlim1,f1,Df1,color='b',fmt='.',label='Soglia = -35mV')
        graph2.errorbar(VAlim2,f2,Df2,color='r',fmt='.',label='Soglia = -70mV')

    graph1.legend();graph2.legend();graph2.semilogy()
    graph1.grid();graph2.grid()

    plt.show()

"""Trigger eventi nei telescopi, ?SCIAMI?"""
def CorTrigD0N(Time,Plot):
    """-Calcolo per le Correzioni dovute al ripple"""
    Dt=[]
    Max85=Max=0
    t1=t2=0
    max=min=0
    for i in range(len(Time)):
        t2=Time[i]

        if t2<t1:
            if (t1-t2<100)&(t1-t2>Max85):
                Max85=t1-t2
            if (t1-t2>200) & (t1>Max):
                Max=t1

            min=t2
        else:
            Dt.append(t2-t1)

        t1=Time[i]

    Dt=np.array(Dt)
    print(f'Max: {Max}')
    print(f'MaxDt stretto: {Max85} .')
    print(f'Delta eventi t: {np.mean(Dt)}+-{np.std(Dt)/np.sqrt(len(Dt))} .')

    if Plot:
        plt.plot(np.arange(len(Time)),Time)
        plt.xlabel('Numero evento letto');plt.ylabel('Tempo [s]')
        plt.show()

    return Max,Max85,np.mean(Dt)
def telescopio(Time,Plot):
    """-Sistema i tempi di rivelazione"""
    Max,Max85,DtEventi=CorTrigD0N(Time,True)

    max=min=0
    t1=t2=0
    max=min=0
    for i in range(len(Time)):
        t2=Time[i]
        max=t2
        if (t2<t1)&(t1-t2<200):
            Time[i:]+=Max85+DtEventi
        if (t2<t1)&(t1-t2>200):
            Time[i:]+=Max+DtEventi

        t1=Time[i]

    if Plot:
        plt.plot(np.arange(len(Time)),Time)
        plt.xlabel('Numero evento letto');plt.ylabel('Tempo [s]')
        plt.show()

    return Time
def Separare(dir,FilePath):
    """-Separa gli eventi letti dal De0Nano"""

    NumTel,Time=np.loadtxt(f'{dir}/{FilePath}',unpack=True)

    TimeTriT1=TimeTriT2=[]          #Canali   1-5
    TimeD12=TimeD13=TimeD23=[]      #Canali   6-7-8
    TimeD45=TimeD46=TimeD56=[]      #Canali   2-3-4

    for i in range(len(NumTel)):
        if (NumTel[i]==1):
            TimeTriT1.append(Time[i])
        if (NumTel[i]==2):
            TimeD45.append(Time[i])
        if (NumTel[i]==3):
            TimeD46.append(Time[i])
        if (NumTel[i]==4):
            TimeD56.append(Time[i])
        if (NumTel[i]==5):
            TimeTriT2.append(Time[i])
        if (NumTel[i]==6):
            TimeD12.append(Time[i])
        if (NumTel[i]==7):
            TimeD13.append(Time[i])
        if (NumTel[i]==8):
            TimeD23.append(Time[i])
    print(len(TimeTriT1))
    TimeTriT1=telescopio(np.array(TimeTriT1),True)
    TimeD12=telescopio(np.array(TimeD12),True)
    TimeD13=telescopio(np.array(TimeD13),True)
    TimeD23=telescopio(np.array(TimeD23),True)

    TimeTriT2=telescopio(np.array(TimeTriT2),True)
    TimeD45=telescopio(np.array(TimeD45),True)
    TimeD46=telescopio(np.array(TimeD46),True)
    TimeD56=telescopio(np.array(TimeD56),True)
    return
def PlotSin(dir,FilePath):
    """-Separa gli eventi letti dal De0Nano"""
    logging.info(f'load data {FilePath}...')
    NumTel,Time=np.loadtxt(f'{dir}/{FilePath}',unpack=True)

    telescopio(np.array(Time),True)

##PROCESS
def Process():
    """main process"""

    #dir='Connessioni8canali_20210424'
    dir='PuntiLavPMT_20210416'
    FilePath='EffPMT3T1S35mV.txt'
    EffPMT(dir,FilePath,True,False,3)



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

    Process()
