import numpy as np
import matplotlib.pyplot as plt
from Time_Of_Flight import TripleRecognize,QuadRecognize
import logging
import itertools

logging.basicConfig(level=logging.INFO)

def RicercaEvento(NumEve):
    """-Ricerca gli indici di inizio e fine dei dati dell'evento NumEve"""
    #togli l'offset per i dati regolari
    RecordLine=NumEve*1037
    StartData=RecordLine+7
    EndData=RecordLine+1037

    return  StartData, EndData

"""-Si può utilizzare questo codice per visualizzare gli eventi che si vuole di una determinata presa dati"""



"""-Cambiare directory file da cercare e evento da plottare"""
ch=3
dir='2020_12_15_doppie'
EventSearch=3



###Process
"""-Apertura file"""
logging.debug('Apro file')
file1=open(f'{dir}\wave0.txt')
file2=open(f'{dir}\wave1.txt')
file3=open(f'{dir}\wave2.txt')
if ch == 4:
    file4=open(f'{dir}\wave3.txt')


StDa,EndDa=RicercaEvento(EventSearch)
"""-Inizio lettura"""
Data1=Data2=Data3=np.array([],dtype=int)
if ch == 4:
    Data4=np.array([],dtype=int)

print('Start!!!')
logging.info(f'Search event {EventSearch}..')
CurrLine=0

        ### CAMBIARE SE CH = 4 ####
for line1,line2,line3 in itertools.zip_longest(file1,file2,file3):
    #print(CurrLine)
    if CurrLine==EndDa:
        break
    if (CurrLine>=StDa) & (CurrLine<EndDa):
        #print(f'{line1}\n{line2}\n{line3}')
        #incrementa la lista di dati dell'evento ricercato
        Data1=np.append(Data1,np.array([int(line1)]))
        Data2=np.append(Data2,np.array([int(line2)]))
        Data3=np.append(Data3,np.array([int(line3)]))
        if ch ==4:
            Data4=np.append(Data4,np.array([int(line4)]))
    CurrLine+=1

file1.close()
file2.close()
file3.close()
if ch==4:
    file4.close()

time=np.arange(len(Data1))


"""-Coincidence"""
logging.info(f'Trovato evento{EventSearch}')
try:
    if ch==3:
        Ts1,C1,Ts2,C2,Ts3,C3=TripleRecognize(Data1,Data2,Data3)
    if ch==4:
        Ts1,C1,Ts2,C2,Ts3,C3,Ts4,C4=QuadRecognize(Data1,Data2,Data3,Data4)
    print(Ts1,Ts2,Ts3);
    if ch==4: print(Ts4)

except:
    print('no Coinc')


"""-plot"""
time=np.arange(len(Data1))

plt.plot(time,Data1,label='PMT1')
plt.plot(time,Data2,label='PMT2')
plt.plot(time,Data3,label='PMT3')

if ch==4:
    plt.plot(time,Data4,label='PMT4')

plt.legend()
plt.show()
