import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging

Cs1=    np.array([168027 , 168023 , 168738 , 169337])
dCs1=  np.array([41 , 32 , 36 , 33])
ErrCs1=np.sqrt(np.std(Cs1)**2+np.mean(dCs1)**2)
Cs1=np.mean(Cs1)

Co11=   np.array([ 290342 , 291295 , 290651 , 291577])
dCo11=  np.array([ 141 , 118 , 93 , 120])
ErrCo11=np.sqrt(np.std(Co11)**2+np.mean(dCo11)**2)
Co11=np.mean(Co11)
Co21=   np.array([ 326215, 326854 , 326133 , 326533  ])
dCo21=  np.array([ 205, 162 , 120 ,180])
ErrCo21=np.sqrt(np.std(Co21)**2+np.mean(dCo21)**2)
Co21=np.mean(Co21)

Ne1=    np.array([ 312010 , 311800 , 311470 , 311690])
dNe1=   np.array([140 , 110 , 100 , 150])
ErrNe1=np.sqrt(np.std(Ne1)**2+np.mean(dNe1)**2)
Ne1=np.mean(Ne1)

Cs2 =   np.array([56157 , 56116 , 55809])
dCs2=   np.array([13 , 12 , 6])
ErrCs2=np.sqrt(np.std(Cs2)**2+np.mean(dCs2)**2)
Cs2=np.mean(Cs2)

Co12=   np.array([98718 , 98952 , 98910 , 98488])
dCo12=  np.array([55 , 46 , 34 , 57])
ErrCo12=np.sqrt(np.std(Co12)**2+np.mean(dCo12)**2)
Co12=np.mean(Co12)

Co22=   np.array([111940 , 111971 , 111994 , 111573])
dCo22=  np.array([83 , 79 , 50 , 97])
ErrCo22=np.sqrt(np.std(Co22)**2+np.mean(dCo22)**2)
Co22=np.mean(Co22)



Ne2=  np.array([107228 , 107195 , 107080 , 106964])
dNe2=   np.array([49 , 37 , 31 , 41])
ErrNe2=np.sqrt(np.std(Ne2)**2+np.mean(dNe2)**2)
Ne2=np.mean(Ne2)

print('PMT1')
Cs=168573
Co1=290972
Ne=311681
Co2=326436
Na=128105
ErrCs=18
ErrCo1=57
ErrNe=62
ErrCo2=81
dNa=22

"""
print('PMT2')
Cs=55932
Co1=98822
Ne=107099
Co2=111916
Na=43103
ErrCs=5.4
ErrCo1=22
ErrNe=20
ErrCo2=35
dNa=7.8
dNa=0
"""



def process(Type):

    ###in ordine: Cesio,Cobalto1,Neon,Cobalto2 [eV]
    ##PMT1
    y  = np.array([Cs,Co1,Ne,Co2])
    dy = np.array([ErrCs,ErrCo1,ErrNe,ErrCo2])
    #Na,dNa=128150,0
    """
    #print(dy/y)
    #print(y)
    #print(dy)

    ##PMT2
    y  = np.array([Cs2,Co12,Ne2,Co22])
    dy = np.array([ErrCs2,ErrCo12,ErrNe2,ErrCo22])
    Na,dNa=43210,0
    print(dy/y)
    print(y)
    print(dy)
    """

    x=np.array([661.6,1173.2,1274.5,1332.5])


    """-best fit"""
    ####   function
    if Type=='Parabola':
        init=[0,0,0]
        def f(x,a,b,c):
            function = a*x**2 + b*x + c
            return function
    if Type=='ParabolaZero':
        init=[0,0]
        def f(x,a,b):
            function = a*x**2 + b*x
            return function
    if Type=='Retta':
        init=[200,0]
        def f(x,a,b):
            return a*x + b
    if Type=='RettaZero':
        init=[200]
        def f(x,a):
            return a*x

    from scipy.optimize import curve_fit
    #   set init value of parameter
    sigma=dy
    #sigma=numpy.sqrt(dy**2+ ((2/14.48)/numpy.sqrt(12))**2)
    w=1/sigma**2        #inverso della varianza
    #   fitting
    pars,covm = curve_fit(f,x,y,init,sigma,absolute_sigma=False)
    #   Calculate the chisquare for the best-fit function
    if Type=='RettaZero':
        chi2 = ((w*(y-f(x,pars[0]))**2)).sum()
    if (Type=='Retta') | (Type=='ParabolaZero'):
        chi2 = ((w*(y-f(x,pars[0],pars[1]))**2)).sum()
    if Type=='Parabola':
        chi2 = ((w*(y-f(x,pars[0],pars[1],pars[2]))**2)).sum()
    ndof=len(x)-len(init)
    #output
    a,da = pars[0], np.sqrt(covm[0,0])
    print(f'a = {pars[0]} +/- {np.sqrt(covm[0,0])}')
    if (Type != 'RettaZero'):
        b,db = pars[1], np.sqrt(covm[1,1])
        print(f'b =  {pars[1]}+/- {np.sqrt(covm[1,1])}')
    if Type=='Parabola':
        c,dc = pars[2], np.sqrt(covm[2,2])
        print(f'c =  {pars[2]}+/- {np.sqrt(covm[2,2])}')
    print(f'chi2 = {chi2}, ndof = {ndof}')

    ##MISURE DI MASSA ELETTRONE
    if Type=='Parabola':
        rdelta  =  np.sqrt(b**2-4*a*(c-Na))
        me  =   (-b+rdelta)/(2*a)
        dsdb= (-1+b/rdelta)/(2*a)
        dsdc= -1/rdelta
        dsdy= +1/rdelta
        dsda= (b-rdelta)/(2*a**2) + (-c+Na)/(rdelta*a)
        dme=np.sqrt(dsdb**2*covm[1,1]+dsdc**2*covm[2,2]+dsda**2*covm[0,0]+2*(dsdb*dsda*covm[0,1])+2*(dsda*dsdc*covm[0,2])+2*(dsdb*dsdc*covm[2,1]))
    if Type=='ParabolaZero':
        rdelta  =  np.sqrt(b**2-4*a*(-Na))
        me=   (-b+rdelta)/(2*a)
        dsdb= (-1+b/rdelta)/(2*a)
        dsdy= +1/rdelta
        dsda= (b-rdelta)/(2*a**2) + (Na)/(rdelta*a)
        dme=np.sqrt(dsdb**2*covm[1,1]+dsda**2*covm[0,0]+2*(dsdb*dsda*covm[0,1]))
    if Type=='Retta':
        me   = (Na-b)/a
        dfdy = dfdb = 1/a
        dfda = (Na-b)/a**2
        dme  = np.sqrt(dfdb**2*db + dfda**2*da + 2*covm[0,1]*dfda*dfdb)
    if Type=='RettaZero':
        me  = Na/a
        dme = np.sqrt( (dNa/Na)**2 + (da/a)**2 )*me
    print(f'me: {me} +- {dme}')



    #Plot
    #plt.figure(1)
    fig, (graph1,graph2)=plt.subplots(2)
    #errorbar Data Plot

    graph1.errorbar(x,y,dy,linestyle='',color='b',marker='.')

    #   plot fit_curve
    xx=np.linspace(min(x),max(x),1000)
    if Type=='RettaZero':
        graph1.plot(x,f(x,pars[0]), color='red')
        r = (y-f(x,pars[0]))
    if (Type=='Retta')|(Type=='ParabolaZero'):
        graph1.plot(x,f(x,pars[0],pars[1]), color='red')
        r = (y-f(x,pars[0],pars[1]))
    if Type=='Parabola':
        graph1.plot(x,f(x,pars[0],pars[1],pars[2]), color='red')
        r = (y-f(x,pars[0],pars[1],pars[2]))


    fig.suptitle(f'Fit con {Type}')
    graph1.set(xlabel='', ylabel='Energia [u.a.]')
    graph2.set(xlabel='Energia [keV]', ylabel='res.')
    graph2.errorbar(x,r,dy,fmt='.');graph2.axhline(0, color='k', linewidth=0.5)

    plt.show()



if __name__=='__main__':

    parser = argparse.ArgumentParser(usage='Letter Counter of a text')

#    parser.add_argument('DataFile', type=str, help = 'path to input file_data')
    parser.add_argument('-info', default=False, action='store_true', help='print info messages')
    parser.add_argument('Type',  help='Indicate the curve fit type', choices=['Retta','Parabola','ParabolaZero','RettaZero'])

    args = parser.parse_args()

    if args.info:
        logging.basicConfig(level=logging.INFO)

    process(args.Type)
