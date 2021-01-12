""" Ho copiato in gran parte il codice di lab 2 e modificato
    per essere più leggibile e riutilizzabile
    """

import numpy
import matplotlib
import pylab
import argparse
import logging


def process(DataFilePath, residui):
    """A simple program to read file data and fit it
    """
    #Import data
    logging.info('Reading data..')
    x,y,dx,dy=pylab.loadtxt(DataFilePath, unpack=True)

    #best fit
    #   function
    def f(x,a,b):
        return a+2/b*x

    from scipy.optimize import curve_fit
    #   set init value of parameter
    init=[20,-15]
    #   set error dy or the propagate error with dx and b
    #sigma=dy
    sigma=numpy.sqrt(dy**2+ ((2/14.48)/numpy.sqrt(12))**2)
    w=1/sigma**2        #inverso della varianza
    #   fitting
    pars,covm = curve_fit(f,x,y,init,sigma,absolute_sigma=False)
    #   Calculate the chisquare for the best-fit function
    chi2 = ((w*(y-f(x,pars[0],pars[1]))**2)).sum()
    ndof=len(x)-len(init)
    #output
    for i in range(len(pars)):
        print(f'{pars[i]:.2f} +/- {covm[i,i]:.2f}')
    print(f'chi2 = {chi2:.2f}, ndof = {ndof}')
    print(f'a = {pars[0]:.2f} +/- {numpy.sqrt(covm[0,0]):.2f}')
    print(f'b =  {pars[1]:.2f}+/- {numpy.sqrt(covm[1,1]):.2f}')
    print(f'norm cov = {covm[0,1]/(numpy.sqrt(covm[0,0]*covm[1,1])):.2f}')


    #Plot
    pylab.figure(1)
    #errorbar Data Plot
    pylab.errorbar(x,y,dy,linestyle='',color='b',marker='.')
    pylab.xlabel('Posizione PMT-3 [cm]');pylab.ylabel('Intervallo di tempo PMT1-2 [ns]');pylab.title('V Luce Barra Scintillante')
    #   plot fit_curve
    xx=numpy.linspace(min(x),max(x),1000)
    pylab.plot(x,f(x,pars[0],pars[1]), color='red')

    if residui is not True:
        pylab.figure(2)
        r = (y-f(x,pars[0],pars[1]))/sigma
        pylab.title('Residui')
        pylab.ylabel('Norm. res.');pylab.xlabel('Posizione PMT-3 [cm]')
        #pylab.ylim((-.9,.9))
        pylab.errorbar(x,r,dy,fmt='.');pylab.axhline(0, color='k', linewidth=0.5)

    #   show
    pylab.show()



if __name__=='__main__':

    parser = argparse.ArgumentParser(usage='Letter Counter of a text')

    parser.add_argument('DataFile', type=str, help = 'path to input file_data')
    parser.add_argument('-info', default=False, action='store_true', help='print info messages')
    parser.add_argument('-residui', default=True, action='store_false', help='do not view the residui plot')

    args = parser.parse_args()

    if args.info:
        logging.basicConfig(level=logging.INFO)

    process(args.DataFile, args.residui)
