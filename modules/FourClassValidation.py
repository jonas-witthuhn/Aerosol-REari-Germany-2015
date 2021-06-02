#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:19:43 2020

@author: walther

This is a translation of the matlab script of Jamie Bright:
https://github.com/JamieMBright/time-series-validation.git
"""
import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

def _warning(
    message,
    *args):
    print(message)

warnings.showwarning = _warning


dtype_general = [
    ( 'Om',     'f8' ), # Observation mean
    ( 'Pm',     'f8' ), # Prediction mean
    ( 'N',      'i4' ), # length of measurements and predictions
    ]
dtype_A = [
    ( 'MBD',    'f8' ),
    ( 'RMSD',   'f8' ),
    ( 'MAD',    'f8' ),
    ( 'SD',     'f8' ),
    ( 'R2',     'f8' ),
    ( 'SBF',    'f8' ),
    ( 'U95',    'f8' ),
    ( 'TS',     'f8' )
    ]
dtype_B = [
    ( 'NSE',    'f8' ),
    ( 'WIA',    'f8' ),
    ( 'LCE',    'f8' )
    ]
dtype_C = [
    ( 'KSI',    'f8' ),
    ( 'OVER',   'f8' ),
    ( 'CPI',    'f8' )
    ]


def normalize_metrics(validation_struct):
    ## normalize in the way that 0 is the best value for all metrics
    for key in validation_struct.dtype.names:
        if key =='MBD':
            # best value is 0 , negative and positive deviation is equal bad
            validation_struct[key]=np.abs(validation_struct[key])
        if key =='R2':
            validation_struct[key]=100.*np.abs(1.- validation_struct[key])

        if key in ['SBF','NSE','WIA','LCE']:
            # original best value was 1 -> now 0
            validation_struct[key]=np.abs(1.- validation_struct[key])
    return validation_struct
            
    
    
def FourClassValidation(Observations,Predictions,class_selection=['A','B','C','D']):

    #   References
    #  
    #   [1] Polo J, Zarzalejo LF, Ramirez L, Espinar B. Iterative filtering of
    #   ground data for qualifying statistical models for solar irradiance
    #   estimation from satellite data. Sol. Energy 2006; 80:2407
    #  
    #   [2] Espinar B, Ramirez L, Drews A, Beyer HG, Zarzalejo LF, Polo J, et
    #   al. Analysis of different comparison parameters applied to solar
    #   radiation data from satellite and German radiometric stations. Sol Energy
    #   2009;83:11825.
    #  
    #   [3] Marsaglia G, Tsang WW, Wang J. Evaluating Kolmogorov's Distribution.
    #   J Stat Softw 2003;8:14
    #  
    #   [4] Taylor KE. Summarizing multiple aspects of model performance in a
    #   single diagram. J Geophys Res 2001;106D:718392.
    #  
    #   [5] Correa CD, Lindstrom P. The mutual information diagram for
    #   uncertainty visualization. Int J Uncertain Quantif 2013;3:187201.
    
    ### checks
    permissible_classes=['A','B','C','D']
    class_selection=np.array(class_selection)
    for c in class_selection:
        if c not in permissible_classes:
            raise ValueError("User defined class '{}' is not a valid class".format(c))
    
    Observations=np.array(Observations)
    Predictions=np.array(Predictions)
    
    
    if np.any(~np.isreal(Observations)) or np.any(~np.isreal(Predictions)):
        raise ValueError("Observations and Predictions must consist of real numbers")
    if Observations.shape != Predictions.shape:
        raise ValueError("Observations and Predictions must be of equal size")
    
    ### collect required datatypes
    dtype=dtype_general.copy()
    for c in class_selection:
        if c=='A':
            dtype.extend(dtype_A)
        elif c=='B':
            dtype.extend(dtype_B)
        elif c=='C':
            dtype.extend(dtype_C)
        else:
            continue # Class D is visualizations, therefore no statistics are calculated

    ### if only one time series, extend dimensions regardless
    Oshape=Observations.shape
    if len(Oshape)==1:
        Observations=Observations[:,np.newaxis]
        Predictions=Predictions[:,np.newaxis]
        
        
    Oshape=Observations.shape
    ### initialize validationstruct
    validation_struct=np.array([tuple(np.zeros(len(dtype)))]*Oshape[-1],
                               dtype=dtype)
    
    ### loop through each time series, assuming each column is a unique site to validate.
    for i in range(Oshape[-1]):
        O=Observations[:,i]
        P=Predictions[:,i]
        
        # remove nans
        ind=~np.isnan(O)*~np.isnan(P)
        O=O[ind]
        P=P[ind]
        
        
        
        # define common usages
        Om=np.mean(O)
        Pm=np.mean(P)
        N=len(O)
        
        validation_struct['Om'][i]=Om
        validation_struct['Pm'][i]=Pm
        validation_struct['N'][i]=N
        
        
        ## Class A - indicators of dispersion
        # These are the indicators that the majority of readers should be most
        # familiar with. They are all expressed here in percent (of Om) rather than
        # in absolute units (W/m2 for irradiances, or MJ/m2 or kWh/ m2 for
        # irradiations) because non-expert stakeholders can much more easily
        # understand percent results. In any case, stating the value of Om in all
        # validation results allows the experts to convert back the percent figures
        # into absolute units if they so desire. Formulas in this section are well
        # established and do not need further references.
        if 'A' in class_selection:
            pOm=100./Om
            # A.1 Mean bias difference (MBD)
            MBD=pOm * np.mean(P-O)
            validation_struct['MBD'][i] =  MBD
            # A.2 Root mean square difference (RMSD)
            RMSD=pOm * np.sqrt(np.mean((P-O)**2))
            validation_struct['RMSD'][i]=  RMSD
            # A.3 Mean absolute difference (MAD)
            MAD=pOm * np.mean(np.abs(P-O))
            validation_struct['MAD'][i] =  MAD
            # A.4 Standard deviation of the residual (SD)
            SD=pOm * np.sqrt( np.sum(N*(P-O)**2)-np.sum((P-O)**2))/float(N)
            validation_struct['SD'][i]  =  SD
            # A.5 Coefficient of determination (R2)
            # A.6 Slope of best-fit line (SBF)
            if np.sum((O-Om)**2)==0:
                warnings.warn(">>FourClassValidation.py - Warning: Observations[%d] are constant - R2 and SBF cannot be used"%(i),UserWarning)
                SBF=np.nan
                R2=np.nan
            else:
                R2= 1.-np.sum((O-P)**2)/ np.sum((O-Om)**2)
                SBF=np.sum((P-Pm)/(O-Om)) / np.sum((O-Om)**2)
            validation_struct['R2'][i]  = R2
            #validation_struct['SBF'][i] = SBF
            # A.7 Uncertainty at 95#
            validation_struct['U95'][i] = 1.96* np.sqrt(validation_struct['SD'][i]**2 + validation_struct['RMSD'][i]**2)
            # A.8 t-statistic (TS)
            if RMSD==MBD:
                TS=0.
            else:
                TS= np.sqrt(float(N-1)*MBD**2/(RMSD**2-MBD**2))
            validation_struct['TS'][i]  = TS
       
        ## Class B - Indicators of overall performance
        # These are indicators that are less common in the solar field than those
        # of Class A. They convey relatively similar information as those of Class
        # A, with the cosmetic advantage that a higher value indicates a better
        # model.
        if 'B' in class_selection:
            if np.sum((O-Om)**2)==0:
                warnings.warn(">>FourClassValidation.py - Warning: Observations[%d] are constant - NSE and LCE cannot be used"%(i),UserWarning)
                NSE=np.nan
                LCE=np.nan 
            else:
                NSE=1.-np.sum((P-O)**2)/np.sum((O-Om)**2)
                LCE=1.-np.sum(np.abs(P-O))/np.sum(np.abs(O-Om))
            # B.1 Nash-Sutcliffe's efficiency (NSE)
            validation_struct['NSE'][i]=NSE
            # B.2 Willmotts's  index of agreement (WIA)
            WIA=1.-np.sum(P-O)**2/np.sum(np.abs(P-Om)+np.abs(O-Om))**2
            validation_struct['WIA'][i]=WIA
            # B.3 Lagates's coefficient of efficiency (LCE)
            validation_struct['LCE'][i]=LCE
            # LCE and NSE vary between 1 for perfect agreement and -inf for complete
            # disagreement, whereas WIE varies only between 1 and 0.
        
        ## Class C
        # The goal is to compare one or more cumulative frequency distribution
        # of modeled data to that of a reference dataset. Can one or more single
        # number provide a measure of the similitude between two or more
        # distributions? Substantial progress in that direction resulted from an
        # initial study by Polo et al. [1],who proposed to use the Kolmogorov
        # Smirnov test when comparing different cumulative distribution functions
        # (CDFs), because of its advantage of being non- parametric and valid for
        # any kind of CDF. Espinar et al. [2] developed the method further, now
        # referring to it as the KolmogorovSmirnov test Integral (KSI)
        if 'C' in class_selection:
            if np.sum((O-Om)**2)==0:
                warnings.warn(">>FourClassValidation.py - Warning: Observations[%d] are constant - KSI, OVER and CPI cannot be used"%(i),UserWarning)
                KSI=np.nan
                CPI=np.nan
                OVER=np.nan
            else:
                # C.1 Kolmogorov-Smirnov test Integral (KSI)
                # irradiance must be binned into x by intervals of n
                # KSI is 0 if the two distributions being compared can be considered
                # identical in a statistical sense.
                xbins=np.linspace(np.nanmin(O),np.nanmax(O),15)
                xmin,xmax= xbins[0],xbins[-1]
                edges=np.array(list(xbins)+[np.inf])
                Od=np.histogram(O,edges)[0]/N
                Pd=np.histogram(P,edges)[0]/N
                # absolute difference between the two normalised distributions
                Dn=np.abs(Od-Pd)
                # pure function of N obtained from [3], though simplified to constant phi(N) \approx 1.63
                Dc = 1.63/N**0.5
                Ac = Dc*(xmax-xmin) 
                fun1 = interp1d(xbins,Dn)
                KSI=(100./Ac)*quad(fun1,xmin,xmax)[0]
                # C.2 Relative frequency of exeedence situations (OVER)
                # The OVER test is the same as the KSI test, but only for those
                # bins that exceed the limit defined by Dc. This is useful when the
                # normalised distribution of modelled data points in specific bins
                # exceeds the critical limit that would make it statistically
                # undistinguisable from the reference distribution.
                # OVER is 0 if the normalised distribution always remains below Dc.
                # OVER can be null indicating that the distribution of the
                # predictions generally respect those of the predictions.
                Dnc=Dn.copy()
                Dnc[Dnc<Dc]=0
                fun2= interp1d(xbins, Dnc)
                OVER=(100./Ac)*quad(fun2,xmin,xmax)[0]
                # C.3 Combined Performance Index (CPI)
                # The interest of CPI is that it combines conventional information
                # about dispersion and biase (through RMSD) with information about
                # distribution likenesses (through KSI and OVER), whilst maintaining a
                # high degree of discrimination between the different models. This
                # feature is of course highly desireable when comparing differnet
                # models of similar performance. This is arguably the most significant
                # statistic to compare different model performance.
                RMSD=(100./Om)*np.sqrt(np.mean((P-O)**2))
                CPI=(KSI+OVER+2.*RMSD)/4.
            validation_struct['KSI'][i]=KSI
            validation_struct['OVER'][i]=OVER
            validation_struct['CPI'][i]=CPI
            
    # Class D
    # This category is completely different from the three previous ones
    # because the goal here is to obtain a visualization rather than summary
    # statistics in the form of a few numbers
    if 'D' in class_selection:
        ## TODO
        
        # The first recommended plots for class D is a Taylor diagram detailed
        # by KE Taylor [4] that combines RMSD, SD and R2 into a single polar
        # diagram. It is ideal for comparing the performance of many different
        # models.
        
        # The second suggestion is a Mutual Information Diagram, which is a
        # revision of the Taylor diagram proposed by [5].
        
        # The box plot is another decent variation for demonstrating
        # performance at different sites.
        pass
    return validation_struct
