'''
Created on 2013-1-25
=====================================The tsunami sediment model3 (function part)================================
#===============================================function dietrichWs=============================================
#Input: 	T- temperature of fluid in deg. C[Array2]
#               D- grain size in mm [Num]
#               rhos- density of sediment in gm/cm**3 [Num]
#               rhow- density of water in gm/cm**3 [Num]
#               Cs- volume concentration of sediment [Num]
#               csf- Corey shape factor, usually taken as 0.7 [Num]
#               P- roundness of grains, usaully taken as 3.5 [Num]
#Output         ws- settling velocity (cm/s) [Array2]
#function       KinVisc(T)
#=================================================Function gelfK================================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array1]
#   		h- water depth (m) [Num]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K- parabolic eddy viscosity [Array1]
#===============================================Function kinvisc================================================
#Input: 	T- temperature of fresh water in degrees celsius [Array2]
#Output:        kinvisc- kinematic viscosity in m**2/s[Num]
#==================================================function linearK=============================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array1]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K-linear eddy viscosity profile [Array1]
#================================================Function parabolick============================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array 1]
#   		h- water depth (m) [Num]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K-parabolic eddy viscosity profile [Array1]
#====================================================Function phi2mm============================================
#Input: 	phi- Grain size in phi [Num]
#Output:        mm-Grain size in mm [Num]
#====================================================function Setvel1===========================================
#Input: 	d- grain diamter (m) [Num]
#		kinvisc- kinematic viscosity (m**2/s) [Array2]
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#Output:        ws- particle settling velocity (m/s) [Array2]
#=================================================Function sw_dens0=============================================
#Input: 	S = salinity    [psu      (PSS-78)] [Array2]
#               T = temperature [degree C (IPTS-68)][Array2]
#Output:        dens0 = density  [kg/m**3] of salt water with properties S,T, [Array2]
#               P=0 (0 db gauge pressure)  [Num]
#===================================================Function sw_smow============================================
Input: 	        T = temperature [degree C (IPTS-68)][Array2]
#Output:        dens = density  [kg/m**3] [Array2]
#================================================function tubesetvel============================================
#Input: 	phi- grain size in phi [Num]
#               rhos- density of sediment (or sphere) in gm/cm**3 [Num]
#               T- temperature of water in settling tube, deg. C [Array2]
#               S- salinity of water in settling tube (added to explore real world) [Array2]
#Output:        ws- settling velocity of grain (cm/s) [Array2]
#Functions:     KinVisc(T) and sw_dens0(S,T)
#================================================function UstarCrit=============================================
#Input: 	d- grain diamter (m) [Num]
#		kinvisc- kinematic velocity (m**2/s) [Array2]
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#               shieldscrit- Shield's parameter for initiation of sediment transport [Num]
#		Shield's parameter= ustarcrit^2/((s-1)*g*d)
#Output:        ustarcrit- critical shear velocity for initiation of sediment transport (m/s) [Array2]
#Function:      KinVisc(T)
#===============================================function zoD====================================================
#Inputs:    D- nominal grain diameter (m) [Num]
#	    ustarc- current shear velocity (m/s) (may be a vector) [Num]
#	    rhow- density of water (kg/m3) [Num]
#           s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#           T- temperature of water in settling tube, deg. C [Array2]
#Output:    zosD- bed roughness for saltating layer (m)  *** must be added to zoN (Nikuradse
#	    bed roughness) to get total roughness  (zo=zoN+zosD) [Array2]
#Function:  KinVisc(T),UstarCrit(D,T,s)
#===================================================function zoWS===============================================
#Inputs:    D- nominal grain diameter (m) [Num]
#	    ustarc- current shear velocity (m/s) (may be a vector) [Num]
#	    rhow- density of water (kg/m3) [Num]
#           s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#           T- temperature of water in settling tube, deg. C [Array2]
#Output:        zos- bed roughness for mobile flat bed (m)  *** must be added to zoN (Nikuradse
#		bed roughness) to get total roughness  (zo=zoN+zos)[Array]
#Function: UstarCrit(D,T,s)
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
outp_path = '../output/'
from pylab import *
from types import *
import csv

#===============================================function dietrichWs===========================================
#This is a fuction that calculates settling velocity from Dietrich'82
#Input: 	T- temperature of fluid in deg. C
#               D- grain size in mm
#               rhos- density of sediment in gm/cm**3
#               rhow- density of water in gm/cm**3
#               Cs- volume concentration of sediment
#               csf- Corey shape factor, usually taken as 0.7
#               P- roundness of grains, usaully takefunction.pyn as 3.5
#Used
#               g- acceleration due to gravity, 980 cm/sec**2
#               k- 1/von karmen's constant, 2.5
#Calculated
#               mu- dynamic viscosity of water (poise)
#               musubs- dynamic viscosity of fluid (poise)
#               rhof- density of fluid in gm/cm^3
#               nu- kinematic viscosity of fluid (poise)
#               dstar- dimensionless particle size
#               wstar- dimensionless settling velocity
#Output         ws- settling velocity (cm/s)
#need function  KinVisc
#originate from Bruce Jaffe, USGS, Agu 6,2003
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def dietrichWs(T,D,rhos,rhow,Cs,csf,P):
    g=9.8  # acceleration due to gravity, cm/sec**2
    k=2.5
    nargin=5
    if(nargin==5):# if csf and P are not set, use defaults
        csf=0.7    # set Corey shape factor of grains to default value
        P=3.5   # set roundness of grains to default value
    lu=KinVisc(T)
    mu=KinVisc(T)*rhow*10000         # convert kinematic viscosity in m**2/s to dynamic viscosity in poise
    #calculate dynamic viscosity for water and sediment
    musubs=mu*(1+k*Cs)
    #calculate fluid density
    rhof=rhow+(rhos-rhow)*Cs
    #calculate nu for the sediment conentration
    nu=musubs/rhof
    #calculate settling velocity, converts grain size from mm to cm
    dstar=(rhos-rhof)*g*(D/10)**3/(rhof*nu**2)             # eq. 6, Dietrich '82 (D82), converting size from mm to cm
    dstarl=log10(dstar)
    r1=-3.76715+1.92944*dstarl-0.09815*dstarl**2-0.00575*dstarl**3+0.00056*dstarl**4 #fitted equation for size and density effects, eq. 9 D82
    omc=1.0-csf
    r2=log10(1-omc/0.85)-omc**2.3*tanh(dstarl-4.6)+0.3*(0.5-csf)*omc**2*(dstarl-4.6) #fitted equation for shape effects. eq. 16 D82
    r3=(0.65-(csf/2.83)*tanh(dstarl-4.6))**(1+(3.5-P)/2.5) #fitted equation for roundness effects, eq. 18 D82
    wstar= r3*10**(r1+r2)                                       # eq. 19, D82
    ws=((rhos-rhof)*g*nu*wstar/rhof)**(1.0/3) # rearranged eq. 5 from D82
    return(ws)

#=================================================Function gelfK=============================================
#This is a function that calculates the parabolic eddy viscosity profile using
#Input: 	z- vector with elevations (m) whsedstats.pyere the eddy viscosity is calculated
#   		h- water depth (m)
#		ustrc- current shear velocity (m/s)
#Output:        K- eddy viscosity
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def gelfK(z,h,ustrc):
    vk=0.41
    K=zeros(len(z))
    for i in range(len(z)):
        K[i]=vk*ustrc*z[i]*exp((-z[i]/h)-3.2*(z[i]/h)**2+2.13*(z[i]/h)**3)
    return K

#===============================================Function kinvisc=============================================
#This is a function that calculates kinematic viscosity of fresh water given temperature
#Uses eq. 3.1.4 from Van Rijn, 1993, "Principles of Sediment Transport in Rivers, Estuaries, and Coastal Seas."
#Input: 	T- temperature of fresh water in degrees celsius
#Output:        kinvisc- kinematic viscosity in m^2/s
#originate from Bruce Jaffe, USGS, Jan 30,2001
#Made into python by Hui Tang, Virginia Tech, Jan 30, 2013 tanghui@vt.edu
def KinVisc(T):
    kinvisc=(1.14-0.031*(T-15)+0.00068*(T-15)**2)*10**-6
    return kinvisc

#==================================================function linearK==========================================
#This is a fuction that calculates the linear eddy viscosity profile
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated
#		ustrc- current shear velocity (m/s)
#originate from Bruce Jaffe, USGS, April 3,2sedstats.py001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def linearK(z,h,ustrc):
    K=zeros(len(z))
    vk=0.41
    for i in range(len(z)):
        K[i]=vk*z[i]*ustrc
    return K

#================================================Function parabolick=========================================
#This is a function that calculates the parabolic eddy viscosity profile
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated
#   		h- water depth (m)
#		ustrc- current shear velocity (m/s)
#Output:        K-parabolic eddy viscosity profile
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def parabolicK(z,h,ustrc):
    K=zeros(len(z))
    vk=0.41
    for i in range(len(z)):
        K[i]=vk*z[i]*(1-z[i]/h)*ustrc
    return K

#====================================================Function phi2mm==========================================
#This is a function that converts phi to mm
#Input: 	phi- Grain size in phi
#Output:        mm-Grain size in mm
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def phi2mm(phi):
    mm=2**(-phi)
    return mm

#====================================================function Setvel1=========================================
#This is a fuction calculate settling velocity of sediment.  Formula from
#Ole Madsen's short course at Coastal Seds '99.  These are fits to data from
#Dietrich '82, van Rijn '90, Cheng '97.
#Input: 	d- grain diamter (m)
#		kinvisc- kinematic viscosity (m**2/s)
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water)
#Output:        ws- particle settling velocity (m/s)
#Variation of Setvel allowing diameters smaller than 0.00006 m
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def Setvel1(d,kinvisc,s):
    g=9.80 #gravitational acceleration
    # calculate sstar
    sstar=(d/(4*kinvisc))*sqrt((s-1)*g*d)
    if(d>= 0.0001 and d <= 0.001 ): #diametesedstats.pyr between 0.1 and 1mm
        wstar=1/(5.40/sstar + 0.92)
    elif(d>= 0.0001 and d <= 0.005): #diameter between 0.1 and 5mm
        wstar=1/(5.21/sstar + 0.88)
    ws=wstar*sqrt((s-1)*g*d) #settling velocity in m/s (?)
    return ws

#=================================================Function sw_dens0==========================================
#This is a function that calculates Density of Sea Water at atmospheric pressure using UNESCO 1983 (EOS 1980) polynomial.
#DISCLAIMER:
#This software is provided "as is" without warranty of any kind.
#REFERENCES:
#Unesco 1983. Algorithms for computation of fundamental properties of seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#UNESCO 1983 p17  Eqn(14) Millero, F.J & Poisson, A. INternational one-atmosphere equation of state for seawater. Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
#Input: 	S = salinity    [psu      (PSS-78)]
#               T = temperature [degree C (IPTS-68)]
#Output:        dens0 = density  [kg/m**3] of salt water with properties S,T,
#               P=0 (0 db gauge pressure)
#originate from Phil Morgan, Sep 05,1992 morgan@ml.csiro.au
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def sw_dens0(S,T):
    nargin=2
    if(nargin!=2):
        print('error:sw_dens0.m: Must pass 2 parameters')
    else:
        #if (type(S)=='float' or type(T)=='float'):
        b0=8.24493e-1
        b1=-4.0899e-3
        b2=7.6438e-5
        b3=-8.2467e-7
        b4=5.3875e-9
        c0=-5.72466e-3
        c1=+1.0227e-4
        c2=-1.6546e-6
        d0=4.8314e-4
        dens1=sw_smow(T)
        dens=dens1+(b0+(b1+(b2+(b3+b4*T)*T)*T)*T)*S+(c0+(c1+c2*T)*T)*S*sqrt(S)+d0*S**2
    return(dens)

#===================================================Function sw_smow==========================================
#This is a function that calculates Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
#DISCLAIMER:
#This software is provided "as is" without warranty of any kind.
#REFERENCES:
#Unesco 1983. Algorithms for computation of fundamental properties of seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#UNESCO 1983 p17  Eqn(14) Millero, F.J & Poisson, A. INternational one-atmosphere equation of state for seawater. Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
#Input: 	T = temperature [degree C (IPTS-68)]
#Output:        dens = density  [kg/m^3]
#originate from Phil Morgan, Sep 05,1992 morgan@ml.csiro.au
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def sw_smow(T):
    #check input arguments
    #test inputs
    nargin=1
    if(nargin!=1):
        print('error:sw_smow.m: Only one input argument allowed')
    else:
        a0=999.842594
        a1=6.793952e-2
        a2=-9.095290e-3
        a3=1.001685e-4
        a4=-1.120083e-6
        a5=6.536332e-9
        dens=a0+(a1+(a2+(a3+(a4+a5*T)*T)*T)*T)*T
    return(dens)

#================================================function tubesetvel=========================================
#This is a function to determine settling velocity in USGS settling
#tubes from grain size output of sedsize (phi).  Uses Gibbs'
#settling formula (from excel spreadsheet provided by John
#Penscil, a contrator who wrote the programs for Mike
#Torresan).  This spreadsheet is GibbsEquations.xls.  I obtained
#a copy in 2003.
#Input: 	phi- grain size in phi
#               rhos- density of sediment (or sphere) in gm/cm**3
#               T- temperature of water in settling tube, deg. C
#               S- salinity of water in settling tube (added to explore real world)
#Used:
#               g- acceleration due to gravity (cm/s**2)
#               rhof- density of fluid in settling tube (gm/cm**3)
#               kv- kinematic viscosity of fluid in settling tube(m/s**2)
#               dynvisc- dynamic viscosity of fluid in settling tube (poise)
#               r- radius of grain (or sphere) in mm
#Output:        ws- settling velocity of grain (m/s)
#only need to input grain size in phi, will use defaults of
#rhos=2.65, T=21, S=0
#Need functions KinVisc and sw_dens0sedstats.py
#originate from Bruce Jaffe, USGS, Aug 1,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def tubesetvel(phi,rhos,T, S):
    g=979.9 #gravitational acceleration in cm/s^2
    nargin=1
    if(nargin==1): #use defaults if only phi size is provided
        rhos=2.65 #quartz density
        T=11.6  #21 deg. C water temperature
        S=0#fresh water in settling tubes
    #calculate dynamic viscosity of water in settling tube (poise)
    kv=KinVisc(T)      # m/s^2
    #kv=1.004e-6
    #print kv
    rhof=sw_dens0(S,T) # kg/m^3,
    #print rhof
    dynvisc=10*kv*rhof #dynamic viscosity, converted to poise (dyne-sec/cm^2)
    #print dynvisc
    rhof=rhof/1000.0#convert fluid density to gm/cm^3 for use in formula
    #reproduce spreadsheet columns and calculate settling velocity
    r=(1.0/(2.0**phi))/20.0
    term= sqrt((9.0*dynvisc**2)+(g*r**2*rhof*(rhos-rhof)*(0.015476+(0.19841*r))))
    ws=(term-(3*dynvisc))/((rhof)*(0.011607+(0.14881*r)))/100.0
    return ws

#================================================function UstarCrit==========================================
#This is a fuction that calculates the critical shear velocity for initiation of sediment transport. Formula from Ole Madsen's short course
#at Coastal Seds '99.
#Input: 	d- grain diamter (m)
#		kinvisc- kinematic velocity (m**2/s)
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water)
#               shieldscrit- Shield's parameter for initiation of sediment transport
#		Shield's parameter= ustarcrit^2/((s-1)*g*d)
#Output:        ustarcrit- critical shear velocity for initiation of sediment transport (m/s)
#Function:      KinVisc(T)
#originate from Bruce Jaffe, USGS, Fed 8,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def UstarCrit(d,kinvisc,s):
    #kinvisc=KinVisc(T)
    if(type(kinvisc)=='float'):
        mK=len(kinvisc)
        nK=len(kinvisc[0][:])
        sstar=zeros(shape=(mK,nK))
        ustarcrit=zeros(shape=(mK,nK))
        for i in range(mK) :
            for j in range(nK):
                g=9.80
                sstar[i][j]=d/(4*kinvisc[i][j])*sqrt((s-1)*g*d)
                if (sstar[i][j]>= 0.8):
                    shieldscrit=0.095*sstar[i][j]**(2.0/3.0) + 0.056*(1-exp((-sstar[i][j]**(3.0/4.0))/20.0))
                else:
                    shieldscrit=0.1*sstar[i][j]**(-2.0/7.0)
                ustarcrit[i][j]=sqrt((s-1)*g*d)*sqrt(shieldscrit) # critical shear velocity in m/s

    else:
        g=9.80 #gravitational acceleration
        #calculate sstar
        sstar=d/(4*kinvisc)*sqrt((s-1)*g*d)
        #print kinvisc
        #calculate Shield parameter for initiation of sediment transport
        if (sstar >= 0.8):
            shieldscrit=0.095*sstar**(-2.0/3.0) + 0.056*(1-exp((-sstar**(3.0/4.0))/20.0))
        else:
            shieldscrit=0.1*sstar**(-2.0/7.0)
        ustarcrit=sqrt((s-1)*g*d)*sqrt(shieldscrit) # critical shear velocity in m/s
    return (ustarcrit)
#===============================================function zoD=================================================
#This is a fuction that calculates zosD, Dietrich's 1982 expression for bottom
#roughness created by saltating layer.
#Inputs:    D- nominal grain diameter (m)
#	    ustarc- current shear velocity (m/s) (may be a vector)
#	    rhow- density of water (kg/m3)
#Used:      alpha- constant, 0.077, from Muddy Creek, Wyoming data
#           (averages about 0.10, can be a max of 0.13 according
#           to Wiberg and Rubin 1989)
#Output:    zosD- bed roughness for saltating layer (m)  *** must be added to zoN (Nikuradse
#	    bed roughness) to get total roughness  (zo=zoN+zosD)
#References:Dietrich, W.D., Flow, boundary shear stress, and sediment transport in a river meander,
#           Ph.D. dissertation, 261 p., Univ. of Wash., Seattle, 1982
#           Wieberg, P. L. and Rubin, D. M., 1989, Bed roughness produced by
#	    saltating sediment. J. Geophys. Res., V. 94, no. C4, pp. 5011-5016.
#originate from Bruce Jaffe, USGS, Fed 8,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu

def zosD(T,D,ustarc,rhow,s):
    alpha=0.077         # constant, set using Muddy Creek, Wyoming data
    s=2.65
    kinvisc=KinVisc(T)
    ustarcrit=UstarCrit(D,T,s)
    alpha=0.077
    taub=rhow*ustarc**2
    mU=len(ustarcrit)
    nU=len(ustarcrit[0][:])
    taucrit=zeros(shape=(mU,nU))
    Tstar=zeros(shape=(mU,nU))
    zosD=zeros(shape=(mU,nU))
    for i in range(mU):
        for j in range(nU):
            taucrit[i][j]=rhow*ustarcrit[i][j]**2
            Tstar[i][j]=taub/taucrit[i][j]
            zosD[i][j]=alpha*D*(0.6*Tstar[i][j]/(1+0.2*Tstar[i][j])) #Dietrich 1982 from Wiberg and Rubin, 1989 equation 3
    return zosD

#===================================================function zoWS===========================================
#This is a fuction that calculates zos, the bed roughness parameter for mobile
#flat bed, using the formula from Wieberg and Rubin, 1989.
#ws stands for "Wiberg and Smith"
#Input: 	D- nominal grain diameter (m)
#		ustarc- current shear velocity (m/s) (may be a vector)
#		ustarcrit- critical shear velocity for initiation of motion (m/s) (may be a vector)
#		rhow- density of water (kg/m3))
#Used:          gammawWS- a constant determined by Wieberg and Rubin to fit the data
#		taub- bottom shear stress (N/m2)
#		taucrit- critical shear stress for initiation of motion (N/m2)
#		Tstar- transport stage, ratio of taub to taucrit
#		delb- average saltation height
#		a1- constant used in fomulation of average saltation height
#		a2- constant used in fomulation of average saltation height
#Output:        zos- bed roughness for mobile flat bed (m)  *** must be added to zoN (Nikuradse
#		bed roughness) to get total roughness  (zo=zoN+zos)
#Reference:     Wieberg, P. L. and Rubin, D. M., 1989, Bed roughness produced by
#		saltating sediment. J. Geophys. Res., V. 94, no. C4, pp. 5011-5016.
#originate from Bruce Jaffe, USGS, Mar 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def zoWS(D,ustarc,ustarcrit,rhow):
    gammaWS = 0.056
    a1=0.68
    a2=0.0204*(log(D*100))**2+0.0220*log(D*100)+0.0709 	#equation 6, D converted to cm, Wiberg and Rubin, 1989
    taub=rhow*ustarc**2
    taucrit=rhow*ustarcrit**2
    Tstar=taub/taucrit
    delb=D*a1*Tstar/(1+a2*Tstar)
    zos=gammaWS*delb   # equation 7 Wiberg and Rubin, 1989
    return zos

from pylab import *
from numpy import *
# from function import *
# from sedstats import *
#============================================Subfunction grading==============================================
#This is a function to calculate grading of sediment based on concentration profiles for different grain size
#classes.
#Input:             sl:Sediment load for each size class (columns) for given z(rows)
#                   phi:size classes in phi(classes=columns of C),row vector
#                   z:elevation (m) above bottom of sediment loads (elevations=row of sl), column vector
#                   intv:interval (m) of sediment texture stats, maybe a set interval or specified(e.g thickness
#                        1cm,2cm,1cm,4cm intv[0.01 0.03 0.04 0.08]),default is 0.01 if not specified
#                   porosity:porosity of deposit,default is 0.35 if not specified
#Function used:     tubesetvel.py, sedstats.py
#Note:set defaults for intervals and porosity if not specified
#originate from Bruce Jaffe 9/29/2003
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def grading(sl,phi,z,intv,porosity,ni):
##============================================Calculate part============================================
    numrow_sl=len(sl)
    numcol_sl=len(sl[0]) #get matrix sl size
    hitsize=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with size of particles hitting bed
    hitload=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with concentration of particles hitting bed
    hitheight=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with height of particles hitting bed
    if((len(z)!=numrow_sl) or (len(phi)!=numcol_sl)):
       print('Number of size classes or number of elevations do not match profiles')
    else:
       z=z.reshape(len(z),1) #z must be a column vector
       #calculate time for sediment to reach reach bottom
       setvel=0.01*tubesetvel(phi,rhos=2.65,T=21.6,S=0)#calculate settling velocity (m/s)
       setvel=setvel.reshape(1,len(setvel))
       numrow_setvel=len(setvel)
       numcol_setvel=len(setvel[0])
       if(numrow_setvel>numcol_setvel):
           setvel=setvel #make sure setvel is a row vector
       settlingtime=zeros(shape=(numrow_sl,numcol_setvel))
       stime=zeros(shape=(numrow_sl*numcol_setvel,3))
       k=0
       #print int(numcol_setvel)
       for i in range(numrow_sl):
           for j in range(numcol_setvel):
               settlingtime[i][j]=z[i]/setvel[0][j]
                      #stime[k][0]=settlingtime[i][j]
       stime[:,0]=settlingtime.reshape(int(numrow_sl*numcol_setvel),)
       #print shape(stime)
#***********************************sort by time it takes to hit bottom*****************************************
       stime1=sorted(stime[:,0])
       ind=argsort(stime[:,0],axis=0)
       for i in range(len(ind)):
           stime[k][1]=int(int(ind[i])/int(numcol_setvel))
           stime[k][2]=int(int(ind[i])%int(numcol_setvel))
           k=k+1
#**********************************find size of sediment hitting the bed****************************************
       #phiind=floor(((stime[:,2]-0.001)/numrow_sl))+1 # find size of sediment
       phiind=stime[:,2] # find size of sediment
       hitsize=zeros(len(phiind))
       for i in range(len(phiind)):
           hitsize[i]=phi[phiind[i]]    # keep track of order sediment size hitting the bed
#***********************************find elevations that sediment started from**********************************
       zind=stime[:,1]
       #print ind
       hitheight=zeros(len(zind))
       for i in range(len(zind)):
           hitheight[i]=z[zind[i]] #keep track of order of height sediment came from
       hitload=zeros(len(zind))
       for i in range(len(zind)):
           hitload[i]=sl[zind[i],phiind[i]]
       thickness=cumsum(hitload/(1-porosity))#cumulative thickness of deposit
       #print sum(thickness)
#********************************************determine interval breaks*******************************************
       thickness2=zeros(len(thickness))
       intbind=zeros(len(thickness))
       indthick=zeros(len(intv))
       if(len(intv)==1):#calculate interval boundaries for set intervals
           for intbound in range(intv,thickness[len(thickness)]-intv/2,intv):
               for i in range(len(thickness)):
                   if(thickness[i]<intv[intbound]):
                       thickness2[j]=thickness[i]
                       j=j+1
               intbind[round(intbound/intv)]=max(thickness2)
           intbind2=zeros(len(intbind)+2)
           intbind2[0]=0
           i=1
           while(i<len(intbind)+1):
               intbind2[i]=intbind[i-1]
               i=i+1
           intbind[i]=len(thickness)#add break at first and last thickness value
       else:#if interval boundaries are specified
           for intbound in range(len(intv)):
               for i in range(len(thickness)):
                   j=0
                   if(thickness[i]<intv[intbound]):
                      thickness2[j]=thickness[i]
                      j=j+1
                      indthick[intbound]=i
               intbind[intbound]=max(thickness2)
           #print indthick
           indthick2=zeros(len(intv)+1)
           intbind2=zeros(len(intbind)+1)
           intbind2[0]=1
           i=1
           while(i<len(intbind)+1):
               intbind2[i]=intbind[i-1]
               i=i+1
           j=0
           while(intbind2[j]!=0):
               j=j+1
           intbind3=zeros(int(j))
           for i in range(int(j)):
               intbind3[i]=intbind2[i]
           intbind3[0]=0
           k=1
           while(k<len(indthick)+1):
               indthick2[k]=indthick[k-1]
               k=k+1
#**************************************calculate stats for each interval*****************************************
       st=zeros(len(intbind3)-1)
       last=zeros(len(intbind3)-1)
       sload=zeros(shape=(len(indthick),len(phi)))
       m1=zeros(len(st))
       m3=zeros(len(st))
       m4=zeros(len(st))
       stdev=zeros(len(st))
       fwmedian=zeros(len(st))
       fwsort=zeros(len(st))
       fwmean=zeros(len(st))
       fwskew=zeros(len(st))
       fwkurt=zeros(len(st))
       p=zeros(shape=(len(st),len(phi)))
       s=zeros(shape=(len(st),len(phi)))
       wpc=zeros(shape=(len(st),len(phi)))
       midint=zeros(len(indthick2)-1)
       for i in range(0,len(intbind3)-1):
           st[i]=intbind3[i]
       j=0
       for i in range(1,len(intbind3)):
           last[j]=intbind3[i]
           j=j+1
       #print last,st
       for j in range(0,len(indthick2)-1):
           for k in range(len(phi)): #get cumulative sediment loads for each phi class
               #sload[k][j]=sum(hitload[int(indthick[j]):int(indthick[j+1])])
               for l in range(int(indthick2[j]),int(indthick2[j+1])):
                    if(phi[stime[l][2]]==phi[k]):
                        sload[j][k]=hitload[l]+sload[j][k]
           sedstats1=sedstats(phi,sload[j])
           m1[j]=sedstats1[0]
           m3[j]=sedstats1[1]
           m4[j]=sedstats1[2]
           stdev[j]=sedstats1[3]
           fwmedian[j]=sedstats1[4]
           fwsort[j]=sedstats1[5]
           fwmean[j]=sedstats1[6]
           fwskew[j]=sedstats1[7]
           fwkurt[j]=sedstats1[8]
           p[j,:]=phi
           s[j,:]=sload[j]
           wpc[j,:]=100*sload[j]/(sum(sload[j])) #weight percent in each phi interval
           midint[j]=(intbind3[j]+intbind3[j+1])/2.0
       #print phi,midint,m1,stdev,m3,m4,fwmedian,fwmean,fwsort,fwskew,fwkurt
#=================================================plot figure part======================================
##============================plot figure for Grain Statistics, 1st Moment for layers===================
##This figure is desinged to show the relationship between location and first hit grain size
       fig1=figure()
       fname1='Grain Statistics, 1st Moment for layers'
       ax1=subplot(111)
       ax1.set_xlabel('1st Moment (phi)')
       ax1.set_ylabel('Midpoint of layer (m)')
       ax1.set_title(fname1)
       ax1.plot(m1,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/m1/'+fname1+' Sample%i.png'%ni,dpi=100)
##============================plot figure for Standard Deviation for layers=============================
##This figure is designed to show the relationship between location and standard deviation
       fig2=figure()
       fname2='Grain Statistics, Standard Deviation for layers'
       ax2=subplot(111)
       ax2.set_xlabel('Standard Deviation (phi)')
       ax2.set_ylabel('Midpoint of layer (m)')
       ax2.set_title(fname2)
       ax2.plot(stdev,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/stdev/'+fname2+' Sample%i.png'%ni,dpi=100)
##=================================plot figure for FW mean for layers===================================
##This figure is designed to show the relationship between location and mean
       fig3=figure()
       fname3='Grain Statistics, FW mean for layers'
       ax3=subplot(111)
       ax3.set_xlabel('FW Mean (phi)')
       ax3.set_ylabel('Midpoint of layer (m)')
       ax3.set_title(fname3)
       ax3.plot(fwmedian,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwmedian/'+fname3+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW sorting for layers==================================
##This figure is designed to show the relationship between location and sorting
       fig4=figure()
       fname4='Grain Statistics, FW sorting for layers'
       ax4=subplot(111)
       ax4.set_xlabel('Sorting (phi)')
       ax4.set_ylabel('Midpoint of layer (m)')
       ax4.set_title(fname4)
       ax4.plot(fwsort,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwsorting/'+fname4+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Skewness for layers=================================
##This figure is designed to show the relationship between location and skewness
       fig5=figure()
       fname5='Grain Statistics, FW Skewness for layers'
       ax5=subplot(111)
       ax5.set_xlabel('Skewness')
       ax5.set_ylabel('Midpoint of layer (m)')
       ax5.set_title(fname5)
       ax5.plot(fwskew,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwskew/'+fname5+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Kurtosis for layers=================================
##This figure is designed to show the reltaionship between location and kurtosis
       fig6=figure()
       fname6='Grain Statistics, FW Kurtosis for layers'
       ax6=subplot(111)
       ax6.set_xlabel('Kurtosis')
       ax6.set_ylabel('Midpoint of layer (m)')
       ax6.set_title(fname6)
       ax6.plot(fwskew,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwkurt/'+fname6+' Sample%i.png'%ni,dpi=100)
##==========================plot figure for Thickness vs. Time sediment hits bottom=====================
##This figures is designed to show the relationship between thickness and time sediment hits bottom
       fig7=figure()
       fname7='Thickness vs. Time sediment hits bottom'
       ax7=subplot(111)
       ax7.set_xscale('log',basex=10)
       ax7.set_xlabel('Time(s)')
       ax7.set_ylabel('Thickness (m)')
       ax7.set_title(fname7)
       ax7.plot(stime1,thickness,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/hittime/'+fname7+' Sample%i.png'%ni,dpi=100)

'''
Created on 2012-12-10
=========================The tsunami water depth and average velocity model==========================
#==================================The calculaet model part==========================================
##===================================handle input parameter==========================================
##================================preprocess the data for equation===================================
#=================================The solve equation part============================================
@author:Hui Tang-tanghui@vt.edu
'''
#====================================================================================================
#====================================================================================================
#============The tsunami water depth and average velocity model based on Andrew L. Moore=============
#====================================================================================================
#====================================================================================================
#This model ia bsed on a paper from Sedimentary Geology 200(2007)336_346
#LANDWARD FINING FROM MULTIPY SOURCES IN A SAND SHEET DEPOSITED BY THE 1929 GRAND BANKS TSUNAMI
#NEWFOUNDLAND
#Andrew L. Moore, Brian G. McAdoo,Alan Ruffman
#Input:         Dl-Largest grain's diameter m
#               Dm-medean grain's diameter m
#               rho-water density
#Used constant: tcr1-Dimensionless critical shear stress 0.045-0.06
#               g-gravity acceleration: 9.8m/s^2
#               rho-Density of fluid: 1050 kg/m^3
#               rhos-Density of sediments: 2770 kg/m^3
#               k-Karmon's constant: 0.41
#               R-Rouse number: 0.55<R<0,63 for fine peak; 1.8<R<2.1 for coaser peak
#Output:        Vi-Smallest possible average velocity m/s
#Hui Tang 05/22/2013
#====================================================================================================
#=================================The calculate model part===========================================
#This part is designed for calculatethe critical shear stress,shear velocity,particle velocity and set
#these equations for this model
from pylab import *
import os
def Moore_function(Dl,Dm,rho):
#========================================= Handle input parameters===================================
#This part is used to read parameters for this model frome document
    separator=':'
    with open('parameter_P16.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
	        if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='Dimensionless_shear_velovity': Tstarcr=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Rouse_number': P=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='Largest_distance': l=float(v)
                elif n=='Gravity_acceleration': g=float(v)
            except IndexError as e:
	        continue
	    except ValueError as e:
	        continue
#======================================Preprocess the data for equations===========================
#This function accumulate the critical shear stress,shear velocity and particle settling velocity
#Take the result as these equation's parameters
    def processdata1(tcr1,D,rhos1,rhos2,g):#get the critical shear stress and shear velocity
        if tcr1>0.045:
            if tcr1<0.06:
                tcr2=tcr1*D*g*(rhos1-rhos2)# get the critical shear stress
                U1=sqrt(tcr2/rhos2)# get the shear velocity
            else:
                    print "please make dimensionless critical shear stress less than 0.06"
        else: print "Please make dimensionless critical shear stress larger than0.045"
        return U1
    def processdata2(R,Vc,u1):#get the particle settling velocity
        if R>0.05:
            if R<0.63:
                ws=R*Vc*u1#get the particle settling velocity
            else:
                    print "Please make Rouse number smaller than 0.63"
        else: print "Please make Rouse number larger than 0.05"
        return ws
    u1=processdata1(Tstarcr,Dl,Rhos,rho,100*g)
    Ws=processdata2(P,vk,u1)
    h=arange(1,10000,1) #Depth of water
    z0=Dm/30 #Roughtness of land surface
    u2=l*Ws/h
    u3=u1*(log(h/z0-(1-z0/h)))/P
#=================================Solve the equation system==================================
#This part is designed for solving equation and get the tsunami's average depth and velocity
    def fn(h,Ws,P,u1,l,z0):
        s=l*Ws/h-u1*(log(h/z0-(1-z0/h)))/P
        return(s)
    def solve_function(s,h,delta=0.0001):
        s1=s
        if s1>0:
            while (s1>0):
                h=h+delta;
                s1=fn(h,Ws,P,u1,l,z0);
            H=h;
        else:
            while(s1<0):
                h=h-delta;
                s1=fn(h,Ws,P,u1,l,z0);
            H=h;
        return H
    h=5
    s=fn(h,Ws,P,u1,l,z0)
    H=solve_function(s,h,delta=0.001) #The average depth of water
    U1=u1*(log(H/z0-(1-z0/H)))/P #The average velocity of flow
    U2=l*Ws/H
    return(U2)

'''
Created on 2012-11-24;
===============================================Output data to CSV document==================================
#======================================================progressBar==========================================
#=================================================output2CSV function one===================================
#=================================================output2CSV function two===================================
@author: Hui Tang tanghui@vt.edu
'''
#===========================================================================================================
#===========================================================================================================
#==================================================Output Data to CSV Script================================
#===========================================================================================================
#===========================================================================================================
#This Script is used to output data to CSV document
from pylab import *
import csv
import sys
# ================================================= progressBar ============================================
# CLASS NAME: DLLInterface
# Author: Larry Bates (lbates@syscononline.com)
# Written: 12/09/2002
# Released under: GNU GENERAL PUBLIC LICENSE
class progressBar:
    def __init__(self, finalcount, title='Progress' ,progresschar='>'):
        import sys
        self.finalcount=finalcount; self.blockcount=0
        if not progresschar: self.block=chr(178)
        else:                self.block=progresschar
        # Get pointer to sys.stdout so I can use the write/flush
        self.f=sys.stdout
        # If the final count is zero, don't start the progress gauge
        if not self.finalcount : return
        self.f.write('\n'+title+'\n[')
        return
    def progress(self, count):
        # Make sure I don't try to go off the end (e.g. >100%)
        count=min(count, self.finalcount)
        # If finalcount is zero, I'm done
        if self.finalcount:
            percentcomplete=int(round(100*count/self.finalcount))
            if percentcomplete < 1: percentcomplete=1
        else:            percentcomplete=100
        #print "percentcomplete=",percentcomplete
        blockcount=int(percentcomplete/2)
        #print "blockcount=",blockcount
        if blockcount > self.blockcount:
            for i in range(self.blockcount,blockcount):
                self.f.write(self.block)
                self.f.flush()
        if percentcomplete >= 100: self.f.write("]")
        self.blockcount=blockcount
        return
# ============================================== output2CSV =============================================
# this function output given data to csv file
def output2CSV(name,data1,data2,text1,text2):
    data=zeros([len(data1),2]);
    data[:,0]=data1; data[:,1]=data2;
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(outp_path+name,'wt') as f:
        writer=csv.writer(f,delimiter=';')
        writer.writerow((text1,text2)) # write description
        for row in data: writer.writerow((row[0],row[1]));pb.progress(n+1);n+=1;
    return data
# ============================================= output2CSV1 =============================================
# this function output given data to csv file
def output2CSV1(name,data1,data2,data3,data4,data5,text1,text2,text3,text4,text5):
    data=zeros([len(data1),5]);
    data[:,0]=data1; data[:,1]=data2;data[:,2]=data3;data[:,3]=data4;data[:,4]=data5
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(outp_path+name,'wt') as f:
        writer=csv.writer(f,delimiter=';')
        writer.writerow((text1,text2,text3,text4,text5)) # write description
        for row in data: writer.writerow((row[0],row[1],row[2],row[3],row[4]));pb.progress(n+1);n+=1;
    return data
# ============================================= output2CSV2 =============================================
# this function output given data to csv file
def output2CSV2(name,data1,data2,data3,text1,text2,text3):
    data=zeros([len(data1),3]);
    data[:,0]=data1; data[:,1]=data2;data[:,2]=data3;
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(name,'wt') as f:
        writer=csv.writer(f,delimiter=';')
        writer.writerow((text1,text2,text3)) # write description
        for row in data: writer.writerow((row[0],row[1],row[2]));pb.progress(n+1);n+=1;
    return data

'''
Created on 2013-5-26;
=====================Joint Tsunami Sediment Speed model pre-processing Script=====================
#==================================Grainsize Distribution Part====================================
#==========================================Water Density Part=====================================
#==========================================Water Depth Part=======================================
@author: Hui Tang, tanghui@vt.edu
'''
#=================================================================================================
#=================================================================================================
#====================Joint Tsunami Sediment Speed Model Pre-processing Script=====================
#=================================================================================================
#=================================================================================================
#This function is designed to preprocess the field data to input the joint model.
#Inputs:            filename-The name of field data document, must be csv documnet contain grainsize
#                            and percentage of each grainsize,separated by ";" and in a row.
#Function used:     sedstat.py,sw_dens0.py
#Output:            Dl: Largest grainsize in phi
#                   Ds: Smallest grainsize in phi
#                   Dm: Medean grainsize in phi
#                   Nc: The number of grainsize class
#                   Rho:Density of sea water
#                   H:Max water depth
from pylab import *
# from function import *
# from sedstats import *
# from ReadCSV import *
# from output2CSV import *
import csv
import os
def preproc(filename):
#===================================================================================================
#===========================================Input Parameter=========================================
#===================================================================================================
#This part is designed for inputing all parameter for all model and major function
    separator=':'
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
	        if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='water_temperature': wtemp=float(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_run_up':Rz=float(v)
                elif n=='slope': m=float(v)
                elif n=='Largest_distance': l=float(v)
                elif n=='Depth': h=float(v)
                elif n=='Sediment_thickness':th=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
#=========================================Grainsize Distribution Part==============================
    filename1='../input/data_process.csv'
    data=readCSV(filename,separator=';')
    phi=data[:,0]
    fr=data[:,1]
    Dl=0
    Ds=data[-1,0]
    Nc=len(data[:,0])
    fr=fr.reshape(len(fr),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    if(len(fr)==1):
            Dm=phi[0]
    else:
            sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
            Dm=sedstat[0]
#==========================================Water Density Part=====================================
    fr=fr.reshape(len(fr))
    rho=sw_dens0(sal,wtemp)
#==========================================Water Depth Part=======================================
    H=h*m*Rz/(m*Rz-l)
    fp=output2CSV1(filename1,phi,100*fr,th,h,l,'phi','fraction','thickness','depth','location')
    return(Dl,Ds,Dm,Nc,rho,H)


# from pylab import *
# from function import *
# from ReadCSV import *
 # from output2CSV import *
# from tools2 import *
# import os
# import csv
#
# filename='data_draw.csv'
# filename1='data_draw1.csv'
# Data=readCSV2(filename,separator=';')
# se=Data[:,0]
# data_x=Data[:,1]
# data_y=Data[:,2]
# for j in range(2):
#     n=0
#     for i in range(len(se)-1):
#         if(se[i]==nan):
#             if(se[i+1]!=nan):
#                 if(i!=len(se)-1):
#                     se[i]=(se[i-1]+se[i+1])/2
#                 else:
#                     se[i]=se[i-1]
#             else:
#                 se[i]=(se[i-1]+se[i+2])/2
#         if(se[i]!=0):
#             n=n+1
# data1_se=zeros(n)
# data1_x=zeros(n)
# data1_y=zeros(n)
# for i in range(n):
#     data1_se[i]=se[i]
#     data1_x[i]=data_x[i]
#     data1_y[i]=data_y[i]
#             #data1_se=append(data1_se,se[i])
#             #data1_x=append(data1_x,data_x[i])
#             #data1_y=append(data1_y,data_y[i])
# fp=output2CSV2(filename1,data1_se,data1_x,data1_y,'Mean_grainsize','x_location','y_location')

'''
Created on 2013-5-24;
===============================================Read data From CSV document==================================
#=================================================readCSV function one======================================
#=================================================readCSV function two======================================
@author: Hui Tang tanghui@vt.edu
'''
#============================================================================================================
#============================================================================================================
#===========================================Read Data From CSV document Script===============================
#============================================================================================================
#============================================================================================================
#This Script is used to read data from csv document
#Writen by Hui Tang 05/24/2013
#==================================================readCSV one===============================================
#This is a function to read data from CSV document,two items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction)
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu
def readCSV(name,separator=';'):
    data=zeros([1,2]); start=True # initialize data array
    with open(name,'rU') as f:
        reader = csv.reader(f) # read CSV file
        for row in reader: # read each line and store them into a n by 2 array
            try: rs=row[0].split(separator); d1=float(rs[0]); d2=float(rs[1]);
            except IndexError: continue # if line is not delimited by separator, continue
            except ValueError: continue # if not numbers, continue
            if (d1!="\x00"):
                if not start: data=append(data,[[d1,d2]],axis=0);
                else: start=False; data[0]=[d1,d2]
    return data
#===============================================readCSV two===================================================
#This is a function to read data from CSV document,five items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction),sediment thickness, water
#                        depth and location
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu
def readCSV1(name,separator=';'):
    data=zeros([1,5]); start=True # initialize data array
    with open(outp_path+name,'rU') as f:
        reader = csv.reader(f) # read CSV file
        for row in reader: # read each line and store them into a n by 2 array
            try: rs=row[0].split(separator); d1=float(rs[0]); d2=float(rs[1]);d3=float(rs[2]);d4=float(rs[3]);d5=float(rs[4])
            except IndexError: continue # if line is not delimited by separator, continue
            except ValueError: continue # if not numbers, continue
            if (d1!="\x00"):
                if not start: data=append(data,[[d1,d2,d3,d4,d5]],axis=0);
                else: start=False; data[0]=[d1,d2,d3,d4,d5]
    return data
#===============================================readCSV three===================================================
#This is a function to read data from CSV document,five items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction),sediment thickness, water
#                        depth and location
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu




'''
Created on 2013-5-24
===========================================The Root Error Caculate Part==================================
#===============================================Input Data===============================================
#==========================================The calculate part============================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=====================================The Root Error Calculate part model================================
#========================================================================================================
#========================================================================================================
#This part is designed to calculate the root square error between field data and model result, if the RSE
# is larger than 10, change the grainsize distribution of sediment source, if less than 10, stop go to nextt
# step directly.
#This part is based on the paper published on Sedimentary Geology 282(2012) 90-109
# Flow Speed estimated by inverse modeling of sandy tsunami deposits: results from the 11 March tsunami on
# the coastal plain near the Sendai Airport, Honshu, Japan
#               Bruce E. Jaffe and Kazuhisa Goto
#Hui Tang 05/24/2013
#Input:		phi-grain size
#           Fr2-Model result Distribution
#           Num-The number of Fraction of sub-intervals in each sample location
#Output:	RSE-Average Root Sqaure Root
#=======================================================================================================
#==========================================Input Data===================================================
#This part is designed to input data to calculate RSE

from pylab import *
import os
# from function import *
# from ReadCSV import *
def calrse(name,phi,Fr2,Num):
    sumerr=zeros(Num)#RSE for each sub-intervals
    for i in range(Num):
        data1=readCSV(name,separator=';')
    	phi=data1[:,0]
        Fr1=data1[:,1]
#=============================================The calculate part========================================
        n=len(Fr1)
        m=len(Fr2)
        err=zeros(n)
        sumerr=zeros(n)
        if(m!=n):
            print('These Dimension of two document are not same, Please check the data structure')
        else:
            for j in range(n):
                err[j]=((Fr1[j]-Fr2[j]))**2
                sumerr[i]=sum(err)
                sumerr[i]=sqrt(sumerr[i]/n)
    RSE=sum(sumerr)
    return(RSE)


from pylab import *
#=============================================function sedstats===========================================
#This function is designed to compute moment measures and Folk and Ward statistics from sediment grain size
#weight % and phi size classes.
#Inputs: phi:grain size (phi) for each size class
#        weight:weight (or weight percent, vol concentration) for size classes
#Calculated: phimdpt:phi midpoints used for moment measures for tsunami samples, [-1.125:.25:10 12]
#Outputs:    m1:mean grain size (phi), 1st moment
#            stdev:standard deviation (phi) similar to sorting, is equal to square root of 2nd moment
#            m3:third moment
#            m4:fourth moment
#            fwmedian:Folk and Ward medium grain size (phi), 50th percentile
#            fwmean:Folk and Ward mean grain size (phi)
#            fwsort:Folk and Ward sorting (phi)
#            fwskew:Folk and Ward skewness
#            fwkurt:Folk and Ward Kurtosis
def sedstats(phi,weight):
    if(len(phi)==1):
        print('Moment measures not valid because there is only one size class')
        #phi=
        #create a cumulative weight with 0% weight in it to allow calculation
    else:
        mpt1=phi[0]+0.5*(phi[0]-phi[1])
        phimdpt=zeros(len(phi))
        phimdpt[0]=mpt1
        phimdpt[1:len(phi)]=phi[1:len(phi)]+0.5*(phi[0:len(phi)-1]-phi[1:len(phi)])
        phimdpt=transpose(phimdpt)
        #p=phimdpt*weight
        #print p
        #print phimdpt
        #print weight
        #calculate the 1st moment (mean grain size)
        m1=sum(phimdpt*weight)/sum(weight)
        dev=phimdpt-m1
        #calculate the standard deviation, similar to sorting, square root of the 2nd moment
        var=sum(weight*(dev**2))/sum(weight)
        stdev=sqrt(var)
        #calculate the third moment
        m3=sum(weight*((dev/stdev)**3))/sum(weight)
        #calculate the fourth moment
        m4=sum(weight*((dev/stdev)**4))/sum(weight)
        #Use phi intervals, not midpoints to calculate Folk and Ward stats- for tsunami samples use
        cw=100*cumsum(weight)/sum(weight)#calculate normalized cumulative weight percent
        cp=[5,16,25,50,75,84,95]
        if(cw[0]>=5):
            print('Folk and Ward Statistics suspect because first size class has >=5% cumulative %')
            cw1=zeros(len(cw)+1)
            cw1[1:len(cw)+1]=cw[:] #create a cumulative weight with 0% weight in it to allow calculation
        #calculate cumulative percents used in Folk and Ward statistics using
        #linear interpolation (could use different interpolation scheme)
        cumphi=zeros(7)
        cw1=zeros(int(len(phi)))
        for i in range(7):
            k=0
            while(cw[k]<=cp[i]):
                cw1[k]=cw[k]
                k=k+1
            lp=argmax(cw1)
            slp=(cp[i]-cw[lp])/(cw[lp+1]-cw[lp])
            cumphi[i]=phi[lp]+slp*(phi[lp+1]-phi[lp])
        #make some names that have meaning
        phi5=cumphi[0]
        phi16=cumphi[1]
        phi25=cumphi[2]
        phi50=cumphi[3]
        phi75=cumphi[4]
        phi84=cumphi[5]
        phi95=cumphi[6]
        #calculate Folk and Ward stats
        fwmedian=phi50
        fwmean=(phi16+phi50+phi84)/3
        fwsort=(phi84-phi16)/4.0+(phi95-phi5)/6.6
        fwskew=(phi16+phi84-2*phi50)/(2*(phi50-phi16))+(phi5+phi95-2*phi50)/(2*(phi95-phi5))
        fwkurt=(phi95-phi5)/(2.44*(phi75-phi25))
        wp=weight*phi
        #print len(weight[0,:])
        #percen=zeros(len(weight[0,:]))
        #print weight[0,39]
        #for j in range(len(weight[0,:])):
            #percen[j]=weight[0,j]
        #percen=percen.reshape(len(percen[:,0]))
        #print percen
        #sumw=sum(weight)
        #wp=zeros(len(phi))
        #percen=zeros(len(phi))
        #for j in range(len(phi)):
            #wp[j]=percen[j]*phi[j]
        fwmean=sum(wp)/100
        sedstat=[m1,m3,m4,stdev,fwmedian,fwsort,fwmean,fwskew,fwkurt]
        return(sedstat)


'''
Created on 2013-5-21
=====================================The source grain size distribution===================================
=====================================The calculate model part=============================================
#=============================================Input parameters============================================
#=====================================Coculate the pecerntage of each grain size==========================
#=====================plot original grain size distibution================================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=======The original tsunami sediment grain size distibution model based normalize distribution==========
#========================================================================================================
#========================================================================================================
#This model is based on idealize grainsize distribution which will be used for an input for joint model,
#it will be set as normalize distribution and can be other dstribution based on the problem.
#Input:             Dl-Largest grain diameter in phi
#                   Dm-Medean grain diameter in phi
#                   Ds-Smallest grain diameter in phi
#                   nc-Number of grain size class
#Output:            Se-Grain size in phi
#                   fr-Percentage of each grain size %
# Hui Tang 2013-5-21
#========================================================================================================
#=============================================Input parameters===========================================
#========================================================================================================
from pylab import *
import os
# from output2CSV import *
# from function import *
# from ReadCSV import *
def source_distribution(Dl,Ds,Dm,nc):
#========================================================================================================
#=====================================Coculate the pecerntage of each grain size=========================
#========================================================================================================
#Conculating the pencertage for different grain size and set different grian size for plot figure in the
#phi scale
    se=linspace(Dl,Ds,nc)
    #Dm=2.8
    filename1='../input/sample_P14abc.csv'
    Field_data=readCSV(filename1,separator=';')
    Field_se=Field_data[:,0]
    fr=exp(-(Field_se-Dm)**2/2*1)/sqrt(2*(3.14))
    fr=fr/sum(fr)
#========================================================================================================
#=====================plot original grain size distibution===============================================
#========================================================================================================
    #fname1='original grain size distibution.csv'
    #fp2=output2CSV(fname1,Field_se,fr,'phi','fraction')
    #f=figure()
    #fname='original grain size distibution'
    #ax=subplot(111)
    #ax.set_xlim(0,7)
    #ax.set_ylim(0,6)
    #ax.set_title(fname)
    #ax.plot(Field_se,100*fr,'-')
    #ax.set_xlabel('Grain size')
    #ax.set_ylabel('Pecentage (%)')
    #grid()
    #show()
    return(Field_se,fr)
def source_distribution1(Dl,Ds,Dm,nc):
#========================================================================================================
#=====================================Coculate the pecerntage of each grain size=========================
#========================================================================================================
#Conculating the pencertage for different grain size and set different grian size for plot figure in the
#phi scale
    se=linspace(Dl,Ds,nc)
    #Dm=2.8
    filename1='../input/sample_P14abc.csv'
    Field_data=readCSV(filename1,separator=';')
    Field_se=Field_data[:,0]
    fr=exp(-(Field_se-Dm)**2/2*1)/sqrt(2*(3.14))
    fr=fr/sum(fr)
#========================================================================================================
#=====================plot original grain size distibution===============================================
#========================================================================================================
    #fname1='original grain size distibution.csv'
    #fp2=output2CSV(fname1,Field_se,fr,'phi','fraction')
    #f=figure()
    #fname='original grain size distibution'
    #ax=subplot(111)
    #ax.set_xlim(0,7)
    #ax.set_ylim(0,6)
    #ax.set_title(fname)
    #ax.plot(Field_se,100*fr,'-')
    #ax.set_xlabel('Grain size')
    #ax.set_ylabel('Pecentage (%)')
    #grid()
    #show()
    return(se,fr)


'''
Created on 2013-1-13
======================================The tsunami sediment model=========================================
#=====================================The calculate model part===========================================
#=============================================Input parameters===========================================
##==============================Conculate the Deposit Velosity===========================================
##================Conculate some different parameters for tsunami sediment model ========================
##=====================================Calculate the thickness at x=0====================================
##==================Calculate the horizontal and vertical sediment run-up along the slope================
##==Calculate different grain size's thickness,the number of grains and area everywhere along the slope==
##=============================calculate percentage of different grian size==============================
#=========================================The plot figure part===========================================
##=================================plot the location of sediments end====================================
##===================plot the different grain size's thickness on the still water level==================
##========================plot different grains thickness on different location==========================
##===================================plot pecentation of different grian size============================
##===================================plot D50 vs different location======================================
##=======================plot thickness vs different location in different grain size====================
##=======================plot total thickness vs different location in different grain size==============
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#===========================The tsunami sediment model based on R L. Soulsby=============================
#========================================================================================================
#========================================================================================================
#This model is based on the paper published in sixth International Symposium on Coastal Sediment Processes
#RECONSTRUCTING TSUNAMI RUN-UP FROM SEDIMENTARY CHARACTERISTICS- A SIMPLE MATHMATICAL MODEL
#             Richard L. Soulsby, David E. Smith, Alan Ruffman
#Input:             V-Average velocity m/s
#                   C0-Initial sediment concentration kg/m^3
#                   se-Sediment grainsize in phi
#                   H-Water depth at still water level
#                   nc-Number of grain size classes
#Constant used:     Rz-Vertical water run-up m
#                   m-slope
#                   Rhos-Sediments density
#                   wtemp-Water tempertature C
#                   sal-Salinity of sea water (psu)
#                   q-sediment deposit rate
#                   nc-Number of grainsize classes
#                   N-Number of sample location
#Fuction used:      Tubesetlevel
#                   inconcerntration
#Output:            th2-Sediment load for each sample location
#                   se1-grain size in phi
#                   fr-Grain size distribution for each sample point
#========================================================================================================
#===========================================The calculate model part=====================================
#This part is designed for calculate the sediment grains deposit velosity, horizontal run-up limits of water,
#horizontal and vertical run-up limits of tsunami sediments, sediment thickness, the number of sediment grains
#and percentage of different grain size in different location along the slope
#========================================================================================================
#=============================================Input parameters===========================================
#========================================================================================================
#This part is designed to input the parameter for this model
from pylab import *
import os
# from output2CSV import *
# from function import *
# from ReadCSV import *
def sousby(se,C0,V,H,nc):
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                if n=='sediment_deposit_rate': q=float(v)
                elif n=='water_run_up':Rz=float(v)
                elif n=='slope': m=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Number_sample': N=int(v)
                elif n=='Filename': filename=str(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='salinity': sal=float(v)
                elif n=='pi': pi=float(v)
                elif n=='Filename_Sousby': filename=str(v)
                elif n=='settling_velocity_type': setType=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
			continue
    C=sum(C0)
    #print C
    se2=phi2mm(se)*0.001
#===================================Calculate the Deposit Velosity=======================================
#Calculating the deposit velosity for different grain size and set different grian size for plot figure in
#the phi scale
    rhow=sw_dens0(sal,wtemp)
    if(setType==1.0):
        ws=tubesetvel(se,Rhos,wtemp,sal)#Settling velocity of each grain size
    elif(setType==2.0):
        ws=100*dietrichWs(wtemp,se2,Rhos,rhow,C/2650.0,0.7,3.5)
    elif(setType==3.0):
        ws[i]=Setvel1(se2,kinvisc,s[i])
    #ws=tubesetvel(se,Rhos,wtemp,sal)#Settling velocity of each grain size
#====================Conculate some different parameters for tsunami sediment model =====================
    Rw=Rz*m #The horizontal run-up limit of tsunami water
    T=Rw/V#run up time
    a=ws*T/H #constant a in model
    Rs=Rw/(1+a*q) #The horizontal run-up limit of tsunami sediment
    #print "Rs"
    #print Rs
    x=linspace(0,350,int(N))
#======================================Calculate the sediment load at x=0=================================
#This part is designed for calculate different grain size's thickness at still water level
    th=zeros(nc)# The sediment load of different grain size sediments in still water level
    th=a*(1+a*q)*C0*H/((1+a)*Rhos)
#=============================Calculate different grain size's sediment load along the slope==================
#This part is designed to calculate sediment load, the area of the sample location, the area of sediment
#will be used to calculate the percentage of different grain size
    i=0 #index for grain size
    j=0 #index for sample location
    th2=zeros(shape=(N,nc)) #The sediment thickness in different sample location
    A=zeros(shape=(N,nc)) #The area of sediment grain in different sample location
    sumA=zeros(N) # The total area of sediment grain in different sample location
    totalth2=zeros(N)# The total sediment thickness in different sample location
    Dt=0.01 #The width of area
    while(j<N): #calculate the thickness, number and area
        i=0
        while(i<nc):
            th2[j,i]=th[i]*(1-x[j]/Rs[i])
            if(th2[j,i]>0):
                A[j,i]=th2[j,i]*Dt
                sumA[j]=sumA[j]+A[j,i]
                totalth2[j]=sumA[j]/Dt
            i=i+1
        j=j+1
#================================calculate percentage of different grian size=================================
#This part is design to use the area to calculate the percentage of different grain size in different sample
#location
    #print th2[9]
    p=0 #index for sample location
    pecentA=zeros(shape=(N,nc)) #The percentage of different grain size
    med=zeros(N)# The D50 of sample location
    while(p<N): # calculate the percentage
        n=0
        sumpecent=0
        while(n<nc):
            pecentA[p,n]=100*A[p,n]/sumA[p]
            sumpecent=sumpecent+pecentA[p,n]
            if(sumpecent<=55):
                med[p]=se[n]
            n=n+1
        p=p+1
#===========================================================================================================
#==================================The plot figure part=====================================================
#This part is designed to plot figure for this model include:total sediment load, mean grain size,flow depth
##==================================Output Data and parameters==========================================
#This part is designed to output and plot Data and parameter
    h=x*(Rz-H)/(Rz*m)-x/m+H# Water depth
    for i in range(N):
        pe=pecentA[i,:]
        filename2="sample%02d"%i+".csv"
        fp2=output2CSV1(filename2,se,pecentA[i,:],totalth2[i],h[i],x[i],'phi','fraction','thickness','depth','location')
        f4=figure()
        filename="Grainsize distribution at different location x=%f"%x[i]+"m"
        ax4=subplot(111)
        ax4.plot(se,pecentA[i,:],'-')
        ax4.set_xlim(0,10)
        ax4.set_ylim(0,40)
        ax4.set_xlabel('grain size($\phi$)')
        ax4.set_ylabel('The percentage of different sediment grains size(%)')
        #savefig(filename+'.png',dpi=100)
        close()
    return(totalth2)

'''
Created on 2013-5-24
========================The Root Error Caculate Part for sediments thickness=============================
#===============================================Input Data===============================================
#==========================================The calculate part============================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=====================================The Root Error Calculate part model================================
#========================================================================================================
#========================================================================================================
#This part is designed to caculate the sediment difference between field data and model result, if the
# difference between two data is larger than 5%, change the average velocity, if less than 5%, stop go to
# next step directly.
#Hui Tang 05/24/2013
#Input:         name-Document name of field data document
#               th2-Model result of suspended sediemnts thickness m
#               N-The number of sample location
#Output:        averthRSE-Average thickness difference between field data and model result
#=======================================================================================================
#==========================================Input Data===================================================
#This part is designed to input data to calculate thickness difference
from pylab import *
import os
# from ReadCSV import *
# from function import *
def thcalrse(name,th2,N):
    rse=zeros(N)#Average thickness difference between field data and model result for each sample location
    data1=readCSV(name,separator=';')
    th1=data1[:,1]
    k=0
    n=len(th1)
    m=len(th2)
    if(n==m):
        for j in range(N):
            if (th1[j]!=0):
                if(th2[j]!=0):
                    rse[j]=(th1[j]-th2[j])*100/th1[j]
                    k=k+1
    else:
        print("These Dimension of two document are not same, Please check the data structure")
    averthRSE=sum(rse)/N
    #print th1
    #print th2
    #print rse
    return(averthRSE)

'''
Created on 2013-5-26
==========================The TsuSedForm Tsunami Sediment formation Model==========================
#========================================Calculate Part============================================
##=========================Calculate SL and z for sediment Formation===============================
##================================Calculate Settling velocity======================================
##=================================Calculate Settling Time=========================================
##=========================Calculate hittingsize and concertration=================================
##=================================Determine Boundary==============================================
##=================================Calculate grainsize distribution================================
##==============================Calculate stats for each interval==================================
#========================================Plot figure part==========================================
##========================Plot grainsize distribution for each interval============================
##===================plot figure for Grain Statistics, 1st Moment for layers=======================
##========================plot figure for Standard Deviation for layers============================
##============================plot figure for FW mean for layers===================================
##==========================plot figure for FW sorting for layers==================================
##==========================plot figure for FW Skewness for layers=================================
##==========================plot figure for FW Kurtosis for layers=================================
##====================plot figure for Thickness vs. Time sediment hits bottom======================
#=========================================OUTPUT Data Part=========================================
@ author Hui Tang- tanghui@vt.edu
'''
#==================================================================================================
#==================================================================================================
#============================The TsuSedForm Tsunami Sediments Formation Model======================
#==================================================================================================
#==================================================================================================
#This model is based on Jaff's TsuSedMod model's grading and static psrts, grading and static parts
#orignate from Jaffe 9/29/2003, Made into python by Hui Tang, Virginia Tech 1/25/2013
#This function is desinged  to form suspended sediments grainsize distribution based on whole grainsize
#distribution in tsunami sediments.
#Input:         name1-Filename of data document from Sousby model
#               h-Water depth m
#               th1-Sediment load m
#               intv-sub-interval bounadry
#               x-Sample location
#Function Used: Tsusemod.py,tubesetvel.py and sedstatics.py
#Output:        se-Suspened sediments grainsize in phi
#               Fr-Suspended seidments percentage for each grain size
from pylab import *
from numpy import *
# from ReadCSV import *
# from function import *
# from sedstats import *
# from TsusedMod import *
# from output2CSV import *
# from ReadCSV import *
def Tsusedform(name1,intv,h,x):
    with open('../input/parameter_P14abc.txt','r') as f:
            for line in f:
                try:
                    s=line.split('=');n=s[0];v=s[1];#remove '\n'
                    if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                    if n=='porosity': porosity=float(v)
                    elif n=='sediment_density': Rhos=float(v)
                    elif n=='Number_sample': N=int(v)
                    elif n=='Filename': filename=str(v)
                    elif n=='water_temperature': wtemp=float(v)
                    elif n=='salinity': sal=float(v)
                    elif n=='pi': pi=float(v)
                    elif n=='Filename_Sousby': filename=str(v)
                except IndexError as e:
                    continue
                except ValueError as e:
			        continue
#========================================Calculate Part============================================
##=========================Calculate SL and z for sediment Formation===============================
    data=Tsusedmod(name1)
    sl=data[0]# Sediment concentration profile
    z=data[1]# Elevation
    phi=data[2]# Grain size
    numrow_sl=len(sl)
    numcol_sl=len(sl[0]) #get matrix sl size
    if((len(z)!=numrow_sl) or (len(phi)!=numcol_sl)):
        print('Number of size classes or number of elevations do not match profiles')
    else:
        z=z.reshape(len(z),1) #z must be a column ector
##================================Calculate Settling velocity======================================
        setvel=tubesetvel(phi,Rhos,wtemp,sal)#calculate settling velocity (m/s)
        setvel=setvel.reshape(1,len(setvel))
        numrow_setvel=len(setvel)
        numcol_setvel=len(setvel[0])
        if(numrow_setvel>numcol_setvel):
            setvel=setvel #make sure setvel is a row vector
##=================================Calculate Settling Time=========================================
        settlingtime=zeros(shape=(numrow_sl,numcol_setvel))#Settling time for each grain size
        stime=zeros(shape=(numrow_sl*numcol_setvel,3))#store for sorting based on settling time
        for i in range(numrow_sl):
            for j in range(numcol_setvel):
                settlingtime[i][j]=z[i]/setvel[0][j]
##=========================Calculate hittingsize and concertration=================================
###==========================Sort by time it takes to hit bottom===================================
        stime[:,0]=settlingtime.reshape(int(numrow_sl*numcol_setvel),)
        ind=argsort(stime[:,0],axis=0)
        k=0
        for i in range(len(ind)):#keep track of grain size and elevation
            stime[k][1]=int(int(ind[i])/int(numcol_setvel))#grain size
            stime[k][2]=int(int(ind[i])%int(numcol_setvel))#location
            k=k+1
        stime1=sorted(stime[:,0])
###===========================find size of sediment hitting the bed================================
        phiind=stime[:,2] # find size of sediment
        hitsize=zeros(len(phiind))#initiate array with size of particles hitting bed
        for i in range(len(phiind)):
            hitsize[i]=phi[phiind[i]]    # keep track of order sediment size hitting the bed
###=========================find elevations that sediment started from=============================
        zind=stime[:,1]# find elevation of sediment
        hitheight=zeros(len(zind)) # initiate array with height of particles hitting bed
        for i in range(len(zind)):
            hitheight[i]=z[zind[i]] #keep track of order of height sediment came from
        hitload=zeros(len(zind)) # initial array with concentration of particles hitting bed
        for i in range(len(zind)):
            hitload[i]=sl[zind[i],phiind[i]]
##=========================Calculate sediment thickness for each hitting grain size=================
        thickness=cumsum(hitload/(1-porosity))#cumulative thickness of deposit
        thickness2=zeros(len(thickness))
##=================================Determine Boundary===============================================
        intbind=zeros(len(thickness))
        indthick=zeros(len(intv))
        if(len(intv)==1):#calculate interval boundaries for set intervals
            for intbound in range(intv,thickness[len(thickness)]-intv/2,intv):
                for i in range(len(thickness)):
                    if(thickness[i]<intv[intbound]):
                        thickness2[j]=thickness[i]
                        j=j+1
                intbind[round(intbound/intv)]=max(thickness2)
            intbind2=zeros(len(intbind)+2)
            intbind2[0]=0
            i=1
            while(i<len(intbind)+1):
                intbind2[i]=intbind[i-1]
                i=i+1
            intbind[i]=len(thickness)#add break at first and last thickness value
        else:#if interval boundaries are specified
            for intbound in range(len(intv)):
                for i in range(len(thickness)):
                    j=0
                    if(thickness[i]<intv[intbound]):
                        thickness2[j]=thickness[i]
                        j=j+1
                        indthick[intbound]=i
                intbind[intbound]=max(thickness2)
            indthick2=zeros(len(intv)+1)
            intbind2=zeros(len(intbind)+1)
            intbind2[0]=1
            i=1
            while(i<len(intbind)+1):
                intbind2[i]=intbind[i-1]
                i=i+1
            j=0
            while(intbind2[j]!=0):
                j=j+1
            intbind3=zeros(int(j))
            for i in range(int(j)):
                intbind3[i]=intbind2[i]
            intbind3[0]=0
            k=1
            while(k<len(indthick)+1):
                indthick2[k]=indthick[k-1]
                k=k+1
##=================================Calculate grainsize distribution================================
        st=zeros(len(intbind3)-1)# start point
        last=zeros(len(intbind3)-1) # end point
        sload=zeros(shape=(len(indthick),len(phi)))# sediment load
        m1=zeros(len(st))# mean grain size
        m3=zeros(len(st))
        m4=zeros(len(st))
        stdev=zeros(len(st))# Standard deviation
        fwmedian=zeros(len(st))# Median grain size 50%
        fwsort=zeros(len(st))# sorting phi
        fwmean=zeros(len(st))# Mean grain size
        fwskew=zeros(len(st))# Skewness
        fwkurt=zeros(len(st))# kurtosis
        p=zeros(shape=(len(st),len(phi)))# grain size
        s=zeros(shape=(len(st),len(phi))) # sediment load
        wpc=zeros(shape=(len(st),len(phi)))
        midint=zeros(len(indthick2)-1)
        for i in range(0,len(intbind3)-1):
            st[i]=intbind3[i]
        j=0
        for i in range(1,len(intbind3)):
            last[j]=intbind3[i]
            j=j+1
        thickness1=zeros(shape=(len(indthick2)-1,len(phi)))
        Fr=zeros(shape=(len(indthick2)-1,len(phi)))#Pecentage for each grain size
        totalth=zeros(len(indthick2)-1)
        for j in range(0,len(indthick2)-1):
            for k in range(len(phi)): #get cumulative sediment loads for each phi class
                #sload[k][j]=sum(hitload[int(indthick[j]):int(indthick[j+1])])
                for l in range(int(indthick2[j]),int(indthick2[j+1])):
                    if(phi[stime[l][2]]==phi[k]):
                        sload[j][k]=hitload[l]+sload[j][k]
                        thickness1[j][k]=sload[j][k]/(1-porosity)
##==============================Calculate stats for each interval==================================
            sedstats1=sedstats(phi,sload[j])
            m1[j]=sedstats1[0]
            m3[j]=sedstats1[1]
            m4[j]=sedstats1[2]
            stdev[j]=sedstats1[3]
            fwmedian[j]=sedstats1[4]
            fwsort[j]=sedstats1[5]
            fwmean[j]=sedstats1[6]
            fwskew[j]=sedstats1[7]
            fwkurt[j]=sedstats1[8]
            p[j,:]=phi
            s[j,:]=sload[j]
            wpc[j,:]=100*sload[j]/(sum(sload[j])) #weight percent in each phi interval
            midint[j]=(intbind3[j]+intbind3[j+1])/2.0
        for j in range(0,len(indthick2)-1):
            for k in range(len(phi)): #get cumulative sediment loads for each phi class
                Fr[j][k]=thickness1[j][k]*100/sum(thickness1[j])
                totalth[j]=sum(thickness1[j])
#========================================Plot figure part==========================================
##========================Plot grainsize distribution for each interval============================
        #This part is designed for plot grainszie distribution and ouput data for each interval
        for i in range(0,int(len(intv))):
            pe=Fr[i,:]
            name=name1.split(".")
            filename=name[0]+"_suspended_sample%02d"%i+".csv"
            fp=output2CSV1(filename,phi,pe,totalth[i],h,x,'phi','fraction','thickness','depth','location')
            f1=figure()
            ax1=subplot(111)
            ax1.plot(phi,pe,'-')
            ax1.set_xlim(0,10)
            ax1.set_ylim(0,40)
            ax1.set_xlabel('grain size($\phi$)')
            ax1.set_ylabel('The percentage of different sediment grains size(%)')
            #savefig(filename+'.png',dpi=100)
            close()
##============================plot figure for Grain Statistics, 1st Moment for layers================
        ##This figure is desinged to show the relationship between location and first hit grain size
        fig1=figure()
        fname1='Grain Statistics, 1st Moment for layers'
        ax1=subplot(111)
        ax1.set_xlabel('1st Moment (phi)')
        ax1.set_ylabel('Midpoint of layer (m)')
        ax1.set_title(fname1)
        ax1.plot(m1,midint,'-')
        close()
        #savefig(fname1+' Sample%i.png'%ni,dpi=100)
##============================plot figure for Standard Deviation for layers=============================
        ##This figure is designed to show the relationship between location and standard deviation
        fig2=figure()
        fname2='Grain Statistics, Standard Deviation for layers'
        ax2=subplot(111)
        ax2.set_xlabel('Standard Deviation (phi)')
        ax2.set_ylabel('Midpoint of layer (m)')
        ax2.set_title(fname2)
        ax2.plot(stdev,midint,'-')
        close()
        #savefig(fname2+' Sample%i.png'%ni,dpi=100)
##=================================plot figure for FW mean for layers===================================
        ##This figure is designed to show the relationship between location and mean
        fig3=figure()
        fname3='Grain Statistics, FW mean for layers'
        ax3=subplot(111)
        ax3.set_xlabel('FW Mean (phi)')
        ax3.set_ylabel('Midpoint of layer (m)')
        ax3.set_title(fname3)
        ax3.plot(fwmedian,midint,'-')
        close()
        #savefig(fname3+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW sorting for layers==================================
        ##This figure is designed to show the relationship between location and sorting
        fig4=figure()
        fname4='Grain Statistics, FW sorting for layers'
        ax4=subplot(111)
        ax4.set_xlabel('Sorting (phi)')
        ax4.set_ylabel('Midpoint of layer (m)')
        ax4.set_title(fname4)
        ax4.plot(fwsort,midint,'-')
        close()
        #savefig(fname4+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Skewness for layers=================================
        ##This figure is designed to show the relationship between location and skewness
        fig5=figure()
        fname5='Grain Statistics, FW Skewness for layers'
        ax5=subplot(111)
        ax5.set_xlabel('Skewness')
        ax5.set_ylabel('Midpoint of layer (m)')
        ax5.set_title(fname5)
        ax5.plot(fwskew,midint,'-')
        close()
        #savefig(fname5+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Kurtosis for layers=================================
        ##This figure is designed to show the reltaionship between location and kurtosis
        fig6=figure()
        fname6='Grain Statistics, FW Kurtosis for layers'
        ax6=subplot(111)
        ax6.set_xlabel('Kurtosis')
        ax6.set_ylabel('Midpoint of layer (m)')
        ax6.set_title(fname6)
        ax6.plot(fwskew,midint,'-')
        close()
        #savefig(fname6+' Sample%i.png'%ni,dpi=100)
##==========================plot figure for Thickness vs. Time sediment hits bottom=====================
        ##This figures is designed to show the relationship between thickness and time sediment hits bottom
        fig7=figure()
        fname7='Thickness vs. Time sediment hits bottom'
        ax7=subplot(111)
        ax7.set_xscale('log',basex=10)
        ax7.set_xlabel('Time(s)')
        ax7.set_ylabel('Thickness (m)')
        ax7.set_title(fname7)
        ax7.plot(stime1,thickness,'-')
        print shape(stime1)
        close()
        #savefig(fname7+' Sample%i.png'%ni,dpi=100)
        return (phi,Fr)

'''
Created on 2013-1-25
========================================The TsuSedMod tsunami sediment transport model==============================
#=============================================The calculate model part==============================================
##==========================================Read data from file for model===========================================
##============================================Calculate water parameter=============================================
###=================================================calculate s=====================================================
###===============================calculate critical shear stress (N/m2) and settling velocity for each=============
###=================================================calculate zo====================================================
###==================================calculate shear velocity needed to suspend deposit=============================
#==============================================subfunction part=====================================================
##============================================Subfunction grading===================================================
###============================================Calculate part=======================================================
###===========================================plot figure part======================================================
####============================plot figure for Grain Statistics, 1st Moment for layers=============================
####===================================plot figure for Standard Deviation for layers================================
####=====================================plot figure for FW mean for layers=========================================
####=====================================plot figure for FW sorting for layers======================================
####=====================================plot figure for FW Skewness for layers=====================================
####======================================plot figure for FW Kurtosis for layers====================================
####==================================plot figure for Thickness vs. Time sediment hits bottom=======================
##=============================================subfunction sedstats=================================================
##=============================================subfunction readModelFile============================================
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
#===================================================================================================================
#===================================================================================================================
#=========================================The tsunami sediment model based on Bruce Jaffe===========================
#===================================================================================================================
#===================================================================================================================
#This model is based on the paper published in Sedimentary Geology
#A SIMPLE MODEL FOR CALCULATING TSUNAMI FLOW SPEED FROM TSUNAMI DEPOSITS
#               Bruce E. Jaffe, Guy Gelfenbuam
#===================================================================================================================
#Fuction:estimates water velocity from tsunami sediment deposit
#Assumes steady uniform flow
#Functions needed: phi2mm.py, sw_dens0.py, sw_smow.py, KinVisc.py, UstarCrit.py, zoWS.py,dietrichWs.py, Setvel1.py,
#                  tubesetvel.py, linearK.py, parabolicK.py, gelfK.py
#Inputs for this model:         th:deposit thickness (m)
#               	        	h:estimate of water depth (m)
#                       		sal:salinity (psu)
#               	        	wtemp:water temperature (deg. C)
#               	        	nclass:number of sed sizes to consider
#                       		Cb:total bed concentration (1-porosity), usually 0.65
#                       		fr:fractional portion of each grain size in bed
#                               phi:grain size in phi units for each size class
#                       		RhoS:sediment density (kg/m^3) for each size class
#                         		Dmn:mean grain diamter(mm);
#                               ustrc:shear velocity, initial guess input and then allowed to change
#                               ustrcAdjFac:Factor that sets how quickly ustrc is adjusted during iterations
#                               nit:number of iterations to run
#                               diststep:how often to iterate size classes versus concentration
#                               nz:number of vertical bins for calculations
#                        		vk:Van Karmen's constant
#                        		g0:the resuspension coefficient  *** standard value 1.4 x 10^-4 from Hill et al., 1989
#                         		concconvfactor:difference between modeled and observed concentration difference allowed
#                               sizeconvfactor:set the difference between modeled and observed size distribution allowed
#                               var:sets eddy viscosity shape; 1=linear, 2=parbolic,3=[Gelfenbaum]
#Used:                          ssfrwanted:suspended size distribution wanted
#                               ssconverge:keeps track of whether amount of sediment in each size class is acceptable
#                               ssoffa:array to track of how far off the suspended sediment concentration is from desired concentration
#                               ssfra:array to track of suspended sediment concentration in each size class
#                               fra:array to track of bed sediment concentration in each size class
#                               offa:array to track how far total suspended loasource_distribution.pyd is off from load required to create deposit
#                               ustrca:array to track how  u*c changes with iteration
#Parameter calculated by import function:
#                               d:grain diameter (set to m)
#                         		rhow:density of water (kg/m^3)sedstat
#                       		kinvisc:kinematic viscosity of water (m^2/s), usually about 10^-6
#                               UcritDmn:critical shear velocity for mean grain size
#                               zotot:Z naut total, from Nickaradse bed roughness and saltation bed roughness
#                         		ws:settling velocity (m/s)
#                        		ucrit:critical shear velocity for initiation of sediment transport
#                               z:log space elevation values where concentrations and velocities are calculated
#                          		K:eddie viscosity (m/s^2)
#                       		zo:Z naut, zero velocity intercept
#parameter calculated in this model:
#                               zoN- Nickaradse grain roughness
#      	                        ustrc:current shear velocity, u*c [m/s], initial value 0.1 m/s
#		                        s:ratio of density of sediment to density of fluid density for each size class
#   	                        thload:sediment load needed to get observed tsunami deposit thickness
#		                        tcrit:critical shear stress (N/m^2)
#		                        taub:bottom shear stress due to currents (N/m^2)
#		                        taustar:bottom shear stress normalized by critical shear stress
#		                        S:excess shear stress
#                               Ca:reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters
#		                        Ki:integral of 1/eddie viscosity
#                               C:Suspended-sediment profile (calculated for each size class)
#                               G:Suspended load (calculated for each size class)
#                               off:normalized difference between suspended load and load required to create deposit
#		                        spd:mean current speed calculated using current shear velocity (m/s)
#                               froude:Froude number calculated using maximum speed
#Options:
#                               'infile':specifies filename of model input file.
#                               'outfile':writes results to comma-separated (.csv) file. If the specified output file already exists,
#                               results will be appended to current file.
#                               'grading':plots results of grading function . The default output depth of this function is 0.01 m.  To
#                                modify the defaut, the user can specify an additional arguement of a vector of depths or a single value.
#Usage examples:
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5;
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('infile','C:\exampleFile.xls');
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',0.01);
#                                   where 0.01 specifies set interval (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',[0.01,0.03,0.04,0.8]);
#                                   where [ ] specifies variable intervals (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('outfile','modelResults.csv');
#Version information:
#written in March 12, 2001 by Bruce Jaffe, based on code written by Chris Sherwood, documented in July 3,2001
#Significant code clean-ups and user friendly input added by Andrew Stevens from March to October, 2007
#modfied on 11/29/07 by Bruce Jaffe, modification of Andrew Stevens' cleaned-up code including: allowing single size class runs
#modified 4/29/08 to output weight percent (Bruce Jaffe) and to fix grading function interval error (Mark Buckley)
#made in python 1/25/2013 by Hui Tang (tanghui@vt.edu)
#==============================================================================================================================================
#==============================================================The caculate model part=========================================================
from pylab import *
# from function import *
# from sedstats import *
# from numpy import *
# from grading import *
# from ReadCSV import *
import csv
import os
def Tsusedmod(name):
##============================================================Read data from file for model====================================================
    print("\n ====================================running==================================== \n")
    #read in model inputs from excel file
    separator=':'
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='number_of_iterations': nit=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='number_of_vertical_bins': nz=float(v)
                elif n=='resuspension_coefficient': g0=float(v)
                elif n=='concentration_convergence_factor': concconvfactor=float(v)
                elif n=='size_convergence_factor': sizeconvfactor=float(v)
                elif n=='shear_velocity': ustrc=float(v)
                elif n=='bed_roughness': zotot=float(v)
                elif n=='settling_velocity_type': setType=float(v)
                elif n=='eddy_viscosity_shape': var=float(v)
                #site info
                elif n=='site_description': t=str(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='bed_concentration': Cb=float(v)
                #data
                elif n=='Number_sample':num=int(v)
                elif n=='Number_grainsize':nclass=int(v)
                elif n=='data_dimension':l=int(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
    filename=name
    data=readCSV1(outp_path+filename,separator=';')
    phi=data[:,0]
    fr1=data[:,1]
    fr=zeros(shape=(len(fr1),1))
    th=data[0,2]
    h=data[0,3]
    x=data[0,4]
    #end data input
    nclass=len(fr1) #number of sed sizes to consider
    fr=fr1.reshape(len(fr1),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    ssfrwanted=fr #susp sed size distribution wanted, keep this for comparison with model results
    ssconverge=ones(shape=(nclass,1))#initiate ss converge to 1 for all size classes (does  converge)
    for i in range(len(fr)):
        if(fr[i]!=0):
            ssconverge[i]=0 #reset ss converge values to 0 (does not converge) for all size classes with sediments
    sizeconverge=sizeconvfactor*ones(shape=(nclass,1))#set percent difference in size distribution allowed, same value for all sizes
    D=phi2mm(phi) #convert grain size to mm for settling velocity calculation using Dietrich's formulae
    d=0.001*D #convert classes to m for calculating settling velocity using other forumlae and for UstarC calculation
    if(nclass==1):
        dMean_phi=phi[0]
    else:
        sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
        dMean_phi=sedstat[0]
        m3a=sedstat[1]
        m4a=sedstat[2]
        stdeva=sedstat[3]
        fwmedian1=sedstat[4]
        fwsort1=sedstat[5]
        fwmean1=sedstat[6]
        fwskew1=sedstat[7]
        fwkurt1=sedstat[8]
    dMean_mm=phi2mm(dMean_phi) # convertmean grain size in phi to mm
    Dmn=0.001*dMean_mm  # convert mean grain size to m for calculating bed roughess
    rhos=Rhos*ones(shape=(nclass,1)) #set sediment density, same value for all sizes
#==========================================================Calculate water parameter=============================================================
    rhow=sw_dens0(sal,wtemp)  #calculate the density of seawater [kg/m^3]
    kinvisc=KinVisc(wtemp)    #calculate the kinematic viscosity of the water
###===============================================================calculate s=====================================================================
    s=rhos/rhow
###===============================calculate critical shear stress (N/m2) and settling velocity for each===========================================
    #size class
    ws=zeros(nclass)
    Ucrit=zeros(nclass)
    tcrit=zeros(nclass)
    if(setType==1.0):
        for i in range(nclass):
            ws[i]=-tubesetvel(phi[i],rhos,wtemp,sal)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==2.0):
        for i in range(nclass):
            ws[i]=-dietrichWs(wtemp,d[i],rhos[i],rhow/1000,0.005,0.7,3.5)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==3.0):
        for i in range(nclass):
            ws[i]=Setvel1(d[i],kinvisc,s[i])
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    thload=th*Cb #suspended sediment load to make deposit of observed thickness, thickness times Cb = 1-porosity
###============================================================calculate zo=======================================================================
    zoN=Dmn/30.0  #Nickaradse bed roughness
    #guess ustrc to make deposit, need this to calculate zos, bed roughness from moving flat bed, saltation
    #Note:Make the initial ustrc an input, also make ustrcAdjFac an input and delelet lines 228 and 229
    ustrcAdjFac=0.01
    UcritDmn=UstarCrit(Dmn,kinvisc,2650/rhow) #need critical shear velocity for mean grain size to calculate bed roughness
    if(zotot!=1):
        zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #zo total (Nickaradse bed roughness and Wiberg/Smith bed roughness from saltation,
###===============================================calculate shear velocity needed to suspend deposit===============================================
    #Begin loop to determine shear velocity needed to suspend deposit
    #This loops used to assume that the bed sediment disribution is the same as the suspended sediment distribution
    #This gives reasonable results when the size of the bed material is number of iterations to run
    diststep=1 #adjusts sediment grain size distribution every 2nd iteration
    ssoffa=zeros(shape=(nit/diststep,nclass)) # initiate array tracking ss deviation from desired concentration for iteration=1:nit
    ssfra=ones(shape=(nit/diststep,nclass)) #initiate array tracking ss concentrations for iteration=1:nit
    fra=zeros(shape=(nit/diststep,nclass)) #initiate array tracking bed sediment concentrations for iteration=1:nit
    offa=zeros(nit)#inititate array to track how far total suspended load is off from desired suspended load
    ustrca=zeros(nit)#inititate array to track how  u*c changes with iteration                                            #ref. Wiberg and Rubin, 1989 JGR
    cconverge=0
    for iteration in range(int(nit)):
        fr=fr/sum(fr) #make sure fractions add up to 1
        if(all(ssconverge) and cconverge):
            break
        taub=rhow*ustrc**2 # bottom stress
        taustar=taub/tcrit
        S=taustar-1     #normalized excess shear stress
        for i in range(len(S)):#returns zero when taub<=tcrit
            if S[i]<0:
                S[i]=0
            else:s[i]=s[i]
        S=S.reshape(len(S),1)
        Ca=zeros(shape=(len(s),1))
        for i in range(len(S)):
            Ca[i]=g0*Cb*fr[i]*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        # resuspension coefficient(g0) for flat moveable bed Madsen 94a "Susp. Sed. Transport on the inner shelf waters during extreme storms"
        # ICCE pg 1849 Madsen, Chisholm, Wright based on Wikramanayake 92 PhD thesis first shot, might want to change this using mean current litterature
        # log-spaced z grid
        z=logspace(log10(zotot),log10(h),num=nz,endpoint=True,base=10.0)
        z=z.reshape(len(z),1)
#*************************************************elevation of suspended load, SL******************************************
        def diff(z):
            diff=zeros(shape=(len(z)-1,1))
            for i in range(0,len(z)-1):
                diff[i]=z[i+1]-z[i]
            return (diff)
        zsl=z[0:nz-1]+diff(z)
        diff1=diff(z)
#**************************************elevation of speed calculated from eddy viscosity***********************************
        zmid=(z[0:nz-1]+z[1:nz])/2
#*******************************************calculate eddy viscosity profile***********************************************
        #var=1 is linear eddy viscosity profile, var=2 is parabolic profile, var=3 is gelf profile
        if var==1.0:
            K=linearK(z,h,ustrc)
        elif var==2.0:
            K=parabolicK(z,h,ustrc)
        elif var==3.0:
            K=gelfK(z,h,ustrc)
        #K=(var==1)*linearK(z,h,ustrc)+(var==2)*parabolicK(z,h,ustrc)+(var==3)*gelfK(z,h,ustrc)
#*******************************************Calculate current and eddie viscosity profiles*********************************
        #Integral of 1/K
        K=K.reshape(len(K),1)
        Ki=cumsum(2.0*diff(z)/(K[0:nz-1]+K[1:nz]))
#*********************************calculate the suspended-sediment profile and sediment load for each class****************
        C=zeros(shape=(nz,nclass))
        sl=zeros(shape=(nz-1,nclass))
        sl1=zeros(shape=(nz,nclass))
        G=zeros(nclass)
        for i in range(nclass):
            #Suspended-sediment profile
            C[:len(Ca[0]),i]=Ca[i][0]
            C[len(Ca[0]):,i]=Ca[i][0]*exp(ws[i]*K[i])
            for j in range(int(nz-2)):
                sl[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
                sl1[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
#******************************************Depth-intergral of sediment profile*********************************************
            G[i] = sum(sl[:,i])
#************************************************** set cconverge to 0*****************************************************
        cconverge=0
        off=sum(G)/thload-1
        offa[iteration]=off
        if(abs(off)*100<concconvfactor):
            cconverge=1
        if(all(ssconverge) and cconverge): #if both conc and size are good enough, get outta here
            break
        ustrca[iteration]=ustrc
        ustrc=ustrc*(1-ustrcAdjFac*off) #change ustrc based on how far off suspended load is from desired value
        if(zotot!=1):
            zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #recalculate bed roughness
#*******check for convergence for size of suspended sediment and adjust bed sediment distribution every diststep***********
        def fix(n):
            if(n>0):
                if(round(n)>n):
                    fix=round(n)-1
                else:fix=round(n)
            if(n<0):
                if(round(n)<n):
                    fix=round(n)+1
            return(n)
        doss=diststep*fix(iteration/diststep)
        if(doss==iteration): #redo distribution every diststep steps
            ssfr=G/sum(G) #calculate susp sed size distribution for model
            ssfra[doss/diststep,:]=ssfr
            fra[doss/diststep,:]=fr.reshape(1,len(fr))
            ssfr=ssfr.reshape(len(ssfr),1)
            #keep track of whether model is converging ** array indexes
            #different from iteration
            for k in range(nclass):
                if(ssfrwanted[k]!=0): # only adjust for classes with material in them
                    ssoff=ssfr[k][0]/ssfrwanted[k]-1
                    ssoffa[doss/diststep,k]=ssoff #change bed grain size distribution based on suspended load distribution
                    fr[k]=fr[k]*(1-0.01*ssoff)
                    if (abs(ssoff)*100<sizeconverge[k]): #check for convergence for each size class
                        ssconverge[k]=1.0 #set to 1 if OK
                if(all(ssconverge) and cconverge):
                    break #if both conc and size are good enough, get outta here
    spd=zeros(len(z)-1)
    diff1=diff(z)
    diff1=diff1.reshape(len(diff1),)
    spd[:]=cumsum(diff1*ustrc*ustrc/gelfK(zmid,h,ustrc))#speed calculated using eddy viscosity profile integration
    if(iteration==nit):
        print('Model may not have converged, use extreme caution with these resuls')
    #stop clock
    predictedssload=sum(G[:])
    maximumspeed=max(spd)
    spd1=0.5*(spd[0:len(spd)-1]+spd[1:len(spd)])
    spd1=spd1.reshape(len(spd1),1)
    diff2=zeros(shape=(len(diff1)-1,1))
    diff2=diff1[1:len(diff1)]
    diff2=diff2.reshape(len(diff2),1)
    avgspeed=sum((spd1) * diff2)/h
    MaxFroude=maximumspeed/sqrt(9.8*h)
    AvgFroude=avgspeed/sqrt(9.8*h)
    if (var==1): ed='Linear'
    elif (var==2): ed='Parabolic'
    elif (var==3): ed='Gelfenbaum'
    print('Deposit thickness:th=%4.3f'%th+'m \n')
    print('Water depth:h=%4.1f'%h+'m \n')
    print('Size classes: %i'%nclass+'\n')
    print('D mean: %.3e'%Dmn+'m \n')
    print('Viscosity Profile:%s'%ed+'\n')
    print('\n Model Results \n')
    print('Shear velocity:%4.2f'%ustrc+'m/s \n')
    print('Bed Roughness:%.3e'%zotot+'\n')
    print('Sediment load needed:%4.3f'%thload+'m \n')
    print('Sediment load predicted:%4.3f'%predictedssload+'m \n')
    print('Max. speed:%5.2f'%maximumspeed+'m/s \n')
    print('Avg. speed:%5.2f'%avgspeed+'m/s \n')
    print('Max. Froude:%4.2f'%MaxFroude+'\n')
    print('Avg. Froude:%4.2f'%AvgFroude+'\n')
    fp=open('result.txt','wb')
    fp.write("Model stats \n")
    fp.write("datestr(now) \n")
    fp.write("iterations:")
    fp.write("%i \n"%iteration)
    fp.write("Input Data \n")
    fp.write("Deposit thickness (m):")
    fp.write("%4.3f \n"%th)
    fp.write("Depth (m):")
    fp.write("%4.1f \n"%h)
    fp.write("size classes:")
    fp.write("%i \n"%nclass)
    fp.write("Mean grain diameter (m):")
    fp.write("%.3e \n"%Dmn)
    fp.write("Eddy viscosity profile:")
    fp.write("%s \n"%ed)
    fp.write("Model Results: \n")
    fp.write("ustarc (shear velocity) (m/s):")
    fp.write("%4.2f \n"%ustrc)
    fp.write("ztot (bed roughness):")
    fp.write("%.3e \n"%zotot)
    fp.write("thload (sediment load need) (m):")
    fp.write("%4.3f \n"%thload)
    fp.write("predictload (sediment load predicted) (m): ")
    fp.write("%4.3f \n"%predictedssload)
    fp.write("max. speed (m/s):")
    fp.write("%5.2f \n"%maximumspeed)
    fp.write("avg. speed (m/s):")
    fp.write("%5.2f \n"%avgspeed)
    fp.write("max. Froude:")
    fp.write("%4.2f \n"%MaxFroude)
    fp.write("avg. Froude:")
    fp.write('%4.2f \n'%AvgFroude)
    fp.write("Detail of structure: \n")
    fp.write("ssfrwanted (suspended size distribution wanted): \n")
    fp.write("%s \n"%ssfrwanted)
    fp.write("ssconverge (keeps track of whether amount of sediment in each size class is acceptable): \n")
    fp.write("%s \n"%ssconverge)
    fp.write("cconverge: \n")
    fp.write("%s \n"%cconverge)
    fp.write("ssfra(array to track of suspended sediment concentration in each size class): \n")
    fp.write("%s \n"%ssfra)
    fp.write("offa(array to track how far total suspended load is off from load required to create deposit): \n")
    fp.write("%s \n"%offa)
    fp.write("C(Suspended-sediment profile (calculated for each size class)): \n")
    fp.write("%s \n"%C)
    fp.write("G(Suspended load (calculated for each size class)): \n")
    fp.write("%s \n"%G)
    fp.write("ws(settling velocity (m/s)): \n")
    fp.write("%s \n"%ws)
    fp.write("SL: \n")
    fp.write("%s \n"%sl)
    fp.write("z(log space elevation values where concentrations and velocities are calculated): \n")
    fp.write("%s \n"%s)
    fp.write("zsl: \n")
    fp.write("%s \n"%zsl)
    fp.write("zmid: \n")
    fp.write("%s \n"%zmid)
    fp.write("ustrca(array to track how  u*c changes with iteration): \n")
    fp.write("%s \n"%ustrca)
    fp.write("S(excess shear stress): \n")
    fp.write("%s \n"%S)
    fp.write("spd(mean current speed calculated using current shear velocity (m/s)): \n")
    fp.write("%s \n"%spd)
    fp.write("results of the structure \n")
    fp.write("zoN(Nickaradse grain roughness): \n")
    fp.write("%s \n"%zoN)
    fp.write("tcrit(critical shear stress (N/m^2)): \n")
    fp.write("%s \n"%tcrit)
    fp.write("taub(bottom shear stress due to currents (N/m^2)): \n")
    fp.write("%s \n"%taub)
    fp.write("taustar(bottom shear stress normalized by critical shear stress): \n")
    fp.write("%s \n"%taustar)
    fp.write("phi(grain size in phi units for each size class): \n")
    fp.write("%s \n"%phi)
    fp.write("Ca(reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters): \n")
    fp.write("%s \n"%Ca)
    fp.write("Ki(integral of 1/eddie viscosity): \n")
    fp.write("%s \n"%Ki)
    fp.write("off(normalized difference between suspended load and load required to create deposit): \n")
    fp.write("%s \n"%off)
    fp.write("\n End")
    fp.write('General results \n')
    fp.write('Predictload (kg/m3) \n')
    fp.write('%f \n'%predictedssload)
    fp.write('maximum speed \n')
    fp.write('%f \n'%maximumspeed)
    fp.write('average speed \n')
    fp.write('%f \n'%avgspeed)
    fp.write('Max froude \n')
    fp.write('%f \n'%MaxFroude)
    fp.write('average froude \n')
    fp.write('%f \n'%AvgFroude)
    fp.close()
    return(sl1,z,phi)


'''
Created on 2013-1-25
========================================The TsuSedMod tsunami sediment transport model==============================
#=============================================The calculate model part==============================================
##==========================================Read data from file for model===========================================
##============================================Calculate water parameter=============================================
###=================================================calculate s=====================================================
###===============================calculate critical shear stress (N/m2) and settling velocity for each=============
###=================================================calculate zo====================================================
###==================================calculate shear velocity needed to suspend deposit=============================
#==============================================subfunction part=====================================================
##============================================Subfunction grading===================================================
###============================================Calculate part=======================================================
###===========================================plot figure part======================================================
####============================plot figure for Grain Statistics, 1st Moment for layers=============================
####===================================plot figure for Standard Deviation for layers================================
####=====================================plot figure for FW mean for layers=========================================
####=====================================plot figure for FW sorting for layers======================================
####=====================================plot figure for FW Skewness for layers=====================================
####======================================plot figure for FW Kurtosis for layers====================================
####==================================plot figure for Thickness vs. Time sediment hits bottom=======================
##=============================================subfunction sedstats=================================================
##=============================================subfunction readModelFile============================================
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
#===================================================================================================================
#===================================================================================================================
#=========================================The tsunami sediment model based on Bruce Jaffe===========================
#===================================================================================================================
#===================================================================================================================
#This model is based on the paper published in Sedimentary Geology
#A SIMPLE MODEL FOR CALCULATING TSUNAMI FLOW SPEED FROM TSUNAMI DEPOSITS
#               Bruce E. Jaffe, Guy Gelfenbuam
#===================================================================================================================
#Fuction:estimates water velocity from tsunami sediment deposit
#Assumes steady uniform flow
#Functions needed: phi2mm.py, sw_dens0.py, sw_smow.py, KinVisc.py, UstarCrit.py, zoWS.py,dietrichWs.py, Setvel1.py,
#                  tubesetvel.py, linearK.py, parabolicK.py, gelfK.py
#Inputs for this model:         th:deposit thickness (m)
#               	        	h:estimate of water depth (m)
#                       		sal:salinity (psu)
#               	        	wtemp:water temperature (deg. C)
#               	        	nclass:number of sed sizes to consider
#                       		Cb:total bed concentration (1-porosity), usually 0.65
#                       		fr:fractional portion of each grain size in bed
#                               phi:grain size in phi units for each size class
#                       		RhoS:sediment density (kg/m^3) for each size class
#                         		Dmn:mean grain diamter(mm);
#                               ustrc:shear velocity, initial guess input and then allowed to change
#                               ustrcAdjFac:Factor that sets how quickly ustrc is adjusted during iterations
#                               nit:number of iterations to run
#                               diststep:how often to iterate size classes versus concentration
#                               nz:number of vertical bins for calculations
#                        		vk:Van Karmen's constant
#                        		g0:the resuspension coefficient  *** standard value 1.4 x 10^-4 from Hill et al., 1989
#                         		concconvfactor:difference between modeled and observed concentration difference allowed
#                               sizeconvfactor:set the difference between modeled and observed size distribution allowed
#                               var:sets eddy viscosity shape; 1=linear, 2=parbolic,3=[Gelfenbaum]
#Used:                          ssfrwanted:suspended size distribution wanted
#                               ssconverge:keeps track of whether amount of sediment in each size class is acceptable
#                               ssoffa:array to track of how far off the suspended sediment concentration is from desired concentration
#                               ssfra:array to track of suspended sediment concentration in each size class
#                               fra:array to track of bed sediment concentration in each size class
#                               offa:array to track how far total suspended loasource_distribution.pyd is off from load required to create deposit
#                               ustrca:array to track how  u*c changes with iteration
#Parameter calculated by import function:
#                               d:grain diameter (set to m)
#                         		rhow:density of water (kg/m^3)sedstat
#                       		kinvisc:kinematic viscosity of water (m^2/s), usually about 10^-6
#                               UcritDmn:critical shear velocity for mean grain size
#                               zotot:Z naut total, from Nickaradse bed roughness and saltation bed roughness
#                         		ws:settling velocity (m/s)
#                        		ucrit:critical shear velocity for initiation of sediment transport
#                               z:log space elevation values where concentrations and velocities are calculated
#                          		K:eddie viscosity (m/s^2)
#                       		zo:Z naut, zero velocity intercept
#parameter calculated in this model:
#                               zoN- Nickaradse grain roughness
#      	                        ustrc:current shear velocity, u*c [m/s], initial value 0.1 m/s
#		                        s:ratio of density of sediment to density of fluid density for each size class
#   	                        thload:sediment load needed to get observed tsunami deposit thickness
#		                        tcrit:critical shear stress (N/m^2)
#		                        taub:bottom shear stress due to currents (N/m^2)
#		                        taustar:bottom shear stress normalized by critical shear stress
#		                        S:excess shear stress
#                               Ca:reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters
#		                        Ki:integral of 1/eddie viscosity
#                               C:Suspended-sediment profile (calculated for each size class)
#                               G:Suspended load (calculated for each size class)
#                               off:normalized difference between suspended load and load required to create deposit
#		                        spd:mean current speed calculated using current shear velocity (m/s)
#                               froude:Froude number calculated using maximum speed
#Options:
#                               'infile':specifies filename of model input file.
#                               'outfile':writes results to comma-separated (.csv) file. If the specified output file already exists,
#                               results will be appended to current file.
#                               'grading':plots results of grading function . The default output depth of this function is 0.01 m.  To
#                                modify the defaut, the user can specify an additional arguement of a vector of depths or a single value.
#Usage examples:
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5;
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('infile','C:\exampleFile.xls');
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',0.01);
#                                   where 0.01 specifies set interval (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',[0.01,0.03,0.04,0.8]);
#                                   where [ ] specifies variable intervals (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('outfile','modelResults.csv');
#Version information:
#written in March 12, 2001 by Bruce Jaffe, based on code written by Chris Sherwood, documented in July 3,2001
#Significant code clean-ups and user friendly input added by Andrew Stevens from March to October, 2007
#modfied on 11/29/07 by Bruce Jaffe, modification of Andrew Stevens' cleaned-up code including: allowing single size class runs
#modified 4/29/08 to output weight percent (Bruce Jaffe) and to fix grading function interval error (Mark Buckley)
#made in python 1/25/2013 by Hui Tang (tanghui@vt.edu)
#==============================================================================================================================================
#==============================================================The caculate model part=========================================================
from pylab import *
# from function import *
# from sedstats import *
# from numpy import *
# from grading import *
# from ReadCSV import *
import csv
import os
def Tsusedmod1(name):
##============================================================Read data from file for model====================================================
    print("\n ====================================running==================================== \n")
    #read in model inputs from excel file
    separator=':'
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='number_of_iterations': nit=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='number_of_vertical_bins': nz=float(v)
                elif n=='resuspension_coefficient': g0=float(v)
                elif n=='concentration_convergence_factor': concconvfactor=float(v)
                elif n=='size_convergence_factor': sizeconvfactor=float(v)
                elif n=='shear_velocity': ustrc=float(v)
                elif n=='bed_roughness': zotot=float(v)
                elif n=='settling_velocity_type': setType=float(v)
                elif n=='eddy_viscosity_shape': var=float(v)
                #site info
                elif n=='site_description': t=str(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='bed_concentration': Cb=float(v)
                #data
                elif n=='Number_sample':num=int(v)
                elif n=='Number_grainsize':nclass=int(v)
                elif n=='data_dimension':l=int(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
    filename=name
    data=readCSV1(filename,separator=';')
    phi=data[:,0]
    fr1=data[:,1]
    fr=zeros(shape=(len(fr1),1))
    th=data[0,2]
    h=data[0,3]
    x=data[0,4]
    #end data input
    nclass=len(fr1) #number of sed sizes to consider
    fr=fr1.reshape(len(fr1),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    ssfrwanted=fr #susp sed size distribution wanted, keep this for comparison with model results
    ssconverge=ones(shape=(nclass,1))#initiate ss converge to 1 for all size classes (does  converge)
    for i in range(len(fr)):
        if(fr[i]!=0):
            ssconverge[i]=0 #reset ss converge values to 0 (does not converge) for all size classes with sediments
    sizeconverge=sizeconvfactor*ones(shape=(nclass,1))#set percent difference in size distribution allowed, same value for all sizes
    D=phi2mm(phi) #convert grain size to mm for settling velocity calculation using Dietrich's formulae
    d=0.001*D #convert classes to m for calculating settling velocity using other forumlae and for UstarC calculation
    if(nclass==1):
        dMean_phi=phi[0]
    else:
        sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
        dMean_phi=sedstat[0]
        m3a=sedstat[1]
        m4a=sedstat[2]
        stdeva=sedstat[3]
        fwmedian1=sedstat[4]
        fwsort1=sedstat[5]
        fwmean1=sedstat[6]
        fwskew1=sedstat[7]
        fwkurt1=sedstat[8]
    dMean_mm=phi2mm(dMean_phi) # convertmean grain size in phi to mm
    Dmn=0.001*dMean_mm  # convert mean grain size to m for calculating bed roughess
    rhos=Rhos*ones(shape=(nclass,1)) #set sediment density, same value for all sizes
#==========================================================Calculate water parameter=============================================================
    rhow=sw_dens0(sal,wtemp)  #calculate the density of seawater [kg/m^3]
    kinvisc=KinVisc(wtemp)    #calculate the kinematic viscosity of the water
###===============================================================calculate s=====================================================================
    s=rhos/rhow
###===============================calculate critical shear stress (N/m2) and settling velocity for each===========================================
    #size class
    ws=zeros(nclass)
    Ucrit=zeros(nclass)
    tcrit=zeros(nclass)
    if(setType==1.0):
        for i in range(nclass):
            ws[i]=-tubesetvel(phi[i],rhos,wtemp,sal)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==2.0):
        for i in range(nclass):
            ws[i]=-dietrichWs(wtemp,d[i],rhos[i],rhow,0.005,0.7,3.5)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==3.0):
        for i in range(nclass):
            ws[i]=Setvel1(d[i],kinvisc,s[i])
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    thload=th*Cb #suspended sediment load to make deposit of observed thickness, thickness times Cb = 1-porosity
###============================================================calculate zo=======================================================================
    zoN=Dmn/30.0  #Nickaradse bed roughness
    #guess ustrc to make deposit, need this to calculate zos, bed roughness from moving flat bed, saltation
    #Note:Make the initial ustrc an input, also make ustrcAdjFac an input and delelet lines 228 and 229
    ustrcAdjFac=0.01
    UcritDmn=UstarCrit(Dmn,kinvisc,2650/rhow) #need critical shear velocity for mean grain size to calculate bed roughness
    if(zotot!=1):
        zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #zo total (Nickaradse bed roughness and Wiberg/Smith bed roughness from saltation,
###===============================================calculate shear velocity needed to suspend deposit===============================================
    #Begin loop to determine shear velocity needed to suspend deposit
    #This loops used to assume that the bed sediment disribution is the same as the suspended sediment distribution
    #This gives reasonable results when the size of the bed material is number of iterations to run
    diststep=1 #adjusts sediment grain size distribution every 2nd iteration
    ssoffa=zeros(shape=(nit/diststep,nclass)) # initiate array tracking ss deviation from desired concentration for iteration=1:nit
    ssfra=ones(shape=(nit/diststep,nclass)) #initiate array tracking ss concentrations for iteration=1:nit
    fra=zeros(shape=(nit/diststep,nclass)) #initiate array tracking bed sediment concentrations for iteration=1:nit
    offa=zeros(nit)#inititate array to track how far total suspended load is off from desired suspended load
    ustrca=zeros(nit)#inititate array to track how  u*c changes with iteration                                            #ref. Wiberg and Rubin, 1989 JGR
    cconverge=0
    for iteration in range(int(nit)):
        fr=fr/sum(fr) #make sure fractions add up to 1
        if(all(ssconverge) and cconverge):
            break
        taub=rhow*ustrc**2 # bottom stress
        taustar=taub/tcrit
        S=taustar-1     #normalized excess shear stress
        for i in range(len(S)):#returns zero when taub<=tcrit
            if S[i]<0:
                S[i]=0
            else:s[i]=s[i]
        S=S.reshape(len(S),1)
        Ca=zeros(shape=(len(s),1))
        for i in range(len(S)):
            Ca[i]=g0*Cb*fr[i]*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        # resuspension coefficient(g0) for flat moveable bed Madsen 94a "Susp. Sed. Transport on the inner shelf waters during extreme storms"
        # ICCE pg 1849 Madsen, Chisholm, Wright based on Wikramanayake 92 PhD thesis first shot, might want to change this using mean current litterature
        # log-spaced z grid
        z=logspace(log10(zotot),log10(h),num=nz,endpoint=True,base=10.0)
        z=z.reshape(len(z),1)
#*************************************************elevation of suspended load, SL******************************************
        def diff(z):
            diff=zeros(shape=(len(z)-1,1))
            for i in range(0,len(z)-1):
                diff[i]=z[i+1]-z[i]
            return (diff)
        zsl=z[0:nz-1]+diff(z)
        diff1=diff(z)
#**************************************elevation of speed calculated from eddy viscosity***********************************
        zmid=(z[0:nz-1]+z[1:nz])/2
        zmid1=zmid.reshape(len(zmid))
#*******************************************calculate eddy viscosity profile***********************************************
        #var=1 is linear eddy viscosity profile, var=2 is parabolic profile, var=3 is gelf profile
        if var==1.0:
            K=linearK(z,h,ustrc)
        elif var==2.0:
            K=parabolicK(z,h,ustrc)
        elif var==3.0:
            K=gelfK(z,h,ustrc)
        #K=(var==1)*linearK(z,h,ustrc)+(var==2)*parabolicK(z,h,ustrc)+(var==3)*gelfK(z,h,ustrc)
#*******************************************Calculate current and eddie viscosity profiles*********************************
        #Integral of 1/K
        K=K.reshape(len(K),1)
        Ki=cumsum(2.0*diff(z)/(K[0:nz-1]+K[1:nz]))
#*********************************calculate the suspended-sediment profile and sediment load for each class****************
        C=zeros(shape=(nz,nclass))
        sl=zeros(shape=(nz-1,nclass))
        sl1=zeros(shape=(nz,nclass))
        G=zeros(nclass)
        for i in range(nclass):
            #Suspended-sediment profile
            C[:len(Ca[0]),i]=Ca[i][0]
            C[len(Ca[0]):,i]=Ca[i][0]*exp(ws[i]*K[i])
            for j in range(int(nz-2)):
                sl[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
                sl1[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
#******************************************Depth-intergral of sediment profile*********************************************
            G[i] = sum(sl[:,i])
#************************************************** set cconverge to 0*****************************************************
        cconverge=0
        off=sum(G)/thload-1
        offa[iteration]=off
        if(abs(off)*100<concconvfactor):
            cconverge=1
        if(all(ssconverge) and cconverge): #if both conc and size are good enough, get outta here
            break
        ustrca[iteration]=ustrc
        ustrc=ustrc*(1-ustrcAdjFac*off) #change ustrc based on how far off suspended load is from desired value
        if(zotot!=1):
            zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #recalculate bed roughness
#*******check for convergence for size of suspended sediment and adjust bed sediment distribution every diststep***********
        def fix(n):
            if(n>0):
                if(round(n)>n):
                    fix=round(n)-1
                else:fix=round(n)
            if(n<0):
                if(round(n)<n):
                    fix=round(n)+1
            return(n)
        doss=diststep*fix(iteration/diststep)
        if(doss==iteration): #redo distribution every diststep steps
            ssfr=G/sum(G) #calculate susp sed size distribution for model
            ssfra[doss/diststep,:]=ssfr
            fra[doss/diststep,:]=fr.reshape(1,len(fr))
            ssfr=ssfr.reshape(len(ssfr),1)
            #keep track of whether model is converging ** array indexes
            #different from iteration
            for k in range(nclass):
                if(ssfrwanted[k]!=0): # only adjust for classes with material in them
                    ssoff=ssfr[k][0]/ssfrwanted[k]-1
                    ssoffa[doss/diststep,k]=ssoff #change bed grain size distribution based on suspended load distribution
                    fr[k]=fr[k]*(1-0.01*ssoff)
                    if (abs(ssoff)*100<sizeconverge[k]): #check for convergence for each size class
                        ssconverge[k]=1.0 #set to 1 if OK
                if(all(ssconverge) and cconverge):
                    break #if both conc and size are good enough, get outta here
    spd=zeros(len(z)-1)
    diff1=diff(z)
    diff1=diff1.reshape(len(diff1),)
    spd[:]=cumsum(diff1*ustrc*ustrc/gelfK(zmid,h,ustrc))#speed calculated using eddy viscosity profile integration
    if(iteration==nit):
        print('Model may not have converged, use extreme caution with these resuls')
    #stop clock
    predictedssload=sum(G[:])
    maximumspeed=max(spd)
    spd1=0.5*(spd[0:len(spd)-1]+spd[1:len(spd)])
    spd1=spd1.reshape(len(spd1),1)
    #print zsl
    #print spd1
    #print shape(spd)
    #print shape(z)
    #fig3=figure()
    #plot(spd,zsl)
    #title("Velocity Ditribution in vertical direction at x=%d"%x+"m")
    #ylabel("Distance from bed (m)")
    #xlabel("Veloctiy (m/s) ")
    #xlim(0,8)
    #savefig(filename+' .png',dpi=100)
    diff2=zeros(shape=(len(diff1)-1,1))
    diff2=diff1[1:len(diff1)]
    diff2=diff2.reshape(len(diff2),1)
    avgspeed=sum((spd1) * diff2)/h
    MaxFroude=maximumspeed/sqrt(9.8*h)
    AvgFroude=avgspeed/sqrt(9.8*h)
    if (var==1): ed='Linear'
    elif (var==2): ed='Parabolic'
    elif (var==3): ed='Gelfenbaum'
    print('Deposit thickness:th=%4.3f'%th+'m \n')
    print('Water depth:h=%4.1f'%h+'m \n')
    print('Size classes: %i'%nclass+'\n')
    print('D mean: %.3e'%Dmn+'m \n')
    print('Viscosity Profile:%s'%ed+'\n')
    print('\n Model Results \n')
    print('Shear velocity:%4.2f'%ustrc+'m/s \n')
    print('Bed Roughness:%.3e'%zotot+'\n')
    print('Sediment load needed:%4.3f'%thload+'m \n')
    print('Sediment load predicted:%4.3f'%predictedssload+'m \n')
    print('Max. speed:%5.2f'%maximumspeed+'m/s \n')
    print('Avg. speed:%5.2f'%avgspeed+'m/s \n')
    print('Max. Froude:%4.2f'%MaxFroude+'\n')
    print('Avg. Froude:%4.2f'%AvgFroude+'\n')
    fp=open('result.txt','wb')
    fp.write("Model stats \n")
    fp.write("datestr(now) \n")
    fp.write("iterations:")
    fp.write("%i \n"%iteration)
    fp.write("Input Data \n")
    fp.write("Deposit thickness (m):")
    fp.write("%4.3f \n"%th)
    fp.write("Depth (m):")
    fp.write("%4.1f \n"%h)
    fp.write("size classes:")
    fp.write("%i \n"%nclass)
    fp.write("Mean grain diameter (m):")
    fp.write("%.3e \n"%Dmn)
    fp.write("Eddy viscosity profile:")
    fp.write("%s \n"%ed)
    fp.write("Model Results: \n")
    fp.write("ustarc (shear velocity) (m/s):")
    fp.write("%4.2f \n"%ustrc)
    fp.write("ztot (bed roughness):")
    fp.write("%.3e \n"%zotot)
    fp.write("thload (sediment load need) (m):")
    fp.write("%4.3f \n"%thload)
    fp.write("predictload (sediment load predicted) (m): ")
    fp.write("%4.3f \n"%predictedssload)
    fp.write("max. speed (m/s):")
    fp.write("%5.2f \n"%maximumspeed)
    fp.write("avg. speed (m/s):")
    fp.write("%5.2f \n"%avgspeed)
    fp.write("max. Froude:")
    fp.write("%4.2f \n"%MaxFroude)
    fp.write("avg. Froude:")
    fp.write('%4.2f \n'%AvgFroude)
    fp.write("Detail of structure: \n")
    fp.write("ssfrwanted (suspended size distribution wanted): \n")
    fp.write("%s \n"%ssfrwanted)
    fp.write("ssconverge (keeps track of whether amount of sediment in each size class is acceptable): \n")
    fp.write("%s \n"%ssconverge)
    fp.write("cconverge: \n")
    fp.write("%s \n"%cconverge)
    fp.write("ssfra(array to track of suspended sediment concentration in each size class): \n")
    fp.write("%s \n"%ssfra)
    fp.write("offa(array to track how far total suspended load is off from load required to create deposit): \n")
    fp.write("%s \n"%offa)
    fp.write("C(Suspended-sediment profile (calculated for each size class)): \n")
    fp.write("%s \n"%C)
    fp.write("G(Suspended load (calculated for each size class)): \n")
    fp.write("%s \n"%G)
    fp.write("ws(settling velocity (m/s)): \n")
    fp.write("%s \n"%ws)
    fp.write("SL: \n")
    fp.write("%s \n"%sl)
    fp.write("z(log space elevation values where concentrations and velocities are calculated): \n")
    fp.write("%s \n"%s)
    fp.write("zsl: \n")
    fp.write("%s \n"%zsl)
    fp.write("zmid: \n")
    fp.write("%s \n"%zmid)
    fp.write("ustrca(array to track how  u*c changes with iteration): \n")
    fp.write("%s \n"%ustrca)
    fp.write("S(excess shear stress): \n")
    fp.write("%s \n"%S)
    fp.write("spd(mean current speed calculated using current shear velocity (m/s)): \n")
    fp.write("%s \n"%spd)
    fp.write("results of the structure \n")
    fp.write("zoN(Nickaradse grain roughness): \n")
    fp.write("%s \n"%zoN)
    fp.write("tcrit(critical shear stress (N/m^2)): \n")
    fp.write("%s \n"%tcrit)
    fp.write("taub(bottom shear stress due to currents (N/m^2)): \n")
    fp.write("%s \n"%taub)
    fp.write("taustar(bottom shear stress normalized by critical shear stress): \n")
    fp.write("%s \n"%taustar)
    fp.write("phi(grain size in phi units for each size class): \n")
    fp.write("%s \n"%phi)
    fp.write("Ca(reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters): \n")
    fp.write("%s \n"%Ca)
    fp.write("Ki(integral of 1/eddie viscosity): \n")
    fp.write("%s \n"%Ki)
    fp.write("off(normalized difference between suspended load and load required to create deposit): \n")
    fp.write("%s \n"%off)
    fp.write("\n End")
    fp.write('General results \n')
    fp.write('Predictload (kg/m3) \n')
    fp.write('%f \n'%predictedssload)
    fp.write('maximum speed \n')
    fp.write('%f \n'%maximumspeed)
    fp.write('average speed \n')
    fp.write('%f \n'%avgspeed)
    fp.write('Max froude \n')
    fp.write('%f \n'%MaxFroude)
    fp.write('average froude \n')
    fp.write('%f \n'%AvgFroude)
    fp.close()
    return(maximumspeed,avgspeed,MaxFroude,AvgFroude)

'''
Created on 2013-1-25
========================================The TsuSedMod tsunami sediment transport model==============================
#=============================================The calculate model part==============================================
##==========================================Read data from file for model===========================================
##============================================Calculate water parameter=============================================
###=================================================calculate s=====================================================
###===============================calculate critical shear stress (N/m2) and settling velocity for each=============
###=================================================calculate zo====================================================
###==================================calculate shear velocity needed to suspend deposit=============================
#==============================================subfunction part=====================================================
##============================================Subfunction grading===================================================
###============================================Calculate part=======================================================
###===========================================plot figure part======================================================
####============================plot figure for Grain Statistics, 1st Moment for layers=============================
####===================================plot figure for Standard Deviation for layers================================
####=====================================plot figure for FW mean for layers=========================================
####=====================================plot figure for FW sorting for layers======================================
####=====================================plot figure for FW Skewness for layers=====================================
####======================================plot figure for FW Kurtosis for layers====================================
####==================================plot figure for Thickness vs. Time sediment hits bottom=======================
##=============================================subfunction sedstats=================================================
##=============================================subfunction readModelFile============================================
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
#===================================================================================================================
#===================================================================================================================
#=========================================The tsunami sediment model based on Bruce Jaffe===========================
#===================================================================================================================
#===================================================================================================================
#This model is based on the paper published in Sedimentary Geology
#A SIMPLE MODEL FOR CALCULATING TSUNAMI FLOW SPEED FROM TSUNAMI DEPOSITS
#               Bruce E. Jaffe, Guy Gelfenbuam
#===================================================================================================================
#Fuction:estimates water velocity from tsunami sediment deposit
#Assumes steady uniform flow
#Functions needed: phi2mm.py, sw_dens0.py, sw_smow.py, KinVisc.py, UstarCrit.py, zoWS.py,dietrichWs.py, Setvel1.py,
#                  tubesetvel.py, linearK.py, parabolicK.py, gelfK.py
#Inputs for this model:         th:deposit thickness (m)
#               	        	h:estimate of water depth (m)
#                       		sal:salinity (psu)
#               	        	wtemp:water temperature (deg. C)
#               	        	nclass:number of sed sizes to consider
#                       		Cb:total bed concentration (1-porosity), usually 0.65
#                       		fr:fractional portion of each grain size in bed
#                               phi:grain size in phi units for each size class
#                       		RhoS:sediment density (kg/m^3) for each size class
#                         		Dmn:mean grain diamter(mm);
#                               ustrc:shear velocity, initial guess input and then allowed to change
#                               ustrcAdjFac:Factor that sets how quickly ustrc is adjusted during iterations
#                               nit:number of iterations to run
#                               diststep:how often to iterate size classes versus concentration
#                               nz:number of vertical bins for calculations
#                        		vk:Van Karmen's constant
#                        		g0:the resuspension coefficient  *** standard value 1.4 x 10^-4 from Hill et al., 1989
#                         		concconvfactor:difference between modeled and observed concentration difference allowed
#                               sizeconvfactor:set the difference between modeled and observed size distribution allowed
#                               var:sets eddy viscosity shape; 1=linear, 2=parbolic,3=[Gelfenbaum]
#Used:                          ssfrwanted:suspended size distribution wanted
#                               ssconverge:keeps track of whether amount of sediment in each size class is acceptable
#                               ssoffa:array to track of how far off the suspended sediment concentration is from desired concentration
#                               ssfra:array to track of suspended sediment concentration in each size class
#                               fra:array to track of bed sediment concentration in each size class
#                               offa:array to track how far total suspended loasource_distribution.pyd is off from load required to create deposit
#                               ustrca:array to track how  u*c changes with iteration
#Parameter calculated by import function:
#                               d:grain diameter (set to m)
#                         		rhow:density of water (kg/m^3)sedstat
#                       		kinvisc:kinematic viscosity of water (m^2/s), usually about 10^-6
#                               UcritDmn:critical shear velocity for mean grain size
#                               zotot:Z naut total, from Nickaradse bed roughness and saltation bed roughness
#                         		ws:settling velocity (m/s)
#                        		ucrit:critical shear velocity for initiation of sediment transport
#                               z:log space elevation values where concentrations and velocities are calculated
#                          		K:eddie viscosity (m/s^2)
#                       		zo:Z naut, zero velocity intercept
#parameter calculated in this model:
#                               zoN- Nickaradse grain roughness
#      	                        ustrc:current shear velocity, u*c [m/s], initial value 0.1 m/s
#		                        s:ratio of density of sediment to density of fluid density for each size class
#   	                        thload:sediment load needed to get observed tsunami deposit thickness
#		                        tcrit:critical shear stress (N/m^2)
#		                        taub:bottom shear stress due to currents (N/m^2)
#		                        taustar:bottom shear stress normalized by critical shear stress
#		                        S:excess shear stress
#                               Ca:reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters
#		                        Ki:integral of 1/eddie viscosity
#                               C:Suspended-sediment profile (calculated for each size class)
#                               G:Suspended load (calculated for each size class)
#                               off:normalized difference between suspended load and load required to create deposit
#		                        spd:mean current speed calculated using current shear velocity (m/s)
#                               froude:Froude number calculated using maximum speed
#Options:
#                               'infile':specifies filename of model input file.
#                               'outfile':writes results to comma-separated (.csv) file. If the specified output file already exists,
#                               results will be appended to current file.
#                               'grading':plots results of grading function . The default output depth of this function is 0.01 m.  To
#                                modify the defaut, the user can specify an additional arguement of a vector of depths or a single value.
#Usage examples:
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5;
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('infile','C:\exampleFile.xls');
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',0.01);
#                                   where 0.01 specifies set interval (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',[0.01,0.03,0.04,0.8]);
#                                   where [ ] specifies variable intervals (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('outfile','modelResults.csv');
#Version information:
#written in March 12, 2001 by Bruce Jaffe, based on code written by Chris Sherwood, documented in July 3,2001
#Significant code clean-ups and user friendly input added by Andrew Stevens from March to October, 2007
#modfied on 11/29/07 by Bruce Jaffe, modification of Andrew Stevens' cleaned-up code including: allowing single size class runs
#modified 4/29/08 to output weight percent (Bruce Jaffe) and to fix grading function interval error (Mark Buckley)
#made in python 1/25/2013 by Hui Tang (tanghui@vt.edu)
#==============================================================================================================================================
#==============================================================The caculate model part=========================================================
from pylab import *
# from function import *
# from sedstats import *
# from numpy import *
# from grading import *
# from ReadCSV import *
import csv
import os
def Tsusedmod2(name):
##============================================================Read data from file for model====================================================
    print("\n ====================================running==================================== \n")
    #read in model inputs from excel file
    separator=':'
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='number_of_iterations': nit=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='number_of_vertical_bins': nz=float(v)
                elif n=='resuspension_coefficient': g0=float(v)
                elif n=='concentration_convergence_factor': concconvfactor=float(v)
                elif n=='size_convergence_factor': sizeconvfactor=float(v)
                elif n=='shear_velocity': ustrc=float(v)
                elif n=='bed_roughness': zotot=float(v)
                elif n=='settling_velocity_type': setType=float(v)
                elif n=='eddy_viscosity_shape': var=float(v)
                #site info
                elif n=='site_description': t=str(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='bed_concentration': Cb=float(v)
                #data
                elif n=='Number_sample':num=int(v)
                elif n=='Number_grainsize':nclass=int(v)
                elif n=='data_dimension':l=int(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
    filename=name
    data=readCSV1(filename,separator=';')
    phi=data[:,0]
    fr1=data[:,1]
    fr=zeros(shape=(len(fr1),1))
    th=data[0,2]
    h=data[0,3]
    x=data[0,4]
    #end data input
    nclass=len(fr1) #number of sed sizes to consider
    fr=fr1.reshape(len(fr1),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    ssfrwanted=fr #susp sed size distribution wanted, keep this for comparison with model results
    ssconverge=ones(shape=(nclass,1))#initiate ss converge to 1 for all size classes (does  converge)
    for i in range(len(fr)):
        if(fr[i]!=0):
            ssconverge[i]=0 #reset ss converge values to 0 (does not converge) for all size classes with sediments
    sizeconverge=sizeconvfactor*ones(shape=(nclass,1))#set percent difference in size distribution allowed, same value for all sizes
    D=phi2mm(phi) #convert grain size to mm for settling velocity calculation using Dietrich's formulae
    d=0.001*D #convert classes to m for calculating settling velocity using other forumlae and for UstarC calculation
    if(nclass==1):
        dMean_phi=phi[0]
    else:
        sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
        dMean_phi=sedstat[0]
        m3a=sedstat[1]
        m4a=sedstat[2]
        stdeva=sedstat[3]
        fwmedian1=sedstat[4]
        fwsort1=sedstat[5]
        fwmean1=sedstat[6]
        fwskew1=sedstat[7]
        fwkurt1=sedstat[8]
    dMean_mm=phi2mm(dMean_phi) # convertmean grain size in phi to mm
    Dmn=0.001*dMean_mm  # convert mean grain size to m for calculating bed roughess
    rhos=Rhos*ones(shape=(nclass,1)) #set sediment density, same value for all sizes
#==========================================================Calculate water parameter=============================================================
    rhow=sw_dens0(sal,wtemp)  #calculate the density of seawater [kg/m^3]
    kinvisc=KinVisc(wtemp)    #calculate the kinematic viscosity of the water
###===============================================================calculate s=====================================================================
    s=rhos/rhow
###===============================calculate critical shear stress (N/m2) and settling velocity for each===========================================
    #size class
    ws=zeros(nclass)
    Ucrit=zeros(nclass)
    tcrit=zeros(nclass)
    if(setType==1.0):
        for i in range(nclass):
            ws[i]=-tubesetvel(phi[i],rhos,wtemp,sal)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==2.0):
        for i in range(nclass):
            ws[i]=-dietrichWs(wtemp,d[i],rhos[i],rhow,0.005,0.7,3.5)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==3.0):
        for i in range(nclass):
            ws[i]=Setvel1(d[i],kinvisc,s[i])
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    thload=th*Cb #suspended sediment load to make deposit of observed thickness, thickness times Cb = 1-porosity
###============================================================calculate zo=======================================================================
    zoN=Dmn/30.0  #Nickaradse bed roughness
    #guess ustrc to make deposit, need this to calculate zos, bed roughness from moving flat bed, saltation
    #Note:Make the initial ustrc an input, also make ustrcAdjFac an input and delelet lines 228 and 229
    ustrcAdjFac=0.01
    UcritDmn=UstarCrit(Dmn,kinvisc,2650/rhow) #need critical shear velocity for mean grain size to calculate bed roughness
    if(zotot!=1):
        zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #zo total (Nickaradse bed roughness and Wiberg/Smith bed roughness from saltation,
###===============================================calculate shear velocity needed to suspend deposit===============================================
    #Begin loop to determine shear velocity needed to suspend deposit
    #This loops used to assume that the bed sediment disribution is the same as the suspended sediment distribution
    #This gives reasonable results when the size of the bed material is number of iterations to run
    diststep=1 #adjusts sediment grain size distribution every 2nd iteration
    ssoffa=zeros(shape=(nit/diststep,nclass)) # initiate array tracking ss deviation from desired concentration for iteration=1:nit
    ssfra=ones(shape=(nit/diststep,nclass)) #initiate array tracking ss concentrations for iteration=1:nit
    fra=zeros(shape=(nit/diststep,nclass)) #initiate array tracking bed sediment concentrations for iteration=1:nit
    offa=zeros(nit)#inititate array to track how far total suspended load is off from desired suspended load
    ustrca=zeros(nit)#inititate array to track how  u*c changes with iteration                                            #ref. Wiberg and Rubin, 1989 JGR
    cconverge=0
    for iteration in range(int(nit)):
        fr=fr/sum(fr) #make sure fractions add up to 1
        if(all(ssconverge) and cconverge):
            break
        taub=rhow*ustrc**2 # bottom stress
        taustar=taub/tcrit
        S=taustar-1     #normalized excess shear stress
        for i in range(len(S)):#returns zero when taub<=tcrit
            if S[i]<0:
                S[i]=0
            else:s[i]=s[i]
        S=S.reshape(len(S),1)
        Ca=zeros(shape=(len(s),1))
        for i in range(len(S)):
            Ca[i]=g0*Cb*fr[i]*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        # resuspension coefficient(g0) for flat moveable bed Madsen 94a "Susp. Sed. Transport on the inner shelf waters during extreme storms"
        # ICCE pg 1849 Madsen, Chisholm, Wright based on Wikramanayake 92 PhD thesis first shot, might want to change this using mean current litterature
        # log-spaced z grid
        z=logspace(log10(zotot),log10(h),num=nz,endpoint=True,base=10.0)
        z=z.reshape(len(z),1)
#*************************************************elevation of suspended load, SL******************************************
        def diff(z):
            diff=zeros(shape=(len(z)-1,1))
            for i in range(0,len(z)-1):
                diff[i]=z[i+1]-z[i]
            return (diff)
        zsl=z[0:nz-1]+diff(z)
        diff1=diff(z)
#**************************************elevation of speed calculated from eddy viscosity***********************************
        zmid=(z[0:nz-1]+z[1:nz])/2
        zmid1=zmid.reshape(len(zmid))
#*******************************************calculate eddy viscosity profile***********************************************
        #var=1 is linear eddy viscosity profile, var=2 is parabolic profile, var=3 is gelf profile
        if var==1.0:
            K=linearK(z,h,ustrc)
        elif var==2.0:
            K=parabolicK(z,h,ustrc)
        elif var==3.0:
            K=gelfK(z,h,ustrc)
        #K=(var==1)*linearK(z,h,ustrc)+(var==2)*parabolicK(z,h,ustrc)+(var==3)*gelfK(z,h,ustrc)
#*******************************************Calculate current and eddie viscosity profiles*********************************
        #Integral of 1/K
        K=K.reshape(len(K),1)
        Ki=cumsum(2.0*diff(z)/(K[0:nz-1]+K[1:nz]))
#*********************************calculate the suspended-sediment profile and sediment load for each class****************
        C=zeros(shape=(nz,nclass))
        sl=zeros(shape=(nz-1,nclass))
        sl1=zeros(shape=(nz,nclass))
        G=zeros(nclass)
        for i in range(nclass):
            #Suspended-sediment profile
            C[:len(Ca[0]),i]=Ca[i][0]
            C[len(Ca[0]):,i]=Ca[i][0]*exp(ws[i]*K[i])
            for j in range(int(nz-2)):
                sl[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
                sl1[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
#******************************************Depth-intergral of sediment profile*********************************************
            G[i] = sum(sl[:,i])
#************************************************** set cconverge to 0*****************************************************
        cconverge=0
        off=sum(G)/thload-1
        offa[iteration]=off
        if(abs(off)*100<concconvfactor):
            cconverge=1
        if(all(ssconverge) and cconverge): #if both conc and size are good enough, get outta here
            break
        ustrca[iteration]=ustrc
        ustrc=ustrc*(1-ustrcAdjFac*off) #change ustrc based on how far off suspended load is from desired value
        if(zotot!=1):
            zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #recalculate bed roughness
#*******check for convergence for size of suspended sediment and adjust bed sediment distribution every diststep***********
        def fix(n):
            if(n>0):
                if(round(n)>n):
                    fix=round(n)-1
                else:fix=round(n)
            if(n<0):
                if(round(n)<n):
                    fix=round(n)+1
            return(n)
        doss=diststep*fix(iteration/diststep)
        if(doss==iteration): #redo distribution every diststep steps
            ssfr=G/sum(G) #calculate susp sed size distribution for model
            ssfra[doss/diststep,:]=ssfr
            fra[doss/diststep,:]=fr.reshape(1,len(fr))
            ssfr=ssfr.reshape(len(ssfr),1)
            #keep track of whether model is converging ** array indexes
            #different from iteration
            for k in range(nclass):
                if(ssfrwanted[k]!=0): # only adjust for classes with material in them
                    ssoff=ssfr[k][0]/ssfrwanted[k]-1
                    ssoffa[doss/diststep,k]=ssoff #change bed grain size distribution based on suspended load distribution
                    fr[k]=fr[k]*(1-0.01*ssoff)
                    if (abs(ssoff)*100<sizeconverge[k]): #check for convergence for each size class
                        ssconverge[k]=1.0 #set to 1 if OK
                if(all(ssconverge) and cconverge):
                    break #if both conc and size are good enough, get outta here
    spd=zeros(len(z)-1)
    diff1=diff(z)
    diff1=diff1.reshape(len(diff1),)
    spd[:]=cumsum(diff1*ustrc*ustrc/gelfK(zmid,h,ustrc))#speed calculated using eddy viscosity profile integration
    if(iteration==nit):
        print('Model may not have converged, use extreme caution with these resuls')
    #stop clock
    predictedssload=sum(G[:])
    maximumspeed=max(spd)
    spd1=0.5*(spd[0:len(spd)-1]+spd[1:len(spd)])
    spd1=spd1.reshape(len(spd1),1)
    diff2=zeros(shape=(len(diff1)-1,1))
    diff2=diff1[1:len(diff1)]
    diff2=diff2.reshape(len(diff2),1)
    avgspeed=sum((spd1) * diff2)/h
    MaxFroude=maximumspeed/sqrt(9.8*h)
    AvgFroude=avgspeed/sqrt(9.8*h)
    if (var==1): ed='Linear'
    elif (var==2): ed='Parabolic'
    elif (var==3): ed='Gelfenbaum'
    speed=zeros(shape=(100,2))
    zsl=zsl.reshape(100,)
    speed[:,0]=zsl
    speed[:,1]=spd
    return(speed)
'''
Created on 2013-5-21
=============================The original tsunami sediment concentration model===========================
#=========================================Input parameters===============================================
#=====================================The calculate model part===========================================
##==============================Caculate critical velocity and critical shear stress=====================
##=================================Caculate the concentration of each grain size=========================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#==============================The original tsunami sediment concentration model=========================
#========================================================================================================
#========================================================================================================
#This model is based on the paper published Proceedings of the 24th International Conferrence on Coastal
#Engeering, 1993
#SUSPENDED SEDIMENT SUSPENSION ON THE INNER SHELF DURING EXTREME STROMS
#              Madsen, O.S, Chisholm, T.A and Wright, L.D
#Hui Tang 05/21/2013
#Input:             Se-Grain size in phi
#                   fr-Percentage of each grain size %
#                   Vi-Smallest possible average velocity m/s
#                   H-water depth at still water level m
#                   Cb-Bed concentration M^3/m^3
#                   rho-Density of sea water
#Used constant:     k-Karman's constant 0.41
#                   rho-Density of fluid: 1050 kg/m^3
#                   rhos-Density of sediments: 2770 kg/m^3
#                   g0-Resuspension coefficient 0.0004
#Function used:     UstarCrit
#Output:            Se-Grain size in phi
#                   C0-initial concentration for sediments source m^3/m^3
#========================================================================================================
#==========================================Input parameters==============================================
#This part is designed to input the parameter for this model
from pylab import *
import os
# from output2CSV import *
# from function import *
def inconcentration(se1,fr,Vi,H,Cb,rho):
    with open('../input/parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                if n=='resuspension_coefficient': g0=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
		    	continue
#========================================================================================================
#======================Deal with the data from original grainsize distribution===========================
#This part is designed to deal with the data from orignal grainsize distribution file, in order to use in
# caculating concentration.
    se=2**(-se1)*0.001
    nc=len(se)
#========================================================================================================
#=====================================The calculate model part===========================================
#This part is designed for calculate the original sediment concentration base on the grain size
##==============================Caculate critical velocity and critical shear stress=====================
    Ucrit=zeros(nc)#Critical shear velocity for initiation of sediment transport
    tcrit=zeros(nc)#Critical shear stress (N/m^2)
    C=zeros(nc)# Original concentration for each grain size
    for i in range(nc):
        Ucrit[i]=UstarCrit(se[i],KinVisc(wtemp),Rhos/rho)# caculate critical shear velocity for initiation of sediment transport
        tcrit[i]=rho*Ucrit[i]**2# Caculate critical shear stress (N/m^2)
    taub=(Vi*vk/(log(30*H/(2**(-1)-2**(-9))-(1-(H*(2**(-1)-2**(-9)))))))**2*rho # bottom stress(?)
    taustar=taub/tcrit
    S=taustar-1     #normalized excess shear stress
    for i in range(nc):#returns zero when taub<=tcrit
        if(S[i]<0):
            S[i]=0
    for i in range(len(S)):
        C[i]=fr[i]*g0*Cb*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        C0=Rhos*C
    return(C0)
