'''
Fennel et al (2006) Nitrogen cycling in the Middle Atlantic Bight: Results from a 
three-dimensional model and implications for the North Atlantic nitrogen budget.
GLOBAL BIOGEOCHEMICAL CYCLES, VOL. 20, GB3007, doi:10.1029/2005GB002456

'''


def load_defaults():
    '''
    This function creates a dictionaries called "par" and "InitCond"
    and pre-loads them with all the default 
    parameters and initial conditions, respectively.
    Also outputs days and dt
    '''
    # Framework
    days = 365 * 3 # Three year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    par['mu0']   = 0.69  
    par['kNO3']  = 0.5    
    par['kNH4']  = 0.5  
    par['alpha'] = 0.125  
    par['gmax']  = 0.6  #Original 0.6
    par['kP']    = 2
    par['mP']    = 0.15    
    par['tau']   = 0.005 
    par['thetaMax']= 0.053
    par['beta']  = 0.75 
    par['lBM']   = 0.01   
    par['lE']    = 0.01
    par['mZ']    = 0.025
    par['rSD']   = 0.3 # Original 0.03
    par['rLD']   = 0.1 # # Original 0.01
    par['nmax']  = 0.05
    par['kI']    = 0.1
    par['I0']    = 0.0095
    par['wP']    = 0.1
    par['wS']    = 0.1
    par['wL']    = 1.
    par['AE_P']  = 0.9  
    par['AE_D']  = 0.2    
    par['AE_Z']  = 0.3  
    par['Bpub']  = 0.43  
    par['Fmax_ref']= 0.03#0.025
    par['GT']    = 0.44
    par['KTempH']= 0.1    
    par['KTempL']= 0.5 
    par['KSaltL']= 0.25
    par['KOxyL'] = 0.02 
    par['KFood'] = 1.    
    par['KRE']   = 0.86
    par['OxyL']  = 17.5
    par['Rm']    = 0.002
    par['SaltL'] = 10.
    par['TempH'] = 25.
    par['TempL'] = -4.
    par['beta']  = 0.12
    par['epsilonP'] = 1.
    par['epsilonD'] = 0.5
    par['epsilonZ'] = 0.3
    par['lamda_nat'] = 0.0 # 0.00137
    par['lamda_harvest'] = 0.0 # 0.001
    
    
    # Physical characteristics of embayment
    par['chi'] = 0.01
    par['X'] = 2000#2000 # Basin length
    par['Y'] = 200#200 # Basin width
    par['Z'] = 10 # Basin depth
    par['V'] = par['X'] * par['Y'] * par['Z']
    par['uwind'] = 0.5
    par['vwind'] = 0.5
    
    # Initial conditions
    InitCond = {}
    InitCond['Phy'] = 0.5
    InitCond['Zoo'] = 0.2
    InitCond['SDet'] = 0.5 
    InitCond['LDet'] = 0.3
    InitCond['NH4'] = 0.1
    InitCond['NO3'] = 5.
    InitCond['Oxy'] = 340. #Oxygen 
    InitCond['Soma'] = 3.0
    InitCond['Gonad'] = 0.0
    InitCond['conc_muss'] = 20.
    InitCond['n_muss'] = InitCond['conc_muss'] * par['V']
    
    return days, dt, par, InitCond
    



    
def run_model(days,dt,InitCond,par,forc):
    '''
    This is your model. Do a brief description.

    '''
    # Import libraries
    import numpy as np
    
    print 'Starting model run with ' + str(InitCond['conc_muss']) + ' mussels/m3 ...'
    
    # Make sure n_muss is correctly estimated
    InitCond['n_muss'] = InitCond['conc_muss'] * par['V']
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    # Create zero-vectors
    Phy = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Zoo = np.zeros((NoSTEPS,),float) # same as above
    SDet = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    LDet = np.zeros((NoSTEPS,),float) # same as above
    NH4 = np.zeros((NoSTEPS,),float) # same as above
    NO3 = np.zeros((NoSTEPS,),float) # same as above
    Oxy = np.zeros((NoSTEPS,),float) # same as above

    
    mu = np.zeros((NoSTEPS,),float) # same as above
    f_I = np.zeros((NoSTEPS,),float) # same as above
    L_NO3 = np.zeros((NoSTEPS,),float) # same as above
    L_NH4 = np.zeros((NoSTEPS,),float) # same as above
    airwater_O2_flux = np.zeros((NoSTEPS,),float) # same as above
    TotN = np.zeros((NoSTEPS,),float) # same as above
    
    # Mussels
    Soma = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Gonad = np.zeros((NoSTEPS,),float) # same as above
    B = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    B_conc = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    L_Temp = np.zeros((NoSTEPS,),float) # same as above
    L_Salt = np.zeros((NoSTEPS,),float) # same as above
    L_Oxy = np.zeros((NoSTEPS,),float) # same as above
    L_Food = np.zeros((NoSTEPS,),float) # same as above
    F = np.zeros((NoSTEPS,),float) # same as above
    A = np.zeros((NoSTEPS,),float) # same as above
    R = np.zeros((NoSTEPS,),float) # same as above
    RE = np.zeros((NoSTEPS,),float) # same as above
    Spawning = np.zeros((NoSTEPS,),float) # same as above
    n_muss = np.zeros((NoSTEPS,),float) # same as above
    CumulativeHarvest = np.zeros((NoSTEPS,),float) # same as above

    
    
    # Initializing with initial conditions
    Phy[0] = InitCond['Phy']
    Zoo[0] = InitCond['Zoo']
    SDet[0] = InitCond['SDet']
    LDet[0] = InitCond['LDet']
    NH4[0] = InitCond['NH4']
    NO3[0] = InitCond['NO3']
    Oxy[0] = InitCond['Oxy']
    
    # Mussels
    Soma[0] = InitCond['Soma']
    Gonad[0] = InitCond['Soma']
    B[0] = InitCond['Soma'] + InitCond['Gonad']
    Spawning[0] = 0.
    n_muss[0] = InitCond['n_muss']
    B_conc[0] = B[0] * n_muss[0] / par['V']
    
    # Forcing
    Salt = forc['Salt']
    Temp = forc['Temp']
    I = forc['I']
    Phy0 = forc['Phy']
    Zoo0 = forc['Zoo']
    SDet0 = forc['SDet']
    LDet0 = forc['LDet']
    NH40 = forc['NH4']
    NO30 = forc['NO3']
    Oxy0 = forc['Oxy']


    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        muMax = par['mu0'] * (1.066 ** Temp[t]) # text
        
        f_I[t] = (par['alpha']*I[t])/(np.sqrt(muMax**2+((par['alpha']**2)*(I[t]**2)))) #Eq5
        
        L_NO3[t] = (NO3[t]/(par['kNO3']+NO3[t])) * (1/(1+(NH4[t]/par['kNH4']))) #Eq3
        
        L_NH4[t] = NH4[t]/(par['kNH4']+NH4[t]) # Eq 4
    
        mu[t] =muMax * f_I[t] * (L_NO3[t] + L_NH4[t]) # Eq2
    
        g = par['gmax'] * ((Phy[t]**2)/(par['kP']+(Phy[t]**2)))
        
        n = par['nmax'] * (1 - max(0,(I[t]-par['I0'])/(par['kI']+I[t]-par['I0'])))
        
        n_O2 = (Oxy[t]/(3.+Oxy[t]))
    
        dPhydt = (mu[t] * Phy[t]) - \
                 (g  * Zoo[t]) - \
                 (par['mP'] * Phy[t]) - \
                 (par['tau']*(SDet[t]+Phy[t])*Phy[t]) - \
                 (par['wP']*Phy[t]/par['Z']) # Eq1
                 
        dZoodt = (g * par['beta'] * Zoo[t]) - \
                 (par['lBM']*Zoo[t]) - \
                 (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                 (par['mZ']*(Zoo[t]**2))#Eq10
        
        dSDetdt = (g * (1-par['beta']) * Zoo[t]) + \
                  (par['mZ']*(Zoo[t]**2)) + \
                  (par['mP'] * Phy[t]) - \
                  (par['tau']*(SDet[t]+Phy[t])*SDet[t]) - \
                  (par['rSD']*SDet[t]) - \
                  (par['wS']*SDet[t]/par['Z'])
                  
        dLDetdt = (par['tau']*((SDet[t]+Phy[t])**2)) - \
                  (par['rLD']*LDet[t]) - \
                  (par['wL']*LDet[t]/par['Z'])
                  
        dNO3dt = -(muMax * f_I[t] * L_NO3[t] * Phy[t]) + \
                  (n * n_O2 * NH4[t])
                 
        dNH4dt = -(muMax * f_I[t] * L_NH4[t] * Phy[t]) - \
                  (n * n_O2 * NH4[t]) + \
                  (par['lBM'] * Zoo[t]) + \
                  (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) + \
                  (par['rSD']*SDet[t]) + \
                  (par['rLD']*LDet[t]) + \
                  (par['wP']*Phy[t]/par['Z']) + \
                  (par['wS']*SDet[t]/par['Z']) + \
                  (par['wL']*LDet[t]/par['Z'])


        # MUSSELS -------------------------------------------------------------
        # Calculate Temperature Limitation
        L_Temp[t] = min(max(0.,1.-np.exp(-par['KTempL']*(Temp[t]-par['TempL']))), \
                     max(0.,1.+((1.-np.exp(par['KTempH']*Temp[t]))/(np.exp(par['KTempH']*par['TempH'])-1.))))
        
        # Calculate Salinity Limitation
        L_Salt[t] = max(0.,1.-np.exp(-par['KSaltL']*(Salt[t]-par['SaltL'])))
        
        # Calculate Oxygen Limitation
        L_Oxy[t] = max(0.,1.-np.exp(-par['KOxyL']*(Oxy[t]-par['OxyL'])))
        
        # Calculate Oxygen Limitation
        L_Food[t] = (Phy[t]+Zoo[t]+SDet[t])/(par['KFood']+Phy[t]+Zoo[t]+SDet[t])
        
        # Calculate Filtration rate
        Fmax  = par['Fmax_ref']*(B[t]**(2./3.))
        
        F[t] = Fmax * L_Temp[t] * L_Salt[t] * L_Oxy[t] * L_Food[t]
        
        A[t] = F[t] * ((par['epsilonP']*par['AE_P']*Phy[t])+ \
                       (par['epsilonZ']*par['AE_Z']*Zoo[t])+ \
                       (par['epsilonD']*par['AE_D']*SDet[t]))
        
        R[t] = (par['Rm']*B[t]) + (par['beta']*A[t])
        
        RE[t] = max(0., (B[t]-par['Bpub'])/(par['KRE'] + B[t] - (2.*par['Bpub'])))
        
        # Spawning
        if n_muss[t] == 0.: 
            Spawning[t] = 0.
            dGonaddt = 0
            dSomadt = 0.
        elif Gonad[t]/B[t] < par['GT']:
            Spawning[t] = 0.
            dGonaddt = (A[t]-R[t]) * RE[t]
            dSomadt =  (A[t]-R[t]) * (1.-RE[t])
            offset = Gonad[t] + (dGonaddt * dt)
            if offset < 0 : # If Gonad is going to be negative... don't apply dynamic allocation
                dGonaddt = 0.
                dSomadt = A[t]-R[t]
        elif Gonad[t]/B[t] >= par['GT']:         
            Spawning[t] = Gonad[t]
            dGonaddt = 0.
            dSomadt = A[t]-R[t]

        #Feedback to NPZD2 model
        # Faeces and Pseudofaeces
        Fae = F[t] * ((par['epsilonP']*(1-par['AE_P'])*Phy[t])+ \
                     (par['epsilonZ']*(1-par['AE_Z'])*Zoo[t])+ \
                     (par['epsilonD']*(1-par['AE_D'])*SDet[t]))
                     
        dLDetdt = dLDetdt + (Fae*(n_muss[t]/par['V']))    
                  
        # Remove eaten Phy, Zoo and SDet from water-column
        dPhydt =  dPhydt-((F[t] *par['epsilonP']*Phy[t])*(n_muss[t]/par['V']))
        dZoodt =  dZoodt-((F[t] *par['epsilonZ']*Zoo[t])*(n_muss[t]/par['V']))
        dSDetdt = dSDetdt -((F[t] *par['epsilonD']*SDet[t])*(n_muss[t]/par['V']))

        # Mussel population
        Harvest = par['lamda_harvest'] * n_muss[t]
        NatMortality = par['lamda_nat'] * n_muss[t]
        dn_mussdt = -NatMortality - Harvest
   
        # Excretion into Ammonia
        dNH4dt = dNH4dt + ((R[t]*n_muss[t]/par['V']) + ((NatMortality*B[t])/par['V']))       
      
        
        
        # Oxygen sub-model =========================================
        
        # Parameters
        OA0 = 2.00907       # Oxygen
        OA1 = 3.22014       # saturation
        OA2 = 4.05010       # coefficients
        OA3 = 4.94457
        OA4 =-0.256847
        OA5 = 3.88767
        OB0 =-0.00624523
        OB1 =-0.00737614
        OB2 =-0.0103410
        OB3 =-0.00817083
        OC0 =-0.000000488682
        rOxNO3= 8.625       # 138/16
        rOxNH4= 6.625       # 106/16
        l2mol = 1000.0/22.9316 # liter to mol
        
        #-----------------------------------------------------------------------
        #  Surface O2 gas exchange.
        #-----------------------------------------------------------------------
        
        #  Compute surface O2 gas exchange.
        cff2=0.31*(24.0/100.0)
        
        #  Compute O2 transfer velocity : u10squared (u10 in m/s)
        u10squ=(par['uwind']*par['uwind'])+(par['vwind']*par['vwind'])
        
        # Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
        SchmidtN_Ox=1953.4-Temp[t]*(128.0-Temp[t]*(3.9918-Temp[t]*0.050091))
        cff3=cff2*u10squ*np.sqrt(660.0/SchmidtN_Ox)        
        
        #  Calculate O2 saturation concentration using Garcia and Gordon
        #  L&O (1992) formula, (EXP(AA) is in ml/l).        
        TS=np.log((298.15-Temp[t])/(273.15+Temp[t]))        
        
        AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+ \
           Salt[t]*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+ \
           OC0*Salt[t]*Salt[t]
        
        # Convert from ml/l to mmol/m3.
        O2satu=l2mol*np.exp(AA)        
        
        # Add in O2 gas exchange.
        O2_Flux = cff3*(O2satu-Oxy[t])
        
        airwater_O2_flux[t] = O2_Flux * (1./par['Z'])
        
        dOxydt = airwater_O2_flux[t]
        
        
        # Production via Photosynthesys
        dOxydt = dOxydt + (muMax * f_I[t] * L_NO3[t] * Phy[t] * rOxNO3) # New production
        dOxydt = dOxydt + (muMax * f_I[t] * L_NH4[t] * Phy[t] * rOxNH4) # Regenerated production
        
        # Respiration
        dOxydt = dOxydt - (((par['lBM']*Zoo[t]) - \
                           (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                           (par['mZ']*(Zoo[t]**2))) * rOxNH4) # Zooplankton
 
        dOxydt = dOxydt - ((n * n_O2 * NH4[t])* rOxNH4 * 2) # Nitrification 
 
        dOxydt = dOxydt - (((par['rSD']*SDet[t])+(par['rLD']*LDet[t])) * rOxNH4) #S and L Detritus remineralization
        
        dOxydt = dOxydt - (((par['wS']*SDet[t]/par['Z']) + \
                          (par['wL']*LDet[t]/par['Z'])) * rOxNH4) #S and L Detritus remineralization in sediments
        
        dOxydt = dOxydt - ((((R[t]*n_muss[t])/par['V']) + (NatMortality/par['V'])) * rOxNH4)  #Mussels

       
        
        
        # Physical Model ======================================================
        dPhydt = dPhydt + (par['chi'] * (Phy0[t] - Phy[t]))
        dZoodt = dZoodt + (par['chi'] * (Zoo0[t] - Zoo[t]))
        dNH4dt = dNH4dt + (par['chi'] * (NH40[t] - NH4[t]))
        dNO3dt = dNO3dt + (par['chi'] * (NO30[t] - NO3[t]))
        dSDetdt = dSDetdt + (par['chi'] * (SDet0[t] - SDet[t]))
        dLDetdt = dLDetdt + (par['chi'] * (LDet0[t] - LDet[t]))
        dOxydt = dOxydt + (par['chi'] * (Oxy0[t] - Oxy[t]))
                
        

        # Update and step ----------------------------------------------------
        Phy[t+1]  = Phy[t]  + (dPhydt * dt)
        Zoo[t+1]  = Zoo[t]  + (dZoodt * dt) + ((Spawning[t]*n_muss[t])/par['V'])
        SDet[t+1] = SDet[t] + (dSDetdt * dt)
        LDet[t+1] = LDet[t] + (dLDetdt * dt)
        Oxy[t+1]  = max(0,Oxy[t] +  (dOxydt * dt))
        NH4[t+1]  = NH4[t]  + (dNH4dt * dt)
        NO3[t+1]  = NO3[t] +  (dNO3dt * dt)
        if NO3[t+1] <= 0.001:
            offset = NO3[t+1]
            NH4[t+1] = NH4[t+1] + offset
            NO3[t+1] = NO3[t+1] - offset
        # Mussels
        Soma[t+1] = Soma[t] + (dSomadt * dt)
        Gonad[t+1] = Gonad[t] + (dGonaddt * dt) - Spawning[t]
        B[t+1] = Soma[t+1] + Gonad[t+1]
        n_muss[t+1] = max(0,n_muss[t] + (dn_mussdt * dt))
        CumulativeHarvest[t+1] = CumulativeHarvest[t] + Harvest
        B_conc[t+1] =  B[t+1] * n_muss[t+1] / par['V']

        # Estimate Total Nitrogen
        TotN[t+1] = Phy[t+1] + Zoo[t+1] + SDet[t+1] + LDet[t+1] + NH4[t+1] + NO3[t+1] + ((B[t+1]*n_muss[t+1])/par['V'])
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['par'] = par
    output['InitCond'] = InitCond
    output['time'] = time
    output['Phy'] = Phy
    output['Zoo'] = Zoo
    output['SDet'] = SDet
    output['LDet'] = LDet
    output['NH4'] = NH4
    output['NO3'] = NO3
    output['Oxy'] = Oxy
    output['I'] = I
    output['mu'] = mu
    output['f_I'] = f_I
    output['L_NO3'] = L_NO3
    output['L_NH4'] = L_NH4
    output['TotN'] = TotN
    output['airwater_O2_flux'] = airwater_O2_flux
    output['Soma'] = Soma
    output['Gonad'] = Gonad
    output['B'] = B
    output['n_muss'] = n_muss
    output['Spawning'] = Spawning
    output['CumulativeHarvest'] = CumulativeHarvest
    output['F'] = F
    output['L_Temp'] = L_Temp
    output['L_Salt'] = L_Salt
    output['L_Oxy'] = L_Oxy
    output['L_Food'] = L_Food
    output['B_conc'] = B_conc

    print "Model run: DONE!!!"
    return output


    
def plot_model(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax, ax2, ax3) = plt.subplots(3,1,figsize=(13,13))
    ax.plot(output['time']/365,output['Phy'],'g-')
    ax.plot(output['time']/365,output['Zoo'],'r-')
    ax.plot(output['time']/365,output['SDet'],'k-')
    ax.plot(output['time']/365,output['LDet'],'k-.')
    ax.plot(output['time']/365,output['NH4'],'m-')
    ax.plot(output['time']/365,output['NO3'],'c-')
    ax.plot(output['time']/365,output['B_conc'],'r.')
    ax.set_ylabel('Nitrogen \n (mmol N m$^{-3}$)')
    ax.set_title('Ecosystem Model - Plankton Ecosystem')
    ax.legend(['Phy','Zoo','SDet','LDet','NH4','NO3','B_conc'])
    
    ax2.plot(output['time']/365,output['Oxy'],'b-')
    ax2.set_ylabel('Oxygen \n (mmol O2 m$^{-3}$)')
    ax2.legend(['f_I','Mu','L_NO3','L_NH4'])
    
    ax3.plot(output['time']/365,output['f_I'],'r-')
    ax3.plot(output['time']/365,output['mu'],'g-')
    ax3.plot(output['time']/365,output['L_NO3'],'b-')
    ax3.plot(output['time']/365,output['L_NH4'],'k-')
    ax3.set_ylabel('Plankton Diagnostics \n (dimensionless)')
    ax3.legend(['f_I','Mu','L_NO3','L_NH4'])


    
    fig2, (ax, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(13,13))
    ax.plot(output['time']/365,output['B'],'r.')
    ax.plot(output['time']/365,output['Gonad'],'b-')
    ax.plot(output['time']/365,output['Soma'],'k-')
    ax.legend(['B','Gonad','Soma'])
    ax.set_title('Ecosystem Model - Mussels')
    
    ax2.plot(output['time']/365,output['n_muss'],'g-')
#    ax2.plot(output['time']/365,output['CumulativeHarvest'],'r-')
    ax2.set_ylabel('Filtration rate \n (L ind$^{-1}$ h$^{-1}$)')
    ax2.legend(['Total Number of \n mussels in bay'])

    ax3.plot(output['time']/365,output['F']*1000/24,'g-')
    ax3.set_ylabel('Filtration rate \n (L ind$^{-1}$ h$^{-1}$)')
    ax3.legend(['F'])

    ax4.plot(output['time']/365,output['L_Temp'],'b-')
    ax4.plot(output['time']/365,output['L_Salt'],'m-')
    ax4.plot(output['time']/365,output['L_Oxy'],'k-')
    ax4.plot(output['time']/365,output['L_Food'],'r-')
    ax4.legend(['L_Temp','L_Salt','L_Oxy','L_Food'])
    ax4.set_ylabel('Mussel Diagnostics \n (dimensionless)')
    ax4.set_xlabel('Time (years)')

    plt.show()
    return
    
def plot_totalNitrogen(output):
    import matplotlib.pyplot as plt
    fig3, (ax) = plt.subplots(1,1)
    ax.plot(output['time']/365,output['TotN'],'y.')
    ax.legend(['TotN'])
    ax.set_title('Ecosystem Model - Total Nitrogen')
    ax.set_ylabel('Nitrogen \n (mmol N)')
    ax.set_xlabel('Time (years)')
    return    
    
    
    
if __name__ == '__main__':
    import load_forcing

    days, dt, par, InitCond = load_defaults()
    forc = load_forcing.get_forcing(dt,days)
    output = run_model(days,dt,InitCond,par,forc)
#    load_forcing.plot_forcing(dt,days,forc)
    plot_model(output)
#    plot_totalNitrogen(output)
