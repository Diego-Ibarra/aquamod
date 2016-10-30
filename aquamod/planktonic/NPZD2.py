
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
    par['mu0']    = 0.69  
    par['kNO3']    = 0.5    
    par['kNH4']    = 0.5  
    par['alpha']    = 0.125  
    par['gmax'] = 0.6  #Original 0.6
    par['kP']       = 0.44
    par['mP']   = 0.15    
    par['tau']   = 0.005 
    par['thetaMax']   = 0.053
    par['beta']    = 0.75 
    par['lBM']    = 0.1    
    par['lE']   = 0.1
    par['mZ']  = 0.25
    par['rSD']    = 0.3 # Original 0.03
    par['rLD'] = 0.1 # # Original 0.01
    par['nmax'] = 0.05
    par['kI'] = 0.1
    par['I0']  = 0.0095
    
    # Initial conditions
    InitCond = {}
    InitCond['Phy'] = 0.2
    InitCond['Zoo'] = 0.1
    InitCond['SDet'] = 1. 
    InitCond['LDet'] = 1.
    InitCond['NH4'] = 0.1
    InitCond['NO3'] = 7.
    InitCond['Temp'] = 6.
    #InitCond['O2'] = 0.5 
    return days, dt, par, InitCond
    



    
def run_model(days,dt,InitCond,par):
    '''
    This is your model. Do a brief description.

    '''
    # Import libraries
    import numpy as np
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    Phy = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Zoo = np.zeros((NoSTEPS,),float) # same as above
    SDet = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    LDet = np.zeros((NoSTEPS,),float) # same as above
    NH4 = np.zeros((NoSTEPS,),float) # same as above
    NO3 = np.zeros((NoSTEPS,),float) # same as above
    I  = np.zeros((NoSTEPS,),float) # same as above
    
    mu = np.zeros((NoSTEPS,),float) # same as above
    f_I = np.zeros((NoSTEPS,),float) # same as above
    L_NO3 = np.zeros((NoSTEPS,),float) # same as above
    L_NH4 = np.zeros((NoSTEPS,),float) # same as above
    TotN = np.zeros((NoSTEPS,),float) # same as above
    
    Temp = np.ones((NoSTEPS,),float) * InitCond['Temp'] #Temperature

    
    # Creating sunlight
    for i in range(len(I)):
        I[i] = 10 * np.sin((2*np.pi*time[i])/1) + \
               8 * np.sin((2*np.pi*time[i])/365)
        # We can't have negative light... so negatives are made zero 
        if I[i] < 0:
            I[i] = 0.0000001
    
    
    # Initializing with initial conditions
    Phy[0] = InitCond['Phy']
    Zoo[0] = InitCond['Zoo']
    SDet[0] = InitCond['SDet']
    LDet[0] = InitCond['LDet']
    NH4[0] = InitCond['NH4']
    NO3[0] = InitCond['NO3']
    
    
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
    # Calculate Temperature Limitation

        muMax = par['mu0'] * (1.066 ** Temp[t]) # text
        
        f_I[t] = (par['alpha']*I[t])/(np.sqrt(muMax**2+((par['alpha']**2)*(I[t]**2)))) #Eq5
        
        L_NO3[t] = max(0,((NO3[t])/(par['kNO3']*NO3[t])) * (1/(1+(NH4[t]/par['kNH4'])))) #Eq3
        
        L_NH4[t] = max(0,NH4[t]/(par['kNH4']+NH4[t])) # Eq 4
    
        mu[t] =muMax * f_I[t] * (L_NO3[t] + L_NH4[t]) # Eq2
    
        g = par['gmax'] * (Phy[t]**2/(par['kP']+(Phy[t]**2)))
        
        n = par['nmax'] * (1 - max(0,(I[t]-par['I0'])/(par['kI']+I[t]-par['I0'])))
    
        dPhydt = (mu[t] * Phy[t]) - \
                 (g  * Zoo[t]) - \
                 (par['mP'] * Phy[t]) - \
                 (par['tau']*(SDet[t]+Phy[t])*Phy[t]) # Eq1
                 
        dZoodt = (g * par['beta'] * Zoo[t]) - \
                 (par['lBM']*Zoo[t]) - \
                 (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                 (par['mZ']*(Zoo[t]**2))#Eq10
                 
        dSDetdt = (g * (1-par['beta']) * Zoo[t]) + \
                  (par['mZ']*(Zoo[t]**2)) + \
                  (par['mP'] * Phy[t]) - \
                  (par['tau']*(SDet[t]+Phy[t])*SDet[t]) - \
                  (par['rSD']*SDet[t])
                  
        dLDetdt = (par['tau']*((SDet[t]+Phy[t])**2)) - \
                  (par['rLD']*LDet[t])
                  
        dNO3dt = -(muMax * f_I[t] * L_NO3[t] * Phy[t]) + \
                  (n * NH4[t])
                 
        dNH4dt = -(muMax * f_I[t] * L_NH4[t] * Phy[t]) - \
                  (n * NH4[t]) + \
                  (par['lBM'] * Zoo[t]) + \
                  (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) + \
                  (par['rSD']*SDet[t]) + \
                  (par['rLD']*LDet[t])
                  
    
        # Update and step ----------------------------------------------------
        Phy[t+1]  = Phy[t]  + (dPhydt * dt)
        Zoo[t+1]  = Zoo[t]  + (dZoodt * dt)
        SDet[t+1] = SDet[t] + (dSDetdt * dt)
        LDet[t+1] = LDet[t] + (dLDetdt * dt)
        NH4[t+1]  = NH4[t]  + (dNH4dt * dt)
        NO3[t+1]  = NO3[t] +  (dNO3dt * dt)
        if NO3[t+1] <= 0.0001:
            offset = NO3[t+1]
            NH4[t+1] = NH4[t+1] + offset
            NO3[t+1] = NO3[t+1] - offset
            
        # Estimate Total Nitrogen
        TotN[t+1] = Phy[t+1] + Zoo[t+1] + SDet[t+1] + LDet[t+1] + NH4[t+1] + NO3[t+1]
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['Phy'] = Phy
    output['Zoo'] = Zoo
    output['SDet'] = SDet
    output['LDet'] = LDet
    output['NH4'] = NH4
    output['NO3'] = NO3
    output['mu'] = mu
    output['f_I'] = f_I
    output['L_NO3'] = L_NO3
    output['L_NH4'] = L_NH4
    output['TotN'] = TotN

    print "Model run: DONE!!!"
    return output


    
def plot_model(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax, ax2) = plt.subplots(2,1)

    ax.plot(output['time']/365,output['mu'],'g-')
    ax.plot(output['time']/365,output['f_I'],'r-')
    ax.plot(output['time']/365,output['L_NO3'],'b-')
    ax.plot(output['time']/365,output['L_NH4'],'k-')
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax.set_title('Fennel et al 2006 Model')
    ax.legend(['Mu','f_I','L_NO3','L_NH4'])
    
    ax2.plot(output['time']/365,output['Phy'],'g-')
    ax2.plot(output['time']/365,output['Zoo'],'r-')
    ax2.plot(output['time']/365,output['SDet'],'k-')
    ax2.plot(output['time']/365,output['LDet'],'k-.')
    ax2.plot(output['time']/365,output['NH4'],'m-')
    ax2.plot(output['time']/365,output['NO3'],'c-')
    ax2.plot(output['time']/365,output['TotN'],'y-')
    ax2.set_xlabel('Time (years)')
    ax2.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax2.set_title('Fennel et al 2006 Model')
    plt.legend(['Phy','Zoo','SDet','LDet','NH4','NO3','TotN'])
    plt.show()
    return
    
if __name__ == '__main__':
    days, dt, par, InitCond = load_defaults()
    output = run_model(days,dt,InitCond,par)
    plot_model(output)