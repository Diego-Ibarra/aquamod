


def load_defaults():
    '''
    This function creates a dictionaries called "par" and "InitCond"
    and pre-loads them with all the default 
    parameters and initial conditions, respectively.
    Also outputs days and dt
    '''
    # Framework
    days = 155 # One year
    dt   = 0.01 # units: days    
    
    # Parameters (defaults)
    par = {}
    par['psi']    = 5.  # Dry weight/fresh weight (units: none)
    par['k']    = 0.12  # Consumption efficiency (units: mgDW^-1 m^2)
    par['Linf']    = 29.9  # Asymptotic length (units: mm)
    par['K']    = 0.013  # Growth coefficient (units: day^-1)
    par['t0']    = 0.06  # Hypothetical age at size 0 (units: day)
    par['Ae']    = 0.8  # Assimilation efficiency (units: day^-1)
    par['T']    = 40.  # Temperature (units: C)
    par['GSI']    = 0.05  # Gonadosomatic index (units: day^-1)
    par['Ew']    = 2*(10**-6)  # Egg weight (units: gDW)
    par['R32']   = 0.2 # This is OUR assumption... not clear in paper
    par['F']  = 2.0 # This is OUR assumption... not clear in paper
    
    # Initial conditions
    InitCond = {}
    InitCond['Bt']  = 1. # Total biomass (units: gDW) 
    InitCond['G']  = 0.1 # Gut content (units: gDW)
    InitCond['Bg'] = 0.03 # Reproductive (non-somatic) biomass (units: gDW)
 
    
    return days, dt, par, InitCond
    



    
def run_model(days,dt,InitCond,par):
    '''
    This is a growth physiology model of a penaeid shrimp based on a publication:
    A.R. Franco, J.G. Ferreira, & A.M. Nobre. (2006). Development of a growth model for penaeid shrimp. Aquaculture, 259, pp. 268–277.
    
    INPUTS:
        days: 155
        dt: time steps (units: days)
        InitCond: Dictionary with all initial conditions
        par: Dictionary with all model parameters
        
    OUTPUTS:
        var1: Bt (total biomass) (gDW)
        var2: G (gonad tissue biomass) (gDW)
        var3: Bs (somatic tissue biomass) (gDW)
    
    
    
    '''
    # Import libraries
    import numpy as np
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    #Make empty arrays
    Bt = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    G = np.zeros((NoSTEPS,),float) # same as above
    M = np.zeros((NoSTEPS,),float) # same as above
    Bg = np.zeros((NoSTEPS,),float) # same as above
    Bs = np.zeros((NoSTEPS,),float) # same as above
#    
    # Initializing with initial conditions
    Bt[0] = InitCond['Bt']
    G[0] = InitCond['G']
    Bg[0] = InitCond['Bg']
    Bs[0] = InitCond['Bt'] - InitCond['Bg']
#    var2[0] = InitCond['dummy_InitialCondition_2']
#    var3[0] = InitCond['dummy_InitialCondition_3']
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        # just as an example
        

        L_T = ((-0.02*(par['T']**2))+(1.44*par['T'])-17.41)/(par['R32']) # Eq. 7
       
        L_Bt = 0.09*(par['psi']*(Bt[t])**0.62) # Eq. 6
        
        Rmax = L_Bt * L_T # Eq. 5
        
        Mw = (0.05 * np.exp(0.07*par['T'])) * ((par['psi']*Bt[t])**-0.185)    
        
        M[t] = Mw * Bt[t]
        
        N =  9.3 * (Bt[t]**2.8)   # 15    
        
        S = N * par['Ew'] # Eq. 16
        
        
        dBgdt = Bt[t] * par['GSI'] #Eq. 14
        
        dGdt = (Rmax * (1-np.exp(-par['k']*par['F']))) - (G[t] * par['Ae']) # Eq. 12
        
        dBtdt = (G[t] * par['Ae']) - M[t]  - S # Eq. 13
        
        # Calculating future time step
        Bg[t+1] = Bg[t] + (dBgdt * dt)
        G[t+1] = G[t] + (dGdt * dt)
        Bt[t+1] = Bt[t] + (dBtdt * dt)
        Bs[t+1] = Bt[t+1] - G[t+1] # somatic tissue biomass (gDW)

    # END of MAIN MODEL LOOP ******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['Bt'] = Bt
    output['G'] = G
    output['M'] = M
    output['Bg'] = Bg
    output['Bs'] = Bs

    print "Model run: DONE!!!"
    
    return output



    
def plot_model(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    # Plotting
    fig, (ax) = plt.subplots(1,1)
    ax.plot(output['time'],output['Bt'],'b-')
    ax.plot(output['time'],output['G'],'g-')
#    ax.plot(output['time'],output['Bg'],'r-')
    ax.plot(output['time'],output['Bs'],'y-')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Biomass (gDW)')
    ax.set_title('Penaeid Shrimp Physiology Model')
    plt.legend(['Bt','G','Bs'])
    plt.show()
    return
    
    
if __name__ == "__main__":
    # load default parameters
    days, dt, par, InitCond = load_defaults()
     # run the model
    output = run_model(days,dt,InitCond,par)
    # plot model
    plot_model(output)
    
    

    