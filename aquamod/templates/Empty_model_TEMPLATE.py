


def load_defaults():
    '''
    This function creates a dictionaries called "par" and "InitCond"
    and pre-loads them with all the default 
    parameters and initial conditions, respectively.
    Also outputs days and dt
    '''
    # Framework
    days = 365 # One year
    dt   = 0.01 # units: days  

    
    # Parameters (defaults)
    par = {}
    par['Ae'] = 0.8
    par['F'] = 1.5
    par['k'] = 0.013
    par['T']    = 26.  # Temperature (units: C)
    par['R32']   = 0.5 #This is OUR assumption... not clear in paper
    par['psi']    = 5.  # Dry weight/fresh weight (units: none)
    
    # Initial conditions
    InitCond = {}
    InitCond['G']  = 0.1 # Variable name (units: d^-1)
    InitCond['Bt']  = 1.
    
    return days, dt, par, InitCond
    



    
def run_model(days,dt,InitCond,par):
    '''
    This is your model. Do a brief description.
    
    INPUTS:
        days: number of days of simulation
        dt: time steps (units: days)
        InitCond: Dictionary with all initial conditions
        par: Dictionary with all model parameters
        
    OUTPUTS:
        var1: name (units)
        var2: name (units)
        var3: name (units)
    
    Don't forget to reference the paper where you got it    
    '''
    # Import libraries
    import numpy as np
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    #Make empty arrays
    G = np.zeros((NoSTEPS,),float) # same as above
    Bt = np.zeros((NoSTEPS,),float) 
#    
    # Initializing with initial conditions
    G[0] = InitCond['G']
    Bt[0] = InitCond['Bt']
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        # just as an example
    
        L_T = ((-0.02*(par['T']**2))+(1.44*par['T'])-17.41)/(par['R32']) # Eq. 7
         
        L_Bt = 0.09*(par['psi']*(Bt[t])**0.62) # Eq. 6

        Rmax = L_Bt * L_T # Eq. 5
        
        dGdt = (Rmax * (1-np.exp(-par['k']*par['F']))) - (G[t] * par['Ae']) # Eq. 12

    # END of MAIN MODEL LOOP ******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['var1'] = var1
    output['var2'] = var2
    output['var3'] = var3

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
    ax.plot(output['time'],output['var1'],'b-')
    ax.plot(output['time'],output['var2'],'g-')
    ax.plot(output['time'],output['var3'],'r-')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('My Units (units here)')
    ax.set_title('TEMPLATE example Model Simulation')
    plt.legend(['var1','var2','var3'])
    plt.show()
    return