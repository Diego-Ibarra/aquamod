'''A box model is used to decribe the growth of mussels, mainly Mytilus edulis, in a small aquaculture site at Upper South Cove
near Lunenburg Nova Scotia. The ecological interactions in the model include 2 competing herbivores, mussels and zooplankton, 
and 2 food sources, phytoplankton and non-plankton seston.

Dowd (1997) On predicting the growth o cultured bivalves. Ecological modelling 104 (1997) 113-131
'''

def load_defaults():
    '''
    This function creates a dictionaries called "par" and "InitCond"
    and pre-loads them with all the default 
    parameters and initial conditions, respectively.
    Also outputs days and dt
    '''
    # Framework
    days = 365 * 5 # One year
    dt   = 0.01 # units: days    
    
    # Parameters (defaults)
    par = {}
    par['I_M']        = 0.1  # Ingestion rate for mussel (units: d^-1) 
    par['R_M']        = 0.01 # Respiration rate for mussel
    par['epsilon_MP'] = 0.9 # Assimilation efficiency for mussels on phytoplankton
    par['epsilon_MS'] = 0.2 # Assimilation efficiency for mussels on seston
    par['mu_M']       = 0.8 # Selection factor for mussels
    par['lambda_M']   = -0.002 # Mortality rate for mussels  (units: d^-1)
    par['kappa_M']    = 1000 # half-saturation constant for mussel ingestion (units: gC m^-3)
    par['Q_MI']       = 0.07 #  Temperature rate constant (units: degree C^-1)
    par['Q_MR']       = 0.07 #  (units:degree C^-1)
    par['b']          = -2 # Allometric exponent
    par['T']          = 10 # temperature (units: degree C)
    par['D_M']        = 0 # Spawning parameter, set to zero for simplicity
    
    # Diego's: I added 3 new parameters =====================================================================================
#    par['C_MI'] = 70
#    par['C_MR1'] = .01
#    par['C_MR2'] = .01
    par['C_MI'] = 40.
    par['C_MR1'] = .01
    par['C_MR2'] = .01
 
   # Initial conditions
    InitCond = {}
    InitCond['M']    = 0.015 # mussel weight
    InitCond['P']    = 0.1 # phytoplanton (units: gC m-3) 
    InitCond['S']    = 1.0 # Seston (units: gC m-3)
   
    
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
    import math
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    #Make empty arrays
    M = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
#    
    # Initializing with initial conditions
    M[0] = InitCond['M']
    P = InitCond['P']
    S = InitCond['S']
 
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):

        # Diego's: I replaced "1s" for parameters par['C_MI'], par['C_MR1'], par['C_MR2'] ===========================================
        f_MI = par['C_MI']*((P+S)/(par['kappa_M']+P+S))*(np.exp(par['Q_MI']*par['T']))*(M[t]**par['b'])
        f_MR = (par['C_MR1']*(M[t]**par['b']))+(par['C_MR2']*(P+S))*(np.exp(par['Q_MR']*par['T']))*(M[t]**par['b']) 
        # Diego's: I was using these prints to see the values of f_MI and f_MR  during the run  
        #print f_MI
        #print f_MR        
        
        # The growth rate of an individual mussel (Eq. 1)
        dMdt = ((par['epsilon_MP']*P)/(P+par['mu_M']*S)+(par['epsilon_MS']*par['mu_M']*S)/(P+par['mu_M']*S))* \
                (f_MI*par['I_M']*M[t])-(f_MR*par['R_M']*M[t])-(par['D_M'])
        
        
        # time stepping
        M[t+1] = M[t] + (dMdt * dt) 
        
    # END of MAIN MODEL LOOP ******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['M'] = M

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
    ax.plot(output['time']/365,output['M'],'b-')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Mussel Dry Weight (gC)')
    ax.set_title('Model predicted growth trajectories for mussels in a coastal inlet near Lunenburg Nova Scotia')
    plt.show()
    return
    
    
# Diego's I added this so that you can run the model by running THIS script, rather than by using the "experiment_run.py" script
if __name__ == '__main__':
    # load default parameters
    days, dt, par, InitCond = load_defaults()
    
     # run the model
    output = run_model(days,dt,InitCond,par)
    
    # plot model
    plot_model(output)