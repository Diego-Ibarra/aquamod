


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

    par['ephy']=2.86
    par['edet']=1.68
    par['eO2']=20
    par['aet']=0.033
    par['aeO']=0.033
    par['Fmax']=48
    par['pj']=0.05
    par['bf']=0.4
    par['br']=0.7
    par['bp']=1.28
    par['Tses']=200.
    par['kf']=0.07
    par['PN']=0.4
    par['c1']=120
    par['c2']=2400
    par['kp1']=0.15
    par['kp2']=0.01
    par['pf0']=0.4
    par['art']=0.768
    par['ar0']=-0.528
    par['ap'] = 0.0057
    par['dp'] = 58.
    
    
    #---------
    par['NPHY']=1.
    par['NDET']=4.
    par['Temp']=1.
    par['SES']=50
    
    # Initial conditions
    InitCond = {}
    InitCond['Oys'] = 0.2

    
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
    Oys = np.zeros((NoSTEPS,),float) # same as above
    sfg= np.zeros((NoSTEPS,),float) # same as above
    assoys= np.zeros((NoSTEPS,),float)
    Ingphy= np.zeros((NoSTEPS,),float)
    Ingdet= np.zeros((NoSTEPS,),float)
    F= np.zeros((NoSTEPS,),float)
    R= np.zeros((NoSTEPS,),float)
#    
    # Initializing with initial conditions
    Oys[0] = InitCond['Oys']

    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        # just as an example
    
        Wbf = (Oys[t]*par['pj'])**par['bf']
    
        F[t]= par['Fmax'] * np.exp(par['kf']*min(0,par['Tses']-par['SES'])) * Wbf

        ct = (F[t]*(par['SES']+((par['NPHY']+par['NDET'])*par['PN']))) / Wbf        
        
        PF = (par['pf0']*(1-np.exp(par['kp1']*min(0,par['c1']-ct)))) + ((1 - par['pf0'])*(1-np.exp(par['kp2']*(min(0,par['c2']-ct)))))
        
        Ingphy[t]= F[t] * par['NPHY'] * (1-PF)

        Ingdet[t]= F[t] * par['NDET'] * (1-PF)

        assoys[t]= (par['aet']*par['Temp']) + par['aeO']
        
        R[t]=((par['art']* par['Temp']) + par['ar0']) * ((Oys[t]*par['pj'])**par['br'])
        
        sfg[t]=assoys[t]*(Ingphy[t]*par['ephy'] + Ingdet[t]*par['edet']) - R[t]*par['eO2']
        
        spawn = (par['ap'] * (Oys[t]**par['bp']))**(1/par['dp'])
        
        dOysdt = sfg[t] - spawn
        
        # Time stepping
        Oys[t+1] = Oys[t] + (dOysdt * dt)

    # END of MAIN MODEL LOOP ******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['Oys'] = Oys
    output['W'] = Oys * par['pj']
    output['sfg'] = sfg
    output['assoys'] = assoys
    

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
    ax.plot(output['time'],output['W'],'b-')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Oyster Dry Weight (g) (units here)')
    ax.set_title('Oyster Model')
    plt.legend(['Oys'])
    plt.show()
    return
    
if __name__ == "__main__":
    days, dt, par, InitCond = load_defaults()
    output = run_model(days,dt,InitCond,par)
    plot_model(output)