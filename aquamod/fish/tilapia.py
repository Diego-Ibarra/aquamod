def load_defaults():
   
    # Framework
    days = 365 * 1 # One year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    par['b'] = 0.62 # efficiency of food assimilation(dimensionless)
    par['a'] = 0.53 # fraction of food assimilated used for feeding catabolism(dimensionless)
    par['h'] = 0.8 #(dimensionless)
    par['s'] = 17.31 #(dimensionless)
    par['m'] = 0.67 #(dimensionless)
    par['n'] = 0.81 #(dimensionless)
    par['PNPP/B'] = 1 # Potential net primary production (PNPP) to standing crop (B) of Tilapia (dimensionless)
    par['kmin'] = 0.00133
    par['j'] = 0.0132
    par['T'] = 35 # Temperature (Degrees C)
    par['Tmin'] = 15 # Minimum Temperature for survival (Degrees C)
    par['DO'] = 0.9 # Dissolved Oxyen mg/l
    par['DOmin'] = 0.3 #Minimum Dissolved Oxygen mg/l
    par['DOcrit'] = 1 # Critical Dissolved Oxygen mg/l
    par['UIA'] = 0.07 # Unionized Ammonia Concentration mg/l
    par['UIAmax'] = 1.4  # maximum Unionized Ammonia Concentrationg/l
    par['UIAcrit'] = 0.06  # critical Unionized Ammonia Concentration mg/l
    

    # Initial conditions
    InitCond = {}
    InitCond['W'] = 13 # (g)
    return days, dt, par, InitCond
    

    
def run_model(days,dt,InitCond,par):
    # Import libraries
    import numpy as np
    
    # Setup the framework 
    NoSTEPS = int(days / dt) # Calculates the number of steps 
    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers 
    
    # Create arrays of zeros
    W = np.zeros((NoSTEPS,),float)
    # Initializing with initial conditions
    W[0] = InitCond['W']
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        
        if par['DO'] > par['DOcrit']:
            delta = 1
        elif par['DOmin'] <= par['DO'] and par['DO'] <= par['DOcrit']:
            delta = (par['DO']-par['DOmin'])/(par['DOcrit']-par['DOmin'])
        elif par['DO'] < par['DOmin']:
            delta = 0
        # Equation 3  (a,b,c) respectively
        if par['UIA'] < par['UIAcrit']:      
            upsilon = 1 # !!!TO DO
        elif par['UIAcrit'] <= par['UIA'] and par['UIA'] <= par['UIAmax']:
            upsilon = (par['UIAmax']-par['UIA'])/(par['UIAmax'] - par['UIAcrit'])
        elif par['UIA'] > par['UIAmax']:
            upsilon = 0
        # Equation 4 (a,b,c) respectively
        k = par['kmin'] * np.exp(par['j']*(par['T'] - par['Tmin']))  # Equation 12
        
        dWdt = ((par['b']*(1 - par['a'])*delta*upsilon*par['h']*(1 - np.exp(-par['s']*par['PNPP/B']))) * (W[t]**par['m'])) - (k*(W[t]**par['n'])) #Equation 11 - Final Growth Rate    
        
        W[t+1] = W[t] + (dWdt*dt)
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['W'] = W

    print "Model run: DONE!!!"
    return output


    
def plot_model(output):
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax) = plt.subplots(1,1)
    ax.plot(output['time'],output['W'],'g-')
    plt.ylabel('Fish Weight (g)')
    plt.xlabel('Time (days)')
    plt.show()
    return
    

if __name__ == "__main__":
    days, dt, par, InitCond = load_defaults()
    output = run_model(days,dt,InitCond,par)
    plot_model(output)