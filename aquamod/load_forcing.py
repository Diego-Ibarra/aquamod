import pandas as pd
import scipy.interpolate as intrp
import numpy as np

def load():
    File1 = 'TempSaltOxy.csv'
    
    File2 = 'Nitrate.csv'
    
    File3 = 'PhyZooSDetLDet.csv'
    
    data = pd.read_csv(File1)
    
    forc = {}
    forc['Temp'] = data['temperature_an'][1:].values.astype(float)
    forc['Salt'] = data['salinity_an'][1:].values.astype(float)
    forc['Oxy'] = data['disOxygen_an'][1:].values.astype(float)
    
    data = pd.read_csv(File2)
    
    forc['NO3'] = data['nitrate_an'][1:].values.astype(float)
    
    data = pd.read_csv(File3)
        
    forc['Phy'] = data['phytoplankton'][1:].values.astype(float)
    forc['Zoo'] = data['zooplankton'][1:].values.astype(float)
    forc['SDet'] = data['SDet'][1:].values.astype(float)
    forc['LDet'] = data['LDet'][1:].values.astype(float)
    forc['NH4'] = data['NH4'][1:].values.astype(float)
        
    return forc

def interp(dt,days,forc):
    # Resize array and Eliminate Gaps by Interpolating in between
 
    forc = load()
    
    interp_forc = {}
    
    for key in forc:
        data = forc[key]
        
        # Copy-paste number of years
        for i in range(0,int(days/365)-1):
            data = np.append(data,data)
            
        months = len(data)
        x = np.arange(0, months)
        f = intrp.interp1d(x, data, kind='linear', fill_value='extrapolate' )

        NoSTEPS = int(days / dt)
        newx = np.linspace(0,days/30,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
        new_data = f(newx)
        
        # Get rid of nans
        nans, x= np.isnan(new_data), lambda z: z.nonzero()[0]
        new_data[nans]= np.interp(x(nans), x(~nans), new_data[~nans])
        interp_forc[key] = new_data
        
    return interp_forc
    
def get_forcing(dt,days):
    forc = load()
    new_forc = interp(dt,days,forc)
    new_forc = add_analytical_light(dt,days,new_forc)
    return new_forc

def add_analytical_light(dt,days,forc):
    
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    I  = np.zeros((NoSTEPS,),float) # same as above
    
    # Creating sunlight
    for i in range(len(I)):
        I[i] = 10 * np.sin((2*np.pi*time[i])/1) + \
               8 * np.sin(((2*np.pi*time[i])/365)-2100)
        # We can't have negative light... so negatives are made zero 
        if I[i] < 0:
            I[i] = 0.0000001
            
    forc['I'] = I
    return forc
    
def plot_forcing(dt,days,forc):
    import matplotlib.pyplot as plt
    
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    fig, (ax, ax2, ax3, ax4, ax5) = plt.subplots(5,1,figsize=(13,13))
    ax.set_title('FORCING - Properties of "outer" box')
    ax.plot(time/365,forc['I'],'r-')
    ax.set_ylabel('Sunlight \n (W m$^{-2}$)')
    ax2.plot(time/365,forc['Oxy'],'b-')
    ax2.set_ylabel('Oxygen \n (mmol O2 m$^{-3}$)')
    ax3.plot(time/365,forc['Phy'],'g-')
    ax3.plot(time/365,forc['Zoo'],'r-')
    ax3.plot(time/365,forc['SDet'],'k-')
    ax3.plot(time/365,forc['LDet'],'k-.')
    ax3.plot(time/365,forc['NH4'],'m-')
    ax3.plot(time/365,forc['NO3'],'b-')
    ax3.set_ylabel('Nitrogen \n (mmol N m$^{-3}$)')
    ax3.legend(['Phy','Zoo','SDet','LDet','NH4','NO3'])
    ax4.plot(time/365,forc['Temp'],'r-')
    ax4.set_ylabel('Temperature \n (oC)')
    ax5.plot(time/365,forc['Salt'],'c-')
    ax5.set_ylabel('Salinity \n (ppt)')
    ax5.set_xlabel('Time (years)')
    return

if __name__ == '__main__':
    dt = 0.01
    days = 365*3
    forc = get_forcing(dt,days)
    plot_forcing(dt,days,forc)
    