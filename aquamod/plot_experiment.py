import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle


def plot_one_variable(output,var='Phy'):  
    fig, (ax) = plt.subplots(1,1,figsize=(13,5))
    ax.set_ylabel(var)
    ax.set_xlabel('Time (years)')
    ax.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    
    levels = sorted([float(i) for i in output.keys()])
    
    legendKeys = []  
    for floatLevel in levels:
        level = str(int(floatLevel))
        legendKeys.append(level)
        ax.plot(output[level]['time']/365,output[level][var],'-')
    
    print "Printing..."
    ax.legend(legendKeys)
    plt.show()
    return


def plot_ecosystem(output):

    levels = sorted([float(i) for i in output.keys()])

    fig, (ax, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7,1,figsize=(13,13))
    ax.set_ylabel('Phy \n (mmol N m$^{-3}$)')
    ax.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax2.set_ylabel('Zoo \n (mmol N m$^{-3}$)')
    ax2.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax3.set_ylabel('NO3 \n (mmol N m$^{-3}$)')
    ax3.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax4.set_ylabel('NH4 \n (mmol N m$^{-3}$)')
    ax4.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax5.set_ylabel('SDet \n (mmol N m$^{-3}$)') 
    ax5.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax6.set_ylabel('LDet \n (mmol N m$^{-3}$)')
    ax6.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])    
    ax7.set_ylabel('Oxy \n (mmol O2 m$^{-3}$)')
    ax7.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    ax7.set_xlabel('Time (years)')
    ax7.set_title('Planktonic Ecosystem')
    
    fig2, (ax9) = plt.subplots(1,1,figsize=(13,5))
    ax9.set_ylabel('B_conc \n (mmol N m$^{-3}$)')
    ax7.set_title('Mussel variables')
    ax9.set_color_cycle([plt.cm.jet(i) for i in np.linspace(0, 1, len(output))])
    
    legendKeys = []   
    for floatLevel in levels:
        level = str(int(floatLevel))
        legendKeys.append(level)
        ax.plot(output[level]['time']/365,output[level]['Phy'],'-')
        ax2.plot(output[level]['time']/365,output[level]['Zoo'],'-')
        ax3.plot(output[level]['time']/365,output[level]['NO3'],'-')
        ax4.plot(output[level]['time']/365,output[level]['NH4'],'-')
        ax5.plot(output[level]['time']/365,output[level]['SDet'],'-')
        ax6.plot(output[level]['time']/365,output[level]['LDet'],'-')
        ax7.plot(output[level]['time']/365,output[level]['Oxy'],'-')
        ax9.plot(output[level]['time']/365,output[level]['B_conc'],'-')
    
    ax.legend(legendKeys)
    ax9.legend(legendKeys)
    plt.show()
    return



def plot_densityVSproduction(output): 
    SeedDensity, Production = [], []
    for level in output:
        SeedDensity.append(float(level))
        Production.append((max(output[level]['B'][-int(365/0.01):])-output[level]['B'][0])*output[level]['n_muss'][0])
    
    
    df = pd.DataFrame({'SeedDensity' : SeedDensity, 'Production' : Production})
    df = df.sort_values(by='SeedDensity')
    
    fig, (ax) = plt.subplots(1,1,figsize=(9,5))
    ax.plot(df['SeedDensity'], df['Production'],'ko-')
    ax.set_xlabel('Mussel density (mussels/m$^{3}$)')
    ax.set_ylabel('Max Production (mmol N)')
    ax.set_title('Density vs Production (max of last year)')
    plt.show()
    return
    
    
def plot_densityVSoxygen(output):
    SeedDensity, MedianOxy, MinOxy, LethalOxy = [], [], [], []
    for level in output:
        SeedDensity.append(float(level))
        MedianOxy.append(np.median(output[level]['Oxy'][-int(365/0.01):]))
        MinOxy.append(min(output[level]['Oxy'][-int(365/0.01):]))
        LethalOxy = 25.
    
    df = pd.DataFrame({'SeedDensity' : SeedDensity,
                       'MedianOxy' : MedianOxy,
                       'MinOxy':MinOxy,
                       'LethalOxy':LethalOxy})
    df = df.sort_values(by='SeedDensity')
    
    fig, (ax) = plt.subplots(1,1,figsize=(9,5))
    ax.plot(df['SeedDensity'], df['MedianOxy'],'ko-')
    ax.plot(df['SeedDensity'], df['MinOxy'],'ro-')
    ax.plot(df['SeedDensity'], df['LethalOxy'],'r-',linewidth=3)
    ax.text(df['SeedDensity'].tail(1)-6, df['LethalOxy'][0]+7, 'Lethal to most fish',color='r',fontsize=14)
    ax.set_xlabel('Mussel density (mussels/m$^{3}$)')
    ax.set_ylabel('Oxygen (mmol O2 m$^{-3}$ )')
    ax.set_title('Density vs Oxygen (median and min of last year)')
    plt.show()
    return



def plot_densityVSammonia(output):
    # NH4
    SeedDensity, MedianNH4, MaxNH4, LethalNH4= [], [], [], []
    for level in output:
        SeedDensity.append(float(level))
        MedianNH4.append(np.median(output[level]['NH4'][-int(365/0.01):]))
        MaxNH4.append(max(output[level]['NH4'][-int(365/0.01):]))
        LethalNH4 = 11. # For Fry Rainbox trout
    
    df = pd.DataFrame({'SeedDensity' : SeedDensity,
                       'MedianNH4' : MedianNH4,
                       'MaxNH4':MaxNH4,
                       'LethalNH4':LethalNH4})
    df = df.sort_values(by='SeedDensity')
    
    fig, (ax) = plt.subplots(1,1,figsize=(9,5))
    ax.plot(df['SeedDensity'], df['MedianNH4'],'ko-')
    ax.plot(df['SeedDensity'], df['MaxNH4'],'ro-')
    ax.plot(df['SeedDensity'], df['LethalNH4'],'r-',linewidth=3)
    ax.text(df['SeedDensity'].tail(1)-9, df['LethalNH4'][0]+0.5, 'Lethal to Raibout trout fry',color='r',fontsize=14)
    ax.set_xlabel('Mussel density (mussels/m$^{3}$)')
    ax.set_ylabel('NH4 (mmol N m$^{-3}$ )')
    ax.set_title('Density vs Ammonia (median and max of last year)')
    plt.show()
    return



def plot_densityVSecosystem(output):
    SeedDensity, MedianPhy, MedianZoo, MedianSDet, MedianLDet,= [], [], [], [], []
    for level in output:
        SeedDensity.append(float(level))
        MedianPhy.append(np.median(output[level]['Phy'][-int(365/0.01):]))
        MedianZoo.append(np.median(output[level]['Zoo'][-int(365/0.01):]))
        MedianSDet.append(np.median(output[level]['SDet'][-int(365/0.01):]))
        MedianLDet.append(np.median(output[level]['LDet'][-int(365/0.01):]))

    df = pd.DataFrame({'SeedDensity' : SeedDensity,
                       'MedianPhy' : MedianPhy,
                       'MedianZoo' : MedianZoo,
                       'MedianSDet' : MedianSDet,
                       'MedianLDet' : MedianLDet,})
    df = df.sort_values(by='SeedDensity')
    
    fig, (ax) = plt.subplots(1,1,figsize=(9,5))
    ax.plot(df['SeedDensity'], df['MedianPhy'],'go-')
    ax.plot(df['SeedDensity'], df['MedianZoo'],'ro-')
    ax.plot(df['SeedDensity'], df['MedianSDet'],'bo-')
    ax.plot(df['SeedDensity'], df['MedianLDet'],'co-')
    ax.set_xlabel('Mussel density (mussels/m$^{3}$)')
    ax.set_ylabel('Nitrogen (mmol N m$^{-3}$ )')
    ax.set_title('Density vs Ecosystem Variables (median of last year)')
    ax.legend(['Phy','Zoo','SDet','LDet'])
    plt.show()
    return

def plot_all(output):
    plot_ecosystem(output)
    plot_densityVSproduction(output)
    plot_densityVSoxygen(output)
    plot_densityVSammonia(output)
    plot_densityVSecosystem(output)
    return

if __name__ == '__main__':
    output = pickle.load( open( 'Exp1_output.p', 'rb' ) )
    plot_all(output)