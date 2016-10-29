import NPZD2_Mussels_Diego_Chi_n_Oxy_Forcing_Sinking as model
import load_forcing
import pickle

days, dt, par, InitCond = model.load_defaults()
forc = load_forcing.get_forcing(dt,days)

MussConc_levels = [0,1,2,5,10,15,20]

output = {}
for level in MussConc_levels:
    InitCond['conc_muss'] = float(level)
    output[str(level)] = model.run_model(days,dt,InitCond,par,forc)
    
pickle.dump( output, open( 'Exp1_output.p', 'wb' ) )

print 'Experiment is DONE! and saved in: model_output.p'