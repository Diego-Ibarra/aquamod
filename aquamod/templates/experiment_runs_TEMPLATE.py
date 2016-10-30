import Empty_model_TEMPLATE as mymodel


# load default parameters
days, dt, par, InitCond = mymodel.load_defaults()

# Change a few things
days = 100
InitCond['dummy_InitialCondition_1']  = 7.5

 # run the model
output = mymodel.run_model(days,dt,InitCond,par)

# plot model
mymodel.plot_model(output)