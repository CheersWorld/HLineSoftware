import main
import globalVarsManager as gvars
import configurationManager as cm

cm.setParameters()
runs = 0
gvars.obs['duration'] = 20
gvars.obs['t_sample'] = .04
for i in range(0, 3):
    gvars.obs['rf_gain'] = i * 10
    for j in range(0, 3):
        gvars.obs['if_gain'] = j * 12
        for k in range(0, 3):
            gvars.obs['bb_gain'] = k * 6
            runs += 1
            print("Run " + str(runs) + " of 27")
            main.observe(gvars.obs, gvars.args.noPlot)            
print(runs)