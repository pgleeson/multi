from  pyneuroml import pynml
import os
import matplotlib.pyplot as plt

# v_clamps_i_seg0_EPSC.12345.dat  v_clamps_i_seg0_IPSC.12345.dat  v_clamps_i_seg1406_EPSC.12345.dat  v_clamps_i_seg1406_IPSC.12345.dat  v_clamps_i_seg2953_EPSC.12345.dat  v_clamps_i_seg2953_IPSC.12345.dat

lims = {'EPSC':(-1e-8,1e-9),'IPSC':(-1e-9,1e-8)}

for x in lims.keys():
    
    print("=====\nPlot for %s"%x)
    xs = []
    ys = []
    labels = []
    for f in os.listdir('./temp'):
        if f.startswith('v_clamps') and x in f and f.endswith('dat'):
            print('Adding %s'%f)
            seed = f.split('.')[1]
            seg = f.split('.')[0].split('_')[3]
            res, i = pynml.reload_standard_dat_file('temp/%s'%f)
            print res.keys()
            xs.append(res['t'])
            ys.append(res[0])
            labels.append('%s %s'%(seg,seed))
            
    pynml.generate_plot(xs,
                        ys,
                        'V clamp currents %s'%x,
                        ylim=lims[x],
                        labels=labels,
                        show_plot_already=False)
    
plt.show()
