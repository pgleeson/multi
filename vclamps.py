from  pyneuroml import pynml
import os
import matplotlib.pyplot as plt
import numpy as np

# v_clamps_i_seg0_EPSC.12345.dat  v_clamps_i_seg0_IPSC.12345.dat  v_clamps_i_seg1406_EPSC.12345.dat  v_clamps_i_seg1406_IPSC.12345.dat  v_clamps_i_seg2953_EPSC.12345.dat  v_clamps_i_seg2953_IPSC.12345.dat

lims = {'EPSC':(-0.16e-8,0.15e-9),'IPSC':(0.9e-8,1.2e-8)}
lims2 = {'EPSC':(-7e-11,0e-10),'IPSC':(0e-8,1.7e-10)}


for x in lims.keys():
    
    print("=====\nPlot for %s"%x)
    xs = []
    ys = []
    labels = []
    all_seg_0 = []
    all_seg_1406 = []
    for f in os.listdir('./temp'):
        if f.startswith('v_clamps') and x in f and f.endswith('dat'):
            print('Adding %s'%f)
            seed = f.split('.')[1]
            seg = f.split('.')[0].split('_')[3]
            res, i = pynml.reload_standard_dat_file('temp/%s'%f)
            print res.keys()
            time_points = res['t']
            xs.append(res['t'])
            ys.append(res[0])
            labels.append('%s %s'%(seg,seed))
            if seg == 'seg0':
                all_seg_0.append(res[0])
            if seg == 'seg1406':
                all_seg_1406.append(res[0])
                
    
    avg_seg0 = np.zeros(len(all_seg_0[0]))
    for l in all_seg_0:
        for i in range(len(l)):
            avg_seg0[i]+=l[i]
            
            
    avg_seg1406 = np.zeros(len(all_seg_1406[0]))
    for l in all_seg_1406:
        for i in range(len(l)):
            avg_seg1406[i]+=l[i]
            
 
    pynml.generate_plot(xs,
                        ys,
                        'V clamp currents %s'%x,
                        ylim=lims[x],
                        labels=labels,
                        show_plot_already=False)
    
    xs = []
    ys = []
    labels = []
    
    xs.append(time_points)
    ys.append(avg_seg0/len(all_seg_0))
    labels.append('%s average'%(seg))
    
    pynml.generate_plot(xs,
                        ys,
                        'V clamp current averages at soma %s'%x,
                        ylim=lims[x],
                        labels=labels,
                        show_plot_already=False)
    
    xs = []
    ys = []
    labels = []
    
    xs.append(time_points)
    ys.append(avg_seg1406/len(all_seg_1406))
    seg = 'seg1406'
    labels.append('%s average'%(seg))
    
    pynml.generate_plot(xs,
                        ys,
                        'V clamp current averages at seg1406 (dendrite) %s'%x,
                        ylim=lims2[x],
                        labels=labels,
                        show_plot_already=False)
    
plt.show()
