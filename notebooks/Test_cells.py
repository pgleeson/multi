
# coding: utf-8

# In[1]:

from pyneuroml import pynml
from matplotlib import pyplot as plt


from plot_cell import plot_cell_firing

show = False
rerun = False

num_processors=12

if rerun:

    #plot_cell_firing('../cells/AllenHH/AllenHH_477127614.cell.nml', num_processors=num_processors, show_plot=show)


    #plot_cell_firing('../cells/AllenHH/AllenHH_476686112.cell.nml', num_processors=num_processors, show_plot=show)


    #plot_cell_firing('../cells/BBP/cNAC187_L23_NBC_9d37c4b1f8_0_0.cell.nml', num_processors=num_processors, show_plot=show)


    plot_cell_firing('../cells/SmithEtAl2013/L23_Retuned_477127614.cell.nml', num_processors=num_processors, show_plot=show)


refs = [('AllenHH_477127614','L23_Retuned_477127614'), ('AllenHH_476686112','cNAC187_L23_NBC_9d37c4b1f8_0_0') ]

for ref in refs:
    
    for ii in ['IF','IV']:
        point_neuron, cols = pynml.reload_standard_dat_file('%s_%s.dat'%(ref[0],ii))
        detailed_neuron, cols = pynml.reload_standard_dat_file('%s_%s.dat'%(ref[1],ii))
        print detailed_neuron

        pynml.generate_plot([point_neuron['t'], detailed_neuron['t']],
                            [point_neuron[0], detailed_neuron[0]],
                            "2 cells: %s, %s"%ref,
                            markers = ['o','o'],
                            labels = [ref[0],ref[1]],
                            show_plot_already=False)
                            
    

plt.show()

