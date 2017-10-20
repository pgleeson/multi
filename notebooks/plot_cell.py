import sys 
from pyneuroml import pynml

from pyneuroml.analysis import generate_current_vs_frequency_curve

def plot_cell_firing(cell_file, num_processors=1,quick=False,show_plot=True, plot2_duration=8000):
    
    cell_id = cell_file.split('/')[-1].split('.')[0]
    
    simulator = 'jNeuroML_NEURON' if num_processors==1 else 'jNeuroML_NetPyNE'
    
    custom_amps_nA =       [-0.1,-0.08,-0.06,-0.04,-0.02,0.0,0.02,0.04,0.06,0.08,0.1,0.12]
    if quick:
        custom_amps_nA =       [-0.08,-0.04,0.0,0.04,0.08,0.12]
        
    
    curve = generate_current_vs_frequency_curve(cell_file, 
                                    cell_id, 
                                    custom_amps_nA =       custom_amps_nA, 
                                    analysis_duration =    1000, 
                                    pre_zero_pulse =       100,
                                    post_zero_pulse =      100,
                                    analysis_delay =       0,
                                    dt =                   0.025,
                                    simulator =            simulator,
                                    num_processors =       num_processors,
                                    plot_voltage_traces =  True,
                                    plot_if =              False,
                                    plot_iv =              False,
                                    temperature =          '34degC',
                                    title_above_plot =      True,
                                    save_voltage_traces_to = "%s_traces.png"%cell_id,
                                    show_plot_already=     show_plot,
                                    verbose = True)
                                    
    # Longer duration, more points 
    
    if not quick:
        curve = generate_current_vs_frequency_curve(cell_file, 
                                    cell_id, 
                                    start_amp_nA =         -0.1, 
                                    end_amp_nA =           0.4, 
                                    step_nA =              0.025, 
                                    analysis_duration =    plot2_duration, 
                                    dt =                   0.025,
                                    pre_zero_pulse =       0,
                                    post_zero_pulse =      0,
                                    analysis_delay =       100,
                                    simulator =            simulator,
                                    num_processors =       num_processors,
                                    plot_voltage_traces =  False,
                                    plot_if =              True,
                                    plot_iv =              True,
                                    temperature =          '34degC',
                                    save_if_figure_to =    "%s_IF.png"%cell_id,
                                    save_iv_figure_to =    "%s_IV.png"%cell_id,
                                    title_above_plot =      True,
                                    save_if_data_to =      "%s_IF.dat"%cell_id,
                                    save_iv_data_to =      "%s_IV.dat"%cell_id,
                                    show_plot_already=     show_plot,
                                    verbose = True)
                                    


if __name__ == '__main__':
    
    if '-summary' in sys.argv:
        
        erefs = ['AllenHH_477127614','L23_Retuned_477127614']#,'L23_Retuned_477127614']
        irefs = ['AllenHH_476686112']
        
        iv_xs = []
        iv_ys = []
        iv_labels = []
        if_xs = []
        if_ys = []
        if_labels = []
        for ref in erefs+irefs:
            
            ivf = '%s_IV.dat'%ref
            data, indices = pynml.reload_standard_dat_file(ivf)
            print("Loaded %s from %s"%(data.keys(),ivf))
            iv_xs.append(data['t'])
            iv_ys.append(data[0])
            iv_labels.append(ref)
            
            iff = '%s_IF.dat'%ref
            data, indices = pynml.reload_standard_dat_file(iff)
            print("Loaded %s from %s"%(data.keys(),iff))
            if_xs.append(data['t'])
            if_ys.append(data[0])
            if_labels.append(ref)
            
            
        pynml.generate_plot(iv_xs,iv_ys, "IV curves", labels=iv_labels, markers=['o']*len(iv_labels), show_plot_already=False)
        pynml.generate_plot(if_xs,if_ys, "IF curves", labels=if_labels, markers=['o']*len(if_labels))
    else:
    
        #plot_cell_firing('../cells/AllenHH/AllenHH_477127614.cell.nml')
        #plot_cell_firing('../cells/AllenHH/AllenHH_476686112.cell.nml')
        #plot_cell_firing('../cells/SmithEtAl2013/L23_Retuned_477127614.cell.nml')
        plot_cell_firing('../cells/Thalamocortical/L23PyrRS.cell.nml',quick=False,num_processors=18,plot2_duration=3000)