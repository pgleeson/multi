
from pyneuroml.analysis import generate_current_vs_frequency_curve

def plot_cell_firing(cell_file, num_processors=1,show_plot=True):
    
    cell_id = cell_file.split('/')[-1].split('.')[0]
    
    simulator = 'jNeuroML_NEURON' if num_processors==1 else 'jNeuroML_NetPyNE'

    curve = generate_current_vs_frequency_curve(cell_file, 
                                    cell_id, 
                                    custom_amps_nA =       [-0.1,-0.08,-0.06,-0.04,-0.02,0.0,0.02,0.04,0.06,0.08,0.1,0.12], 
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
                                    show_plot_already=     show_plot)
                                    
    # Longer duration, more points 
    
    curve = generate_current_vs_frequency_curve(cell_file, 
                                    cell_id, 
                                    start_amp_nA =         -0.1, 
                                    end_amp_nA =           0.4, 
                                    step_nA =              0.025, 
                                    analysis_duration =    2000, 
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
                                    show_plot_already=     show_plot)