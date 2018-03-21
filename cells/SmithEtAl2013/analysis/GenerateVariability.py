

import opencortex.core as oc
from pyelectro import analysis
from pyneuroml import pynml


def generate(cell_id, duration, reference, 
             format='hdf5',
             num_cells = 10,
             simulator=None):

    #Insyn = int(Ensyn * 0.2)
    #bInsyn = int(bEnsyn * 0.2)
    
    cell_file = '../%s.cell.nml'%cell_id
    if '/' in cell_id:
        cell_id=cell_id.split('/')[-1]

    nml_doc, network = oc.generate_network(reference, temperature='35degC')

    oc.include_neuroml2_cell_and_channels(nml_doc,cell_file,cell_id)


    pop = oc.add_population_in_rectangular_region(network,
                                        'L23_pop',
                                        cell_id,
                                        num_cells,
                                        0,0,0,
                                        100,100,100)
    
          
    to_plot = {'Some_voltages':[]}
    to_save = {'%s_%s_voltages.dat'%(reference,cell_id):[]}
    
    interesting_seg_ids = [0,200,1000,2000,2500,2949] # [soma, .. some dends .. , axon]
    interesting_seg_ids = [0] # [soma, .. some dends .. , axon]

    pg0 = oc.add_pulse_generator(nml_doc,
                           id="pg0",
                           delay="200ms",
                           duration="600ms",
                           amplitude="0.4nA")
                           
    pg1 = oc.add_pulse_generator(nml_doc,
                           id="pg1",
                           delay="100ms",
                           duration="400ms",
                           amplitude="0.02nA")
    

    oc.add_targeted_inputs_to_population(network, 
                                         "PG_fixed",
                                         pop, 
                                         pg0.id, 
                                         segment_group='soma_group',
                                         number_per_cell = 1,
                                         all_cells=True)
                                         
    oc.add_targeted_inputs_to_population(network, 
                                         "PG_variable",
                                         pop, 
                                         'cond0',             # from ../../../ConductanceClamp.xml
                                         segment_group='soma_group',
                                         number_per_cell = 1,
                                         all_cells=True,
                                         weights = '1-2*random()')
                   
    for i in range(num_cells):
        
        for seg_id in interesting_seg_ids:
            to_plot.values()[0].append('%s/%s/%s/%s/v'%(pop.id, i, pop.component,seg_id))
            to_save.values()[0].append('%s/%s/%s/%s/v'%(pop.id, i, pop.component,seg_id))
            
                                

    nml_file_name = '%s.net.nml'%network.id + ('.h5' if format=='hdf5' else '')
    
    target_dir='./'
    
    oc.save_network(nml_doc, 
                    nml_file_name, 
                    validate=False, 
                    use_subfolder=True,
                    target_dir=target_dir,
                    format=format)



    lems_file_name, lems_sim = oc.generate_lems_simulation(nml_doc,
                                network, 
                                nml_file_name, 
                                duration, 
                                dt = 0.025,
                                target_dir=target_dir,
                                include_extra_lems_files=['../../../ConductanceClamp.xml'],
                                gen_plots_for_all_v = False,
                                plot_all_segments = False,
                                gen_plots_for_quantities = to_plot,   #  Dict with displays vs lists of quantity paths
                                gen_saves_for_all_v = False,
                                save_all_segments = False,
                                gen_saves_for_quantities = to_save,   #  Dict with file names vs lists of quantity paths
                                verbose = True)
                                
    if simulator:
                
        print ("Running %s for %sms in %s"%(lems_file_name, duration, simulator))

        traces, events = oc.simulate_network(lems_file_name,
                 simulator,
                 max_memory='4000M',
                 nogui=True,
                 load_saved_data=True,
                 reload_events=True,
                 plot=False,
                 verbose=True,
                 num_processors=min(num_cells,16))
                 
                                
                                
 
    
if __name__ == "__main__":
    
    cell_id = 'L23_NoHotSpot'
    cell_id = 'singleCompAllChans'
    cell_id = '../BBP/cADpyr229_L23_PC_5ecbf9b163_0_0'
    reference = "L23_Variability"
    duration = 1200
    
    simulator='jNeuroML_NEURON'
    #simulator='jNeuroML_NetPyNE'
    #simulator=None
        
    generate(cell_id, 
             duration=duration, 
             reference=reference, 
             num_cells = 10,
             format ='xml',
             simulator=simulator)
