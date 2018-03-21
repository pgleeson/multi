

import opencortex.core as oc


def generate(cell_id, duration, reference, 
             Bee=1,
             Ensyn = 10, 
             Erates = [50,100],
             st_onset = 0,
             st_duration = 1e9,
             format='hdf5'):

    #Insyn = int(Ensyn * 0.2)
    #bInsyn = int(bEnsyn * 0.2)
    
    cell_file = '../%s.cell.nml'%cell_id

    nml_doc, network = oc.generate_network(reference, temperature='35degC')

    oc.include_neuroml2_cell_and_channels(nml_doc,cell_file,cell_id)
    
    
    synAmpaEE = oc.add_exp_one_syn(nml_doc, id="synAmpaEE", gbase="%snS"%Bee,
                             erev="0mV", tau_decay="1ms")




    pop = oc.add_population_in_rectangular_region(network,
                                        'L23_pop',
                                        cell_id,
                                        len(Erates),
                                        0,0,0,
                                        100,100,100)
    
          
    target_group = 'dendrite_group'
    target_group = 'soma_group'
    
    to_plot = {'Some_voltages':[]}
    to_save = {'%s_voltages.dat'%cell_id:[]}
    
    interesting_seg_ids = [0,200,1000,2000,2500,2949] # [soma, .. some dends .. , axon]
    interesting_seg_ids = [0] # [soma, .. some dends .. , axon]

    
    for i,r in enumerate(Erates):
        
        syn0 = oc.add_transient_poisson_firing_synapse(nml_doc,
                                           id="%s_stim_%s"%(synAmpaEE.id,r),
                                           average_rate="%s Hz"%r,
                                           synapse_id=synAmpaEE.id,
                                           delay='%s ms'%st_onset,
                                           duration='%s ms'%st_duration)
                                           
        oc.add_targeted_inputs_to_population(network, 
                                             "Esyn_%s"%r,
                                             pop, 
                                             syn0.id, 
                                             segment_group=target_group,
                                             number_per_cell = Ensyn,
                                             all_cells=False,
                                             only_cells=[i])
                                             
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



    oc.generate_lems_simulation(nml_doc,
                                network, 
                                nml_file_name, 
                                duration, 
                                dt = 0.025,
                                target_dir=target_dir,
                                gen_plots_for_all_v = False,
                                plot_all_segments = False,
                                gen_plots_for_quantities = to_plot,   #  Dict with displays vs lists of quantity paths
                                gen_saves_for_all_v = False,
                                save_all_segments = False,
                                gen_saves_for_quantities = to_save,   #  Dict with file names vs lists of quantity paths
                                verbose = True)
 
    
if __name__ == "__main__":
    
    cell_id = 'L23_NoHotSpot'
    cell_id = 'singleCompAllChans'
    reference = "L23_FF"
    duration = 600
    Erates = range(50,500,50)
        
    generate(cell_id, duration=duration, reference=reference, Erates=Erates, format ='xml')
