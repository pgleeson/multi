

import opencortex.core as oc
from pyelectro import analysis
from pyneuroml import pynml


def get_rate_from_trace(times, volts):

    analysis_var={'peak_delta':0,'baseline':0,'dvdt_threshold':0, 'peak_threshold':0}

    try:
        analysis_data=analysis.IClampAnalysis(volts,
                                           times,
                                           analysis_var,
                                           start_analysis=0,
                                           end_analysis=times[-1],
                                           smooth_data=False,
                                           show_smoothed_data=False,
                                           verbose=True)

        analysed = analysis_data.analyse()

        import pprint; pp = pprint.PrettyPrinter()
        #pp.pprint(analysed)

        return analysed['mean_spike_frequency']
    
    except Exception as e:
        print("Problem getting rate: %s"%e)
        return 0

def generate(cell_id, duration, reference, 
             Bee=1,
             Ensyn = 10, 
             Erates = [50,100],
             st_onset = 0,
             st_duration = 1e9,
             format='hdf5',
             simulator=None,
             num_processors=1,
             target_group='soma_group',
             temperature='35degC'):

    #Insyn = int(Ensyn * 0.2)
    #bInsyn = int(bEnsyn * 0.2)
    
    cell_file = '../%s.cell.nml'%cell_id
    if '/' in cell_id:
        cell_id=cell_id.split('/')[-1]

    nml_doc, network = oc.generate_network(reference, temperature=temperature)

    oc.include_neuroml2_cell_and_channels(nml_doc,cell_file,cell_id)
    
    
    synAmpaEE = oc.add_exp_one_syn(nml_doc, id="synAmpaEE", gbase="%snS"%Bee,
                             erev="0mV", tau_decay="1ms")


    pop = oc.add_population_in_rectangular_region(network,
                                        'L23_pop',
                                        cell_id,
                                        len(Erates),
                                        0,0,0,
                                        100,100,100)
    
          
    
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



    lems_file_name, lems_sim = oc.generate_lems_simulation(nml_doc,
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
                 num_processors=num_processors)
                 
        rates = {}
        tt = [t*1000 for t in traces['t']]
        for tk in traces.keys():
            if tk!='t':
                rate = get_rate_from_trace(tt,[v*1000 for v in traces[tk]])
                print("Cell %s has rate %s Hz"%(tk,rate))
                i = int(tk.split('/')[1])
                rates[Erates[i]]=rate

        print Erates
        print rates

        ax = pynml.generate_plot([Erates],             
                                 [ [rates[r] for r in Erates] ],                
                                 "FF plots",               
                                 xaxis = 'Input frequency (Hz)',        
                                 yaxis = 'Firing frequency (Hz)',  
                                 markers=['o'],
                                 show_plot_already=True)     # Save figure
                                
        file_name = '%s.%s.%ssyns.%s.rates'%(cell_id,target_group,Ensyn,temperature)     
        f = open(file_name,'w')
        for r in Erates:
            f.write('%s\t%s\n'%(r,rates[r]))
        f.close()
        
        print("Finished! Saved rates data to %s"%file_name)
 
    
if __name__ == "__main__":
    
    cell_id = 'L23_NoHotSpot'
    #cell_id = '../BBP/cADpyr229_L23_PC_5ecbf9b163_0_0'
    #cell_id = 'singleCompAllChans'
    reference = "L23_FF"
    
    target_group = 'dendrite_group'
    target_group = 'soma_group'
    #target_group = 'soma_0'
    
    quick = True
    quick = False
    
    if quick:
        duration = 600
        Erates = range(5,1000,100)
    else:
        duration = 1000
        Erates = range(1,100,5)
    
    simulator='jNeuroML_NEURON'
    simulator='jNeuroML_NetPyNE'
    #simulator=None
    
    num_processors = min(16,len(Erates))
        
    generate(cell_id, 
             duration=duration, 
             reference=reference,
             Ensyn = 100,
             Erates=Erates, 
             format ='xml',
             simulator=simulator,
             num_processors=num_processors,
             target_group=target_group)
