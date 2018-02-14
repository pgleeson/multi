'''
Generates a NeuroML 2 file with many types of cells, populations and inputs
for testing purposes
'''

import opencortex
import opencortex.core as oc

import opencortex.utils.color as occ

import sys
import math

import numpy as np
import pylab as pl
from pyneuroml import pynml

from pyelectro import analysis
import pprint

pp = pprint.PrettyPrinter(indent=4)

min_pop_size = 1

exc_color = occ.L23_PRINCIPAL_CELL
exc2_color = occ.L23_PRINCIPAL_CELL_2
inh_color = occ.L23_INTERNEURON
inh_color = occ.L23_INTERNEURON_2

exc_color = '0 0 1'
exc2_color = '0 1 0'
inh_color = '1 0 0'
inh2_color = '1 0 1'


# transitent time to discard the data (ms)
Ttrans = 150.
# simulation time before perturbation (ms)
Tblank= 500.
# simulation time of perturbation (ms)
Tstim = 500.


def scale_pop_size(baseline, scale):
    return max(min_pop_size, int(baseline*scale))


def generate(scale_populations = 1,
             percentage_exc_detailed=0,
             exc2_cell = 'SmithEtAl2013/L23_Retuned_477127614',
             percentage_inh_detailed=0,
             scalex=1,
             scaley=1,
             scalez=1,
             Bee = .1,
             Bei = .1,
             Bie = -.2,
             Bii = -.2,
             Be_bkg = .1,
             Be_stim = .1,
             r_bkg = 0,
             r_stim = 0,
             percent_inh_pert=0.75,
             connections=True,
             exc_target_dendrites=False,
             inh_target_dendrites=False,
             duration = 1000,
             dt = 0.025,
             global_delay = .1,
             max_in_pop_to_plot_and_save = 10,
             format='xml',
             suffix='',
             run_in_simulator = None,
             num_processors = 1,
             target_dir='./temp/',
             exc_clamp=None):       # exc_clamp is work in progress...
                 
    reference = ("ISN_net%s"%(suffix)).replace('.','_')
                    

    num_exc = scale_pop_size(80,scale_populations)
    num_exc2  = int(math.ceil(num_exc*percentage_exc_detailed/100.0))
    num_exc -= num_exc2
    
    num_inh = scale_pop_size(20,scale_populations)
    num_inh2  = int(math.ceil(num_inh*percentage_inh_detailed/100.0))
    num_inh -= num_inh2
    
    nml_doc, network = oc.generate_network(reference, network_seed=1234)
    
    #exc_cell_id = 'AllenHH_480351780'
    exc_cell_id = 'AllenHH_477127614'
    #exc_cell_id = 'HH_477127614'
    exc_type = exc_cell_id.split('_')[0]
    oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/%s/%s.cell.nml'%(exc_type,exc_cell_id), exc_cell_id)
    
    
    #inh_cell_id = 'AllenHH_485058595'
    inh_cell_id = 'AllenHH_476686112'
    #inh_cell_id = 'HH_476686112'
    inh_type = exc_cell_id.split('_')[0]
    oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/%s/%s.cell.nml'%(inh_type,inh_cell_id), inh_cell_id)

    if percentage_exc_detailed>0:
        exc2_cell_id = exc2_cell.split('/')[1]
        exc2_cell_dir = exc2_cell.split('/')[0]
        oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/%s/%s.cell.nml'%(exc2_cell_dir,exc2_cell_id), exc2_cell_id)

    if percentage_inh_detailed>0:
        inh2_cell_id = 'cNAC187_L23_NBC_9d37c4b1f8_0_0'
        oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/BBP/%s.cell.nml'%inh2_cell_id, inh2_cell_id)
    

    xDim = 700*scalex
    yDim = 200*scaley
    zDim = 700*scalez

    xs = -1*xDim/2
    ys = -1*yDim/2
    zs = -1*zDim/2

    #####   Synapses
    

    synAmpaEE = oc.add_exp_one_syn(nml_doc, id="synAmpaEE", gbase="%snS"%Bee,
                             erev="0mV", tau_decay="1ms")
    synAmpaEI = oc.add_exp_one_syn(nml_doc, id="synAmpaEI", gbase="%snS"%Bei,
                             erev="0mV", tau_decay="1ms")

    synGabaIE = oc.add_exp_one_syn(nml_doc, id="synGabaIE", gbase="%snS"%abs(Bie),
                             erev="-80mV", tau_decay="1ms")
    synGabaII = oc.add_exp_one_syn(nml_doc, id="synGabaII", gbase="%snS"%abs(Bii),
                             erev="-80mV", tau_decay="1ms")

    synAmpaBkg = oc.add_exp_one_syn(nml_doc, id="synAmpaBkg", gbase="%snS"%Be_bkg,
                             erev="0mV", tau_decay="1ms")
    synAmpaStim = oc.add_exp_one_syn(nml_doc, id="synAmpaStim", gbase="%snS"%Be_stim,
                             erev="0mV", tau_decay="1ms")
                             
    #####   Input types


    tpfsA = oc.add_transient_poisson_firing_synapse(nml_doc,
                                       id="tpsfA",
                                       average_rate="%s Hz"%r_bkg,
                                       delay = '0ms', 
                                       duration = '%sms'%(Ttrans+Tblank),
                                       synapse_id=synAmpaBkg.id)

    tpfsB = oc.add_transient_poisson_firing_synapse(nml_doc,
                                       id="tpsfB",
                                       average_rate="%s Hz"%r_bkg,
                                       delay = '%sms'%(Ttrans+Tblank),
                                       duration = '%sms'%(Tstim),
                                       synapse_id=synAmpaBkg.id)

    tpfsC = oc.add_transient_poisson_firing_synapse(nml_doc,
                                       id="tpsfC",
                                       average_rate="%s Hz"%(r_bkg+r_stim),
                                       delay = '%sms'%(Ttrans+Tblank),
                                       duration = '%sms'%(Tstim),
                                       synapse_id=synAmpaStim.id)


    #####   Populations

    popExc = oc.add_population_in_rectangular_region(network,
                                                  'popExc',
                                                  exc_cell_id,
                                                  num_exc,
                                                  xs,ys,zs,
                                                  xDim,yDim,zDim,
                                                  color=exc_color)                 
    allExc = [popExc]

    if num_exc2>0:
        popExc2 = oc.add_population_in_rectangular_region(network,
                                                  'popExc2',
                                                  exc2_cell_id,
                                                  num_exc2,
                                                  xs,ys,zs,
                                                  xDim,yDim,zDim,
                                                  color=exc2_color)
                                                  
        allExc.append(popExc2)

    popInh = oc.add_population_in_rectangular_region(network,
                                                  'popInh',
                                                  inh_cell_id,
                                                  num_inh,
                                                  xs,ys,zs,
                                                  xDim,yDim,zDim,
                                                  color=inh_color)           
    allInh = [popInh]
    
    if num_inh2>0:
        popInh2 = oc.add_population_in_rectangular_region(network,
                                                  'popInh2',
                                                  inh2_cell_id,
                                                  num_inh2,
                                                  xs,ys,zs,
                                                  xDim,yDim,zDim,
                                                  color=inh2_color)
                                                  
        allInh.append(popInh2)


    #####   Projections

    if connections:

        exc_exc_conn_prob = 0.15
        exc_inh_conn_prob = 0.15
        inh_exc_conn_prob = 1
        inh_inh_conn_prob = 1
        
        for popEpr in allExc:
            
            for popEpo in allExc:
                proj = add_projection(network, "projEE",
                                      popEpr, popEpo,
                                      synAmpaEE.id, exc_exc_conn_prob, 
                                      global_delay,
                                      exc_target_dendrites)
                                                
            for popIpo in allInh:
                proj = add_projection(network, "projEI",
                                      popEpr, popIpo,
                                      synAmpaEI.id, exc_inh_conn_prob, 
                                      global_delay,
                                      exc_target_dendrites)

            
        for popIpr in allInh:
            
            for popEpo in allExc:
                proj = add_projection(network, "projIE",
                                      popIpr, popEpo,
                                      synGabaIE.id, inh_exc_conn_prob, 
                                      global_delay,
                                      inh_target_dendrites)
        
            for popIpo in allInh:
                proj = add_projection(network, "projII",
                                      popIpr, popIpo,
                                      synGabaII.id, inh_inh_conn_prob, 
                                      global_delay,
                                      inh_target_dendrites)

                                        

    #####   Inputs

    for pop in allExc:
        oc.add_inputs_to_population(network, "Stim_pre_%s"%pop.id,
                                    pop, tpfsA.id,
                                    all_cells=True)
    for pop in allInh:
        oc.add_inputs_to_population(network, "Stim_pre_%s"%pop.id,
                                    pop, tpfsA.id,
                                    all_cells=True)
                                 
    for pop in allExc:
        oc.add_inputs_to_population(network, "Stim_E_%s"%pop.id,
                                    pop, tpfsB.id,
                                    all_cells=True)   
                        
    for pop in allInh:
        
        num_inh_pert = int(pop.get_size()*percent_inh_pert)
           
        oc.add_inputs_to_population(network, "Stim_I_pert_%s"%pop.id,
                                    pop, tpfsC.id,
                                    all_cells=False,
                                    only_cells=range(0,num_inh_pert))   
                                    
        oc.add_inputs_to_population(network, "Stim_I_nonpert_%s"%pop.id,
                                    pop, tpfsB.id,
                                    all_cells=False,
                                    only_cells=range(num_inh_pert,pop.get_size()))  
    


    # Work in progress...
    # General idea: clamp one (or more) exc cell at rev pot of inh syn and see only exc inputs
    #
    if exc_clamp:
        
        vc = oc.add_voltage_clamp_triple(nml_doc, id="exc_clamp", 
                             delay='0ms', 
                             duration='%sms'%duration, 
                             conditioning_voltage=synGaba1.erev,
                             testing_voltage=synGaba1.erev,
                             return_voltage=synGaba1.erev, 
                             simple_series_resistance="1e5ohm",
                             active = "1")
                             
        for pop in exc_clamp:
            oc.add_inputs_to_population(network, "exc_clamp_%s"%pop,
                                        network.get_by_id(pop), vc.id,
                                        all_cells=False,
                                        only_cells=exc_clamp[pop])
                


    #####   Save NeuroML and LEMS Simulation files      
    

    nml_file_name = '%s.net.%s'%(network.id,'nml.h5' if format == 'hdf5' else 'nml')
    oc.save_network(nml_doc, 
                    nml_file_name, 
                    validate=(format=='xml'),
                    format = format,
                    target_dir=target_dir)
        

    if format=='xml':
        
        
        save_v = {}
        plot_v = {}
        
        if num_exc>0:
            exc_traces = '%s_%s_v.dat'%(network.id,popExc.id)
            save_v[exc_traces] = []
            plot_v[popExc.id] = []
            
        if num_inh>0:
            inh_traces = '%s_%s_v.dat'%(network.id,popInh.id)
            save_v[inh_traces] = []
            plot_v[popInh.id] = []
            
        if num_exc2>0:
            exc2_traces = '%s_%s_v.dat'%(network.id,popExc2.id)
            save_v[exc2_traces] = []
            plot_v[popExc2.id] = []
            
        if num_inh2>0:
            inh2_traces = '%s_%s_v.dat'%(network.id,popInh2.id)
            save_v[inh2_traces] = []
            plot_v[popInh2.id] = []
        
        
        for i in range(min(max_in_pop_to_plot_and_save,num_exc)):
            plot_v[popExc.id].append("%s/%i/%s/v"%(popExc.id,i,popExc.component))
            save_v[exc_traces].append("%s/%i/%s/v"%(popExc.id,i,popExc.component))
            
        for i in range(min(max_in_pop_to_plot_and_save,num_exc2)):
            plot_v[popExc2.id].append("%s/%i/%s/v"%(popExc2.id,i,popExc2.component))
            save_v[exc2_traces].append("%s/%i/%s/v"%(popExc2.id,i,popExc2.component))
            
        for i in range(min(max_in_pop_to_plot_and_save,num_inh)):
            plot_v[popInh.id].append("%s/%i/%s/v"%(popInh.id,i,popInh.component))
            save_v[inh_traces].append("%s/%i/%s/v"%(popInh.id,i,popInh.component))
            
        for i in range(min(max_in_pop_to_plot_and_save,num_inh2)):
            plot_v[popInh2.id].append("%s/%i/%s/v"%(popInh2.id,i,popInh2.component))
            save_v[inh2_traces].append("%s/%i/%s/v"%(popInh2.id,i,popInh2.component))
            
        gen_spike_saves_for_all_somas = run_in_simulator!='jNeuroML_NetPyNE'
            
        lems_file_name = oc.generate_lems_simulation(nml_doc, network, 
                                target_dir+nml_file_name, 
                                duration =      duration, 
                                dt =            dt,
                                gen_plots_for_all_v = False,
                                gen_plots_for_quantities = plot_v,
                                gen_saves_for_all_v = False,
                                gen_saves_for_quantities = save_v,
                                gen_spike_saves_for_all_somas = gen_spike_saves_for_all_somas,
                                target_dir=target_dir)
                                
        
        if run_in_simulator:
            
            print ("Running %s for %sms in %s"%(lems_file_name, duration, run_in_simulator))
            
            traces, events = oc.simulate_network(lems_file_name,
                     run_in_simulator,
                     max_memory='4000M',
                     nogui=True,
                     load_saved_data=True,
                     reload_events=True,
                     plot=False,
                     verbose=False,
                     num_processors=num_processors)
                     
                     
            print("Reloaded traces: %s"%traces.keys())
            #print("Reloaded events: %s"%events.keys())
            
            use_events_for_rates = False
            
            exc_rate = 0
            inh_rate = 0
            
            if use_events_for_rates:
                if (run_in_simulator=='jNeuroML_NetPyNE'):
                    raise('Saving of spikes (and so calculation of rates) not yet supported in jNeuroML_NetPyNE')
                for ek in events.keys():
                    rate = 1000 * len(events[ek])/float(duration)
                    print("Cell %s has rate %s Hz"%(ek,rate))
                    if 'popExc' in ek:
                        exc_rate += rate/num_exc
                    if 'popInh' in ek:
                        inh_rate += rate/num_inh
            
            else:
                tot_exc_rate = 0 
                exc_cells = 0
                tot_inh_rate = 0 
                inh_cells = 0
                tt = [t*1000 for t in traces['t']]
                for tk in traces.keys():
                    if tk!='t':
                        rate = get_rate_from_trace(tt,[v*1000 for v in traces[tk]])
                        print("Cell %s has rate %s Hz"%(tk,rate))
                        if 'popExc' in tk:
                            tot_exc_rate += rate
                            exc_cells+=1
                        if 'popInh' in tk:
                            tot_inh_rate += rate
                            inh_cells+=1
                            
                exc_rate = tot_exc_rate/exc_cells
                inh_rate = tot_inh_rate/inh_cells
                    
                    
                    
            print("Run %s: Exc rate: %s Hz; Inh rate %s Hz"%(reference,exc_rate, inh_rate))
                     
            return exc_rate, inh_rate, traces
        
    else:
        lems_file_name = None
                                
    return nml_doc, nml_file_name, lems_file_name
                        
                        
def add_projection(network, 
                   proj_id,
                   pop_pre, 
                   pop_post,
                   syn_id,
                   conn_prob,
                   delay,
                   target_dendrites):


    if pop_post.size > pop_pre.size:
        num_connections = pop_pre.size * conn_prob
        targeting_mode='convergent'
    else:
        num_connections = pop_post.size * conn_prob
        targeting_mode='divergent'

    post_segment_group = 'soma_group'
    
    if '2' in pop_post.id and target_dendrites:
        post_segment_group = 'dendrite_group'
        

    proj = oc.add_targeted_projection(network,
                                    proj_id,
                                    pop_pre,
                                    pop_post,
                                    targeting_mode=targeting_mode,
                                    synapse_list=[syn_id],
                                    pre_segment_group = 'soma_group',
                                    post_segment_group = post_segment_group,
                                    number_conns_per_cell=num_connections,
                                    delays_dict = {syn_id:delay})
    return proj
       
def get_rate_from_trace(times, volts):

    analysis_var={'peak_delta':0,'baseline':0,'dvdt_threshold':0, 'peak_threshold':0}

    try:
        analysis_data=analysis.IClampAnalysis(volts,
                                           times,
                                           analysis_var,
                                           start_analysis=0,
                                           end_analysis=times[-1],
                                           smooth_data=False,
                                           show_smoothed_data=False)

        analysed = analysis_data.analyse()

        pp.pprint(analysed)

        return analysed['mean_spike_frequency']
    
    except:
        return 0

                         
def _plot_(X, g_rng, i_rng, sbplt=111, ttl=[]):
    ax = pl.subplot(sbplt)
    pl.title(ttl)
    pl.imshow(X, origin='lower', interpolation='none')
    pl.xlabel('Ratio inh/exc')
    pl.ylabel('Input (Hz)')
    ax.set_xticks(range(0,len(g_rng))); ax.set_xticklabels(g_rng)
    ax.set_yticks(range(0,len(i_rng))); ax.set_yticklabels(i_rng)
    pl.colorbar()


if __name__ == '__main__':
    
    if '-all' in sys.argv:
        generate()
        
        generate(scale_populations = 5,
             scalex=2,
             scalez=2,
             connections=False)
             
             
    elif '-paramSweep' in sys.argv:     
        
        run_in_simulator='jNeuroML_NEURON'
        run_in_simulator='jNeuroML_NetPyNE'
        num_processors = 18
        
        duration = 1000
        dt = 0.025
        scale_populations = 0.2
        percentage_exc_detailed = 0.01
        percentage_exc_detailed = 100
        #percentage_exc_detailed = 0
        percentage_inh_detailed = 0
        #percentage_inh_detailed = 100
        
        quick = False
        quick = True
        
        g_rng = np.arange(.5, 4.5, .5)
        i_rng = np.arange(50, 400, 50)
        trace_highlight = [(2,150)]
        
        if quick:
            g_rng = [2]
            #g_rng = [2,3,4]
            i_rng = [350]
            #i_rng = [250,300]
            #i_rng = [100,150,200]
            trace_highlight = [(g_rng[0],i_rng[0])]
            
            duration = 1000
            scale_populations = .1
            percentage_exc_detailed = 100
            #percentage_exc_detailed = 0.01
            #percentage_exc_detailed = 0
            percentage_inh_detailed = 100
            #percentage_inh_detailed = 0
            run_in_simulator='jNeuroML_NEURON'
            run_in_simulator='jNeuroML_NetPyNE'
            num_processors = 12


        Rexc = np.zeros((len(g_rng), len(i_rng)))
        Rinh = np.zeros((len(g_rng), len(i_rng)))
        
        desc = '%s_%s_%sms_%se2_%si2'%(run_in_simulator,scale_populations, duration,percentage_exc_detailed,percentage_inh_detailed)
        
        count=1
        for i1, g in enumerate(g_rng):
            for i2, i in enumerate(i_rng):
                print("====================================")
                highlight = False
                for h in trace_highlight:
                    if h[0]==g and h[1]==i:
                        highlight = True
                        
                print(" Run %s of %s: scale=%s; g=%s; i=%s (highlighting: %s)"%(count, len(g_rng)*len(i_rng), scale_populations, g, i, highlight))
                
                info = generate(scale_populations = scale_populations,
                    scalex=2,
                    scalez=2,
                    duration = duration,
                    dt = dt,
                    max_in_pop_to_plot_and_save = 5,
                    percentage_exc_detailed = percentage_exc_detailed,
                    percentage_inh_detailed = percentage_inh_detailed,
                    global_delay = 2,
                    ratio_inh_exc = g,
                    input_rate=i,
                    run_in_simulator=run_in_simulator,
                    num_processors = num_processors )
                    
                Rexc[i1,i2] = info[0]
                Rinh[i1,i2] = info[1]

                if highlight:
                    traces = info[2]
                    all_t = []
                    all_v = []
                    colours = []
                    tr_shade_e=1
                    tr_shade_e2=1
                    tr_shade_i=1
                    tr_shade_i2=1
                    trace_keys_ordered = []              
                    for vs in traces.keys():
                        if 'Exc2' in vs or 'Inh2' in vs:
                            trace_keys_ordered.insert(0,vs)
                        else:
                            trace_keys_ordered.append(vs)
                            
                    for vs in trace_keys_ordered:
                        if vs!='t':
                            all_v.append([v*1000.0 for v in traces[vs]])
                            all_t.append([t*1000.0 for t in traces['t']])
                            if 'Exc2' in vs:
                                colours.append((1-tr_shade_e2,1,1-tr_shade_e2))
                                tr_shade_e2*=0.8
                            elif 'Exc' in vs:
                                colours.append((1-tr_shade_e,1-tr_shade_e,1))
                                tr_shade_e*=0.8
                            elif 'Inh2' in vs:
                                colours.append((1,1-tr_shade_i,1))
                                tr_shade_e*=0.8
                            else:
                                colours.append((1,1-tr_shade_i,1-tr_shade_i))
                                tr_shade_i*=0.8
                                
                    
                    print(colours)
                    pynml.generate_plot(all_t, all_v, 
                                        "Sim g=%s, i=%s"%(g,i),
                                        colors=colours,
                                        show_plot_already=False,
                                        xaxis = 'Time (ms)',            # x axis legend
                                        yaxis = 'Membrane potential (mV)',   # y axis legend
                                        save_figure_to='%s_traces.png'%desc)
                count+=1
                    
                

        fig = pl.figure(figsize=(16,8))
        info = "%s: scale %s, %s ms"%(run_in_simulator,scale_populations, duration)

        fig.canvas.set_window_title(info)
        pl.suptitle(info)

        _plot_(Rexc.T, g_rng, i_rng, 221, 'Rates Exc (Hz)')
        _plot_(Rinh.T, g_rng, i_rng, 222, 'Rates Inh (Hz)')
        
        

        pl.subplots_adjust(wspace=.3, hspace=.3)


        pl.savefig("%s_rates.png"%desc, bbox_inches='tight')
        print("Finished: "+info)
        pl.show()
        
    elif '-standard' in sys.argv:   
        
        opencortex.set_verbose()
        '''
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="A",
                 percentage_exc_detailed=0,
                 target_dir='./NeuroML2/')
     
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="B",
                 percentage_exc_detailed=0.1,
                 exc2_cell='Thalamocortical/L23PyrRS',
                 exc_target_dendrites=True,
                 target_dir='./NeuroML2/')
     
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="C",
                 percentage_exc_detailed=0.1,
                 exc2_cell='SmithEtAl2013/L23_NoHotSpot',
                 exc_target_dendrites=True,
                 target_dir='./NeuroML2/')'''
     
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="D",
                 percentage_exc_detailed=0.1,
                 exc2_cell='BBP/cADpyr229_L23_PC_5ecbf9b163_0_0',
                 exc_target_dendrites=True,
                 target_dir='./NeuroML2/')
        '''
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="B2",
                 exc_target_dendrites=True,
                 inh_target_dendrites=True,
                 percentage_exc_detailed=0.1,
                 target_dir='./NeuroML2/')
        
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="C",
                 percentage_exc_detailed=100,
                 target_dir='./NeuroML2/')
        
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scale_populations=1,
                 suffix="D",
                 percentage_exc_detailed=100,
                 percentage_inh_detailed=100,
                 target_dir='./NeuroML2/')'''

                     
        
    else:
        # background rate (sp/s)
        r_bkg = 10000.-400.
        # rate of perturbation (sp/s)
        r_stim = -400.
        
        percent_inh_pert = 0.75
        
        Be = .1
        Bi = -.2*1
        Be_bkg = .1
        Be_stim = .1
        
        scale_up = 5
        Be*=scale_up
        Bi*=scale_up
        Be_bkg*=scale_up
        Be_stim*=scale_up
        
        Bee = Be
        Bei = Be
        Bie = Bi
        Bii = Bi
        
        if '-test' in sys.argv:   
            #r_bkg = 20
            r_stim = -2000
            
            percent_inh_pert = .75
            scale_populations = 2
            
        
        generate(Bee = Bee,
                 Bei = Bei,
                 Bie = Bie,
                 Bii = Bii,
                 Be_bkg = Be_bkg,
                 Be_stim = Be_stim,
                 r_bkg = r_bkg,
                 r_stim = r_stim,
                 percent_inh_pert=percent_inh_pert,
                 duration = 1000,
                 dt = 0.025,
                 scale_populations=scale_populations,
                 format='xml',
                 percentage_exc_detailed=0,
                 target_dir='./temp/')