'''
Generates a NeuroML 2 file with many types of cells, populations and inputs
for testing purposes
'''

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

min_pop_size = 3

exc_color = occ.L23_PRINCIPAL_CELL
exc2_color = occ.L23_PRINCIPAL_CELL
inh_color = occ.L23_INTERNEURON

exc_color = '0 0 1'
exc2_color = '0 1 0'
inh_color = '1 0 0'


def scale_pop_size(baseline, scale):
    return max(min_pop_size, int(baseline*scale))


def generate(scalePops = 1,
             percentage_exc_detailed=0,
             scalex=1,
             scaley=1,
             scalez=1,
             ratio_inh_exc=2,
             connections=True,
             duration = 1000,
             input_rate = 150,
             global_delay = 2,
             max_in_pop_to_plot_and_save = 5,
             format='xml',
             suffix='',
             run_in_simulator = None,
             num_processors = 1,
             target_dir='./temp/'):
                 
    reference = ("Multiscale__g%s__i%s%s"%(ratio_inh_exc,input_rate,suffix)).replace('.','_')
                    

    num_exc = scale_pop_size(80,scalePops)
    num_exc2  = int(math.ceil(num_exc*percentage_exc_detailed/100.0))
    num_exc -= num_exc2
    num_inh = scale_pop_size(40,scalePops)
    
    nml_doc, network = oc.generate_network(reference)
    
    #exc_cell_id = 'AllenHH_479704527'
    #oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/AllenHH/AllenHH_479704527.cell.nml', exc_cell_id)
    exc_cell_id = 'AllenHH_480351780'
    oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/AllenHH/AllenHH_480351780.cell.nml', exc_cell_id)
    
    
    inh_cell_id = 'AllenHH_485058595'
    oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/AllenHH/AllenHH_485058595.cell.nml', inh_cell_id)

    #oc.include_opencortex_cell(nml_doc, 'AllenInstituteCellTypesDB_HH/HH_477127614.cell.nml')
    #oc.include_opencortex_cell(nml_doc, 'AllenInstituteCellTypesDB_HH/HH_476686112.cell.nml')
    if percentage_exc_detailed>0:
        exc2_cell_id = 'L23_Retuned'
        oc.include_neuroml2_cell_and_channels(nml_doc, 'cells/SmithEtAl2013/L23_Retuned.cell.nml', exc2_cell_id)
        #oc.include_opencortex_cell(nml_doc, 'L23Pyr_SmithEtAl2013/L23_NoHotSpot.cell.nml')
    

    xDim = 700*scalex
    yDim = 200*scaley
    zDim = 700*scalez

    xs = -1*xDim/2
    ys = -1*yDim/2
    zs = -1*zDim/2

    #####   Synapses
    
    exc_syn_nS = 1.

    synAmpa1 = oc.add_exp_two_syn(nml_doc, id="synAmpa1", gbase="%snS"%exc_syn_nS,
                             erev="0mV", tau_rise="0.5ms", tau_decay="5ms")

    synGaba1 = oc.add_exp_two_syn(nml_doc, id="synGaba1", gbase="%snS"%(exc_syn_nS*ratio_inh_exc),
                             erev="-80mV", tau_rise="1ms", tau_decay="20ms")

    #####   Input types


    pfs1 = oc.add_poisson_firing_synapse(nml_doc,
                                       id="psf1",
                                       average_rate="%s Hz"%input_rate,
                                       synapse_id=synAmpa1.id)


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


    #####   Projections

    if connections:

        for pop1 in allExc:
            
            for pop2 in allExc:
                proj = oc.add_probabilistic_projection(network, "proj0",
                                                pop1, pop2,
                                                synAmpa1.id, 0.5, delay = global_delay)

            proj = oc.add_probabilistic_projection(network, "proj1",
                                            pop1, popInh,
                                            synAmpa1.id, 0.7, delay = global_delay)

            proj = oc.add_probabilistic_projection(network, "proj2",
                                            popInh, pop1,
                                            synGaba1.id, 0.7, delay = global_delay)

        proj = oc.add_probabilistic_projection(network, "proj3",
                                        popInh, popInh,
                                        synGaba1.id, 0.5, delay = global_delay)

                                        

    #####   Inputs

    for pop in allExc:
        oc.add_inputs_to_population(network, "Stim_%s"%pop.id,
                                    pop, pfs1.id,
                                    all_cells=True)



    #####   Save NeuroML and LEMS Simulation files      
    

    nml_file_name = '%s.net.%s'%(network.id,'nml.h5' if format == 'hdf5' else 'nml')
    oc.save_network(nml_doc, 
                    nml_file_name, 
                    validate=(format=='xml'),
                    format = format,
                    target_dir=target_dir)
        

    if format=='xml':
        
        plot_v = {popExc.id:[],popInh.id:[]}
        
        exc_traces = '%s_%s_v.dat'%(network.id,popExc.id)
        inh_traces = '%s_%s_v.dat'%(network.id,popInh.id)
        save_v = {exc_traces:[], inh_traces:[]}
        
        if num_exc2>0:
            exc2_traces = '%s_%s_v.dat'%(network.id,popExc2.id)
            save_v = {exc_traces:[], inh_traces:[], exc2_traces:[]}
            plot_v = {popExc.id:[],popExc2.id:[],popInh.id:[]}
        
        
        for i in range(min(max_in_pop_to_plot_and_save,num_exc)):
            plot_v[popExc.id].append("%s/%i/%s/v"%(popExc.id,i,popExc.component))
            save_v[exc_traces].append("%s/%i/%s/v"%(popExc.id,i,popExc.component))
            
        for i in range(min(max_in_pop_to_plot_and_save,num_exc2)):
            plot_v[popExc2.id].append("%s/%i/%s/v"%(popExc2.id,i,popExc2.component))
            save_v[exc2_traces].append("%s/%i/%s/v"%(popExc2.id,i,popExc2.component))
            
        for i in range(min(max_in_pop_to_plot_and_save,num_inh)):
            plot_v[popInh.id].append("%s/%i/%s/v"%(popInh.id,i,popInh.component))
            save_v[inh_traces].append("%s/%i/%s/v"%(popInh.id,i,popInh.component))
            
        gen_spike_saves_for_all_somas = run_in_simulator!='jNeuroML_NetPyNE'
            
        lems_file_name = oc.generate_lems_simulation(nml_doc, network, 
                                target_dir+nml_file_name, 
                                duration =      duration, 
                                dt =            0.025,
                                gen_plots_for_all_v = False,
                                gen_plots_for_quantities = plot_v,
                                gen_saves_for_all_v = False,
                                gen_saves_for_quantities = save_v,
                                gen_spike_saves_for_all_somas = gen_spike_saves_for_all_somas,
                                target_dir=target_dir)
                                
        
        if run_in_simulator:
            
            print ("Running %s in %s"%(lems_file_name, run_in_simulator))
            
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
        
        generate(scalePops = 5,
             scalex=2,
             scalez=2,
             connections=False)
             
             
    elif '-paramSweep' in sys.argv:     
        
        run_in_simulator='jNeuroML_NEURON'
        run_in_simulator='jNeuroML_NetPyNE'
        num_processors = 18
        
        duration = 500
        scalePops = .5
        percentage_exc_detailed = 0
        
        quick = False
        quick = True
        
        g_rng = np.arange(.5, 4.5, .5)
        i_rng = np.arange(50, 400, 50)
        trace_highlight = [(2,150)]
        
        if quick:
            g_rng = [2]
            #g_rng = [2,3,4]
            i_rng = [250]
            #i_rng = [250,300]
            #i_rng = [100,150,200]
            trace_highlight = [(g_rng[0],i_rng[0])]
            
            duration = 1000
            scalePops = .25
            percentage_exc_detailed = 100


        Rexc = np.zeros((len(g_rng), len(i_rng)))
        Rinh = np.zeros((len(g_rng), len(i_rng)))
        
        desc = '%s_%s_%sms_%se2'%(run_in_simulator,scalePops, duration,percentage_exc_detailed)
        
        count=1
        for i1, g in enumerate(g_rng):
            for i2, i in enumerate(i_rng):
                print("====================================")
                highlight = False
                for h in trace_highlight:
                    if h[0]==g and h[1]==i:
                        highlight = True
                print(" Run %s of %s: scale=%s; g=%s; i=%s (highlighting: %s)"%(count, len(g_rng)*len(i_rng), scalePops, g, i, highlight))
                info = generate(scalePops = scalePops,
                    scalex=2,
                    scalez=2,
                    duration = duration,
                    max_in_pop_to_plot_and_save = 5,
                    percentage_exc_detailed = percentage_exc_detailed,
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
                    for vs in traces.keys():
                        if vs!='t':
                            all_v.append([v*1000.0 for v in traces[vs]])
                            all_t.append([t*1000.0 for t in traces['t']])
                            if 'Exc2' in vs:
                                colours.append((1-tr_shade_e2,1,1-tr_shade_e2))
                                tr_shade_e2*=0.8
                            elif 'Exc' in vs:
                                colours.append((1-tr_shade_e,1-tr_shade_e,1))
                                tr_shade_e*=0.8
                            else:
                                colours.append((1,1-tr_shade_i,1-tr_shade_i))
                                tr_shade_i*=0.8
                                
                    
                    print colours
                    pynml.generate_plot(all_t, all_v, 
                                        "Sim g=%s, i=%s"%(g,i),
                                        colors=colours,
                                        show_plot_already=False,
                                        xaxis = 'Time (ms)',            # x axis legend
                                        yaxis = 'Membrane potential (mV)',   # y axis legend
                                        save_figure_to='%s_traces.png'%desc)
                count+=1
                    
                

        fig = pl.figure(figsize=(16,8))
        info = "%s: scale %s, %s ms"%(run_in_simulator,scalePops, duration)

        fig.canvas.set_window_title(info)
        pl.suptitle(info)

        _plot_(Rexc.T, g_rng, i_rng, 221, 'Rates Exc (Hz)')
        _plot_(Rinh.T, g_rng, i_rng, 222, 'Rates Inh (Hz)')
        
        

        pl.subplots_adjust(wspace=.3, hspace=.3)


        pl.savefig("%s_rates.png"%desc, bbox_inches='tight')
        print("Finished: "+info)
        pl.show()
        
    elif '-test' in sys.argv:   
        
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scalePops=1,
                 suffix="A",
                 percentage_exc_detailed=0,
                 target_dir='./NeuroML2/')
        
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scalePops=1,
                 suffix="B",
                 percentage_exc_detailed=0.1,
                 target_dir='./NeuroML2/')
        
        generate(ratio_inh_exc=1.5,
                 duration = 1000,
                 input_rate = 250,
                 scalePops=1,
                 suffix="C",
                 percentage_exc_detailed=100,
                 target_dir='./NeuroML2/')

    else:
        generate(ratio_inh_exc=1.5,
                 duration = 500,
                 input_rate = 250,
                 scalePops=1,
                 percentage_exc_detailed=0.1,
                 target_dir='./NeuroML2/')