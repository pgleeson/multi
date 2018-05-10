
files = ['temp/Sim_ISN_net.popExc.spikes','temp/Sim_ISN_net.popInh.spikes']

out = open('ISN-nest-EI-0.gdf','w')

for fn in files:
    f = open(fn)
    
    for l in f.readlines():
        w  = l.split('\t')
        i = int(w[0])
        t = float(w[1])*1000
        if 'Inh' in fn:
            i+=800
        
        out.write('%s   %s\n'%(t, i))
        
out.close()
        
        
