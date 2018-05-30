import os

def convert(base_dir='./temp/'):
    files = ['%sSim_ISN_net.popExc.spikes'%base_dir,'%sSim_ISN_net.popInh.spikes'%base_dir]

    out_file = 'ISN-nest-EI-0.gdf'
    out = open(out_file,'w')

    for fn in files:
        f = open(fn)

        for l in f.readlines():
            w  = l.split('\t')
            i = int(w[0])
            t = float(w[1])*1000
            if 'Inh' in fn:
                i+=800

            out.write('%s   %s\n'%(t, i))
            
        print("> Converted: %s"%os.path.abspath(fn))
        
    print("> Saved:     %s"%os.path.abspath(out_file))
        
    out.close()


if __name__ == '__main__':
    
    convert()
        
        
