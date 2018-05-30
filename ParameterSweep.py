from ISN_tuned import run_one
from analysis_perturbation_pynn import analyse
import sys
from to_gdf import convert



class ParameterSweep():

    def __init__(self, vary, fixed={}):

        print("Initialising ParameterSweep with %s, %s" % (vary, fixed))
        self.fixed = fixed
        self.vary = vary
        self.complete = 0
        self.total_todo = 1
        for v in vary:
            self.total_todo *= len(vary[v])


    def _rem_key(self, d, key):
        r = dict(d)
        del r[key]
        return r


    def _run_instance(self, ** kwargs):

        print('============================================================= \n     Instance (%s/%s): %s' % (self.complete, self.total_todo, kwargs))
        '''     '''
        run_one( ** kwargs)
        convert(base_dir='./')
        analyse(kwargs['fraction_to_stim'], 1000, '%s.png' % kwargs['reference'], show=False)


    def _sweep(self, v, f, reference=''):

        #print("  VAR: <%s>\n  FIX <%s>"%(v,f))
        keys = list(v)

        if len(keys) > 1:
            vals = v[keys[0]]
            others = v
            others = self._rem_key(others, keys[0])
            for val in vals:
                all_params = f
                all_params[keys[0]] = val

                self._sweep(others, all_params, reference='%s-%s%s' % (reference, keys[0], val))

        else:
            vals = v[keys[0]]
            for val in vals:
                all_params = f
                all_params[keys[0]] = val
                r = '%s_%s%s' % (reference, keys[0], val)
                all_params['reference'] = 'REFb%s%s' % (self.complete, r)
                self._run_instance( ** all_params)
                self.complete += 1


    def run(self):

        print("Running")
        self._sweep(self.vary, self.fixed)


if __name__ == '__main__':

    if not sys.version_info[0] == 3:
        print('Run with Python 3...')
        quit()

    fixed = {'dt':0.025, 'fraction_to_stim':0.9, 'simulation_seed':100}

    quick = False
    #quick=True


    vary = {'Bee':[.55], 'Bei':[.5]}
    vary = {'Bee':[.5, .55, .6], 'Bei':[.5, .55], 'Bii':[.8, 1], 'Bie':[.8, 1, 1.2]}
    vary['r_bkg_ExtExc'] = [3600, 3800, 4000, 4200]
    vary['r_bkg_ExtExc'] = [4000, 4100, 4200]
    vary['r_bkg_ExtExc'] = [4000]
    vary['r_bkg_ExtInh'] = [3600, 3800, 4000]
    vary['r_bkg_ExtInh'] = [3600, 3700, 3800]
    vary['r_bkg_ExtInh'] = [3600]
    vary['r_stim'] = [-150, -160, -180, -200]
    vary['r_stim'] = [-180]

    if quick:
        vary = {'Bee':[.5], 'Bei':[.6]}
        vary['r_bkg_ExtExc'] = [4000]
        vary['r_bkg_ExtInh'] = [3500]
        vary['r_stim'] = [-160, -180]

    ps = ParameterSweep(vary, fixed)

    ps.run()
    '''
    vary['i'] = [1,2,3]

    ps = ParameterSweep(vary,fixed)

    ps.run()'''