from __future__ import division
import numpy as np
from scipy.io import loadmat, savemat
from scipy.stats import scoreatpercentile
import matplotlib.pyplot as plt
import os
from os.path import join, splitext, isfile, basename

from pyhsmm.models import WeakLimitStickyHDPHMM
from pyhsmm.basic.distributions import Gaussian
from pyhsmm.util.stats import cov

from config import getdirs


class Trace(object):
    '''
    A container object for trace experiments.
    '''

    def __init__(self,matfilepath,start=0,end=None):
        self.matfilepath = matfilepath
        self.start = int(start)
        self.end = int(end) if end is not None else None
        self.name = _name_from_matfilepath(matfilepath)

        # load matfile
        mat = loadmat(matfilepath)

        # add model_data
        self.model_data = \
            np.hstack((mat['unsmoothedGrI'].T,mat['unsmoothedRedI'].T)
                      )[self.start:self.end]

        # we start out with no model
        self.model = None

        # add everything from the matlab dict to this instance
        for k, v in mat.items():
            if not k.startswith('__'):
                self.__dict__[k] = np.squeeze(v)

    ### model fitting

    def run_gibbs(self,niter):
        if self.model is None:
            self.reset_model()

        for itr in xrange(niter):
            self.model.resample_model()

        self.model.states_list[0].Viterbi()

    def reset_model(self):
        self.model = self._make_model()
        self.model.add_data(self.model_data)

    def _make_model(self):
        data = self.model_data
        # The parameters are:
        # Gaussian observation distributions (ellipses in red-green intensity space)
        # mu_0 and sigma parameterize our prior belief about the means and sigma of each state
        # nu_0 expresses your confidence in the prior--it's the number of data points that
        # you claim got you these prior parameters. Nu_0 has to be strictly bigger than the
        # number of dimensions (2, in our case). You could do 2.01.
        # The nominal covariance is sigma_0/nu_0, so hence the 3 in sigma_0.
        # kappa_0: Uncertainty in the mean should be related to uncertainty in the covariance.
        # kappa_0 is an amplitude for that. Smaller number means other states' means will be
        # further away.
        obs_hypparams = dict(
            mu_0=data.mean(0),
            sigma_0=3.*cov(data),
            nu_0=3.,
            kappa_0=0.5,
        )

        # In the function call below:
        # (1) alpha and gamma bias how many states there are. We're telling it to expect
        # one state (conservative)
        # (2) kappa controls the self-transition bias. Bigger number means becomes more expensive
        # for states to non-self-transition (that is, change to a different state).
        model = WeakLimitStickyHDPHMM(
            alpha=1.,gamma=1.,init_state_distn='uniform',
            kappa=500.,
            obs_distns=[Gaussian(**obs_hypparams) for _ in range(10)],
        )

        return model

    ### model fit results

    @property
    def model_stateseq(self):
        return self.model.stateseqs[0]

    @property
    def model_durations(self):
        return self.model.durations[0]

    @property
    def model_stateseq_norep(self):
        return self.model.stateseqs_norep[0]

    @property
    def model_redgreenseq(self):
        # construct corresponding mean sequence
        m = self.model
        return np.array([m.obs_distns[state].mu for state in self.model_stateseq])

    @property
    def model_fretseq(self):
        g, r = self.model_redgreenseq.T
        return r/(r+g)

    ### plotting

    def model_fretplot(self):
        raw_g, raw_r = self.unsmoothedGrI, self.unsmoothedRedI
        raw_fretseq = raw_r/(raw_r+raw_g)

        t_model = np.arange(self.start, self.end if self.end is not None else self.unsmoothedGrI.shape[0])
        g, r = self.model_redgreenseq.T
        fretseq = r/(r+g)

        fig, (ax1, ax2) = plt.subplots(2,1,figsize=(12,6), sharex=True)

        plt.subplot(ax1)
        # plt.title('intensities')
        fig.suptitle(self.name)
        # plot raw values
        plt.plot(raw_g,'g-',alpha=0.5)
        plt.plot(raw_r,'r-',alpha=0.5)
        # plot model values
        plt.plot(t_model,g,'g--',linewidth=2)
        plt.plot(t_model,r,'r--',linewidth=2)
        plt.xlabel('Frame')
        plt.ylabel('Intensity (a.u.)')

        plt.subplot(ax2)
        # plt.title('fret values')
        # plot raw values
        plt.plot(raw_fretseq,'k-',alpha=0.5)
        # plot model values
        plt.plot(t_model,fretseq,'b-')
        # set ylim to handle outliers
        plt.ylim(scoreatpercentile(raw_fretseq,5)-0.2,scoreatpercentile(raw_fretseq,95)+0.2)
        plt.xlabel('Frame')
        plt.ylabel('FRET')

        plt.xlim(0,raw_g.shape[0])
        # fig.suptitle(self.name)

    ### saving to matfile

    def save_matfile(self,matfiledir=None):
        if matfiledir is None:
            mkdir('Results')
            matfilepath = join('Results',self.name+'Results')
        else:
                matfilepath = join(matfiledir,self.name+'_Results')
                
        tosave = ['start','unsmoothedGrI','unsmoothedRedI',
                'rawGrI','rawRedI','unsmoothedFRET','FRET','RedI','GrI',
                'model_data','fps','t_Inj']
        mdict = {key:self.__dict__[key] for key in tosave}
        mdict.update({
            'model_stateseq':self.model_stateseq,
            'model_stateseq_norep':self.model_stateseq_norep,
            'model_durations':self.model_durations,
            'model_redgreenseq':self.model_redgreenseq,
            'model_fretseq':self.model_fretseq,
            })
        savemat(matfilepath,mdict)


def load_experiment(expname,datadir=None,resultsdir=None):
    '''
    Given an experiment name, returns a dict where keys are trace names and
    values are corresponding Traces.

    Reads the experiment's 'goodtraces.txt' file, if it exists, and only
    includes in the returned dict traces whose lines do not start with '#'.

    # Example #

    traces = load_experiment('wtWB')
    model.add_data(traces['Spot1_186_141208'].model_data)
    '''

    default_datadir, default_resultsdir = getdirs()
    datadir = datadir if datadir is not None else default_datadir
    resultsdir = resultsdir if resultsdir is not None else default_resultsdir

    goodfile = join(resultsdir,expname,'goodtraces.txt')

    if isfile(goodfile):
        # parse goodtraces.txt, create corresponding Traces
        with open(goodfile,'r') as infile:
            lines = infile.readlines()

        traces = {}
        for l in lines:
            if not l.startswith('#'):
                # This line gets the filename, splitting off anything after a whitespace
                matfilepath = join(datadir,l.split()[0])
                name = _name_from_matfilepath(matfilepath)
                kwds = dict(pair.split('=') for pair in l.split()[1:])
                traces[name] = Trace(matfilepath, **kwds)

        return traces

    else:
        raise ValueError('Did not find goodtraces file: %s' % goodfile)


def load_trace(datadir,tracename):
    '''
    Given a trace name, returns the corresponding Trace object.

    # Example #

    trace = load_trace('Spot1_186_141208')
    model.add_data(trace.model_data)
    '''

    # Turn the _ in the original name back into an /
    tempname = tracename.split('_Spot',1)
    tracefilename = 'Spot' + tempname[1] + '.mat'
    datasubdir = join(datadir,tempname[0])
    for root, dirnames, filenames in os.walk(datasubdir):
        for filename in filenames:
            if filename == tracefilename:
                return Trace(join(root,filename))


def _name_from_matfilepath(f):
    assert not f.endswith('/')
    rest, filename = os.path.split(f)
    _, dirname = os.path.split(rest)
    spot_id, _ = splitext(filename)
    return dirname + '_' + spot_id


def mkdir(path):
    # from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    import errno
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

