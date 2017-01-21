from __future__ import division

import sys
import cPickle as pickle

import matplotlib.pyplot as plt
from os.path import join

from fret import load_experiment
import config

# Haven't decided whether I want this reproducibility here yet or not
# np.random.seed(0)

# Load all the data

# Addition here to work with the new loading scheme:

datadir, resultsdir = config.getdirs()
modelfitsdir = join(resultsdir,sys.argv[1])
figdir = join(resultsdir,sys.argv[1],'ResultsFigs')

alltraces = load_experiment(sys.argv[1])

# Update Steph 5/2015: Allowing an option to re-analyze only one trace in this set,
# by passing a second argument
if len(sys.argv) not in (2,3):
    usage()
    sys.exit(1)
elif len(sys.argv) == 2:
    # For each trace in alltraces, fit (separately) a model with Gibbs
    for name, trace in alltraces.items():
        trace.run_gibbs(500)
        print 'Done with %s!' % name
else:
    name = sys.argv[2]
    trace = alltraces[name]
    trace.run_gibbs(500)
    print 'Done with %s!' % name

# save the model-fit traces for later, just in case
with open(join(modelfitsdir, 'modelfits.pkl'), 'w') as outfile:
    pickle.dump(alltraces, outfile, protocol=-1)

# Now save everything as a mat file for Steph's code
if len(sys.argv) == 2:
    for name, trace in alltraces.items():
        trace.model_fretplot()
        plt.savefig(join(figdir,'model_fretplot_' + name + '.png'))
        plt.close()
        trace.save_matfile(modelfitsdir)
else:
    trace.model_fretplot()
    plt.savefig(join(figdir,'model_fretplot_' + name + '.png'))
    plt.close()
    trace.save_matfile(modelfitsdir)
