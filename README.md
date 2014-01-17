#Single-Molecule FRET Analysis Code#

Requires Matlab with the following toolboxes installed:
* Image Processing Toolbox
* Curve Fitting Toolbox
* Optimization Toolbox

##Getting Started##

Open the ```smFRETsetup.m``` file and change parameters as necessary. For example, you will probably want to change ```defaultdatadir``` and ```codedir```.  You will also want to follow the instructions in ```smFRETsetup``` under Display Defaults to make sure figures will show up well on your screen.

To analyze the sample data in this repository, run

```matlab
smFRET('Nucleosomes4000xdil')
```

from the command line.  Pass '1' as the second optional input to run in debug mode, which will provide some additional figures.

You will first be asked if you want to create a channel mapping, or load an old one.  You can perform a channel mapping with any of the three provided bead maps; however, the code will fail if you use BadBeads (working on that).

If you choose to perform a channel mapping, press enter and navigate to one of the bead folders in SampleData.  Follow the instructions to adjust bead-detection thresholds and find beads.

After you're finished, or if you decided to use the existing map in the GoodBeads directory, you will then be asked to navigate to the directory with data to be analyzed (which in this case is RealData).

By the way, any time it asks you to press "anything else", it will probably crash if you press the space bar.  Choose a different anything else!

##To Do##
* Channel Mapping: my CalcChannelMap versus Matlab's fitgeotrans.  Fitgeotrans does better when you overlay the images, but mine does better in terms of error in calculated spot center versus detected spot center, especially near the edge between channels.
* Image smoothing before spot finding?
* Figure out image scaling, loading as uint16 vs double
* Refining spot centers--Gaussian fit or pseudo-fit as in Ha lab code
* Calculating intensities and subtracting background better as in Ha lab code 


