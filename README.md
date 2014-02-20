#Single-Molecule FRET Analysis Code#

Welcome to Stephanie Johnson's smFRET code repository!  While you are welcome to download the repository as-is and run it off the shelf (at your own risk, I make no guarantees this code is "right"), this code is written as modularly as possible, so that it can be adapted as needed for your own experiment. For example, `FindSpots` should find spots in any image you give it (as a function of user-specified parameters, such as the sizes of the spots, how close together they can be, etc). You should be able to write your own FindSpots function and as long as it takes an image and returns spot centers, it will interface with the rest of the code. A wrapper function that calls various portions of the analysis in proper order, called `smFRET.m`, is included but is highly tailored to my setup.

The general outline of the analysis workflow is:
(1) Calculate a map that correlates pixels in the acceptor channel image to pixels in the donor channel image, as these will never be perfectly aligned in an smFRET setup.  Or, load an old one.
(2) Find fluorescent spots in a movie or set of movies and allow the user to scroll through intensity-versus-time traces for each spot, and then save good traces for further analysis later.

More information (including detailed derivations of, for example, the channel mapping linear algebra) can be found in the (currently non-existent) smFRET Users's Manual pdf that accompanies this repository.

I would be very grateful if you would [open an issue](http://www.youtube.com/watch?v=TJlYiMp8FuY) if you find bugs in the code or find it doesn't work well with your data.

##System Requirements##

Requires Matlab with the following toolboxes installed:
* Image Processing Toolbox
* Curve Fitting Toolbox
* Optimization Toolbox

This code was written for analyzing prism-based TIRF images; it may not work well for objective-based setups (for example, you may have to alter or re-write `FindSpots`).

The wrapper function that is tailored to our experimental setup assumes the following about the data:
* We acquire data with the open-source UCSF software [Micro-Manager](http://www.micro-manager.org), using its Multi-D acquisition tool. Micro-Manager then stores the data as follows: You give Micro-Manager a root file name, and then every time you hit acquire it creates a new directory called rootname_1, rootname_2, etc.  Each of these directories contains one tiff per frame from the acquired movie, with filename img_000000000_000.tif, img_000000001_000.tif, etc. In each movie's directory Micro-Manager also saves a metadata.txt file with information about frame rate, etc.
* The personalized function `LoadUManagerTifs` assumes this file and directory structure; if you use something other than Micro-Manager to acquire data, you will have to write your own function to load movie files into Matlab for analysis.
* Channel mapping is accomplished by means of images of fluorescent beads, which are very bright in the donor channel and somewhat bright in the acceptor channel.  The wrapper function assumes these beads are in folders that contain the word "Bead".  Although the wrapper function does not assume there are movies in these folders, I've found it's better to average 10-20 frames before finding bead positions.
* While you can set certain setup-specific parameters in `smFRETsetup.m`, such as whether the donor channel is the right half or left half (or top or bottom) of the image, you will have to write your own wrapper function for certain extensions such as images obtained on a three-color FRET setup.  

##Getting Started##

Open the `smFRETsetup.m` file and change parameters as necessary. For example, you will probably want to change `defaultdatadir` and `codedir`.  You will also want to follow the instructions under Display Defaults to make sure figures will show up well on your screen.

Sample data that can be analyzed by my code can be found in my [smFRETdata](https://github.com/stephlj/smFRETdata) repository. To analyze the sample data in that  repository, run

```matlab
smFRET('Nucleosomes4000xdil')
```

from the command line.  Pass '1' as the second optional input to run in debug mode, which will provide some additional figures.

You will first be asked if you want to create a channel mapping, or load an old one.  You can perform a channel mapping with any of the three provided bead maps; however, the code will fail if you use BadBeads (working on that).

If you choose to perform a channel mapping, press enter and navigate to one of the bead folders in SampleData.  Follow the instructions to adjust bead-detection thresholds and find beads.

After you're finished, or if you decided to use the existing map in the GoodBeads directory, you will then be asked to navigate to the directory with data to be analyzed (which in this case is RealData).

By the way, any time it asks you to press "anything else", it will probably crash if you press the space bar.  Choose a different anything else!

##Acknowledgements##

I am extremely grateful to [Matt Johnson](http://www.mit.edu/~mattjj/) and [Luke Breuer](http://luke.breuer.com) for extensive help developing this software.  This project would have been dead in the water without the significant intellectual and time contributions that they made! I am also indebted to TJ Ha's lab at UIUC, which has made what I consider to be "industry-standard" [smFRET code](http://cplc.illinois.edu/software/) freely available on their website. In this repository I make frequent allusions to the Ha lab code, noting where my code differs from theirs and where I have followed their lead. Finally I am grateful for financial support through NIH grants to my PI, Dr. Geeta Narlikar.

##Licensing Information##

Coming soon!

Copyright (c) 2014 Stephanie Johnson, Ph.D.

University of California, San Francisco

##To Do##
* Image smoothing before spot finding? Seems to be fine without it ... 
* Refining intensity calculation--Gaussian fit or pseudo-fit as in Ha lab code
* Whole image background subtraction as in Ha lab code 
* Have the option to re-find spots every 10 frames

