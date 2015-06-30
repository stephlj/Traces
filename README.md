#Traces (single molecule FRET analysis code)#

Traces is a software suite designed to calculate FRET-versus-time traces from a standard prism-based TIRF single-molecule FRET (smFRET) microscopy setup.  The most recent version of this suite can be found at [GitHub.com](https://github.com/stephlj/smFRETcode).  It can be run as a stand-alone analysis suite, but it is written as modularly as possible in the hopes that it can be adapted fairly easily to other microscopy setups and acquisition software.  

We have attempted to avoid as much as possible "buried" parameters, collecting them instead in the TracesSetup.m file.  TracesSetup.m is the only function a user should have to edit regularly (though some initial modifications may need to be made to, for example, the image loading function; see the accompanying manual for more details).

The general outline of the analysis workflow is: (1) Calculate a map that correlates pixels in the acceptor channel image to pixels in the donor channel image, as these will never be perfectly aligned in an smFRET setup. Or, load an old one. (2) Find fluorescent spots in a movie or set of movies, calculate intensity-versus-time traces for each spot, and then allow the user to save good traces for downstream analyses.

More information (including detailed derivations of, for example, the channel mapping linear algebra) can be found in the manual that accompanies this repository.

Please [open an issue](http://www.youtube.com/watch?v=TJlYiMp8FuY) if you find bugs in Traces. 

##Why use Traces?##

Traces developed out of three main problems we were facing analyzing smFRET data in our lab: (1) The need for more customization options for our particular setup than were offered by other software options; (2) the need for more automation, additional computational and manual checks to enrich for good data, and the ability to save and rerun select parts of the analysis, in order to handle the large data sets we generate; and (3) the need to be able to "trace" dye intensity-versus-time data back to the images that generated those data. This last was so important to us, especially for distinguishing complicated but real signals from experimental artifacts, that we named our software "Traces"--we want the user to always be able to "trace back" to the original images, to the analysis parameters used to get intensity-versus-time data out of those images, etc.  (Plus, the point of this software is to generate intensity-versus-time trajectories, which are sometimes called "traces".) 

More generally, just as controls and sanity checks are essential for rigorously interpreting experimental results, similar controls and checks are crucial for verifying that any code's output is "right". Therefore Traces offers the user visualizations, plots, and sanity checks after every manipulation of the original data set, so that the user can assess how Traces is doing, instead of treating the process of taking a set of raw camera images and generating FRET-versus-time traces as a black box.

##System Requirements##

Requires Matlab with the following toolbox installed:
* Image Processing Toolbox

and will use the following if available:
* Curve Fitting Toolbox
* Optimization Toolbox

##Data Format Requirements##

This software assumes that images from the camera were acquired using the open-source UCSF software [Micro-Manager](http://www.micro-manager.org), using its Multi-D acquisition tool. You will have to modify or re-write `LoadUManagerTifs` and `GetInfoFromMetaData` if you use different acquisition software.

More information can be found in the accompanying manual.

##Sample Data##

Sample data that can be analyzed by Traces can be found at [smFRETdata](https://github.com/stephlj/smFRETdata). To analyze the sample data in that  repository, run

```matlab
Traces('HighFRET25pM')
```

from the command line and then follow the instructions in the manual.

##Acknowledgements##

This project was a collaborative effort with [Matt Johnson](http://www.themattjohnson.com/) and [Luke Breuer](http://luke.breuer.com), who made significant intellectual and time contributions to its development.  I am also indebted to TJ Ha's lab at UIUC, which has made what I consider to be "industry-standard" [smFRET code](http://cplc.illinois.edu/software/) freely available on their website. In this repository I make frequent allusions to the Ha lab code, noting where my code differs from theirs and where I have followed their lead. Thanks also to the UCSF SMUG crew (Margaret Elvekrog, Thomas Noriega, Megan Mayerle, Sarah Ledoux, Daniel Elnatan) and the Narlikar lab, especially John Leonard, for helpful discussions about code as well as experimental development, and to Megan Mayerle for the DNAs used for the sample data. Finally I am grateful for financial support through NIH grants to my PI, Dr. Geeta Narlikar, and a Leukemia and Lymphoma Society Career Development Program fellowship.

##Licensing and Copyright Information##

Software and manual copyright (C)  2014 Stephanie Johnson, University of California, San Francisco.

 The MIT License (MIT)
 
 Copyright (c) 2014 Stephanie Johnson, University of California, San Francisco
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.