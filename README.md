#Traces (single molecule FRET analysis code)#

Traces is a software suite designed to calculate FRET-versus-time traces from a standard prism-based TIRF single molecule FRET (smFRET) microscopy setup.  The most current version of this suite can be found at [GitHub.com](https://github.com/stephlj/smFRETcode).  It can be run as a stand-alone analysis suite, but it is written as modularly as possible in the hopes that it can be adapted fairly easily for other microscopy setups and image acquisition software.  Some setup-specific parameters can be changed in the smFRETsetup file, but you may need to modify at least some parts of this code (for example, the file loading function, unless images from your camera are saved the same way ours are; see the accompanying manual for more details). 

The general outline of the analysis workflow is: (1) Calculate a map that correlates pixels in the acceptor channel image to pixels in the donor channel image, as these will never be perfectly aligned in an smFRET setup. Or, load an old one. (2) Find fluorescent spots in a movie or set of movies and allow the user to scroll through intensity-versus-time traces for each spot, and then save good traces for further analysis later.

More information (including detailed derivations of, for example, the channel mapping linear algebra) can be found in the documentation/manual  that accompanies this repository.

I would be  grateful if you would [open an issue](http://www.youtube.com/watch?v=TJlYiMp8FuY) if you find bugs in the code or find it doesn't work well with your data. Also, a port of parts or all of this software to Python (rather than Matlab), or some other freely available language, would be a great help to the scientific community.

##Why use Traces?##

Traces developed out of three main problems we were facing analyzing smFRET data in our lab: (1) The need for more customization options for our particular microscope than were offered by other software options; (2) ; and (3) the need to be able to "trace" dye intensity-versus-time data back to the images that generated those data. This last was so important to us, especially for distinguishing complicated but real signals from experimental artifacts, that we named our software "Traces"--we want the user to always be able to "trace back" to the original images, to the analysis parameters used to get intensity-versus-time data out of those images, etc.  (Plus, the point of this software is to generate intensity-versus-time trajectories, which are sometimes called "traces".)

##System Requirements##

Requires Matlab with the following toolbox installed:
* Image Processing Toolbox

and will use the following if available:
* Curve Fitting Toolbox
* Optimization Toolbox

This code was written for analyzing prism-based TIRF images; it may not work well for objective-based setups (for example, you may have to alter or re-write `FindSpots`).

##Data Format Requirements##

This software assumes that images from the camera were acquired using the open-source UCSF software [Micro-Manager](http://www.micro-manager.org), using its Multi-D acquisition tool. You will have to modify or re-write `LoadUManagerTifs` and `GetInfoFromMetaData` if you use different acquisition software.

More information can be found in the accompanying manual.

##Sample Data##

Sample data that can be analyzed by my code can be found in my [smFRETdata](https://github.com/stephlj/smFRETdata) repository. To analyze the sample data in that  repository, run

```matlab
smFRET('HighFRET25pM')
```

from the command line and then follow the instructions in the manual.

##Acknowledgements##

This project was a collaborative effort with [Matt Johnson](http://www.themattjohnson.com/) and [Luke Breuer](http://luke.breuer.com), who made significant intellectual and time contributions to its development.  I am also indebted to TJ Ha's lab at UIUC, which has made what I consider to be "industry-standard" [smFRET code](http://cplc.illinois.edu/software/) freely available on their website. In this repository I make frequent allusions to the Ha lab code, noting where my code differs from theirs and where I have followed their lead. Thanks also to the UCSF SMUG crew (Margaret Elvekrog, Thomas Noriega, Megan Mayerle, Sarah Ledoux, Daniel Elnatan) and the Narlikar lab, especially John Leonard, for helpful discussions about code as well as experimental development, and to Sarah Ledoux for the DNAs used for the sample data. Finally I am grateful for financial support through NIH grants to my PI, Dr. Geeta Narlikar, and a Leukemia and Lymphoma Society Career Development Program fellowship.

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