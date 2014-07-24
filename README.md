#Single-Molecule FRET Analysis Code#

This is a software suite designed to calculate FRET-versus-time traces from a standard prism TIRF-based single-molecule FRET microscopy setup.  The most current version of this suite can be found at [GitHub.com](https://github.com/stephlj/smFRETcode).  It can be run as a stand-alone analysis suite, but it is written as modularly as possible in the hopes that it can be adapted fairly easily for other microscopy setups and acquisition software.  A limited number of setup-specific parameters can be changed in the smFRETsetup file, but it is highly likely that you will need to modify at least some parts of this code (for example, the file loading function, unless images from your camera are saved the same way ours are; see the accompanying manual for more details). The included manual discusses some additional options for customizing this suite to suit your particular needs.

The general outline of the analysis workflow is: (1) Calculate a map that correlates pixels in the acceptor channel image to pixels in the donor channel image, as these will never be perfectly aligned in an smFRET setup. Or, load an old one. (2) Find fluorescent spots in a movie or set of movies and allow the user to scroll through intensity-versus-time traces for each spot, and then save good traces for further analysis later.

More information (including detailed derivations of, for example, the channel mapping linear algebra) can be found in the documentation/manual  that accompanies this repository.

I would be very grateful if you would [open an issue](http://www.youtube.com/watch?v=TJlYiMp8FuY) if you find bugs in the code or find it doesn't work well with your data. Also, a port of parts or all of this software to Python (rather than Matlab), or some other freely available language, would be a great help to the scientific community.

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

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE.txt file that accompanies this software; it can also be found at the [GNU website](http://www.gnu.org/licenses/).
