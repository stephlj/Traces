# Traces: single molecule FRET analysis code

Traces calculates FRET-versus-time trajectories (“traces”) from a standard prism-based TIRF single-molecule FRET (smFRET) microscopy setup. It can be run as a stand-alone analysis suite, but it is written as modularly as possible with the hope that it can be adapted fairly easily to other microscopy setups and acquisition software.  

We have attempted to avoid as much as possible "buried" parameters, collecting them instead in the TracesSetup.m file.  TracesSetup.m is the only function a user should have to edit regularly (though some initial modifications may need to be made to, for example, the image loading function; see the accompanying manual for more details).

The general outline of the analysis workflow is: (1) Calculate a map that correlates pixels in the acceptor channel image to pixels in the donor channel image, as these will never be perfectly aligned in an smFRET setup. Or, load an old one. (2) Find fluorescent spots in a movie or set of movies, calculate intensity-versus-time traces for each spot, and then allow the user to save good traces for downstream analyses.

Traces also includes code for downstream analyses of FRET-versus-time trajectories, based on fitting hidden Markov models (HMMs) with Gaussian emissions to the data to quantify "dwells" or "pauses". We have also developed a companion package called [Slopey](https://github.com/stephlj/slopey) for quantifying fast, non-instantaneous transitions between dwells.

More information (including detailed derivations of, for example, the channel mapping linear algebra) can be found in the manual that accompanies this repository.

Please [open an issue](http://www.youtube.com/watch?v=TJlYiMp8FuY) if you find bugs in Traces. 

## Why Traces?

Traces developed out of three main needs we were facing in analyzing smFRET data in our lab that were not met by other available software options: the need (1) for more customization options for our particular setup; (2) for more automation, additional computational and manual checks to enrich for good data with less work on the part of the user, and the ability to save and rerun select parts of the analysis, in order to handle the large data sets we generate; and (3) to be able to "trace" intensity-versus-time data back to the images that generated those data. 

This last point was so important to us, especially for distinguishing complicated but real signals from experimental artifacts, that we named our software "Traces": we want the user to always be able to "trace back" to the original images, to the analysis parameters used to get intensity-versus-time data out of those images, etc. For example, Traces allows the user to play sections of the movie, or to examine averages of sections of the movie, that generated the intensity-versus-time data being displayed for analysis.

More generally, just as controls and sanity checks are essential for rigorously interpreting experimental results, similar controls and checks are crucial for verifying that any software's output is "right". Therefore Traces offers the user visualizations, plots, and sanity checks after every manipulation of the data, so that the user can assess how Traces is doing, instead of treating the process of taking a set of raw camera images and generating FRET-versus-time traces as a black box.

Several other smFRET analysis software packages exist that span a spectrum from bare-bones to highly automated, with Traces falling somewhere in the middle. At one end of the spectrum is the pioneering [package](http://cplc.illinois.edu/software/) developed by TJ Ha’s lab at UIUC. Algorithmically, Traces is very similar to the Ha lab code; however, Traces does not treat the analysis pipeline from images to intensity-versus-time trajectories as a black box. Traces also offers a rudimentary GUI and some additional features such as the ability to re-localize FRET pairs, if the calculated channel map was sub-optimal. At the other end of the spectrum are packages such as [iSMS](http://isms.au.dk/) and [SPARTAN](http://www.scottcblanchardlab.com/#!software/i7met), which include GUIs and many automated features, and can be run without as much programming knowledge on the part of the user. However, these packages tend to be more difficult to customize and may work best for setups that are similar to the setup for which they were originally designed. For example, iSMS finds FRET pairs using a nearest-neighbor algorithm, rather than a channel mapping procedure, which does not work well with the high spot density and single laser illumination in our setup. 

If none of the available software options work well for you, it is our hope that the algorithmic details in the Traces manual and its modularity will help you to develop your own software. For example, you can replace the background correction function `CalcBkgnd.m` with one that uses a different method to compute a background to subtract (or one the subtracts no background at all). As long as your function has a compatible input/output structure as `CalcBkgnd.m,` you can still take advantage of the rest of the Traces suite.

## System Requirements

Requires Matlab with the following toolbox installed:
* Image Processing Toolbox

and will use the following if available:
* Curve Fitting Toolbox
* Optimization Toolbox

## Data Format Requirements

Traces can load .tif and .pma image files. The functions `LoadUManagerTifs` and `GetInfoFromMetaData` assume a file storage structure like that generated by the open-source UCSF software [Micro-Manager](http://www.micro-manager.org), using its Multi-D acquisition tool. You will have to modify or re-write either the loading functions or your data structure if you use different acquisition software.

More information can be found in the accompanying manual.

## Sample Data

Sample data that can be analyzed by Traces can be found at [TracesSampleData](https://github.com/stephlj/TracesSampleData). To analyze the sample data in that repository, first copy the sample TracesSetup.m file from the SampleFiles subdirectory to the directory that contains Traces.m. Edit the parameters in TracesSetup.m (see comments in that file as well as additional details in the manual). Then run

```matlab
Traces('SNF2h51nMATP1mM')
```

from the command line and follow the instructions in the manual. The TracesSampleData repository also has analysis files to accompany the raw data.

## Citing Traces

If you use Traces to analyze data for a publication, please cite it!

Johnson S.L., Johnson M.J., Breuer L.A., and Narlikar G.J. (2018) Traces: single molecule FRET analysis code. DOI: 10.5281/zenodo.1211233.

## Acknowledgements

This project was a collaborative effort with [Matt Johnson](http://www.themattjohnson.com/) and [Luke Breuer](http://luke.breuer.com), who made significant intellectual and time contributions to its development.  I am also indebted to the freely available [smFRET analysis software](http://cplc.illinois.edu/software/) from TJ Ha's lab at UIUC. In this repository I make frequent allusions to the Ha lab code, noting where my code differs from theirs and where I have followed their lead. Thanks also to the UCSF SMUG crew (Margaret Elvekrog, Thomas Noriega, Megan Mayerle, Sarah Ledoux, Daniel Elnatan) and the Narlikar lab, especially John Leonard, for helpful discussions about code as well as experimental development. Finally I am grateful for financial support through NIH grants to my PI, Dr. Geeta Narlikar, and a Leukemia and Lymphoma Society Career Development Program fellowship.

## Licensing and Copyright Information

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