## Least-squares horizons with local slopes and multi-grid correlations

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
to implement the 2D and 3D horizon extraction methods,
discussed in our Geophysics paper 
[Least-squares horizons with local slopes and multi-grid correlations]
(http://www.jsg.utexas.edu/wu/files/wu2018LeastSquaresHorizons.pdf).

If you find this work helpful in your research, please cite:

    @article{doi:wu2018least,
        author = {Xinming Wu and Sergey Fomel},
        title = {Least-squares horizons with local slopes and multigrid correlations},
        journal = {GEOPHYSICS},
        volume = {83},
        number = {4},
        pages = {IM29-IM40},
        year = {2018},
        doi = {10.1190/geo2017-0830.1},
        URL = {https://doi.org/10.1190/geo2017-0830.1},
    }

This software depends on that in the [Mines Java Toolkit
(JTK)](https://github.com/dhale/jtk/). If you want to do more than browse the
source code, you must first download and build the Mines JTK using
[Gradle](http://www.gradle.org). The build process for software in
this repository is the same.

Like the Mines JTK, this is a toolkit for computer programmers. It is not a
complete system for seismic interpretation of geologic faults. Others
(including commercial software companies) have built such systems using
earlier versions of one or more of the software tools provided in this
repository.

Although based on [research first published in 2012]
(http://inside.mines.edu/~dhale/research.html)
most of the software in this repository was written in the summer of 2014
as Dave prepared for the upcoming lecture tour. This software will change
during the course of that tour in response to questions and discussions
inspired by the lecture.

### Summary

Here are brief descriptions of key components:

#### HorizonExtractor2
Implements 3 methods for extracting 2D horizon curves:
1) predictive horizons with local slopes only;
2) least-squares horizons with local slopes only;
3) least-squares horizons with both local slopes and multi-grid correlations.

#### HorizonExtractor3
Implements 2 methods for extracting 3D horizon surfaces:
1) least-squares horizons with local slopes only;
2) least-squares horizons with both local slopes and multi-grid correlations.

#### demo2
type ./jy demo2 
to run the three 2D horizon extraction methods 
in 3 different examples 

#### demo3
type ./jy demo3 
to run the two 3D horizon extraction methods 
in 2 examples 

---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
