## Least-squares horizons with local slopes and multi-grid correlations

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
for 2D and 3D horizon extraction that is discussed in our Geophysics paper 
[Least-squares horizons with local slopes and multi-grid correlations](http://www.jsg.utexas.edu/wu/files/wu2018LeastSquaresHorizons.pdf).

If you find this work helpful in your research, please cite:

    @article{doi:wu2018least,
        author = {Xinming Wu and Sergey Fomel},
        title = {Least-squares horizons with local slopes and multi-grid correlations},
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
complete system for seismic interpretation of geologic horizons. Others
(including commercial software companies) have built such systems using
earlier versions of one or more of the software tools provided in this
repository.

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

please email me xinming.wu@beg.utexas.edu to ask for the 3D datasets

## Examples

2D and 3D examples published in the [paper]((http://www.jsg.utexas.edu/wu/files/wu2018LeastSquaresHorizons.pdf))
### 2D examples

#### 1) Teapot Dome (provided by the Rocky Mountain Oilfield Testing Center)
<p align="left">
  <img src="png/2d/tpd/tpdp.png" width="295px" height="500px"/>
  <img src="png/2d/tpd/tpds.png" width="295px" height="500px"/>
  <img src="png/2d/tpd/tpdm.png" width="295px" height="500px"/>
</p>

#### 2) Netherlands off-shore F3 (provided by the Dutch Government through TNO and dGB Earth Sciences)
<p align="left">
  <img src="png/2d/f3d/f3dp.png" width="295px" height="200px"/>
  <img src="png/2d/f3d/f3ds.png" width="295px" height="200px"/>
  <img src="png/2d/f3d/f3dm.png" width="295px" height="200px"/>
</p>

#### 3) Curt (provided by Australian government)
<p align="left">
  <img src="png/2d/curt/curtp.png" width="695px" height="250px"/>
</p>
<p align="left">
  <img src="png/2d/curt/curts.png" width="695px" height="250px"/>
</p>
<p align="left">
  <img src="png/2d/curt/curtm.png" width="695px" height="250px"/>
</p>

---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
