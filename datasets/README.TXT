This folder contains the video sequences used in the paper. They are stored in a .mat format and the first
frame represents the background. All frames were processed to have a certain resolution and to contain only
grayscale images. 

-----------------------------------------------------------------------
Each file contains the following variables:

* imStack: 3-dimensional array, where imStack(:,:,f) represents frame f
* nFrames: number of frames
* siz:     resolution of each frame
-----------------------------------------------------------------------


-----------------------------------------------------------------------
Sources of the sequences:

***************
Changedetection
***************

Link: 
  http://www.changedetection.net/

Publication: 
  N. Goyette, P.-M. Jodoin, F. Porikli, J. Konrad, and P. Ishwar
  changedetection.net: A new change detection benchmark dataset
  Proc. IEEE Workshop on Change Detection (CDW-2012) at CVPR
  2012 

Datasets:
  Name                   |  Frames   |  Resolution 
  baselineHighway        |  1-500    |  120x160
  baselinePETS2006       |  1-1200   |  144x180
  thermalPark            |  20-600   |  144x176
  cameraJitterTraffic    |  1-250    |  120x160
  dynamicBackgroundCanoe |  820-1100 |  160x120


*********
Stuttgart
*********

Link:
  http://www.vis.uni-stuttgart.de/index.php?id=sabs

Publication: 
  S. Brutzer, B. Höferlin
  Evaluation of Background Subtraction Techniques for Video Surveillance
  Computer Vision and Pattern Recognition (CVPR) pp. 1937-1944
  2011 

Dataset:
  (from SABS-Basic.rar)
  Name                        |  Frames    |  Resolution 
  Stuttgart_144x192_0011_0310 |  0011-0310 |  144x192


**************
Project Caviar
**************

Link: 
  http://groups.inf.ed.ac.uk/vision/CAVIAR/CAVIARDATA1/

Publication: 
  EC Funded CAVIAR project/IST 2001 37540

Datasets:
  (from OneShopOneWait1front.tar.gz)
  Name                          |  Frames    |  Resolution 
  VascoDaGama_144x192_0344_1235 |  0344-1235 |  288x384 (original)


*******************
Video Trace Library
*******************

Link:
  http://trace.eas.asu.edu/index.html

Publication:
  P. Seeling, M. Reisslein
  Video Traffic Characteristics of Modern Encoding Standards: 
  H.264/AVC with SVC and MVC Extensions and H.265/HEVC
  The Scientfic World Journal, Vol. 2014, Article ID 189481, pp. 1-16
  2014

Dataset:
  Name        |  Frames    |  Resolution 
  hall128x128 |  18-300    |  128x128


**************
GarrettWarnell
**************

Link:
  http://garrettwarnell.com/ARCS-1.0.zip
  (for the original, see view 5 of 
  ftp://ftp.cs.rdg.ac.uk/pub/PETS2009/Crowd_PETS09_dataset/a_data/Crowd_PETS09/S2_L1.tar.bz2

Publication:
  G. Warnell, S. Bhattacharya, R. Chellappa, T. Basar
  Adaptive-Rate Compressive Sensing via Side Information 
  http://arxiv.org/abs/1401.0583
  2014

Dataset:
  Name        |  Frames    |  Resolution 
  PETS09_S2L1 |  1-171     |  116x116
-----------------------------------------------------------------------
