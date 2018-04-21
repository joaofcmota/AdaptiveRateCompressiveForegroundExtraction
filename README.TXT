This software provides Matlab code to execute the experiments described in

[1] J. F. C. Mota, N. Deligiannis, A. C. Sankaranarayanan, V. Cevher, 
    and M. R. D. Rodrigues
    *Adaptive-Rate Reconstruction of Time-Varying Signals with Application 
    in Compressive Foreground Extraction*
    accepted, IEEE Transactions on Signal Processing
    http://arxiv.org/abs/1503.03231
    2015

and

[2] J. F. C. Mota, N. Deligiannis, A. C. Sankaranarayanan, V. Cevher, 
    and M. R. D. Rodrigues
    *Dynamic Sparse State Estimation Using L1-L1 Minimization: Adaptive-Rate 
    Measurement Bounds, Algorithms and Applications*
    IEEE International Conference on Acoustics, Speech, and Signal Processing 
    (ICASSP)
    2015

If you want to cite this code, please use [1] and/or [2]. 


==================================================================================
Organization:


  **************
  createFigures:   Code for replicating the figures in [1]. The figures in the 
  **************   were postprocessed in LaTeX, so the visuals are not exactly as 
                   in the paper. Some experiments do not generate figures, only 
                   the data.

      * FFTMeasurements: Generates data to create Figure 4 in [1].

      * GaussianMeasurements: Generates data to create Figures 3(a), 3(b), 3(e),
                              and 3(f) in [1].

      * GaussianMeasurementsNoise: Generates data to create Figures 3(c) and 3(d).



  ***********
  Algorithms:   Code to solve L1-L1 minimization. There are two implementations:
  ***********   

      * ADMM:   (folder basisPursuitPlusL1) this is the solver used in the 
                experiments in

                [3] J. F. C. Mota, N. Deligiannis, M. R. D. Rodrigues
                    Compressed Sensing with Prior Information: Optimal Strategies, 
                    Geometry, and Bounds
                    in revision, IEEE Transactions on Information Theory
                    arXiv: http://arxiv.org/abs/1408.5250

      * DECOPT: (folder PrimalDualFramework-TranDinhCevher) this is the solver 
                proposed in 

                [4] Q. Tran-Dinh, V. Cevher
                    Constrained convex minimization via model-based excessive gap
                    Proceedings of the annual conference on Neural Information 
                    Processing Systems Foundation (NIPS) 
                    2014

                [5] Q. Tran-Dinh, V. Cevher 
                    A Primal-Dual Algorithmic Framework for Constrained Convex 
                    Minimization
                    LIONS Tech. Report EPFL-REPORT-199844 
                    2014

                The code for the original solver can be found in:

                http://lions.epfl.ch/decopt/

                We slightly modified the code to make it applicable to L1-L1 
                minimization.


      --------------------------------------------------------------------------

      * Modified-CS: (folder otherApproaches/modifiedCS-Vaswani) this solves the
                     Modified-CS problem in

                     [6] N. Vaswani, W. Lu
                         Modified-CS: Modifying Compressive Sensing for Problems
                         With Partially Known Support 
                         IEEE Transactions on Signal Processing, Vol. 58, No. 9, 
                         2010.


  *********
  datasets:   Contains .mat files with the video sequences used in [1]. See the
  *********   README.TXT file in that folder for the sources of the sequences.

==================================================================================



==================================================================================
Copyright (C) 2015  Joao Mota

This program is free software: you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program.  If not, see <http://www.gnu.org/licenses/>.
==================================================================================

