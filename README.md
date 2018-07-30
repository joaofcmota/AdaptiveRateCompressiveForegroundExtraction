# Adaptive Rate Reconstruction of Time Varying Signals With Applications in Compressive Foreground Extraction 

This code replicates the experiments in 

1. **[Adaptive-Rate Reconstruction of Time-Varying Signals with Application in
   Compressive Foreground
   Extraction](http://dx.doi.org/10.1109/TSP.2016.2544744)**.  
  J. F. C. Mota, N. Deligiannis, A. C. Sankaranarayanan, V. Cevher, M. R. D. Rodrigues.  
  IEEE Transactions on Signal Processing, Vol. 64, No. 14, pp. 3651-3666, 2016.  
  [link](http://dx.doi.org/10.1109/TSP.2016.2544744),
  [arXiv](http://arxiv.org/abs/1503.03231)

and

2. **[Dynamic Sparse State Estimation Using L1-L1 Minimization: Adaptive-Rate Measurement Bounds, Algorithms and Applications]( http://dx.doi.org/10.1109/ICASSP.2015.7178588 )**.  
  J. F. C. Mota, N. Deligiannis, A. C. Sankaranarayanan, V. Cevher, M. R. D. Rodrigues. 
  IEEE International Conference on Acoustics, Speech, and Signal Processing
  (ICASSP), Brisbane, 2015. 
  [link]( http://dx.doi.org/10.1109/ICASSP.2015.7178588 )


## Organization


* **createFigures:**
  Replicates the figures in [1]. As all figures were postprocessed
  in LaTeX, the visuals are not exactly as in the paper. Some experiments only
  generate data, not figures.

    * FFTMeasurements: Generates data to create Figure 4 in [1].

    * GaussianMeasurements: Generates data to create Figures 3(a), 3(b), 3(e), and 3(f) in [1].

    * GaussianMeasurementsNoise: Generates data to create Figures 3(c) and 3(d).


* **Algorithms:**
  Code to solve the L1-L1 minimization problem:

  <a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\|x&space;-&space;\overline{x}\|_1&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\|x&space;-&space;\overline{x}\|_1&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b&space;\end{array}" title="\begin{array}[t]{ll} \underset{x}{\text{minimize}} & \|x\|_1 + \|x - \overline{x}\|_1 \\ \text{subject to} & Ax = b. \end{array}" /></a>
  
  There are two implementations:

  * basisPursuitPlusL1/ADMM:   
    Solver used in the experiments in

    **[Compressed Sensing with Prior Information: Strategies, Geometry, and
    Bounds](https://doi.org/10.1109/TIT.2017.2695614)**.  
    J. F. C. Mota, N. Deligiannis, M. R. D. Rodrigues.  
    IEEE Transactions on Information Theory, Vol. 63, No. 7, pp. 4472-4496, 2017.  
    [link](https://doi.org/10.1109/TIT.2017.2695614), 
    [arXiv](http://arxiv.org/abs/1408.5250), 
    [Github page](https://github.com/joaofcmota/cs-with-prior-information)

  * PrimalDualFramework-TranDinhCevher/DECOPT:  
    Solver proposed in 

    **[Constrained convex minimization via model-based excessive gap](
    http://papers.nips.cc/paper/5494-constrained-convex-minimization-via-model-based-excessive-gap)**.  
    Q. Tran-Dinh, V. Cevher.  
    Proceedings of the annual conference on Neural Information Processing
    Systems Foundation (NIPS), 2014.  
    [link](http://papers.nips.cc/paper/5494-constrained-convex-minimization-via-model-based-excessive-gap),
    see also [code webpage](http://lions.epfl.ch/decopt/)

    The original code was slightly modified to be applicable to L1-L1 minimization.


  * otherApproaches/modifiedCS-Vaswani/Modified-CS:  
    Solver for the *Modified-CS* problem:

    <a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x_{T^c}\|_1\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b\,,&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x_{T^c}\|_1\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b\,,&space;\end{array}" title="\begin{array}[t]{ll} \underset{x}{\text{minimize}} & \|x_{T^c}\|_1\\ \text{subject to} & Ax = b\,, \end{array}" /></a>

    which was proposed in

    **[Modified-CS: Modifying Compressive Sensing for Problems With Partially Known Support](
    https://ieeexplore.ieee.org/abstract/document/5471173/)**.  
    N. Vaswani, W. Lu.  
    IEEE Transactions on Signal Processing, Vol. 58, No. 9, 2010.  
    [link](https://ieeexplore.ieee.org/abstract/document/5471173/)

  * datasets: 
    Contains `.mat` files with the video sequences used in [1]. 
    For further information about the datasets check the `README.TXT` file in
    that folder.


---

License: [ GPLv3 ]( https://www.gnu.org/licenses/gpl-3.0.en.html )
