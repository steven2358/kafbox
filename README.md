Kernel Adaptive Filtering Toolbox
=================================

A Matlab benchmarking toolbox for kernel adaptive filtering.

The Kernel Adaptive Filtering Toolbox contains implementations of kernel adaptive filtering algorithms and tools to compare their performance.

Maintainer: [Steven Van Vaerenbergh](http://gtas.unican.es/people/steven) (steven at gtas dot dicom dot unican dot es)  
Contributors: [Miguel Lazaro-Gredilla](http://www.tsc.uc3m.es/~miguel), [Sohan Seth](http://www.sohanseth.com/)  
Official web: https://sourceforge.net/projects/kafbox  

This toolbox is a collaborative effort: every developer wishing to contribute code or suggestions can do so through https://github.com/steven2358/kafbox

Directories included in the toolbox
-----------------------------------

`data/` - data sets

`demo/` - demos and test files

`lib/` - algorithm libraries and utilities

Setup
-----

Run install.m

Octave / Matlab pre-2008a
-------------------------
This toolbox uses the `classdef` command which is not supported in Matlab pre-2008a and not yet in Octave. The older 0.x versions of this toolbox do not use `classdef` and can therefore be used with all versions of Matlab and Octave. http://sourceforge.net/projects/kafbox/files/

Usage
-----
Each kernel adaptive filtering algorithm is implemented as a Matlab class. To use one, first define its options:
```matlab
options = struct('nu',1E-4,'kerneltype','gauss','kernelpar',32);
```
Next, create an instance of the filter:
```matlab
kaf = aldkrls(options);
```
One iteration of training is performed by feeding one input-output data pair to the filter:
```matlab
kaf = kaf.train(x,y);
```
The outputs for one or more test inputs are evaluated as follows:
```matlab
Y_test = kaf.evaluate(X_test);
```

Example: time-series prediction
-------------------------------
Code from `demo/demo_prediction.m`
```matlab
% Demo: 1-step ahead prediction on Lorenz attractor time-series data
[X,Y] = kafbox_data(struct('file','lorenz.dat','embedding',6));

% make a kernel adaptive filter object of class aldkrls with options: 
% ALD threshold 1E-4, Gaussian kernel, and kernel width 32
kaf = aldkrls(struct('nu',1E-4,'kerneltype','gauss','kernelpar',32));

%% RUN ALGORITHM
N = size(X,1);
Y_est = zeros(N,1);
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end % progress indicator, 10 dots
    Y_est(i) = kaf.evaluate(X(i,:)); % predict the next output
    kaf = kaf.train(X(i,:),Y(i)); % train with one input-output pair
end
fprintf('\n');
SE = (Y-Y_est).^2; % test error

%% OUTPUT
fprintf('MSE after first 1000 samples: %.2fdB\n\n',10*log10(mean(SE(1001:end))));
```
Result:

    MSE after first 1000 samples: -40.17dB

Included algorithms
-------------------
- Approximate Linear Dependency Kernel Recursive Least-Squares (ALD-KRLS), as proposed in Y. Engel, S. Mannor, and R. Meir. "The kernel recursive least-squares algorithm", IEEE Transactions on Signal Processing, volume 52, no. 8, pages 2275-2285, 2004.
- Sliding-Window Kernel Recursive Least-Squares (SW-KRLS), as proposed in S. Van Vaerenbergh, J. Via, and I. Santamaria. "A sliding-window kernel RLS algorithm and its application to nonlinear channel identification", 2006 IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP), Toulouse, France, 2006.
- Naive Online Regularized Risk Minimization Algorithm (NORMA), as proposed in J. Kivinen, A. Smola and C. Williamson. "Online Learning with Kernels", IEEE Transactions on Signal Processing, volume 52, no. 8, pages 2165-2176, 2004.
- Kernel Least-Mean-Square (KLMS), as proposed in W. Liu, P.P. Pokharel, and J.C. Principe, "The Kernel Least-Mean-Square Algorithm," IEEE Transactions on Signal Processing, vol.56, no.2, pp.543-554, Feb. 2008.
- Fixed-Budget Kernel Recursive Least-Squares (FB-KRLS), as proposed in S. Van Vaerenbergh, I. Santamaria, W. Liu and J. C. Principe, "Fixed-Budget Kernel Recursive Least-Squares", 2010 IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP 2010), Dallas, Texas, U.S.A., March 2010.
- Kernel Recursive Least-Squares Tracker (KRLS-T), as proposed in S. Van Vaerenbergh, M. Lazaro-Gredilla, and I. Santamaria, "Kernel Recursive Least-Squares Tracker for Time-Varying Regression," Neural Networks and Learning Systems, IEEE Transactions on , vol.23, no.8, pp.1313-1326, Aug. 2012.
- Quantized Kernel Least Mean Squares (QKLMS), as proposed in Chen B., Zhao S., Zhu P., Principe J.C. "Quantized Kernel Least Mean Square Algorithm," IEEE Transactions on Neural Networks and Learning Systems, vol.23, no.1, Jan. 2012, pages 22-32.
- Random Fourier Fourier Feature Kernel Least Mean Squares (RFF-KLMS), as proposed in Abhishek Singh, Narendra Ahuja and Pierre Moulin, "Online Learning With Kernels: Overcoming The Growing Sum Problem", 2012 IEEE International Workshop on Machine Learning For Signal Processing.
- Extended Kernel Recursive Least Squares (EX-KRLS), as proposed in W. Liu and I. Park and Y. Wang and J.C. Principe, "Extended kernel recursive least squares algorithm", IEEE Transactions on Signal Processing, volume 57, number 10, pp. 3801-3814, oct. 2009.
- Gaussian-Process based estimation of the parameters of KRLS-T, as proposed in Steven Van Vaerenbergh, Ignacio Santamaria, and Miguel Lazaro-Gredilla, "Estimation of the forgetting factor in kernel recursive least squares," 2012 IEEE International Workshop on Machine Learning for Signal Processing (MLSP), 2012.
- Kernel Affine Projection algorithm with Coherence Criterion, as proposed in C. Richard, J.C.M. Bermudez, P. Honeine, "Online Prediction of Time Series Data With Kernels," IEEE Transactions on Signal Processing, vol.57, no.3, pp.1058,1067, March 2009.
- Kernel Normalized Least-Mean-Square algorithm with Coherence Criterion, as proposed in C. Richard, J.C.M. Bermudez, P. Honeine, "Online Prediction of Time Series Data With Kernels," IEEE Transactions on Signal Processing, vol.57, no.3, pp.1058,1067, March 2009.
- Recursive Least-Squares algorithm with exponential weighting (RLS), as described in S. Haykin, "Adaptive Filtering Theory (3rd Ed.)", Prentice Hall, Chapter 13.
- Multikernel Normalized Least Mean Square algorithm with Coherence-based Sparsification (MKNLMS-CS), as proposed in M. Yukawa, "Multikernel Adaptive Filtering", IEEE Transactions on Signal Processing, vol.60, no.9, pp.4672-4682, Sept. 2012.
- Parallel HYperslab Projection along Affine SubSpace (PHYPASS) algorithm, as described in M. Takizawa and M. Yukawa, "An Efficient Data-Reusing Kernel Adaptive Filtering Algorithm Based on Parallel Hyperslab Projection Along Affine Subspace," 2013 IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP), pp.3557-3561, May 2013.

Contributing
------------
If you wish to contribute to the toolbox, please [fork it on GitHub](https://github.com/steven2358/kafbox), push your change to a named branch, then send me a pull request.

License
-------
This source code is released under the FreeBSD License.
