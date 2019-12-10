# TwoColorUnmixing
Code for correlation analysis of two-color imaging data with spectral unmixing

Coded and tested with Matlab R2018b and Windows 10.0
Required toolboxes:
Statistics and Machine Learning Toolbox

No installation is required; just navigate Matlab's working directory to this folder.
Install time: 1 min
Expected run time depends on selected parameter settings, >1 min per 1000 shuffle instances on a typical computer.

The purpose of this code is to demonstrate the process by which null distributions for sample correlations were generated in (Fig 2e-g) of our manuscript describing jYCaMP. Two-color imaging results in tradeoffs between collection efficiency and crosstalk. To have reasonable collection efficiencies, it is normal to use filter sets that result in small amounts of crosstalk (~10%) between emitters, for example when imaging GFP and RFP. Like several studies before us, we wanted to use correlation between two functional indicators to detect interacting pre- and post-synaptic compartments. However, to perform appropriate statistical tests, we need to account for the effect of bleedthrough on sample correlations.

One common approach to compensate for bleedthrough is linear unmixing, which is appropriate for Gaussian-distributed intensities. However, in functional imaging photon counts are Poisson-distributed, and unmixing (either least-squares or maximum likelihood under the Poisson model) results in biased sample correlations. Least squares spectral unmixing of Poisson-distributed measurements results in negative bias in computed sample correlations. Maximum likelihood unmixing of Poisson-distributed measurements results in positive bias; we did not use maximum-likelihood unmixing in our study. We found that we could compensate for either bias and recover accurate null distributions (to a precision much finer than measured effects) by a process of bootstrap sampling followed by median subtraction.

To validate that approach, this code performs the sampling process for simulated ground truth data with a variety of distributions. It reports the 2-sample KS statistic that was used as a test statistic in our study.
