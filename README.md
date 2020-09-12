# Speech Enhancement Using Wiener Filtering and PSR Phase Reconstruction

The algorithm estimates the magnitude spectrum of the underlying speech signal using Wiener filtering and the phase spectrum of the signal. The proposed method for phase reconstruction takes advantage of certain properties of the pitch-synchronous representation (PSR) of harmonic signals which make estimation of the phase spectrum faster and more accurate. The computation of the pitch-synchronous STFT requires fundamental frequency detection as a preliminary step.

The [document](https://github.com/geetkhatri/speech-enhancement-psr/blob/master/Speech%20Enhancement%20Using%20Wiener%20Filtering%20and%20PSR%20Phase%20Reconstruction.pdf) explains how the algorithm works and the motivation behind it.

### References

The MATLAB code for the implementation of the algorithm makes use of:

- Mike Brookes. VOICEBOX: A speech processing toolbox for MATLAB (http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html).
- Pascal Scalart. Wiener filter for Noise Reduction and speech enhancement (https://www.mathworks.com/matlabcentral/fileexchange/24462-wiener-filter-for-noise-reduction-and-speech-enhancement), MATLAB Central File Exchange.
