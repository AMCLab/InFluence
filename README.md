# InFluence

InFluence converts unmodulated simulated TEM images to modulated images with a user defined flux.

The Stobbs factor is the discrepancy between simulated images and experimental images. Its origin is instrument related, and can therefore be rectified by including the effect of the imaging system on the observed image intensities. Electrons scatter in the sensor layer of direct detectors; the lateral spread of the electrons can be far greater than the pixel size, resulting in the detection of electrons in pixels neighbouring the pixel of incidence. InFluence simulates modulation of the expected intensities due to the lateral spread of electrons by employing a Monte Carlo single scattering model to calculate the trajectories of the electrons in the sensor.

InFluence has been designed around the Medipix3 detector operating in single pixel mode, but we expect that it can be applied to other detectors that have a similar design (single sensor layer sitting on pixel readout circuitry).

Usage of the scripts is given in the Usage_Guide.txt document. Example unmodulated images can be found in the ExampleData folder. The slanted edge is oriented at approximately 0.1 rad to the x-axis.

Note that these scripts will be refactored in the future to be more modular for easy integration with other simulation packages, so usage will change.

InFluence: An Open-Source Python Package to Model Images Captured with Direct Electron Detectors, 
Microscopy and Microanalysis, Volume 29, Issue 4, August 2023, Pages 1380–1401,
https://doi.org/10.1093/micmic/ozad064
