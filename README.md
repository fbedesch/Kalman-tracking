# Kalman-tracking
This code contains a standalone version of the tracking code included in DELPHES.
The calculation of the track covariance matrix can be done either with standard full covariance matrix of the measurements, including multiple scattering effects, or with a Kalman filter approach.
Running instructions:

root

root> .L examples/LoadTrk.c+

root>LoadTrk()

root> KalmanCheck("./Geometries/GeoIDEA_NewDCH.txt")

This program will generate plots using both the standard version or the Kalman version of the code.

root> PlotGeo("./Geometries/GeoIDEA_NewDCH.txt")

will make a simple plot of the geometry and display the material as a function of the polar angle or its cosine.
