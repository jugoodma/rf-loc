# 798 Proj

Using the `spotfi` for angle of arrival of `ofdm` signal packets.
Angle measurements will be taken by a moving car, at a calculated delta distance.
Using this, we can triangulate the position of the transmitter using angle-side-angle from geometry.

* `parameters.m` contains `ofdm` and general code parameters
* `calibration.m` transmits a sine wave to test if your antennas are calibrated
* `capture.m` transmits an `ofdm` packet continuously (with a very small buffer), and saves the recording to a `.dat` file
* `match.m` analyzes a `.dat` recording and computes AoA (and ToF) estimates per packet, as well as per symbol in the packet, using the channel estimate `H`
