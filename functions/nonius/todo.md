# General

- [ ] make a script for calibration/validation in the TMS setup (to be ran before each run)
- [ ] it would be good to have the option to launch a GUI to visualize the different runs. But that's too much I fear.
- [ ] write the docs and examples
- [ ] replace the shadedErrorBar function

# New features

- [ ] compute vergence and version, and compare with theoretical vergence (this could also be used during screening)
- [ ] manual confirmation of outlier rejection (mark or unmark as outlier)

# Improvements

- [ ] arg checks
- [ ] unit tests
- [ ] assert time window within recording
- [ ] make sure there are no spaces in epochName
- [ ] in the heatmap function could plot a visual angle reference axis using concentric circles overlaid in the image
- [ ] ? assert the movement delay to be constant across the eyes ?
- [ ] ? add option to pass design as input to estimateCalibrationParams and then call setCalibrationDesign within the former method
- [ ] denoise and outlier rejection for calibration data (like removing blinks and stuff). Now I have a dataset to try that.
- [ ] apply trf to camera coordinate measuremnts? need to rescale the calibration grid to match the measurements.

# Bugs


