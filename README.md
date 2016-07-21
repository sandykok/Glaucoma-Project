# Glaucoma-Project

This project objectives are
a. To localize and detect optic disk.
b. To segment and extract the optic cup and optic disk.
c. To find cup to disk ratio.
d. To classify the images based on the features extracted using SVM classifier

Algorithm

1.  Optical Disc size estimation,
2.  Optic Disk localization,
	a. Background normalization,
	b. Template matching,
	c. Directional matched filtering,
3.  Optic Disk segmentation.
	a. Image Preprocessing,
		i. Blood vessel removal,
		ii. Bright region removal,
	b. Level Set Model,
	c. Least-Square Ellipse Fitting
4. Optic Cup segmentation.
5.Cup to Disk Calculation.
6.Classification.
	a. Feature Extraction.
	b. Classification using Cup to Disk.
	c. SVM Classification.

Results

Specificity - 94.28%
Sensitivity - 95%
Accuracy - 94.61%
