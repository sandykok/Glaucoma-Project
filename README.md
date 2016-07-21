# Glaucoma-Project

This project objectives are
To localize and detect optic disk.
To segment and extract the optic cup and optic disk.
To find cup to disk ratio.
To classify the images based on the features extracted using SVM classifier

Algorithm

Modules
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

