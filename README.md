# Tadpole Behavior Automation
This is a Matlab scrip and funtion for use in automating the avoidance behavior analysis of Xenopus tadpoles. To run this script it is necessary to download the BatchRunTAD9000.m and the TadFunctionTest.m. Currently the Tad_Det_Assign_9001.m is the test file kept up to date with the TadFunctionTest.m that is used for testing and confirmation purposes.

### Program Goals:

This script was designed to be used for expedited and automated analysis of visual response behavior of Xenopus tadpoles to a moving dot stimulus. The program uses blob filtering methods to find the locations of the tadpoles then assigns each tadpoles location based on distance to the next location. The location and radius of each oncoming dot is then extracted from the video. Using the location of the tadpoles and the location and radius of the dots the data is analyzed to find the frames in which a tadpole comes in contact with a dot. Based on set criteria from previous analysis protocols avoidances, in relation to the dot, of the tadpole and and encounters with a dot are counted. This provides useful data used to analyze recovery behavior. 

### Pre-Run Dependencies:

- https://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m--extrema2-m
- https://www.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker

To run the BatchRunTAD9000.m script or the TadFunctionTest.m function the two above dependencies must be installed off the Matlab File Exchange. The Extrema2 package is used for finding the locations of the local maximums for detection of tadpoles, and the Simple Tracker package is used for the Munkres assignment algorithm. These packages must be installed on a computer running Matlab version R2017_a or later with the Image Processing Toolbox installed and at least 16 Gb of memory.

### Running Program:

Once installed the BatchRunTAD9000.m should be opened in the same folder as the TadFunctionTest.m. From the BatchRunTAD9000.m script hit run and select the folder containing the videos taken during the stimulus trials. This will begin the video analysis process creating folders for each video run saving location and avoidance data throughout the process. Videos with detections of less than half of the detections in the first frame analyzed will be skipped and shown as an error. 

### Situational Adjustments:
- Depth of Gaussian filter
- Threshold value for detections
- Tadpole speed variability and, X and Y noise
- Acceptable distance for assignment based on detections
- Area to crop out of video for dot detection
- Distance between tadpole gut and eyes
- Angle tolerance for avoidance determination

### Screenshots:
![tadpoleoriginal](https://user-images.githubusercontent.com/41453184/43158936-dd0b34a2-8f35-11e8-88a6-6bb53058d6ca.png)
The original image from the behavior videos before any image processing is shown above. 

![tadpoledotsandassignments](https://user-images.githubusercontent.com/41453184/43158847-9c6ff4d2-8f35-11e8-905a-e43e728ef397.png)
After tadpole detection and assignment and, dot detection the overlay of the extracted location data can be seen visually above. 
