# Tadpole Behavior Automation
This is a Matlab scrip and function for use in automating the avoidance behavior analysis of Xenopus tadpoles. To run this script it is necessary to download the BatchRunTAD9000.m and the TadFunctionTest.m. Currently the Tad_Det_Assign_9001.m is the test file kept up to date with the TadFunctionTest.m that is used for testing and confirmation purposes.

### Program Goals:

This script was designed to be used for expedited and automated analysis of visual response behavior of Xenopus tadpoles to a moving dot stimulus. The program uses blob filtering methods to find the locations of the tadpoles then assigns each tadpoles location based on distance to the next location. The location and radius of each oncoming dot is then extracted from the video. Using the location of the tadpoles and the location and radius of the dots the data is analyzed to find the frames in which a tadpole comes in contact with a dot. Based on set criteria from previous analysis protocols avoidances, in relation to the dot, of the tadpole and encounters with a dot are counted. This provides useful data used to analyze recovery behavior. 

### Pre-Run Dependencies:

- https://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m--extrema2-m
- https://www.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker

To run the BatchRunTAD9000.m script or the TadFunctionTest.m function the two above dependencies must be installed off the Matlab File Exchange. The Extrema2 package is used for finding the locations of the local maximums for detection of tadpoles, and the Simple Tracker package is used for the Munkres assignment algorithm. These packages must be installed on a computer running Matlab version R2017_a or later with the Image Processing Toolbox installed and at least 16 Gb of memory.

### Running Program:

Once installed the BatchRunTAD9000.m should be opened in the same folder as the TadFunctionTest.m. From the BatchRunTAD9000.m script hit run and select the folder containing the videos taken during the stimulus trials. This will begin the video analysis process creating folders for each video run saving location and avoidance data throughout the process. Videos with detections of less than half of the detections in the first frame analyzed will be skipped and shown as an error. 

### Situational Adjustments:
- Size of LoG filter
- Threshold value for detections
  - Increasing the threshold value will cause more noise from the original pictures to be discarded below the specified threshold value. This will leave the deeper points, tadpole gut, found by the LoG filter as maximums allowing the program to search for all local maxima and save the maxima as X and Y coordinates. If there are too many detections that aren’t tadpoles this can be increased incrementally until a reasonable value is reached to eliminate excess detections. 
  
- Tadpole speed variability and, X and Y noise
  - The speed variability is taken as an estimate of the tadpoles speed from frame to frame given in pixels/frame. The X and Y noise are also taken as estimates of how much potential there is for tadpole movement variability in the detections. These three values should not have to be changed as the default values are accurate. 
  
- Area to crop out of video for dot detection
  - These values are given as the area to crop out of the video when processing the images for finding the locations of the dots. These values are hard coded into the script and cannot be changed via the initial prompt so to adjust the cropping must be done directly in the function. The values are xleft, xright, ytop, and ybot. Increasing xleft and decreasing xright result in the picture getting smaller in the x direction. Increasing ytop and decreasing ybot result in the picture getting smaller in the y direction. 
  
- Distance between tadpole gut and eyes
  - The program uses the gut of the tadpole to find the location because it is the most visible part of the tadpole. This being said estimated coordinates of the tadpoles field of view needed to be found in order to assure that the locations of the dots crossed over the tadpoles eyes. After testing it was determined that the average distance between tadpoles gut, at this stage of growth, and eyes was about 15 pixels. Using this as a radius points in a semicircle around the head and eyes are found and taken as new coordinates for the tadpoles. Adjusting this value could be done either by increasing or decreasing this radius for a tadpole at either a later or earlier stage of growth respectively. 
  
- Angle tolerance for avoidance determination
  - Counting an avoidance event is dependent on how much a tadpoles swimming direction angle changes within eight frames of crossing the location of a dot. To quantify this turn upper and lower bounds for a change in swimming direction angle were set at 110° and 70° respectively. This means that if a tadpole is swimming, in the video, to the right an avoidance is counted if the tadpoles swimming angle changes greater than or equal to 70° of the angle it was previously going. Now, for a tadpole moving to the left an avoidance is counted if the tadpoles swimming angle is less than or equal to 110° of the angle it was previously at. The values for upper and lower avoidance determination can be changed by increasing the upper avoidance angle and decreasing the lower avoidance angle for a smaller turn to be counted as an avoidance. To increase the size of the turn required to count an avoidance the upper and lower angles can be decreased and increased up to 95° and 85° respectively.

### Screenshots:
![tadpoleoriginal](https://user-images.githubusercontent.com/41453184/43158936-dd0b34a2-8f35-11e8-88a6-6bb53058d6ca.png)
The original image from the behavior videos before any image processing is shown above. 

![tadpoledotsandassignments](https://user-images.githubusercontent.com/41453184/43158847-9c6ff4d2-8f35-11e8-905a-e43e728ef397.png)
After tadpole detection and assignment and, dot detection the overlay of the extracted location data can be seen visually above. 
