# Canny Edge Detection
This project demonstrates the use of a custom implementation of the Canny Edge Detection technique to identify potential tumors in MRI brain images.

# Algorithm
- Input Image Preprocessing: Load the MRI brain image and convert it to grayscale for analysis.
- Gaussian Smoothing: Apply a Gaussian filter to remove noise and smooth the image (sigma = 3).
- Gradient Computation: Compute image gradients (img_dx, img_dy) using Gaussian derivatives to obtain the magnitude of edges.
- Non-Max Suppression: Retain local maxima in the gradient magnitude to thin the edges.
- Double Thresholding: Classify pixels as strong edges, weak edges, or non-edges based on the high and low threshold ratios (highThresholdRatio = 0.275, lowThresholdRatio = 0.25).
- Edge Tracking by Hysteresis: Connect weak edges to strong edges recursively to form complete boundaries.
- Tumor Edge Detection: The final output highlights tumor-like regions in the MRI brain image.

# Output Visualizations
The algorithm provides the following outputs at each stage:
- Smoothed image after Gaussian filtering.
- Gradient magnitude visualization.
- Non-Max Suppression result.
- Double thresholding result.
- Final edge-detected MRI image with highlighted tumor regions.
