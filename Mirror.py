import numpy as np
import matplotlib.pyplot as plt
import collections
import math
import matplotlib
from Fits import FitsFile
from Segment import Segment
from scipy.ndimage import gaussian_filter

class Mirror:
    """
    Class Mirror:
        This class represents a mirror in a MIRA Image. It provides functionalities for analyzing mirror segments, cleaning the mirror data, and visualizing mirror properties.

        Attributes:
            radian_to_arcsec (float): Conversion factor from radians to arcseconds.
            sqrt_of_3_div_2 (float): Square root of 3 divided by 2.
            sensor_offset (float): Offset value for sensor.
            segment_side_a (float): Length of one side of a mirror segment.
            denom (float): Denominator value used for calculations.
            d2_arcsec_scale (float): Scaling factor for converting radians to arcseconds.
            n_segments (int): Number of mirror segments.
            n_actuators (int): Number of actuators.
            n_sensors (int): Number of sensors.
            n_edges (int): Number of edges.
            file (FitsFile): Object representing a FITS file containing mirror data.
            center (tuple): Coordinates of the center of the mirror.
            segments (list): List of Segment objects representing mirror segments.
            plate_scale (float): Plate scale of the mirror.
            scale (float): Scaling factor based on plate scale and other parameters.
            size (int): Size of the boxes around each mirror segment.

        Methods:
            getSegmentCoordinates(segNr): Returns the coordinates of a given segment number.
            getTxM(ang): Returns the transformation matrix for a given angle.
            rotate_coordinates(segNr, ang): Rotates the coordinates of a segment by a given angle.
            clean_mirror(debug=False, min_distance=10): Cleans the mirror segments based on distance criteria.
            visualizeCentroiding(x, y, size=6): Visualizes centroiding by marking the coordinates and search area on an image.
            plot_mirror(): Plots the mirror segments on an image.
            plot_blurred_image(sigma=20): Plots the blurred image of the mirror data.
            find_center(plot=True, debug=False): Finds the center of the mirror based on image data.
            find_mirror_segments(size=6, sigma=10, debug=True, plot=True, attempt=1): Finds mirror segments in the image.
            find_new_center(data): Finds the new center based on image data.
            goodFit(points): Checks if a set of mirror segment points constitute a good fit.
            construct_new_mirror(debug=False): Constructs a new mirror based on existing segments.
            remove_duplicate_segments(): Removes duplicate mirror segments.
            move_point(segment, angle_degrees, distance): Moves a mirror segment point by a given angle and distance.
            close_angles(angle1, angle2, threshold=3): Checks if two angles are close within a threshold.
            calculate_center(): Calculates the center of the mirror based on existing segments.
            recenter(box): Recenters the mirror segments within a specified box size.
            remove_less_than_3_neighbors(scale): Removes mirror segments with less than three neighbors within a specified scale.
            calc_rotation_angle(point1, point2): Calculates the rotation angle between two points.
            dist_between(point1, point2): Calculates the distance between two points.
    """

    def __init__(self, path):
        # Constants
        self.radian_to_arcsec = 3600.0 * 180.0 / math.pi
        self.sqrt_of_3_div_2 = math.sqrt(3) / 2
        self.sensor_offset = 0.055
        self.segment_side_a = 0.9
        self.denom = self.sensor_offset * self.segment_side_a * self.sqrt_of_3_div_2
        self.d2_arcsec_scale = 1E-9 * self.radian_to_arcsec / self.denom
        self.n_segments = 36
        self.n_actuators = 108
        self.n_sensors = 168
        self.n_edges = 84
        # Variables from Image
        self.file = FitsFile(path)
        self.rot_angle = self.file.rot_angle #if the rotation angle is slightly off we want to know the new and original angles
        self.center = (len(self.file.data[0]) // 2,(len(self.file.data) // 2)) # center of the mirror
        self.segments = [] # array of segments
        self.plate_scale = self.file.instrument_data[self.file.instrument]['PlateScale'] # Plate Scale 
        self.scale = (1E-9 * self.radian_to_arcsec / self.plate_scale / self.denom * self.file.pmfm)
        self.seg_size = int(1.2 / (self.plate_scale+.03)) # Size of the boxes around each mirror segment
        if self.file.instrument == 'KPF' or self.file.instrument == 'KCWI':
            self.seg_size += 3
        

    def getSegmentCoordinates(self, segNr):
        """
        Calculates the coordinates of a segment.
        Parameters:
        - segNr (int): Segment number ranging from 0 to 35.
        Returns:
        - list: A list containing the x and y coordinates of the segment.
        """
        
        segmPosXY = [
            [1, 1],  [0, 1],   [-1, 1],  [-1, -1], [0, -1],  # 1-5
            [1, -1], [2, 1],   [1, 2],   [0, 2],   [-1, 2],  # 6-10
            [-2, 1], [-2, 0],  [-2, -1], [-1, -2], [0, -2],  # 11-15
            [1, -2], [2, -1],  [2, 0],   [3, 2],   [2, 2],  # 16-20
            [1, 3],  [0, 3],   [-1, 3],  [-2, 2],  [-3, 2],  # 21-25
            [-3, 1], [-3, -1], [-3, -2], [-2, -2], [-1, -3],  # 26-30
            [0, -3], [1, -3],  [2, -2],  [3, -2],  [3, -1], [3, 1]]  # 31-36
        result = []       
        segmentSideH = self.segment_side_a * self.sqrt_of_3_div_2
        offset = -segmentSideH
        ix = segmPosXY[segNr][0]
        iy = segmPosXY[segNr][1]
        
        if (ix & 1) == 0: # Adjust offset for even row segments
            offset = 0 # Mirror pattern and the reason there is an offset--mirror centers are not exactly 2 mirror radii away
        
        result.append(ix * (self.segment_side_a + self.segment_side_a / 2)) # Calculate x coordinate
        
        if iy < 0:
            offset = -offset # Adjust offset for negative y values
        
        result.append(iy * 2 * segmentSideH + offset) # Calculate y coordinate
        return result


    def getTxM(self, ang):
        """
        Returns the 2x2 transformation matrix for a given angle of rotation.
        Parameters:
        - ang (float): Angle of rotation in degrees.
        Returns:
        - list: 2x2 transformation matrix representing the rotation.
        """
        
        ca = math.cos(math.radians(ang))
        sa = math.sin(math.radians(ang))
        txM = [[ca, -sa], [sa, ca]]
        return txM


    def rotate_coordinates(self, segNr, ang):
        """
        Rotates the coordinates of a segment rotated around the origin by a given angle.
        Parameters:
        - segNr (int): Index of the segment.
        - ang (float): Angle of rotation in degrees.
        Returns:
        - list: List containing the rotated coordinates of the segment.
        """
        
        txM = self.getTxM(ang)
        result = self.getSegmentCoordinates(segNr) # Get the coordinates of the segment
        xa, ya = result[0], result[1]
        result[0] = float(txM[0][0] * xa + txM[0][1] * ya) # Multiply the coordinates by the rotation matrix 
        result[1] = float(-txM[1][0] * xa - txM[1][1] * ya)
        return result
    
    
    def clean_mirror(self, debug=False, min_distance=10):
        """
        Removes segments that are closer than a specified minimum distance.
        Parameters:
        - debug (bool, optional): If True, prints debug information. Defaults to False.
        - min_distance (int, optional): The minimum distance between segments to be maintained. Defaults to 10.
        """
        
        min_distance = abs(self.scale)
        if debug:
            print(f"Before Clean: {len(self.segments)}")
            print(f"Removing Segments < {min_distance}")
        
        self.segments = sorted(self.segments)
        dist = len(self.segments) - 1
        removes = []
        for point in self.segments:
            proceed = True
            for rem in removes: 
                if point == rem: # Do not remove both points that are close to eachother
                    proceed = False
            if proceed:
                for other_point in self.segments:
                    d = self.dist_between(point, other_point)
                    if d < min_distance and d != -1:
                        removes.append(other_point)
        for point in removes:
            try:
                self.segments.remove(point)
            except Exception as e:
                if debug:
                    print(e)
        if debug:
            print(f"After Clean: {len(self.segments)}")
           
        
    #-----------------------------------------------------#
    #---------------Visualization Functions---------------#
    #-----------------------------------------------------#
    def visualizeCentroiding(self, x, y, size=6):
        """
        Visualizes centroiding by displaying an image with a marked search area around the centroid.
        Parameters:
        - x (int): x-coordinate of the centroid.
        - y (int): y-coordinate of the centroid.
        - size (int, optional): Size of the search area around the centroid. Defaults to 6.
        """
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
        image_data = self.file.data[y - size:y + size,x - size:x + size]
        idata = gaussian_filter(image_data,2)
        ax2.imshow(idata)
        ax1.scatter(x,y, s = 1)
        r = matplotlib.patches.Rectangle((x-size,y-size),size*2,size*2, fill = False)
        ax1.add_patch(r)
        ax1.imshow(image_data)
        ax1.set_title("Coordinate & Search Area Marked on the Image")
        ax2.set_title("Image Zoomed to the Box")
        
        
    def plot_mirror(self):
        """
        Plots the mirror segments that have been found on an image.
        """
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
        rows, cols = np.shape(self.file.data)
        rows, cols = rows / 2, cols / 2
        scale = int(self.scale)
        for segment in self.segments:
            ax2.scatter(segment.x, segment.y, color = 'red', s = 2)
        ax2.imshow(self.file.data, cmap = 'gray')
        ax1.imshow(self.file.data, cmap = 'gray')
        plt.title('Image')
        plt.show()
        
        
    def plot_blurred_image(self, sigma=20):
        """
        Plots the original image and the image with a specified blur factor. 
        """
        
        blurred_image_data = gaussian_filter(self.file.data, sigma)# Apply Gaussian filter for blurring | Adjust sigma for the desired blur strength
        plt.figure(figsize=(15, 7)) # Plot the original and blurred images side by side
        plt.subplot(1, 2, 1) # Plot the original image
        plt.imshow(self.file.data, cmap='gray')
        plt.title('Original Image')
        plt.subplot(1, 2, 2) # Plot the blurred image
        plt.imshow(blurred_image_data, cmap='gray')
        plt.title('Blurred Image')
        plt.show() # Display the plot

        
    #-----------------------------------------------------#
    #--------------Initial Mirror Acquisition-------------#
    #-----------------------------------------------------#
    def find_center(self, plot=True, debug=False):
        """
        Finds the center of the mirror image.
        Parameters:
        - plot (bool): Whether to plot the center-finding process.
        - debug (bool): Whether to print detailed information about the process.
        """

        sigma1, sigma2, sigma3 = 40,50,40 # Set the amount of blur to apply to each of the three subsets of image data
        # Adjust these parameters based on the size of the mirror image and amount of background noise
        # More noise= more blur needed
        dimX, dimY = len(self.file.data), len(self.file.data[0])
        box1 = 100
        box2 = int(min(16/self.plate_scale, (dimX / 2) - 25)) # Side length of the last two boxes used to fit the center
        box3 = int(min(14/self.plate_scale, (dimX / 2) - 80)) # These boxes are based on the size of the mirror
        if self.file.instrument == 'KCWI':
            sigma1, sigma2, sigma3 = 120,80,20
            box1 = 1
            box3 = int(min(3/self.plate_scale, (dimX / 2) - 80))  
        if self.file.instrument == 'KPF': # Kpf needs a different blur method due to its gradient
            sigma2, sigma3 = 80,15
            box3 = int(min(3/self.plate_scale, (dimX / 2) - 80))  
        if debug:
            print(f"First Blur Size: {box2} | Second Blur Size: {box3}")
            print(f"This is a {dimX} x {dimY} image")
        #------------------------------------------------------------------------------#
        #------------------------------Find Center Try 1-------------------------------#
        #------------------------------------------------------------------------------#
        #------------All find center attempts follow the same methodology--------------#
        #------------1. Take a subset of the total image and blur it-------------------#
        #------------2. Find the brightest (or dimmest) spot and make that the center--#
        #------------3. Return the new center or keep going with a smaller subset------#
        #------------------------------------------------------------------------------#
        
        if dimX < 400 or dimY < 400:
            box1 = 50
            sigma2 = 25
            sigma3 = 30
        data_new = self.file.data[box1:-box1, box1:-box1] # 1. Take the image subset (these get smaller as we get closer to the actual center
        data_new = gaussian_filter(data_new, sigma1)  # 1. Blur the image (datanew will sum to find the brightest 'y')
        data_t = data_new.transpose() # data_t will sum to find the brightest 'x'
        # 2. Find the brightest spot
        sumx, sumy = [], []
        for r in range(min(len(data_new), len(data_t))):
            sumx.append(sum(data_t[r]))
            sumy.append(sum(data_new[r]))
        midy, midx = self.center[0], self.center[1]
        max1X = max(sumx, default=midx) 
        max1Y = max(sumy, default=midy)
        # 2.(cont) Make that the center
        xNew = sumx.index(max1X) + box1
        yNew = sumy.index(max1Y) + box1
        if dimX < 400 or dimY < 400:
            threshold = box1
            if xNew > (midx + threshold):
                xNew= midx + threshold
            elif xNew < (midx - threshold):
                xNew = midx - threshold
            if yNew > (midy + threshold):
                yNew = midy + threshold
            elif yNew < (midy - threshold):
                yNew = midy - threshold
        if debug:
            print(f"Center1: {xNew}, {yNew}")
        #------------------------------------------------------------------------------#
        #------------------------------Find Center Try 2-------------------------------#
        #------------------------------------------------------------------------------#
        xNew1, yNew1 = xNew, yNew
        if self.file.instrument != "KCWI":
            if yNew - box2 < 0:
                box2 = yNew
                # print(f"Edge: {box2}")
            if xNew - box2 < 0:
                box2 = min(yNew, xNew)
                # print(f"Edge: {box2}")
            if yNew + box2 > dimY:
                box2 = dimY - yNew
                # print(f"Edge: {box2}")
            if xNew + box2 > dimX:
                box2 = min(dimY - yNew, dimX - xNew)
                # print(f"Edge: {box2}")
            data_new1 = self.file.data[yNew - box2: yNew + box2, xNew - box2: xNew + box2]
            # if len(data_new1) == 0: 
            #     # This happens when the new center coordinates are too close to the edge
            #     # and the box thats created goes off of the image. In that case just use the center of the image
            #     # because most of the mirror images should be relatively close to the center of the image
            #     yNew, xNew = int(len(self.file.data) / 2), int(len(self.file.data[0]) / 2)
            #     data_new1 = self.file.data[yNew - box2: yNew + box2, xNew - box2: xNew + box2]
            data_new1 = gaussian_filter(data_new1, sigma2)
            data_t = data_new1.transpose()
            # Find the maximum of the smaller subset of the image
            sumx, sumy = [], []
            try:
                for r in range(len(data_new1)):
                    sumx.append(sum(data_t[r]))
                    sumy.append(sum(data_new1[r]))
                midy, midx = int(len(data_new1 / 2)), int(len(data_new1[0] / 2)) 
                max1X = max(sumx, default=midx)
                max1Y = max(sumy, default=midy)
                xNew1 = sumx.index(max1X) + xNew - box2
                yNew1 = sumy.index(max1Y) + yNew - box2
                if debug:
                    print(f"Center2: {xNew1}, {yNew1}")
            except:
                xNew1 = xNew
                yNew1 = yNew
                
            threshold = box2
            if xNew1 > (xNew + threshold):
                xNew1= xNew + threshold
            elif xNew1 < (xNew - threshold):
                xNew1 = xNew - threshold
            if yNew1 > (yNew + threshold):
                yNew1 = yNew + threshold
            elif yNew1 < (yNew - threshold):
                yNew1 = yNew - threshold
        #------------------------------------------------------------------------------#
        #------------------------------Find Center Try 3-------------------------------#
        #------------------------------------------------------------------------------#
        if yNew1 - box3 < 0:
            box3 = yNew1
        if xNew1 - box3 < 0:
            box3 = min(yNew1, xNew1)
        if yNew1 + box3 > dimY:
            box3 = dimY - yNew1
        if xNew1 + box3 > dimX:
            box3 = min(dimY - yNew1, dimX - xNew1)
        data_new2 = self.file.data[yNew1 - box3: yNew1 + box3, xNew1 - box3: xNew1 + box3]
        data_new2 = gaussian_filter(data_new2, sigma3)
        data_t = data_new2.transpose()
        # Find the final max
        sumx, sumy = [], []
        for r in range(len(data_new2)):
            sumx.append(sum(data_t[r]))
            sumy.append(sum(data_new2[r]))
        midy, midx = int(len(data_new2) / 2), int(len(data_new2[0]) / 2) 
        if self.file.instrument == 'KPF' or self.file.instrument == 'KCWI':
            if debug:
                print("KPF/KCWI CENTER")
            # Since the center of kpf and kcwi images are big enough to have a dim spot we can find a more accurate
            # center by finding the darkest point on the center of the image
            # this corrects for kpf gradient, which is when one side of the mirror is brighter and skews the center
            miny, minx = np.unravel_index(np.argmin(data_new2, axis=None), data_new2.shape)
            xNew2 = minx + xNew1 - box3
            yNew2 = miny + yNew1 - box3
        else:
            max1X = max(sumx, default=midx)
            max1Y = max(sumy, default=midy)
            xNew2 = sumx.index(max1X) + xNew1 - box3
            yNew2 = sumy.index(max1Y) + yNew1 - box3
        threshold = abs(int(self.scale)) # We can only move 1 mirror segment in any direction (this avoids overcorrecting)
        if xNew2 > (xNew1 + threshold):
            xNew2 = xNew1 + threshold
        elif xNew2 < (xNew1 - threshold):
            xNew2 = xNew1 - threshold
        if yNew2 > (yNew1 + threshold):
            yNew2 = yNew1 + threshold
        elif yNew2 < (yNew1 - threshold):
            yNew2 = yNew1 - threshold
        if debug:
            print(f"Center3: {xNew2}, {yNew2}")
            print(f"Sigmas: {sigma1, sigma2, sigma3}")
        if plot:
            fig, (ax5, ax6, ax7, ax8) = plt.subplots(1, 4, figsize=(18, 9))
            ax5.imshow(self.file.data)
            ax5.scatter(xNew, yNew, s=20, color='red')
            ax5.scatter(xNew1, yNew1, s=20, color='orange')
            ax5.scatter(xNew2, yNew2, s=20, color='green')
            ax6.scatter(xNew - 100, yNew - 100, s=20, color='red')
            ax7.scatter(xNew1 - xNew + box2, yNew1 - yNew + box2, s=20, color='orange')
            ax8.scatter(xNew2 - xNew1 + box3, yNew2 - yNew1 + box3, s=20, color='green')
            ax6.imshow(data_new)
            try:
                ax7.imshow(data_new1) 
            except: 
                pass
            ax8.imshow(data_new2)
        self.center = xNew2, yNew2
        
      
    def find_mirror_segments(self, size=6, debug=False, plot=True, attempt=1):
        """
        Finds the mirror segments and fits them to the image.
        Parameters:
        - size (int): Size of the mirror segments.
        - debug (bool): Whether to print detailed information.
        - plot (bool): Whether to plot the fitting process.
        - attempt (int): Attempt number.
        """
        try:
            self.clean_mirror
        except:
            pass
        if len(self.segments) != 36:
            xCenter, yCenter = self.center[0], self.center[1] # There is an option to move the center around, save the new centers as temp variables
            data = gaussian_filter(self.file.data, 2)

            if plot: # If we are plotting the segments then initialize the subplots.
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9)) # Axis 1 contains the original point guesses of the mirror | Axis 2 contains the fitted mirror segments    

            if debug: 
                print(self.file) # Describes the details of the file that we are fitting mirror segments to

            size = self.seg_size # Create a soze variable that we can change
            radius = self.scale

            if attempt == 1: # 
                angles = [self.file.rot_angle, -self.file.rot_angle]
            else:
                angles = [self.file.rot_angle, self.file.rot_angle - 10, self.file.rot_angle - 20, self.file.rot_angle - 30, self.file.rot_angle - 40, self.file.rot_angle - 50]
            m = int(math.sqrt(math.sqrt(self.scale**2)))
            centers = [(0,0)]
            all_points, bestPoints = [],[]
            best_range = 9999
            fit = False # Fit is false until we find a good fit for the mirror segments

            for center in centers:
                xCenter, yCenter = (xCenter + center[0]), (yCenter + center[1])

                for angle in angles: # Loop through all of the rotation angles
                    coordinates = []

                    for i in range(self.n_segments):  # Loop through all of the mirror segments
                        size = self.seg_size
                        coords = self.rotate_coordinates(i, angle) # After rotating the coordinates, center them around the actual mirror center and scale them by the 'scale' value
                        coords[0] = int(coords[0] * radius) + xCenter
                        coords[1] = int(coords[1] * radius) + yCenter
                        point = coords
                        x = point[0]
                        y = point[1]
                        mirror = self.file.data[y - size:y + size,x - size:x + size] # Take a subset of the image
                        if len(mirror) == 0: # Only happens when the mirror segment box goes off the image, which means our image is messed 
                            xC, yC = 0,0
                        else:
                            xC, yC = self.find_new_center(mirror) # Find the new point estimate for that mirror segment
                        # Update x and y accordingly
                        x, y = (x + xC - size), (y + yC - size)
                        # Shrink the search box (we should be somewhere on the mirror now)
                        size -= 1
                        # Redo the fitting process for a better fit on the mirror segment
                        mirror = self.file.data[y - size:y + size,x - size:x + size]
                        if len(mirror) == 0:
                            xC, yC = 0,0
                        else:
                            xC, yC = self.find_new_center(mirror)
                        x, y = (x + xC - size), (y + yC - size)
                        # Save the updated coordinates
                        s = Segment(x,y)
                        coordinates.append(s)
                    # Check to see if the updated 36 segments are a good fit
                    all_points.append(coordinates)
                    fit, range_new = self.goodFit(coordinates)
                    # If they are a good fit then update the range and update the coordinates (coordinates are the best center of each of the 36 segments)
                    if fit:
                        bestPoints = coordinates
                        best_range = range_new
                        self.rot_angle = angle
                        break
                    # If they are not a good fit but still a better fit then update our 'best' fit. 
                    # Even if none of the coordinates are a good fit we still want to save the best fit (might only be 1 or 2 segments away)
                    elif range_new < best_range:
                        best_range = range_new
                        bestPoints = coordinates

            self.segments = bestPoints
            self.clean_mirror(debug=debug)
            if debug:
                print(f"Segments Found: {len(self.segments)}")
                print(f"Actual Rotation Angle: {angle} degrees")
            if plot:
                for i in range(36):
                    orig_coords = self.rotate_coordinates(i,angle)
                    orig_coords[0] = int(orig_coords[0] * radius) + xCenter
                    orig_coords[1] = int(orig_coords[1] * radius) + yCenter
                    orig_point = orig_coords
                    x = orig_point[0]
                    y = orig_point[1]
                    r = matplotlib.patches.Rectangle((x-self.seg_size,y-self.seg_size),self.seg_size*2,self.seg_size*2, fill = False, color = "red")
                    ax1.add_patch(r)

                for i in range(len(self.segments)): 
                    point = self.segments[i]
                    ax2.scatter(point.x,point.y, marker = 'o', s = 20, color = "red", facecolors='none')

                ax2.scatter(self.center[0], self.center[1], c='green', marker='x', s=30)
                ax1.set_title(f"Original Mirror Segment Guess")
                i = ax1.imshow(self.file.data, cmap = 'gray')
                ax1.scatter(int(len(self.file.data[0]) / 2), int(len(self.file.data) / 2), s=30, c='green')
                fig.colorbar(i)
                im = ax2.imshow(self.file.data, cmap = 'gray')
                fig.colorbar(im)
                ax2.set_title(f'{self.file.instrument}: {self.file.file[12:20]} | Best Mirror Fit')
                plt.show()

    
    def find_new_center(self, data):
        """
        Finds the new center coordinates based on the maximum sums of the given data.
        Parameters:
        - data (numpy.ndarray): Input data array.
        Returns:
        - tuple: Tuple containing the new center coordinates.
        """
        if len(data) == 0:
            return 0,0
        data_t = data.transpose() # Transpose the data to work with columns as rows and vice versa
        sumx, sumy = [], []
        
        for r in range(min(len(data), len(data_t))): # Calculate the sums along columns and rows
            sumx.append(sum(data_t[r]))
            sumy.append(sum(data[r]))
        
        midx, midy = int(len(data) / 2), int(len(data[0]) / 2) # Calculate midpoints for default values
        
        max1X = max(sumx, default=midx) # Find the maximum sums along columns and rows
        max1Y = max(sumy, default=midy)
        try:
            xAdj = sumx.index(max(sumx)) # Find the indices of the maximum sums
            yAdj = sumy.index(max(sumy))
        except: 
            return midx, midy
        return xAdj, yAdj


    def goodFit(self, points):
        copy = self.segments
        self.segments = points
        self.clean_mirror()
        fit, range_new = (len(self.segments) == 36), len(self.segments)
        self.segments = copy
        return fit, range_new

    #-----------------------------------------------------#
    #------Get the last couple segments--if missing-------#
    #-----------------------------------------------------#
    def construct_new_mirror(self, debug=False):
        if len(self.segments) > 2 and len(self.segments) != 36 and self.file.data is not None: 
            # self.find_center_with_points()    
            if debug:
                print(f"Starting Find Center from current {len(self.segments)} points")
            # Takes a list of mirror points found and returns the center of those points
            radii, distances = [], []
            for point1 in self.segments:
                for point2 in self.segments:
                    dist = self.dist_between(point1, point2)
                    radii.append(dist)
                    distances.append(round(dist,0))
            distances.sort()
            try:
                distances.remove(-1)
            except:
                pass
            counter = dict(collections.Counter(distances))
            sorted(counter.keys())
            counter[-1] = 0
            min_scale = next((key for key, value in counter.items() if value >= 10), self.scale)
            radii = [r for r in radii if r < min_scale * 1.25]
            radii = [r for r in radii if r > min_scale]
            scale = np.mean(radii)
            if debug:
                print(f"Scale of the mirror: {scale}")
            # Now for each point, visit its neighbors and search for a point
            new_points = []
            for segment in self.segments:
                new_points.append(self.move_point(segment, self.file.rot_angle+30, scale))
                new_points.append(self.move_point(segment, self.file.rot_angle+90, scale))
                new_points.append(self.move_point(segment, self.file.rot_angle-30, scale))
                new_points.append(self.move_point(segment, self.file.rot_angle+150, scale))
                new_points.append(self.move_point(segment, self.file.rot_angle-90, scale))
                new_points.append(self.move_point(segment, self.file.rot_angle+210, scale))
            for point in new_points:
                self.segments.append(point)
            if debug:
                self.plot_mirror()
            self.recenter(box = self.seg_size * 1.5)
            self.recenter(box = self.seg_size * 1.2)
            self.recenter(box = self.seg_size * 1)
            self.recenter(box = self.seg_size * 0.8)
            self.remove_duplicate_segments()
            if debug:
                self.plot_mirror()
            self.remove_less_than_3_neighbors(scale+5)
            self.remove_less_than_3_neighbors(scale+5)
            self.clean_mirror()
            
            if debug:
                self.plot_mirror()
                print(f"End with {len(self.segments)} segments found")
                for seg in self.segments:
                    print(seg.x, seg.y)
        
        
    def remove_duplicate_segments(self):
        objects_dict = {}
        for seg in self.segments:
            x = seg.x
            y = seg.y
            key = (x, y)
            if key in objects_dict:
                objects_dict[key].append(seg)
            else:
                objects_dict[key] = [seg]
        duplicates = [objs[0] for objs in objects_dict.values() if len(objs) > 1]
        self.segments = duplicates

        
    def move_point(self, segment, angle_degrees, distance):
        try:
            angle_radians = math.radians(angle_degrees)
            new_x = segment.x + distance * math.cos(angle_radians)
            new_y = segment.y + distance * math.sin(angle_radians)
            s = Segment(int(new_x), int(new_y))
        except:
            print(segment.x, segment.y, angle_degrees, distance)
        return s
            
        
    def close_angles(self, angle1, angle2, threshold=3): 
        return abs(angle1 - angle2) <= threshold        

    
    def calculate_center(self):
        if len(self.segments) == 36:
            x,y = [],[]
            for segment in self.segments:
                x.append(segment.x)
                y.append(segment.y)
            self.center = int(sum(x) / len(x)), int(sum(y) / len(y))
    
    
    def recenter(self, box):
        box = int(box)
        for i in range(len(self.segments)):
            seg = self.segments[i]
            mirror = self.file.data[seg.y - box:seg.y + box,seg.x - box:seg.x + box]
            newx, newy = self.find_new_center(mirror)
            s = Segment(newx + seg.x - self.seg_size, newy + seg.y - self.seg_size)
            self.segments[i] = s
        
        
    def remove_less_than_3_neighbors(self, scale):
        new_points = []
        counts = []
        for point1 in self.segments:
            count = 0
            for point2 in self.segments:
                if point1 != point2 and self.dist_between(point1, point2) <= scale:
                    count += 1
            if count >= 3:
                new_points.append(point1)
            counts.append(count)
        # print(counts)
        self.segments = new_points


    def calc_rotation_angle(self, point1, point2):
        delta_x = point2.x - point1.x
        delta_y = point2.y - point1.y
        angle_rad = math.atan2(delta_y, delta_x)
        angle_deg = math.degrees(angle_rad)
        return angle_rad, angle_deg


    def dist_between(self, point1, point2):
        delta_x = point2.x - point1.x
        delta_y = point2.y - point1.y
        if delta_x + delta_y == 0:
            return -1
        return math.sqrt(delta_x**2 + delta_y**2)