import matplotlib.pyplot as plt
import math
import sys
import argparse
import glob
from Fits import FitsFile
from Segment import Segment
from Mirror import Mirror

class MiraImage:
    def __init__(self, file):
        print(file)
        self.file = FitsFile(file)
        self.mirror = Mirror(file) #stores angle and center

        
    def find_angle(self, centroid, point):
        #find the rotation angle of a point around a center
        dx = point[0] - centroid[0]
        dy = point[1] - centroid[1]
        angle_rad = math.atan2(dy, dx)
        angle_deg = math.degrees(angle_rad)
        return angle_deg if angle_deg >= 0 else angle_deg + 360


    def create_mirror(self, center, angle, radius):
        #create the coordinates of what the mirror should look like with a given center, angle and radius
        coordinates = []
        for i in range(self.mirror.n_segments):
            coords = self.mirror.rotate_coordinates(i, angle)
            coords[0] = int((coords[0] * radius) + center[0])
            coords[1] = int((coords[1] * radius) + center[1])
            coordinates.append((coords[0],coords[1]))
        return coordinates

    
    def evaluate_segments(self, plot=True):
        self.mirror.calculate_center()
        if len(self.mirror.segments) == 36:
            correct_mirror = self.create_mirror(self.mirror.center, math.radians(self.file.rot_angle), radius = self.mirror.scale)
            if plot:
                fig, (ax, ax2) = plt.subplots(1, 2, figsize=(20, 10))
                for i in range(self.mirror.n_segments):
                    ax.imshow(self.file.data)
                    ax.set_title("Actual Centroids (red) vs Expected Centroids (green)")
                    ax2.imshow(self.file.data)
                    ax2.scatter(self.mirror.segments[i].x, self.mirror.segments[i].y, color = 'green', s = 30, facecolors = 'None')
                    # ax2.annotate('%s' % (i+1), xy=(self.mirror.segments[i].x,self.mirror.segments[i].y), textcoords='data', color = 'white')
                    ax.scatter(correct_mirror[i][0], correct_mirror[i][1], color = 'green', s = 20, facecolors ='None')
                    ax.scatter(self.mirror.center[0], self.mirror.center[1], color = 'green', s = 30, marker = '+')
                    ax2.scatter(self.mirror.center[0], self.mirror.center[1], color = 'green', s = 30, marker = '+')
                    ax.scatter(self.mirror.segments[i].x, self.mirror.segments[i].y, color = 'red', s = 5, marker = 'x')
                    # ax.annotate('%s' % (i+1), xy=(self.mirror.segments[i].x,self.mirror.segments[i].y), textcoords='data', color = 'red')
                    # ax.annotate('%s' % (i+1), xy=correct_mirror[i], textcoords='data', color = 'green')
            #Zernike calculation
            segs = []
            for seg in self.mirror.segments:
                segs.append((seg.x, seg.y))
            segs = [Segment(u, v) for u, v in segs]
            cen = [Segment(u, v) for u, v in correct_mirror]
            return (self.zernikes(segs, cen))
        return "Not All Segments Found"
    

    def Su(self, cen):
        return sum(point.x for point in cen)
    def Sv(self, cen):
        return sum(point.y for point in cen)
    def Sxx(self, segs):
        return sum(point.x * point.x for point in segs)
    def Syy(self, segs):
        return sum(point.y * point.y for point in segs)
    def Sux(self, segs, cen):
        return sum(point.x * seg.x for point, seg in zip(cen, segs))
    def Suy(self, segs, cen):
        return sum(point.x * seg.y for point, seg in zip(cen, segs))
    def Svx(self, segs, cen):
        return sum(point.y * seg.x for point, seg in zip(cen, segs))
    def Svy(self, segs, cen):
        return sum(point.y * seg.y for point, seg in zip(cen, segs))
    def Suxx(self, segs, cen):
        return sum(point.x * seg.x * seg.x for point, seg in zip(cen, segs))
    def Suyy(self, segs, cen):
        return sum(point.x * seg.y * seg.y for point, seg in zip(cen, segs))
    def Svxx(self, segs, cen):
        return sum(point.y * seg.x * seg.x for point, seg in zip(cen, segs))
    def Svyy(self, segs, cen):
        return sum(point.y * seg.y * seg.y for point, seg in zip(cen, segs))
    def Suxy(self, segs, cen):
        return sum(point.x * seg.x * seg.y for point, seg in zip(cen, segs))
    def Svxy(self, segs, cen):
        return sum(point.y * seg.x * seg.y for point, seg in zip(cen, segs))
    def Sxxxx(self, segs):
        return sum(point.x**4 for point in segs)
    

    def zernikes(self, segs, cen):
        #Segment center coords expected--
        #image centroids
        RADIUS = 5.47 #Radius of the primary mirror in meters

        ARCSEC_PER_RAD = 206264.806241  
        # Powers of the RADIUS
        r1 = RADIUS
        r2 = RADIUS * RADIUS
        r4 = RADIUS * RADIUS * RADIUS * RADIUS

        # RADIUS/2 and a radian-conversion factor
        r1O2s = r1 / (2. * ARCSEC_PER_RAD)

        # The "a"s: from the mirror only
        a11 = 36.
        a22 = 36.
        a33 = 8. * self.Sxx(segs) / r2
        a44 = 4. * self.Sxx(segs) / r2
        a55 = 4. * self.Sxx(segs) / r2
        a66 = 20. * self.Sxxxx(segs) / r4 - 8. * self.Sxx(segs) / r2 + 24.
        a77 = 20. * self.Sxxxx(segs) / r4 - 8. * self.Sxx(segs) / r2 + 24.
        a61 = 2. * self.Sxx(segs) / r2 - 12.
        a72 = 2. * self.Sxx(segs) / r2 - 12.
        a16 = 6. * a61
        a27 = 6. * a72

        # The "B"s calculated from mirror and image data
        B1 = r1O2s * self.Sv(cen)
        B2 = r1O2s * self.Su(cen)

        B3 = r1O2s * (self.Sux(segs, cen) / r1 + self.Svy(segs, cen) / r1)
        B4 = r1O2s * (self.Suy(segs, cen) / r1 + self.Svx(segs, cen) / r1)
        B5 = r1O2s * (self.Sux(segs, cen) / r1 - self.Svy(segs, cen) / r1)

        B6 = r1O2s * ((self.Suxy(segs, cen) / r2) +
                      (self.Svxx(segs, cen) / (2. * r2)) +
                      (3. * self.Svyy(segs, cen) / (2. * r2)) -
                      (self.Sv(cen) / 3.))

        B7 = r1O2s * (3. * self.Suxx(segs, cen) / (2. * r2) +
                      self.Suyy(segs, cen) / (2. * r2) -
                      self.Su(cen) / 3. +
                      self.Svxy(segs, cen) / r2)

        # The Zernikes
        zc1m1 = (a66 * B1 - a16 * B6) / (a11 * a66 - a61 * a16) #Defocus--Light rays coming from a single point do not converge at the focal point, which creates a blurry image
        zc1p1 = (a77 * B2 - a27 * B7) / (a22 * a77 - a72 * a27) #Defocus
        zc20 = B3 / a33  #Astigmatism--Light rays do not converge due to curvature of the optical surface, which creates a stretched or distorted image
        zc2m1 = B4 / a44 #Astigmatism
        zc2p1 = B5 / a55 #Astigmatism
        zc3m1 = (a11 * B6 - a61 * B1) / (a11 * a66 - a61 * a16) #Coma--Affects clarity of stars in the image. more pronounced towards the edge of the image
        zc3p1 = (a22 * B7 - a72 * B2) / (a22 * a77 - a72 * a27) #Coma

        return zc1m1, zc1p1, zc20, zc2m1, zc2p1, zc3m1, zc3p1

    
def parseArguments(in_args):
    description = "Acquires MIRA Mirror Segments"
    usage = "\n{} python MiraImage.py [-f] filename OR [-d] date [-v] and [-p] to see results \n".format(in_args[0])
    epilog = ""
    parser = argparse.ArgumentParser(description = description, usage = usage, epilog = epilog)
    parser.add_argument("-f", "--fileName", dest = "file", type = str, help = "Path to the MIRA file to be analyzed", default = None)
    parser.add_argument("-d", "--date", dest = "date", type = str, help = "Date to be analyzed (format: day,month,year)", default = None)
    parser.add_argument("-v", "--verbose", dest = "debug", type = bool, help = "Debugger output? (True/False)", default = False)
    parser.add_argument("-p", "--plot", dest = "plot", type = bool, help = "Plot Segments? (True/False)", default = False)
    args = None
    try:
        args = parser.parse_args(in_args[1:])
    except Exception as e:
        print(e)
        parser.print_help()
        sys.exit(0) 
    return args 


def getFiles(d, m, y):
    day, month, year = "{:02d}".format(d), "{:02d}".format(m), "{:02d}".format(y)
    path = ''
    pattern = '[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].fits'
    fits_files = glob.glob(f'{path}/{pattern}')
    fits_files.sort()
    return fits_files


if __name__ == "__main__":
    files = []
    args = parseArguments(sys.argv)
    if args.date is not None:
        d,m,y = args.date.split(",")
        files = getFiles(int(d),int(m),int(y))
    if args.file is not None:
        files.append(args.file)
    print(files)
    for file in files:
        m = MiraImage(file)
        m.mirror.find_center(plot = args.plot)
        m.mirror.find_mirror_segments(debug = args.debug)
        m.mirror.construct_new_mirror()
        m.mirror.find_mirror_segments(attempt = 2, debug = args.debug)
        z = m.evaluate_segments()
        print(f"File {file} analyzed")
        print(f"Center: {m.mirror.center} Rotation Angle: {m.mirror.rot_angle}")
        print(f"{len(m.mirror.segments)} segments found")
        print(f"Zernikes: {z}")
        