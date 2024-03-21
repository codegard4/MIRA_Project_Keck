from astropy.io import fits
import numpy as np
class FitsFile:
    '''
    A Fits File stores the following attributes of a MIRA fits file
    - file: the string path to the file
    - data: the image data from the file
    - flip: whether the image is flipped (1) or not (0)
    - rot_posn: the rotator position from the file header
    - rot_angle: the rotation angle of the mirror in the image
    - pmfm: the PMFM of the image
    - instrument: the instrument used in the MIRA image
    - instrument_data: additional information about each of the instruments used on Keck (as of 1/1/24)
    '''
    
    def __init__(self, file):
        self.initialize_instrument_data()
        self.file = file
        
        self.pull_header()
        self.pull_data()
    
        
    def pull_data(self):    
        # Open the FITS file
        with fits.open(self.file) as hdul:
            # Get the data from the first image extension 
            data = hdul[0].data
        #Close the FITS file
        hdul.close()    
        if self.flip == 1:
            data = data[::-1]
            data = data[:,::-1] #flip around the x and y axes
        self.data = data
        self.shape = np.array(data).shape
        self.scaledArray()
        
        
    def scaledArray(self):
        """
        Scales the input array to the range [0, 1].
        Parameters:
        - data (numpy.ndarray): Input array to be scaled.
        Returns:
        - numpy.ndarray: Scaled array with values in the range [0, 1].
        """
        # Find the minimum and maximum values in the input array
        min_val = np.min(self.data)
        max_val = np.max(self.data)
        # Scale the input array to the range [0, 1]
        self.data = (self.data - min_val) / (max_val - min_val)
        
        
    def initialize_instrument_data(self):
        self.instrument_data = {
            'UNKNOWN': {'Name': 'UNKNOWN', 'RotAngle': 0.0, 'Flip': 0, 'PlateScale': 0.153, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'HIRES': {'Name': 'HIRES', 'RotAngle': 0.0, 'Flip': 0, 'PlateScale': 0.086, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'LRISBLUE': {'Name': 'LRISB', 'RotAngle': 0.0, 'Flip': 0, 'PlateScale': 0.135, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'LRISpBLUE': {'Name': 'LRISB', 'RotAngle': 0.0, 'Flip': 0, 'PlateScale': 0.135, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'LRISRED': {'Name': 'LRISRED', 'RotAngle': 90.0, 'Flip': 1, 'PlateScale': 0.209, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'LRIS': {'Name': 'LRIS', 'RotAngle': 90.0, 'Flip': 1, 'PlateScale': 0.209, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'LWS': {'Name': 'LWS', 'RotAngle': 64.0, 'Flip': 0, 'PlateScale': 0.08, 'Secondary': 25, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'ACAM': {'Name': 'ACAM', 'RotAngle': 0.0, 'Flip': 1, 'PlateScale': 0.1325, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'KCWI': {'Name': 'KCWI', 'RotAngle': -90.61, 'Flip': 0, 'PlateScale': 0.03037, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'NIRES': {'Name': 'NIRES', 'RotAngle': 180.0, 'Flip': 1, 'PlateScale': 0.12, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'NIRSPEC': {'Name': 'NIRSPEC', 'RotAngle': 275.0, 'Flip': 0, 'PlateScale': 0.155, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'NIRC2': {'Name': 'NIRC2', 'RotAngle': 0.0, 'Flip': 1, 'PlateScale': 1.0, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'NIRC': {'Name': 'NIRC', 'RotAngle': 180.0, 'Flip': 1, 'PlateScale': 0.15, 'Secondary': 25, 'StackCutoff': 0.15, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'ESI': {'Name': 'ESI', 'RotAngle': 97.2, 'Flip': 0, 'PlateScale': 0.138, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'DEIMOS': {'Name': 'DEIMOS', 'RotAngle': 179.0, 'Flip': 1, 'PlateScale': 0.103, 'Secondary': 15, 'StackCutoff': 1.0, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'MOSFIRE': {'Name': 'MOSFIRE', 'RotAngle': -0.244, 'Flip': 1, 'PlateScale': 0.18, 'Secondary': 15, 'StackCutoff': 0.1, 'ComaCutoff': 0.01, 'FocusCutoff': 0.01},
            'SSC1': {'Name': 'SSC1', 'RotAngle': 218.5, 'Flip': 1, 'PlateScale': 0.105, 'Secondary': 15, 'StackCutoff': 0.15, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'SSC2': {'Name': 'SSC2', 'RotAngle': 265.0, 'Flip': 1, 'PlateScale': 0.105, 'Secondary': 15, 'StackCutoff': 0.0, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'SSC': {'Name': 'SSC', 'RotAngle': 218.5, 'Flip': 1, 'PlateScale': 0.105, 'Secondary': 15, 'StackCutoff': 0.15, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15},
            'KPF': {'Name': 'KPF', 'RotAngle': -90.0, 'Flip': 1, 'PlateScale': 0.056, 'Secondary': 15, 'StackCutoff': 0.15, 'ComaCutoff': 0.15, 'FocusCutoff': 0.15}
        }
    # KCWI Plate Scale Original--0.03637
    # DEIMOS Plate Scale Original--0.119
    # ESI Plate Scale Original--0.153
    def pull_header(self):
        #If image does not contain the name of the instrument then get the instrument from the camera name
        camname_to_instrument = {'camname' : 'Instrument','mosfire' : 'MOSFIRE','hiresSlit' : 'HIRES','lrisSlit' : 'LRIS','kpf' : 'KPF','lrisslit': 'LRIS', 'ESI: Echellette Spectrograph and Imager': 'ESI', 'DEIMOS: real science mosaic CCD subsystem with PowerPC in VME crate': 'DEIMOS'}
        rot_angle = 0
        # Open the FITS file
        try:
            hdulist = fits.open(self.file)
        except: 
            print(f"This is an invalid file name")
        # Get the primary header
        header = hdulist[0].header
        # Attempt to extract instrument information
        try: # Try to get instrument information from 'INSTRUME' parameter
            self.instrument = header['INSTRUME']
        except KeyError:
            try: # Try to get instrument information from 'CAMNAME' parameter   
                camname = header['CAMNAME']
                try:
                    self.instrument = camname_to_instrument[camname]
                except KeyError: # If instrument is still not found, set to 'UNKNOWN'
                    self.instrument = 'UNKNOWN'
            except KeyError: #Cannot get file name from the header or the camname parameter
                try: # Try to get instrument information from file name
                    file_split = self.file.split("_")[0]
                    self.instrument = camname_to_instrument[file_split]
                except KeyError:
                    print("No way to determine instrument")
                    self.instrument = 'UNKNOWN'
        try: # Try to turn the instrument name into a stored instrument
            self.instrument = camname_to_instrument[self.instrument]
        except: 
            pass
        try:
            rot_angle = self.instrument_data[self.instrument]['RotAngle']
        except:
            pass
        try:
            self.flip = self.instrument_data[self.instrument]['Flip']
        except: 
            self.flip = 0
        try:
            self.pmfm = header['PMFM']
        except KeyError:
            self.pmfm = 300 
        try: 
            self.rot_posn = header['ROTPPOSN']
        except KeyError: 
            self.rot_posn = 0
        try:
            self.plate_scale = self.instrument_data[self.instrument]['PlateScale']
        except:
            self.plate_scale = 0.153
        # Close the FITS file
        hdulist.close()
        self.rot_angle = (rot_angle+self.rot_posn) 
        
        
    def set_rot_angle(self, angle):
        self.rot_angle = angle
        
        
    def __str__(self):
        s = f"FitsFile: {self.file}" + '\n' + f"{self.instrument} | PMFM: {self.pmfm} | Rotation Angle: {self.rot_angle} | Plate Scale: {self.plate_scale}"
        return s
        
        
if __name__ == '__main__':
    path = ""
    fileName = ""
    f = FitsFile(f'{path}{fileName}')
    print(f)