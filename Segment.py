class Segment:
    '''
    A class representing a mirror segment.

    Attributes:
    - number (int): Segment number
    - x (float): x-coordinate of the center of the segment.
    - y (float): y-coordinate of the center of the segment.
    - brightness (float): Brightness value of the center of the segment, calculated based on the image the segment is on.
    - diameter (float): Diameter of the segment.

    Methods:
    - set_brightness(value): Set the brightness value of the segment.
    - set_diameter(diam): Set the diameter of the segment.
    - set_number(num): Set the segment number.
    - calculate_diameter(data, threshold=0.5): Calculate the diameter of the segment based on the provided image data and a threshold.

    Magic Methods:
    - __lt__(Segment): Less than comparison based on x-coordinate.
    - __gt__(Segment): Greater than comparison based on x-coordinate.
    - __str__(): String representation of the segment.
    '''
    
    def __init__(self, x, y):
        self.number = None #Segment number on the mirror
        self.x = x # x coord of the center of the segment
        self.y = y # y coord of the center of the segment
        self.brightness = 0 # brightness value of the center of the segment that we can calculate given the image the segment is on
        self.diameter = None # diameter of the segment
        
    def set_brightness(value):
        self.brightness = value
        
        
    def set_diameter(diam):
        self.diameter = diam
        
        
    def set_number(num):
        self.nummber = num
        
        
    def __lt__(self, Segment):
        return self.x < Segment.x
    
    
    def __gt__(self, Segment):
        return self.x > Segment.x
        
        
    def calculate_diameter(self, data, threshold=0.5): #Very rough but simple way to calculate the diameter of a segment
        peak_value = self.brightness
        width = 0
        x_low = self.x
        x_high = self.x
        y_low = self.y
        y_high = self.y
        while data[x_low, self.y] >= (peak_value * threshold) and x_low >= 0:
            x_low -= 1
            width += 1
        while data[x_high, self.y] >= (peak_value * threshold) and x_high < data.shape[0] - 1:
            x_high += 1
            width += 1
        while data[self.x, y_low] >= (peak_value * threshold) and y_low >= 0:
            y_low -= 1
            width += 1
        while data[self.x, y_high] >= (peak_value * threshold) and y_high < data.shape[1] - 1:
            y_high += 1
            width += 1
        self.diameter = (width-4) / 2
        
        
    def __str__(self):
        s = f"Segment: ({self.x}, {self.y})" + '\n' + f"Brightness: {self.brightness} | Diameter: {self.diameter} | Number: {self.number}"
        return s
        
        
    if __name__ == '__main__':
        s = Segment(1,1)
        s.set_brightness(10)
        s.set_diameter(3.0)
        s.set_number(1)
        print(s)