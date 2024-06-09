from PIL import Image
import math
import argparse
# ...
class myClass:
    def __init__(self, fileName):
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-gray.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-smallgap.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-smoothcolor.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-checkers.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-depth.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-elements.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-sRGB.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-gammabox.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-perspective.txt'
        ##self.file_path = '/Users/ziwei/Documents/UIUC-MCS/CS418/rasterizer-files/rast-frustum.txt'
        self.file_path = fileName
        self.png = []
        self.position4 = []
        self.color3 = []
        self.drawArraysTriangles = []
        self.image = Image.new("RGBA", (0, 0), (0,0,0,0))
        self.depth = False
        self.pixel_pos = {}
        self.elements = []
        self.sRGB = False
        self.hyp = False
    def process_input(self):
        with open(self.file_path, encoding = "unicode-escape") as f1:
            for line in f1:
                if (len(line.strip()) == 0):
                    continue
                words = line.strip().split()
                keyword = words[0]
                if keyword == 'png':
                    width = int(words[1])
                    height = int(words[2])
                    filename = words[3]
                    self.png = [width, height, filename]
                    self.image = Image.new("RGBA", (width, height), (0,0,0,0))
                elif keyword == 'sRGB':
                    self.sRGB = True
                elif keyword =='hyp':
                    self.hyp = True
                elif keyword == 'position' and words[1] == '4':
                    self.position4 = [(float(words[i]), float(words[i+1]),  float(words[i+2]),  float(words[i+3])) for i in range(2, len(words), 4)]
                elif keyword == 'color' and words[1] == '3':
                    color = []
                    for i in range(2, len(words), 3):
                        if (self.sRGB):
                            r = float(words[i])
                            g = float(words[i+1])
                            b = float(words[i+2])
                            a = 1
                            color.append([r,g,b,a])
                        else:
                            r = float(words[i]) * 255
                            g = float(words[i+1]) * 255
                            b = float(words[i+2]) * 255
                            a = 255
                            color.append([r,g,b,a])
                    self.color3 = color
                elif keyword == 'drawArraysTriangles':
                    self.drawArraysTriangle(int(words[1]), int(words[2]))
                elif keyword == 'depth':
                    self.depth = True
                elif keyword == 'drawElementsTriangles':
                    self.drawElementsTriangle(int(words[1]), int(words[2]))
                elif keyword == 'elements':
                    for i in range(1, len(words), 3):
                        self.elements.append([int(words[i]), int(words[i+1]), int(words[i+2])])
                    ##print(len(self.elements))
                    ##print(len(self.position4))
                    ##self.element
                    
        self.image.save(self.png[2])


    def drawElementsTriangle(self, count, offset):
        view_port_position4 = [
            [(x / w + 1) * self.png[0] / 2, (y / w + 1) *  self.png[1] / 2, z, w] for x, y, z, w in self.position4
        ]
        for c in range(int(count/3)):
            
            element = self.elements[c]
            pos1 = view_port_position4[element[offset + 0]]
            pos2 = view_port_position4[element[offset + 1]]
            pos3 = view_port_position4[element[offset + 2]]
            pos1_index = view_port_position4.index(pos1)
            pos2_index = view_port_position4.index(pos2)
            pos3_index = view_port_position4.index(pos3)
            color1 = self.color3[pos1_index]
            color2 = self.color3[pos2_index]
            color3 = self.color3[pos3_index]
            if self.depth == False:
                self.pixel_pos = {}
            ##print("pos", self.pixel_pos)
            self.scanline_algorithm(pos1 + color1, pos2 + color2, pos3 + color3)
            self.draw_image()
        
    def drawArraysTriangle(self, first, total_pos):
        view_port_position4 = [[(x / w + 1) * self.png[0] / 2, (y / w + 1) *  self.png[1] / 2, z, w] for x, y, z, w in self.position4]
        ##print("1", self.position4, self.color3)
        if (self.hyp == True):
            view_port_position4 = [
                [x, y, z/w, w] for x, y, z, w in view_port_position4
            ]
            ##print("2", view_port_position4)
        count = int(total_pos/3)
        for c in range(count):
            pos1 = view_port_position4[first + c * 3 + 0]
            pos2 = view_port_position4[first + c * 3 + 1]
            pos3 = view_port_position4[first + c * 3 + 2]
            color1 = self.color3[first + c * 3 + 0]
            color2 = self.color3[first + c * 3 + 1]
            color3 = self.color3[first + c * 3 + 2]
            ##print("initial", pos1 + color1, pos2 + color2, pos3 + color3)
            if (self.hyp == True):
                color1 = [col/pos1[3] for col in color1]
                color2 = [col/pos2[3] for col in color2]
                color3 = [col/pos3[3] for col in color3]

                pos1[3] = 1/pos1[3] ##1/w
                pos2[3] = 1/pos2[3] ##1/w
                pos3[3] = 1/pos3[3] ##1/w
                ##print(color1, color2, color3)
            if self.depth == False:
                self.pixel_pos = {}
                
            ##print("initial", pos1 + color1, pos2 + color2, pos3 + color3)
            self.scanline_algorithm(pos1 + color1, pos2 + color2, pos3 + color3)
            ##print(self.pixel_pos)
            self.draw_image()
        
    def draw_image(self):
        if self.depth:
            draw_pixels = []
            for key, value in self.pixel_pos.items():
                if (self.hyp):
                    sorted_value = sorted(value, key=lambda x: x[2])
                    pixel = sorted_value[0]
                    draw_pixels.append(pixel)
                    ##print("3", pixel)
                    
                    for i in range(4, len(pixel)):
                       pixel[i] = pixel[i]/pixel[3] ##rgb/w'
                    pixel[3] = 1/pixel[3] ##1/w'

                    ##print("4", pixel)
                    if (self.sRGB):
                        for i in range(4, len(pixel) - 1):
                            if pixel[i] <= 0.0031308: pixel[i] = pixel[i] * 12.92
                            else: pixel[i] = 1.055 * ((pixel[i]) ** (1/2.4)) - 0.055
                        ##print(pixel)
##                    for i in range(4, len(pixel) - 1):
##                       pixel[i] = pixel[i]/pixel[3] ##rgb/w'
##                    pixel[3] = 1/pixel[3] ##1/w'
                    self.image.im.putpixel((int(pixel[0]), int(pixel[1])), (int(pixel[4] * 255), int(pixel[5] * 255), int(pixel[6] * 255), int(pixel[7] * 255)))
                else:
                    sorted_value = sorted(value, key=lambda x: x[2])
                    ##print(sorted_value)
                    
                    pixel = sorted_value[0]     
                    self.image.im.putpixel((int(pixel[0]), int(pixel[1])), (int(pixel[4]), int(pixel[5]), int(pixel[6]), int(pixel[7])))
                
        else:
            ##print(self.pixel_pos)
            for key, value in self.pixel_pos.items():
                pixel = value[0]
                if (pixel[0] >= self.png[0] or pixel[1] >= self.png[1]): continue
                if (self.sRGB):
                    for i in range(4, len(pixel) - 1):
                        if pixel[i] <= 0.0031308: pixel[i] = pixel[i] * 12.92
                        else: pixel[i] = 1.055 * ((pixel[i]) ** (1/2.4)) - 0.055
                ##if (self.hyp):
                    self.image.im.putpixel((int(pixel[0]), int(pixel[1])), (int(pixel[4] * 255), int(pixel[5] * 255), int(pixel[6] * 255), int(pixel[7] * 255)))
                ##print(pixel)
                else:
                    self.image.im.putpixel((int(pixel[0]), int(pixel[1])), (int(pixel[4]), int(pixel[5]), int(pixel[6]), int(pixel[7])))
    
    def dda_algorithm(self, a, b, is_x_axis):
        e = 0
        s = []
        p = []
        d = []

        if is_x_axis:
            if (a[0][0] == b[0][0]): return
            if (a[0][0] > b[0][0]): a, b = b, a
            d = (b[0][0] - a[0][0])
            e = math.ceil(a[0][0]) - a[0][0]
        else:      
            if (a[0][1] == b[0][1]): return 
            if (a[0][1] > b[0][1]): a, b = b, a
            d = (b[0][1] - a[0][1])
            e = math.ceil(a[0][1]) - a[0][1]
         
        delta = [element1 - element2 for element1, element2 in  zip(b[1], a[1])]
        s = [element/d for element in delta]
        o = [e * element for element in s]
        p = [p1 + p2 for p1, p2 in zip(a[1], o)]
        
        self.dda_find_all_points(p, b, s, is_x_axis)
        return (((p[:2]), p), ((s[:2]), s))
    
    def dda_find_all_points(self, p, b, s, is_x_axis):
        if is_x_axis:
            ##self.pixel_pos[tuple(p[:2])] = p
            ##print('lx points', p[0], b)
            while(p[0] < b[0][0]):
                if tuple(p[:2]) in self.pixel_pos:
                    ##print("duplicate")
                    self.pixel_pos[tuple(p[:2])].append(p)
                else:
                    self.pixel_pos[tuple(p[:2])] = [p]
                p = [p1 + s1 for p1, s1 in zip(p, s)]               
                ##self.pixel_pos[tuple(p[:2])] = p
            
    def sRGB_Conversion(self, color):
        srgb_color = []
        for element in color:
##            if element <= 0.0031308: element = element * 12.92
##            else: element = 1.055 * ((element) ** (1/2.4)) - 0.055
##            if element <= 0.04045: element = element/12.92
##            else: element = ((element + 0.055)/1.055) ** 2.4
            srgb_color.append(element)
        ##srgb_color.append(1)
        return srgb_color
    def scanline_algorithm(self, p, q, r):
        vectors = [p, q, r]    
        vertex = {}
        t = {}
        m = {}
        b = {}
        
        p_prime, s_prime, p_, s_ = {}, [], {}, []
        
        for vec in vectors:
            key = tuple(vec[:2])  
            if key in vertex:
                vertex[key] = vec
            else:
                vertex[key] = vec
        sorted_dict = dict(sorted(vertex.items(), key=lambda item: item[0][1]))
        converted_list = list(sorted_dict.items())
        ##print("tmb", converted_list)
        if (len(converted_list) < 3): return
        t = converted_list[0]
        m = converted_list[1]
        b = converted_list[2]
        
        tb = self.dda_algorithm(t, b, False)
        p_vec_prime = tb[0]
        s_vec_prime = tb[1]
        tm = self.dda_algorithm(t, m, False)
        if (tm is not None):
            p_vec = tm[0]
            s_vec = tm[1]

            while (p_vec[0][1] < m[0][1]):
                self.dda_algorithm(p_vec, p_vec_prime, True)
                p_vec_prime = ([p + s for p, s in zip(p_vec_prime[0], s_vec_prime[0])], [p + s for p, s in zip(p_vec_prime[1], s_vec_prime[1])])
                p_vec = ([p + s for p, s in zip(p_vec[0], s_vec[0])], [p + s for p, s in zip(p_vec[1], s_vec[1])])

            
        
        mb = self.dda_algorithm(m, b, False)
        if (mb is not None):
            p_vec = mb[0]
            s_vec = mb[1]

            while (p_vec[0][1] < b[0][1]):
                self.dda_algorithm(p_vec, p_vec_prime, True)
                p_vec_prime = ([p + s for p, s in zip(p_vec_prime[0], s_vec_prime[0])], [p + s for p, s in zip(p_vec_prime[1], s_vec_prime[1])])
                p_vec = ([p + s for p, s in zip(p_vec[0], s_vec[0])], [p + s for p, s in zip(p_vec[1], s_vec[1])])              
        
        
if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Read file.")
     parser.add_argument('filename', help='Name of the file')
     args = parser.parse_args()
     filename = args.filename
     draw = myClass(filename)
     draw.process_input()

##    
##draw = myClass()
##draw.process_input()
