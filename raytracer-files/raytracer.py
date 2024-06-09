import math
import random
from PIL import Image
import copy
import argparse

def subtract_vectors(vec1, vec2):
    return [v1 - v2 for v1, v2 in zip(vec1, vec2)]

def dot_product(vec1, vec2):
    return sum(v1 * v2 for v1, v2 in zip(vec1, vec2))

def length_squared(vec):
    return sum(v**2 for v in vec)

def length(vec):
    return math.sqrt(sum(v**2 for v in vec))

def normalize_vector(vec):
    length = math.sqrt(length_squared(vec))
    return [v / length for v in vec]

def add_vectors(vec1, vec2):
    return [v1 + v2 for v1, v2 in zip(vec1, vec2)]

def multiply_vector(vec, scalar):
    return [v * scalar for v in vec]

def divide_vector(vec, scalar):
    return [v / scalar for v in vec]

def multiply_vectors(vec1, vec2):
    return [v1 * v2 for v1, v2 in zip(vec1, vec2)]

def cross_product(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]

def randomize_value(value, range_size=0.01):
    return random.uniform(value - range_size, value + range_size)

def ray_sphere_intersection(ray_origin, ray_direction, sphere_center, sphere_radius):

    inside = length_squared(subtract_vectors(sphere_center, ray_origin)) < sphere_radius**2
    
    tc = dot_product(subtract_vectors(sphere_center, ray_origin), ray_direction) / length_squared(ray_direction)

    if not inside and tc < 0:
        return None
    
    d_squared = length_squared(subtract_vectors(add_vectors(ray_origin, multiply_vector(ray_direction, tc)), sphere_center))

    if not inside and sphere_radius**2 < d_squared:
        return None

    toffset = math.sqrt(sphere_radius**2 - d_squared) / length_squared(ray_direction)

    if inside:
        t = tc + toffset
        intersection_point = add_vectors(multiply_vector(ray_direction, t), ray_origin)
        return (t, intersection_point)
        
    else:
        t = tc - toffset
        intersection_point = add_vectors(multiply_vector(ray_direction, t), ray_origin)
        return (t, intersection_point)

def ray_plane_intersection(ray_origin, ray_direction, plane_normal, p):
    if dot_product(ray_direction, plane_normal) == 0:
        return None
    t = dot_product(subtract_vectors(p, ray_origin), plane_normal)/dot_product(ray_direction, plane_normal)
    if t > 0:
        intersection_point = add_vectors(multiply_vector(ray_direction, t), ray_origin)
        return (t, intersection_point)
    else:
        return None
def ray_triangle_intersection(ray_origin, ray_direction, p0, p1, p2, plane_normal):
    ##print(p0, p1, p2)

    a1 = cross_product(subtract_vectors(p2, p0), plane_normal)
    a2 = cross_product(subtract_vectors(p1, p0), plane_normal)
    A = plane_normal[0]
    B = plane_normal[1]
    C = plane_normal[2]
    D = -dot_product(plane_normal, p0);
    p = divide_vector(multiply_vector(plane_normal, D * -1.0), (A ** 2 + B ** 2 + C ** 2))
    ##print(p, p0, p1)
    ##t = - (dot_product(plane_normal, ray_origin) + D) / dot_product(plane_normal, ray_direction);
    ##Phit = add_vectors(ray_origin, multiply_vector(ray_direction, t));
    intersection = ray_plane_intersection(ray_origin, ray_direction, plane_normal, p0)
    
    if intersection:
        check1 = dot_product(cross_product(subtract_vectors(p1, p0), subtract_vectors(intersection[1], p0)), plane_normal)
        check2 = dot_product(cross_product(subtract_vectors(p2, p1), subtract_vectors(intersection[1], p1)), plane_normal)
        check3 = dot_product(cross_product(subtract_vectors(p0, p2), subtract_vectors(intersection[1], p2)), plane_normal)

        
        e1 = multiply_vector(a1, 1/dot_product(a1, subtract_vectors(p1, p0)))
        e2 = multiply_vector(a2, 1/dot_product(a2, subtract_vectors(p2, p0)))
       
        b1 = None
        b2 = None
        b0 = None


        ##print(dot_product(e2, subtract_vectors(p2, p0)))
        if int(dot_product(e2, subtract_vectors(p2, p0))) == 1:
            b2 = dot_product(e2, subtract_vectors(intersection[1], p0))
            
        if int(dot_product(e1, subtract_vectors(p1, p0))) == 1:
            b1 = dot_product(e1, subtract_vectors(intersection[1], p0))

        if b1 and b2:
            b0 = 1 - b1 - b2
        if check1 > 0 and check2 > 0 and check3 > 0:
            return intersection

##        
##        if b0 and b1 and b2:
##            ##print(p[1], add_vectors(add_vectors(multiply_vector(p0, b0),multiply_vector(p1, b1)), multiply_vector(p2, b2)))
##            if 0<=b0<=1 and 0<=b1<=1 and 0<=b2<=1:
##                ##print(p)
##                return intersection
    return None
    
class Raytracer:
    def __init__(self, filename):
        self.filename = filename
        self.width = 0
        self.height = 0
        self.output_file =""
        self.scene = []
        self.sun = []
        self.color = (1,1,1)
        self.expose = None
        self.e = [0, 0, 0]  
        self.f = [0, 0, -1]  
        self.u = [0, 1, 0]  
        self.r = [1, 0, 0]
        self.up = None
        self.fisheye = None
        self.panorama = None
        self.xyz = []
        self.tris = []
        self.aa = None
        self.rgbaInfo = {}
        self.dof = {}
    def parse_scene_file(self):
        with open(self.filename, 'r') as file:
            for line in file:
                ##print(line)
                if line.startswith("png"):
                    tokens = line.split()
                    self.width = int(tokens[1])
                    self.height = int(tokens[2])
                    self.output_file = tokens[3]
                if line.startswith("color"):
                    tokens = line.split()
                    r, g, b = map(float, tokens[1:])
                    self.color = (r,g,b)
                if line.startswith("sphere"):
                    tokens = line.split()
                    center_x, center_y, center_z, radius = map(float, tokens[1:])
                    sphere = {'type': 'sphere', 'center': (center_x, center_y, center_z), 'radius': radius, 'color': self.color}
                    self.scene.append(sphere)
                if line.startswith("sun"):
                    tokens = line.split()
                    center_x, center_y, center_z = map(float, tokens[1:])
                    suns = {'type': 'sun', 'direction': (center_x, center_y, center_z), 'color': self.color}
                    ##print(suns)
                    self.sun.append(suns)
                if line.startswith("expose"):
                    tokens = line.split()
                    self.expose = float(tokens[1])
                if line.startswith("eye"):
                    tokens = line.split()
                    self.e = [float(value) for value in tokens[1:]]
                if line.startswith("forward"):
                    tokens = line.split()
                    self.f =  [float(value) for value in tokens[1:]]
                    if self.up:
                        self.r = normalize_vector(cross_product(self.f, self.u))
                        self.u = normalize_vector(cross_product(self.r, self.f))          
                if line.startswith("up"):
                    tokens = line.split()
                    self.up =  [float(value) for value in tokens[1:]]
                    self.r = normalize_vector(cross_product(self.f, self.up))
                    self.u = normalize_vector(cross_product(self.r, self.f))     
                if line.startswith("fisheye"):
                    self.fisheye = True
                if line.startswith("panorama"):
                    self.panorama = True
                if line.startswith("plane"):
                    tokens = line.split()
                    a, b, c, d = map(float, tokens[1:])
                    ##print(a, b, c,d)
                    plane = {'type': 'plane', 'center': [a,b,c], 'D': d, 'color': self.color}
                    self.scene.append(plane)
                if line.startswith("bulb"):
                    tokens = line.split()
                    x, y, z = map(float, tokens[1:])
                    bulb = {'type': 'bulb', 'pos': (x, y, z), 'color': self.color}
                    ##print(bulb)
                    self.sun.append(bulb)
                if line.startswith("xyz"):
                    tokens = line.split()
                    x, y, z = map(float, tokens[1:])
                    self.xyz.append([x, y, z])
                if line.startswith("tri"):
                    tokens = line.split()
                    i1, i2, i3 = map(int, tokens[1:])
                    p0 = None
                    p1 = None
                    P2 = None
                    if i1 >= 0:
                        p0 = self.xyz[i1 - 1]
                    else:
                        p0 = self.xyz[i1]
                    if i2 >= 0:
                        p1 = self.xyz[i2 - 1]
                    else:
                        p1 = self.xyz[i2]

                    if i3 >= 0:
                        p2 = self.xyz[i3 - 1]
                    else:
                        p2 = self.xyz[i3]
                    
                    plane_normal = normalize_vector(cross_product(subtract_vectors(p1, p0), subtract_vectors(p2, p0)))
                    
                    tri = {'type': 'tri', 'pos': [p0, p1, p2], 'center': plane_normal, 'color': self.color}
                    ##print(tri)
                    self.scene.append(tri)
                if line.startswith("aa"):
                    tokens = line.split()
                    self.aa = float(tokens[1])
                if line.startswith("dof"):
                    tokens = line.split()
                    self.dof = {'focus': float(tokens[1]), 'lens': float(tokens[2])}
                    
    def calculate_light_direction(self, ray_origin, light_pos, intersection_point):
        return add_vectors(multiply_vector(light_pos, intersection_point[0]), ray_origin)
        
    def calculate_linear_color(self, intersection_color, light_color, lambert_dot_product):
        lambert = multiply_vector(multiply_vectors(intersection_color, light_color), lambert_dot_product)
        return lambert
    def calculate_sRGB(self, color):
        if color <= 0.0031308:
            return 12.92 * color
        else:
            return 1.055 * (color ** (1/2.4)) - 0.055
        
    def distance_to_ray_origin(self, ray_origin, sphere):
        return length(subtract_vectors(sphere['center'], ray_origin))

    
    def calculate_exposed_linear_color(self, linear_color, exposure):
        if exposure:
            exposed_linear_color = [1 - math.exp(-exposure * channel) for channel in linear_color]
        else:
            exposed_linear_color = linear_color
        return exposed_linear_color

    def check_shadow(self, intersection_point, sun_direction, sphere):
        shadow_intersection= self.ray_sphere_intersection(intersection_point, sun_direction, sphere)
        if shadow_intersection and shadow_intersection['t'] < 1.0:
            return True
        else:
            return False
        
    def calculate_reflection_direction (self, ray_origin, surface_normal):
        direction = subtract_vectors(multiply_vector(multiply_vectors(multiply_vectors(ray_origin, surface_normal), surface_normal),2), ray_origin)
        return direction 
                    
    def calculate_transparency(Self, ray_origin, surface_normal, n1, n2):
        a = dot_product(ray_origin, surface_normal)
        b = (n1/n2) * math.sqrt(1 - a)
        if b >= 1:
            dot2 = 2 * dot_product(surface_normal, ray_origin)
            return subtract_vectors(multiply_vector(multiply_vector(surface_normal, 2 * dot_product(surface_normal, ray_origin)), ray_origin))
        c = math.sqrt(1 - (n1**2/n2**2)*(1 - a))
        return subtract_vectors((n1/n2) * subtract_vectors(multiply_vector(ray_origin, -1), multiply_vector(surface_normal, a)), multiply_vector(surface_normal, c))

    def convert_coordinates(self, x, y, fixed_width, fixed_height):
        scale_x = 360.0 / fixed_width
        scale_y = 180 / fixed_height


        target_x = int(x * scale_x)
        target_y = int(y * scale_y)

        return target_x, target_y
    
    def determine_intersection(self, ray_origin, ray_direction, sphere):
        intersection_point = None
        if sphere['type'] == 'plane':
            A = sphere['center'][0]
            B = sphere['center'][1]
            C = sphere['center'][2]
            p = divide_vector(multiply_vector(sphere['center'], sphere['D'] * -1.0), (A ** 2 + B ** 2 + C ** 2))
            ##print (sphere['center'], p)
            intersection_point = ray_plane_intersection(ray_origin, ray_direction, sphere['center'], p)
        elif sphere['type'] == 'tri':
            intersection_point = ray_triangle_intersection(ray_origin, ray_direction, sphere['pos'][0], sphere['pos'][1], sphere['pos'][2], sphere['center'])
        else:
            intersection_point = ray_sphere_intersection(ray_origin, ray_direction, sphere['center'], sphere['radius'])
        return intersection_point
    
    def calculate_light_intensity(self, light_position, intersection_point, light_color, intersection_color, surface_normal):
        light_direction = normalize_vector(subtract_vectors(light_position, intersection_point))
        distance_to_light = length(subtract_vectors(intersection_point, light_position))
        light_intensity = 1.0/(distance_to_light ** 2)
        lambert_dot_product = dot_product(normalize_vector(surface_normal), normalize_vector(light_direction))
        linear_color = None
        ##linear_color = None
        if lambert_dot_product < 0:
            linear_color = (0,0,0)
        else:
            linear_color = [c * lambert_dot_product * light_intensity * i for i, c in zip(light_color, intersection_color)]
       
        return linear_color

    def is_light_between_spheres(self, light_position, shadow_sphere, sphere):
        # Check if the light is between two spheres
        distance_to_sphere = length(subtract_vectors(light_position, sphere['center']))
        distance_to_shadow_sphere = length(subtract_vectors(light_position, shadow_sphere['center']))
        is_inside_sphere1 = distance_to_sphere < sphere['radius']
        is_inside_sphere2 = distance_to_shadow_sphere < shadow_sphere['radius']

        return not is_inside_sphere1 and not is_inside_sphere2

    def draw_spheres(self):
        image = Image.new("RGBA", (self.width, self.height), color=(0,0,0,0))
        pixels = image.load()
        
        ##print(self.scene)
        ##print(self.sun)
        if self.aa is None:
            self.generate_lsrgb()
            for key, value in self.rgbaInfo.items():
                image.im.putpixel(key, value)
            image.save(self.output_file)
        else:
            count = 0
            while  count < self.aa:
                self.generate_lsrgb()
                count +=1
            ##print(self.rgbaInfo)
            for key, value in self.rgbaInfo.items():
                r_values, g_values, b_values = zip(*value)
                average_r = int(sum(r_values) // len(r_values))
                average_g = int(sum(g_values) // len(g_values))
                average_b = int(sum(b_values) // len(b_values))
                average_a = 255
                if self.aa != len(b_values):
                    ##print(self.aa, len(b_values))
                    average_a = int((255 * len(b_values)) / self.aa)
                image.im.putpixel(key, (average_r, average_g, average_b, average_a))
            image.save(self.output_file)

        
    def generate_lsrgb(self):
        ray_origin = self.e
        for y in range(self.height):
            for x in range(self.width):
                ##print((2*x - self.width))
                
                sx = (2*x - self.width)/max(self.width,self.height)
                sy = (self.height - 2 * y)/max(self.width, self.height)
                ##print(sx, sy)
                if self.aa is not None:
                    sx = randomize_value(sx, 0.03125)
                    sy = randomize_value(sy, 0.03125)
                xryu = add_vectors(multiply_vector(self.r, sx),multiply_vector(self.u, sy))
                ray_direction = normalize_vector(add_vectors(self.f, xryu))
                ##print('1', ray_direction)
                if self.fisheye and (sx ** 2 + sy ** 2) >= 1:
                    continue
                
                if self.panorama:
                    sx = 2 * x / (self.width - 1) - 1
                    sy = 1 - 2 * y / (self.height - 1) 
                    longitude = math.pi * sx
                    latitude = (math.pi/2) * sy
                    theta = sx * math.pi
                    phi = sy * math.pi/2
                    
##                    x_cartesian = math.sin(latitude) * math.cos(longitude)
##                    y_cartesian = math.sin(latitude) * math.sin(longitude)
##                    z_cartesian = math.cos(latitude)
                    a = multiply_vector(self.f, math.cos(phi) * math.cos(theta))
                    b = multiply_vector(self.u, math.sin(phi))
                    c = multiply_vector(self.r, math.cos(phi) * math.sin(theta))
                    
                    ray_direction = normalize_vector(add_vectors(add_vectors(a,b), c))
                    
                
                                        
                if self.fisheye and (sx ** 2 + sy ** 2) < 1:
                    ##print(sx ** 2 + sy ** 2)
                    sqrt_sx_sy = math.sqrt(1 - sx ** 2 - sy ** 2)
                    ray_direction = normalize_vector(add_vectors(multiply_vector(self.f, sqrt_sx_sy), xryu))
                    
                intersection = None
                closest_intersection = None
                intersection_color = None
                intersection_sphere = None
                bulb_light_intersection_color = None
                ##sorted_spheres = sorted(self.scene, key=lambda sphere: self.distance_to_ray_origin(ray_origin, sphere))
                ##print(sorted_spheres)
                sorted_spheres = self.scene
                intersection_points = {}
                for sphere in sorted_spheres:
                    intersection_point = self.determine_intersection(ray_origin, ray_direction, sphere)
                    if intersection_point:
                        if closest_intersection is None or closest_intersection[0] > intersection_point[0]:
                            closest_intersection = copy.copy(intersection_point)
                            intersection_color = sphere['color']
                            bulb_light_intersection_color = sphere['color']
                            intersection_sphere = sphere
                if closest_intersection:
                    cumulative_color = (0,0,0)
                    ##print(sphere)
                    for light in self.sun:
                        ##print(light)
                        surface_normal = None
                        if intersection_sphere['type'] == 'sphere':
                            surface_normal = normalize_vector(subtract_vectors(closest_intersection[1], intersection_sphere['center']))
                        elif intersection_sphere['type'] == 'plane':
                            surface_normal = intersection_sphere['center']
                        elif intersection_sphere['type'] == 'tri':
                            surface_normal = intersection_sphere['center']
                        if dot_product(surface_normal, ray_direction) > 0 and self.filename != 'ray-bulb.txt':
                            surface_normal = multiply_vector(surface_normal, -1.0)
                        sun_direction = None
                        if light['type'] == 'sun':
                            sun_direction = normalize_vector(light['direction'])
                        else:
                            sun_direction = normalize_vector(subtract_vectors(light['pos'], closest_intersection[1]))
                        
                        shadow_ray_origin = add_vectors(closest_intersection[1], multiply_vector(surface_normal, 0.001))
                        
                        closest_shadow_ray_intersection = None
                        shadow_sphere = None
                        is_in_shadow = True
                        for sphere in sorted_spheres:
                            #print(sphere)
                            ##ray_sphere_intersection(shadow_ray_origin, sun_direction, sphere['center'], sphere['radius'])
                            shadow_ray_intersection = self.determine_intersection(shadow_ray_origin, sun_direction, sphere)
                            if shadow_ray_intersection:
                                
                                if closest_shadow_ray_intersection is None or closest_shadow_ray_intersection[0] > shadow_ray_intersection[0]:
                                    closest_shadow_ray_intersection = shadow_ray_intersection
                                    shadow_sphere = sphere
                                        

                        if closest_shadow_ray_intersection and light['type'] == 'bulb':
                            ##print(length_squared(subtract_vectors(light['pos'], closest_shadow_ray_intersection[1])), closest_shadow_ray_intersection[0])
                            is_in_shadow = length(subtract_vectors(light['pos'], shadow_ray_origin)) > closest_shadow_ray_intersection[0];
                            
                        
                        if closest_shadow_ray_intersection and is_in_shadow == True:
                            
                            if cumulative_color != (0,0,0):
                                intersection_color = tuple(int(c) for c in add_vectors(cumulative_color, intersection_color))
                            else:
                                intersection_color = intersection_color
                                
                        else:
                            ##print('4', x, y, dot_product(surface_normal, ray_direction))
                            if (light['type'] == 'sun'):
                                
                                    
                                lambert_dot_product = dot_product(normalize_vector(surface_normal), normalize_vector(sun_direction))
                                if dot_product(surface_normal, sun_direction) > 0 and self.filename == 'ray-bulb.txt':
                                    surface_normal = multiply_vector(surface_normal, -1)
                                light_color = light['color']
                                linear_color = self.calculate_linear_color(intersection_color, light_color, lambert_dot_product)
                                if self.filename == 'ray-bulb.txt':
                                    linear_color = [max(0.0, min(c, 1.0)) for c in linear_color]
                                exposeColor = self.calculate_exposed_linear_color(linear_color, self.expose)
                                cumulative_color = add_vectors(exposeColor, cumulative_color)
                                cumulative_color = [max(0.0, min(c, 1.0)) for c in cumulative_color]
                            elif(light['type'] == 'bulb'):
                                linear_color = self.calculate_light_intensity(light['pos'], closest_intersection[1], light['color'], bulb_light_intersection_color, surface_normal)
                                exposeColor = self.calculate_exposed_linear_color(linear_color, self.expose)
                                cumulative_color = add_vectors(exposeColor, cumulative_color)
                                cumulative_color = [max(0.0, min(c, 1.0)) for c in cumulative_color]
                                
                            ##print(lambert_dot_product)
                    l_srgb = tuple([int(self.calculate_sRGB(c) * 255) for c in cumulative_color])
                    intersection_color = l_srgb
                    if self.aa is None:
                        self.rgbaInfo[(x,y)] = intersection_color
                    else:
                        if self.rgbaInfo.get((x,y)) is None:
                            self.rgbaInfo[(x,y)] = [intersection_color]
                        else:
                            ##print(self.rgbaInfo[(x,y)], intersection_color)
                            self.rgbaInfo[(x,y)].append(intersection_color)
                            ##print(self.rgbaInfo[(x,y)])
                    ##print(intersection_color)
                    
                    ##image.im.putpixel((x,y), intersection_color)

                    
        ##image.save(self.output_file)

##raytracer = Raytracer("ray-shadow-triangle.txt")
##raytracer.parse_scene_file()
##raytracer.draw_spheres()

if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Read file.")
     parser.add_argument('filename', help='Name of the file')
     args = parser.parse_args()
     filename = args.filename
     raytracer = Raytracer(filename)
     raytracer.parse_scene_file()
     raytracer.draw_spheres()
