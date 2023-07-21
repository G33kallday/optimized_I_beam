from engineeringProps import *

#section properties of rectangles
class Rectangle:
    def __init__(self,width,height,centroid):
        self.width=width
        self.height=height
        self.centroid=centroid

        self.area = width*height
        self.first_moment = self.area * centroid
        self.inertia = 1/12 * width * height ** 3

#properties of I section
class Section:
    def __init__(self,x,moment,top_flange_width,top_flange_thickness,bot_flange_width,bot_flange_thickness,web_thickness,web_height):
        #x coordinate of section
        self.x=x
        self.moment=moment
        self.top_flange_width=top_flange_width
        self.top_flange_thickness=top_flange_thickness
        self.bot_flange_width=bot_flange_width
        self.bot_flange_thickness=bot_flange_thickness
        self.web_thickness=web_thickness
        self.web_height=web_height

        n_ratio = flange_modulus/web_modulus
        self.top_flange_eff_width = (top_flange_width-web_thickness) * n_ratio
        self.bot_flange_eff_width = (bot_flange_width - web_thickness) * n_ratio

        self.top_flange = Rectangle(self.top_flange_eff_width,top_flange_thickness,web_height-top_flange_thickness/2)
        self.web = Rectangle(web_thickness,web_height,web_height/2)
        self.bot_flange = Rectangle(self.bot_flange_eff_width,bot_flange_thickness,bot_flange_thickness/2)

        self.rectangles = [self.top_flange,self.web,self.bot_flange]

        self.sum_first_moment = 0
        self.area = 0
        for rectangle in self.rectangles:
            self.sum_first_moment+= rectangle.first_moment
            self.area += rectangle.area
        self.centroid = self.sum_first_moment/self.area

        trans_inertia = 0
        for rectangle in self.rectangles:
            trans_inertia += rectangle.inertia + rectangle.area * (rectangle.centroid - self.centroid) ** 2
        self.trans_inertia = trans_inertia

        self.top_flange_stress = -self.moment * (web_height - self.centroid) / trans_inertia * n_ratio
        self.bot_flange_stress = self.moment * self.centroid / trans_inertia * n_ratio