#beam object and related functions

from tabulate import tabulate
from engineeringProps import *
from scipy.optimize import root,root_scalar
from scipy import linspace,integrate
from section import *
from geometry import *
import numpy as np

#bending force along beam. Assume loading is a UDL the length of the loading plate in the center of span
#Inputs: x coordinate along beam (mm), load applied by hydraulic press (N)
#Output: bending force at x (kNmm)
def moment(x,load):
    unloaded_length = (span-plate_length)/2
    if x<(span-plate_length)/2:
        return load*x/2
    elif x>(span+plate_length)/2:
        return load*((span+plate_length)/2-x)/2+load*unloaded_length/2
    else:
        return (load*unloaded_length)/2 + 0.5 * (load/2-load/plate_length*(x-unloaded_length)+load/2)*(x-unloaded_length)

#fucntion used to solve where the web is insuffient in bending and flanges need to start
#Outputs: difference between the stress in web and the web strength as a 2 element array where the 1st element is the
# stress difference at the top and 2nd element is the stress difference at the bottom

def start_of_flange_function(x,load,web_thickness,web_height):
    section = Section(x,moment(x,load),0,top_flange_thickness,0,bot_flange_thickness,web_thickness,web_height)
    net_top_stress = section.top_flange_stress - web_comp_strength
    net_bot_stress = section.bot_flange_stress - web_ten_strength
    return min(-net_top_stress,net_bot_stress)

#function used to solve the flange widths
#Outputs: difference between the stress at the top and bottom of the flanges and the flange strengths as a 2 element array
def flange_function(flange_widths,x,load,web_thickness,web_height):
    if flange_widths[0] >= 0 and flange_widths[1] >= 0:
        section = Section(x,moment(x,load),flange_widths[0],top_flange_thickness,flange_widths[1],bot_flange_thickness,web_thickness,web_height)
        net_top_stress = section.top_flange_stress - comp_strength
        net_bot_stress = section.bot_flange_stress - ten_strength

        return([-net_top_stress,net_bot_stress])

    else:
        return([-100,-100])

#approximates the mass of the beam by using Simpson intregration of the areas of the defined sections
#Inputes: Beam object
#Outputs: mass (g)
def beam_mass(beam):
    x_arr=[]
    y_arr=[]
    for section in beam.sections:
        x_arr.append(section.x)
        y_arr.append(section.area)

    volume = 2 * integrate.simpson(y_arr,x_arr)
    return volume*density

#creates a beam object. Beam will be sized as to resist the applied load
#Input: load, design load of beam
#Output: beam object
class Beam:
    def __init__(self,load):
        self.load = load
        # solve for web thickness
        self.shear_force = 0.5 * self.load
        self.web_thickness = 3 / 2 * self.shear_force / (web_height * shear_strength)

        self.start_of_flange = root_scalar(start_of_flange_function, method='secant', x0=0, x1=1,
                                           args=(load, self.web_thickness, web_height)).root

        #number of sections solved for in beam
        self.number_of_sections = number_of_sections
        #x coordinate of sections. Beam is symetrical so only half the beam is modelled
        x_arr = np.linspace(self.start_of_flange, span / 2, self.number_of_sections-1)

        self.sections = []
        #section at the start of the beam
        self.sections.append(Section(x=0, moment=moment(0, self.load),
                            top_flange_width=0, top_flange_thickness=top_flange_thickness,
                            bot_flange_width=0, bot_flange_thickness=bot_flange_thickness,
                            web_thickness=self.web_thickness, web_height=web_height))
        #section where the flanges start
        self.sections.append(Section(x=self.start_of_flange, moment=moment(self.start_of_flange, load),
                            top_flange_width=0, top_flange_thickness=top_flange_thickness,
                            bot_flange_width=0, bot_flange_thickness=bot_flange_thickness,
                            web_thickness=self.web_thickness, web_height=web_height))

        #remaining sections. Solve for flange widths by solving the roots of the flange functions.
        #this could be improved by solving for the section that is strong enough with the minimum area. It is currently
        #unknown how unoptimal the current result is
        for x_coor in x_arr[1:]:
            flange_widths = root(flange_function, [0, 0], args=(x_coor, load, self.web_thickness, web_height)).x
            self.sections.append(Section(x=x_coor, moment=moment(x_coor, load),
                                top_flange_width=flange_widths[0], top_flange_thickness=top_flange_thickness,
                                bot_flange_width=flange_widths[1], bot_flange_thickness=bot_flange_thickness,
                                web_thickness=self.web_thickness, web_height=web_height))
        self.mass = beam_mass(self)


#outputs a table of the relevant beam properties
def beam_tabler(beam):

    table=[[beam.web_thickness,web_height]]
    headers = ["web thickness","web height"]
    print(tabulate(table,headers,tablefmt="simple_grid",floatfmt=".1f"))

    table = [[top_flange_thickness, bot_flange_thickness]]
    headers = ["TF thickness", "BF thickness"]
    print(tabulate(table, headers, tablefmt="simple_grid", floatfmt=".1f"))

    table=[]
    for section in beam.sections:
        table.append([section.x,section.top_flange_width,section.bot_flange_width])
    headers = ["x coor",'TF width',"BF width"]
    print(tabulate(table,headers,tablefmt="simple_grid",floatfmt=".1f"))
    print("Beam mass = ",round(beam.mass),"g")
    print("Column mass = ",round(column_mass(beam.load)),"g")
    print("Max load = ",round(beam.load),"N")

#appromixate column mass
def column_mass(load):
    return abs(load/comp_strength * column_height * density)

#function used to solve for a beam weighing 6kg (-1kg for unaccounted for details)
def beam_load_function(load):
    beam = Beam(load)
    return 6000 - beam_mass(beam) - column_mass(load)-1000
