from section import *
from scipy.optimize import root,root_scalar
import numpy as np
from beam import *


if __name__ == '__main__':

    # solution = root(flange_function,[200,200])
    # print(solution.x)
    # print(start_of_flange_function(20))
    # for section in beam:
    #     print(round(section.x),round(section.top_flange_width),round(section.bot_flange_width))
    #
    # print(beam_load_function(1))
    # print(beam_load_function(100000))

    # beam_a = Beam(1)
    # beam_tabler(beam_a)




    #6kg bridge
    max_load = root_scalar(beam_load_function,bracket=[1,100000],method='brentq').root
    beam = Beam(max_load)
    beam_tabler(beam)

    # beam_tabler(beam)
    # print(beam_load_function(0))
    # print(beam_load_function(1000000))