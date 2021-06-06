# PYTHON PACKAGES
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy
import utils

class Environment:
    def __init__(self, name='Xenon129 Environment'):
        self.name = name
        self.B0 = 0.  # DC z magnetic field             [Tesla]
        self.Bnoise = 0.  # z magnetic field noise          [Tesla]

        # System parameters
        self.wd_x = 0.  # x drive frequency               [rad  s^-1]
        self.Ad_x = 0.  # x drive amplitude               [s^-1]
        self.wd_y = 0.  # y drive frequency               [rad  s^-1]
        self.Ad_y = 0.  # y drive amplitude               [s^-1]

        # world rotation
        self.wr = 0.  # world rotation around z axis    [rad  s^1]

    def set_state(self, wr=0, B0=0, Bnoise=0, wd_x=0, Ad_x=0, wd_y=0, Ad_y=0):
        self.wr = wr
        self.B0 = B0
        self.Bnoise = Bnoise
        self.wd_x = wd_x
        self.Ad_x = Ad_x
        self.wd_y = wd_y
        self.Ad_y = Ad_y

    def display_params(self):
        print('===================================================================')
        print(f'| {self.name}:')
        print(f'| ----------')
        print(f'| B0:                     {self.B0}')
        print(f'| B_noise:                {self.Bnoise}')
        print(f'| {phy.OMEGA}d_x:                   {self.wd_x}')
        print(f'| {phy.BIG_OMEGA}d_x:                   {self.Ad_x}')
        print(f'| {phy.OMEGA}d_y:                   {self.wd_y}')
        print(f'| {phy.BIG_OMEGA}d_y:                   {self.Ad_y}')
        print(f'| {phy.BIG_OMEGA}r:                     {self.wr}')
        print('===================================================================')
