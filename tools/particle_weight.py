#!/usr/bin/env python

# Calculate the particles weight

import math

particle = {"x": 0,
            "y": 5}
landmark = {"x": 2, 
            "y": 1}

sigma_x = 0.3
sigma_y = 0.3

normalizer = 1 / (2 * math.pi * sigma_x * sigma_y)

exponent = - (math.pow((particle["x"] - landmark["x"]), 2) / (2 * math.pow(sigma_x, 2)) +
              math.pow((particle["y"] - landmark["y"]), 2) / (2 * math.pow(sigma_y, 2)))
p_weight = normalizer * math.exp(exponent)

print("particle weigth=", p_weight)