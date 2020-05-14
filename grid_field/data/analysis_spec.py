import matplotlib as mp
import pandas as pd

posx = pd.read_csv("XPOS.csv")
posy = pd.read_csv("YPOS.csv")
posz = pd.read_csv("ZPOS.csv")
px   = pd.read_csv("MOMENTUM_X.csv")
py   = pd.read_csv("MOMENTUM_Y.csv")
pz   = pd.read_csv("MOMENTUM_Z.csv")
time = pd.read_csv("TIME.csv")
energy = pd.read_csv("ENERGY.csv")


time