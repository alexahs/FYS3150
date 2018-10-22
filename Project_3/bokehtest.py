from bokeh.plotting import figure, show
import numpy as np

x = np.linspace(-2*np.pi, 2*np.pi)
y = np.sin(x)

p = figure()
p.line(x, y)
show(p)
