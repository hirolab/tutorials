# %%
"""from IPython import get_ipython

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
get_ipython().run_line_magic('matplotlib', 'qt')
"""
# %%
import numpy as np 
from kalman_filter import SimpleUAV3D

kf = SimpleUAV3D()
kf.X[3] = -1

kf.predict(np.array([[1], [0], [0]]))
print(kf.X, kf.P)
    
kf.correct(np.array([[0.9], [0], [0]]))
kf.X, kf.P 

kf.plot_state()

# %% Animation growing covariance
import numpy as np 
from kalman_filter import SimpleUAV3D
import matplotlib.pyplot as plt 

T = 0.1
kf = SimpleUAV3D(T=T, Q=(np.diag([T**2/0.1, T**2/2, T**2/2, T, T, T])*10)**2)

ax = kf.plot_state()


for t in np.arange(0, 3, T):
    kf.predict(np.array([[0.5], [0], [0]]))
    kf.plot_state(ax)
    plt.pause(0.1)
    
