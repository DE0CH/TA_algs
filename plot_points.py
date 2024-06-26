# %%
with open('points.txt') as f:
    points = []
    for line in f:
        if line.strip():
            points.append(list(map(float, line.split())))
            d = len(points[0])

n = len(points)

# %%
with open('progress.txt') as f:
    progress = []
    lines = f.readlines()
    lines = lines[:-1]
    for line in lines:
        if line.strip():
            progress.append(list(map(int, line.split())))

for id_ in range(d):
    for i in range(n):
        k = progress[i][id_] - 1
        progress[i][id_] = points[k][id_]

for i in range(n):
    progress[i].append(i/n)

# %%
progress_t = []
for id_ in range(d+1):
    a = []
    for i in range(n):
        a.append(progress[i][id_])
    progress_t.append(a)

# %%
import ipyvolume as ipv
import numpy as np
ipv.pylab.style.set_style_dark()

#%%
x, y, z = progress_t
x = np.array(x)
y = np.array(y)
z = np.array(z)

# %%
fig = ipv.figure()
ipv.plot(x, y, z, marker='sphere')
ipv.show()
