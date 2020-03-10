import plotly.graph_objects as go
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


# https://stackoverflow.com/questions/46571624/sorting-points-from-distance-to-a-given-point-x-y-here-in-my-case-x-0-y-o
def e_dist(a, b, metric='euclidean'):
    """Distance calculation for 1D, 2D and 3D points using einsum
    : a, b   - list, tuple, array in 1,2 or 3D form
    : metric - euclidean ('e','eu'...), sqeuclidean ('s','sq'...),
    :-----------------------------------------------------------------------
    """
    a = np.asarray(a)
    b = np.atleast_2d(b)
    a_dim = a.ndim
    b_dim = b.ndim
    if a_dim == 1:
        a = a.reshape(1, 1, a.shape[0])
    if a_dim >= 2:
        a = a.reshape(np.prod(a.shape[:-1]), 1, a.shape[-1])
    if b_dim > 2:
        b = b.reshape(np.prod(b.shape[:-1]), b.shape[-1])
    diff = a - b
    dist_arr = np.einsum('ijk,ijk->ij', diff, diff)
    if metric[:1] == 'e':
        dist_arr = np.sqrt(dist_arr)
    dist_arr = np.squeeze(dist_arr)
    return dist_arr


scaler = MinMaxScaler(feature_range=(0.5, 1.0))

points = 500
x = np.random.rand(points) * 100
y = np.random.rand(points) * 100

a = np.array([x, y]).T
idx = np.argsort(e_dist(a, [50, 50]))
closest = a[idx]
z = sorted(scaler.fit_transform(list(np.random.normal(0.7, 0.2, (502, 1)))).flatten(), reverse=True)[1:-1]

df = pd.DataFrame({"x": closest[:, 0], "y": closest[:, 1], "c": z})

x = df["x"]
y = df["y"]
z = df["c"]

fig = go.Figure(
    data=[
        go.Mesh3d(
            x=x, y=y, z=z,
            colorscale="Jet",
            colorbar_title="MCC",
            intensity=z,
            intensitymode="vertex",
            showscale=False)])

fig.add_trace(
    go.Scatter3d(
        x=x, y=y, z=z,
        marker=dict(size=1, color="black"),
        mode="markers",
        error_z={"type": "data",
                 "symmetric": False,
                 "array": len(z) * [0.05],
                 "arrayminus": len(z) * [0.0],
                 "thickness": 0.5}))

fig.update_layout(
    template="plotly_white",
    scene=dict(
        zaxis=dict(range=[0.0, 1.0]),
        xaxis_title='Encoding X',
        yaxis_title='Encoding Y',
        zaxis_title='MCC (+/-std)'),
    # width=700,
    margin=dict(r=20, b=10, l=10, t=10))

# fig.show()
fig.write_image("encoding_landscape.svg")
