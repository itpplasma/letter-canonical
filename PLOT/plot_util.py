import numpy as np
import matplotlib.pyplot as plt

shape = (75, 64, 100)
shape_spl = (74, 63, 49)

phi = np.linspace(0, 2 * np.pi, shape[0])

def doplot(data, title, rname="R"):
    # Plot 2D slice as contour
    try:
        plt.figure()
        plt.contour(data[:, :, 38])
        plt.colorbar()
        plt.xlabel(r"$\varphi$")
        plt.ylabel("Z")
        plt.title(title)
    except:
        pass

    try:
        plt.figure()
        plt.contour(data[:, 1, :])
        plt.colorbar()
        plt.xlabel(rname)
        plt.ylabel("Z")
        plt.title(title)
    except:
        pass

    try:
        plt.figure()
        plt.contour(data[1, :, :])
        plt.colorbar()
        plt.xlabel(rname)
        plt.ylabel(r"$\varphi$")
        plt.title(title)
    except:
        pass

    try:
        plt.figure()
        plt.plot(data[30, 30, :])
        plt.xlabel(rname)
        plt.title(title)
    except:
        pass

def doplot3d(data, title, rname="R"):
    import plotly.graph_objects as go

    X, Y, Z = np.mgrid[0:1:shape[0]*1j, 0:1:shape[1]*1j, 0:1:shape[2]*1j]

    # Plot 3D isosurfaces
    fig = go.Figure(data=go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=data.flatten(),
        isomin=data.min(),
        isomax=data.max(),
        opacity=0.3,
        surface_count=10,
    ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(title="Z"),
            yaxis=dict(title="φ"),
            zaxis=dict(title=rname),
        )
    )

    fig.show()
