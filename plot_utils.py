""""
plot_utils.py



Plotting Utilities (assumes Plotly)

Functions
---------
makeSphere                  : return coordinates to draw a sphere in 3D
makeCircle                  : return coordinates to draw a circle in 2D
quickPlotTraj               : quickly plot a list of trajectories
addTrajToPlot               : add trajectory data to plot
plotLowThrustControlHistory : plot low-thrust control history
plotPoincareMap             : plot a basic 2-D Poincare Map
savePlotlyFig               : save a plotly figure
"""

import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np


def makeSphere(x, y, z, r, resolution=20):    
    """
    Return the coordinates for plotting a sphere centered at (x,y,z) with radius r
    
    Parameters
    ----------
    x : float 
        sphere center x-coordinate
    y : float 
        sphere center y-coordinate
    z : float 
        sphere center z-coordinate
    r : float 
        sphere radius
    resolution : int 
        sphere resolution

    Returns
    -------
    (X,Y,Z) : tuple
        sphere x,y,z coordinates
    """
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = r * np.cos(u)*np.sin(v) + x
    Y = r * np.sin(u)*np.sin(v) + y
    Z = r * np.cos(v) + z
    return (X, Y, Z)


def makeCircle(x, y, r):
    """ 
    Return coordinates for plotting a circle centered at (x,y) with radius r
    
    Parameters
    ----------
    x : float 
        circle center x-coordinate
    y : float 
        circle center y-coordinate
    r : float 
        circle radius

    Returns
    -------
    (X,Y) : tuple 
        circle x,y coordiantes
    """
    t = np.linspace(0, 2*np.pi, 1000)
    X = x + r*np.cos(t)
    Y = y + r*np.sin(t)
    return (X,Y)


def quickPlotTraj(trajList: list, ode=None, plot_P1=False, plot_P2=True, 
                  hide_legend=False, toScaleP1P2=False, p1Rad:float=-1.0,
                  p2Rad:float=-1.0, plot_lpoints={})->go.Figure:
    """
    Quickly plot a list of trajectories

    Parameters
    ----------
    trajList : list 
        trajectory plot list. note: expects a list of np.array
        that ASSET typically returns after integration, solving,
        or optimizing
    ode : asset_asrl.OptimalControl.OdeBase 
        ode class to plot celestial bodies or lagrange point locations
    plot_P1 : bool
        flag to indicate whether or not to plot P1
    plot_P2 : bool
        flag to indicate whether or not to plot P2
    hide_legend : bool
        flag to indicate whether or not to 
        include the plot legend
    toScaleP1P2 : bool
        flag to plot "to scale" P1 or P2
    p1Rad : float
        P1 radius (non-dim)
    p2Rad : float
        P2 radius (non-dim)
    plot_lpoints : dict
        dict containing Lagrange Points to plot,
        e.g., plot_lpoints={'L2' : True}

    Returns
    -------
    fig : ploty.graph_obects.Figure
        plotly Figure
    """
    fig = go.Figure()

    for i, t in enumerate(trajList):
        t = np.array(t).T
        fig.add_trace(go.Scatter3d(x=t[0],y=t[1],z=t[2],mode='lines',name=f'Traj {i}'))
    
    if ode is not None and plot_P1:
        try:
            if toScaleP1P2 and p1Rad > 0.0:
                x,y,z = makeSphere(ode.P1[0], 0.0, 0.0, r=p1Rad)
                fig.add_trace(go.Surface(x=x,y=y,z=z,colorscale=[[0,'blue'],[1,'blue']],showscale=False,name='P1'))
            else:
                fig.add_trace(go.Scatter3d(x=[ode.P1[0]],y=[ode.P1[1]],z=[ode.P1[2]],mode='markers',marker=dict(size=10),name='P1'))
        except:
            print('ODE class does not contain P1 location.')
    
    if ode is not None and plot_P2:
        try:
            if toScaleP1P2 and p2Rad > 0.0:
                x,y,z = makeSphere(ode.P2[0], 0.0, 0.0, r=p2Rad)
                fig.add_trace(go.Surface(x=x,y=y,z=z,colorscale=[[0,'grey'],[1,'grey']],showscale=False,name='P2'))
            else:
                fig.add_trace(go.Scatter3d(x=[ode.P2[0]],y=[ode.P2[1]],z=[ode.P2[2]],mode='markers',marker=dict(size=10),name='P2'))
        except:
            print('ODE class does not contain P2 location or does not have LStar data')
    
    if ode is not None and any(plot_lpoints.values()):
        try:
            for k,v in plot_lpoints.items():
                if v:
                    lpoint_loc = getattr(ode, k)
                    fig.add_trace(go.Scatter3d(x=[lpoint_loc[0]],y=[lpoint_loc[1]],z=[lpoint_loc[2]],mode='markers',marker=dict(size=7),name=k))
        except:
            print('ODE class does not contain Lagrange Point locations.')

    if hide_legend:
        fig.update_layout(showlegend=False)

    fig.update_layout(scene=dict(aspectmode='data'))

    return fig


def addTrajToPlot(trajList:list, fig:go.Figure, width:int=1, color:str='blue')->go.Figure:
    """
    Add trajectory data to plot

    Parameters
    ----------
    traj : list
        list of trajectories to plot
    fig : go.Figure
        fig to add traj data to
    width : int
        trajectory line width
    color : str
        trajectory line color

    Returns
    -------
    fig : go.Figure
        fig to add traj data to
    """

    for i, t in enumerate(trajList):
        t = np.array(t).T
        fig.add_trace(go.Scatter3d(x=t[0],y=t[1],z=t[2],mode='lines',line=dict(color=color,width=width),name=f'Traj {i}'))
    
    return fig


def plotLowThrustControlHistory(traj, unorm=False, tuidx=[6,7,8,9])->go.Figure:
    """
    Plot low-thrust control history

    Parameters
    ----------
    traj : list of array
        Trajectory control history to plot.
    unorm : bool
        Flag to indicate whether or not to plot the control norm 
        (default False).
    tuidx : list
        Time and control variable indexes of state vector, assumes standard 
        CR3BP low-thrust problem of ASSET format where time is 0-based index 6 
        and control is 7,8,9.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Figure of plot
    """

    time_idx  = tuidx[0]
    ctrl_idx  = tuidx[1:]

    fig = go.Figure()

    if unorm: u_norm =[np.linalg.norm(t[ctrl_idx]) for t in traj]

    plotter_traj = np.array(traj).T

    if unorm: 
        fig.add_trace(go.Scattergl(x=plotter_traj[time_idx], y=u_norm,mode='lines',line=dict(color='dodgerblue',width=2),name=r'$||U||$'))
    else:
        fig.add_trace(go.Scattergl(x=plotter_traj[time_idx],y=plotter_traj[ctrl_idx[0]],mode='lines',line=dict(color='red',width=2),name=r'$U_x$'))
        fig.add_trace(go.Scattergl(x=plotter_traj[time_idx],y=plotter_traj[ctrl_idx[1]],mode='lines',line=dict(color='green',width=2),name=r'$U_y$'))
        fig.add_trace(go.Scattergl(x=plotter_traj[time_idx],y=plotter_traj[ctrl_idx[2]],mode='lines',line=dict(color='blue',width=2),name=r'$U_z$'))
    
    fig.update_layout(xaxis=dict(showline=True,mirror=True,linecolor='black',gridcolor='lightgrey'),
                      yaxis=dict(showline=True,mirror=True,linecolor='black',gridcolor='lightgrey'),
                      legend=dict(x=0.01,y=0.99,bgcolor='rgba(0,0,0,0)'),
                      plot_bgcolor='rgba(0,0,0,0)',
                      font=dict(family='CMU Serif',color='black',size=20))
    
    return fig


def plotPoincareMap(ode, crossings:list, xIdx:int, 
                    yIdx:int, xRange:list=[], yRange:list=[], 
                    plotTitle='', cmapIdx=-1, hidelegend=False,
                    scale:bool=False)->go.Figure:

    """
    Plot a basic 2-D Poincare Map

    Parameters
    ----------
    ode : ASSET OC.ODEBase
        ODE
    crossings : list of dict
        list containing hyperplane crossing data
    xIdx : int 
        index of state variable to plot on x-axis, e.g., 0 = x
    yIdx : int
        index of state variable to plot on y-axis, e.g., 3 = vx
    xRange : list, optional
        x-axis limits of plot
    yRange : list, optional
        y-axis limits of plot
    plotTitle : str, optional
        plot title
    cmapIdx : int, optional
        index of state variable to create color map
    scale : bool, optional
        flag to dimensionalize non-dim values
    
    Returns
    -------
    fig : go.Figure
        plotly figure of 2-D poincare map
    """

    match xIdx:
        case 0:
            xlabel = 'x [nd]'
        case 1:
            xlabel = 'y [nd]'
        case 2:
            xlabel = 'z [nd]'
        case 3:
            xlabel = 'vx [nd]'
        case 4:
            xlabel = 'vy [nd]'
        case 5:
            xlabel = 'vz [nd]'
        case _:
            raise ValueError('State variable index not valid.')
        
    match yIdx:
        case 0:
            ylabel = 'x [nd]'
        case 1:
            ylabel = 'y [nd]'
        case 2:
            ylabel = 'z [nd]'
        case 3:
            ylabel = 'vx [nd]'
        case 4:
            ylabel = 'vy [nd]'
        case 5:
            ylabel = 'vz [nd]'
        case _:
            raise ValueError('State variable index not valid.')

    CMAP_ON = False
    if cmapIdx >= 0:
        CMAP_ON = True

        match cmapIdx:
            case 0:
                zlabel = 'x [nd]'
            case 1:
                zlabel = 'y [nd]'
            case 2:
                zlabel = 'z [nd]'
            case 3:
                zlabel = 'vx [nd]'
            case 4:
                zlabel = 'vy [nd]'
            case 5:
                zlabel = 'vz [nd]'
            case _:
                raise ValueError('State variable index not valid.')
            
    fig = go.Figure() 

    if scale:
        xscale = ode.vstar if 'v' in xlabel else ode.lstar
        yscale = ode.vstar if 'v' in ylabel else ode.lstar
    if scale and CMAP_ON:
        zscale = ode.vstar if 'v' in zlabel else ode.lstar

    for crossing_dict in crossings:

        event_x, event_y = [], []
        if CMAP_ON:
            event_z = []

        try:
            data   = crossing_dict['data']
            name   = crossing_dict['name']
            mcolor = crossing_dict['color']
        except KeyError as e:
            raise e('Crossings dict incorrectly defined.')

        for c in data:
            event_x.append(c[xIdx])
            event_y.append(c[yIdx])
            if CMAP_ON:
                event_z.append(c[cmapIdx])

        if scale:
            event_x *= xscale
            event_y *= yscale
            if CMAP_ON:
                event_z *= zscale

        if CMAP_ON:
            fig.add_trace(go.Scattergl(x=event_x, y=event_y, mode='markers',marker=dict(size=3, color=event_z,colorscale=px.colors.sequential.Viridis),name=name))
        else:
            fig.add_trace(go.Scattergl(x=event_x, y=event_y, mode='markers',marker=dict(size=3, color=mcolor),name=name))

        
    fig.update_layout(xaxis=dict(showline=True,showgrid=True,mirror=True,linecolor='black',gridcolor='lightgrey',title=xlabel),
                      yaxis=dict(showline=True,showgrid=True,mirror=True,linecolor='black',gridcolor='lightgrey',title=ylabel,scaleanchor='x',scaleratio=1),
                      font=dict(
                          family='CMU Serif',
                          size=15,
                          color='black'),
                      plot_bgcolor='rgba(0,0,0,0)',
                      title=dict(text=plotTitle))

    if xRange:
        fig.update_layout(xaxis=dict(range=[xRange[0], xRange[1]]))

    if yRange:
        fig.update_layout(yaxis=dict(range=[yRange[0], yRange[1]]))

    if hidelegend:
        fig.update_layout(showlegend=False)

    return fig


def savePlotlyFig(fig:go.Figure, path, fmt)->None:
    """
    Save a plotly fig

    Parameters
    ----------
    fig : go.Figure
        figure to save
    path : str or Path
        output path
    fmt : str
        output format, e.g., 'png' or 'svg'
    """

    if fmt == 'html':
        fig.write_html(path,include_plotlyjs='cdn',auto_open=False)
    else:
        fig.write_image(path, fmt, engine='kaleido')

    return


