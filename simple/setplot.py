
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = linspace(0,2,101)
        yl = 0.6 + 0.1*xl
        plot(xl,yl,'k')
        yl = 0.4 + 0.3*(xl-0.5) - 0.2*(xl-0.5)**2
        plot(xl,yl,'k')
    

    def sigmatr(current_data):
        q = current_data.q
        return q[0,:,:] + q[1,:,:]

    def div(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            d = zeros(u.shape)
            return d
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        ux = (u[I+1,:][:,J] - u[I-1,:][:,J]) / (2*dx)
        vy = (v[:,J+1][I,:] - v[:,J-1][I,:]) / (2*dy)
        dint = ux + vy
        
        #zx = zeros((mx-2,1))
        #zy = zeros((1,my))
        #d = vstack((zy, hstack((zx, ux+vy, zx)), zy))
        
        d0 = dint[:,0]
        d1 = dint[:,-1]
        d2 = vstack((d0, dint.T, d1)).T
        d0 = d2[0,:]
        d1 = d2[-1,:]
        d = vstack((d0,d2,d1))      
        return d

    def logu(current_data):
        from numpy import log,exp,abs

        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        c=log(abs(u)+exp(1.0e-16))
        return c

    def logv(current_data):
        from numpy import log,exp,abs

        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        c=log(abs(v)+exp(1.0e-16))
        return c

    def curl(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            c = zeros(u.shape)
            return c
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        vx = (v[I+1,:][:,J] - v[I-1,:][:,J]) / (2*dx)
        uy = (u[:,J+1][I,:] - u[:,J-1][I,:]) / (2*dy)
        cint = vx - uy

        c0 = cint[:,0]
        c1 = cint[:,-1]
        c2 = vstack((c0, cint.T, c1)).T
        c0 = c2[0,:]
        c1 = c2[-1,:]
        c = vstack((c0,c2,c1))      
        return c

    # Figure for trace(sigma) and sigma_12 side by side
    plotfigure = plotdata.new_plotfigure(name='P and S waves', figno=11)
    plotfigure.kwargs = {'figsize':(15,10)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.5,.35,.35])' # 'subplot(221)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'u'
    plotaxes.scaled = True
    #plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = sigmatr
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -3.0e-3
    plotitem.pcolor_cmax = 3.0e-3
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.5,.5,.35,.35])' # 'subplot(222)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Stress'
    plotaxes.scaled = True
    #plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 4
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -3.0e-3
    plotitem.pcolor_cmax = 3.0e-3
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.35,.35])' # 'subplot(223)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'x disp'
    plotaxes.scaled = True
    #plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 5
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -3.0e-3
    plotitem.pcolor_cmax = 3.0e-3
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    # Figure for curl:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.5,.1,.35,.35])' # 'subplot(224)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'y disp'
    plotaxes.scaled = True
    #plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 6
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -3.0e-3
    plotitem.pcolor_cmax = 3.0e-3
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]



    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [1,0,0]
    plotitem.amr_patchedges_show = [1]


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Horizontal velocity'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Vertical velocity'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.plotstyle = 'b-'



    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
