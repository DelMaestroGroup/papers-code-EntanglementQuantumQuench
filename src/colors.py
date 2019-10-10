# A set of colors that can be easily loaded for plotting

pastel = ["#688EAF", "#FC991D", "#7DEB74", "#FA6781", "#8B981D", "#BB7548",
          "#AD8FE4", "#96E4AA", "#D669B0", "#E1C947", "#A78200", "#7C9FE4",
          "#957DA6", "#75BF38", "#C3B059", "#51C17A", "#79AEBB", "#2790AC",
          "#688ECE", "#749DB7"]*4

ggplot = ["#E24A33", "#348ABD", "#988ED5", "#777777", "#FBC15E", "#8EBA42",
          "#FFB5B8"]*4

cbqual = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
          "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"]*4

resonant = ["#2078B5", "#FF7F0F", "#2CA12C", "#D72827", "#9467BE", "#8C574B",
            "#E478C2", "#808080", "#BCBE20", "#17BED0", "#AEC8E9", "#FFBC79", 
            "#98E08B", "#FF9896", "#C6B1D6", "#C59D94", "#F8B7D3", "#C8C8C8", 
           "#DCDC8E", "#9EDAE6"]*2

guzik = ['#818181','#6198c9','#c0e7f6','#f1f293','#fec06a','#ed5732','#57af5e','#9a6ead']*2
guzikd = {'grey':'#818181','blue':'#6198c9', 'lightblue':'#c0e7f6',
          'yellow':'#f1f293', 'orange':'#fec06a', 'red':'#ed5732', 
          'green':'#57af5e', 'purple':'#9a6ead'}

#-------------------------------------------------------------------------------
def get_cycle_colors():
    '''Get the current color cycle.'''

    import matplotlib.pyplot as plt
    return plt.rcParams['axes.prop_cycle'].by_key()['color']

#-------------------------------------------------------------------------------
def get_linear_colors(cmap,num_colors,reverse=False):
    '''Return num_colors colors in hex from the colormap cmap.'''
    
    from matplotlib import cm
    from matplotlib import colors as mplcolors
    import numpy as np

    cmap = cm.get_cmap(cmap)

    colors_ = []
    for n in np.linspace(0,1.0,num_colors):
        colors_.append(mplcolors.to_hex(cmap(n)))

    if reverse:
        colors_ = colors_[::-1]
    return colors_

#-------------------------------------------------------------------------------
def hex_to_rgb(value,transmit=None, full=False):
    '''Convert a hex color to rgb tuple.'''
    value = value.lstrip('#')
    lv = len(value)
    step = int(lv/3) 
    scale = 1.0/255.0
    if full:
        scale = 1
    col = tuple(scale*int(value[i:i+step], 16) for i in range(0, lv, step))

    if not transmit:
        return col
    else:
        return col + (transmit,)

#-------------------------------------------------------------------------------
def rgb_to_hex(value):
    '''Convert a rgb tuple to a hex color string.'''
    if value[0] < 1:
        scale = 255
    else:
        scale = 1
    rgb = [int(scale*k) for k in value]
    return '#%02x%02x%02x' % (rgb[0],rgb[1],rgb[2]) 

#-------------------------------------------------------------------------------
# https://www.viget.com/articles/equating-color-and-transparency
def get_alpha_hex(value,alpha):
    '''Convert a hex color to an equivalent non-transparent version.'''

    #first we get the rgb
    rgb = hex_to_rgb(value)

    # apply the transparency
    target = [alpha*k + (0.999-alpha) for k in rgb] 

    return rgb_to_hex(target)
