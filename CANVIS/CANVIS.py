import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches


from scipy.stats import norm
import math
from optparse import OptionParser
import svgutils.transform as sg
import sys
#import cairosvg
import warnings
import os

def vararg_callback(option, opt_str, value, parser):
    """Function that allows for a variable number of arguments at the command line"""
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        if arg[:2] == "--" and len(arg) > 2:
            break
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

def Read_Input(locus_fname, zscore_names, ld_fnames, annotation_fname, specific_annotations, interval, bp_head):
    """Function that reads in all your data files"""
    zscore_data = pd.read_csv(locus_fname, delim_whitespace=True)
    zscores = zscore_data[zscore_names]
    location = zscore_data[bp_head]
    pos_prob = zscore_data['Posterior_Prob']

    if interval is not None: # user input an interval
        a = int(interval[0])
        b = int(interval[1])
    elif interval is None:  # user did not input an interval; set interval to whole interval
        a = np.amin(location)
        b = np.amax(location)
    if a < location[0] or a > location[len(location)-1]:
        # user input out of range interval; set interval to whole interval
        warnings.warn('Specified interval is out of range; left bound set to first valid location')
        a = np.amin(location)
    if b > location[len(location) - 1] or b < location[0]:
        warnings.warn('Specified interval is out of range; right bound set to last valid location')
        b = np.amax(location)

    indices = np.where((location >= a) & (location <= b))
    N = indices[0][0]
    M = indices[-1][-1]
    lds = []
    if ld_fnames is not None:
        for ld_fname in ld_fnames:
            ld = pd.read_csv(ld_fname, header=None, delim_whitespace=True)
            ld_matrix = ld.as_matrix()
            # calculate index for location form location
            ld_matrix = ld_matrix[N:M, N:M]
            ld = pd.DataFrame(data=ld_matrix)
            lds.append(ld)
    else:
        lds = None
    if annotation_fname is not None:
        annotation_data = pd.read_csv(annotation_fname, delim_whitespace=True)
        if specific_annotations is not None:
            annotations = annotation_data[specific_annotations]
        else: # only data; no names
            header = pd.read_csv(annotation_fname, delim_whitespace=True, header=None)
            header = header.values.tolist()
            specific_annotations = header[0]
            annotations = annotation_data[specific_annotations]
        annotations = annotations.as_matrix()
        annotations = annotations[N:M]
    else: # no data or names
        annotations = None
    zscores = zscores.as_matrix()
    zscores = zscores[N:M, :]
    pos_prob = pos_prob.as_matrix()
    pos_prob = pos_prob[N:M]
    location = location.as_matrix()
    location = location[N:M]

    return [zscores,pos_prob,location, lds, annotations, specific_annotations]

def Zscore_to_Pvalue(zscore):
    """Function that converts zscores to pvalues"""
    abs_zscore = np.absolute(zscore)
    pvalue = -1 * (norm.logsf(abs_zscore) / math.log(10))
    return pvalue

# Find the top SNP and return the vector of SNPs relative to it

def Find_Top_SNP(zscore_vect, correlation_matrix):
    correlation_matrix = correlation_matrix.as_matrix()
    # use r^2
    correlation_matrix = np.square(correlation_matrix)
    zscore_vect = np.absolute(zscore_vect)
    top_SNP = zscore_vect.argmax() # returns index
    # get column corresponding to top SNP
    top_vect = correlation_matrix[:][top_SNP]
    return top_vect, top_SNP

# Zscores Plot

def Plot_Statistic_Value(position, zscore, zscore_names, greyscale, lds):
    """function that plots pvalues from given zscores"""

    zscore_tuple = []
    for i in range(0, len(zscore_names)):
        fig = plt.figure(figsize=(6, 3.25))
        sub = fig.add_subplot(1,1,1, facecolor='white')
        plt.xlim(np.amin(position), np.amax(position) + 1)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.ylabel('-log10(pvalue)', fontsize=10)
        z = zscore[:, i]
        pvalue = Zscore_to_Pvalue(z)

        if lds is not None:
            if i < len(lds): # exists a corresponding LD
                correlation_matrix = lds[i]
                [top_vect, top_SNP] = Find_Top_SNP(z, correlation_matrix)

            else: # no corresponding LD, so use previously calculated one
                warnings.warn("Warning: no corresponding LD matrix for zscore. Plot is made using previous LD matrix.")
                n = len(lds) - 1
                correlation_matrix = lds[n]
                [top_vect, top_SNP] = Find_Top_SNP(z, correlation_matrix)

            if greyscale == 'y':
                sub.scatter(position, pvalue, c=top_vect, cmap='Greys', zorder=1, clip_on=False)
            else:
                sub.scatter(position, pvalue, c=top_vect, cmap='jet', zorder=1, clip_on=False)
            x = position[top_SNP]
            y = pvalue[top_SNP]
            sub.plot(x, y, marker='D', color='black', zorder=2, clip_on=False)
        else:
            if greyscale == "y":
                sub.scatter(position, pvalue, color='#6B6B6B', clip_on=False)
            else:
                sub.scatter(position, pvalue, color='#D64541', clip_on=False)

        plt.gca().set_ylim(bottom=0)
        #add threshold line at 5*10-8
        x = [np.amin(position), np.amax(position) + 1]
        y = [-1*math.log(5*10**-8)/(math.log(10)), -1*math.log(5*10**-8)/(math.log(10))]
        plt.plot(x,y,'gray', linestyle='dashed', clip_on=False)
        label = mpatches.Patch(color='#FFFFFF', label=zscore_names[i])
        legend = plt.legend(handles=[label])
        for label in legend.get_texts():
            label.set_fontsize('small')
        value_plot = fig

        bar = plt.figure()

        if lds is not None:
            # add color bar
            min_value = np.amin(top_vect)
            max_value = np.amax(top_vect)
            fig = plt.figure(figsize=(3, 1.0))
            ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
            if greyscale == 'y':
                cmap = mpl.cm.binary
            else:
                cmap = mpl.cm.jet
            norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
            mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
            bar = fig

        zscore_tuple.append((value_plot,bar))

    return zscore_tuple

def Plot_Position_Value(position, pos_prob, threshold, greyscale):
    """Function that plots z-scores, posterior probabilites, other features """
    if greyscale == "y":
        plot_color = '#BEBEBE'
        set_color = '#000000'
    else:
        plot_color = '#2980b9'
        set_color = '#D91E18'
    [credible_loc, credible_prob] = Credible_Set(position, pos_prob, threshold)
    fig = plt.figure(figsize=(6, 3.25))
    sub1 = fig.add_subplot(1,1,1, facecolor='white')
    plt.xlim(np.amin(position), np.amax(position)+1)
    plt.ylabel('Posterior probabilities', fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlabel('Location', fontsize=10)
    sub1.scatter(position, pos_prob, color=plot_color, clip_on=False)
    if threshold != 0:
        sub1.scatter(credible_loc, credible_prob, color=set_color, label='Credible Set', clip_on=False)
        title = "Credible Set: " + str(threshold*100) + "%"
        credible_set = mpatches.Patch(color=set_color, label=title)
        legend = plt.legend(handles=[credible_set])
        for label in legend.get_texts():
            label.set_fontsize(10)
    plt.gca().set_ylim(bottom=0)
    value_plots = fig
    return value_plots

def Credible_Set(position, pos_prob, threshold):
    """Function that finds the credible set according to a set threshold"""
    total = sum(pos_prob)
    bounds = threshold*total
    #make into tuples
    tuple_vec = []
    for i in range(0, len(position)):
        tup = (position[i], pos_prob[i])
        tuple_vec.append(tup)
    #order tuple from largest to smallest
    tuple_vec = sorted(tuple_vec, key=lambda x: x[1], reverse=True)
    credible_set_value = []
    credible_set_loc = []
    total = 0
    for tup in tuple_vec:
        total += tup[1]
        credible_set_loc.append(tup[0])
        credible_set_value.append(tup[1])
        if total > bounds:
            break
    return credible_set_loc, credible_set_value

def Plot_Heatmap(lds, greyscale, large_ld):
    """Function that plots heatmap of LD matrix"""
    ld_arr = []
    for correlation_matrix in lds:
        n = correlation_matrix.shape
        correlation_matrix = np.square(correlation_matrix)
        if n[0] > 350 and large_ld == 'n':
            warnings.warn('LD matrix is too large and will not be produced. To override, add "--L y"')
            heatmap = None
            bar = None
            ld_arr.append((heatmap, bar))
            break

        if n[0] > 350 and large_ld == 'y':
            warnings.warn('LD matrix is too large but will be produced due to override flag')

        fig = plt.figure(figsize=(3.25, 3.25))
        sns.set(style="white")
        correlation = correlation_matrix.corr()
        mask = np.zeros_like(correlation, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True
        if greyscale == "y":
            cmap = sns.light_palette("black", as_cmap=True)
        else:
            cmap = None
        sns.heatmap(correlation, mask=mask, cmap=cmap, square=True,
                    linewidths=0, cbar=False, xticklabels=False, yticklabels=False, ax=None)
        heatmap = fig

        matrix = correlation_matrix.as_matrix()
        min_value = np.amin(matrix)
        max_value = np.amax(matrix)
        fig = plt.figure(figsize=(3, 1.0))
        ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        if greyscale == 'y':
            cmap = mpl.cm.binary
        else:
            cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
        mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
        bar = fig
        ld_arr.append((heatmap, bar))
    return ld_arr

def Plot_Annotations(annotation_names, annotation_vectors, greyscale):
    """Plot the annotations with labels"""
    annotation_tuple = []
    for i in range(0, len(annotation_names)):
        annotation = annotation_vectors[:,i]
        colors = []
        if greyscale == "y":
            for a in annotation:
                if a == 1:
                    colors.append('#000000')
                else:
                    colors.append('#FFFFFF')
        else:
            color_array = ['#2980b9']
            for a in annotation:
                if a == 1:
                    colors.append(color_array[0])
                else:
                    colors.append('#FFFFFF')
        fig = plt.figure(figsize=(5, .75))
        ax2 = fig.add_axes([0.05, 0.8, 0.9, 0.15])
        cmap = mpl.colors.ListedColormap(colors)
        cmap.set_over('0.25')
        cmap.set_under('0.75')
        bounds = range(1, len(annotation)+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        annotation_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional',
                                                    orientation='horizontal')
        annotation_plot.set_label(annotation_names[i], fontsize=8)
        annotation_plot.set_ticks([])
        annotation_plot = fig
        annotation_tuple.append(annotation_plot)
    return annotation_tuple

def Assemble_Figure(zscore_plots, value_plots, heatmaps, annotation_plot, output, horizontal):
    """Assemble everything together and return svg and pdf of final figure"""
    DPI = 300
    size_prob_plot = 215
    size_stat_plot = 225
    size_annotation_plot = 30
    size_heatmap = 200

    if heatmaps == None:
        horizontal = 'n'
    elif len(heatmaps)>1:
        horizontal='y'

    if horizontal == 'y':
        size_width = "9in"
        size_height = '9in'
    else:
        size_width = "5in"
        size_height = '11in'

    fig = sg.SVGFigure(size_width, size_height)
    value_plots.savefig('value_plots.svg', format='svg', dpi=DPI, transparent=True)
    value_plots = sg.fromfile('value_plots.svg')
    plot1 = value_plots.getroot()
    if annotation_plot is not None:
        len_ann_plot = (len(annotation_plot))
    else:
        len_ann_plot = 0
    plot1.moveto(0, 0)
    fig.append(plot1)
    # plot heatmap(s)
    len_annotation_plot = size_annotation_plot * (len_ann_plot + 1)
    if heatmaps is not None:
        heatmap_count = 0
        for heatmap in heatmaps:
            plot4 = heatmap[0]
            try:
                plot4.savefig('heatmap.svg', format='svg', dpi=DPI, transparent=True)
                plot4 = sg.fromfile('heatmap.svg')
                plot4 = plot4.getroot()

                #transform and add heatmap figure; must be added first for correct layering
                if horizontal=='y':
                    y_scale = len_annotation_plot + size_prob_plot +size_heatmap*heatmap_count+ 75
                    plot4.moveto(375,y_scale, scale=1.40)
                    plot4.rotate(-45, 0, 0)
                    fig.append(plot4)
                    x_move = 510
            except Exception as e:
                print('generate heatmap error')
                print(e)

            else:
                y_scale = size_stat_plot*len(zscore_plots) + (size_heatmap+1) * heatmap_count + len_annotation_plot + size_prob_plot + 110
                plot4.moveto(0,y_scale, scale=1.40)
                plot4.rotate(-45, 0, 0)
                fig.append(plot4)
                x_move=110

            heatmap_count = heatmap_count + 1
        colorbar_h = heatmap[1]
        try:
            colorbar_h.savefig('colorbar_h.svg', format='svg', dpi=DPI, transparent=True)
            colorbar_h = sg.fromfile('colorbar_h.svg')
            colorbar_h = colorbar_h.getroot()
            colorbar_h.moveto(x_move, y_scale + size_heatmap)
            fig.append(colorbar_h)
        except Exception as e:
            print('generate colorbar_h error')
            print(e)


    if annotation_plot is not None:
        # transform and add annotations plots
        index = 0
        for plot in annotation_plot:
            plot.savefig('annotation_plot.svg', format='svg', dpi=DPI, transparent=True)
            plot = sg.fromfile('annotation_plot.svg')
            plot3 = plot.getroot()
            y_move = size_prob_plot + size_annotation_plot * (index + 1)
            plot3.moveto(30, y_move, scale=1.05)
            index += 1
            fig.append(plot3)

    #transform and add zscore plots
    index = 1

    for plot in zscore_plots:
        plot2 = plot[0]
        plot2.savefig('stats_plot.svg', format='svg', dpi=DPI, transparent=True)
        plot2 = sg.fromfile('stats_plot.svg')
        plot2 = plot2.getroot()
        y_move = size_stat_plot * index + len_annotation_plot
        index += 1
        plot2.moveto(0, y_move)
        fig.append(plot2)


    # extract colorbar
    y_move = size_stat_plot * len(zscore_plots) + len_annotation_plot + size_prob_plot
    plot = zscore_plots[0]
    colorbar = plot[1]
    colorbar.savefig('colorbar.svg', format='svg', dpi=DPI, transparent=True)
    colorbar = sg.fromfile('colorbar.svg')
    colorbar = colorbar.getroot()
    colorbar.moveto(100, y_move +40)
    fig.append(colorbar)

    #export final figure as a svg and pdf
    #svgfile = "canvis.svg"
    svgfile = output + '.svg'
    fig.save(svgfile)


    """ Uncomment if want to convert to PDF. Note: must have CarioSVG libraries installed

    pdffile = output + ".pdf"
    cairosvg.svg2pdf(url=svgfile, write_to=pdffile)

    """

    html_file = open("canvis.html",'w+')
    html_str = """
    <img src="canvis.svg" >

    """
    html_file.write(html_str)
    html_file.close()


def main():

    # Parse the command line data
    parser = OptionParser()
    parser.add_option("-l", "--locus_name", dest="locus_name")
    parser.add_option("-z", "--zscores", dest="zscores", action='callback', callback=vararg_callback)
    parser.add_option("-a", "--annotations", dest="annotations")
    parser.add_option("-s", "--specific_annotations", dest="specific_annotations", action='callback', callback=vararg_callback)
    parser.add_option("-r", "--ld_name", dest="ld_name", action='callback', callback=vararg_callback)
    parser.add_option("-t", "--threshold", dest="threshold", default=0)
    parser.add_option("-g", "--greyscale", dest="greyscale", default='n')
    parser.add_option("-o", "--output", dest="output", default='fig_final')
    parser.add_option("-i", "--interval", dest="interval", nargs=2)
    parser.add_option("-L", "--large_ld", dest="large_ld", default='n')
    parser.add_option("-H", "--horizontal", dest="horizontal", default='n')
    parser.add_option("-b", "--bp_head", default="pos", dest="bp_head")

    # extract options
    (options, args) = parser.parse_args()
    locus_name = options.locus_name
    zscore_names = options.zscores
    ld_name = options.ld_name
    annotations = options.annotations
    annotation_names = options.specific_annotations
    threshold = options.threshold
    threshold = int(threshold)
    if threshold < 0 or threshold > 100:
        warnings.warn('Specified threshold is not valid; threshold is set to 0')
        threshold = 0
    else:
        threshold = (threshold)*.01
    greyscale = options.greyscale
    output = options.output
    interval = options.interval
    large_ld = options.large_ld
    horizontal = options.horizontal
    bp_head = options.bp_head

    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position)
        --zscores [-z] specific zscores to be plotted
        --annotations [-a]  specify annotation file name
        --specific_annotations [-s] annotations to be plotted
        --ld_name [r] specify the ld_matrix file name
        --threshold [-t] threshold for credible set [default: 0]
        --greyscale [-g] sets colorscheme to greyscale [default: n]
        --output [-o] desired name of output file
        --interval [-i] designated interval [default: all locations]
        --large_ld [-L] overrides to produce large LD despite large size
        """

    #check if required flags are presnt
    if(locus_name == None or zscore_names == None):
        sys.exit(usage)

    [zscores, pos_prob, location, ld, annotations, annotation_names] = Read_Input(locus_name, zscore_names,
                                                                                  ld_name, annotations, annotation_names, interval, bp_head)
    zscore_plots = Plot_Statistic_Value(location, zscores, zscore_names, greyscale, ld)
    value_plots = Plot_Position_Value(location, pos_prob, threshold, greyscale)

    if ld is not None:
        heatmap = Plot_Heatmap(ld, greyscale, large_ld)
        
    else:
        heatmap = None

    if annotations is not None:
        annotation_plot = Plot_Annotations(annotation_names, annotations, greyscale)
    else:
        annotation_plot = None

    Assemble_Figure(zscore_plots, value_plots, heatmap, annotation_plot, output, horizontal)

    #remove extraneous files
    if heatmap is not None:
        try:
            os.remove('heatmap.svg')
            os.remove('colorbar_h.svg')
        except:
            pass
    os.remove('stats_plot.svg')
    if annotation_plot is not None:
        os.remove('annotation_plot.svg')
    #os.remove("canvis.svg")
    os.remove('value_plots.svg')

if __name__ == "__main__":
    main()
