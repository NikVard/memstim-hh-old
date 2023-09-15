# -*- coding: utf-8 -*-

# OS stuff
import os
import sys
import warnings
from pathlib import Path
from tqdm import tqdm

# Computational stuff
import numpy as np
import matplotlib as mplb
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib import ticker
import matplotlib.patheffects as patheffects

# TensorPAC
from tensorpac import Pac
from tensorpac.utils import PSD

# Other scripts and my stuff
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = Path(script_dir).parent
sys.path.insert(0, os.path.abspath(parent_dir))

import parameters
from src.freq_analysis import *
from src.figure_plots_parameters import *

# ILLUSTRATOR STUFF
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.titlesize'] = fsize_titles
plt.rcParams['axes.labelsize'] = fsize_xylabels

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
})

# Arial font everywhere
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'

# Sizebar func
def add_sizebar(ax, xlocs, ylocs, bcolor, text, textx, texty, fsize, rot, ha, va):
    """ Add a sizebar to the provided axis """
    ax.plot(xlocs, ylocs, ls='-', c=bcolor, linewidth=1., rasterized=False, clip_on=False)

    # add the text
    if type(text) is list:
        # bottom text
        ax.text(x=textx[0], y=texty[0], s=text[0], rotation=rot[0], fontsize=fsize, va=va, ha=ha, clip_on=False)

        # top text
        ax.text(x=textx[1], y=texty[1], s=text[1], rotation=rot[1], fontsize=fsize, va=va, ha=ha, clip_on=False)
    else:
        ax.text(x=textx, y=texty, s=text, rotation=rot, fontsize=fsize, va=va, ha=ha, clip_on=False)


# get the root directory for the simulations
root_cluster = os.path.join(parent_dir, 'results', 'fig4_kN2')
# print(root_cluster)

# get a list of the directories with oscillator amplitudes
osc_amplitude_dirs = sorted(next(os.walk(root_cluster))[1])
# print(osc_amplitude_dirs)

# set the directory with the precomputed results
root_precomp = os.path.join(parent_dir, 'figures', 'fig4_kN2', 'data')

# Timing parameters
second = 1
ms = 1e-3
duration = 10*second
stim_onset = 1800*ms
dt = 0.1*ms
fs = int(1*second/dt)
tv = np.linspace(0, duration, fs*duration)/ms
winsize_FR = 5*ms
overlap_FR = 0.9
winstep_FR = winsize_FR*round(1-overlap_FR,4)
fs_FR = int(1/winstep_FR)

# Area names and sizes
areas = [['EC_pyCAN', 'EC_inh'], ['DG_py', 'DG_inh'], ['CA3_pyCAN', 'CA3_inh'], ['CA1_pyCAN', 'CA1_inh']]
area_labels = ['EC', 'DG', 'CA3', 'CA1']
N_tot = [[10000, 1000], [10000, 100], [1000, 100], [10000, 1000]]

# Raster downsampling
N_scaling = 100
N_gap = 10

# Color selection
c_inh = '#bf616a'
c_exc = '#5e81ac'

# Firing rates plotting gap
rates_gap = 250 # Hz

# get the oscillator amplitude values
# osc_amps = np.arange(0.01, 0.23, 0.01)
osc_amps = []
for val in osc_amplitude_dirs:
    osc_amps.append(float(val.split('_')[1]))
osc_amps.sort()
osc_amps = np.array(osc_amps)

# kN values = np.arange(1, 11, 1)
kN_vals = np.arange(0.,55.,5.)

# make a meshgrid for plotting heatmaps
X, Y = np.meshgrid(osc_amps, kN_vals)
MI_heatmap_E = np.load(os.path.join(root_precomp, 'MI_E_heatmap.npy'))
MI_heatmap_I = np.load(os.path.join(root_precomp, 'MI_I_heatmap.npy'))

# Define points of interest % list of tuples: (label, [osc, stim] nA)
points_of_interest = []
points_of_interest.append((0.20, 5.0)) # A - no activity % post-stim
# points_of_interest.append((0.20, 15.0)) # B - transient activity
points_of_interest.append((0.20, 35.0)) # C - restoration of physiological activity
points_of_interest_labels = ['A', 'B']

# theta / gamma power matrices
theta_heatmap_E = np.load(os.path.join(root_precomp, 'theta_E_heatmap.npy'))
theta_heatmap_I = np.load(os.path.join(root_precomp, 'theta_I_heatmap.npy'))
gamma_heatmap_E = np.load(os.path.join(root_precomp, 'gamma_E_heatmap.npy'))
gamma_heatmap_I = np.load(os.path.join(root_precomp, 'gamma_I_heatmap.npy'))

# ILLUSTRATOR STUFF
mplb.rcParams['pdf.fonttype'] = 42
mplb.rcParams['ps.fonttype'] = 42

# Make figure outline
fig_width = 7.5
fig_height = 8.75
fig = plt.figure(figsize=(fig_width, fig_height), tight_layout=False)

# Add gridspecs
gs_outer = fig.add_gridspec(1, 2) # if you add wspace=<%> remove the tight_layout!!
gs_left = GridSpecFromSubplotSpec(len(points_of_interest_labels), 1, hspace=0.2, subplot_spec=gs_outer[0])
gs_right = GridSpecFromSubplotSpec(3, 1, hspace=0.25, subplot_spec=gs_outer[1])

# Organize the axes
axs_left = []
axs_right = []
axs_rasters = []
axs_FRs = []

for elem in gs_right:
    ax = fig.add_subplot(elem)
    axs_right.append(ax)

# Parameters (backup) filename
fname = 'parameters_bak.json'
# t_stim = data['stimulation']['onset']*1000
t_stim = 1800.0 # ms

# iterate over oscillator amplitudes
axs_idxs = []
for osc_amp_dir in tqdm(osc_amplitude_dirs, desc='Plotting [A] w/ points of interest'):
    curr_osc_amp = float(osc_amp_dir.split('_')[1])
    idx_osc_amp = np.where(osc_amps == curr_osc_amp) # index for heatmap

    # one step in
    curr_osc_dir = os.path.join(root_cluster, osc_amp_dir)
    print('osc: ', osc_amp_dir)

    # two steps in - stim 7.0nA @ 1800ms
    for kN_dir in next(os.walk(os.path.join(curr_osc_dir, 'data', '7.0_nA', '0.00_1800.0_ms')))[1]:

        # one more step in - kN
        curr_kN_dir = os.path.join(curr_osc_dir, 'data', '7.0_nA', '0.00_1800.0_ms', kN_dir)

        # load parameters backup data
        try:
            data = parameters.load(os.path.join(curr_kN_dir, fname))
        except Exception as e:
            print('[!]' + "Error code " + str(e.errno) + ": " + e.strerror + ' | Parameters file not found!')
            continue;

        curr_kN_val = data['Kuramoto']['kN']
        idx_kN_val = np.where(kN_vals == curr_kN_val)[0] # index for heatmap

        # Set the current vals #-> used to decide if we're plotting traces
        curr_vals = (curr_osc_amp, curr_kN_val)

        # Check if we need to plot the traces
        if curr_vals in points_of_interest:
            print(curr_vals)
            # where is it found, what's the label of the point
            label_idx = points_of_interest.index(curr_vals)
            print("[>] Found label ", points_of_interest_labels[label_idx])

            # load the data (spikes) for the CA1 E-group
            curr_data_dir = os.path.join(curr_kN_dir, 'data')
            curr_spikes_dir = os.path.join(curr_data_dir, 'spikes')
            # print('[L5]-----data: ', curr_data_dir)
            # print('[L6]------spikes: ', curr_spikes_dir)

            print(curr_data_dir)

            # load the rhythm
            rhythm = np.loadtxt(os.path.join(curr_data_dir, 'order_param_mon_rhythm.txt'))

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning, append=1)

                if len([entry for entry in os.listdir(curr_spikes_dir) if os.path.isfile(os.path.join(curr_spikes_dir, entry))]) == 16:
                    i_exc = np.loadtxt(os.path.join(curr_spikes_dir, 'CA1_pyCAN_spikemon_i.txt'))
                    t_exc = np.loadtxt(os.path.join(curr_spikes_dir, 'CA1_pyCAN_spikemon_t.txt'))
                    i_inh = np.loadtxt(os.path.join(curr_spikes_dir, 'CA1_inh_spikemon_i.txt'))
                    t_inh = np.loadtxt(os.path.join(curr_spikes_dir, 'CA1_inh_spikemon_t.txt'))
                else:
                    print("[!] Warning: files missing!", "osc_amp: ", curr_osc_amp, " kN: ", curr_kN_val, " dir: ", curr_data_dir)
                    continue

            # Fix the timings -> from ms to sec
            i_exc = i_exc.astype(int)
            t_exc = t_exc*ms
            i_inh = i_inh.astype(int)
            t_inh = t_inh*ms

            # sort based on index number (lower to higher)
            idx_sort_exc = np.argsort(i_exc)
            idx_sort_inh = np.argsort(i_inh)
            i_exc = i_exc[idx_sort_exc]
            t_exc = t_exc[idx_sort_exc]
            i_inh = i_inh[idx_sort_inh]
            t_inh = t_inh[idx_sort_inh]

            # set number of neurons {CA1}
            N_exc = N_tot[3][0]
            N_inh = N_tot[3][1]

            # select some neurons randomly, subscaling
            # exc_mixed = np.random.permutation(np.arange(N_exc))[:N_scaling]
            # inh_mixed = np.random.permutation(np.arange(N_inh))[:N_scaling]
            exc_mixed = np.arange(0, N_exc+1, int(N_exc/N_scaling))
            inh_mixed = np.arange(0, N_inh+1, int(N_inh/N_scaling))

            idx_exc = np.in1d(i_exc, exc_mixed)
            idx_inh = np.in1d(i_inh, inh_mixed)

            i_exc_sub = i_exc[idx_exc]
            t_exc_sub = t_exc[idx_exc]
            i_inh_sub = i_inh[idx_inh]
            t_inh_sub = t_inh[idx_inh]

            # assign new neuron count numbers
            cnt = 0
            i_exc_sub_new = np.copy(i_exc_sub)
            for ii in exc_mixed:
                idx_tmp = np.where(i_exc_sub == ii)
                # print('changing ', ii, 'to ', cnt)
                i_exc_sub_new[idx_tmp] = cnt
                cnt += 1
            i_exc_sub = i_exc_sub_new

            # cnt = 0
            cnt += N_gap
            i_inh_sub_new = np.copy(i_inh_sub)
            for ii in inh_mixed:
                idx_tmp = np.where(i_inh_sub == ii)
                # print('changing ', ii, 'to ', cnt)
                i_inh_sub_new[idx_tmp] = cnt
                cnt += 1
            i_inh_sub = i_inh_sub_new

            # get the firing rate (spikes-to-rates)
            tv_inh_FR, FR_inh, fs_FR2 = my_FR(spikes=t_inh, duration=duration, window_size=winsize_FR, overlap=overlap_FR)
            tv_exc_FR, FR_exc, _ = my_FR(spikes=t_exc, duration=duration, window_size=winsize_FR, overlap=overlap_FR)

            # Normalize the FRs
            FR_inh_norm = (FR_inh/winsize_FR)/N_tot[3][1]
            FR_exc_norm = (FR_exc/winsize_FR)/N_tot[3][0]

            # avoid divisions by zero
            # FR_inh_norm[FR_inh_norm<=1.e-5] = 1.e-5
            # FR_exc_norm[FR_exc_norm<=1.e-5] = 1.e-5

            # subplots rasters | FRs
            gs_curr = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_left[label_idx], hspace=0)
            # ax_rhythm = fig.add_subplot(gs_left[label_idx])
            # ax_rhythm.plot(tv*ms, rhythm)

            # plot the rasters for E-I
            ax_rasters_curr = fig.add_subplot(gs_curr[0])
            ax_rasters_curr.scatter(t_inh_sub, i_inh_sub, s=1.25, linewidth=1., marker='.', c=c_inh, edgecolors='none', alpha=1., rasterized=True)
            ax_rasters_curr.scatter(t_exc_sub, i_exc_sub, s=1.25, linewidth=1., marker='.', c=c_exc, edgecolors='none', alpha=1., rasterized=True)

            # plot the FRs for E-I
            ax_FRs_curr = fig.add_subplot(gs_curr[1])
            # ax_FRs_curr.plot(tv_inh_FR, FR_inh_norm+200)
            # ax_FRs_curr.plot(tv_exc_FR, FR_exc_norm, c=)
            ax_FRs_curr.plot(tv_inh_FR, FR_inh_norm+rates_gap, ls='-', linewidth=1., c=c_inh, label='inh', zorder=10, rasterized=False)
            ax_FRs_curr.plot(tv_exc_FR, FR_exc_norm, ls='-', linewidth=1., c=c_exc, label='exc', zorder=10, rasterized=False)

            # ax_rasters_curr.text(x=0.5, y=0.5, s=points_of_interest_labels[label_idx], transform=ax_rasters_curr.transAxes, fontsize=12)

            # add the axes to the lists
            axs_idxs.append(label_idx)
            axs_rasters.append(ax_rasters_curr)
            axs_FRs.append(ax_FRs_curr)


# Done with the iterations
print("Done.")

# Share x axes in FRs
for ax, label in zip(axs_rasters, [points_of_interest_labels[idx] for idx in axs_idxs]):
    # Hide some spines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Text to mark box
    # ax.text(x=-0.02, y=0.85, transform=ax.transAxes, weight='bold', s=label, fontsize=fsize_misc, ha='center', color='red', bbox=dict(boxstyle='square', edgecolor='red', facecolor='none'), clip_on=False)

    # Set xlims
    ax.set_xlim([(t_stim-500)*ms, (t_stim+2500)*ms])
    ax.set_ylim([0, cnt+1])

    # Set y-labels
    ax.set_ylabel('CA1 Neurons', fontsize=fsize_xylabels, labelpad=20)

    # ax.set_title(label)

    # Set the ticks
    # ax.xaxis.set_major_locator(ticker.FixedLocator(ax_rate_majors))
    # ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    # Add the title
    ax.set_title(label+'.', loc='left', weight='bold', fontsize=fsize_panels)


ax_rate_majors = np.arange(0., duration, .5) #[0.5, 1.0, 1.5...]
for ax in axs_FRs:
    # Hide some spines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set xticks
    ax.xaxis.set_major_locator(ticker.FixedLocator(ax_rate_majors))

    # Set xlims
    ax.set_xlim([(t_stim-500)*ms, (t_stim+2500)*ms])
    ax.set_ylim([-1, 600])

    # Set y-labels
    ax.set_ylabel('CA1 Firing Rates', fontsize=fsize_xylabels, labelpad=20)

    # Set the ticks
    # ax.xaxis.set_major_locator(ticker.FixedLocator(ax_rate_majors))
    # ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

# Last axis shows time
# axs_FRs[-1].spines['bottom'].set_visible(True)

# Set xlabel
# axs_FRs[-1].set_xlabel('Time [s]', fontsize=fsize_xylabels)

# Add sizebars for x-y axes
xlims_sz = [t_stim*ms-0.550, t_stim*ms-0.300]
ylims_sz = [-75, 25]
add_sizebar(axs_FRs[np.argmax(axs_idxs)], xlims_sz, [ylims_sz[0]]*2, 'black', '250 ms', fsize=fsize_xylabels, rot=0, textx=np.mean(xlims_sz), texty=ylims_sz[0]-50, ha='center', va='top')
add_sizebar(axs_FRs[np.argmax(axs_idxs)], [xlims_sz[0]]*2, ylims_sz, 'black', '100 Hz', fsize=fsize_xylabels, rot=90, textx=xlims_sz[0]-0.05, texty=np.mean(ylims_sz), ha='right', va='center')

# Sort the axes from top to bottom
axs_idxs_sorted_I = np.argsort(axs_idxs)
axs_rasters_sorted = np.array(axs_rasters)[axs_idxs_sorted_I]
axs_FRs_sorted = np.array(axs_FRs)[axs_idxs_sorted_I]

# Stimulation lines + marker
# First raster
axs_rasters_sorted[0].axvline(x=t_stim*ms, ymin=-0.1, ymax=1.0, color='gray', alpha=0.75, ls='-', linewidth=0.75, zorder=10, rasterized=False, clip_on=True)

# Rasters + FRs
for ax in axs_rasters_sorted:
    ax.axvline(x=t_stim*ms, ymin=-0.1, ymax=1.1, color='gray', alpha=0.75, ls='-', linewidth=0.75, zorder=10, rasterized=False, clip_on=False)

for ax in axs_FRs_sorted:
    ax.axvline(x=t_stim*ms, ymin=-0.1, ymax=1.1, color='gray', alpha=0.75, ls='-', linewidth=0.75, zorder=10, rasterized=False, clip_on=False)

# Last FR
axs_FRs_sorted[-1].axvline(x=t_stim*ms, ymin=0, ymax=1.1, color='gray', alpha=0.75, ls='-', linewidth=0.75, zorder=10, rasterized=False, clip_on=True)

# Don't forget to mark the stimulation onset!
# add marker for stimulation onset
axs_rasters[np.argmin(axs_idxs)].scatter(x=t_stim*ms, y=cnt+30, s=15, linewidth=1., marker='v', c='gray', edgecolors=None, alpha=1, rasterized=False, clip_on=False)

# set vmin/vmax for plotting
# vmin = 1e-12
# vmax = 1.
vmin_MI = min(MI_heatmap_I.min(), MI_heatmap_E.min())
vmax_MI = max(MI_heatmap_I.max(), MI_heatmap_E.max())
vmin_theta = min(theta_heatmap_I.min(), theta_heatmap_E.min())
vmax_theta = max(theta_heatmap_I.max(), theta_heatmap_E.max())
vmin_gamma = min(gamma_heatmap_I.min(), gamma_heatmap_E.min())
vmax_gamma = max(gamma_heatmap_I.max(), gamma_heatmap_E.max())

# normalize the colors
norm_MI = mplb.colors.Normalize(vmin=vmin_MI, vmax=vmax_MI)
norm_theta = mplb.colors.Normalize(vmin=vmin_theta, vmax=vmax_theta)
norm_gamma = mplb.colors.Normalize(vmin=vmin_gamma, vmax=vmax_gamma)

# Plot the heatmaps
# im1 = axs_right[0].contourf(X+X.min()/2, Y+Y.min()/2, MI_heatmap_I.T, levels=5, norm=norm_MI)
# im2 = axs_right[0].contourf(X+X.min()/2, Y+Y.min()/2, MI_heatmap_E.T, levels=5, norm=norm_MI)
im2 = axs_right[2].pcolormesh(X+X.min()/2, Y+Y.min()/2, MI_heatmap_E.T, shading='auto', rasterized=True)
# cs2 = axs_right[2].contour(X+X.min()/2, Y+Y.min()/2, MI_heatmap_E.T, levels=[0.15*MI_heatmap_E.max()], colors='white', origin='lower')

# Theta power
# im3 = axs_theta[0].contourf(X+X.min()/2, Y+Y.min()/2, theta_heatmap_I.T, levels=5, norm=norm_theta)
# im4 = axs_right[1].contourf(X+X.min()/2, Y+Y.min()/2, theta_heatmap_E.T, levels=5, norm=norm_theta)
im4 = axs_right[0].pcolormesh(X+X.min()/2, Y+Y.min()/2, theta_heatmap_E.T, cmap='inferno', shading='auto', rasterized=True)
# cs4 = axs_right[0].contour(X+X.min()/2, Y+Y.min()/2, theta_heatmap_E.T, levels=[0.15*theta_heatmap_E.max()], colors='white', origin='lower')

# Gamma power
# im5 = axs_right[2].contourf(X+X.min()/2, Y+Y.min()/2, gamma_heatmap_I.T, levels=5, norm=norm_gamma)
# im6 = axs_right[2].contourf(X+X.min()/2, Y+Y.min()/2, gamma_heatmap_E.T, levels=5, norm=norm_gamma)
im6 = axs_right[1].pcolormesh(X+X.min()/2, Y+Y.min()/2, gamma_heatmap_E.T, cmap='inferno', shading='auto', rasterized=True)
# cs6 = axs_right[1].contour(X+X.min()/2, Y+Y.min()/2, gamma_heatmap_E.T, levels=[0.15*gamma_heatmap_E.max()], colors='white')

# Mark points of interest using rectangles
rect_list = []

# Add text annotations
text_c = [(points_of_interest[0][0]-0.05, points_of_interest[0][1]-0.0),
          (points_of_interest[1][0]-0.05, points_of_interest[1][1]-0.0),
          # (points_of_interest[2][0]-0.05, points_of_interest[2][1]-0.0)
         ]
for point, label, text_coords in zip(points_of_interest, points_of_interest_labels, text_c):
    # Add annotations
    x1,y1 = point; x2,y2 = text_coords; ratio = 0.85
    xn = ratio*x1 + (1-ratio)*x2
    yn = ratio*y1 + (1-ratio)*y2

    # Point of interest
    axs_right[0].plot([x1], [y1], marker='o', mec='black', mfc='white',  markersize=5)

    # Annotation
    PE = [patheffects.withStroke(linewidth=2, foreground='black', capstyle="round")]
    anno = axs_right[0].annotate(label,
                          xy=(xn, yn),
                          xytext=text_coords,
                          arrowprops=dict(arrowstyle='-', path_effects=PE, edgecolor='white', linewidth=2),
                          path_effects=PE,
                          # bbox=dict(boxstyle='square', edgecolor='white', facecolor='none'),
                          weight='bold', fontsize=fsize_titles, color='white',
                          horizontalalignment='center',
                          verticalalignment='center')
    # Line adjustment
    anno.arrow_patch.set_path_effects([
        patheffects.Stroke(linewidth=3, foreground="k"),
        patheffects.Normal()])

# axs_right[0].annotate(label,
#                       xy=text_coords,
#                       xytext=text_coords,
#                       bbox=dict(boxstyle='square', edgecolor='white', facecolor='none'),
#                       weight='bold', fontsize=fsize_misc, color='white',
#                       horizontalalignment='center',
#                       verticalalignment='top')
# axs_right[0].annotate("",
#                       xy=(xn, yn),
#                       xytext=text_coords,
#                       arrowprops=dict(facecolor='none', edgecolor='white', width=3, headlength=6, headwidth=6),
#                       bbox=dict(boxstyle='square', edgecolor='white', facecolor='none'),
#                       weight='bold', fontsize=fsize_misc, color='white',
#                       horizontalalignment='center',
#                       verticalalignment='top')
# axs_right[0].text(x=xn, y=yn-0.5, s='(x={x:.2f},y={y:d})'.format(x=x1, y=int(y1)), fontsize=fsize_legends, ha='left', va='top', color='red')

# # Place the rectangles
# for rect, textloc, textlabel in zip(rect_list, text_c, points_of_interest_labels):
#     axs_right[0].add_patch(rect)
#     axs_right[0].annotate(textlabel, textloc,
#                         color='red', weight='bold', fontsize=fsize_legends,
#                         ha='center', va='center')

# Set xlims of images
for ax in axs_right[0:]:
    ax.set_xlim((0., 0.21))

# Add points, mark other axes
for ax in axs_right[1:]:
    for point in points_of_interest:
        # Point of interest
        x1,y1 = point
        ax.plot([x1], [y1], marker='o', mec='black', mfc='white',  markersize=5)

# Set x-y ticks properly
for ax in axs_right:
    ax.set_xticks([])
    ax.set_yticks(kN_vals[1::2])
    ax.tick_params(axis='both', which='both', labelsize=fsize_ticks)
ax.set_xticks(osc_amps[1::4])

# Set x-y labels
axs_right[0].set_ylabel('Synchronization (k/N)', fontsize=fsize_xylabels)
axs_right[1].set_ylabel('Synchronization (k/N)', fontsize=fsize_xylabels)
axs_right[2].set_ylabel('Synchronization (k/N)', fontsize=fsize_xylabels)
axs_right[2].set_xlabel('Osc. amp. [nA]', fontsize=fsize_xylabels)

# Set titles
axs_right[0].set_title('Theta band power', fontsize=fsize_titles)
axs_right[1].set_title('Gamma band power', fontsize=fsize_titles)
axs_right[2].set_title('Phase-Amplitude Coupling', fontsize=fsize_titles)

# Add the colorbars
# plt.colorbar(im2, label="Modulation Index", ax=axs_right[0])
# plt.colorbar(im4, label="Theta Power", ax=axs_right[1])
# plt.colorbar(im6, label="Gamma Power", ax=axs_right[2])
cbar_theta = plt.colorbar(im4, ax=axs_right[0])
cbar_gamma = plt.colorbar(im6, ax=axs_right[1])
cbar_MI = plt.colorbar(im2, ax=axs_right[2])

cbar_theta.ax.set_title('a.u.', fontsize=fsize_xylabels)
cbar_gamma.ax.set_title('a.u.', fontsize=fsize_xylabels)
cbar_MI.ax.set_title('a.u.', fontsize=fsize_xylabels)

# Add contour lines
# cbar_MI.add_lines(cs2)
# cbar_theta.add_lines(cs4)
# cbar_gamma.add_lines(cs6)

# Tick parameters
cbar_MI.ax.tick_params(labelsize=fsize_ticks)
cbar_theta.ax.tick_params(labelsize=fsize_ticks)
cbar_gamma.ax.tick_params(labelsize=fsize_ticks)

# Set the panel labels
# axs_rasters[0].set_title('A.', loc='left', weight='bold', fontsize=fsize_panels)
# axs_rasters[1].set_title('B.', loc='left', weight='bold', fontsize=fsize_panels)
# axs_rasters[2].set_title('C.', loc='left', weight='bold', fontsize=fsize_panels)
axs_right[0].set_title('C.', loc='left', weight='bold', fontsize=fsize_panels)

# Try the Gridspec tight_layout approach
gs_outer.tight_layout(fig)

# Save the figure
print('[+] Saving the figure...')
fig.savefig(os.path.join(parent_dir, 'figures', 'fig5', 'fig5_Supp4_kN_full.png'), transparent=True, dpi=300, format='png', bbox_inches='tight')
fig.savefig(os.path.join(parent_dir, 'figures', 'fig5', 'fig5_Supp4_kN_full.pdf'), transparent=True, dpi=300, format='pdf', bbox_inches='tight')

# Show the images
plt.show()

# Exit regularly
sys.exit(0)
