import matplotlib.pyplot as plt
import numpy as np
import glymur
import os
import glob
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import imageio

x = glob.glob('**/*.jp2', recursive=True)
x = sorted(x)
nbins = 7
img = np.zeros([7, 1024, 1024])
sub_titles=['LogT = 5.7-6.1','LogT = 6.1-6.3', 'LogT = 6.3-6.5', 'LogT = 6.5-6.7', 'LogT = 6.7-6.9', 'LogT = 6.9-7.1']

for j in range(0, int(len(x) / 7)):
    fig = plt.figure(figsize=(12.5, 16))
    fig.suptitle(x[j*7][-52:-48]+'/'+x[j*7][-47:-45]+'/'+x[j*7][-44:-42]+' : '+x[j*7][-40:-38]+'H')
    for i in range(0, 7):
        jp2 = glymur.Jp2k(x[j * 7 + i])
        img[i, :, :] = jp2[:]
    img[1, :, :] = (img[0, :, :] + img[1, :, :]) / 2
    for i in range(1, 7):
        ax=fig.add_subplot(3, 2, i)
        plt.imshow(img[i, :, :], vmin=0, vmax=255, cmap='inferno')
        ax.set_title(sub_titles[i-1])
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cmap = cm.ScalarMappable(
        norm=mcolors.Normalize(0, 255),
        cmap=plt.get_cmap('inferno'))
    cmap.set_array([])
    fig.colorbar(cmap, cax=cbar_ax)
    plt.savefig(str(j).zfill(4) + '.png')
    plt.close()

filenames = sorted(glob.glob('*.png'))
with imageio.get_writer('movie.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
