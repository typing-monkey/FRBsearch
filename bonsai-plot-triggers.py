#!/usr/bin/env python
#
# A quick-and-dirty plot script for the triggers.hdf5 file.
# Currently very primitive, feel free to improve it!


import os
import sys

if len(sys.argv) != 2:
    print >>sys.stderr, 'usage: bonsai-plot-triggers.py <infile.hdf5>'
    sys.exit(2)

input_filename = sys.argv[1]

if not os.path.exists(input_filename):
    print >>sys.stderr, '%s: file not found' % input_filename
    sys.exit(1)

if input_filename.endswith('.h5'):
    input_stem = input_filename[:-3]
elif input_filename.endswith('.hdf5'):
    input_stem = input_filename[:-5]
else:
    print >>sys.stderr, 'bonsai-plot-triggers.py: this script is simpleminded and requires the input filename to end in .h5 or .hdf5'
    sys.exit(1)


####################################################################################################


import PIL.Image
import numpy as np


def write_png(filename, arr, weights=None, transpose=False, ytop_to_bottom=False):
    """This is a quick-and-dirty plotting routine that I cut-and-paste everywhere."""

    arr = np.array(arr, dtype=np.float)
    assert arr.ndim == 2

    if weights is None:
        weights = np.ones(arr.shape, dtype=np.float)
    else:
        weights = np.array(weights, dtype=np.float)
        assert weights.shape == arr.shape
    
    if not transpose:
        arr = np.transpose(arr)
        weights = np.transpose(weights) if (weights is not None) else None
    
    if not ytop_to_bottom:
        arr = arr[::-1]
        weights = weights[::-1] if (weights is not None) else None

    (wmin, wmax) = (np.min(weights), np.max(weights))
    if wmin < 0:
        raise RuntimeError('write_png: negative weights are currently treated as an error')

    # A corner case..
    if wmax == 0.0:
        print >>sys.stderr, '%s: array was completely masked, writing all-black image' % filename
        rgb = np.zeros((arr.shape[0], arr.shape[1], 3), dtype=np.uint8)
        img = PIL.Image.fromarray(rgb)
        img.save(filename)
        return

    # normalize weights to [0,1]
    weights = weights/wmax

    # weighted mean and rms
    mean = np.sum(weights*arr) / np.sum(weights)
    var = np.sum((weights*(arr-mean))**2) / np.sum(weights**2)
    rms = np.sqrt(var) if (var > 0.0) else 1.0    # The "1.0" handles the corner case of a constant array.

    # color in range [0,1].
    color = 0.5 + 0.16*(arr-mean)/rms    # factor 0.16 preserves convention from some old code
    color = np.maximum(color, 0.0)
    color = np.minimum(color, 0.999999)  # 0.99999 instead of 1.0, to make roundoff-robust
    
    # rgb in range [0,1]
    red = 256. * color * weights
    blue = 256. * (1-color) * weights

    rgb = np.zeros((arr.shape[0],arr.shape[1],3), dtype=np.uint8)
    rgb[:,:,0] = red
    rgb[:,:,2] = blue

    img = PIL.Image.fromarray(rgb)
    img.save(filename)
    print >>sys.stderr, 'wrote %s' % filename


####################################################################################################


import h5py

f = h5py.File(input_filename, 'r')
ntrees = f.attrs['NTREES']
warnflag = False

for itree in xrange(ntrees):
    trigger_dataset = f['TREE%d' % itree]['TRIGGERS']
    (ndm, nsm, nbeta, nt) = trigger_dataset.shape
    
    is_monster_plot = (ndm > 1024) or (nt > 2048)

    if is_monster_plot and not warnflag:
        print >>sys.stderr, "Warning: bonsai-plot-triggers.py is very primitive, and splitting large files into multiple plots isn't implemented yet"
        warnflag = True

    for ism in xrange(nsm):
        for ibeta in xrange(nbeta):
            id_str = ''
            if ntrees > 1:
                id_str += ('_tree%d' % itree)
            if nsm > 1:
                id_str += ('_sm%d' % ism)
            if nbeta > 1:
                id_str += ('_beta%d' % ibeta)

            if is_monster_plot:
                print >>sys.stderr, 'Attempting to write %d-by-%d monster plot' % (ndm, nt)

            output_filename = '%s%s.png' % (input_stem, id_str)
            write_png(output_filename, trigger_dataset[:,ism,ibeta,:], transpose=True)
