#!/usr/bin/env pvbatch

# state file generated using paraview version 5.8.0

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

import argparse
import numpy as np

import paratools

parser = argparse.ArgumentParser(description="Renders volume fraction chi.")
parser.add_argument('files', nargs='+', help="list of data files 'data_*.vtk'")
parser.add_argument('--force',
                    action="store_true",
                    help="Overwrite existing files")
parser.add_argument('--resy',
                    type=int,
                    default=900,
                    help="Output image resolution")
args = parser.parse_args()

sources_ft = []
timearrays = []

files_chi = paratools.ReplaceFilename(args.files, "data_{}.xmf")
source_chi = XDMFReader(FileNames=files_chi)
source_chi.CellArrayStatus = ['chi']
(source_chi, ), (timearray, ) = paratools.ApplyForceTime([source_chi])
sources_ft.append(source_chi)
timearrays.append(timearray)

renderView1 = CreateView('RenderView')
renderView1.OrientationAxesVisibility = 0
renderView1.UseLight = 0
box = paratools.GetBoundingBox(source_chi)
boxsize = box[1] - box[0]
boxc = (box[0] + box[1]) * 0.5
viewsize = np.round(boxsize[:2] * args.resy / boxsize[1]).astype(int)
renderView1.ViewSize = viewsize
renderView1.CameraPosition = [boxc[0], boxc[1], 100]
renderView1.CameraFocalPoint = [boxc[0], boxc[1], 0]
renderView1.CameraParallelScale = boxsize[1] * 0.5
renderView1.CameraParallelProjection = 1
renderView1.Background = [1] * 3

celltopoint_chi = CellDatatoPointData(Input=source_chi)
celltopoint_chi.CellDataArraytoprocess = ['chi']

chiLUT = GetColorTransferFunction('chi')
chiLUT.RGBPoints = [
    0.0, 0.403922, 0.0, 0.121569, 0.0588252296408856, 0.576932, 0.055363,
    0.14925, 0.12157207458920345, 0.728489, 0.155017, 0.197386,
    0.1843189195375213, 0.81707, 0.33218, 0.281046, 0.24706626450054098,
    0.894579, 0.503806, 0.399769, 0.3098131094488588, 0.960323, 0.66782,
    0.536332, 0.37255995439717665, 0.982468, 0.800692, 0.706113,
    0.43530679934549443, 0.994925, 0.908651, 0.857901, 0.4980538593001341,
    0.999846, 0.997232, 0.995694, 0.5608009892568321, 0.926105, 0.926105,
    0.926105, 0.6235478342051499, 0.843368, 0.843368, 0.843368,
    0.6862946791534676, 0.749865, 0.749865, 0.749865, 0.7490415241017854,
    0.631373, 0.631373, 0.631373, 0.8117883690501032, 0.502653, 0.502653,
    0.502653, 0.874535714013123, 0.359939, 0.359939, 0.359939,
    0.9372825589614407, 0.227451, 0.227451, 0.227451, 1.0000294039097586,
    0.101961, 0.101961, 0.101961
]
chiLUT.Discretize = 0
chiLUT.ColorSpace = 'Lab'
chiLUT.ScalarRangeInitialized = 1.0

celltopoint_chiDisplay = Show(celltopoint_chi, renderView1)
celltopoint_chiDisplay.Representation = 'Surface'
celltopoint_chiDisplay.ColorArrayName = ['POINTS', 'chi']
celltopoint_chiDisplay.LookupTable = chiLUT

steps = paratools.GetSteps(args.files)
paratools.SaveAnimation(steps,
                        renderView1,
                        sources_ft,
                        timearrays,
                        force=args.force)
