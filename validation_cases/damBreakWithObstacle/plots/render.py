import sys
import os.path
import numpy as np
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def screenshot(fr, out, res):
  for i in fr:
    cnt=int(round(i*50))
    fname = 'plots/frames/'+str(out)+'.'+str('{:03d}'.format(cnt))+'.png'
    if os.path.exists(fname) == False:
        animationScene1.AnimationTime = round(i,2)
        renderView1.ViewTime = round(i,2)
        print('Plotting frame number {}'.format(cnt))
        SaveScreenshot(fname, renderView1, ImageResolution=res,TransparentBackground=1)

# create a new 'OpenFOAMReader'
damBreakWithObstaclefoam = OpenFOAMReader(FileName=sys.argv[1]+'.foam')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on damBreakWithObstaclefoam
damBreakWithObstaclefoam.SkipZeroTime = 0

# Properties modified on damBreakWithObstaclefoam
damBreakWithObstaclefoam.CellArrays = ['alpha.water']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1526, 789]

# show data in view
damBreakWithObstaclefoamDisplay = Show(damBreakWithObstaclefoam, renderView1)
# trace defaults for the display properties.
damBreakWithObstaclefoamDisplay.Representation = 'Surface'
damBreakWithObstaclefoamDisplay.DiffuseColor = [0.996078431372549, 0.996078431372549, 0.996078431372549]
# Properties modified on damBreakWithObstaclefoamDisplay
damBreakWithObstaclefoamDisplay.BackfaceRepresentation = 'Cull Backface'

# Properties modified on damBreakWithObstaclefoamDisplay
damBreakWithObstaclefoamDisplay.BackfaceRepresentation = 'Cull Frontface'
# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
#Left view
#renderView1.CameraPosition = [-0.69, 0.83, 2.42]
#renderView1.CameraFocalPoint = [1.31, 0.037509902814557905, 0.27]
#renderView1.CameraViewUp = [0.11, 0.96, -0.25]
#renderView1.CameraParallelScale = 1.7

#Right view
#renderView1.CameraPosition = [4.13, 0.94, 2.43]
#renderView1.CameraFocalPoint = [1.88, 0.16, 0.34]
#renderView1.CameraViewUp = [-0.17, 0.97, -0.18]
#renderView1.CameraParallelScale = 1.76

#Center view
renderView1.CameraPosition = [1.65, 1.89, 3.88]
renderView1.CameraFocalPoint = [1.62, 0.43, 0.33]
renderView1.CameraViewUp = [0, 0.92, -0.38]
renderView1.CameraParallelScale = 1.76

# update the view to ensure updated data information
renderView1.Update()

renderView1.OrientationAxesVisibility = 0

# create a new 'Clip'
clip1 = Clip(Input=damBreakWithObstaclefoam)

# Properties modified on clip1
clip1.ClipType = 'Scalar'
clip1.Value = 0.5

# show data in view
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.Representation = 'Surface'

# hide data in view
#Hide(damBreakWithObstaclefoam, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# change solid color
clip1Display.DiffuseColor = [0.20392156862745098, 0.2823529411764706, 1.0]

clip1Display.Triangulate = 1

animationScene1.GoToFirst()

# save screenshot
res=[1530, 790]
frames = np.arange(float(sys.argv[2]), float(sys.argv[3]), 0.02)
# save screenshot
screenshot(frames, 'final', res)

