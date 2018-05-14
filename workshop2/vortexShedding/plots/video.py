#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
vortexfoam = OpenFOAMReader(FileName='foam.foam')

# Properties modified on vortexfoam
vortexfoam.CaseType = 'Decomposed Case'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1600, 802]

# show data in view
vortexfoamDisplay = Show(vortexfoam, renderView1)
# trace defaults for the display properties.
vortexfoamDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
vortexfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set scalar coloring
ColorBy(vortexfoamDisplay, ('POINTS', 'U', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
vortexfoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
vortexfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')

#change interaction mode for render view
renderView1.InteractionMode = '2D'

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()
uLUT.ApplyPreset('Cool to Warm (Extended)', True)
# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)

# change scalar bar placement
uLUTColorBar.Orientation = 'Horizontal'
uLUTColorBar.WindowLocation = 'LowerCenter'

# Properties modified on uLUTColorBar
uLUTColorBar.TitleFontFamily = 'Times'
uLUTColorBar.TitleItalic = 1
uLUTColorBar.TitleFontSize = 7
uLUTColorBar.LabelFontFamily = 'Times'
uLUTColorBar.LabelItalic = 1
uLUTColorBar.LabelFontSize = 7
uLUTColorBar.AutomaticLabelFormat = 0
uLUTColorBar.LabelFormat = '%-#1.1f'
uLUTColorBar.RangeLabelFormat = '%-#1.1f'
uLUTColorBar.ScalarBarThickness = 5
uLUTColorBar.ScalarBarLength = 0.3

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
nf= int(float(sys.argv[1]) / float(sys.argv[2]))
GetAnimationScene().NumberOfFrames = nf
GetAnimationScene().EndTime = float(sys.argv[1])

#horizontal pipe
GetAnimationScene().GoToFirst()
renderView1.CameraPosition = [8.08, 0, 59.21]
renderView1.CameraFocalPoint = [8.08, 0, 0.0]
renderView1.CameraParallelScale = 4.41

# save animation
SaveAnimation('plots/frames/vortex.png', renderView1, ImageResolution=[1600, 802], FrameWindow=[0, nf])
  
# save snapshot  
animationScene1.AnimationTime = 250
renderView1.ViewTime = 250
SaveScreenshot('plots/frames/250.png', renderView1, ImageResolution=[1600, 802], TransparentBackground=1)
