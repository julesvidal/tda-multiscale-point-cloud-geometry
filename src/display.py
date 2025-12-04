# state file generated using paraview version 5.11.1-1-g70075740e
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [728, 817]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [4.016390623170057, 2.8706062811019133, 5.590519113423502]
renderView1.CameraViewUp = [-0.19467038376336326, 0.9223433148376388, -0.33374578837753355]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.9437987771755119
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'Grid Axes 3D Actor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitle = 'Birth'
renderView1.AxesGrid.YTitle = 'Death'
renderView1.AxesGrid.ZTitle = '$-f$'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontSize = 24
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontSize = 24
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontSize = 24
renderView1.AxesGrid.FacesToRender = 21
renderView1.AxesGrid.GridColor = [0.282352941176471, 0.27843137254902, 0.27843137254902]
renderView1.AxesGrid.ShowGrid = 1
renderView1.AxesGrid.AxesToLabel = 14
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontSize = 24
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontSize = 24
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontSize = 24
renderView1.AxesGrid.XAxisNotation = 'Fixed'
renderView1.AxesGrid.XAxisPrecision = 0
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0]
renderView1.AxesGrid.YAxisNotation = 'Fixed'
renderView1.AxesGrid.YAxisPrecision = 0
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [0.0]
renderView1.AxesGrid.ZAxisPrecision = 0
renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.UseCustomBounds = 1
renderView1.AxesGrid.CustomBounds = [-7.0, 7.0, -7.0, 7.0, 1.0, 10.0]
renderView1.AxesGrid.DataBoundsScaleFactor = 1.0

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [728, 825]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'Grid Axes 3D Actor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [-7.09992253780365, -9.332368202507496, -5.6673361621797085]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-7.09992253780365, -9.332368202507496, 46.119549029881426]
renderView2.CameraFocalPoint = [-7.09992253780365, -9.332368202507496, -5.6673361621797085]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 11.160257318558584
renderView2.LegendGrid = 'Legend Grid Actor'
renderView2.PolarGrid = 'Polar Grid Actor'
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# init the 'Grid Axes 3D Actor' selected for 'AxesGrid'
renderView2.AxesGrid.XTitle = 'Birth'
renderView2.AxesGrid.YTitle = 'Death'
renderView2.AxesGrid.ZTitle = '$-f$'
renderView2.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.XTitleFontSize = 24
renderView2.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.YTitleFontSize = 24
renderView2.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.ZTitleFontSize = 24
renderView2.AxesGrid.FacesToRender = 21
renderView2.AxesGrid.GridColor = [0.282352941176471, 0.27843137254902, 0.27843137254902]
renderView2.AxesGrid.ShowGrid = 1
renderView2.AxesGrid.AxesToLabel = 14
renderView2.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.XLabelFontSize = 24
renderView2.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.YLabelFontSize = 24
renderView2.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.ZLabelFontSize = 24
renderView2.AxesGrid.XAxisNotation = 'Fixed'
renderView2.AxesGrid.XAxisPrecision = 0
renderView2.AxesGrid.XAxisUseCustomLabels = 1
renderView2.AxesGrid.XAxisLabels = [0.0]
renderView2.AxesGrid.YAxisNotation = 'Fixed'
renderView2.AxesGrid.YAxisPrecision = 0
renderView2.AxesGrid.YAxisUseCustomLabels = 1
renderView2.AxesGrid.YAxisLabels = [0.0]
renderView2.AxesGrid.ZAxisPrecision = 0
renderView2.AxesGrid.ZAxisUseCustomLabels = 1
renderView2.AxesGrid.UseCustomBounds = 1
renderView2.AxesGrid.CustomBounds = [-7.0, 7.0, -7.0, 7.0, 1.0, 10.0]
renderView2.AxesGrid.DataBoundsScaleFactor = 1.0

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)
layout1.SetSize(1457, 825)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML MultiBlock Data Reader'
diagramvtm = XMLMultiBlockDataReader(registrationName='diagram.vtm', FileName=['NAME_diagram.vtm'])
diagramvtm.CellArrayStatus = ['PairIdentifier', 'PairType', 'Persistence', 'Birth', 'IsFinite', 'score', 'chosenId', 'scaleRange']
diagramvtm.PointArrayStatus = ['ttkVertexScalarField', 'CriticalType', 'Coordinates']
diagramvtm.TimeArray = 'None'

# create a new 'Extract Block'
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=diagramvtm)
extractBlock1.Assembly = 'Hierarchy'
extractBlock1.Selectors = ['/Root/Block1']

# create a new 'Merge Blocks'
mergeBlocks3 = MergeBlocks(registrationName='MergeBlocks3', Input=extractBlock1)
mergeBlocks3.MergePoints = 0

# create a new 'XML MultiBlock Data Reader'
primitivesvtm = XMLMultiBlockDataReader(registrationName='primitives.vtm', FileName=['NAME_primitives.vtm'])
primitivesvtm.PointArrayStatus = ['Normals', 'size', 'isBaryals', 'meanNormals', 'baryOffset', 'primitiveId', 'chosenId', 'projCoords', 'projNormals']
primitivesvtm.TimeArray = 'None'

# create a new 'Calculator'
flag = Calculator(registrationName='flag', Input=mergeBlocks3)
flag.AttributeType = 'Cell Data'
flag.ResultArrayName = 'flag'
flag.Function = 'chosenId'

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=flag)
calculator2.ResultArrayName = 'z'
calculator2.Function = 'coordsZ'

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(registrationName='PythonCalculator1', Input=calculator2)
pythonCalculator1.Expression = 'z.GetRange()'

# create a new 'Calculator'
calculator4 = Calculator(registrationName='Calculator4', Input=pythonCalculator1)
calculator4.ResultArrayName = 'y'
calculator4.Function = 'coordsY'

# create a new 'Python Calculator'
pythonCalculator2 = PythonCalculator(registrationName='PythonCalculator2', Input=calculator4)
pythonCalculator2.Expression = 'y.GetRange()'
pythonCalculator2.ArrayName = 'yRange'

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=pythonCalculator2)
calculator3.CoordinateResults = 1
calculator3.Function = 'coordsX * iHat + coordsY*kHat + (yRange_X+coordsZ - (yRange_Y-yRange_X)/3)*jHat'

# create a new 'Append Datasets'
appendDatasets1 = AppendDatasets(registrationName='AppendDatasets1', Input=[calculator3, pythonCalculator2])

# create a new 'Threshold'
dISCARDED = Threshold(registrationName='DISCARDED', Input=appendDatasets1)
dISCARDED.Scalars = ['CELLS', 'flag']
dISCARDED.LowerThreshold = -1.0
dISCARDED.UpperThreshold = 4.0
dISCARDED.ThresholdMethod = 'Below Lower Threshold'

# create a new 'Extract Edges'
extractEdges3 = ExtractEdges(registrationName='ExtractEdges3', Input=dISCARDED)

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(registrationName='CellDatatoPointData2', Input=extractEdges3)

# create a new 'Calculator'
calculator6 = Calculator(registrationName='Calculator6', Input=cellDatatoPointData2)
calculator6.CoordinateResults = 1
calculator6.Function = 'coords + kHat*(-10-flag)'

# create a new 'Threshold'
cHOSEN = Threshold(registrationName='CHOSEN', Input=appendDatasets1)
cHOSEN.Scalars = ['CELLS', 'flag']
cHOSEN.LowerThreshold = -2.0
cHOSEN.ThresholdMethod = 'Above Upper Threshold'

# create a new 'Extract Edges'
extractEdges2 = ExtractEdges(registrationName='ExtractEdges2', Input=cHOSEN)

# create a new 'TTK IdentifierRandomizer'
tTKIdentifierRandomizer1 = TTKIdentifierRandomizer(registrationName='TTKIdentifierRandomizer1', Input=cHOSEN)
tTKIdentifierRandomizer1.ScalarField = ['CELLS', 'chosenId']
tTKIdentifierRandomizer1.RandomSeed = 1
tTKIdentifierRandomizer1.CompactRange = 1

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=extractEdges2)

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=extractSurface1)
cellDatatoPointData1.CellDataArraytoprocess = ['Persistence', 'birth', 'chosenId', 'death', 'flag', 'maxScale', 'minScale', 'productPersistence', 'scalePersistence', 'score']

# create a new 'Mask Points'
maskPoints1 = MaskPoints(registrationName='MaskPoints1', Input=cellDatatoPointData1)
maskPoints1.OnRatio = 0
maskPoints1.MaximumNumberofPoints = 50000000
maskPoints1.GenerateVertices = 1
maskPoints1.SingleVertexPerCell = 1

# create a new 'Tube'
tube1 = Tube(registrationName='Tube1', Input=extractSurface1)
tube1.Scalars = ['POINTS', '']
tube1.Vectors = ['POINTS', '1']
tube1.Radius = 0.13

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from diagramvtm
diagramvtmDisplay = Show(diagramvtm, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'vtkBlockColors'
vtkBlockColorsTF2D = GetTransferFunction2D('vtkBlockColors')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.AutomaticRescaleRangeMode = "Update on 'Apply'"
vtkBlockColorsLUT.TransferFunction2D = vtkBlockColorsTF2D
vtkBlockColorsLUT.RGBPoints = [255.0, 0.23137254902, 0.298039215686, 0.752941176471, 306.05321335669095, 0.865, 0.865, 0.865, 357.1064267133819, 0.705882352941, 0.0156862745098, 0.149019607843]
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')
vtkBlockColorsPWF.Points = [0.0, 1.0, 0.5, 0.0, 33.28009843826294, 1.0, 0.5, 0.0, 45.0, 1.0, 0.5, 0.0]

# trace defaults for the display properties.
diagramvtmDisplay.Representation = 'Surface'
diagramvtmDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
diagramvtmDisplay.LookupTable = vtkBlockColorsLUT
diagramvtmDisplay.Specular = 1.0
diagramvtmDisplay.SelectNormalArray = 'None'
diagramvtmDisplay.SelectTangentArray = 'None'
diagramvtmDisplay.SelectTCoordArray = 'None'
diagramvtmDisplay.TextureTransform = 'Transform2'
diagramvtmDisplay.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
diagramvtmDisplay.OSPRayScaleArray = 'Coordinates'
diagramvtmDisplay.OSPRayScaleFunction = 'Piecewise Function'
diagramvtmDisplay.Assembly = 'Hierarchy'
diagramvtmDisplay.SelectedBlockSelectors = ['']
diagramvtmDisplay.SelectOrientationVectors = 'Coordinates'
diagramvtmDisplay.ScaleFactor = 1.1334692459730833
diagramvtmDisplay.SelectScaleArray = 'Coordinates'
diagramvtmDisplay.GlyphType = 'Arrow'
diagramvtmDisplay.GlyphTableIndexArray = 'Coordinates'
diagramvtmDisplay.GaussianRadius = 0.05667346229865416
diagramvtmDisplay.SetScaleArray = ['POINTS', 'Coordinates']
diagramvtmDisplay.ScaleTransferFunction = 'Piecewise Function'
diagramvtmDisplay.OpacityArray = ['POINTS', 'Coordinates']
diagramvtmDisplay.OpacityTransferFunction = 'Piecewise Function'
diagramvtmDisplay.DataAxesGrid = 'Grid Axes Representation'
diagramvtmDisplay.PolarAxes = 'Polar Axes Representation'
diagramvtmDisplay.ScalarOpacityFunction = vtkBlockColorsPWF
diagramvtmDisplay.ScalarOpacityUnitDistance = 0.2076315595217025
diagramvtmDisplay.OpacityArrayName = ['POINTS', 'Coordinates']
diagramvtmDisplay.SelectInputVectors = ['POINTS', 'Coordinates']
diagramvtmDisplay.WriteLog = ''

# show data from primitivesvtm
primitivesvtmDisplay = Show(primitivesvtm, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'chosenId'
chosenIdTF2D = GetTransferFunction2D('chosenId')

# get color transfer function/color map for 'chosenId'
chosenIdLUT = GetColorTransferFunction('chosenId')
chosenIdLUT.AutomaticRescaleRangeMode = "Update on 'Apply'"
chosenIdLUT.TransferFunction2D = chosenIdTF2D
chosenIdLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 11.44, 0.0, 0.0, 0.360784313725, 22.799999999999997, 0.0, 1.0, 1.0, 34.32, 0.0, 0.501960784314, 0.0, 45.67999999999999, 1.0, 1.0, 0.0, 57.11999999999999, 1.0, 0.380392156863, 0.0, 68.56, 0.419607843137, 0.0, 0.0, 80.0, 0.878431372549, 0.301960784314, 0.301960784314]
chosenIdLUT.ColorSpace = 'RGB'
chosenIdLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chosenId'
chosenIdPWF = GetOpacityTransferFunction('chosenId')
chosenIdPWF.Points = [0.0, 1.0, 0.5, 0.0, 59.16461944580078, 1.0, 0.5, 0.0, 80.0, 1.0, 0.5, 0.0]
chosenIdPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
primitivesvtmDisplay.Representation = 'Surface'
primitivesvtmDisplay.ColorArrayName = ['POINTS', 'chosenId']
primitivesvtmDisplay.LookupTable = chosenIdLUT
primitivesvtmDisplay.PointSize = 5.0
primitivesvtmDisplay.Specular = 1.0
primitivesvtmDisplay.SelectNormalArray = 'Normals'
primitivesvtmDisplay.SelectTangentArray = 'None'
primitivesvtmDisplay.SelectTCoordArray = 'None'
primitivesvtmDisplay.TextureTransform = 'Transform2'
primitivesvtmDisplay.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
primitivesvtmDisplay.OSPRayScaleArray = 'Normals'
primitivesvtmDisplay.OSPRayScaleFunction = 'Piecewise Function'
primitivesvtmDisplay.Assembly = 'Hierarchy'
primitivesvtmDisplay.SelectedBlockSelectors = ['']
primitivesvtmDisplay.SelectOrientationVectors = 'None'
primitivesvtmDisplay.ScaleFactor = 0.2
primitivesvtmDisplay.SelectScaleArray = 'None'
primitivesvtmDisplay.GlyphType = 'Arrow'
primitivesvtmDisplay.GlyphTableIndexArray = 'None'
primitivesvtmDisplay.GaussianRadius = 0.01
primitivesvtmDisplay.SetScaleArray = ['POINTS', 'Normals']
primitivesvtmDisplay.ScaleTransferFunction = 'Piecewise Function'
primitivesvtmDisplay.OpacityArray = ['POINTS', 'Normals']
primitivesvtmDisplay.OpacityTransferFunction = 'Piecewise Function'
primitivesvtmDisplay.DataAxesGrid = 'Grid Axes Representation'
primitivesvtmDisplay.PolarAxes = 'Polar Axes Representation'
primitivesvtmDisplay.ScalarOpacityFunction = chosenIdPWF
primitivesvtmDisplay.ScalarOpacityUnitDistance = 0.12600950653072795
primitivesvtmDisplay.OpacityArrayName = ['POINTS', 'Normals']
primitivesvtmDisplay.SelectInputVectors = ['POINTS', 'Normals']
primitivesvtmDisplay.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
primitivesvtmDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
primitivesvtmDisplay.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
primitivesvtmDisplay.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(diagramvtm, renderView1)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', 'chosenId']
cellDatatoPointData1Display.LookupTable = chosenIdLUT
cellDatatoPointData1Display.InterpolateScalarsBeforeMapping = 0
cellDatatoPointData1Display.LineWidth = 2.0
cellDatatoPointData1Display.Specular = 1.0
cellDatatoPointData1Display.SelectNormalArray = 'None'
cellDatatoPointData1Display.SelectTangentArray = 'None'
cellDatatoPointData1Display.SelectTCoordArray = 'None'
cellDatatoPointData1Display.TextureTransform = 'Transform2'
cellDatatoPointData1Display.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
cellDatatoPointData1Display.OSPRayScaleArray = 'flag'
cellDatatoPointData1Display.OSPRayScaleFunction = 'Piecewise Function'
cellDatatoPointData1Display.Assembly = ''
cellDatatoPointData1Display.SelectedBlockSelectors = ['']
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 1.8664757233290403
cellDatatoPointData1Display.SelectScaleArray = 'flag'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'flag'
cellDatatoPointData1Display.GaussianRadius = 0.09332378616645201
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'flag']
cellDatatoPointData1Display.ScaleTransferFunction = 'Piecewise Function'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'flag']
cellDatatoPointData1Display.OpacityTransferFunction = 'Piecewise Function'
cellDatatoPointData1Display.DataAxesGrid = 'Grid Axes Representation'
cellDatatoPointData1Display.PolarAxes = 'Polar Axes Representation'
cellDatatoPointData1Display.SelectInputVectors = [None, '']
cellDatatoPointData1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
cellDatatoPointData1Display.OSPRayScaleFunction.Points = [255.0, 0.0, 0.5, 0.0, 357.1064267133819, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 80.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 80.0, 1.0, 0.5, 0.0]

# init the 'Polar Axes Representation' selected for 'PolarAxes'
cellDatatoPointData1Display.PolarAxes.EnableOverallColor = 0
cellDatatoPointData1Display.PolarAxes.ArcTickMatchesRadialAxes = 0

# show data from tube1
tube1Display = Show(tube1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = ['CELLS', 'chosenId']
tube1Display.LookupTable = chosenIdLUT
tube1Display.InterpolateScalarsBeforeMapping = 0
tube1Display.LineWidth = 2.0
tube1Display.Specular = 1.0
tube1Display.SelectNormalArray = 'TubeNormals'
tube1Display.SelectTangentArray = 'None'
tube1Display.SelectTCoordArray = 'None'
tube1Display.TextureTransform = 'Transform2'
tube1Display.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
tube1Display.OSPRayScaleArray = 'y'
tube1Display.OSPRayScaleFunction = 'Piecewise Function'
tube1Display.Assembly = ''
tube1Display.SelectedBlockSelectors = ['']
tube1Display.SelectOrientationVectors = 'None'
tube1Display.ScaleFactor = 1.8924756392836573
tube1Display.SelectScaleArray = 'y'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'y'
tube1Display.GaussianRadius = 0.09462378196418286
tube1Display.SetScaleArray = ['POINTS', 'y']
tube1Display.ScaleTransferFunction = 'Piecewise Function'
tube1Display.OpacityArray = ['POINTS', 'y']
tube1Display.OpacityTransferFunction = 'Piecewise Function'
tube1Display.DataAxesGrid = 'Grid Axes Representation'
tube1Display.PolarAxes = 'Polar Axes Representation'
tube1Display.SelectInputVectors = ['POINTS', 'TubeNormals']
tube1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
tube1Display.OSPRayScaleFunction.Points = [255.0, 0.0, 0.5, 0.0, 357.1064267133819, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
tube1Display.ScaleTransferFunction.Points = [-11.33468246459961, 0.0, 0.5, 0.0, 9.995131222240161e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
tube1Display.OpacityTransferFunction.Points = [-11.33468246459961, 0.0, 0.5, 0.0, 9.995131222240161e-06, 1.0, 0.5, 0.0]

# init the 'Polar Axes Representation' selected for 'PolarAxes'
tube1Display.PolarAxes.EnableOverallColor = 0
tube1Display.PolarAxes.ArcTickMatchesRadialAxes = 0

# show data from maskPoints1
maskPoints1Display = Show(maskPoints1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
maskPoints1Display.Representation = 'Surface'
maskPoints1Display.ColorArrayName = ['POINTS', 'chosenId']
maskPoints1Display.LookupTable = chosenIdLUT
maskPoints1Display.InterpolateScalarsBeforeMapping = 0
maskPoints1Display.PointSize = 20.0
maskPoints1Display.LineWidth = 2.0
maskPoints1Display.RenderPointsAsSpheres = 1
maskPoints1Display.Specular = 1.0
maskPoints1Display.SelectNormalArray = 'None'
maskPoints1Display.SelectTangentArray = 'None'
maskPoints1Display.SelectTCoordArray = 'None'
maskPoints1Display.TextureTransform = 'Transform2'
maskPoints1Display.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
maskPoints1Display.OSPRayScaleArray = 'flag'
maskPoints1Display.OSPRayScaleFunction = 'Piecewise Function'
maskPoints1Display.Assembly = ''
maskPoints1Display.SelectedBlockSelectors = ['']
maskPoints1Display.SelectOrientationVectors = 'None'
maskPoints1Display.ScaleFactor = 1.8664757233290403
maskPoints1Display.SelectScaleArray = 'flag'
maskPoints1Display.GlyphType = 'Arrow'
maskPoints1Display.GlyphTableIndexArray = 'flag'
maskPoints1Display.GaussianRadius = 0.09332378616645201
maskPoints1Display.SetScaleArray = ['POINTS', 'flag']
maskPoints1Display.ScaleTransferFunction = 'Piecewise Function'
maskPoints1Display.OpacityArray = ['POINTS', 'flag']
maskPoints1Display.OpacityTransferFunction = 'Piecewise Function'
maskPoints1Display.DataAxesGrid = 'Grid Axes Representation'
maskPoints1Display.PolarAxes = 'Polar Axes Representation'
maskPoints1Display.SelectInputVectors = [None, '']
maskPoints1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
maskPoints1Display.OSPRayScaleFunction.Points = [255.0, 0.0, 0.5, 0.0, 357.1064267133819, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
maskPoints1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 80.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
maskPoints1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 80.0, 1.0, 0.5, 0.0]

# init the 'Polar Axes Representation' selected for 'PolarAxes'
maskPoints1Display.PolarAxes.EnableOverallColor = 0
maskPoints1Display.PolarAxes.ArcTickMatchesRadialAxes = 0

# show data from calculator6
calculator6Display = Show(calculator6, renderView2, 'GeometryRepresentation')

# get 2D transfer function for 'flag'
flagTF2D = GetTransferFunction2D('flag')

# get color transfer function/color map for 'flag'
flagLUT = GetColorTransferFunction('flag')
flagLUT.InterpretValuesAsCategories = 1
flagLUT.AnnotationsInitialized = 1
flagLUT.AutomaticRescaleRangeMode = "Update on 'Apply'"
flagLUT.TransferFunction2D = flagTF2D
flagLUT.RGBPoints = [-4.0, 0.3, 0.3, 0.3, -2.958230972290039, 0.5254901960784314, 0.5254901960784314, 0.5254901960784314, -2.0, 0.8, 0.8, 0.8]
flagLUT.UseOpacityControlPointsFreehandDrawing = 1
flagLUT.ColorSpace = 'RGB'
flagLUT.NanColor = [1.0, 0.0, 0.0]
flagLUT.Annotations = ['-4', '-4', '-3', '-3', '-2', '-1', '-1', '-1']
flagLUT.ActiveAnnotatedValues = ['-4', '-2']
flagLUT.IndexedColors = [0.5490196078431373, 0.5490196078431373, 0.5490196078431373, 0.0, 0.09411764705882353, 0.07058823529411765, 0.8823529411764706, 0.8823529411764706, 0.8823529411764706, 0.8823529411764706, 0.8823529411764706, 0.8823529411764706]

# trace defaults for the display properties.
calculator6Display.Representation = 'Surface'
calculator6Display.ColorArrayName = ['POINTS', 'flag']
calculator6Display.LookupTable = flagLUT
calculator6Display.InterpolateScalarsBeforeMapping = 0
calculator6Display.LineWidth = 2.0
calculator6Display.Specular = 1.0
calculator6Display.SelectNormalArray = 'None'
calculator6Display.SelectTangentArray = 'None'
calculator6Display.SelectTCoordArray = 'None'
calculator6Display.TextureTransform = 'Transform2'
calculator6Display.EdgeColor = [0.23921568627450981, 0.23921568627450981, 0.23921568627450981]
calculator6Display.OSPRayScaleArray = 'flag'
calculator6Display.OSPRayScaleFunction = 'Piecewise Function'
calculator6Display.Assembly = ''
calculator6Display.SelectedBlockSelectors = ['']
calculator6Display.SelectOrientationVectors = 'None'
calculator6Display.ScaleFactor = 1.773208659887314
calculator6Display.SelectScaleArray = 'flag'
calculator6Display.GlyphType = 'Arrow'
calculator6Display.GlyphTableIndexArray = 'flag'
calculator6Display.GaussianRadius = 0.0886604329943657
calculator6Display.SetScaleArray = ['POINTS', 'flag']
calculator6Display.ScaleTransferFunction = 'Piecewise Function'
calculator6Display.OpacityArray = ['POINTS', 'flag']
calculator6Display.OpacityTransferFunction = 'Piecewise Function'
calculator6Display.DataAxesGrid = 'Grid Axes Representation'
calculator6Display.PolarAxes = 'Polar Axes Representation'
calculator6Display.SelectInputVectors = [None, '']
calculator6Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
calculator6Display.OSPRayScaleFunction.Points = [255.0, 0.0, 0.5, 0.0, 357.1064267133819, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
calculator6Display.ScaleTransferFunction.Points = [-4.0, 0.0, 0.5, 0.0, -2.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
calculator6Display.OpacityTransferFunction.Points = [-4.0, 0.0, 0.5, 0.0, -2.0, 1.0, 0.5, 0.0]

# init the 'Polar Axes Representation' selected for 'PolarAxes'
calculator6Display.PolarAxes.EnableOverallColor = 0
calculator6Display.PolarAxes.ArcTickMatchesRadialAxes = 0

# hide data in view
Hide(cellDatatoPointData1, renderView2)

# hide data in view
Hide(maskPoints1, renderView2)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'flag'
flagPWF = GetOpacityTransferFunction('flag')
flagPWF.Points = [0.0, 1.0, 0.5, 0.0, 33.28009843826294, 1.0, 0.5, 0.0, 45.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = [renderView1, renderView2]
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0

# ----------------------------------------------------------------
# restore active source
SetActiveSource(calculator6)

renderView1.ResetCamera(False, 0.9)
renderView1.AxesGrid.Visibility = 0
renderView2.ResetCamera(False, 0.9)
renderView2.AxesGrid.Visibility = 0

SetActiveSource(tube1)
SetActiveSource(primitivesvtm)
