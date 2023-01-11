import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

#
# BoneModelExtractor
#

class BoneModelExtractor(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BoneModelExtractor"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["BoneModelExtractor"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Saima Safdar"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
This module takes image stacks and produces nrrd file, auto-segmentation file, models (ply files) and labels for all the bones based on the Computed tomography (CT).
See more information in <a href="https://github.com/organization/projectname#BoneModelExtractor">module documentation</a>.
"""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)

#
# Register sample data sets in Sample Data module
#

def registerSampleData():
  """
  Add data sets to Sample Data module.
  """
  # It is always recommended to provide sample data for users to make it easy to try the module,
  # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

  import SampleData
  iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

  # To ensure that the source code repository remains small (can be downloaded and installed quickly)
  # it is recommended to store data sets that are larger than a few MB in a Github release.

  # BoneModelExtractor1
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='BoneModelExtractor',
    sampleName='BoneModelExtractor1',
    # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
    # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
    thumbnailFileName=os.path.join(iconsPath, 'BoneModelExtractor1.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
    fileNames='BoneModelExtractor1.nrrd',
    # Checksum to ensure file integrity. Can be computed by this command:
    #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
    checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
    # This node name will be used when the data set is loaded
    nodeNames='BoneModelExtractor1'
  )

  # BoneModelExtractor2
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='BoneModelExtractor',
    sampleName='BoneModelExtractor2',
    thumbnailFileName=os.path.join(iconsPath, 'BoneModelExtractor2.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
    fileNames='BoneModelExtractor2.nrrd',
    checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
    # This node name will be used when the data set is loaded
    nodeNames='BoneModelExtractor2'
  )

#
# BoneModelExtractorWidget
#

class BoneModelExtractorWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/BoneModelExtractor.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = BoneModelExtractorLogic()

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    

    # Buttons
    self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    if not self._parameterNode.GetNodeReference("InputVolume"):
      firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
      if firstVolumeNode:
        self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Update node selectors and sliders
    self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputVolume"))
    

    # Update buttons states and tooltips
    if self._parameterNode.GetNodeReference("InputVolume"):
      self.ui.applyButton.toolTip = "Compute output volume"
      self.ui.applyButton.enabled = True
    else:
      self.ui.applyButton.toolTip = "Select input volume nodes"
      self.ui.applyButton.enabled = False

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("InputVolume", self.ui.inputSelector.currentNodeID)
    
    self._parameterNode.EndModify(wasModified)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    try:

      # Compute output
      
         
      self.logic.process(self.ui.inputSelector.currentNode(), self.ui.iterationsVal.value, self.ui.smoothingModel.value, self.ui.modelInitial.text, self.ui.directoryName.currentPath)

    except Exception as e:
      slicer.util.errorDisplay("Failed to compute results: "+str(e))
      import traceback
      traceback.print_exc()


#
# BoneModelExtractorLogic
#

class BoneModelExtractorLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("Threshold"):
      parameterNode.SetParameter("Threshold", "100.0")
    if not parameterNode.GetParameter("Invert"):
      parameterNode.SetParameter("Invert", "false")

  def process(self, inputVolume,iterations, modelSmoothing, modelInitial, directoryName):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param thresholding: auto-threshold method selection
    :param iterations: number of iterations for the shrink wrap
    :param model smoothing: number to smooth the models (double)
    :param Model initial: name for the model
    :param Save models: directory to save all model and data
    """

    if not inputVolume:
      raise ValueError("Input is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

  
    #slicer morph 
    #imagestacks
    print(iterations)
    print(type(iterations))
   
    
    #segment editor. 
    #threshold
    # Create segmentation
    segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
    segmentationNode.CreateDefaultDisplayNodes() # only needed for display
    segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(inputVolume)
    
    # Create temporary segment editor to get access to effects
    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
    segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setMasterVolumeNode(inputVolume)
    
    # Create a new segment in segemnt editor effects
    addedSegmentID = segmentationNode.GetSegmentation().AddEmptySegment(modelInitial)
    segmentEditorNode.SetSelectedSegmentID(addedSegmentID)
    
    # thresholding
    #getting bone threshold value min and maximum
    import vtkITK
    thresholdCalculator = vtkITK.vtkITKImageThresholdCalculator()
    thresholdCalculator.SetInputData(inputVolume.GetImageData())
    thresholdCalculator.SetMethodToOtsu()
    thresholdCalculator.Update()
    boneThresholdValue = thresholdCalculator.GetThreshold()
    volumeScalarRange = inputVolume.GetImageData().GetScalarRange()
    logging.info("Volume minimum = {0}, maximum = {1}, bone threshold = {2}".format(volumeScalarRange[0], volumeScalarRange[1], boneThresholdValue))
    
    #applying threshold effect
    segmentEditorWidget.setActiveEffectByName("Threshold")
    effect = segmentEditorWidget.activeEffect()
    effect.setParameter("MinimumThreshold",str(boneThresholdValue))
    effect.setParameter("MaximumThreshold",str(volumeScalarRange[1]))
    #effect.setParameter("MinimumThreshold",0.0)#11130.40)#9386.40)
    #effect.setParameter("MaximumThreshold",0)#31049.00)#21807.00)
    #effect.setParameter("AutoThresholdMethod", "OTSU")
    #effect.setParameter("AutoThresholdMode","SET_UPPER")
    effect.self().onApply()
    
   
    #segmentEditorNode.SetOverwriteMode(slicer.vtkMRMLSegmentEditorNode.OverwriteNone)
    #set smoothing factor for the 3d view
    # Create Closed Surface Representation and set the smoothing of the closed surface representation
    segmentationNode.GetSegmentation().SetConversionParameter("Oversampling factor", "1.5")
    segmentationNode.GetSegmentation().SetConversionParameter("Joint smoothing", "0.00")
    segmentationNode.GetSegmentation().SetConversionParameter("Smoothing factor",str(modelSmoothing))
    segmentationNode.GetSegmentation().SetConversionParameter("Decimation factor", "0.00")
    segmentationNode.CreateClosedSurfaceRepresentation()
    
    print(directoryName)
    print(str(modelInitial))
    
    #remove small islands
    segmentEditorWidget.setActiveEffectByName("Islands")
    effect = segmentEditorWidget.activeEffect()
    operationName = 'REMOVE_SMALL_ISLANDS'
    minsize = 1000
    effect.setParameter("Operation", operationName)
    effect.setParameter("MinimumSize",minsize)
    effect.self().onApply()

    #split island set island min > 1000 to select the bones  
    operationName = 'SPLIT_ISLANDS_TO_SEGMENTS'
    minsize = 1000
    effect.setParameter("Operation", operationName)
    effect.setParameter("MinimumSize",minsize)
    effect.self().onApply()
    
    # Compute centroids of each bone and place the fiducial for labeling in 3d and image view
    import SegmentStatistics
    segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
    segStatLogic.getParameterNode().SetParameter("Segmentation", segmentationNode.GetID())
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.centroid_ras.enabled", str(True))
    segStatLogic.computeStatistics()
    stats = segStatLogic.getStatistics()
    
    #use shrink wrapping to wrap the bone and get the full surface 
    segmentEditorWidget.setActiveEffectByName("Wrap Solidify")
    effect = segmentEditorWidget.activeEffect()
    
    # Place a markup point in each centroid
    markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    #locking the fiducial as these are only used for labeling
    markupsNode.SetLocked(1)
    markupsNode.CreateDefaultDisplayNodes()
    for segmentId in stats['SegmentIDs']:
        centroid_ras = stats[segmentId,"LabelmapSegmentStatisticsPlugin.centroid_ras"]
        segmentName = segmentationNode.GetSegmentation().GetSegment(segmentId).GetName()
        markupsNode.AddFiducialFromArray(centroid_ras, segmentName)
        segmentEditorWidget.setCurrentSegmentID(segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(segmentName))
        effect.setParameter("region", "outerSurface")
        effect.setParameter("regionSegmentID", "segment")
        #effect.setParameter("REGION_OUTER_SURFACE", "outerSurface")
        #effect.setParameter("REGION_LARGEST_CAVITY", "largestCavity")
        effect.setParameter("carveHolesInOuterSurface", True)
        effect.setParameter("carveHolesInOuterSurfaceDiameter", 10)
        effect.setParameter("smoothingFactor", modelSmoothing)
        effect.setParameter("shrinkwrapIterations", iterations)
        effect.setParameter("outputType", "segment")
        #effect.setParameter("outputModelNode", segmentName)
        effect.self().onApply()
        #f.write(str(centroid_ras))
        #f.write("\n")
    #create sements 
    #extract each label and assign it to segment make visible the label on segment
    #run wrap solidify
    
    
    #export it to model ply file
    #segmentationNode.SetName("S")
    #slicer.modules.segmentations.logic().ExportSegmentsClosedSurfaceRepresentationToFiles("/home/saima/2021_5_28/slicer work/modelFiles/", segmentationNode, None,"PLY", True, 1.0, False )
    # Export segments to models
    #shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    #import SurfaceToolbox
    #logic = SurfaceToolbox.SurfaceToolboxLogic()
    #parameterNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode")
    #logic.setDefaultParameters(parameterNode)
    
    shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    exportFolderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), segmentationNode.GetName())
    slicer.modules.segmentations.logic().ExportAllSegmentsToModels(segmentationNode, exportFolderItemId)
    segmentModels = vtk.vtkCollection()
    shNode.GetDataNodesInBranch(exportFolderItemId, segmentModels)
    
    modelTreeID = shNode.GetItemChildWithName(shNode.GetSceneItemID(), segmentationNode.GetName() ) 
    modelTreeChildrenIDs = vtk.vtkIdList()
    shNode.GetItemChildren(modelTreeID, modelTreeChildrenIDs)
    for i in range(modelTreeChildrenIDs.GetNumberOfIds()):
      print("modelChildID   : ", modelTreeChildrenIDs.GetId(i)) 
      print("modelChildName : ", shNode.GetItemName(modelTreeChildrenIDs.GetId(i)))
      #print("modelChildNode  : ", shNode.GetItemDataNode(modelTreeChildrenIDs.GetId(i)))
      #modelName = shNode.GetItemName(modelTreeChildrenIDs.GetId(i))
      #model = shNode.GetItemDataNode(modelTreeChildrenIDs.GetId(i))
      #print(modelName)
      #print(model)
      modelNode = segmentModels.GetItemAsObject(i)
      #parameterNode.SetNodeReferenceID("inputModel", modelNode.GetID())
      #parameterNode.SetNodeReferenceID("outputModel", modelNode.GetID())
      #parameterNode.SetParameter("normals", "true")
      #logic.applyFilters(parameterNode)
      print(modelNode)
      slicer.modules.models.logic().SaveModel(directoryName+"/"+shNode.GetItemName(modelTreeChildrenIDs.GetId(i))+".PLY", modelNode)
      

    """models=slicer.util.getNodesByClass('vtkMRMLModelNode')
    i=1
    for model in models:   
      if not model.GetHideFromEditors():
        n = model.GetName() #n get the name of the model fromt he user
        n = str(modelInitial)
        slicer.modules.models.logic().SaveModel(directoryName+"/"+n+str(i)+".PLY", model)
        i=i+1 """
    slicer.util.saveNode(slicer.mrmlScene.GetNodeByID(segmentationNode.GetID()), directoryName+"/"+str(segmentationNode.GetName())+".seg.nrrd")   
    slicer.util.saveNode(inputVolume, directoryName+"/Data.nrrd")
    slicer.util.saveNode(markupsNode, directoryName+"/Labels.fcsv")    
   
    #slicer.mrmlScene.RemoveNode(shNode)
    segmentEditorWidget = None
    slicer.mrmlScene.RemoveNode(segmentEditorNode)
   
    
    

    

    stopTime = time.time()
    logging.info('Processing completed in {0:.2f} seconds'.format(stopTime-startTime))

#
# BoneModelExtractorTest
#

class BoneModelExtractorTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_BoneModelExtractor1()

  def test_BoneModelExtractor1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    # Get/create input data

    import SampleData
    registerSampleData()
    inputVolume = SampleData.downloadSample('BoneModelExtractor1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 100

    # Test the module logic

    logic = BoneModelExtractorLogic()

    # Test algorithm with non-inverted threshold
    logic.process(inputVolume, outputVolume, threshold, True)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputVolume, outputVolume, threshold, False)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
