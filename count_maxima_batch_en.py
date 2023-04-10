# Counting the number of fluorescence maxima per cell
# 
# This FIJI macro allows counting of intracellular fluorescent condensates observed in time-lapse microscopy of yeast cells.
# As input for the timel-lapse image data, a hyperstack file that contains a pair of phase contrast (Ch #1) and fluorescence (Ch #2) z-stack images at each frame is required.
# This macro is optimized for images accuired by using a 60x objective lens.

# author: Masak Takaine

from ij import IJ, ImagePlus, Prefs
from ij import WindowManager as WM
from ij.process import ImageStatistics as IS
options = IS.MEAN | IS.AREA | IS.STD_DEV
from ij.gui import Roi
from ij.plugin.frame import RoiManager as RM
from ij.measure import ResultsTable
from ij.measure import Measurements as ms
#options_is = IS.MEAN | IS.MEDIAN | IS.MODE
#options_imp = ImagePlus.MEAN | ImagePlus.MEDIAN | ImagePlus.MODE
options_ms = ms.MEAN | ms.MEDIAN | ms.MODE

import os
from os import path
from ij.plugin import ChannelSplitter
from ij.plugin import Duplicator
from ij.plugin import ZProjector as zp
from ij.plugin.filter import BackgroundSubtracter
createBackground = False
lightBackground = False
useParaboloid = False
doPresmooth = False
correctCorners = False

from ij.plugin.filter import GaussianBlur
from ij.plugin.filter import Analyzer
from ij.plugin.filter import ParticleAnalyzer as PA
options_pa = PA.SHOW_NONE + PA.CLEAR_WORKSHEET+PA.ADD_TO_MANAGER # +PA.EXCLUDE_EDGE_PARTICLES + PA.SHOW_RESULTS+PA.SHOW_OUTLINES

from ij.plugin.filter import MaximumFinder
outputType = MaximumFinder.COUNT
excludeOnEdges = False

#@ String(label="Date of experiments, e.g., 2020-08-22") edate1
#@ File (label="Choose source Folder", style="directory") dirS0
#@ File (label="Choose destination Folder", style="directory") dirD0
#@ int(label="Prominence #1 for FindMaxima function") prom1
#@ int(label="Prominence #2 for FindMaxima function") prom2
#@ int(label="Prominence #3 for FindMaxima function") prom3
#@ int(label="Prominence #3 for FindMaxima function") prom4
#@ String(label="Hide/Show the active image? The Show slows the analysis.", choices={"hide","show"}, style="radioButtonHorizontal") sbm

proms = []
proms.append(prom1)
proms.append(prom2)
proms.append(prom3)
proms.append(prom4)

# Segmentation of phase contrast images.
def ch1fr1_segment(imp, radius, sigma, minsize_um): # BackgroundSubtraction radius, GaussianBlur sigma, minimum size (µm^2) to analyze
	filename = imp.getTitle().split(".")[0]  # file name without a extension
	cal = imp.getCalibration() # Returns this image's calibration
	ch1 = ChannelSplitter.split(imp)[0]
	ch1fr1 = Duplicator().run(ch1, 1, 1, 1, 11, 1, 1)
	ch1fr1_mip = zp.run(ch1fr1,"max all")
	
	bs = BackgroundSubtracter()
#	radius = 25
	bs.rollingBallBackground(ch1fr1_mip.getProcessor(), radius, createBackground, lightBackground, useParaboloid, doPresmooth, correctCorners)
	
	gb = GaussianBlur()
#	sigma = 1
	gb.blurGaussian(ch1fr1_mip.getProcessor(), sigma)
	
	IJ.setAutoThreshold(ch1fr1_mip, "Default dark")
	Prefs.blackBackground = True
	
	IJ.run(ch1fr1_mip, "Convert to Mask", "")
	IJ.run(ch1fr1_mip, "Open", "")
	IJ.run(ch1fr1_mip, "Fill Holes", "")
	IJ.run(ch1fr1_mip, "Adjustable Watershed", "tolerance=3")
	rt = ResultsTable()
	
	MAXSIZE = 1000000
	MINSIZE = minsize_um * (1/cal.pixelWidth) ** 2 # equals to 5 µm^2 in pixel unit
	
	p = PA(options_pa, PA.AREA + PA.MEAN + PA.PERIMETER + PA.FERET + PA.SHAPE_DESCRIPTORS, rt, MINSIZE, MAXSIZE)
	p.analyze(ch1fr1_mip)
	
	rm = RM.getRoiManager()
	rm.runCommand(ch1fr1_mip,"Show None")
	rm.runCommand(ch1fr1_mip,"Show All with labels")
	
	nResults = rt.size()

	for i in range(0, nResults): 
		rt.setValue("file", i, filename)
		rt.setValue("radius", i, radius)
		rt.setValue("sigma", i, sigma)
		rt.setValue("minsize_um", i, minsize_um)
		rt.setValue("MINSIZE", i, MINSIZE)

	ch1fr1_mip.setTitle("mask-"+filename)

	return ch1fr1_mip, rt, nResults, filename, radius, sigma, MINSIZE
	# These return values are gathered into a tapple. 

def preprocess(imp, last_frame, radius, sigma):
	ch1 = ChannelSplitter.split(imp)[0] 
#	last_frame = 15
	ch1sub = Duplicator().run(ch1, 1, 1, 1, 11, 1, last_frame)
#	radius = 25
	roll_condition = "rolling="+str(radius)+" stack" 
	IJ.run(ch1sub, "Subtract Background...", roll_condition)
		
	ch1sub_mip = zp.run(ch1sub, "max all")
	IJ.run(ch1sub_mip, "Subtract Background...", roll_condition)
		
	stats = ch1sub_mip.getStatistics(options_ms)
	bc_option = "correction=[Simple Ratio] background=" + str(stats.median)
	IJ.run(ch1sub_mip, "Bleach Correction", bc_option)
	
	winlist = WM.getImageTitles() # Acquired a namelist of windows open
	l_dup = [t for t in winlist if t.startswith("DUP")]
	IJ.selectWindow(l_dup[0])					
	ch1sub_mip_bc = IJ.getImage()					
	ch1sub_mip_bc.hide()					
	
#	sigma = 0.5
# "stack" allows application to all slices.
	IJ.run(ch1sub_mip_bc, "Gaussian Blur...", "sigma="+str(sigma)+" stack")  # A white space is required before "stack"
	return ch1sub_mip_bc, last_frame, radius, sigma
		
def findmaxima_stack(imp, prom):
	stack = imp.getStack()
	mf = MaximumFinder()
	# iterate over all frames
	for frame_no in range(1, imp.getNFrames()+ 1):
	    stack_index = imp.getStackIndex(1, 1, frame_no)
	    ip = stack.getProcessor(stack_index)
	    mf.findMaxima(ip, prom, outputType, excludeOnEdges)
	    
	rt = ResultsTable.getResultsTable()
	for k in range(0, rt.size()):
		rt.setValue("total_cell_num",k, total_cell_num)
		rt.setValue("count_per_cell",k, rt.getValue("Count", k)/total_cell_num)
		rt.setValue("prominence",k, prom)
		rt.setValue("date", k, edate1)
		rt.setValue("file", k, filename)
	return rt, prom

# Save the Result Table in csv format
def save_result_table(directory,result_table, prefix, filename):
	title = prefix+filename
	resultfile = os.path.join(directory, title + ".csv") 
	result_table.saveAs(resultfile)

# Save the mask image in tiff format
def save_mask_as_tif(directory, mask):
    title = mask.getTitle()
    outputfile = os.path.join(directory, title)
    IJ.saveAsTiff(mask, outputfile)

################### Main code
# Make directories
edate = " "+edate1 ## Insert a blank to prevent automatic modification on Excel.
dirSG = os.path.join(str(dirD0), edate1+"_segmentation")
if not os.path.exists(dirSG): 
	os.mkdir(dirSG)								
dirMX = os.path.join(str(dirD0), edate1+"_maxima_counts")
if not os.path.exists(dirMX):
	os.mkdir(dirMX)
dirDR = os.path.join(str(dirD0), edate1+"_drawings") # Create a folder for mask images
if not os.path.exists(dirDR):
	os.mkdir(dirDR)

# Acquire a list of files in the directory
filelist = os.listdir(str(dirS0))
# List comprehension, extract nd2 files.
nd2_files = [f for f in filelist if f.split(".")[-1] == "nd2"]

for nd2_file in nd2_files:
	current_file_path = os.path.join(str(dirS0), nd2_file)
	imp = IJ.openImage(current_file_path) 		# Create a ImagePlus object, assign it into imp
	r_seg = ch1fr1_segment(imp, 30, 1, 5.0)		
	
	mask = r_seg[0]			# the mask image
	rt_seg = r_seg[1]		# the Result Table
	total_cell_num = r_seg[2] # the number of cells detected
	filename = r_seg[3] 	# the file name without extension
	
	# Add a column of "date" to the table
	for m in range(0, rt_seg.size()): 
		rt_seg.setValue("date", m, edate)
		
	save_result_table(str(dirSG), rt_seg, "cell_segments_", filename)
	save_mask_as_tif(str(dirDR), mask)

	prep = preprocess(imp,12,25,0.5)
	imp2 = prep[0]
	last_frame = prep[1] 
	radius = prep[2] 
	sigma = prep[3]
	
	for prom in proms:
		out = findmaxima_stack(imp2, prom)
		rt = out[0]
		for l in range(0, rt.size()):
			rt.setValue("last_frame", l, last_frame)
			rt.setValue("radius", l, radius)
			rt.setValue("sigma", l, sigma)
		prefix = "maxima_counts_prom" + str(prom) + "_"
		save_result_table(str(dirMX), rt, prefix, filename)
		IJ.run("Clear Results")

print "Done. \n"
IJ.run("Clear Results")
rm = RM.getRoiManager()
rm.reset()
IJ.run("Close All")
