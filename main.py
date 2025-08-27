import pandas as pd
from dataclasses import dataclass
import string
import numpy as np
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import pickle
import re
from alive_progress import alive_bar
import multiprocessing as mtp
from multiprocessing import Queue
from multiprocessing import Pool
from multiprocessing import Process
import time
import traceback
import tkinter as tk
from tkinter import filedialog
from tkinter import *
from tkinter.ttk import *
from bs4 import BeautifulSoup
import requests
import urllib.parse
import warnings
import time
import io
from io import StringIO
from PIL import ImageGrab, ImageFont, Image, ImageDraw
import mygene
from copy import deepcopy
import openpyxl
import umap
import math
from sklearn.preprocessing import StandardScaler
import seaborn
#from sksurv.nonparametric import kaplan_meier_estimator
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.datasets import load_waltons

# use a dictionary and have a text file with defined values
# TODO: This code is almost a thousand lines long, try to split it up into multiple files... <-- this first I guess

# TODO: heatmap of fischer transform for ulcerated vs non ulcerated
# TODO: More functionality! Class to interpret TCGA barcode, include sequencing and methylated regions, etc.

# TODO: train a classifier to segregate the N/A ulceration data
TCGA_PATIENT_ID_START = 0
TCGA_PATIENT_ID_END = 12

MIR_PREDICTION_ALGORITHM = 'MIR_PREDICTION_ALGORITHM'

class TCGA_ExpressionData:
    def __init__(self):
        self.gene_ID = ''
        #self.TCGA_Barcode: string = ''
        self.ExpressionValues = {}

    def CheckExpressionType(self, exp):
        if exp in self.ExpressionValues.keys():
            return True
        return False
'''
how do we classify data for each patient?
   * overview
       ** mortality
       ** dx's
       ** melanoma type
       ** melanoma site(s)
    * Clinical follow ups
'''

class ClinicalData:
    def __init__(self, alivestatus_a, dx_a, pfi_a, survival_a):
        self.alivestatus = alivestatus_a
        self.dx = dx_a # how to normalize so its consistent accross datasets
        self.dx_site = dx_a
        self.pfi = pfi_a
        self.survival = survival_a

        #self.follow_ups <-- deal with this later


class TCGA_Patient:
    def __init__(self, ID_a):
        self.ID = ID_a
        self.ExpressionData = {}
        self.ClinicalFeatures = {}

        #sequencing data too

    def AppendClinicalFeature(self, feature, data):
        if feature not in self.ClinicalFeatures.keys():
            self.ClinicalFeatures[feature] = data
        
    # stores the expression of a certain gene for this patient
    def StoreExpression(self, df, gene_ID):
        #get all columns with name, and then search for our gene of interest, if its not there, error and move on
        ''' Its a bit ugly, but it works, the  boolean array selects the first row so that we can have
            the type of expression data for the dictionary'''
        genedata = df.filter(regex=self.ID)
        if (genedata.empty):
            return
        
        genedata = genedata.loc[(df.iloc[:,0] == gene_ID) | np.array(([True] + [False] * (len(df)-1)))]
        if genedata.isnull().values.any() or len(genedata.index) < 2:
            return
            ''' TODO: temporary,ideally you would want
                to still store the other values
                THIS IS TEMPORARY PLS CHANGE EVENTUALLY'''
        data = TCGA_ExpressionData()
        data.gene_ID = gene_ID
        for val in range(len(genedata.columns)):
            temp = genedata.iloc[:,val]
            data.ExpressionValues[temp.iloc[0]] = temp.iloc[-1]
            
        self.ExpressionData[gene_ID] = data
    
    def GetExpressionValue(self, gene_ID, expression):
        if self.CheckGeneDict(gene_ID, expression):
            return self.ExpressionData[gene_ID].ExpressionValues[expression]
        return None

    def CheckGeneDict(self, gene_ID, expressiontype):
        if gene_ID in self.ExpressionData.keys():
            if self.ExpressionData[gene_ID].CheckExpressionType(expressiontype):
                return True
        return False

    def GetClinicalFeature(self, feature):
        if feature in self.ClinicalFeatures.keys():
            return self.ClinicalFeatures[feature]

        return None
TCGA_Patients = []
GeneList = [] # list of all the stored genes
ProcessedDataBuffer = []
GeneDict = {}

def CheckPatientDuplicates(TGCA_Barcode):
    # checks for patient duplicates in array
    for value in TCGA_Patients:
        if (value.ID == GetPatientID(TGCA_Barcode)):
            return True
    return False

# create TCGA patients off of a dataset
# TODO: other functions to add or remove patients
def CreateTCGAPatients(df):
    # for now just use second row from top for names
    #namelist = df.loc[[0]].squeeze()
    #namelist = namelist.iloc[3:]
    namelist = df.columns[2:]
    
    for value in namelist:
        if (not CheckPatientDuplicates(value)):
            TCGA_Patients.append( TCGA_Patient( GetPatientID(value) ) )

def GetPatient(ID):
    for value in TCGA_Patients:
        if (value.ID == ID):
            return value
    return None # TODO: this will cause an error eventually, keep sharp eyes

def GetPatientID(TCGA_Barcode):
    return TCGA_Barcode[TCGA_PATIENT_ID_START:TCGA_PATIENT_ID_END]    

def StoreGeneForAllPatients(df, gene_ID):
    for value in TCGA_Patients:
        value.StoreExpression(df, gene_ID)
    
    GeneList.append(gene_ID)

# I don't like this method but atp I don't have a choice :/
def __CreateGeneListFromPatients():
    global GeneList

    GeneList.clear()

    with alive_bar(len(TCGA_Patients)) as bar:
        for pt in TCGA_Patients:
            for val in pt.ExpressionData.keys():
                if val not in GeneList:
                    GeneList.append(val)
            bar()

# this needs to be done patient to patient because we modify the patient class
PatientJobQueue = []
TotalJobs = 0
JobsDone = 0
UpdatedJobs = 0

GlobalLock = mtp.Lock()
TargetDataFrame = None

def GetNextJob():
    return PatientJobQueue.get()

def UpdateJobsProgress():
    GlobalLock.acquire()
    UpdatedJobs += 1
    GlobalLock.release()

def ProcessJob(ID):
    try:
        print("Starting process")
    
        print(ID)
        patient = GetPatient( ID )
        for value in df.iloc[1:,0]:
            patient.StoreExpression(df, value)
    except Exception as e:
        traceback.print_exc()
    
 #   UpdateJobsProgress()

def UpdateProgressBar(barfunc):
    GlobalLock.acquire()
    while UpdatedJobs > 0:
        barfunc()
        JobsDone += 1
        UpdatedJobs -= 1
    GlobalLock.release()

def ProgressBarProcess(): # todo: this needs to be a pipe
    print("progress bar starting")
    with alive_bar(TotalJobs) as bar:
        while TotalJobs > JobsDone:
            time.sleep(0.4)
            UpdateProgressBar(bar)

NumOfProcesses = 4
def StoreEverything(df):
    with alive_bar(len(df.iloc[1:,0])) as bar:
        for value in df.iloc[1:,0]:
            StoreGeneForAllPatients(df,value)
            bar()

def SetDesiredProcesses(num):
    NumOfProcesses = num

def PrepareJobQueue(df):
    global TotalJobs
    global JobsDone
    global UpdatedJobs
    global TargetDataFrame
    
    for value in df.columns:
        PatientJobQueue.append( GetPatientID(value) )
    print("Job Queue progress check")
    
    TotalJobs = len(PatientJobQueue)
    TargetDataFrame = df

def PrintPatients():
    for value in TCGA_Patients:
        print(value.ID)


def RestorePatientDatas():
    global TCGA_Patients

    TCGA_Patients.clear()

    root_win = Tk()
    root_win.withdraw()

    filenames = filedialog.askopenfilenames()

    tempTCGA = []
    tempGeneList = []
    tempProcessedData = []
    
    for file in filenames:
        with open(file, 'rb') as f:
            tempTCGA, tempGeneList, tempProcessedData = pickle.load(f)
            TCGA_Patients += tempTCGA

def ImportTCGAData():
    filename = filedialog.askopenfilename()
    
    df = pd.read_csv(filename)
    df.columns = df.iloc[0]
    df = df[1:]
    
    return df

def IsGeneStored(gene_ID):
    #returns true if a gene is in gene list
    if gene_ID in GeneList:
        return True
    
    return False

def SaveLegacy():
    root_win = tk.Tk() # TODO: fix so only need to create root window once
    root_win.withdraw()

    filename = filedialog.askopenfilename()
    with open(filename, 'wb') as f:
        pickle.dump((TCGA_Patients, GeneList, ProcessedDataBuffer), f)

def RestoreLegacy():
    global GeneList
    global ProcessedDataBuffer
    global TCGA_Patients
    
    GeneList.clear()
    TCGA_Patients.clear()
    ProcessedDataBuffer.clear()

    root_win = Tk()
    root_win.withdraw()

    filename = filedialog.askopenfilename()
    
    f = open(filename, 'rb')
    TCGA_Patients, GeneList, ProcessedDataBuffer = pickle.load(f)

# TODO: change to work as zip file
def SavePatientData():
    root_win = tk.Tk() # TODO: fix so only need to create root window once
    root_win.withdraw()

    filename = filedialog.asksaveasfilename()
    
    with open(filename, 'wb') as f:
        pickle.dump((TCGA_Patients, GeneList), f)

def RestorePatientData():
    global TCGA_Patients
    
    GeneList.clear()
    TCGA_Patients.clear()
    ProcessedDataBuffer.clear()

    root_win = Tk()
    root_win.withdraw()

    filename = filedialog.askopenfilename()
    
    f = open(filename, 'rb')
    TCGA_Patients = pickle.load(f)

def SaveProcessedData():
    root_win = tk.Tk() # TODO: fix so only need to create root window once
    root_win.withdraw()

    filename = filedialog.asksaveasfilename()
    
    with open(filename, 'wb') as f:
        pickle.dump((GeneList, ProcessedDataBuffer), f)

def RestoreProcessedData():
    global ProcessedDataBuffer
    global GeneList
    
    GeneList.clear()
    ProcessedDataBuffer.clear()

    root_win = Tk()
    root_win.withdraw()

    filename = filedialog.askopenfilename()
    
    f = open(filename, 'rb')
    GeneList, ProcessedDataBuffer = pickle.load(f)


def RestoreGeneDict():
    global GeneDict

    root_win = Tk()
    root_win.withdraw()
    
    filename = filedialog.askopenfilename()

    f = open(filename, 'rb')
    GeneDict = pickle.load(f)

def SaveGeneDict():
    global GeneDict

    root_win = Tk()
    root_win.withdraw()

    filename = filedialog.asksaveasfilename()

    f = open(filename, 'wb')
    pickle.write(GeneDict, f)

def SaveTCGAPatients():
    global TCGA_Patients
    
    root_win = Tk()
    root_win.withdraw()
    
    filename = filedialog.asksaveasfilename()

    f = open(filename, 'wb')
    pickle.dump(TCGA_Patients, f)

def RestoreTCGAPatients():
    global TCGA_Patients
    
    root_win = Tk()
    root_win.withdraw()
    
    filename = filedialog.askopenfilename()

    f = open(filename, 'rb')
    TCGA_Patients = pickle.load(f)

class SaveDataClass:
    def __init__(self, genelist, processedbuff, patients):
        self.genelist = genelist
        self.processedbuff = processedbuff
        self.patients = patients

# for development in interpreter
def GetSaveData():
    global GeneList
    global ProcessedDataBuffer
    global TCGA_Patients

    return deepcopy( SaveDataClass(GeneList, ProcessedDataBuffer, TCGA_Patients) )

def RestoreSaveData(data):
    global GeneList
    global ProcessedDataBuffer
    global TCGA_Patients

    GeneList = data.genelist
    ProcessedDataBuffer = data.processedbuff
    TCGA_Patients = data.patients

# workspace is a compressed cluster of files that contains data on patients, processed regression data, and other things that can be retrieved using the component system class
#def OpenWorkSpace():
    # import gene list and 

# there is a lot of data to save and restore, so this component system is needed to make things easier to manage

# functor class for save and restore functions
@dataclass
class SaveRestore:
    def __init(self, savefn, restorefn):
        # save and restore should be functions
        self.savefn = savefn
        self.restorefn = restorefn

#ComponentList = {'Patient Data' : SaveRestore(), 'Gene List', 'Processed Data', 'miR Targets'}
class ComponentSystem:
    def __init__(self):
        self.Components = []

#    def Restore(self, component):
        
    
# costly, maybe include this in save file
def GetAllGenes():
    tempgenelist = []
    with alive_bar(len(TCGA_Patients)) as bar:
        for value in TCGA_Patients:
            for genekey in value.ExpressionData: # find a way to vectorize this
                if genekey not in tempgenelist:
                    tempgenelist.append(genekey)
            bar()
    return tempgenelist

def CreateGeneList():
    global GeneList
    GeneList = GetAllGenes()

def CreateGeneListFromDF(df):
    for val in df.iloc[1:0,0]:
        if val not in GeneList:
            GeneList.append(val)

def TranslateGenes(genelist):
    # converts gene names into readable gene names
    for i in range(len(genelist)):
        genelist[i] = GeneDict[genelist[i]]

class ExpressionData:
    def __init__(self, Igene, Iexpressiontype, Dgene, Dexpressiontype):
        self.IGene = Igene
        self.IGeneEtype = Iexpressiontype
        self.DGene = Dgene
        self.DGeneEtype = Dexpressiontype
        self.IValues = []
        self.DValues = []

        self.CorCoeff = 0.0
        self.pval = 0.0
        self.slope = 0.0
        self.intercept = 0.0

        # number of matches against miR prediction algorithms, etc.
        self.metadata = {}
        self.metadata[MIR_PREDICTION_ALGORITHM] = []

    def GetMetadata(self, key):
        if key in self.metadata.keys():
            return self.metadata[key]

        return None

    def ProcessData(self):
        # go through each patient, calculate the log of the independent and dependent genes and store them
        for value in TCGA_Patients:
            '''#todo: figure out a better way to do this (group genes based off of the way their expression is quantified perhaps?)
            if value.CheckGeneDict(self.IGene, self.IGeneEtype) and value.CheckGeneDict(self.DGene, self.DGeneEtype):
                Itemp = float(value.GetExpressionValue(self.IGene, self.IGeneEtype))
                Dtemp = float(value.GetExpressionValue(self.DGene, self.DGeneEtype))
                if not Itemp == 0 and not Dtemp == 0:
                    self.IValues.append(np.log2(Itemp))
                    self.DValues.append(np.log2(Dtemp))'''
            self.ProcessPatient(value)

        if (len(self.IValues) == 0 or len(self.DValues) == 0):
            return
        
        tempresult = scipy.stats.linregress(self.IValues, self.DValues)
        self.CorCoeff = tempresult.rvalue
        self.pval = tempresult.pvalue
        self.slope = tempresult.slope
        self.intercept = tempresult.intercept

    def ProcessDataForClinicalFeature(self, feature, val):
        for value in filter(lambda x: x.GetClinicalFeature(feature) == val, TCGA_Patients):
            self.ProcessPatient(value)

        if (len(self.IValues) == 0 or len(self.DValues) == 0):
            return

        tempresult = scipy.stats.linregress(self.IValues, self.DValues)
        self.CorCoeff = tempresult.rvalue
        self.pval = tempresult.pvalue
        self.slope = tempresult.slope
        self.intercept = tempresult.intercept

        self.metadata[feature] = val
        
    def ProcessPatient(self, patient):
        #todo: figure out a better way to do this (group genes based off of the way their expression is quantified perhaps?)
        if patient.CheckGeneDict(self.IGene, self.IGeneEtype) and patient.CheckGeneDict(self.DGene, self.DGeneEtype):
            Itemp = float(patient.GetExpressionValue(self.IGene, self.IGeneEtype))
            Dtemp = float(patient.GetExpressionValue(self.DGene, self.DGeneEtype))
            if not Itemp == 0 and not Dtemp == 0:
                self.IValues.append(np.log2(Itemp))
                self.DValues.append(np.log2(Dtemp))
                
    
    def PlotData(self):
        fig, ax = plt.subplots()
        
        ax.scatter(self.IValues, self.DValues)

        xlowbound = min(self.IValues)
        xhighbound = max(self.IValues)
        
        xseq = np.linspace(xlowbound, xhighbound, num=100)
        ax.plot(xseq, self.intercept + self.slope * xseq, color="k", lw=2.5)
        plt.xlabel("Log2 of expression of " + self.IGene)
        plt.ylabel("Log2 of expression of " + self.DGene)
        plt.subplots_adjust(top=0.833)

        fig.text(0.3, 0.9, "Correlation Coefficient: " + str(self.CorCoeff))
        fig.text(0.3, 0.87, "P-Value: " + str(self.pval))
        
        fig.show()
        fig.suptitle(str(self.IGene) + " vs. " + str(self.DGene))
        print('---- ' + str(self.IGene) + " vs. " + str(self.DGene) + ' ----')
        print("Correlation Coefficient: " + str(self.CorCoeff))
        print("p-value: " + str(self.pval))
        print("Regression line slope: " + str(self.slope))
        print("Regression line intercept: " + str(self.intercept))

    def PlotDataForEmbedded(self, fig, ax):
        plt.xlabel("Log2 of expression of " + self.IGene)
        plt.ylabel("Log2 of expression of " + self.DGene)
        
        ax.scatter(self.IValues, self.DValues)

        if len(self.IValues) == 0 or len(self.DValues) == 0:
            return
        
        xlowbound = min(self.IValues)
        xhighbound = max(self.IValues)

        fig.suptitle(str(self.IGene) + " vs. " + str(self.DGene))
        xseq = np.linspace(xlowbound, xhighbound, num=100)
        ax.plot(xseq, self.intercept + self.slope * xseq, color="k", lw=2.5)
        plt.subplots_adjust(top=0.833)

        
    def GetSubplots(self):
        fig, ax = plt.subplots()
        
        ax.scatter(self.IValues, self.DValues)
        
        xlowbound = min(self.IValues)
        xhighbound = max(self.IValues)
        
        xseq = np.linspace(xlowbound, xhighbound, num=100)
        ax.plot(xseq, self.intercept + self.slope * xseq, color="k", lw=2.5)
        plt.xlabel("Log2 of expression of " + self.IGene)
        plt.ylabel("Log2 of expression of " + self.DGene)
        plt.subplots_adjust(top=0.833)

        return (fig,ax)

    # getters for sort functions and the like
    @staticmethod
    def GetCC(graph):
        return graph.CorCoeff

    @staticmethod
    def GetPval(graph):
        return graph.pval

    @staticmethod
    def GetSlope(graph):
        return graph.slope

    @staticmethod
    def GetIntercept(graph):
        return graph.intercept

def FastGraphGenes(Igene, Iexpressiontype, Dgene, Dexpressiontype):
    graph = ExpressionData(Igene, Iexpressiontype, Dgene, Dexpressiontype)
    graph.ProcessData()
    graph.PlotData()

# TODO: this is very much not optimal, this functionality should be housed within the ExpressionDataContainer class
def _ExpressionDataVisualizerNavigateRight(container, data, root):
    if data[0] == len(container.Data)-1:
        return
    data[0] = data[0] + 1 # index of expression data, list box data structure, figure, ax, toolbar
        
    #embed matplotlib window
    data[3].clear()
    container.Data[data[0]].PlotDataForEmbedded(data[2], data[3])
    data[2].canvas.draw()

    data[1].delete(0, END)

    data[1].insert(END, "Correlation Coefficient: " + str(container.Data[data[0]].CorCoeff))
    data[1].insert(END, "p-value: " + str(container.Data[data[0]].pval))
    data[1].insert(END, "Regression line slope: " + str(container.Data[data[0]].slope))
    data[1].insert(END, "Regression line intercept: " + str(container.Data[data[0]].intercept))
        
    for key, value in container.Data[data[0]].metadata.items():
        data[1].insert(END, key + ": " + str(value))

    data[1].insert(END, str(data[0]+1) + "/" + str(len(container.Data)))
    
    root.update()
        
def _ExpressionDataVisualizerNavigateLeft(container, data, root):
    if data[0] == 0:
        return
    data[0] = data[0] - 1 # index of expression data, list box data structure, figure, ax, toolbar
        
    #embed matplotlib window
    data[3].clear()
    container.Data[data[0]].PlotDataForEmbedded(data[2], data[3])
    data[2].canvas.draw()
    
    data[1].delete(0, END)

    data[1].insert(END, "Correlation Coefficient: " + str(container.Data[data[0]].CorCoeff))
    data[1].insert(END, "p-value: " + str(container.Data[data[0]].pval))
    data[1].insert(END, "Regression line slope: " + str(container.Data[data[0]].slope))
    data[1].insert(END, "Regression line intercept: " + str(container.Data[data[0]].intercept))
        
    for key, value in container.Data[data[0]].metadata.items():
        data[1].insert(END, key + ": " + str(value))
    
    data[1].insert(END, str(data[0]+1) + "/" + str(len(container.Data)))
    
    root.update()

def Util_WrapText(text, maxlen, font, fontsize):
    font = ImageFont.truetype(font, fontsize)
    index=0
    returnVal=''
    while True:
        temp = text.find(' ', index, len(text)-1)
        if temp == -1:
            if (index <len(text)):
                temp = len(text)-1
                # there is probably one more word
            else:
                break
        
        word = text[index:temp+1]

        tempword = returnVal + word
        tempnl = tempword.rfind('\n')+1
        if font.getlength(tempword[tempnl:]) > maxlen:
            returnVal = returnVal + '\n' + word
        else:
            returnVal = returnVal + word

        index = temp+1

    return returnVal

def _SaveWindow(window, fig, data, filepath):
    # going to have to paste graph onto new image, then write the rest of the metadata onto the library
    FONT = "arial.ttf"
    FONTSIZE = 10
    
    if filepath is None:
        return

    #draw call for tostring_rgb()
    fig.canvas.draw()
    
    buf = fig.canvas.tostring_argb()
    ncols, nrows = fig.canvas.get_width_height()
    img_array = np.fromstring(buf, dtype=np.uint8).reshape(nrows, ncols, 4)

    # save figure size
    figsize = fig.get_size_inches()
    # we want 700wx400h
    fig.set_size_inches(( 700/fig.dpi, 400/fig.dpi ))
    
    CCstr = Util_WrapText("Pearson Correlation Coefficient: " + str(data.CorCoeff), 550, FONT, FONTSIZE)
    pvalstr = Util_WrapText("P-Value: " + str(data.pval), 550, FONT, FONTSIZE)
    slopestr = Util_WrapText("Slope: " + str(data.slope), 550, FONT, FONTSIZE)
    interceptstr = Util_WrapText("Intercept: " + str(data.intercept), 550, FONT, FONTSIZE)

    # ugly ugly ugly way to copy over the figure, but any other method I tried yeilded bit jobbled garbage
    save_image = Image.new("RGBA", (700, 400+500) , "#ffffffff")
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi = fig.dpi)
    buf.seek(0)
    mplimage = deepcopy(Image.open(buf))
    buf.close()
    
    save_image.paste(mplimage, (0,0))

    drawtool = ImageDraw.Draw(save_image) # paste figure onto new image that we will be saving

    font = ImageFont.truetype(FONT, size=FONTSIZE)

    heightoffset = 30
    left, top, right, bottom = drawtool.textbbox((0,0), CCstr, font=font)
    drawtool.text((20, 400 + heightoffset), CCstr, fill=(0,0,0,255))
    heightoffset = heightoffset + bottom+10
    
    left, top, right, bottom = drawtool.textbbox((0,0), pvalstr, font=font)
    drawtool.text((20, 400 + heightoffset), pvalstr, fill=(0,0,0,255))
    heightoffset = heightoffset+ bottom+10
    
    left, top, right, bottom = drawtool.textbbox((0,0), slopestr, font=font)
    drawtool.text((20, 400 + heightoffset), slopestr, fill=(0,0,0,255))
    heightoffset = heightoffset + bottom+10
    
    left, top, right, bottom = drawtool.textbbox((0,0), interceptstr, font=font)
    drawtool.text((20, 400 + heightoffset), interceptstr, fill=(0,0,0,255))
    heightoffset = heightoffset + bottom+10
    
    for key,value in data.metadata.items():
        tempstr = Util_WrapText(str(key) + " : " + str(value), 550, FONT, FONTSIZE)

        left, top, right, bottom = drawtool.textbbox((0,0), tempstr, font=font)
        drawtool.text((20, 400 + heightoffset), tempstr, fill=(0,0,0,255))
        heightoffset = heightoffset + bottom+10

    save_image.save(filepath)
    
    fig.set_size_inches(figsize) # set size back to what we saved

    
def _QuitDataViewer(root):
    root.quit()
    root.destroy()

# contatiner for data so we can have sort functions that arent global and the like
class ExpressionDataContainer:
    def __init__(self):
        self.Data = []
        self.tkwindow = None
        self.metadata = []

    def __init__(self, data):
        self.Data = data
        self.metadata = []
        
    def ImportData(self, data):
        self.Data = data
        self.metadata = data.metadata

    # NOTE: processed data buffer NEEDS to be copied for new data container so we don't accidentally
    # mutate the raw data that we very much want to keep constant so it isn't accidentally
    # saved to a file and the original values arent overwritten
    @staticmethod
    def CreateForClinicalFeature(feature,val):
        return ExpressionDataContainter(filter(lambda x: x.GetMetadata(feature) == val, ProcessedDataBuffer.copy()))

    @staticmethod
    def CreateForGeneList(Igenes, Dgenes):
        return ExpressionDatacontainer(filter(lambda x: x.IGene in Igenes or x.DGene in Dgenes, ProcessedDataBuffer.copy()))
    
    @staticmethod
    def CreateFromDefault():
        return ExpressionDataContainer(ProcessedDataBuffer.copy())
        
    # returns container of sorted data from lowest to highest, opposite if inv is true
    def sort(self, Fn, inv=False):
        if not inv:
            self.Data.sort(key=Fn)
        else:
            self.Data.sort(key=Fn)
            self.Data.reverse()

    def copy(self):
        copy = ExpressionData()
        temp = self.Data.copy()
        copy.ImportData(temp)
        return copy

    def filter(self, Fn):
        # returns a copy
        return ExpressionDataContainer(list(filter(Fn, self.Data)))

    @staticmethod
    def _FilterSmallData(data):
        # filter out data with less than 10 values
        if len(data.IValues) < 10:
            return False
        return True

    def FilterSmallData(self):
        return self.filter(ExpressionDataContainer._FilterSmallData)

    @staticmethod
    def _FilterNSPVal(data):
        #filter out data with a pvalue greater than .05
        if data.pval > .05:
            return False
        return True

    def FilterNSPVal(self):
        return self.filter(ExpressionDataContainer._FilterNSPVal)
    
    @staticmethod
    def _FilterLowCC(data):
        #filter out data with a correlation coefficient -.3 < x <.3
        if abs(data.CorCoeff) < .3:
            return False

        return True

    def FilterLowCC(self):
        return self.filter(ExpressionDataContainer._FilterLowCC)

    def FilterForGenes(Igenes, Dgenes):
        return self.filter(lambda x: x.IGene in Igenes or x.DGene in Dgenes)
    
    def SearchForGene(self, gene):
        for val in self.Data:
            if val.DGene == gene:
                return val

        return None


    #def GetMeanOfMetric(metricfn):

    #def GetMedianOfMetric(metricfn):

    #def GetSDOfMetric(metricfn):

    # assume predarr is sorted by which most likely targets are first, and least likely are last
    def MatchMiRTargetPredData(self, name, predarr):
        for val in self.Data:
            if val.DGene in predarr:
                val.metadata[MIR_PREDICTION_ALGORITHM].append(name)
                val.metadata[str(name + " priority")] = predarr.index(val.DGene)


    def FilterByPredAlgorithm(self):
        return self.filter(lambda x: len(x.metadata[MIR_PREDICTION_ALGORITHM]) > 0)

    def _FilterAndSaveItems(self, child_win, selection):
        directory = filedialog.askdirectory(title="Select Directory", mustexist=True)
        print(directory)
        
        # temp figure for subplots
        fig,ax = plt.subplots()
        
        for val in selection.curselection():
            # get index, which are the characters after the " "
            graphname = selection.get(val)
            str_idx = graphname.rfind(" ")+1

            data_idx = int(graphname[str_idx:])
            graphdata = self.Data[data_idx]

            print(data_idx)
            
            # create file path
            if directory[-1] != '/':
                directory += '/'

            # remove tricky characters and build filename
            temptrans = str.maketrans("", "", "?!.")
            
            filepath = directory + graphname[0:str_idx-1].translate(temptrans) + ".png"
            
            graphdata.PlotDataForEmbedded(fig, ax)

            print(filepath)
            
            _SaveWindow(child_win, fig, graphdata, filepath)
    

    def _SaveVisualDialog(self):
        if self.tkwindow is None:
            return
        
        # create child window
        child_win = Toplevel(self.tkwindow)
        child_win.title("Save Dialog")
        child_win.geometry("600x300")
        
        # create check buttons for each data
        data_strs = [data.IGene + " vs. " + data.DGene + " " + str(index) for index,data in enumerate(self.Data)]
        
        # create scrollbar
        scrollbar = Scrollbar(child_win, orient='vertical')
        scrollbar.pack(side=RIGHT, fill=Y)
        
        listbox = Listbox(child_win, selectmode='multiple', yscrollcommand = scrollbar.set)
        listbox.pack(padx=10,pady=10,expand=YES,fill="both")

        for data_str in data_strs:
            listbox.insert(END, data_str)
        
        
        SaveButton = Button(child_win, text="Save", command=lambda: self._FilterAndSaveItems(child_win, listbox))
        #CancelButton = Button(child_win, text="Cancel", command=lambda: child_win.destroy())
        
        SaveButton.pack(side=BOTTOM, fill=X)
        #CancelButton.pack(side=BOTTOM, fill=X)
        
        child_win.protocol("WM_DELETE_WINDOW", lambda: child_win.destroy())
        
        # display everything
        child_win.mainloop()
        
        

    def VisualizeData(self):
        # click through each data set, which displays a scatterplot with a linear regression, and any other metadata
        root_win = Tk()
        root_win.title("Expression Data")

        self.tkwindow = root_win
        #root_win.geometry("700x1600")
        
        '''becuase there is no explicit referencing/copying in python all variables that are changed
            by event function must be in an array or dictionary so they can be changed... >:('''
        data = [0] # index of expression data, list box data structure, figure, ax, toolbar

        #first embed matplotlib window
        fig, ax = self.Data[0].GetSubplots()
        fig.suptitle(str(self.Data[0].IGene) + " vs. " + str(self.Data[0].DGene))
        canvas = FigureCanvasTkAgg(fig, master=root_win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, root_win)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        scrollbar = Scrollbar(root_win)
        scrollbar.pack( side = RIGHT )

        valuelist = Listbox(root_win, yscrollcommand = scrollbar.set, width=100 )
        valuelist.insert(END, "Correlation Coefficient: " + str(self.Data[0].CorCoeff))
        valuelist.insert(END, "p-value: " + str(self.Data[0].pval))
        valuelist.insert(END, "Regression line slope: " + str(self.Data[0].slope))
        valuelist.insert(END, "Regression line intercept: " + str(self.Data[0].intercept))
        
        for key, value in self.Data[0].metadata.items():
           valuelist.insert(END, key + ": " + str(value))
        
        valuelist.insert(END, str(data[0]+1) + "/" + str(len(self.Data)))
        
        valuelist.pack( side = LEFT, fill = BOTH )
        scrollbar.config( command = valuelist.yview )

        data.append(valuelist)
        data.append(fig)
        data.append(ax)
        data.append(toolbar)

        prevbutton = Button(root_win, text="next", command=lambda: _ExpressionDataVisualizerNavigateRight(self, data, root_win))
        nextbutton = Button(root_win, text="prev", command=lambda: _ExpressionDataVisualizerNavigateLeft(self, data, root_win))
        savebutton = Button(root_win, text="save", command=lambda: self._SaveVisualDialog())
        
        prevbutton.pack(side=RIGHT)
        nextbutton.pack(side=RIGHT)
        savebutton.pack(side=RIGHT)
        
        root_win.protocol("WM_DELETE_WINDOW", lambda: _QuitDataViewer(root_win))
        
        tk.mainloop()
        
def ProcessAllGenesForIGene(IGene, Iexpressiontype, Dexpressiontype):
    global ProcessedDataBuffer
    ProcessedDataBuffer.clear()
    with alive_bar(len(GeneList)) as bar:
        for val in GeneList:
            if IGene == val: continue
            tempgraph = ExpressionData(IGene, Iexpressiontype, val, Dexpressiontype)
            tempgraph.ProcessData()
            ProcessedDataBuffer.append(tempgraph)
            bar()

# returns an array of processed genes seperated by the category in the clinical feature
def ProcessGeneForClinicalFeature(IGene, Iexpressiontype, Dexpressiontype, feature):
    global ProcessedDataBuffer
    ProcessedDataBuffer.clear()
    # first iterate through patient data, gather possible clinical feature values
    featurestates = []
    for pt in TCGA_Patients:
        temp = pt.GetClinicalFeature(feature)
        if temp is None:
            continue

        if type(temp) is not str:
            continue #sometimes pandas will be dumb and replace NA with NAN even when I want NA to mean 'not applicable'
        
        if temp not in featurestates:
            featurestates.append(temp)
    
    # for each possible value, compute expression data of that group
    with alive_bar(len(GeneList) * len(featurestates)) as bar:
        for featureval in featurestates:
            for val in GeneList:
                if IGene is val: continue
                tempgraph = ExpressionData(IGene, Iexpressiontype, val, Dexpressiontype)
                tempgraph.ProcessDataForClinicalFeature(feature, featureval)
                ProcessedDataBuffer.append(tempgraph)
                bar()
    
    
''' TODO: CITE THESE DATABASES '''

# one might want to use some prediction algorithms to narrow the search for miRNA's of interest, so thats what were gonna do
#_PicTarRequestURL = 'https://pictar.mdc-berlin.de/cgi-bin/PicTar_vertebrate.cgi'
'''def _PicTarSubmitRequest(miRNA):
    #TODO: gotta captialize the r in miR or else the request won't recognize it
    
    
    requestpayload = {'clade' : 'vertebrate',
                      'dataset' : 'target predictions for all human microRNAs based on conservation in mammals (human, chimp, mouse, rat, dog)',
                      'name2' : miRNA,
                      'name1' : '',
                      'action' : 'Search for targets of a miRNA'}
    
    r = requests.post(_PicTarRequestURL, data=requestpayload)

    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")
        
    return r.text
    '''
# returns a list of the predicted genes, REQUIRES GeneDict AND GeneList TO BE LOADED INTO MEMORY
def _PicTarGetPredictionSet(miRNA):
    requestpayload = {'clade' : 'vertebrate',
                      'dataset' : 'target predictions for all human microRNAs based on conservation in mammals (human, chimp, mouse, rat, dog)',
                      'name2' : miRNA,
                      'name1' : '',
                      'action' : 'Search for targets of a miRNA'}
    
    r = requests.post(_PicTarRequestURL, data=requestpayload)

    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")
    
    soup = BeautifulSoup(r.text, 'html.parser')
    table = soup.find('table').find_next('table')
    RetList = []
    
    for val in table.contents[3:]:
        #iterate through each entry in the table
        genestr = val.contents[-4].string # this way of indexing is some vodoo stuff i just kept subtracting 1 until it gave me what I wanted lol
        
        openpos = genestr.find('(')
        closepos = genestr.find(')')
        if openpos == -1 or closepos == -1:
            continue
        
        while closepos < len(genestr):
            if openpos > closepos:
                closepos = genestr.find(')', closepos)
                if closepos == -1:
                    temp = None
                    break;

            temp = genestr[openpos+1 : closepos]

            if temp in GeneList:
                break
            
            if temp in GeneDict.keys():
                if GeneDict[temp] in GeneList:
                    temp = GeneDict[temp]
                    break
            
            openpos = genestr.find('(', openpos+1)
            closepos = genestr.find(')', closepos+1)

            if openpos == -1 or closepos == -1:
                temp = None
                break
        
        if temp is None: continue
        
        # the gene is usually in parenthases, so use that as a way to tell for now
        RetList.append(temp)
        
    return RetList

# structure that defines standard group of operations that can be done on a webpage
@dataclass
class PredictionSet:
    def __init__(self, name, predfunc):
        self.name = name
        self.predfunc = predfunc
        
_MirDBRequestURL = 'https://mirdb.org/cgi-bin/search.cgi'
_MirDBBaseURL = 'https://mirdb.org'
'''def _MirDBSubmitRequest(miRNA):
    requestpayload = {'species' : 'Human',
                      'searchBox' : miRNA,
                      'submitButton' : 'Go',
                      'searchType' : 'miRNA'}
    
    r = requests.post(_MirDBRequestURL, data=requestpayload)

    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

    return r.text
'''
def cleanhtml(text):
  cleantext = re.sub(re.compile('<.*?>'), '', text)
  return cleantext

def _MirDBParseWindowFunc(selection, combobox, root):
    selection.append(combobox.get())
    root.quit()

def _MirDBGetPredictionSet(miRNA):
    requestpayload = {'species' : 'Human',
                      'searchBox' : miRNA,
                      'submitButton' : 'Go',
                      'searchType' : 'miRNA'}
    
    r = requests.post(_MirDBRequestURL, data=requestpayload)

    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

    soup = BeautifulSoup(r.text, 'html.parser')
    # first check if it directed us to results or another page letting us specify by checking if there is an h3
    if soup.find('h3'):
       # yep, multiple entries, make a window popup and ask the user which one they want
       
        temp = soup.find(href=re.compile('searchType=miRNA'))
        mirlist = temp.find_all('h4')
        for i in range(len(mirlist)):
            mirlist[i] = cleanhtml(mirlist[i].string)
        
        selection = [] # have to use an array because python doesn't have explicit references :(
        root_win = Tk()

        label = Label(root_win, text='Multiple genes have been returned from search, please select desired gene')
        label.pack()
        combobox = Combobox(root_win, state="readonly", values=mirlist)
        combobox.pack()
        button = Button(root_win, text='OK', command=lambda: _MirDBParseWindowFunc(selection, combobox, root_win) )
        button.pack()
        
        root_win.mainloop()

        root_win.destroy()
        # now make the request using the selection
        urltag = temp.find(string=selection[0]).parent.parent
        r = requests.get(_MirDBBaseURL + urltag['href'])
        if not r.status_code == 200:
            warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

        soup = BeautifulSoup(r.text, 'html.parser')

    table = soup.find('table', id='table1')
    RetList = []
    
    for val in table.contents[2:]:
        #iterate through each entry in the table
        # every other value is '\n', but its not a big deal
        if val.string == '\n':
            continue
        
        if not val.find('a', target='_blank'):
            continue
        
        temp = val.find('a', target='_blank').string
        temp = temp.strip()

        if temp in GeneList:
            RetList.append(temp)

        if temp in GeneDict.keys():
            if GeneDict[temp] in GeneList:        
                RetList.append(GeneDict[temp])
        
    return RetList

_TargetScanRequestURL = 'https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi?'
_TargetScanBaseURL = 'https://www.targetscan.org'
'''def _TargetScanSubmitRequest(miRNA):
    payload = {'species' : 'Human',
               'gid' : '',
               'mir_sc' : '',
               'mir_nc' : '',
               'mir_vnc' : '',
               'mirg' : miRNA}

    # GET request >:(
    r = requests.get(_TargetScanRequestURL + urllib.parse.urlencode(payload))
    
    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

    return r.text
'''
def _TargetScanGetPredictionSet(miRNA):
    payload = {'species' : 'Human',
               'gid' : '',
               'mir_sc' : '',
               'mir_nc' : '',
               'mir_vnc' : '',
               'mirg' : miRNA}

    # GET request >:(
    r = requests.get(_TargetScanRequestURL + urllib.parse.urlencode(payload))
    
    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")
    
    soup = BeautifulSoup(r.text, 'html.parser')
    # first check if it directed us to results or another page letting us specify by checking if there is an h3
    if soup.find(string=re.compile('matches multiple families in our miRNA database')):
        # yep, multiple entries, make a window popup and ask the user which one they want
        templist = soup.find('h3').find('table').find_all('a')
        mirlist = []
        for val in templist:
            mirlist.append(val.string)
        
        selection = [] # have to use an array because python doesn't have explicit references :(
        root_win = Tk()

        label = Label(root_win, text='Multiple genes have been returned from search, please select desired gene')
        label.pack()
        combobox = Combobox(root_win, state="readonly", values=mirlist)
        combobox.pack()
        button = Button(root_win, text='OK', command=lambda: _MirDBParseWindowFunc(selection, combobox, root_win) )
        button.pack()
        
        root_win.mainloop()

        root_win.destroy()
        # now make the request using the selection
        urltag = soup.find('h3').find('a',string=selection[0])
        
        r = requests.get(_TargetScanBaseURL + urltag['href'])
        if not r.status_code == 200:
            warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")
        soup = BeautifulSoup(r.text, 'html.parser')

    genetags = soup.find('table', id='restable').find_all('a',target='new')
    RetList = []

    for val in genetags:
        #iterate through each entry in the table
        # every other value is '\n', but its not a big deal
        if val.string == '\n':
            continue

        temp = val.string

        if temp in GeneList:
            RetList.append(temp)

        if temp in GeneDict.keys():
            if GeneDict[temp] in GeneList:
                RetList.append(GeneDict[temp])
        
    return RetList

_MirWalkRequestURL = 'http://mirwalk.umm.uni-heidelberg.de/search/?'
'''def _MirWalkSubmitRequest(miRNA):
    payload = {'species' : 'human', 'gene' : '', 'mirna' : miRNA}

    r = requests.get(_MirWalkRequestURL + urllib.parse.urlencode(payload))

    if not r.status_code == 200:
        warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

    return r.text
'''
def _MirWalkGetPredictionSet(miRNA):
    with request.Session() as session:
        payload = {'species' : 'human', 'gene' : '', 'mirna' : miRNA}
        r = requests.get(_MirWalkRequestURL + urllib.parse.urlencode(payload))

        if not r.status_code == 200:
            warnings.warn("Status code of the returned data is not OK! Data may be incorrect or corrupt!")

        # create session up here
        soup = BeautifulSoup(r.text, 'html.parser')

        if soup.find(string=re.compile('was not found in our database')):
            templist = soup.find(class_='table-container').find_all('tr')
            mirlist = []
            saveurl = {}
            for val in templist[1:]:
                temp = val.find_all('a')[1]
                mirlist.append(temp.string)
                saveurl[mirlist[-1]] = temp['href']
                
            selection = [] # have to use an array because python doesn't have explicit references :(
            root_win = Tk()
                
            label = Label(root_win, text='Multiple genes have been returned from search, please select desired gene')
            label.pack()
            combobox = Combobox(root_win, state="readonly", values=mirlist)
            combobox.pack()
            button = Button(root_win, text='OK', command=lambda: _MirDBParseWindowFunc(selection, combobox, root_win) )
            button.pack()
                
            root_win.mainloop()
                
            root_win.destroy()
            # now make the request using the selection
            
        r = session.get(saveurl[selection[0]])
        data = b''
        with session.get('http://mirwalk.umm.uni-heidelberg.de/export/csv/', stream=True) as response:
            with alive_bar(int(response.headers.get('Content-Length'))/4096) as bar:
                for chunk in response.iter_content(chunk_size=4096):
                    if chunk:
                        data = data + chunk
                        bar()
        
        
        df = pd.read_csv(StringIO(data.encode('UTF-8')))

        df = df[ df['bindingp'] > 0.95]
        
        return df['mirnaid']

class MirandaGeneData:
    def __init__(self):
        self.MirName = ""
        self.TargetName = ""
        self.TotalScore = ""
        self.TotalEnergy = ""

def ReadMirandaOutFile():
    file = filedialog.askopenfile(mode="r")

    if file is None:
        return None
    
    data = file.read().replace('\n', '')
    returnVal = []
    
    barpos=0
    fstart = 0
    fend = len(data)-1
    with alive_bar(fend) as bar:
        while True:
            mirdata = MirandaGeneData()
            
            index = data.find('>>', fstart, fend)
            if index == -1:
                while barpos<fend:
                    bar()
                    barpos = barpos+1
                    
                break
            
            index = index+2 # set position to start of next name

            #temp = data.find('\t', index, fend)
            mirdata.MirName = data[index:temp]

            index = temp+1
            temp = data.find('\t', index, fend)
            mirdata.TargetName = data[index:temp]

            index = temp+1
            temp = data.find('\t', index, fend)
            mirdata.TotalScore = data[index:temp]

            index = temp+1
            temp = data.find('\t', index, fend)
            mirdata.TotalEnergy = data[index:temp]

            returnVal.append(mirdata)

            fstart = temp
            
            while barpos<fstart:
                bar()
                barpos = barpos+1
    
    # sort based on score energy
    returnVal.sort(key=lambda x : x.TotalScore, reverse=True)
    return returnVal

RefseqConversionTable = None
def LoadRefseqConversionTable():
    global RefseqConversionTable
    
    file = filedialog.askopenfilename()
    RefseqConversionTable = pd.read_csv(file)

def GetGeneIDFromRefseq(refseq):
    if RefseqConversionTable is None:
        print("Please load conversion table")
        return None
    try:
        return RefseqConversionTable.loc[RefseqConversionTable['name'] == refseq]['name2'].values[0]
    except:
        print("it seems that either the refseq ID is either not in the table, or the table is not formatted correctly")
        return None

def _MirandaGetPredictionSet():
    mirandadata = ReadMirandaOutFile()
    strdata = []
    
    if mirandadata is None:
        return None

    if RefseqConversionTable is None:
        print("please load refseq conversion table")
        return None

    for element in mirandadata:
        element.TargetName = GetGeneIDFromRefseq(element.TargetName)
        strdata.append(element.TargetName)

    return mirandadata, strdata

genecopoia_URL = "https://www.genecopoeia.com/product/search/result.php?pageNum_Recordset1=0&field=11&tax_id=10090&key=mmu-miR-196b-5p&prt=16&totalRows_Recordset1=252"
genecopoia_baseURL = "https://www.genecopoeia.com"

def _AuxGenecopiaFindPageIndex(pages, index):
    for page in pages:
        try:
            if (int(page.text) == index): return page
        except ValueError as ve:
            continue

# these targets are for mouse, no simple way to navigate site so just provide URL of target query results
def _GetGenecopiaMMTargets():
    if GeneDict == {}:
        print("The gene dictionary is not loaded!")
        return
    
    filename = filedialog.askopenfilename()
    cookies = None
    returnVal = []

    with open(filename, 'rb') as f:
        cookies = pickle.load(f)

    c_jar = requests.cookies.RequestsCookieJar()
    
    for cookie in cookies:
        c_jar.set(cookie['name'], cookie['value'], domain=cookie['domain'], path=cookie['path'])
    
    r = requests.get(genecopoia_URL,cookies=c_jar)

    soup = BeautifulSoup(r.text, 'html.parser')

    pages = soup.find(class_='pagination').find_all('li')
    
    # find total number of pages
    numpages = int(pages[-1].text)
    
    for i in range(2, numpages+1):
    # for each page (including starting page)
        
        #   get body of results
        table = soup.find(class_="table_org dcf-table dcf-table-responsive dcf-table-bordered dcf-table-striped dcf-w-100%").find_all('tr')
        pages = soup.find(class_='pagination').find_all('li')

        # skip first row
        for idx in range(1, len(table)):
            # get all RNA names given (data-label: "Symbol"), or position 3 in .contents
            gene = table[idx].contents[3].text.upper()
            
            if gene in GeneDict.keys():
                if GeneDict[gene] in GeneList:        
                    returnVal.append(GeneDict[gene])
        
            
        #   navigate to next page
        link = genecopoia_baseURL + _AuxGenecopiaFindPageIndex(pages, i).a['href']
        
        #   repeat
        soup = BeautifulSoup(requests.get(link, cookies=c_jar).text, 'html.parser')


    return returnVal

# keep in mind these values will be saved as strings
def LoadFirebrowseClinicalFeatures():
    path = filedialog.askopenfilename()
    
    # is the file excel or tsv?
    df = pd.read_csv(path, sep='\t', keep_default_na=False)
    # create dataframe from excel file (or tsv)
    global TCGA_Patients

    # for each patient (TCGA-##-####)
    #   append each feature to clinical feature list yes, no, N/A <-- frustrating :(
    for pt in TCGA_Patients:
        try:
            temprow = df[pt.ID.lower()]
        except:
            continue

        # set colum headers to their descriptors
        temprow.index = df['bcr_patient_barcode']

        for value in temprow.index:
            pt.AppendClinicalFeature(value, temprow[value])

def UlcerationGraph(Igene, stratification_metric=None):
    # crete matplotib bar graph
    fig, ax = plt.subplots()

    ulceration_stats = ['Ulcerated', 'Non-ulcerated']

    colors = [(0, 0, 1, 1), (1, 0, 0, 1)]
    if stratification_metric is None:
        values = [ [],[] ]
        for pt in TCGA_Patients:
            if pt.ClinicalFeatures['melanoma_ulceration_indicator'] == 'yes':
                if pt.CheckGeneDict(Igene, 'reads_per_million_miRNA_mapped'):
                    values[0].append( float(pt.GetExpressionValue(Igene, 'reads_per_million_miRNA_mapped')) )
            elif pt.ClinicalFeatures['melanoma_ulceration_indicator'] == 'no':
                if pt.CheckGeneDict(Igene, 'reads_per_million_miRNA_mapped'):
                    values[1].append( float(pt.GetExpressionValue(Igene, 'reads_per_million_miRNA_mapped')) )
        
        err = [scipy.stats.sem(values[0]), scipy.stats.sem(values[1])]

        counts = [sum(values[0]) / len(values[0]), sum(values[1]) / len(values[1])]
        
        # normalize to fold change
        nzf_rdpm = counts[1]

        counts[0] = np.log2(counts[1] / nzf_rdpm)
        counts[1] = np.log2(1.0)

        ttestbefore = scipy.stats.ttest_ind(values[0], values[1])
        
        for index,val in enumerate(values[0]):
            values[0][index] = np.log2(val/nzf_rdpm)

        for index,val in enumerate(values[1]):
            values[1][index] = np.log2(val/nzf_rdpm)

        # calculate standard error of mean
        err = [scipy.stats.sem(values[0]), scipy.stats.sem(values[1])]
        
        ax.bar(ulceration_stats, counts, width=0.8)
        plt.errorbar(ulceration_stats, [counts[0], counts[1]], yerr=err, fmt="o", color="g")
        ax.set_ylabel('reads_per_million_miRNA_mapped')
        ax.set_title('196b expression in ulcerated vs non ulcerated Melanoma')

        for i in range(len(ulceration_stats)):
             ax.scatter(i + np.random.random(len(values[i])) * .8 - .8 / 2, values[i], color=colors[i])

        ttestafter = scipy.stats.ttest_ind(values[0], values[1])
    else:
        print('')
    
    plt.show()
# batch normalization
def CreateUlcerationWorkbook(mirs, readtypes):
    savename = filedialog.asksaveasfilename(defaultextension='.xlsx', filetypes=(("excel file", "*.xlsx"),("All Files", "*.*") ))
    
    #create excel workbook
    wb = openpyxl.Workbook()
    # create ulcerated and non ulcerated sheets
    ulceratedws = wb.create_sheet('Ulcerated')
    nonulceratedws = wb.create_sheet('Non-Ulcerated')
    
    # this is kinda backwards... but it should work
    # interate through row of each column for each patient and segregate data
    for index,row in enumerate(ulceratedws.iter_rows(min_row=1, max_row=len(TCGA_Patients)+2, max_col=len(mirs) * len(readtypes)+1)):
        if index == 0:# title row 1, fill with gene names
            for gindex,gene in enumerate(mirs):
                row[gindex * len(readtypes) + 1].value = gene
            continue
        elif index == 1:
            for i in range(1,len(row)):
                row[i].value = readtypes[(i-1) % len(readtypes)]
            continue
        
        pt = TCGA_Patients[index-2]
        row[0].value = pt.ID # set index
        
        if pt.GetClinicalFeature('melanoma_ulceration_indicator') == 'yes':
            for gindex,gene in enumerate(mirs):
                for eindex,expression in enumerate(readtypes):
                    exprval = pt.GetExpressionValue(gene, expression)
                    if exprval is not None:
                        row[gindex*len(readtypes) + eindex + 1].value = str(np.log2(float(pt.GetExpressionValue(gene, expression))))
    for index,row in enumerate(nonulceratedws.iter_rows(min_row=1, max_row=len(TCGA_Patients)+2, max_col=len(mirs) * len(readtypes)+1)):
        if index == 0:# title row 1, fill with gene names
            for gindex,gene in enumerate(mirs):
                row[gindex * len(readtypes) + 1].value = gene
            continue
        elif index == 1:
            for i in range(1,len(row)):
                row[i].value = readtypes[(i-1) % len(readtypes)]
            continue
        
        pt = TCGA_Patients[index-2]
        row[0].value = pt.ID # set index
        
        if pt.GetClinicalFeature('melanoma_ulceration_indicator') == 'no':
            for gindex,gene in enumerate(mirs):
                for eindex,expression in enumerate(readtypes):
                    exprval = pt.GetExpressionValue(gene, expression)
                    if exprval is not None:
                        row[gindex*len(readtypes) + eindex + 1].value = str(np.log2(float(pt.GetExpressionValue(gene, expression))))

    wb.save(savename)

# visualize slope data against prediction algorithms too!
def GenerateUlcerationHeatmap(expressioncontainer):
    # set up dataframe such that column headers are "ulcerated", "non-ulcerated", "N/A"
    # and row headers are gene names
    global GeneList
    global ProcessedDataBuffer

    clinicalfeaturemap = {'no' : 'non-ulcerated', 'yes' : 'ulcerated'}
    
    ulcerationdata = {'ulcerated' : [0.] * len(GeneList), 'non-ulcerated' : [0.] * len(GeneList)}
    df = pd.DataFrame(ulcerationdata, index=GeneList.copy())

    # for each entry in the processed data
        # seggregate the slope values into the dataframe
    for val in ProcessedDataBuffer:
        if len(val.IValues) < 10 or val.DGene not in df.index:
            continue

        cf = val.GetMetadata('melanoma_ulceration_indicator')

        if cf not in clinicalfeaturemap.keys():
            continue

        if val.pval > 0.1 or len(val.IValues) < 30:
            continue
        
        df.at[val.DGene, clinicalfeaturemap[cf]] = val.slope

    # scale to unit variance (SD)
    #scaled_data = StandardScaler().fit_transform(df.transpose())
    #df = pd.DataFrame(scaled_data.transpose(), index = df.index, columns = df.columns)
    # Normalize each value to the non ulcerated and take the log fold change
    for index,row in df.iterrows():
        '''tempU = row['ulcerated']
        tempNU = row['non-ulcerated']
        try:
            row['ulcerated'] = np.log2(tempU / tempNU)
            row['non-ulcerated'] = 1
        except:
            df.drop(index)
        if np.isnan(row['ulcerated']):
            df = df.drop(index)'''
        if abs(row['ulcerated'] - row['non-ulcerated']) < 0.4:
            df = df.drop(index)

    print(df.index.tolist())
    # visualize with heatmap
    # (only ulcerated because we know non ulcerated is 1)
    seaborn.clustermap(df, z_score=None, standard_scale=None, yticklabels=True)
    plt.show()


tm_strings = {'DSS': ['DSS_cr', 'DSS.time.cr'], 'PFI': ['PFI.cr', 'PFI.time.cr']}
def SurvivalCurve196(percentile_low, percentile_high, timemetric):
    global TCGA_Patients
    stagefeaturename = 'pathologic_stage'
    stages = ['stage 0', 'stage i', 'stage ia', 'stage ib', 'stage ic', 'i/ii nos', 'stage ii', 'stage iia', 'stage iib', 'stage iic', 'stage iii', 'stage iiia', 'stage iiib', 'stage iiic', 'stage iv']
    
    # load survival data
    survival_df = pd.read_csv(filedialog.askopenfilename(defaultextension='.csv', filetypes=(("comma seperated values file", "*.csv"),("All Files", "*.*") )))
    
    # build matrix from survival data
    # get SKCM
    skcm_df = survival_df[ survival_df['type'] == 'SKCM']

    # get mir-196b expression for percentile seggregation
    expression_196 = [float(val.GetExpressionValue('hsa-mir-196b', 'reads_per_million_miRNA_mapped')) for val in TCGA_Patients if val is not None and val.GetExpressionValue('hsa-mir-196b', 'reads_per_million_miRNA_mapped') is not None]

    # load into data matrix
    percentile_196 = {x : [] for x in stages}
    pt_status = {x : [] for x in stages}
    pfi_time = {x : [] for x in stages}

    val_high = np.percentile(expression_196, percentile_high)
    val_low = np.percentile(expression_196, percentile_low)    

    len(percentile_196)
    
    for val in skcm_df['bcr_patient_barcode']:
        pt = GetPatient(val)

        if pt is None:
            continue

        stage = pt.GetClinicalFeature('pathologic_stage')

        if stage not in stages:
            continue
        
        masked = skcm_df[skcm_df['bcr_patient_barcode'] == val]
        tempstatus = masked[tm_strings[timemetric][0]].array[0]
        temptime = masked[tm_strings[timemetric][1]].array[0]

        if tempstatus == '#N/A':
            continue
        elif temptime == '#N/A' or temptime == 0:
            continue

        tempexp = pt.GetExpressionValue('hsa-mir-196b', 'reads_per_million_miRNA_mapped')
        if tempexp is None:
            continue
        
        if float(tempexp) > val_high:
            percentile_196[stage].append('hi')
        elif float(tempexp) < val_low:
            percentile_196[stage].append('lo')
        else:
            continue

        if tempstatus == 1:
            pt_status[stage].append(1)
        else:
            pt_status[stage].append(0)

        pfi_time[stage].append(float(temptime))
        if np.isnan(temptime): # for some reason a NaN slips through when all of the entries are N/A, temporary workaround
            percentile_196[stage].pop()
            pt_status[stage].pop()
            pfi_time[stage].pop()

    skcm_dfs = []
    
    for stage in stages:
        skcm_dfs.append(pd.DataFrame({'mir-196' : percentile_196[stage],
                            'Status' : pt_status[stage],
                            'Time': pfi_time[stage]}))
    
    
    '''maskhi = skcm_df['mir-196'] == 'hi'
    masklo = skcm_df['mir-196'] == 'lo'
    
    # load into learner of choice
    time_cell, survival_prob_cell, conf_int = kaplan_meier_curve(
        skcm_df['Status'][maskhi], skcm_df['Time'][maskhi])
    
    plt.step(time_cell, survival_prob_cell, where="post", label=f"hi (n = {maskhi.sum()})")
    
    time_cell, survival_prob_cell, conf_int = kaplan_meier_curve(
        skcm_df['Status'][masklo], skcm_df['Time'][masklo])
    
    plt.step(time_cell, survival_prob_cell, where="post", label=f"hi (n = {maskhi.sum()})")


    plt.ylim(0, 1)
    plt.ylabel("est. probability of survival")
    plt.xlabel("time")
    plt.legend(loc="best")'''
    fig,axes = plt.subplots(ncols=3,nrows=5)
    kmf = [KaplanMeierFitter()] * len(skcm_dfs)
    print("=" * 30)
    print(f"SURVIVAL CURVE {percentile_low}-{percentile_high}")
    print('='*30)
    for idx,df in enumerate(skcm_dfs):
        hi_df = df[df['mir-196'] == 'hi']
        lo_df = df[df['mir-196'] == 'lo']
        
        if hi_df.empty or lo_df.empty:
            continue
        
        kmf[idx].fit(durations=hi_df['Time'], event_observed=hi_df['Status'], label=f'hi n={len(hi_df.index)}')
        kmf[idx].plot_survival_function(ax=axes[idx // 3, idx % 3])
        kmf[idx].fit(durations=lo_df['Time'], event_observed=lo_df['Status'], label=f'lo n={len(lo_df.index)}')
        kmf[idx].plot_survival_function(ax=axes[idx // 3, idx % 3])

        results = logrank_test(durations_A=hi_df['Time'], durations_B=lo_df['Time'], event_observed_A=hi_df['Status'], event_observed_B=lo_df['Status'])
        print(f"\nlog rank of group: {stages[idx]}")
        results.print_summary()
        print('=' * 30)
        
        axes[idx // 3, idx % 3].set_title(f'Kaplan-Meier Estimates by Group, {stages[idx]}')
        axes[idx // 3, idx % 3].set_xlabel('Time')
        axes[idx // 3, idx % 3].set_ylabel('Survival Probability')
        axes[idx // 3, idx % 3].legend()

    plt.subplots_adjust(hspace=1.5)
    
    fig.show()
    
# Pa = A(A^T * A)^(-1) * A^T
# eventually do PCA
#def ProjectPatientData():
    # set up colum vectors of patient gene expression


#def UMAPPatientData(expressiontype):
    # create pt data matrix using pandas dataframes
#    for pt in TCGA_Patients:
        

        
    # remove columns with zero or leave it

#TODO: target 
def _psRNATargetGetPredictionSet():
    df = pd.read_csv(filedialog.askopenfilename(defaultextension='.txt', filetypes=(("psRNATarget text output", "*.txt"),("All Files", "*.*"))), sep='\t')
    temp = df['Target_Acc.'].tolist()

    returnval = []
    for val in temp:
        # first word is the gene ID
        genestr = temp[val.find('|')+1:]

        if val in GeneDict.keys():
            returnval.append(GeneDict[val])

    return returnval;

def _MirTarBaseGetPredictionSet():
    df = pd.read_csv( filedialog.askopenfilename(defaultextension='.csv', filetypes=(("comma seperated values file", "*.csv"),("All Files", "*.*"))) )
    temp = df['Target'].tolist()
    global GeneDict
    returnval = []
    for val in temp:
        if type(val) is not str:
            continue

        if val in GeneDict.keys():
            returnval.append(GeneDict[val])

    return returnval

def GenerateExcelOfPredTar(mirna):
    TargetScan = _TargetScanGetPredictionSet(mirna)
    mir_db = _MirDBGetPredictionSet(mirna)
    mirtar = _MirTarBaseGetPredictionSet()
    psrnatar = _psRNATargetGetPredictionSet()

    algorithms = {'TargetScan' : TargetScan, 'mirDB' : mir_db, 'MirTarBase' : mirtar, 'psRNATarget' : psrnatar}
    global GeneDict
    genes = GeneDict.values()

    lists = [{}, {}, {}, {}]
    
    for gene in genes:
        numalgs = 0
        agg = 0
        for key,value in algorithms.items():
            if gene in value:
                numalgs += 1
                agg += value.index(gene)
        if gene == 'TSPAN12|23554':
            print("numalgs: %i agg: %i", numalgs, agg)
        
        if not numalgs == 0:
           lists[numalgs-1][gene] = agg

    wb = openpyxl.Workbook()
    for idx,genedict in enumerate(lists):
        genedict = dict(sorted(genedict.items(), key=lambda item: item[1]))
        ws = wb.create_sheet(f"{idx+1} prediction algorithms")
        for key in genedict.keys():
            ws.append([key])

    wb.save(filedialog.asksaveasfilename(defaultextension='.xlsx', filetypes=(("excel file", "*.xlsx"),("All Files", "*.*") )))

    
# dictionary to match webstites to parser algorithms
'''Websites = {'MirDB' : WebOperationGroup(_MirDBSubmitRequest, _MirDBParseResponse) , 'PicTar' : WebOperationGroup(_PicTarSubmitRequest, _PicTarParseResponse),
            'TargetScan' : WebOperationGroup(_TargetScanSubmitRequest, _TargetScanParseResponse),
            'miRWalk' : WebOperationGroup(_MirWalkSubmitRequest, _MirWalkParseResponse)}

# TODO: change this system so no nested calls, that was dumb
class WebsiteInterface:
    def __init__(self, ID):
        self.ID = ID

    def GetPredictionForGene(self, gene):
        # eventually this may need to take metadata into account
        temp = Websites[self.ID]
        return temp.ParseResponse(temp.SubmitRequest)

    @staticmethod
    def GetPredictionForGene(gene, ID):
        return Websites[ID](gene)'''
    

#Restore('patientsave.pickle')
#mrna_df = ImportTCGAData('TCGA_SKCM_mRNA_9.29.23.csv')
#StoreEverything(mrna_df)
#Save('newpatientsave.pickle'
