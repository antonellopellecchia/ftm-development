#!/usr/bin/python3

import os, re
from array import array
import ROOT as rt
import numpy as np
import pandas as pd

from physlibs.root import root_style_ftm

class GainCurve:
    def __init__(self, path):
        self.rootFile = rt.TFile(path)
        self.gainGraph = self.rootFile.Get('Graph')
        self.gainFit = self.gainGraph.GetFunction('f')
    
    def ls(self):
        self.rootFile.ls()

    def GetGain(self, amplVoltage):
        return self.gainFit.Eval(amplVoltage)
    
class PrimaryCurrentCurve:
    def __init__(self, path, gainCurve):
        self.rootFile = rt.TFile(path)
        self.currentGraph = self.rootFile.Get('Graph')
        self.currentFit = self.currentGraph.GetFunction('powerScanFit')
        self.gain = gainCurve.GetGain(300) # power scan is done at 300 V
    
    def ls(self):
        self.rootFile.ls()

    def GetPrimaryCurrent(self, sourcePower):
        return self.currentFit.Eval(sourcePower)/self.gain

class CurrentMeasurement:
    def __init__(self, currents, parameters=None, tstep=1):
        self.currents = currents
        self.parameters = parameters
        self.tstep = tstep

    def FromFile(currentPath, parameters=None, sep='\t', tstep=1, multiplier=1e9):
        with open(currentPath) as f:
            currentsStrList = f.read().split(sep)
            currentsFloatList = [ float(s) for s in currentsStrList ]
        currents = np.array(currentsFloatList)*multiplier
        return CurrentMeasurement(currents, parameters, tstep)

    def __str__(self): return f'{self.mean:.2e} +/- {self.error:.2e}'

    def __mul__(self, constant):
        return CurrentMeasurement(self.currents*constant, self.parameters, self.tstep)
    
    def __rmul__(self, constant):
        return CurrentMeasurement(self.currents*constant, self.parameters, self.tstep)

    def __abs__(self):
        return CurrentMeasurement(abs(self.currents), self.parameters, self.tstep)

    def __add__(self, other):
        return CurrentMeasurement(self.currents+other.currents, self.parameters, self.tstep)

    @property
    def title(self):
        return '_'.join([ str(s) for s in self.parameters.values()])

    @property
    def npoints(self):
        return self.currents.size
    
    @property
    def mean(self): return self.currents.mean()

    @property
    def rms(self): return self.currents.std()

    @property
    def error(self): return self.rms/self.npoints**0.5

    @property
    def plot(self):
        try: return self._plot
        except AttributeError:
            times = np.linspace(0, self.npoints*self.tstep, self.npoints)
            self._plot = rt.TGraph()
            for i,(t,c) in enumerate(zip(times,self.currents)):
                self._plot.SetPoint(i, t, c)
            self._plot.SetTitle(';Time (s);Current (nA)')
            return self._plot

class CurrentMeasurementOnOff:
    def __init__(self, measurementOn, measurementOff):
        self.measurementOn = measurementOn
        self.measurementOff = measurementOff

    def __str__(self):
        return f'{self.measurementOn.parameters} {self.mean:.2e} +/- {self.error:.2e}'

    def __repr__(self):
        return str(self)

    def __abs__(self):
        return CurrentMeasurementOnOff(abs(self.measurementOn), abs(self.measurementOff))

    def __add__(self, other):
        return CurrentMeasurementOnOff(self.measurementOn+other.measurementOn, self.measurementOff+other.measurementOff)

    def __radd__(self, other):
        if other==0: return self
        return self+other

    @property
    def title(self):
        return self.measurementOn.title

    @property
    def mean(self):
        return self.measurementOn.mean-self.measurementOff.mean
        
    @property
    def error(self):
        return (self.measurementOn.error**2+self.measurementOff.error**2)**0.5

    @property
    def parameters(self):
        if self.measurementOn.parameters==self.measurementOff.parameters:
            return self.measurementOn.parameters
        else: raise ValueError('Parameters with source on and off should be the identical')
    
    @property
    def plot(self):
        try: return self._plot
        except AttributeError:
            self._plot = rt.TMultiGraph()
            legend = rt.TLegend(.18, .95, .6, .99)
            legend.SetNColumns(2)
            colors = [rt.kRed,rt.kBlue]
            l = ['on', 'off']
            for j,g in enumerate([self.measurementOn.plot, self.measurementOff.plot]):
                g.SetLineColor(colors[j])
                self._plot.Add(g, 'pl')
                legend.AddEntry(g, 'Source '+l[j], 'pl')
            self._plot.SetTitle(';Time (s);Current (nA)')
            return self._plot, legend

class CurrentScan:
    def __init__(self, currents, parnames, xvariable, setup):
        self.currents = currents
        self.parnames = parnames
        self.xvariable = xvariable
        self.setup = setup

    def FromDirectory(directory, pattern, parnames, xvariable, runParameters=None):
        ''' get run parameters from setup file '''
        if runParameters==None: runParameters = f'{directory}/setup.txt'
        setup = dict(pd.read_csv(runParameters, sep='\t'))
        for key,val in setup.items():
            try: setup[key] = float(val)
            except ValueError: setup[key] = val

        ''' read current measurement files '''
        currents = list()
        inputFiles = sorted(os.listdir(directory))

        filesoff, fileson = list(), list()
        for inputFile in inputFiles: # scan files in directory to find current data
            match = re.match(pattern, inputFile)
            if match: # file contains currents data
                if '_off' in inputFile: filesoff.append(inputFile)
                elif '_on' in inputFile: fileson.append(inputFile)

        # create CurrentMeasurementOnOff for each current data file
        for i,(fileoff,fileon) in enumerate(zip(filesoff,fileson)):
            matchOff = re.match(pattern, fileoff)
            matchOn = re.match(pattern, fileon)
            title = '.'.join(fileon.split('.')[:-1])
            if not fileoff.replace('_off', '')==fileon.replace('_on', ''): pass

            parameters = dict()
            for j,parname in enumerate(parnames):
                par = matchOn.group(j+1)
                try: par = float(par)
                except ValueError: pass
                parameters[parname] = par
            measurementOn = CurrentMeasurement.FromFile(f'{directory}/{fileon}', parameters)
            measurementOff = CurrentMeasurement.FromFile(f'{directory}/{fileoff}', parameters)
            currents.append(CurrentMeasurementOnOff(measurementOn, measurementOff))

        return CurrentScan(currents, parnames, xvariable, setup)

    def __iter__(self):
        return iter(self.currents)

    def __str__(self):
        return str(list(self))

    '''def __add__(self, other):
        currentSum = [ c1+c2 for c1,c2 in zip(self, other) ]
        return CurrentScan(currentSum, self.parnames, self.xvariable, self.setup)

    def __radd__(self, other): return self+other'''

    def __abs__(self):
        absCurrents = [ abs(c) for c in self.currents ]
        return CurrentScan(absCurrents, self.parnames, self.xvariable, self.setup)

    def __len__(self): return len(self.currents)

    def __getitem__(self, i): return self.currents[i]

    @property
    def xvalues(self):
        return [ current.parameters[self.xvariable] for current in self ]

    @property
    def plot(self):
        #try: return self._plot
        #except AttributeError:
        plot = rt.TGraphErrors()
        for i,current in enumerate(self):
            plot.SetPoint(i, float(current.parameters[self.xvariable]), current.mean)
            plot.SetPointError(i, 0, current.error)
        return plot

    def Get(self, x):
        return self.currents[self.xvalues.index(x)]

    def GetParameter(self, key): return self.setup[key]

    def SetScalex(self, multiplier):
        for current in self: current.parameters[self.xvariable] *= multiplier

    def SetOffsetx(self, offset):
        for current in self: current.parameters[self.xvariable] += offset

    def SaveCurrentPlots(self, saveFolder):
        for current,x in zip(self,self.xvalues):
            currentTimeCanvas = rt.TCanvas(f'CurrentTimeCanvas{current.title}', '', 1000, 600)
            plot, legend = current.plot
            plot.Draw('a')
            legend.Draw()
            latex = rt.TLatex()
            latex.SetTextSize(.035)
            latex.SetTextAlign(32)
            latex.DrawLatexNDC(.95, .97, f'{x} V')

            currentTimeCanvas.SaveAs(f'{saveFolder}/{current.title}.eps')

class MultiElectrodeScan:
    def __init__(self, electrodeScans, electrodeNames):
        self.scans = electrodeScans
        self.names = electrodeNames

    def EnableTotal(self):
        self.scans.append(self.total)
        self.names.append('total')

    @property
    def total(self):
        # sum of all the scans
        xs = set(self.scans[0].xvalues)
        for scan in self.scans: xs = xs & set(scan.xvalues)
        currents = [
            sum(scan.Get(x) for scan in self.scans)
            for x in xs
        ]
        return CurrentScan(currents, self.scans[0].parnames, self.scans[0].xvariable, self.scans[0].setup)

    @property
    def plot(self):
        try: return self._plot
        except AttributeError:
            self._plot = rt.TMultiGraph()
            self._plot.SetTitle(';Amplification voltage (V);Electrode current (nA)')
            plots = [ scan.plot for scan in self.scans ]
            colors = [ rt.kRed+2, rt.kBlue+2, rt.kGreen+2, rt.kOrange+1 ]
            for plot,color,name in zip(plots,colors,self.names):
                plot.SetMarkerColor(color)
                self._plot.Add(plot, 'p')
                self.legend.AddEntry(plot, name, 'p')
            return self._plot