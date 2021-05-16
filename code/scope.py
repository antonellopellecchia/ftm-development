#!/usr/bin/python3

import os
import lecroyparser
import re
from array import array
import ROOT as rt
import numpy as np
import pandas as pd
from scipy import fftpack
from scipy import signal
from tqdm import tqdm

#from physlibs.root import root_style_ftm as root_style
from physlibs.root import functions

class ScopeData:
    def __init__(self, x, y):
        self.x, self.y = x, y
    
    def __repr__(self):
        return str(self.x)+str(self.y)

class ScopeSignal:
    def __init__(self, x, y, name, negative=False, scopeImpedence=50):
        '''if y is None: # reading from file
            self.scopeFile = x
            self.scopeData = None'''
        self.name = name
        if x is None and y is None: self.scopeData = None
        else: self.scopeData = ScopeData(x, y*(-1)**negative)
        self.scopeImpedence = scopeImpedence
        #self._graph = None
        self.noiseHisto = None
        self.negative = negative
    
    def FromFile(filePath, name, negative=False):
        # return an empty signal with file path inside, to read it later
        signal = ScopeSignal(None, None, name, negative)
        signal.scopeFile = filePath
        return signal

    def __sub__(self, signal):
        return ScopeSignal(self.x, self.y-self.y, str(self), self.scopeImpedence)

    def __add__(self, signal):
        return ScopeSignal(self.x, self.y+signal.scopeData.y, str(self), self.scopeImpedence)

    def __truediv__(self, constant):
        return ScopeSignal(self.x, self.y/constant, str(self), self.scopeImpedence)

    def __str__(self):
        return self.GetName()
    
    def __del__(self):
        try:
            if self._graph is not None: del self._graph
            if self.scopeData is not None: del self.scopeData
        except AttributeError: pass

    @property
    def x(self):
        if self.scopeData is None: self.ReadSignal()
        return self.scopeData.x

    @x.setter
    def x(self, x):
        self.scopeData.x = x

    @property
    def y(self):
        if self.scopeData is None: self.ReadSignal()
        return self.scopeData.y

    @y.setter
    def y(self, y):
        self.scopeData.y = y

    def ReadSignal(self):
        if self.scopeData is not None: return self.scopeData
        self.scopeData = lecroyparser.ScopeData(self.scopeFile, sparse = 10000)
        self.x, self.y = np.array(self.x, dtype=float), np.array(self.y, dtype=float)
        self.x *= 1e9
        self.y *= 1e3
        if self.negative: self.y *= -1

    def GetSamplingTime(self):
        try: return self.samplingTime
        except AttributeError:
            self.ReadSignal()
            deltax = [ x-self.x[i-1] for i,x in enumerate(self.x[1:])]
            self.samplingTime = sum(deltax)/len(deltax) # sampling time in ns
        return self.samplingTime
    
    '''def CreateGraph(self):
        self._graph = rt.TGraph(len(self.x), array('f', self.x), array('f', self.y))
        self._graph.SetTitle(';Time (ns);Amplitude (mV)')
        return self._graph'''
    
    def GetName(self):
        return self.scopeFile.split('/')[-1].replace('.trc', '')
    
    @property
    def graph(self, fcut=None, tstart=None, tend=None):
        try: return self._graph
        except AttributeError: pass
        self._graph = rt.TGraph(len(self.x), self.x, self.y)
        self._graph.SetTitle(';Time (ns);Amplitude (mV)')
        return self._graph
        if fcut is not None: return self.GetFilteredGraph(fcut)
        if tstart is None and tend is None:
            if not self.scopeData: self.ReadSignal()
            if not self._graph: self.CreateGraph()
            return self._graph
        x = self.x[self.x>tstart]
        y = self.y[self.x>tstart]

    def GetNoiseList(self, tmax=-5): # return list with points before tmax
        return self.y[self.x<tmax]

    def GetNoiseHisto(self, tmax):
        if self.noiseHisto: return self.noiseHisto
        imax = 0
        for t in self.x:
            if t>tmax: break
            imax += 1
        noiseMin = min(self.y[:imax])
        noiseMax = max(self.y[:imax])
        self.noiseHisto = rt.TH1F('noiseHisto_'+self.GetName(), ';Amplitude (mV);Counts', 10, noiseMin-0.5*abs(noiseMin), noiseMax+0.5*abs(noiseMax))
        for i,v in enumerate(self.y):
            if i>imax: break
            self.noiseHisto.Fill(v)
        return self.noiseHisto

    def SetAttributes(self, attributes):
        self.attributes = attributes

    def GetAttributes(self):
        return self.attributes

    def GetAttributesString(self):
        return '_'.join([ ''.join(it) for it in self.GetAttributes().items() ])
    
    def GetCharge(self):
        return self.GetChargeBetween()

    def GetChargeBetween(self, tLow=None, tHigh=None, fcut=0., baseline=0, findBaseline=True):
        if findBaseline: baseline = self.baseline
        if fcut>0: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.graph
        if tLow is None: tLow = self.x[0]
        if tHigh is None: tHigh = self.x[-1]
        iLow = 0
        iHigh = signalGraph.GetN()-1
        while self.x[iLow]<tLow: iLow += 1
        while self.x[iHigh]>=tHigh: iHigh -= 1
        if iHigh<=iLow: raise ValueError('tLow should be lower than tHigh')
        integral = 0
        for i in range(iLow, iHigh):
            '''tdelta = signalGraph.GetPointX(i+1)-signalGraph.GetPointX(i)
            ysum = signalGraph.GetPointY(i)-baseline+signalGraph.GetPointY(i+1)-baseline'''
            tdelta = self.x[i+1]-self.x[i]
            ysum = self.y[i]-baseline+self.y[i+1]-baseline
            integral += 0.5*ysum*tdelta
        return integral/self.scopeImpedence
    
    def GetChargeRoot(self):
        return self.GetGraph().Integral()/self.scopeImpedence

    def GetChargeBetweenRoot(self, tLow, tHigh=None, fcut=0.):
        if fcut>0: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()
        if tHigh is None: tHigh = self.x[-1]
        iLow = 0
        iHigh = signalGraph.GetN()-1
        while self.x[iLow]<tLow: iLow += 1
        while self.x[iHigh]>=tHigh: iHigh -= 1
        #print(iLow, iHigh, self.GetGraph().Integral(iLow, iHigh)/self.scopeImpedence)
        return self.GetGraph().Integral(iLow, iHigh)/self.scopeImpedence

    def GetArrivalTimeByThreshold(self, threshold):
        self.ReadSignal()
        #signalTuples = np.array((self.x, self.y)).T
        condition = lambda v: v>=threshold
        for (time,voltage) in zip(self.x, self.y):
            if condition(voltage): return time
        raise ValueError('Signal does not cross threshold')

    def GetArrivalTimeByLinearInterpolation(self, threshold, edge, fcut=0.):
        if fcut>0.: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()

        time0 = self.GetArrivalTimeByLinearFit(edge, fcut)
        tstep = 5e-3 # 1 ps precision
        
        v = signalGraph.Eval(time0)
        if edge=='positive': condition = lambda v: v<=threshold
        elif edge=='negative': condition = lambda v: v>=threshold
        while condition(v):
            time0 += tstep
            v = signalGraph.Eval(time0)
        return time0
    
    def GetRisingEdgeFit(self, edge, fcut=0):
        # get approximate 10-90% rise times
        if edge=='positive': ampl = self.GetAmplitudeMax()
        else: ampl = self.GetAmplitudeMin()
        #time30 = self.GetArrivalTimeByThreshold(amplitudeMax*0.3, edge)
        time90 = self.GetArrivalTimeByThreshold(ampl*0.9, edge)
        # fit signal rising edge linearly
        self.pol1 = rt.TF1('pol', 'pol1', time90-1., time90)
        if fcut>0: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()
        signalGraph.Fit(self.pol1, 'R0')
        return self.pol1

    def GetArrivalTimeByLinearFit(self, edge, fcut=0.):
        # calculate intersection of rising edge with 0
        self.pol1 = self.GetRisingEdgeFit(edge, fcut)
        return -self.pol1.GetParameter(0)/self.pol1.GetParameter(1)

    def GetArrivalTimeBySigmoidFitOld(self, edge, fcut=0.): # too complicated, to be simplified
        if edge=='positive': ampl = self.GetAmplitudeMax()
        else: ampl = self.GetAmplitudeMin()
        ttop = self.GetArrivalTimeByThreshold(ampl*0.9, edge)
        tmid = self.GetArrivalTimeByThreshold(ampl*0.5, edge)
        #tbot = self.x[0]
        tbot = tmid-5.
        self.pol1 = self.GetRisingEdgeFit(edge, fcut)
        self.sigmoid = functions.GetSigmoid(tbot, ttop, ampl, self.pol1.GetParameter(1), tmid, 0)
        self.sigmoid.SetLineColor(rt.kRed+2)
        if fcut>0.: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()
        fitResult = signalGraph.Fit(self.sigmoid, 'RS')
        status = int(fitResult)
        if status>0: raise ValueError('Fit failed')
        return self.sigmoid.GetParameter(2)

    def GetArrivalTimeBySigmoidFit(self, findStart=False, tstart=-2, fraction=0.2, fcut=None, status=False, returnChi2=False):
        #try: baseline = self.GetBaseline(self.x[-1]-10)
        try: baseline = self.baseline
        except ValueError as e:
            if status: return 0, -1
            else: raise e
        peak = self.GetAmplitudeMax()
        tend = self.GetArrivalTimeByThreshold(baseline+(peak-baseline)*0.99)
        if findStart: tstart = self.GetArrivalTimeByThreshold(baseline+(peak-baseline)*0.20)-10
        self.sigmoid = functions.GetSigmoid(tstart, tend, peak, (peak-baseline)/(tend-tstart)*50, 0.5*(tstart+tend), baseline)
        self.sigmoid.SetLineColor(rt.kRed+2)
        if fcut is None: signalGraph = self.GetGraph()
        else: signalGraph = self.GetFilteredGraph(fcut)
        fitResult = signalGraph.Fit(self.sigmoid, 'RS')
        if not status and int(fitResult)>0: raise ValueError('Fit failed')
        #return self.sigmoid.GetParameter(2)
        if int(fitResult)==0 and self.sigmoid.GetNDF()>0:
            chi2 = self.sigmoid.GetChisquare()/self.sigmoid.GetNDF()
            if chi2>1: fitResult = 1
        else:
            chi2 = 0
            fitResult = 1
        slope, t0 = self.sigmoid.GetParameter(1), self.sigmoid.GetParameter(2)
        try: sat = t0 - 1/slope*np.log(1/fraction-1) # by default return time at 20% peak
        except ZeroDivisionError:
            return 0, 4, 0
        if returnChi2: return sat, int(fitResult), chi2
        if status: return sat, int(fitResult)
        else: return sat

    def GetSignalCFD(self, att=0.6, delay=1):
        t0 = -20
        x = self.x[self.x>t0]
        y = self.y[self.x>t0]

        indexDelay = int(delay/(x[1]-x[0]))
        yAttenuated = -y*att
        yDelayed = np.roll(y, indexDelay)

        yCFD = yAttenuated+yDelayed
        return ScopeSignal(x, yCFD, name=self.GetName()+'_cfd', scopeImpedence=self.scopeImpedence)

    def GetArrivalTimeCFDByInterpolation(self, att=0.6, delay=1):
        signalCFD = self.GetSignalCFD(att, delay)
        graphCFD = signalCFD.GetGraph()
        t0 = signalCFD.scopeData.x[0]
        t1 = signalCFD.scopeData.x[-1]
        tstep = 1e-3 # ps interpolation precision
        for t in np.arange(t0, t1, tstep):
            if graphCFD.Eval(t0)>0: return t
        #for i,t in enumerate(signalCFD.scopeData.x):
        #    if t<=0 and signalCFD.scopeData.x[i-1]>0: return t
        raise ValueError('Unable to find zero-crossing time with CFD')

    def GetArrivalTimeCFDByLinearFit(self, att=0.7, delay=1):
        signalCFD = self.GetSignalCFD(att, delay)
        #graphCFD = signalCFD.GetGraph()
        t0 = signalCFD.GetTimeMax()
        t1 = signalCFD.GetTimeMin()
        pol1 = rt.TF1('pol', 'pol1', t0, t1)
        signalCFD.GetGraph().Fit(pol1, 'RQ')
        return -pol1.GetParameter(0)/pol1.GetParameter(1)
    
    def GetRisingSlope(self, edge): # return dV/dt
        return self.GetRisingEdgeFit(edge).GetParameter(1)
    
    def GetRiseTime(self, edge):
        fit = self.GetRisingEdgeFit(edge)
        if edge=='negative': return (0.9*self.GetAmplitudeMin()-0.1*self.GetAmplitudeMin())/fit.GetParameter(1)
        else: return (0.9*self.GetAmplitudeMax()-0.1*self.GetAmplitudeMin())/fit.GetParameter(1)

    def GetSigmaNoise(self, tmax=-5):
        #return self.GetNoiseHisto(tmax).GetRMS()
        return self.GetNoiseList(tmax).std()

    @property
    def baseline(self, tmax=-5):
        try: return self._baseline
        except AttributeError: pass
        self._baseline = self.GetNoiseList(tmax).mean()
        print(self._baseline, self.GetNoiseList(tmax))
        return self._baseline
        #return self.GetNoiseHisto(tmax).GetMean()

    def GetAmplitudeMax(self):
        return max(self.y)

    def GetAmplitudeMin(self):
        return min(self.y)

    def GetTimeMax(self):
        amplMax = self.GetAmplitudeMax()
        return self.x[self.y==amplMax][0]

    def GetTimeMin(self):
        amplMin = self.GetAmplitudeMin()
        return self.x[self.y==amplMin][0]
    
    def GetSpectrumSNRatio(self, fcut=None):
        # return power spectrum ratio between signal and non-signal time windows
        npoints = int(50/self.GetSamplingTime()) # must be equal for signal and noise
        xfNoise, yfNoise = self.GetPowerSpectrum(55, npoints=npoints, fcut=fcut)
        xfSignal, yfSignal = self.GetPowerSpectrum(-10, npoints=npoints, fcut=fcut)
        snratio = yfSignal-yfNoise
        return sum(snratio)/len(snratio)
    
    def IsNoise(self, snCut=0.5):
        return self.GetSpectrumSNRatio(fcut=1)<snCut

    ''' utilities for fourier transform '''
    def GetFFT(self, tmin, tmax=None, npoints=None):
        t,y = self.x[self.x>=tmin], self.y[self.x>=tmin]
        if tmax is not None: t,y = t[t<=tmax],y[t<=tmax]
        elif npoints is not None: t,y = t[:npoints],y[:npoints]
        else: raise ValueError('Need to specify either tmax or npoints')
        timestep = self.GetSamplingTime()
        yf = fftpack.fft(y)
        xf = fftpack.fftfreq(y.size, d=timestep)
        return xf, yf

    def GetPowerSpectrum(self, tmin, tmax=None, fcut=None, npoints=None):
        xf, yfft = self.GetFFT(tmin, tmax=tmax, npoints=npoints)
        xf, yf = xf[xf>0], np.log10(np.abs(yfft)**2)[xf>0]
        if fcut is not None: xf, yf = xf[xf<fcut], yf[xf<fcut]
        return xf, yf

    def GetSpectrumGraph(self, tmin, tmax):
        self.fftGraph = rt.TGraph()
        xf, yf = self.GetPowerSpectrum(tmin, tmax)
        for i,(x,y) in enumerate(zip(xf, yf)):
            self.fftGraph.SetPoint(i+1, x, y)
        self.fftGraph.SetTitle(';Frequency (GHz);Amplitude (dB)')
        #self.fftGraph.GetXaxis().SetRangeUser(0., xf[-1])
        return self.fftGraph

    def GetSignalOverNoiseGraph(self, fcut=None):
        self.snGraph = rt.TGraph()
        npoints = int(50/self.GetSamplingTime()) # must be equal for signal and noise
        xfNoise, yfNoise = self.GetPowerSpectrum(55, npoints=npoints, fcut=fcut)
        xfSignal, yfSignal = self.GetPowerSpectrum(-10, npoints=npoints, fcut=fcut)
        for i,(x,y) in enumerate(zip(xfSignal, yfSignal)):
            self.snGraph.SetPoint(i+1, x, y-yfNoise[i])
        self.snGraph.SetTitle(';Frequency (GHz);Amplitude (dB)')
        #self.snGraph.GetXaxis().SetRangeUser(0., xf[-1])
        return self.snGraph

    def GetFilteredSignal(self, fcut):
        xf, yf = self.GetFFT(self.x[0], self.x[-1])
        yf[np.abs(xf)>fcut] = 0
        yfiltered = np.real(fftpack.ifft(yf))
        filteredSignal = ScopeSignal(self.x, yfiltered, self.GetName(), scopeImpedence=self.scopeImpedence)
        return filteredSignal

    def GetFilteredGraph(self, fcut):
        try:
            if self.fcut==fcut: return self.filteredGraph
        except AttributeError: pass
        self.fcut = fcut
        return self.GetFilteredSignal(fcut).GetGraph()
        #xf, yf = self.GetFFT(self.x[0], self.x[-1])
        ##xf = fftpack.fftfreq(y.size, d=timestep)
        ##yf = fftpack.fft(self.y)
        #yf[np.abs(xf)>fcut] = 0
        ##print(xf)
        ##print(yfft)
        #yfiltered = fftpack.ifft(yf)
        #yfilteredFloat = np.array([ float(y) for y in yfiltered ])
        ##print(yfilteredFloat)
        #self.filteredGraph = rt.TGraph(len(self.x), self.x, yfilteredFloat)
        #self.filteredGraph.SetTitle(';Time (ns);Amplitude (mV)')
        #return self.filteredGraph
        
    def GetButterworthFilteredSignal(self, fcut, order=4):
        #try: if self.fcut==fcut: return self.filteredGraph
        #except AttributeError: pass
        #self.fcut = fcut
        butA, butB = signal.butter(order, fcut, 'low', analog=True)
        yfiltered = signal.filtfilt(butA, butB, self.y)
        #yfilteredFloat = np.array([ float(y) for y in yfiltered ])
        #print(yfilteredFloat)
        self.filteredGraph = rt.TGraph(len(self.x), self.x, yfiltered)
        self.filteredGraph.SetTitle(';Time (ns);Amplitude (mV)')
        return self.filteredGraph
    
    def SetCanvas(self, canvas):
        self.canvas = canvas
    
    def Draw(self, opt=''):
        signalGraph = self.GetGraph()
        self.canvas.cd()
        signalGraph.SetLineWidth(1)
        signalGraph.Draw(opt)
        latex = rt.TLatex()
        latex.SetTextAlign(32)
        latex.SetTextSize(.035)
        latex.DrawLatexNDC(.95, .97, 'Signal charge %1.2f pC'%(self.GetCharge()))
    
    def SaveAs(self, fname):
        signalCanvas = rt.TCanvas('signalCanvas'+self.name, '', 800, 600)
        graph = self.GetGraph()
        graph.Draw('AL')
        try: os.makedirs(os.path.dirname(fname))
        except FileExistsError: pass
        signalCanvas.SaveAs(fname)


class ScopeEvent:
    # stores one or more signals at a single trigger
    def __init__(self, charge=None, amplitude=None, overThreshold=None, preamplifier=None):
        if charge is not None: self._detectorCharge = charge
        if overThreshold is not None: self._isOverThreshold = overThreshold
        if amplitude is not None: self._detectorAmplitude = amplitude
        self.preamplifier = preamplifier

    def FromSignals(signals, triggerSignals, detectorSignals, threshold):
        event = ScopeEvent()
        event.signals = signals # list of signals, both trigger and detector
        event.triggerSignals = triggerSignals # index of trigger signals in event list
        event.detectorSignals = detectorSignals  # index of detector signals in event list
        event.threshold = threshold
        return event

    def __repr__(self):
        return self.name

    @property
    def name(self):
        return self.signals[0].name

    @property
    def triggerSignal(self):
        return self.signals[self.triggerSignals[0]]

    @property
    def detectorSignal(self):
        return self.signals[self.detectorSignals[0]]

    @property
    def detectorCharge(self):
        try: return self._detectorCharge
        except AttributeError: pass
        self._detectorCharge = self.detectorSignal.GetCharge()
        return self._detectorCharge

    @property
    def collectedCharge(self):
        if self.preamplifier: return self.detectorAmplitude/self.preamplifier
        else: return self.detectorCharge

    @property
    def detectorAmplitude(self):
        try: return self._detectorAmplitude
        except AttributeError: pass
        self._detectorAmplitude = self.detectorSignal.GetAmplitudeMax()
        return self._detectorAmplitude

    @property
    def isOverThreshold(self):
        try: return self._isOverThreshold
        except AttributeError: pass
        self._isOverThreshold = self.detectorSignal.GetAmplitudeMax()>self.threshold
        return self._isOverThreshold

    def SaveAs(self, fname):
        # plot trigger and detector signals side to side
        signalCanvas = rt.TCanvas('signalCanvas'+self.name, '', 1200, 600)
        signalCanvas.Divide(2,1)
        signalCanvas.cd(1)
        
        multiGraphs = [rt.TMultiGraph(),rt.TMultiGraph()]
        for i,g in enumerate(multiGraphs):
            signalCanvas.cd(i+1)
            for j in [self.triggerSignals,self.detectorSignals][i]:
                signal = self.signals[j]
                #print(signal.graph, type(signal.graph))
                multiGraphs[i].Add(signal.graph, 'l')
                #print('Processing signal', j)
            multiGraphs[i].Draw('a')

        thresholdLine = rt.TLine(-200, self.threshold, 200, self.threshold)
        thresholdLine.SetLineColor(rt.kRed)
        thresholdLine.SetLineStyle(2)
        thresholdLine.Draw()

        baseLine = rt.TLine(-200, self.signals[-1].baseline, 200, self.signals[-1].baseline)
        baseLine.SetLineColor(rt.kGreen)
        baseLine.SetLineStyle(2)
        baseLine.Draw()

        latex = rt.TLatex()
        latex.SetTextSize(0.03)
        latex.DrawLatexNDC(.18, .87, f'Amplitude {self.detectorAmplitude:1.2f} mV')
        latex.DrawLatexNDC(.18, .82, f'Charge {self.detectorCharge:1.2f} pC')
        latex.DrawLatexNDC(.18, .77, f'Baseline {self.signals[-1].baseline:1.2f} mV')

        '''signalCanvas.cd(1)
        self.signals[0].graph.Draw('AL')
        signalCanvas.cd(2)
        self.signals[1].graph.Draw('AL')'''

        try: os.makedirs(os.path.dirname(fname))
        except FileExistsError: pass
        signalCanvas.SaveAs(fname)


class DataTaking:
    # groups many signals acquired at fixed setup conditions
    def __init__(self, scopeEvents, setup):
        self.scopeEvents = scopeEvents
        #self.signalsTrigger = signalsTrigger
        #self.signalsDetector = signalsDetector
        self.setup = setup
    
    def FromDirectory(path, setup, chtrigger, chdetector, negative):
        signalsTrigger = list()
        signalsDetector = list()
        
        threshold = float(setup['THRESHOLD'])

        files = sorted(os.listdir(path))
        for f in files:
            fileNameMatch = re.match(r'C(\d)--Trace--(\d+).trc', f)
            if not fileNameMatch: continue
            channel, traceNumber = int(fileNameMatch.group(1)), fileNameMatch.group(2)
            name = f'C{channel}-{traceNumber}'

            otherChannels = [chtrigger,chdetector]
            if not channel in otherChannels: continue # exclude files that are not signals

            # keep only files that appear coupled with other channels
            # to avoid remainings of previous acquisition
            # should not happen if datataking was done correctly
            otherChannels.remove(channel)
            exists = True
            for otherChannel in otherChannels:
                if not os.path.exists('%s/C%d--Trace--%s.trc'%(path, otherChannel, traceNumber)): exists = False

            # separate files between trigger and signal
            if not exists: continue
            elif channel==chtrigger: signalsTrigger.append(ScopeSignal.FromFile(f'{path}/{f}', name))
            elif channel==chdetector: signalsDetector.append(ScopeSignal.FromFile(f'{path}/{f}', name, negative))
            else: continue

        scopeEvents = list()
        for signalTrigger,signalDetector in zip(signalsTrigger, signalsDetector):
            scopeEvents.append(ScopeEvent.FromSignals([signalTrigger,signalDetector], [0], [1], threshold))
        dataTaking = DataTaking(scopeEvents, setup)
        return dataTaking

    def FromNtuples(path, setup):
        #rootFile = rt.TFile(path)
        treeDf = rt.RDataFrame('DataTakingTree', path)
        eventTree = treeDf.AsNumpy()
        scopeEvents = list()

        if setup['READOUT']=='ortec': preamplifier = 4 # Ortec sensitivity in mV/fC
        else: preamplifier = None
        
        for charge,amplitude,overThreshold in zip(eventTree['detectorCharge'],eventTree['detectorAmplitude'],eventTree['isOverThreshold']):
            scopeEvents.append(ScopeEvent(charge,amplitude,overThreshold,preamplifier))
        return DataTaking(scopeEvents, setup)
    
    def __iter__(self):
        return self.scopeEvents.__iter__()

    def __str__(self):
        return str(self.setup)

    def __len__(self):
        return len(self.scopeEvents)

    @property
    def name(self):
        return '_'.join([ str(s) for s in self.setup.values()])

    def ToNtuples(self, path, events=-1):
        rootFile = rt.TFile(path, 'RECREATE')
        tree = rt.TTree('DataTakingTree', '')
        
        detectorCharge = np.zeros(1, dtype=float)
        tree.Branch('detectorCharge', detectorCharge, 'detectorCharge/D')
        isOverThreshold = np.zeros(1, dtype=float)
        tree.Branch('isOverThreshold', isOverThreshold, 'isOverThreshold/D')
        detectorAmplitude = np.zeros(1, dtype=float)
        tree.Branch('detectorAmplitude', detectorAmplitude, 'detectorAmplitude/D')

        for event in tqdm(self.scopeEvents[:events]):
            detectorCharge[0] = event.detectorCharge
            detectorAmplitude[0] = event.detectorAmplitude
            isOverThreshold[0] = event.isOverThreshold
            tree.Fill()
        
        rootFile.Write()
        rootFile.Close()

    def GetEfficiency(self):
        try: return self._efficiency
        except AttributeError: pass
        nevents = len(self)
        eventsOverThreshold = [ event for event in self if event.isOverThreshold ]
        nsignals = len(eventsOverThreshold)
        self._efficiency = nsignals/nevents
        return self._efficiency

    def GetAverageCharge(self):
        return np.array([event.collectedCharge for event in self.scopeEvents]).mean()

    def DrawWaveforms(self, path, nevents=None):
        if nevents is None: nevents=len(self.scopeEvents)
        for scopeEvent in self.scopeEvents[:nevents]:
            #print(scopeEvent)
            scopeEvent.SaveAs(f'{path}/{scopeEvent.name}.png')


class Measurement:
    # groups multiple datatakings at different setup conditions (e.g. ampl. field)
    def __init__(self, dataTakings, setup):
        self.dataTakings = dataTakings
        self.setup = setup

    def FromDirectory(path, parameters, chtrigger, chdetector, negative):
        # initialize from datatakings in directory
        subfolders = sorted(os.listdir(path))
        subfolders = [ d for d in subfolders if os.path.isdir(f'{path}/{d}') ]

        setupFile = f'{path}/setup.txt'
        measurementSetupDict = pd.read_csv(setupFile, sep='\t').to_dict()
        measurementSetup = {key:arg[0] for (key,arg) in measurementSetupDict.items()}

        dataTakings = list()
        for d in subfolders:
            dataTakingSetup = {key:arg for (key,arg) in zip(parameters,d.split('_'))}
            dataTakingSetup.update(measurementSetup)
            dataTakings.append(DataTaking.FromDirectory(f'{path}/{d}', dataTakingSetup, chtrigger, chdetector, negative))
        return Measurement(dataTakings, measurementSetup)

    def FromNtuples(path, parameters):
        # initialize from datatakings in directory
        rootFiles = sorted(os.listdir(path))
        rootFiles = [f for f in rootFiles if f[-5:]=='.root']

        setupFile = f'{path}/setup.txt'
        measurementSetupDict = pd.read_csv(setupFile, sep='\t').to_dict()
        measurementSetup = {key:arg[0] for (key,arg) in measurementSetupDict.items()}

        dataTakings = list()
        for rootFile in rootFiles:
            dataTakingSetup = {key:arg for (key,arg) in zip(parameters,rootFile[:-5].split('_'))}
            dataTakingSetup.update(measurementSetup)
            dataTakings.append(DataTaking.FromNtuples(f'{path}/{rootFile}', dataTakingSetup))
        return Measurement(dataTakings, measurementSetup)

    def ToNtuples(self, path, events=-1):
        for i,dataTaking in enumerate(self):
            print('Point', i, dataTaking)
            dataTaking.ToNtuples(f'{path}/{dataTaking.name}.root', events)
    
    def __iter__(self):
        return self.dataTakings.__iter__()
    
    def GetEfficiencyPlot(self):
        efficiencyPlot = rt.TGraphErrors()
        opticalDensity = float(self.setup['OD'])
        for i,dataTaking in tqdm(enumerate(self)):
            laserPower = float(dataTaking.setup['power'])*51/100*10**(-opticalDensity)
            averageCharge = dataTaking.GetAverageCharge()
            #efficiencyPlot.SetPoint(i, laserPower, dataTaking.GetEfficiency())
            efficiencyPlot.SetPoint(i, averageCharge, dataTaking.GetEfficiency())
        efficiencyPlot.SetTitle(';Collected charge (fC);Efficiency')
        return efficiencyPlot

    def GetEfficiencyPlotPrimaries(self, calibrationCurve, calibration='gain'):
        efficiencyPlot = rt.TGraphErrors()
        opticalDensity = float(self.setup['OD'])
        for i,dataTaking in tqdm(enumerate(self)):
            if calibration=='gain': # retrieve primary charge from gain curve
                gain = calibrationCurve.GetGain(self.setup['AMPLIFICATION'])
                averageCharge = dataTaking.GetAverageCharge() # in fC
                primaryCharge = averageCharge/gain
                primaryElectrons = primaryCharge*1e-15/1.6e-19
            elif calibration=='current': # retrieve primary charge from primary current scan
                laserPower = float(dataTaking.setup['power'])*51/100*10**(-opticalDensity)
                primaryCharge = calibrationCurve.GetPrimaryCurrent(laserPower)/100
                primaryElectrons = primaryCharge*1e-9/1.6e-19
            efficiencyPlot.SetPoint(i, primaryElectrons, dataTaking.GetEfficiency())
        efficiencyPlot.SetTitle(';Primary electrons;Efficiency')

        return efficiencyPlot