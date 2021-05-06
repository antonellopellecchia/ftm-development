#!/usr/bin/python3

import os
import lecroyparser
import re
from array import array
import ROOT as rt
import numpy as np
from scipy import fftpack
from scipy import signal

from physlibs.root import root_style
from physlibs.root import functions

class ScopeData:
    def __init__(self, x, y):
        self.x, self.y = x, y

class ScopeSignal:
    def __init__(self, x=None, y=None, name=None, scopeImpedence=50):
        if y is None: # reading from file
            self.scopeFile = x
            self.scopeData = None
        else: # reading from x and y data
            self.scopeFile = name
            self.scopeData = ScopeData(x, y)
        self.scopeImpedence = scopeImpedence
        self.signalGraph = None
        self.noiseHisto = None

    def __sub__(self, signal):
        return ScopeSignal(self.scopeData.x, self.scopeData.y-self.scopeData.y, str(self), self.scopeImpedence)

    def __add__(self, signal):
        return ScopeSignal(self.scopeData.x, self.scopeData.y+signal.scopeData.y, str(self), self.scopeImpedence)

    def __truediv__(self, constant):
        return ScopeSignal(self.scopeData.x, self.scopeData.y/constant, str(self), self.scopeImpedence)

    def __str__(self):
        return self.GetName()
    
    def __del__(self):
        try:
            if self.signalGraph is not None: self.signalGraph.Delete()
            if self.scopeData is not None: del self.scopeData
        except AttributeError: pass

    def ReadSignal(self, negative=False):
        if self.scopeData is not None: return self.scopeData 
        self.negative = negative
        self.scopeData = lecroyparser.ScopeData(self.scopeFile)
        self.scopeData.x *= 1e9
        self.scopeData.y *= 1e3
        if negative: self.scopeData.y *= -1

    def GetSamplingTime(self):
        try: return self.samplingTime
        except AttributeError:
            self.ReadSignal()
            deltax = [ x-self.scopeData.x[i-1] for i,x in enumerate(self.scopeData.x[1:])]
            self.samplingTime = sum(deltax)/len(deltax) # sampling time in ns
        return self.samplingTime
    
    def CreateGraph(self):
        self.signalGraph = rt.TGraph(len(self.scopeData.x), array('f', self.scopeData.x), array('f', self.scopeData.y))
        self.signalGraph.SetTitle(';Time (ns);Amplitude (mV)')
        return self.signalGraph
    
    def GetName(self):
        return self.scopeFile.split('/')[-1].replace('.trc', '')
    
    def GetGraph(self, fcut=None, tstart=None, tend=None):
        if fcut is not None: return self.GetFilteredGraph(fcut)
        if tstart is None and tend is None:
            if not self.scopeData: self.ReadSignal()
            if not self.signalGraph: self.CreateGraph()
            return self.signalGraph
        x = self.scopeData.x[self.scopeData.x>tstart]
        y = self.scopeData.y[self.scopeData.x>tstart]
        signalGraph = rt.TGraph(len(x), array('f', x), array('f', y))
        signalGraph.SetTitle(';Time (ns);Amplitude (mV)')
        return signalGraph

    def GetNoiseHisto(self, tmax):
        if self.noiseHisto: return self.noiseHisto
        imax = 0
        for t in self.scopeData.x:
            if t>tmax: break
            imax += 1
        noiseMin = min(self.scopeData.y[:imax])
        noiseMax = max(self.scopeData.y[:imax])
        self.noiseHisto = rt.TH1F('noiseHisto_'+self.GetName(), ';Amplitude (mV);Counts', 10, noiseMin-0.5*abs(noiseMin), noiseMax+0.5*abs(noiseMax))
        for i,v in enumerate(self.scopeData.y):
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
        if findBaseline: baseline = self.GetBaseline()
        if fcut>0: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()
        if tLow is None: tLow = self.scopeData.x[0]
        if tHigh is None: tHigh = self.scopeData.x[-1]
        iLow = 0
        iHigh = signalGraph.GetN()-1
        while self.scopeData.x[iLow]<tLow: iLow += 1
        while self.scopeData.x[iHigh]>=tHigh: iHigh -= 1
        if iHigh<=iLow: raise ValueError('tLow should be lower than tHigh')
        integral = 0
        for i in range(iLow, iHigh):
            tdelta = signalGraph.GetPointX(i+1)-signalGraph.GetPointX(i)
            ysum = signalGraph.GetPointY(i)-baseline+signalGraph.GetPointY(i+1)-baseline
            integral += 0.5*ysum*tdelta
        return integral/self.scopeImpedence
    
    def GetChargeRoot(self):
        return self.GetGraph().Integral()/self.scopeImpedence

    def GetChargeBetweenRoot(self, tLow, tHigh=None, fcut=0.):
        if fcut>0: signalGraph = self.GetFilteredGraph(fcut)
        else: signalGraph = self.GetGraph()
        if tHigh is None: tHigh = self.scopeData.x[-1]
        iLow = 0
        iHigh = signalGraph.GetN()-1
        while self.scopeData.x[iLow]<tLow: iLow += 1
        while self.scopeData.x[iHigh]>=tHigh: iHigh -= 1
        #print(iLow, iHigh, self.GetGraph().Integral(iLow, iHigh)/self.scopeImpedence)
        return self.GetGraph().Integral(iLow, iHigh)/self.scopeImpedence

    def GetArrivalTimeByThreshold(self, threshold):
        self.ReadSignal()
        #signalTuples = np.array((self.scopeData.x, self.scopeData.y)).T
        condition = lambda v: v>=threshold
        for (time,voltage) in zip(self.scopeData.x, self.scopeData.y):
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
        #tbot = self.scopeData.x[0]
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
        #try: baseline = self.GetBaseline(self.scopeData.x[-1]-10)
        try: baseline = self.GetBaseline()
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
        x = self.scopeData.x[self.scopeData.x>t0]
        y = self.scopeData.y[self.scopeData.x>t0]

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

    def GetSigmaNoise(self, tmax):
        return self.GetNoiseHisto(tmax).GetRMS()

    def GetBaseline(self, tmax=-5):
        return self.GetNoiseHisto(tmax).GetMean()

    def GetAmplitudeMax(self):
        return max(self.scopeData.y)

    def GetAmplitudeMin(self):
        return min(self.scopeData.y)

    def GetTimeMax(self):
        amplMax = self.GetAmplitudeMax()
        return self.scopeData.x[self.scopeData.y==amplMax][0]

    def GetTimeMin(self):
        amplMin = self.GetAmplitudeMin()
        return self.scopeData.x[self.scopeData.y==amplMin][0]
    
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
        t,y = self.scopeData.x[self.scopeData.x>=tmin], self.scopeData.y[self.scopeData.x>=tmin]
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
        xf, yf = self.GetFFT(self.scopeData.x[0], self.scopeData.x[-1])
        yf[np.abs(xf)>fcut] = 0
        yfiltered = np.real(fftpack.ifft(yf))
        filteredSignal = ScopeSignal(self.scopeData.x, yfiltered, self.GetName(), scopeImpedence=self.scopeImpedence)
        return filteredSignal

    def GetFilteredGraph(self, fcut):
        try:
            if self.fcut==fcut: return self.filteredGraph
        except AttributeError: pass
        self.fcut = fcut
        return self.GetFilteredSignal(fcut).GetGraph()
        #xf, yf = self.GetFFT(self.scopeData.x[0], self.scopeData.x[-1])
        ##xf = fftpack.fftfreq(y.size, d=timestep)
        ##yf = fftpack.fft(self.scopeData.y)
        #yf[np.abs(xf)>fcut] = 0
        ##print(xf)
        ##print(yfft)
        #yfiltered = fftpack.ifft(yf)
        #yfilteredFloat = np.array([ float(y) for y in yfiltered ])
        ##print(yfilteredFloat)
        #self.filteredGraph = rt.TGraph(len(self.scopeData.x), self.scopeData.x, yfilteredFloat)
        #self.filteredGraph.SetTitle(';Time (ns);Amplitude (mV)')
        #return self.filteredGraph
        
    def GetButterworthFilteredSignal(self, fcut, order=4):
        #try: if self.fcut==fcut: return self.filteredGraph
        #except AttributeError: pass
        #self.fcut = fcut
        butA, butB = signal.butter(order, fcut, 'low', analog=True)
        yfiltered = signal.filtfilt(butA, butB, self.scopeData.y)
        #yfilteredFloat = np.array([ float(y) for y in yfiltered ])
        #print(yfilteredFloat)
        self.filteredGraph = rt.TGraph(len(self.scopeData.x), self.scopeData.x, yfiltered)
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

class ScopeSignalSet:
    def __init__(self, signalDirectory):
        self.signalDirectory = signalDirectory
    
    def GetSignals(self, pattern, labels):
        scopeFiles = sorted(os.listdir(self.signalDirectory))
        scopeSignals = list()
        skipped = list()
        for scopeFile in scopeFiles:
            scopeFilePath = self.signalDirectory+'/'+scopeFile
            attributes = dict()
            m = re.match(pattern, scopeFile)
            if not m:
                skipped.append(scopeFile)
                continue
            for i,attr in enumerate(m.groups()):
                attributes[labels[i]] = attr
            scopeSignal = ScopeSignal(scopeFilePath)
            scopeSignal.SetAttributes(attributes)
            scopeSignals.append(scopeSignal)
        self.scopeSignals = scopeSignals
        self.skipped = skipped
        return scopeSignals, skipped