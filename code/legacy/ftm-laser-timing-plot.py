#!/usr/bin/python3

import os, sys
import numpy as np
import ROOT as rt
import pandas as pd
import argparse
import time, datetime
import re

from ftm_analysis import scope
from physlibs.root import root_style_ftm
from physlibs.root import functions

operations = {
    '<': lambda x,t: x<t,
    '>': lambda x,t: x>t,
    '=': lambda x,t: x==t
}

def applyCut(event, variable, operation, value):
    return operations[operation](event[variable],value)

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input')
    ap.add_argument('--output')
    ap.add_argument('--cut', type=str)
    ap.add_argument('--time', nargs='+', type=float)
    ap.add_argument('--triggercharge', nargs='+', type=float)
    ap.add_argument('--triggerjitter')
    ap.add_argument('--charge', nargs='+', type=float)
    ap.add_argument('--amplitude', nargs='+', type=float)
    ap.add_argument('--risetime', nargs='+', type=float)
    ap.add_argument('--psd', nargs='+', type=float)
    ap.add_argument('--pedestal', action='store_true')
    ap.add_argument('--setup')
    ap.add_argument('--label', type=str)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', help='Ignored', action='store_true')
    options = ap.parse_args(sys.argv[1:])

    outdir = options.output

    if options.setup:
        setupCsv = pd.read_csv(options.setup, header=0, index_col=0).T
        power = setupCsv['power'][0]
        bias = setupCsv['bias'][0]

    for d in [outdir]:#, outdir+'/TimeByAmplitude', outdir+'/TimeByCharge']:
        try: os.makedirs(d)
        except FileExistsError: pass

    treeFile = rt.TFile(options.input, 'READ')
    eventTree = treeFile.Get('EventTree')
    eventTree.Print()
    eventData, eventColumns = eventTree.AsMatrix(return_labels=True)
    eventDataFrame = pd.DataFrame(data=eventData, columns=eventColumns)

    #isEventFilter = eventDataFrame['isEvent']==True
    #eventDataFrame = eventDataFrame[isEventFilter]
    print(eventDataFrame)

    results = dict() # save time jitter, avg charge etc.
    
    if options.risetime:
        ''' plot rise time histogram '''
        riseTimeBot = options.risetime[0]
        riseTimeTop = options.risetime[1]
        riseTimeBins = int(options.risetime[2])
        
        riseTimeHisto = rt.TH1F('RiseTimeHisto', ';Signal rise time (ns);', riseTimeBins, riseTimeBot, riseTimeTop)

        riseTimes = eventDataFrame['riseTime']
        for riseTime in riseTimes:
            riseTimeHisto.Fill(riseTime)

        riseTimeSpectrumCanvas = rt.TCanvas('RiseTimeSpectrumCanvas', '', 800, 600)
        riseTimeHisto.Draw()
        riseTimeSpectrumCanvas.SaveAs('%s/RiseTimeSpectrum.eps'%(outdir))

    if options.charge:
        ''' plot charge spectrum '''
        chargeBot = options.charge[0]
        chargeTop = options.charge[1]
        chargeBins = int(options.charge[2])
        
        psdBot = options.psd[0]
        psdTop = options.psd[1]
        psdBins = int(options.psd[2])

        chargeHisto = rt.TH1F('ChargeHisto', ';Collected charge (pC);', chargeBins, chargeBot, chargeTop)
        chargeOutHisto = rt.TH1F('ChargeOutHisto', ';Collected charge (pC);', chargeBins, chargeBot, chargeTop)
        psdHisto = rt.TH1F('PsdHisto', ';PSD;', psdBins, psdBot, psdTop)
        chargePsdHisto = rt.TH2F('ChargePsdHisto', ';Signal charge (pC);PSD;', chargeBins, chargeBot, chargeTop, psdBins, psdBot, psdTop)

        chargesLong = eventDataFrame['chargeLong']
        chargesOut = eventDataFrame['chargeOut']
        PSDs = eventDataFrame['psd']
        for chargeLong,chargeOut,psd in zip(chargesLong,chargesOut,PSDs):
            chargeHisto.Fill(chargeLong)
            chargeOutHisto.Fill(chargeOut)
            psdHisto.Fill(psd)
            chargePsdHisto.Fill(chargeLong, psd)

        chargeSpectrumCanvas = rt.TCanvas('ChargeSpectrumCanvas', '', 800, 600)
        
        gaus = rt.TF1('g', 'gaus(0)')
        chargeOutHisto.Fit(gaus)
        
        chargeHisto.SetLineColor(rt.kBlue+2)
        chargeHisto.Draw()

        chargeHisto.Fit(gaus)
        averageCharge = gaus.GetParameter(1)
        averageChargeError = gaus.GetParError(1)

        chargeSpectrumCanvas.Update()
        chargeSpectrumCanvas.Draw()
        fitBox = chargeHisto.FindObject('stats')
        fitBox.SetTextColor(rt.kRed+2)
        fitBox.SetX1NDC(0.7)
        fitBox.SetY1NDC(0.62)
        fitBox.SetX2NDC(0.95)
        fitBox.SetY2NDC(0.92)

        '''latex = rt.TLatex()
        latex.SetTextSize(.029)
        latex.SetTextAlign(12)
        latex.DrawLatexNDC(.2, .9, 'Average collected charge %1.2f pC'%(averageCharge))'''

        #chargeOutHisto.SetLineColor(rt.kBlue)
        #chargeOutHisto.Draw('same')
        
        if options.label: root_style_ftm.labelRight(chargeSpectrumCanvas, options.label)
        root_style_ftm.labelFtm(chargeSpectrumCanvas)
        #chargeSpectrumCanvas.SetLogy()
        chargeSpectrumCanvas.SaveAs('%s/ChargeSpectrum.eps'%(outdir))

        psdCanvas = rt.TCanvas('PsdCanvas', '', 800, 600)
        psdHisto.Draw()
        psdCanvas.SaveAs('%s/psd.eps'%(outdir))

        chargePsdCanvas = rt.TCanvas('ChargePsdCanvas', '', 800, 600)
        chargePsdHisto.Draw('colz')
        chargePsdCanvas.SaveAs('%s/ChargePSD.eps'%(outdir))

        results['charge'] = averageCharge
        results['errCharge'] = chargeHisto.GetRMS()/chargeHisto.GetEntries()**0.5
    
    if options.amplitude:
        ''' plot amplitude histogram '''
        amplitudeBot = options.amplitude[0]
        amplitudeTop = options.amplitude[1]
        amplitudeBins = int(options.amplitude[2])
        
        amplitudeHisto = rt.TH1F('AmplitudeHisto', ';Signal amplitude (mV);', amplitudeBins, amplitudeBot, amplitudeTop)
        amplitudes = eventDataFrame['amplitude']
        for amplitude in amplitudes:
            amplitudeHisto.Fill(amplitude)

        amplitudeSpectrumCanvas = rt.TCanvas('AmplitudeSpectrumCanvas', '', 800, 600)
        amplitudeHisto.Draw()
        amplitudeSpectrumCanvas.SaveAs('%s/AmplitudeSpectrum.eps'%(outdir))

        results['amplitude'] = amplitudeHisto.GetMean()
        results['errAmplitude'] = amplitudeHisto.GetRMS()/amplitudeHisto.GetEntries()**0.5

    if options.triggercharge:
        ''' determine APD time resolution from previous fit '''
        triggerChargeBot = options.triggercharge[0]
        triggerChargeTop = options.triggercharge[1]
        triggerChargeBins = int(options.triggercharge[2])
        
        triggerChargeHisto = rt.TH1F('TriggerChargeHisto', ';Collected charge (pC);', triggerChargeBins, triggerChargeBot, triggerChargeTop)

        chargesTrigger = eventDataFrame['chargeTrigger']
        for charge in chargesTrigger:
            triggerChargeHisto.Fill(charge)

        triggerChargeSpectrumCanvas = rt.TCanvas('TriggerChargeSpectrumCanvas', '', 800, 600)
        triggerChargeHisto.Draw()
        fit = rt.TF1('g', 'gaus(0)')
        triggerChargeHisto.Fit(fit)
        triggerChargeSpectrumCanvas.Update()
        triggerChargeSpectrumCanvas.Draw()
        fitBox = triggerChargeHisto.FindObject('stats')
        fitBox.SetTextColor(rt.kRed+2)
        fitBox.SetX1NDC(0.25)
        fitBox.SetY1NDC(0.45)
        fitBox.SetX2NDC(0.55)
        fitBox.SetY2NDC(0.85)
        triggerChargeSpectrumCanvas.SaveAs('%s/TriggerChargeSpectrum.eps'%(outdir))

        apdCharge = fit.GetParameter(1)
        print(apdCharge)
        errApdCharge = fit.GetParError(1)
        
        jitterFile = rt.TFile(options.triggerjitter)
        jitterFile.ls()
        jitterFunction = jitterFile.Get('h')
        apdJitter = jitterFunction.Eval(apdCharge)*1e-3
    
    if options.time:
        ''' plot arrival time spectrum '''
        timeBot = options.time[0]
        timeTop = options.time[1]
        timeBins = int(options.time[2])

        timeHisto = rt.TH1F('TimeHisto', ';Signal time (ns);', timeBins, timeBot, timeTop)

        timesTrigger = eventDataFrame['timeTrigger']
        timesSignal = eventDataFrame['timeSignal']
        for timeTrigger,timeSignal in zip(timesTrigger,timesSignal):
            timeHisto.Fill(timeSignal-timeTrigger)

        timeCanvas = rt.TCanvas('TimeCanvas', '', 800, 600)

        ampl1, mean1, sigma1 = timeHisto.GetEntries()/10, timeHisto.GetMean(), timeHisto.GetRMS()/3
        ampl2, mean2, sigma2 = timeHisto.GetEntries()/100, timeHisto.GetMean(), timeHisto.GetRMS()*1.5
        #ampl3, mean3, sigma3 = timeHisto.GetEntries()/100, timeHisto.GetMean(), timeHisto.GetRMS()*2
        fit = functions.Get2Gauss(timeBot, timeTop, ampl1, mean1, sigma1, ampl2, mean2, sigma2) # 2gauss
        #fit = functions.Get3Gauss(timeBot, timeTop, ampl1, mean1, sigma1, ampl2, mean2, sigma2, ampl3, mean3, sigma3)
        #fit = functions.Get3Gauss(timeBot, timeTop, 300, -1.2, 70e-3, 300, -0.8, 70e-3, 150, 0.15, 100e-3)
        #fit.SetLineColor(rt.kRed)
        timeHisto.Fit(fit)
        timeJitter = fit.GetParameter(2)
        errTimeJitter = fit.GetParError(2)
        timeHisto.SetFillColor(rt.kYellow)
        timeHisto.Draw()

        gaus1 = rt.TF1('gaus1', 'gaus(0)', timeBot, timeTop)
        gaus1.SetParameters(fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2))
        gaus2 = rt.TF1('gaus2', 'gaus(0)', timeBot, timeTop)
        gaus2.SetParameters(fit.GetParameter(3), fit.GetParameter(4), fit.GetParameter(5))
        #gaus3 = rt.TF1('gaus3', 'gaus(0)', timeBot, timeTop)
        #gaus3.SetParameters(fit.GetParameter(6), fit.GetParameter(7), fit.GetParameter(8))
        gaus1.SetLineStyle(7)
        gaus2.SetLineStyle(7)
        #gaus3.SetLineStyle(7)
        gaus1.Draw('same')
        gaus2.Draw('same')
        #gaus3.Draw('same')

        timeCanvas.Update()
        timeCanvas.Draw()
        fitBox = timeHisto.FindObject('stats')
        fitBox.SetTextColor(rt.kRed+2)
        fitBox.SetX1NDC(0.65)
        fitBox.SetY1NDC(0.25)
        fitBox.SetX2NDC(0.9)
        fitBox.SetY2NDC(0.70)

        latex = rt.TLatex()
        latex.SetTextSize(.035)
        latex.SetTextAlign(32)
        latex.DrawLatexNDC(.9, .88, 'Total time jitter %1.2f ps'%(timeJitter*1e3))
        print('Total time jitter %1.2f ps'%(timeJitter*1e3))

        if options.triggercharge:
            ftmJitter = (timeJitter**2-apdJitter**2)**0.5
            print('APD time jitter %1.2f ps'%(apdJitter*1e3))
            print('MCP-PMT time jitter %1.2f ps'%(ftmJitter*1e3))
            latex.DrawLatexNDC(.9, .82, 'APD time jitter %1.2f ps'%(apdJitter*1e3))
            latex.DrawLatexNDC(.9, .76, 'MCP-PMT time jitter %1.2f ps'%(ftmJitter*1e3))
        else: ftmJitter = timeJitter

        if options.label: root_style_ftm.labelRight(timeCanvas, options.label)
        root_style_ftm.labelFtm(timeCanvas)

        timeCanvas.SaveAs('%s/TimeJitter.eps'%(outdir))

        results['timeJitter'] = timeJitter
        results['errTimeJitter'] = errTimeJitter

    if options.charge and options.time:
        ''' plot psd vs arrival time 2D distribution '''
        timeChargeHisto = rt.TH2F('TimeChargeHisto', ';Signal time (ns);Charge;', timeBins, timeBot, timeTop, chargeBins, chargeBot, chargeTop)
        timePsdHisto = rt.TH2F('TimePsdHisto', ';Signal time (ns);PSD;', timeBins, timeBot, timeTop, psdBins, psdBot, psdTop)
        for timeTrigger,timeSignal,charge,psd in zip(timesTrigger,timesSignal,chargesLong,PSDs):
            t0 = timeSignal-timeTrigger
            timeChargeHisto.Fill(t0, charge)
            timePsdHisto.Fill(t0, psd)

        timeChargeCanvas = rt.TCanvas('TimeChargeCanvas', '', 800, 600)
        timeChargeHisto.Draw('colz')
        timeChargeCanvas.SaveAs('%s/TimeCharge.eps'%(outdir))
    
        timePsdCanvas = rt.TCanvas('TimePsdCanvas', '', 800, 600)
        timePsdHisto.Draw('colz')
        timePsdCanvas.SaveAs('%s/TimePsd.eps'%(outdir))
    
    # save file with results
    resultsdf = pd.DataFrame.from_dict(results, orient='index')
    resultsdf.to_csv('%s/results.csv'%(outdir))

if __name__=='__main__': main()