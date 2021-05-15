#!/usr/bin/python3

import os, sys
import numpy as np
import ROOT as rt
import pandas as pd
import argparse
import time, datetime
import re

import scope, femtoammeter
from physlibs.root import root_style_ftm as root_style
from physlibs.root import functions

rt.gROOT.SetBatch(rt.kTRUE)
rt.gErrorIgnoreLevel = rt.kWarning

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input')
    ap.add_argument('--output')
    ap.add_argument('--inputFormat')
    ap.add_argument('--outputFormat')
    ap.add_argument('--events', type=int, const=None)
    ap.add_argument('--chtrigger', type=int)
    ap.add_argument('--chdetector', type=int)
    ap.add_argument('--gainCurve')
    ap.add_argument('--primaryCurrent')
    ap.add_argument('--polarity')
    ap.add_argument('--setup')
    ap.add_argument('--draw')
    options = ap.parse_args(sys.argv[1:])

    options.input = options.input
    options.output = options.output
    try: os.makedirs(options.output)
    except FileExistsError: pass

    if not options.polarity: polarity = 'negative'
    else: polarity = options.polarity
    negative = polarity=='negative'

    parameters = ['drift', 'anode', 'power']
    #threshold = 0.7

    if options.gainCurve:
        # use gain to find primary current
        gainCurve = femtoammeter.GainCurve(options.gainCurve)
        if options.primaryCurrent:
            ''' use primary current scan to find primary charge
            this method is independent of the signal amplitude
            and should give you the same result if gain is the same
            measured in signal mode and current mode '''
            primaryCurrentCurve = femtoammeter.PrimaryCurrentCurve(options.primaryCurrent, gainCurve)

    if options.inputFormat=='raw':
        print('Input is raw...')
        scopeMeasurement = scope.Measurement.FromDirectory(options.input, parameters, options.chtrigger, options.chdetector, negative)
    elif options.inputFormat=='ntuples':
        scopeMeasurement = scope.Measurement.FromNtuples(options.input, parameters)

    if options.outputFormat=='ntuples':
        print('Converting to ntuples...')
        scopeMeasurement.ToNtuples(options.output, options.events)
    elif options.outputFormat=='draw':
        for scopeDataTaking in scopeMeasurement:
            print('Processing', scopeDataTaking)
            scopeDataTaking.DrawWaveforms(f'{options.output}/{scopeDataTaking.name}', nevents=options.events)
    elif options.outputFormat=='results':
        print('Saving efficiency plot...')
        # save efficiency as function of primary charge if gain is available,
        # instead as function of collected charge:
        if options.primaryCurrent and options.gainCurve:
            # if we have both, find primary charge from primary current:
            efficiencyPlot = scopeMeasurement.GetEfficiencyPlotPrimaries(primaryCurrentCurve, calibration='current')
        elif options.gainCurve:
            # if we have only gain curve, find primary charge from signal:
            efficiencyPlot = scopeMeasurement.GetEfficiencyPlotPrimaries(gainCurve, calibration='gain')
        else: # otherwise do not try to calculate primary charge:
            efficiencyPlot = scopeMeasurement.GetEfficiencyPlot()

        xmin, xmax = efficiencyPlot.GetXaxis().GetXmin(), efficiencyPlot.GetXaxis().GetXmax()
        efficiencySigmoid = functions.GetSigmoid(xmin, xmax, 1, 1, 0.5*(xmin+xmax), 0)
        efficiencySigmoid.SetLineColor(rt.kGreen+2)
        efficiencyPlot.Fit(efficiencySigmoid)
        
        efficiencyCanvas = rt.TCanvas('EfficiencyCanvas', '', 800, 600)
        efficiencyPlot.Draw('AP')

        efficiencyCanvas.Update()
        efficiencyCanvas.Draw()
        fitBox = efficiencyPlot.FindObject('stats')
        fitBox.SetTextColor(rt.kGreen+2)
        fitBox.SetLineColor(rt.kGreen+2)
        fitBox.SetFillColor(0)
        fitBox.SetFillStyle(1001)
        fitBox.SetBorderSize(1)
        x1,y1,x2,y2 = 0.6,0.15,0.9,0.35
        fitBox.SetX1NDC(x1)
        fitBox.SetY1NDC(y1)
        fitBox.SetX2NDC(x2)
        fitBox.SetY2NDC(y2)

        root_style.labelFtm(efficiencyCanvas)
        root_style.labelRight(efficiencyCanvas, 'FTM small-size')

        infoText = rt.TPaveText(0.19, 0.77, 0.5, 0.93, 'BL NDC')
        infoText.SetFillColorAlpha(0, 0)
        infoText.SetTextAlign(13)
        infoText.SetBorderSize(0)
        infoText.SetTextSize(.03)

        print(scopeMeasurement.setup)
        
        try:
            foilType = scopeMeasurement.setup['FOIL']
            infoText.AddText(f'{foilType} foil')
        except KeyError: pass
        try:
            mixture = scopeMeasurement.setup['GAS']
            infoText.AddText(f'Gas mixture {mixture}')
        except KeyError: pass
        try:
            driftField = scopeMeasurement.setup['DRIFT']
            infoText.AddText(f'Drift field {driftField} kV/cm')
        except KeyError: pass
        try:
            gain = scopeMeasurement.setup['AMPLIFICATION']
            infoText.AddText(f'Gain {gainCurve.GetGain(gain)/1e4:1.1f} #times 10^{{4}}')
        except KeyError: pass
        infoText.Draw()

        efficiencyCanvas.SaveAs(f'{options.output}/Efficiency.eps')
        print('Efficiency plot saved.')

    return

    if options.cache:
        cacheFile = rt.TFile(options.cache, 'READ')
        cacheTree = cacheFile.Get('EventTree')
        eventData, eventColumns = cacheTree.AsMatrix(return_labels=True)
        cacheDataFrame = pd.DataFrame(data=eventData, columns=eventColumns)

    if options.draw:
        for d in [options.draw+'/noise', options.draw+'/triggered', options.draw+'/signal', options.draw+'/pick']:
            try: os.makedirs(d)
            except FileExistsError: pass

    eventFile = rt.TFile('%s/signals.root'%(options.output), 'RECREATE')

    files = sorted(os.listdir(options.input))
    signalFilesTrigger = list()
    signalFilesDetector = list()
    for file in files:
        fileNameMatch = re.match(r'C(\d)--Trace--(\d+).trc', file)
        if not fileNameMatch: continue
        channel = int(fileNameMatch.group(1))
        traceNumber = fileNameMatch.group(2)
        otherChannels = [options.chtrigger,options.chdetector]
        if options.pick and not traceNumber==options.pick: continue 
        if not channel in otherChannels: continue
        otherChannels.remove(channel)
        exists = True
        for otherChannel in otherChannels:
            if not os.path.exists('%s/C%d--Trace--%s.trc'%(options.input, otherChannel, traceNumber)): exists = False
        if not exists: continue
        elif channel==options.chtrigger: signalFilesTrigger.append(file)
        elif channel==options.chdetector: signalFilesDetector.append(file)
        else: continue
    if options.sample:
        signalFilesTrigger = signalFilesTrigger[:options.sample]
        signalFilesDetector = signalFilesDetector[:options.sample]
    totalEvents = len(signalFilesDetector)
    signalFilesTrigger = np.array(signalFilesTrigger)
    signalFilesDetector = np.array(signalFilesDetector)

    eventTree = rt.TTree('EventTree', '')

    signalId = np.zeros(1, dtype='int')
    isEvent = np.zeros(1, dtype='bool')
    timeTrigger = np.zeros(1, dtype='float32')
    timeSignal = np.zeros(1, dtype='float32')
    chargeShort = np.zeros(1, dtype='float32')
    chargeLong = np.zeros(1, dtype='float32')
    chargeOut = np.zeros(1, dtype='float32')
    chargeTrigger = np.zeros(1, dtype='float32')
    psd = np.zeros(1, dtype='float32')
    riseTime = np.zeros(1, dtype='float32')
    amplitude = np.zeros(1, dtype='float32')
    sigmoidChi2 = np.zeros(1, dtype='float32')
    
    eventTree.Branch('signalId', signalId, 'signalId/I')
    eventTree.Branch('isEvent', isEvent, 'isEvent/I')
    eventTree.Branch('timeTrigger', timeTrigger, 'timeTrigger/F')
    eventTree.Branch('timeSignal', timeSignal, 'timeSignal/F')
    eventTree.Branch('chargeShort', chargeShort, 'chargeShort/F')
    eventTree.Branch('chargeLong', chargeLong, 'chargeLong/F')
    eventTree.Branch('chargeOut', chargeOut, 'chargeOut/F')
    eventTree.Branch('chargeTrigger', chargeTrigger, 'chargeTrigger/F')
    eventTree.Branch('psd', psd, 'psd/F')
    eventTree.Branch('riseTime', riseTime, 'riseTime/F')
    eventTree.Branch('amplitude', amplitude, 'amplitude/F')
    eventTree.Branch('sigmoidChi2', sigmoidChi2, 'sigmoidChi2/F')


    for signalNumber,(signalFileTrigger,signalFileDetector) in enumerate(zip(signalFilesTrigger,signalFilesDetector)):
        
        if options.verbose: print(signalNumber,signalFileTrigger,signalFileDetector)

        if options.cache:
            eventFilter = cacheDataFrame['signalId']==signalNumber
            if cacheDataFrame[eventFilter].size==0: continue # noise events not contained in cache tree
            if not bool(list(cacheDataFrame[eventFilter]['isEvent'])[0]): continue # triggered events with no MCP signal

        signalId[0] = signalNumber

        filePathTrigger = '%s/%s'%(options.input, signalFileTrigger)
        scopeSignalTrigger = scope.ScopeSignal(filePathTrigger)
        filePathDetector = '%s/%s'%(options.input, signalFileDetector)
        scopeSignalDetector = scope.ScopeSignal(filePathDetector)
        try:
            scopeSignalTrigger.ReadSignal()
            scopeSignalDetector.ReadSignal(negative=negative)
            scopeSignalDetector = scopeSignalDetector.GetFilteredSignal(0.2)
        except:
            print('Skipping', scopeSignalTrigger, scopeSignalDetector)
            continue

        try:
            timeTrigger[0] = scopeSignalTrigger.GetArrivalTimeBySigmoidFit(tstart=-2)
            noise = False
        except ValueError:
            timeTrigger[0] = 0
            noise = True
            noise = False

        if not noise:
            amplitude[0] = scopeSignalDetector.GetAmplitudeMax()-scopeSignalDetector.GetBaseline(-10)
            # 1.5 mV signal threshold
            #voltageThreshold = 1.7
            #if scopeSignalDetector.GetAmplitudeMax()<voltageThreshold: # APD triggered, but no MCP signal
            #    fitstatus = -1
            #    isEvent[0] = False
            #    timeSignal[0] = 0
            #else: # APD triggered and good MCP signal
            #frequencyCut = 0.2
            timeSignal[0], fitstatus, sigmoidChi2[0] = scopeSignalDetector.GetArrivalTimeBySigmoidFit(findStart=True, status=True, returnChi2=True)

            try:
                chargeTrigger[0] = scopeSignalTrigger.GetChargeBetween(timeTrigger[0]-2, timeTrigger[0]+20)
                chargeShort[0] = scopeSignalDetector.GetChargeBetween(timeSignal[0]-50, timeSignal[0]+5)
                #chargeLong[0] = scopeSignalDetector.GetChargeBetween(timeSignal[0]-50, timeSignal[0]+100)
                chargeLong[0] = scopeSignalDetector.GetChargeBetween(0, 200)
                chargeOut[0] = scopeSignalDetector.GetChargeBetween(150, 200)
                psd[0] = (chargeLong[0]-chargeShort[0])/chargeLong[0]
            except IndexError:
                isEvent[0] = False
                chargeTrigger[0] = 0
                chargeShort[0] = 0
                chargeLong[0] = 0
                chargeOut[0] = 0
                continue

            # discard both fit failed and events with chi2/ndf>1:
            if fitstatus==0:
                isEvent[0] = True
                eventTree.Fill()
            else:
                timeSignal[0] = 0
                isEvent[0] = False
        else: # APD not triggered, noise
            fitstatus = -1
            isEvent[0] = False


        if options.draw:
            sat = [timeTrigger[0], timeSignal[0]]
            tbot = [sat[0]-5, sat[1]-30]
            ttop = [sat[0]+5, sat[1]+150]
            signalCanvas = rt.TCanvas('signalCanvas'+signalFileDetector[2:], '', 1600, 600)
            signalCanvas.Divide(2,1)
            for i,scopeSignal in enumerate([scopeSignalTrigger,scopeSignalDetector]):
                signalCanvas.cd(i+1)
                #graph = scopeSignal.GetGraph(fcut=0.5)
                graph = scopeSignal.GetGraph()
                #if i==0: graph = scopeSignal.GetGraph()
                #else: graph = scopeSignal.GetFilteredSignal(frequencyCut)
                graph.GetXaxis().SetRangeUser(tbot[i], ttop[i])
                graph.Draw('AL')
                latex = rt.TLatex()
                latex.SetTextSize(.035)
                if i==1:
                    if isEvent[0]:
                        latex.DrawLatexNDC(.6, .85, 'Time %1.2f'%(sat[i]))
                        y0 = graph.GetFunction('sigmoid').Eval(sat[i])
                        line1 = rt.TLine(sat[i], -1, sat[i], amplitude[0])
                        line1.SetLineStyle(7)
                        line1.SetLineColor(rt.kRed)
                        line1.SetLineWidth(1)
                        line1.Draw()
                        line2 = rt.TLine(-sat[i]-20, y0, sat[i]+20, y0)
                        line2.SetLineStyle(7)
                        line2.SetLineColor(rt.kRed)
                        line2.SetLineWidth(1)
                        line2.Draw()
                    else: latex.DrawLatexNDC(.6, .85, 'Noise')
                    latex.DrawLatexNDC(.6, .8, 'Fit status %d'%(fitstatus))
                    latex.DrawLatexNDC(.6, .75, 'Charge %1.2f pC'%(chargeLong[0]))
                    latex.DrawLatexNDC(.6, .7, 'Amplitude %1.2f mV'%(amplitude[0]))
            #if options.pick: signalCanvas.SaveAs('%s/pick/%s.eps'%(options.draw, signalFileDetector[2:]))
            if noise: signalCanvas.SaveAs('%s/noise/%s.png'%(options.draw, signalFileDetector[2:]))
            elif not isEvent[0]: signalCanvas.SaveAs('%s/triggered/%s.png'%(options.draw, signalFileDetector[2:]))
            else: signalCanvas.SaveAs('%s/signal/%s.png'%(options.draw, signalFileDetector[2:]))
            #signalCanvas.SaveAs('%s/signal/%s.png'%(options.draw, signalFileDetector[2:]))

    eventTree.Print()
    eventFile.Write()
    eventFile.Close()

if __name__=='__main__': main()