#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

import ROOT as rt
from physlibs.root import root_style_ftm

def readCurrentFile(path):
    with open(path) as f:
        currentsStrList = f.read().split('\t')
        currentsFloatList = [ float(s) for s in currentsStrList ]
        return np.array(currentsFloatList)

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--source')
    ap.add_argument('--drift')
    ap.add_argument('--anode')
    ap.add_argument('--ground')
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--energy', type=float)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', action='store_true', help='ROOT batch mode')
    options = ap.parse_args(sys.argv[1:])

    for d in ['.', 'currents']:
        try: os.makedirs(options.out+'/'+d)
        except: pass

    if not options.source or options.source=='laser':
        source='laser'
        if options.energy: fullEnergy = options.energy
        else: fullEnergy = 51
        laserRate = 100
    elif options.source=='xray': source='xray'

    currentFolders = [options.drift, options.anode, options.ground]
    electrodeCurrents, errElectrodeCurrents = dict(), dict()
    anodeVoltages, errAnodeVoltages = dict(), dict()
    electrodes = ['drift', 'anode', 'ground']
    for electrode in electrodes:
        electrodeCurrents[electrode], errElectrodeCurrents[electrode] = list(), list()
        anodeVoltages[electrode], errAnodeVoltages[electrode] = list(), list()

    for iFolder,electrodeFolder in enumerate(currentFolders):
        electrode = electrodes[iFolder]
        currentSetupFile = pd.read_csv(electrodeFolder+'/setup.txt', sep='\t')
        try: driftField = currentSetupFile['DRIFT'][0]
        except KeyError: driftField = 0
        if source=='laser':
            opticalDensity = currentSetupFile['OD'][0]
            attenuation = currentSetupFile['ATTENUATION'][0]
            energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
            print('Laser energy', energy, 'uJ')
        elif source=='xray':
            energy = currentSetupFile['CURRENT'][0]
            print('X-ray current', energy, 'uA')

        inputFiles = os.listdir(electrodeFolder)
        currentFilesOff = list()
        currentFilesOn = list()
        for inputFile in inputFiles:
            match = re.match(r'\d+_(\w+).csv', inputFile)
            if match:
                if match.group(1)=='off': currentFilesOff.append(inputFile)
                elif match.group(1)=='on': currentFilesOn.append(inputFile)
        getVoltage = lambda x: int(re.match(r'(\d+)_(\w+).csv',x).group(1))
        currentFilesOff = sorted(currentFilesOff, key=getVoltage)
        currentFilesOn = sorted(currentFilesOn, key=getVoltage)

        for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
            matchOff = re.match(r'(\d+)_\w+.csv', currentFileOff)
            matchOn = re.match(r'(\d+)_\w+.csv', currentFileOn)
            title = currentFileOn[:-4]
            if not matchOff.group(1)==matchOn.group(1): pass
            voltage = float(matchOn.group(1))
            #if voltage>300: continue

            currentPointsOff = 1e9*readCurrentFile(electrodeFolder+'/'+currentFileOff)
            currentPointsOn = 1e9*readCurrentFile(electrodeFolder+'/'+currentFileOn)
            times = np.linspace(0, len(currentPointsOn), len(currentPointsOn))

            currentTimeCanvas = rt.TCanvas('CurrentTimeCanvas'+str(i), '', 1000, 600)
            currentTimeMultiGraph = rt.TMultiGraph()
            currentTimeOnPlot = rt.TGraph(len(currentPointsOn), times, currentPointsOn)
            currentTimeOnPlot.SetLineColor(rt.kBlue)
            currentTimeOffPlot = rt.TGraph(len(currentPointsOff), times, currentPointsOff)
            currentTimeOffPlot.SetLineColor(rt.kRed)
            currentTimeMultiGraph.Add(currentTimeOnPlot, 'pl')
            currentTimeMultiGraph.Add(currentTimeOffPlot, 'pl')
            currentTimeMultiGraph.SetTitle(';Time (s);Current (nA)')
            currentTimeMultiGraph.Draw('a')
            currentTimeCanvas.SaveAs('%s/currents/%s/%s.eps'%(options.out,electrode,title))

            currentAverage = currentPointsOn.mean()-currentPointsOff.mean()
            currentError = ( currentPointsOn.std()**2/currentPointsOn.size + currentPointsOff.std()**2/currentPointsOff.size ) **0.5
            electrodeCurrents[electrode].append(currentAverage)
            errElectrodeCurrents[electrode].append(currentError)
            if not voltage in anodeVoltages[electrode]:
                anodeVoltages[electrode].append(voltage)
                errAnodeVoltages[electrode].append(0)

    currentGraphs = list()
    currentGraph = rt.TMultiGraph()
    legend = rt.TLegend(0.2, 0.15, 0.7, 0.2)
    legend.SetNColumns(4)

    for i,electrode in enumerate(electrodes):
        electrodeCurrents[electrode] = np.array(electrodeCurrents[electrode])
        errElectrodeCurrents[electrode] = np.array(errElectrodeCurrents[electrode])
        anodeVoltages[electrode], errAnodeVoltages[electrode] = np.array(anodeVoltages[electrode]), np.array(errAnodeVoltages[electrode])
        currents, errCurrents = electrodeCurrents[electrode], errElectrodeCurrents[electrode]
        currentGraphs.append( rt.TGraphErrors(len(currents), anodeVoltages[electrode], currents, errAnodeVoltages[electrode], errCurrents) )
        g = currentGraphs[-1]
        g.SetMarkerColor(rt.kRed+i*2)
        legend.AddEntry(g, electrode, 'p')
        currentGraph.Add(g, 'p')

    #electrodes.append('total')
    #electrodeCurrents['total'] = electrodeCurrents['drift']+electrodeCurrents['anode']+electrodeCurrents['ground']
    #errElectrodeCurrents['total'] = errElectrodeCurrents['drift']+errElectrodeCurrents['anode']+errElectrodeCurrents['ground']

    #currentsSum = electrodeCurrents['drift']+electrodeCurrents['anode']+electrodeCurrents['ground']
    #errCurrentsSum = errElectrodeCurrents['drift']+errElectrodeCurrents['anode']+errElectrodeCurrents['ground']
    #currentSumGraph = rt.TGraphErrors(len(currents), anodeVoltages, currentsSum, errAnodeVoltages, errCurrentsSum)
    #legend.AddEntry()
    #currentGraph.Add(currentSumGraph)


    totalCurrentGraph = rt.TGraph()
    voltages = anodeVoltages['drift']
    totalCurrents = list()
    for i,voltage in enumerate(voltages):
        totalCurrent = 0
        for j,electrode in enumerate(electrodes): totalCurrent += currentGraphs[j].Eval(voltage)
        totalCurrents.append(totalCurrent)
        totalCurrentGraph.SetPoint(i, float(voltage), totalCurrent)
    totalCurrentGraph.SetMarkerColor(rt.kBlue+3)
    totalCurrentCanvas = rt.TCanvas('TotalCurrentCanvas', '', 800, 600)
    totalCurrentGraph.SetTitle(';Amplification voltage (V);Total current (nA)')
    totalCurrentGraph.Draw('AP')
    totalCurrentCanvas.SetGrid()
    totalCurrentCanvas.SaveAs('%s/TotalCurrent.eps'%(options.out))

    currentCanvas = rt.TCanvas('CurrentCanvas', '', 800, 600)
    legend.AddEntry(totalCurrentGraph, 'total', 'p')
    currentGraph.Add(totalCurrentGraph, 'p')
    currentGraph.SetTitle(';Amplification voltage (V);Electrode current (nA)')
    currentGraph.Draw('a')
    currentCanvas.SetGrid()
    legend.Draw()

    root_style_ftm.labelFtm(currentCanvas)
    if options.label: root_style_ftm.labelRight(currentCanvas, options.label)

    currentCanvas.SaveAs(options.out+'/Currents.eps')

    primaryCurrent = 0
    for i,electrode in enumerate(electrodes):
        electrodeCurrentCanvas = rt.TCanvas('%sCurrentCanvas'%(electrode.capitalize()), '', 800, 600)
        graph = rt.TGraphErrors(len(electrodeCurrents[electrode]), anodeVoltages[electrode], electrodeCurrents[electrode], errAnodeVoltages[electrode], errElectrodeCurrents[electrode])
        graph.SetTitle(';Amplification voltage (V);%s current (nA)'%(electrode.capitalize()))
        graph.Draw('AP')
        graph.GetYaxis().SetTitleOffset(1.5)
        if electrode=='ground' or electrode=='anode':
            color = rt.kRed*(electrode=='ground')+rt.kBlue*(electrode=='anode')+2
            if electrode=='ground': x1,y1,x2,y2 = 0.2,0.15,0.5,0.35
            if electrode=='anode': x1,y1,x2,y2 = 0.2,0.7,0.5,0.9

            fit = rt.TF1('f', '[0]+[1]*exp([2]*x)', 10, 250)
            fit.SetParNames('i_{0}', 'k', '#alpha')
            fit.SetParameters(-0.002, -0.01, 0.1)
            fit.SetLineColor(color)
            graph.Fit(fit, 'R')

            electrodeCurrentCanvas.Update()
            electrodeCurrentCanvas.Draw()
            fitBox = graph.FindObject('stats')
            fitBox.SetTextColor(color)
            fitBox.SetLineColor(color)
            fitBox.SetFillColor(0)
            fitBox.SetFillStyle(1001)
            fitBox.SetBorderSize(1)
            fitBox.SetX1NDC(x1)
            fitBox.SetY1NDC(y1)
            fitBox.SetX2NDC(x2)
            fitBox.SetY2NDC(y2)

            primaryCurrent += fit.GetParameter(0)

        root_style_ftm.labelFtm(electrodeCurrentCanvas)
        if options.label: root_style_ftm.labelRight(electrodeCurrentCanvas, options.label)

        #graph.GetYaxis().SetRangeUser(1e-8, 1e2)
        #electrodeCurrentCanvas.SetLogy()
        electrodeCurrentCanvas.SetGrid()
        electrodeCurrentCanvas.SaveAs('%s/%sCurrent.eps'%(options.out, electrode.capitalize()))    
    print('Primary current', primaryCurrent, 'nA')

    gainCanvas = rt.TCanvas('GainCanvas', '', 800, 600)
    gain = electrodeCurrents['ground']/primaryCurrent
    errGain = errElectrodeCurrents['ground']/primaryCurrent
    gainPlot = rt.TGraphErrors(len(gain), anodeVoltages['ground'], gain, errAnodeVoltages['ground'], errGain)
    gainPlot.SetTitle(';Amplification voltage (V);Effective gas gain')
    gainPlot.GetYaxis().SetTitleOffset(1.5)
    color = rt.kRed+2
    gainPlot.SetLineColor(color)
    gainPlot.SetMarkerStyle(4)
    gainPlot.Draw('AP')

    fit = rt.TF1('f', '[0]*exp([1]*x) + [2]*x + [3]', 200, 500)
    fit = rt.TF1('f', 'expo(0)+pol1(2)', 20, 500)
    fit = rt.TF1('f', 'expo(0)', 250, 460)
    #fit.SetParNames('i_{0}', 'k', '#alpha')
    #fit.SetParameters(-5, 0.023, 0.01, .1)
    fit.SetLineColor(color)
    gainPlot.Fit(fit, 'R')

    '''fit = rt.TF1('f', '[0]*exp([1]*x) + [2]+[3]*x', 20, 460)
    fit.SetParameters(0.4317, 0.002653, -2.993, 0.01793)
    fit.SetParNames('a', 'b', '#alpha', '#beta')
    fit.SetLineColor(color)
    gainPlot.Fit(fit, 'R')'''

    fit.SetLineStyle(7)
    fit.SetRange(250, 600)
    fit.Draw('same')

    gainCanvas.Update()
    gainCanvas.Draw()
    fitBox = gainPlot.FindObject('stats')
    fitBox.SetTextColor(color)
    fitBox.SetLineColor(color)
    fitBox.SetFillColor(0)
    fitBox.SetFillStyle(1001)
    fitBox.SetBorderSize(1)
    x1,y1,x2,y2 = 0.2,0.75,0.5,0.9
    fitBox.SetX1NDC(x1)
    fitBox.SetY1NDC(y1)
    fitBox.SetX2NDC(x2)
    fitBox.SetY2NDC(y2)

    gainCanvas.SetLogy()
    gainCanvas.SetGrid()
    root_style_ftm.labelFtm(gainCanvas)
    if options.label: root_style_ftm.labelRight(gainCanvas, options.label)

    gainPlot.SaveAs('%s/Gain.root'%(options.out))
    gainCanvas.SaveAs('%s/Gain.eps'%(options.out))

if __name__=='__main__': main()
