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
    ap.add_argument('--currents', nargs='+')
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--energy', type=float)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', action='store_true', help='ROOT batch mode')
    options = ap.parse_args(sys.argv[1:])

    for d in ['', '/currents/drift', '/currents/anode', '/currents/ground']:
        try: os.makedirs(options.out+d)
        except: pass

    if not options.source or options.source=='laser':
        source='laser'
        if options.energy: fullEnergy = options.energy
        else: fullEnergy = 51
        laserRate = 100
    elif options.source=='xray': source='xray'

    if options.currents:
        electrodeCurrents, errElectrodeCurrents = dict(), dict()
        electrodes = ['drift', 'anode', 'ground']
        scale = {'drift':1e3, 'anode':1e3, 'ground':1e9}
        anodeVoltages, errAnodeVoltages = list(), list()
        powers, errPowers = list(), list()
        for electrode in electrodes: electrodeCurrents[electrode], errElectrodeCurrents[electrode] = list(), list()

        for iFolder,currentFolder in enumerate(options.currents):
            
            currentSetupFile = pd.read_csv(currentFolder+'/setup.txt', sep='\t')
            try: driftField = currentSetupFile['DRIFT'][0]
            except KeyError: driftField = 0
            if source=='laser':
                opticalDensity = currentSetupFile['OD'][0]
                try:
                    attenuation = currentSetupFile['ATTENUATION'][0]
                    energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
                    print('Laser energy', energy, 'uJ')
                except KeyError: pass
            elif source=='xray':
                energy = currentSetupFile['CURRENT'][0]
                print('X-ray current', energy, 'uA')

            inputFiles = os.listdir(currentFolder)

            # decide which variable is changing
            scanAnode,scanDrift,scanPower=list(),list(),list()
            for inputFile in inputFiles:
                match = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', inputFile)
                if match:
                    fdrift,fanode,fpower = match.group(1),match.group(2),match.group(3)
                    if not fanode in scanAnode: scanAnode.append(fanode)
                    if not fdrift in scanDrift: scanDrift.append(fdrift)
                    if not fpower in scanPower: scanPower.append(fpower)
            if len(scanAnode)>1: scanType='anode'
            elif len(scanDrift)>1: scanType='drift'
            elif len(scanPower)>1: scanType='power'
            print(scanType, 'scan')

            for electrode in electrodes:
                currentFilesOff = list()
                currentFilesOn = list()
                for inputFile in inputFiles:
                    match = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', inputFile)
                    if match:
                        fdrift,fanode,fpower,fstatus,felectrode = match.group(1),match.group(2),match.group(3),match.group(4),match.group(5)
                        if not felectrode==electrode: continue
                        if fstatus=='off': currentFilesOff.append(inputFile)
                        elif fstatus=='on': currentFilesOn.append(inputFile)
                if scanType=='anode':
                    getx = lambda x: int(re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv',x).group(2))
                elif scanType=='power':
                    getx = lambda x: int(re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv',x).group(3))
                currentFilesOff = sorted(currentFilesOff, key=getx)
                currentFilesOn = sorted(currentFilesOn, key=getx)

                for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
                    matchOff = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOff)
                    matchOn = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOn)
                    title = currentFileOn[:-4]
                    if not matchOff.group(2)==matchOn.group(2): pass
                    voltage = float(matchOn.group(2))
                    power = float(matchOn.group(3))

                    currentPointsOff = scale[electrode]*readCurrentFile(currentFolder+'/'+currentFileOff)
                    currentPointsOn = scale[electrode]*readCurrentFile(currentFolder+'/'+currentFileOn)
                    times = np.linspace(0, len(currentPointsOn), len(currentPointsOn))

                    currentTimeCanvas = rt.TCanvas('CurrentTimeCanvas'+str(i), '', 1000, 600)
                    currentTimeMultiGraph = rt.TMultiGraph()
                    graphs = [
                        rt.TGraph(len(currentPointsOn), times, currentPointsOn),
                        rt.TGraph(len(currentPointsOff), times, currentPointsOff)
                    ]
                    colors = [rt.kRed,rt.kBlue]
                    legend = rt.TLegend(.18, .95, .6, .99)
                    legend.SetNColumns(2)
                    l = ['on', 'off']
                    for j,g in enumerate(graphs):
                        g.SetLineColor(colors[j])
                        #g.SetMarkerColor(colors[j])
                        currentTimeMultiGraph.Add(g, 'pl')
                        legend.AddEntry(g, 'Source '+l[j], 'pl')
                    currentTimeMultiGraph.SetTitle(';Time (s);Current (nA)')
                    currentTimeMultiGraph.Draw('a')
                    legend.Draw()
                    latex = rt.TLatex()
                    latex.SetTextSize(.035)
                    latex.SetTextAlign(32)
                    latex.DrawLatexNDC(.95, .97, 'Amplification voltage %d V'%(voltage))
                    currentTimeCanvas.SaveAs('%s/currents/%s/%s.eps'%(options.out,electrode,title))

                    currentAverage = currentPointsOn.mean()-currentPointsOff.mean()
                    currentError = ( (currentPointsOn.std()/np.sqrt(currentPointsOn.size))**2 + (currentPointsOff.std()/np.sqrt(currentPointsOff.size))**2 ) **0.5
                    #currentError = ( (currentPointsOn.std())**2 + (currentPointsOff.std())**2 ) **0.5
                    if not options.drift and len(options.currents)>1:
                        currentAverage /= energy
                        currentError /= energy
                    electrodeCurrents[electrode].append(currentAverage)
                    errElectrodeCurrents[electrode].append(currentError)
                    if not voltage in anodeVoltages:
                        anodeVoltages.append(voltage)
                        errAnodeVoltages.append(0)
                    if not power in powers:
                        powers.append(power)
                        errPowers.append(0)

        currentCanvas = rt.TCanvas('CurrentCanvas', '', 800, 600)
        currentGraphs = list()
        currentGraph = rt.TMultiGraph()

        legend = rt.TLegend(0.2, 0.85, 0.65, 0.88)
        legend.SetNColumns(4)

        for electrode in electrodes:
            electrodeCurrents[electrode] = np.array(electrodeCurrents[electrode])
            errElectrodeCurrents[electrode] = np.array(errElectrodeCurrents[electrode])
        electrodes.append('total')
        electrodeCurrents['total'] = electrodeCurrents['drift']+electrodeCurrents['anode']+electrodeCurrents['ground']
        errElectrodeCurrents['total'] = errElectrodeCurrents['drift']+errElectrodeCurrents['anode']+errElectrodeCurrents['ground']
        anodeVoltages, errAnodeVoltages = np.array(anodeVoltages), np.array(errAnodeVoltages)
        powers, errPowers = np.array(powers), np.array(errPowers)
        for i,electrode in enumerate(electrodes):
            currents, errCurrents = electrodeCurrents[electrode], errElectrodeCurrents[electrode]
            if scanType=='anode':
                currentGraphs.append( rt.TGraphErrors(len(currents), anodeVoltages, currents, errAnodeVoltages, errCurrents) )
            if scanType=='power':
                currentGraphs.append( rt.TGraphErrors(len(currents), powers, currents, errPowers, errCurrents) )
            g = currentGraphs[-1]
            g.SetMarkerColor(rt.kRed+i)
            legend.AddEntry(g, electrode, 'p')
            currentGraph.Add(g, 'p')

        #currentsSum = electrodeCurrents['drift']+electrodeCurrents['anode']+electrodeCurrents['ground']
        #errCurrentsSum = errElectrodeCurrents['drift']+errElectrodeCurrents['anode']+errElectrodeCurrents['ground']
        #currentSumGraph = rt.TGraphErrors(len(currents), anodeVoltages, currentsSum, errAnodeVoltages, errCurrentsSum)
        #legend.AddEntry()
        #currentGraph.Add(currentSumGraph)

        if scanType=='anode': currentGraph.SetTitle(';Amplification voltage (V);Electrode current (nA)')
        if scanType=='power': currentGraph.SetTitle(';Laser pulse energy (#muJ);Electrode current (nA)')
        currentGraph.Draw('a')
        currentCanvas.SetGrid()
        legend.Draw()
        currentCanvas.SaveAs(options.out+'/Currents.eps')

        for electrode in electrodes:
            electrodeCurrentCanvas = rt.TCanvas('electrodeCurrentCanvas', '', 800, 600)
            if scanType=='anode':
                electrodeCurrentGraph = rt.TGraphErrors(len(currents), anodeVoltages, abs(electrodeCurrents[electrode]), errAnodeVoltages, errElectrodeCurrents[electrode])
                electrodeCurrentGraph.SetTitle(';Amplification voltage (V);%s current (nA)'%(electrode.capitalize()))
            if scanType=='power':
                electrodeCurrentGraph = rt.TGraphErrors(len(powers), powers, abs(electrodeCurrents[electrode]), errPowers, errElectrodeCurrents[electrode])
                electrodeCurrentGraph.SetTitle(';Laser pulse energy (#muJ);%s current (nA)'%(electrode.capitalize()))
            if electrode=='ground':
                fit = rt.TF1('f', 'pol3(0)', 20, 250)
                #electrodeCurrentGraph.Fit(fit, 'R')
            electrodeCurrentGraph.Draw('AP')
            #electrodeCurrentGraph.GetYaxis().SetRangeUser(1e-4, 1e2)
            if scanType=='anode': electrodeCurrentCanvas.SetLogy()
            electrodeCurrentCanvas.SetGrid()
            electrodeCurrentCanvas.SaveAs('%s/%sCurrent.eps'%(options.out, electrode.capitalize()))

if __name__=='__main__': main()
