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

    try: os.makedirs(options.out)
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
        for electrode in electrodes: electrodeCurrents[electrode], errElectrodeCurrents[electrode] = list(), list()

        for iFolder,currentFolder in enumerate(options.currents):
            
            currentSetupFile = pd.read_csv(currentFolder+'/setup.txt', sep='\t')
            try: anodeVoltage = currentSetupFile['VANODE'][0]
            except KeyError: anodeVoltage = 0
            if source=='laser':
                opticalDensity = currentSetupFile['OD'][0]
                attenuation = currentSetupFile['ATTENUATION'][0]
                energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
                print('Laser energy', energy, 'uJ')
            elif source=='xray':
                energy = currentSetupFile['CURRENT'][0]
                print('X-ray current', energy, 'uA')

            inputFiles = os.listdir(currentFolder)
            for electrode in electrodes:
                currentFilesOff = list()
                currentFilesOn = list()
                for inputFile in inputFiles:
                    match = re.match(r'\d+_(\w+)_(\w+).csv', inputFile)
                    if match:
                        if not match.group(2)==electrode: continue
                        if match.group(1)=='off': currentFilesOff.append(inputFile)
                        elif match.group(1)=='on': currentFilesOn.append(inputFile)
                getVoltage = lambda x: int(re.match(r'(\d+)_(\w+)_(\w+).csv',x).group(1))
                currentFilesOff = sorted(currentFilesOff, key=getVoltage)
                currentFilesOn = sorted(currentFilesOn, key=getVoltage)

                for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
                    matchOff = re.match(r'(\d+)_\w+_\w+.csv', currentFileOff)
                    matchOn = re.match(r'(\d+)_\w+_\w+.csv', currentFileOn)
                    if not matchOff.group(1)==matchOn.group(1): pass
                    voltage = float(matchOn.group(1))

                    currentPointsOff = scale[electrode]*readCurrentFile(currentFolder+'/'+currentFileOff)
                    currentPointsOn = scale[electrode]*readCurrentFile(currentFolder+'/'+currentFileOn)

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
        for i,electrode in enumerate(electrodes):
            currents, errCurrents = electrodeCurrents[electrode], errElectrodeCurrents[electrode]
            currentGraphs.append( rt.TGraphErrors(len(currents), anodeVoltages, currents, errAnodeVoltages, errCurrents) )
            g = currentGraphs[-1]
            g.SetMarkerColor(rt.kRed+i)
            legend.AddEntry(g, electrode, 'p')
            currentGraph.Add(g, 'p')

        #currentsSum = electrodeCurrents['drift']+electrodeCurrents['anode']+electrodeCurrents['ground']
        #errCurrentsSum = errElectrodeCurrents['drift']+errElectrodeCurrents['anode']+errElectrodeCurrents['ground']
        #currentSumGraph = rt.TGraphErrors(len(currents), anodeVoltages, currentsSum, errAnodeVoltages, errCurrentsSum)
        #legend.AddEntry()
        #currentGraph.Add(currentSumGraph)

        currentGraph.SetTitle(';Amplification voltage (V);Electrode current (nA)')
        currentGraph.Draw('a')
        currentCanvas.SetGrid()
        legend.Draw()
        currentCanvas.SaveAs(options.out+'/Currents.eps')

        electrode = 'ground'
        errElectrodeCurrents[electrode] /= max(abs(electrodeCurrents[electrode]))
        electrodeCurrents[electrode] /= max(abs(electrodeCurrents[electrode]))
        anodeVoltages /= 1e3*0.5
        electrodeCurrentCanvas = rt.TCanvas('electrodeCurrentCanvas', '', 800, 600)
        electrodeCurrentGraph = rt.TGraphErrors(len(anodeVoltages), anodeVoltages, abs(electrodeCurrents[electrode]), errAnodeVoltages, errElectrodeCurrents[electrode])
        electrodeCurrentGraph.SetTitle(';Drift field (kV/cm);Relative effective gain')
        electrodeCurrentGraph.SetMarkerColor(rt.kBlue+3)
        electrodeCurrentGraph.SetLineColor(rt.kBlue+3)
        electrodeCurrentGraph.SetMarkerStyle(24)
        electrodeCurrentGraph.Draw('ALP')
        #electrodeCurrentGraph.GetYaxis().SetRangeUser()
        #electrodeCurrentCanvas.SetLogy()
        electrodeCurrentCanvas.SetGrid()

        root_style_ftm.labelFtm(electrodeCurrentCanvas)
        if options.label: root_style_ftm.labelRight(electrodeCurrentCanvas, options.label)

        latex = rt.TLatex()
        latex.SetTextSize(0.035)
        latex.DrawLatexNDC(.6, .25, 'Ar:CO_{2} 70:30')
        latex.DrawLatexNDC(.6, .2, 'Amplification voltage %d V'%(anodeVoltage))

        electrodeCurrentCanvas.SaveAs('%s/GainDriftScan.eps'%(options.out))

if __name__=='__main__': main()