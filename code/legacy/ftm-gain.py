#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

def readCurrentFile(path):
    with open(path) as f:
        currentsStrList = f.read().split('\t')
        currentsFloatList = [ float(s) for s in currentsStrList ]
        return np.array(currentsFloatList)

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--gnd', type=str, nargs='+')
    ap.add_argument('--source', type=str)
    ap.add_argument('--drift', type=str)
    ap.add_argument('--drift_dark', type=str)
    ap.add_argument('--gain', action='store_true')
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--energy', type=float)
    ap.add_argument('--backend', type=str)
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

    if not options.backend or options.backend=='root':
        backend = 'root'
        import ROOT as rt
        from physlibs.root import root_style_ftm
    elif options.backend=='matplotlib':
        backend = 'matplotlib'
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        plt.style.use('classic')

    if options.drift:
        driftSetupFile = pd.read_csv(options.drift+'/setup.txt', sep='\t')
        if source=='laser': opticalDensity = driftSetupFile['OD'][0]
        driftVoltage = driftSetupFile['DRIFT'][0]

        inputFiles = os.listdir(options.drift)
        currentFilesOff = list()
        currentFilesOn = list()
        for inputFile in inputFiles:
            ''' need one measurement off for each measurement on
             due to changes in the offset current over long times: '''
            match = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', inputFile)
            if match and match.group(5)=='ground':
                if match.group(4)=='off': currentFilesOff.append(inputFile)
                elif match.group(4)=='on': currentFilesOn.append(inputFile)
        getPower = lambda x: int(re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv',x).group(3))
        currentFilesOff = sorted(currentFilesOff, key=getPower)
        currentFilesOn = sorted(currentFilesOn, key=getPower)

        driftCurrentGraph = rt.TGraphErrors()
        for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
            matchOff = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOff)
            matchOn = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOn)
            if not matchOff.group(3)==matchOn.group(3): pass

            if source=='laser':
                attenuation = float(matchOn.group(3))
                energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
            elif source=='xray':
                xrayCurrent = float(matchOn.group(3))

            currentPointsOff = -1e12*readCurrentFile(options.drift+'/'+currentFileOff)
            currentPointsOn = -1e12*readCurrentFile(options.drift+'/'+currentFileOn)

            currentAverage = currentPointsOn.mean()-currentPointsOff.mean()
            currentError = ( (currentPointsOn.std()/np.sqrt(currentPointsOn.size))**2 + currentPointsOff.std()**2/np.sqrt(currentPointsOff.size)) **0.5

            if source=='laser': driftCurrentGraph.SetPoint(i, energy, currentAverage)
            elif source=='xray': driftCurrentGraph.SetPoint(i, xrayCurrent, currentAverage)
            driftCurrentGraph.SetPointError(i, 0, currentError)

        driftCurrentCanvas = rt.TCanvas('DriftCurrentCanvas', '', 800, 600)
        if source=='laser': driftCurrentGraph.SetTitle(';Laser pulse energy (#muJ);Drift current (pA)')
        elif source=='xray': driftCurrentGraph.SetTitle(';X-ray current (#muA);Drift current (pA)')
        driftCurrentGraph.Draw('AP')
        if source=='laser': fit = rt.TF1('primaryCurrent', 'pol2', 0, 100)
        elif source=='xray': fit = rt.TF1('primaryCurrent', 'pol1', 0, 100)
        #driftCurrentGraph.Fit(fit)
        driftCurrentCanvas.SetGrid()
        #primaryCurrentFunction = driftCurrentGraph.GetFunction('primaryCurrent')
        #primaryCurrentFunction.SetLineStyle(7)

        driftCurrentCanvas.SaveAs(options.out+'/DriftCurrent.eps')

    '''if options.drift_dark:
        driftSetupFile = pd.read_csv(options.drift_dark+'/setup.txt', sep='\t')
        if source=='laser': opticalDensity = driftSetupFile['OD'][0]
        driftVoltage = driftSetupFile['DRIFT'][0]

        inputFiles = os.listdir(options.drift_dark)
        currentFilesOff = list()
        currentFilesOn = list()
        for inputFile in inputFiles:
            # need one measurement off for each measurement on
            # due to changes in the offset current over long times:
            match = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', inputFile)
            if match and match.group(5)=='ground':
                if match.group(4)=='off': currentFilesOff.append(inputFile)
                elif match.group(4)=='on': currentFilesOn.append(inputFile)
        getPower = lambda x: int(re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv',x).group(3))
        currentFilesOff = sorted(currentFilesOff, key=getPower)
        currentFilesOn = sorted(currentFilesOn, key=getPower)
        if len(currentFilesOff)==0: currentFilesOff = currentFilesOn

        driftCurrentGraphDark = rt.TGraphErrors()
        for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
            matchOff = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOff)
            matchOn = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOn)
            if not matchOff.group(1)==matchOn.group(1): pass

            if source=='laser':
                attenuation = float(matchOn.group(1))
                energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
            elif source=='xray':
                xrayCurrent = float(matchOn.group(1))

            currentPointsOff = -1e12*readCurrentFile(options.drift_dark+'/'+currentFileOff)
            currentPointsOn = -1e12*readCurrentFile(options.drift_dark+'/'+currentFileOn)

            if matchOff.group(2)==matchOn.group(2): currentAverage = currentPointsOn.mean()
            else: currentAverage = currentPointsOn.mean()-currentPointsOff.mean()
            currentError = ( (currentPointsOn.std()/np.sqrt(currentPointsOn.size))**2 + currentPointsOff.std()**2/np.sqrt(currentPointsOff.size)) **0.5

            if source=='laser': driftCurrentGraphDark.SetPoint(i, energy, currentAverage)
            elif source=='xray': driftCurrentGraphDark.SetPoint(i, xrayCurrent, currentAverage)
            driftCurrentGraphDark.SetPointError(i, 0, currentError)

        driftCurrentCanvasDark = rt.TCanvas('DriftCurrentCanvasDark', '', 800, 600)
        if source=='laser': driftCurrentGraphDark.SetTitle(';Laser pulse energy (#muJ);Drift current (pA)')
        elif source=='xray': driftCurrentGraphDark.SetTitle(';X-ray current (#muA);Drift current (pA)')
        driftCurrentCanvasDark.SetGrid()
        driftCurrentGraphDark.Draw('AP')
        #exp = rt.TF1('primaryCurrentDark', '[0] + [1]*(1-[2]*exp([3]*x))', 0, 100)
        #exp.SetParameters(0, 50, 1, 1)
        if options.drift:
            f = rt.TF1('primaryCurrentDark', 'pol5', 0, 100)
            driftCurrentGraphDark.Fit(f)
            primaryCurrentFunctionDark = driftCurrentGraphDark.GetFunction('primaryCurrentDark')
            primaryCurrentFunctionDark.SetLineStyle(7)
            primaryCurrentFunction = rt.TF1('primaryCurrentOn', 'primaryCurrent-primaryCurrentDark', 0, 50)
        driftCurrentCanvasDark.SaveAs(options.out+'/DriftCurrentDark.eps')

        if options.drift:
            driftCurrentCanvasOn = rt.TCanvas('DriftCurrentCanvasOn', '', 800, 600)
            primaryCurrentFunction.Draw()
            driftCurrentCanvasOn.SaveAs(options.out+'/DriftCurrentOn.eps')'''

    if options.gnd:
        gainArray, errGainArray = list(), list()
        gainVoltageArray, errGainVoltageArray = list(), list()

        gndCurrentArray, errGndCurrentArray = list(), list()
        voltageArray, errVoltageArray = list(), list()

        for iFolder,gndFolder in enumerate(options.gnd):
            gndSetupFile = pd.read_csv(gndFolder+'/setup.txt', sep='\t')
            driftField = float(gndSetupFile['DRIFT'][0])
            if source=='laser':
                opticalDensity = gndSetupFile['OD'][0]
                attenuation = gndSetupFile['ATTENUATION'][0]
                energy = fullEnergy*10**(-opticalDensity)*(1e-2*attenuation)
                print(opticalDensity, attenuation, fullEnergy)
                print('Laser energy', energy, 'uJ')
            elif source=='xray':
                energy = gndSetupFile['CURRENT'][0]
                print('X-ray current', energy, 'uA')
            if options.drift:
                #primaryCurrent = primaryCurrentFunction.Eval(energy)
                primaryCurrent = driftCurrentGraph.Eval(energy)
                print('Primary current', primaryCurrent, 'pA')

            inputFiles = os.listdir(gndFolder)
            currentFilesOff = list()
            currentFilesOn = list()
            for inputFile in inputFiles:
                match = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', inputFile)
                if match and match.group(5)=='ground':
                    if match.group(4)=='off': currentFilesOff.append(inputFile)
                    elif match.group(4)=='on': currentFilesOn.append(inputFile)
            getVoltage = lambda x: int(re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv',x).group(1))
            currentFilesOff = sorted(currentFilesOff, key=getVoltage)
            currentFilesOn = sorted(currentFilesOn, key=getVoltage)

            for i,(currentFileOff,currentFileOn) in enumerate(zip(currentFilesOff,currentFilesOn)):
                matchOff = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOff)
                matchOn = re.match(r'(\d+)_(\d+)_(\d+)_(\w+)_(\w+).csv', currentFileOn)
                if not matchOff.group(2)==matchOn.group(2): pass
                voltage = float(matchOn.group(2))

                currentPointsOff = -1e12*readCurrentFile(gndFolder+'/'+currentFileOff)
                currentPointsOn = -1e12*readCurrentFile(gndFolder+'/'+currentFileOn)

                currentAverage = currentPointsOn.mean()-currentPointsOff.mean()
                currentError = ( currentPointsOn.std()**2/currentPointsOn.size + currentPointsOff.std()**2/currentPointsOff.size ) **0.5
                #currentError = ( (currentPointsOn.std())**2 + (currentPointsOff.std())**2 ) **0.5
                if not options.drift and len(options.gnd)>1:
                    currentAverage /= energy
                    currentError /= energy
                gndCurrentArray.append(currentAverage)
                errGndCurrentArray.append(currentError)
                voltageArray.append(voltage)
                errVoltageArray.append(0)

                if options.drift:
                    gain = currentAverage/primaryCurrent
                    gainError = currentError/primaryCurrent
                    gainArray.append(gain)
                    errGainArray.append(gainError)
                    gainVoltageArray.append(voltage)
                    errGainVoltageArray.append(0)
                    #gainPlot.SetPoint(iGainPoint, voltage, gain)
                    #gainPlot.SetPointError(iGainPoint, 0, gainError)
                    #iGainPoint += 1

        gndCurrentArray, errGndCurrentArray = np.array(gndCurrentArray), np.array(errGndCurrentArray)
        voltageArray, errVoltageArray = np.array(voltageArray), np.array(errVoltageArray)
        if backend=='root':
            gndCurrentCanvas = rt.TCanvas('GndCurrentCanvas', '', 800, 600)
            gndCurrentGraph = rt.TGraphErrors(len(gndCurrentArray), voltageArray, gndCurrentArray, errVoltageArray, errGndCurrentArray)
            gndCurrentGraph.SetTitle(';Amplification voltage (V);Ground current (pA)')
            gndCurrentGraph.Draw('AP')
            gndCurrentCanvas.SetLogy()
            #gndCurrentGraph.GetYaxis().SetRangeUser(1e-2, 1e4)
            gndCurrentCanvas.SetGrid()
            gndCurrentCanvas.SaveAs(options.out+'/GndCurrent.eps')
        elif backend=='matplotlib':
            fig = plt.figure()
            plt.errorbar(voltageArray, gndCurrentArray, yerr=errGndCurrentArray, fmt='o-')
            plt.xlabel('Amplification voltage (V)')
            plt.ylabel('Ground current (pA)')
            plt.yscale('log')
            fig.savefig(options.out+'/GndCurrent.pdf')

        if options.gain:
            if not options.drift:
                gainArray = np.array(gndCurrentArray)/min(gndCurrentArray)
                errGainArray = np.array(errGndCurrentArray)/min(gndCurrentArray)
                gainVoltageArray, errGainVoltageArray = voltageArray, errVoltageArray
            gainArray, errGainArray = np.array(gainArray), np.array(errGainArray)
            gainVoltageArray, errGainVoltageArray = np.array(gainVoltageArray), np.array(errGainVoltageArray)
            if backend=='root':
                gainCanvas = rt.TCanvas('gainCanvas', '', 800, 600)
                gainPlot = rt.TGraphErrors(len(gainArray), gainVoltageArray, gainArray, errGainVoltageArray, errGainArray)
                gainPlot.SetTitle(';Amplification voltage (V);Effective gain')
                gainPlot.Draw('AP')
                color = rt.kRed+2
                #gainPlot.SetLineColor(color)
                gainPlot.SetMarkerStyle(24)
                gainPlot.GetXaxis().SetLimits(0, 600)
                gainPlot.GetYaxis().SetRangeUser(5e-2, 3.5e3)
                #gainPlot.Fit(expo, 'R')
                #gainPlot.GetYaxis().SetNdivisions(103)
                gainCanvas.SetLogy()
                gainCanvas.SetGrid()
                
                #fit = rt.TF1('exp', 'expo(0)', 200, 600)
                fit = rt.TF1('f', '[0]*exp([1]*x) + [2]+[3]*x', 20, 460)
                fit.SetParameters(0.4317, 0.002653, -2.993, 0.01793)
                fit.SetParNames('a', 'b', '#alpha', '#beta')
                fit.SetLineColor(color)
                gainPlot.Fit(fit, 'R')

                fit.SetLineStyle(7)
                fit.SetRange(0, 600)
                fit.Draw('same')
                print(fit.Eval(580), 'gain at 580 V')

                gainCanvas.Update()
                gainCanvas.Draw()
                fitBox = gainPlot.FindObject('stats')
                fitBox.SetTextColor(color)
                fitBox.SetLineColor(color)
                fitBox.SetFillColor(0)
                fitBox.SetFillStyle(1001)
                fitBox.SetBorderSize(1)
                x1,y1,x2,y2 = 0.2,0.7,0.5,0.9
                fitBox.SetX1NDC(x1)
                fitBox.SetY1NDC(y1)
                fitBox.SetX2NDC(x2)
                fitBox.SetY2NDC(y2)

                root_style_ftm.labelFtm(gainCanvas)
                if options.label: root_style_ftm.labelRight(gainCanvas, options.label)

                latex = rt.TLatex()
                latex.SetTextSize(0.03)
                latex.DrawLatexNDC(0.21, 0.62, 'Ar:CO_{2} 70:30 - Drift field %d kV/cm'%(driftField))

                gainPlot.SaveAs(options.out+'/Gain.root')
                gainCanvas.SaveAs(options.out+'/Gain.eps')
            elif backend=='matplotlib': pass

if __name__=='__main__': main()
