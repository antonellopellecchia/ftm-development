#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

import ROOT as rt
from physlibs.root import root_style_ftm

from ftm_analysis import femtoammeter

rt.gErrorIgnoreLevel = rt.kWarning

'''def readCurrentFile(path):
    with open(path) as f:
        currentsStrList = f.read().split('\t')
        currentsFloatList = [ float(s) for s in currentsStrList ]
        return np.array(currentsFloatList)'''

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--source')
    ap.add_argument('--drift')
    ap.add_argument('--anode')
    ap.add_argument('--ground')
    ap.add_argument('--low_ground')
    ap.add_argument('--low_anode')
    ap.add_argument('--low_drift')
    ap.add_argument('--power_scan')
    ap.add_argument('--cut', type=float)
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--energy', type=float)
    ap.add_argument('--gain', action='store_true')
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', action='store_true', help='ROOT batch mode')
    options = ap.parse_args(sys.argv[1:])
    optionsdict = dict(options._get_kwargs())

    driftGap = 0.5 # cm

    for d in ['.', 'currents']:
        try: os.makedirs(f'{options.out}/{d}')
        except: pass

    if not options.source or options.source=='laser':
        source='laser'
        if options.energy: fullEnergy = options.energy
        else: fullEnergy = 51
        laserRate = 100
    elif options.source=='xray': source='xray'

    ''' Ground current vs source power scan '''
    if options.power_scan:
        outdir = f'{options.out}/currents/power_scan'
        try: os.makedirs(outdir)
        except FileExistsError: pass

        pattern = r'(\d+)_(\d+)_(\d+)_(\w+)_ammeter.csv'
        parnames = [ 'drift', 'amplification', 'energy', 'source' ]
        scanVariable = 'energy'

        powerScanCurrents = femtoammeter.CurrentScan.FromDirectory(options.power_scan, pattern, parnames, scanVariable)
        opticalDensity = powerScanCurrents.GetParameter('OD')
        powerScanCurrents.SetScalex(1/100*10**(-opticalDensity)*fullEnergy)
        powerScanCurrents = abs(powerScanCurrents)
        powerScanCurrents.SaveCurrentPlots(outdir)

        powerScanCanvas = rt.TCanvas(f'PowerScanCanvas', '', 800, 600)
        powerScanGraph = powerScanCurrents.plot
        powerScanGraph.SetTitle(f';Source power (#muJ);Ground current (nA)')
        powerScanGraph.GetYaxis().SetTitleOffset(1.6)
        powerScanGraph.Draw('ap')

        powerScanFit = rt.TF1('powerScanFit', '[0]*pow(x,[1])', 0, 100)
        powerScanFit.SetParameters(0.1, 1)
        powerScanFit.SetParNames('k', '#alpha')
        powerScanGraph.Fit(powerScanFit)
        powerScanExponent = powerScanFit.GetParameter(1)

        root_style_ftm.labelFtm(powerScanCanvas)
        if options.label: root_style_ftm.labelRight(powerScanCanvas, options.label)
        powerScanCanvas.SetGrid()
        powerScanCanvas.SaveAs(f'{options.out}/PowerScan.eps')

    ''' Ground current vs amplification voltage scan '''
    if options.ground:
        outdir = f'{options.out}/currents/anode_scan'
        try: os.makedirs(outdir)
        except FileExistsError: pass

        pattern = r'(\d+)_(\d+)_(\d+)_(\w+)_ammeter.csv'
        parnames = [ 'drift', 'amplification', 'energy', 'source' ]
        scanVariable = 'amplification'

        groundCurrentScan = femtoammeter.CurrentScan.FromDirectory(options.ground, pattern, parnames, scanVariable)
        groundCurrentScan = abs(groundCurrentScan)
        groundCurrentScan.SaveCurrentPlots(outdir)

        groundCurrentCanvas = rt.TCanvas(f'GroundCurrentCanvas', '', 800, 600)
        groundCurrentGraph = groundCurrentScan.plot
        groundCurrentGraph.SetTitle(f';Amplification voltage (V);Ground current (nA)')
        groundCurrentGraph.GetYaxis().SetTitleOffset(1.6)
        groundCurrentGraph.Draw('ap')

        root_style_ftm.labelFtm(groundCurrentCanvas)
        if options.label: root_style_ftm.labelRight(groundCurrentCanvas, options.label)

        groundCurrentCanvas.SetGrid()
        groundCurrentCanvas.SetLogy()
        groundCurrentCanvas.SaveAs(f'{options.out}/GroundCurrent.eps')

    ''' Ground and anode current scans at low amplification fields '''
    lowFieldCurrentScans = dict()
    primaryCurrent, primaryCurrentError2 = 0, 0
    for electrode,directory in zip(['ground', 'anode', 'drift'],[options.low_ground, options.low_anode, options.low_drift]):
        if directory:
            outdir = f'{options.out}/currents/anode_scan_low/{electrode}'
            try: os.makedirs(outdir)
            except FileExistsError: pass

            pattern = r'(\d+)_(\d+)_(\d+)_(\w+)_ammeter.csv'
            parnames = [ 'drift', 'amplification', 'energy', 'source' ]
            scanVariable = 'amplification'

            lowFieldCurrentScans[electrode] = femtoammeter.CurrentScan.FromDirectory(directory, pattern, parnames, scanVariable)
            if electrode == 'drift':
                driftOffset = lowFieldCurrentScans[electrode].GetParameter('DRIFT') * driftGap * 1e3
                lowFieldCurrentScans[electrode].SetOffsetx(-driftOffset)
            lowFieldCurrentScans[electrode].SaveCurrentPlots(outdir)

            electrodeCurrentCanvas = rt.TCanvas(f'{electrode.capitalize()}CurrentCanvas', '', 800, 600)
            electrodeCurrentGraph = lowFieldCurrentScans[electrode].plot
            electrodeCurrentGraph.SetTitle(f';Amplification voltage (V);{electrode.capitalize()} current (nA)')
            electrodeCurrentGraph.GetYaxis().SetTitleOffset(1.6)
            electrodeCurrentGraph.Draw('ap')

            if electrode in ['ground', 'anode']:
                color = rt.kRed*(electrode=='ground')+rt.kBlue*(electrode=='anode')+2
                if electrode=='ground': x1,y1,x2,y2 = 0.2,0.15,0.5,0.35
                if electrode=='anode': x1,y1,x2,y2 = 0.2,0.7,0.5,0.9

                fit = rt.TF1('f', '[0]+[1]*exp([2]*x)', 10, 250)
                fit.SetParNames('i_{0}', 'k', '#alpha')
                fit.SetParameters(-0.002, -0.01, 0.1)
                fit.SetLineColor(color)
                electrodeCurrentGraph.Fit(fit, 'R')

                electrodeCurrentCanvas.Update()
                electrodeCurrentCanvas.Draw()
                fitBox = electrodeCurrentGraph.FindObject('stats')
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
                primaryCurrentError2 += fit.GetParError(0)**2

            root_style_ftm.labelFtm(electrodeCurrentCanvas)
            if options.label: root_style_ftm.labelRight(electrodeCurrentCanvas, options.label)

            electrodeCurrentCanvas.SetGrid()
            electrodeCurrentCanvas.SaveAs(f'{options.out}/LowField{electrode.capitalize()}Current.eps')

    primaryCurrent = abs(primaryCurrent)
    primaryCurrentError = primaryCurrentError2**0.5
    print(f'Primary current {primaryCurrent:1.2e} +/- {primaryCurrentError:1.2e} nA')

    ''' Gain calculation from currents '''
    if options.ground and options.low_anode and options.low_ground and options.power_scan:
        # we can calculate the gain: need primary current, power scan exponent and ground current
        primaryCurrentExtrapolated = primaryCurrent * (10/80)**powerScanExponent
        primaryCurrentExtrapolatedError = primaryCurrentError * (10/80)**powerScanExponent
        primaryCurrentErrorExtrapolatedRelative = primaryCurrentExtrapolatedError/primaryCurrentExtrapolated

        gainPlot = rt.TGraphErrors()
        gainPlot.SetName('Graph')
        for i,current in enumerate(groundCurrentScan):
            vampl = float(current.parameters['amplification'])
            gain = current.mean/primaryCurrentExtrapolated
            gainError = gain * ((current.error/current.mean)**2 + primaryCurrentErrorExtrapolatedRelative**2)**0.5
            gainPlot.SetPoint(i, vampl, gain)
            gainPlot.SetPointError(i, 0, gainError)

        gainCanvas = rt.TCanvas('GainCanvas', '', 800, 600)
        gainPlot.SetTitle(';Amplification voltage (V);Effective gas gain')
        gainPlot.GetYaxis().SetTitleOffset(1.5)
        gainPlot.GetYaxis().SetRangeUser(0.1, 1e3)
        color = rt.kRed+2
        gainPlot.SetLineColor(color)
        gainPlot.SetMarkerStyle(24)
        gainPlot.Draw('AP')

        #fit = rt.TF1('f', 'expo(0)+pol1(2)', 20, 500)
        fit = rt.TF1('f', 'expo(0)', 300, 600)
        fit.SetParNames('a', 'b', '#alpha', '#beta')
        fit.SetLineColor(color)
        gainPlot.Fit(fit, 'R')
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
        x1,y1,x2,y2 = 0.2,0.7,0.5,0.9
        fitBox.SetX1NDC(x1)
        fitBox.SetY1NDC(y1)
        fitBox.SetX2NDC(x2)
        fitBox.SetY2NDC(y2)

        gainCanvas.SetLogy()
        gainCanvas.SetGrid()
        root_style_ftm.labelFtm(gainCanvas)
        if options.label: root_style_ftm.labelRight(gainCanvas, options.label)

        legend = rt.TLegend(0.6, 0.15, 0.9, 0.28)
        legend.SetTextSize(0.033)
        legend.SetBorderSize(1)
        try: mixture = str(groundCurrentScan.GetParameter('GAS')[0])
        except AttributeError: mixture = 'Ar:CO_{2} 70:30'
        legend.AddEntry(gainPlot, mixture, 'p')
        legend.Draw()        

        gainPlot.SaveAs('%s/Gain.root'%(options.out))
        gainCanvas.SaveAs('%s/Gain.eps'%(options.out))

    ''' Make plot with all measurements together '''
    if options.low_anode and options.low_ground and options.low_drift:
        electrodes = ['ground', 'anode', 'drift']
        scans = [ lowFieldCurrentScans[el] for el in electrodes]
        multiScan = femtoammeter.MultiElectrodeScan(scans, electrodes)
        multiScan.EnableTotal()

        currentCanvas = rt.TCanvas('CurrentCanvas', '', 800, 600)
        legend = rt.TLegend(0.2, 0.85, 0.65, 0.88)
        legend.SetNColumns(4)
        multiScan.legend = legend
        multiScanPlot = multiScan.plot
        multiScanPlot.Draw('a')
        legend.Draw()

        currentCanvas.SetGrid()
        currentCanvas.SaveAs(f'{options.out}/Currents.eps')

if __name__=='__main__': main()