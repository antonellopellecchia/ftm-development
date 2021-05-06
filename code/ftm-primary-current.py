#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

import ROOT as rt
from physlibs.root import root_style_ftm

import femtoammeter

rt.gErrorIgnoreLevel = rt.kWarning

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--drift', nargs='+')
    ap.add_argument('--mfcScales', nargs='+', type=float)
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    options = ap.parse_args(sys.argv[1:])
    optionsdict = dict(options._get_kwargs())

    driftGap = 0.5 # cm

    for d in ['.', 'currents']:
        try: os.makedirs(f'{options.out}/{d}')
        except FileExistsError: pass


    pattern = r'(\d+)_(\d+)_(\d+)_(\w+)_ammeter.csv'
    parnames = [ 'drift', 'amplification', 'energy', 'source' ]
    scanVariable = 'amplification'

    ''' Drift current scans at low amplification fields '''
    gasFractionList = list()
    primaryCurrents, primaryCurrentErrors = list(), list()
    for ipoint,directory in enumerate(options.drift):
        try: driftScan = femtoammeter.CurrentScan.FromDirectory(directory, pattern, parnames, scanVariable)
        except FileNotFoundError:
            print(f'{directory} not found.')
            break
        driftOffset = driftScan.GetParameter('DRIFT') * driftGap * 1e3
        driftScan.SetOffsetx(-driftOffset)

        gas = driftScan.GetParameter('GAS')[0]
        gasDir = gas.replace(' ', '_').replace(':', '_')
        try: os.makedirs(f'{options.out}/{gasDir}')
        except FileExistsError: pass
        
        gasFractions = gas.split(' ')[-1].split(':')
        gasFractions = [float(f) for f in gasFractions]
        gasFractions = [f/100*options.mfcScales[i] for i,f in enumerate(gasFractions)]
        gas1Fraction, gas2Fraction = gasFractions[0]/sum(gasFractions),gasFractions[1]/sum(gasFractions)
        gasFractionList.append(gas2Fraction)

        driftCurrentCanvas = rt.TCanvas(f'DriftCurrentCanvas', '', 800, 600)
        driftCurrentGraph = driftScan.plot
        driftCurrentGraph.SetTitle(f';Amplification voltage (V);Drift current (nA)')
        driftCurrentGraph.GetYaxis().SetTitleOffset(1.6)
        driftCurrentGraph.Draw('ap')

        color = rt.kGreen+2
        x1,y1,x2,y2 = 0.2,0.7,0.5,0.9

        fit = rt.TF1('f', '[0]+[1]*exp([2]*x)', 10, 250)
        fit.SetParNames('i_{0}', 'k', '#alpha')
        fit.SetParameters(-0.002, -0.01, 0.1)
        fit.SetLineColor(color)
        driftCurrentGraph.Fit(fit, 'R')

        driftCurrentCanvas.Update()
        driftCurrentCanvas.Draw()
        fitBox = driftCurrentGraph.FindObject('stats')
        fitBox.SetTextColor(color)
        fitBox.SetLineColor(color)
        fitBox.SetFillColor(0)
        fitBox.SetFillStyle(1001)
        fitBox.SetBorderSize(1)
        fitBox.SetX1NDC(x1)
        fitBox.SetY1NDC(y1)
        fitBox.SetX2NDC(x2)
        fitBox.SetY2NDC(y2)
            
        # primary current calculated from drift current fit:
        primaryCurrents.append(fit.GetParameter(0))
        primaryCurrentErrors.append(fit.GetParError(0))
        print(f'Primary current {primaryCurrents[-1]:1.2e} +/- {primaryCurrentErrors[-1]:1.2e} nA')

        root_style_ftm.labelFtm(driftCurrentCanvas)
        if options.label: root_style_ftm.labelRight(driftCurrentCanvas, options.label)

        driftCurrentCanvas.SetGrid()
        driftCurrentCanvas.SaveAs(f'{options.out}/{gasDir}/LowFieldDriftCurrent.eps')

    primaryCurrents = np.array(primaryCurrents)*1e3
    primaryCurrentErrors = np.array(primaryCurrentErrors)*1e3
    gasFractionList = np.array(gasFractionList)*100
    gasFractionErrorList = np.zeros(len(gasFractionList))

    primaryCurrentCanvas = rt.TCanvas('PrimaryCurrentCanvas', '', 800, 600)
    primaryCurrentPlot = rt.TGraphErrors(
        len(primaryCurrents), gasFractionList, primaryCurrents/primaryCurrents.mean(), gasFractionErrorList, primaryCurrentErrors/primaryCurrents)
    primaryCurrentPlot.SetTitle(';iC_{4}H_{10} fraction (%);Primary ionization current (pA)')
    primaryCurrentPlot.Draw('AP')
    primaryCurrentCanvas.SaveAs(f'{options.out}/PrimaryCurrent.eps')

if __name__=='__main__': main()