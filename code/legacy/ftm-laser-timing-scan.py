#!/usr/bin/python3

import os, sys
import numpy as np
import ROOT as rt
import pandas as pd
import argparse
import time, datetime
import re

from physlibs.root import root_style_ftm
from physlibs.root import functions

gainDatasheet = pd.read_csv('/home/anto/Documents/Dottorato/FTM/MCP-PMT/cosmics/gain-datasheet.csv', header=None, index_col=0).T

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input', nargs='+')
    ap.add_argument('--setup', nargs='+')
    ap.add_argument('--output')
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', help='Ignored', action='store_true')
    options = ap.parse_args(sys.argv[1:])

    try: os.makedirs(options.output)
    except FileExistsError: pass

    chargePlot = rt.TGraphErrors()
    primaryChargePlot = rt.TGraphErrors()
    timeJitterPlot = rt.TGraphErrors()

    laserFullEnergy = 51

    for i,(resultFile,setupFile) in enumerate(zip(options.input,options.setup)):
        setupdf = pd.read_csv(setupFile, sep='\t')
        gain = float(setupdf['GAIN'][0])
        attenuation = float(setupdf['ATT'][0])
        opticalDensity = float(setupdf['OD'][0])
        laserEnergy = laserFullEnergy*attenuation/100*10**(-opticalDensity)

        resultdf = pd.read_csv(resultFile, header=0, index_col=0).T
        charge = resultdf['charge'][0]
        errCharge = resultdf['errCharge'][0]
        primaryCharge = charge/gain/1.6e-7
        errPrimaryCharge = errCharge/gain/1.6e-7
        timeJitter = resultdf['timeJitter'][0]
        errTimeJitter = resultdf['errTimeJitter'][0]

        chargePlot.SetPoint(i, laserEnergy, charge)
        primaryChargePlot.SetPoint(i, laserEnergy, primaryCharge)
        chargePlot.SetPointError(i, 0., errCharge)
        primaryChargePlot.SetPointError(i, 0., errPrimaryCharge)
        timeJitterPlot.SetPoint(i, primaryCharge, timeJitter)
        timeJitterPlot.SetPointError(i, errPrimaryCharge, errTimeJitter)

    chargeCanvas = rt.TCanvas('ChargeCanvas', '', 800, 600)
    chargePlot.SetTitle(';Laser pulse energy (#muJ);Average collected charge (pC)')
    chargePlot.Draw('AP')
    pol2 = rt.TF1('p2', '[0]+[1]*x+pow(x,[2])', 0, 5)
    pol2.SetParameters(0, 1, 2)
    chargePlot.Fit(pol2)
    #chargeCanvas.SetLogy()
    chargeCanvas.SaveAs('%s/Charge.eps'%(options.output))

    primaryChargeCanvas = rt.TCanvas('PrimaryChargeCanvas', '', 800, 600)
    primaryChargePlot.SetTitle(';Laser pulse energy (#muJ);Average primary charge (electrons)')
    primaryChargePlot.Draw('AP')
    pol2 = rt.TF1('p2', '[0]+[1]*x+pow(x,[2])', 0, 5)
    pol2.SetParameters(0, 1, 2)
    rt.gStyle.SetOptFit(0)
    primaryChargePlot.Fit(pol2)
    #primaryChargeCanvas.SetLogy()
    primaryChargeCanvas.SaveAs('%s/PrimaryCharge.eps'%(options.output))

    timeJitterCanvas = rt.TCanvas('TimeJitterCanvas', '', 800, 600)
    timeJitterCanvas.SetLeftMargin(.12)
    timeJitterCanvas.SetRightMargin(.12)
    timeJitterPlot.SetTitle(';Number of primary electrons;Time resolution (ns)')
    timeJitterPlot.Draw('AP')
    f = rt.TF1('f', '[0]/pow(x,[1])+[2]', 0, 1)
    f.SetParNames('#sigma_{0}', '#alpha', 'offset')
    f.SetParameters(5.0, 0.5)
    timeJitterPlot.Fit(f)
    timeJitterCanvas.SaveAs('%s/TimeJitter.eps'%(options.output))


if __name__=='__main__': main()