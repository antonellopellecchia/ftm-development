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

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--calibration')
    ap.add_argument('--out')
    ap.add_argument('--flow', type=float)
    ap.add_argument('-b', action='store_true', help='ROOT batch mode')
    options = ap.parse_args(sys.argv[1:])

    try: os.makedirs(options.out)
    except FileExistsError: pass

    calibrationDf = pd.read_excel(options.calibration, skiprows=[1])
    mfcFlow = np.array(calibrationDf['MFC flow'], dtype=float)
    gasFlow = np.array(calibrationDf['Flow'])

    calibrationPlot = rt.TGraph(len(gasFlow), gasFlow, mfcFlow)
    calibrationPlot.SetTitle(';Gas flow (ccn/min);MFC measured flow (ccn/min)')

    fit = rt.TF1('fit', '[0]*(1-exp(-[1]*x))', 0, 80)
    fit.SetParameters(50, 10)
    calibrationPlot.Fit('pol4')
    calibrationPlot.SaveAs(f'{options.out}/Calibration.root')

    calibrationCanvas = rt.TCanvas('CalibrationCanvas', '', 800, 600)
    calibrationPlot.Draw('AP')
    calibrationCanvas.SaveAs(f'{options.out}/Calibration.eps')

    setFlow = calibrationPlot.Eval(options.flow)
    print(f'Set MFC flow to {setFlow:1.2f} ccn/min to get a real flow of {options.flow:1.2f} ccn/min')

if __name__=='__main__': main()