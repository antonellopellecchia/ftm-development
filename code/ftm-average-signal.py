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

from tqdm import tqdm

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input')
    ap.add_argument('--output')
    ap.add_argument('--polarity')
    ap.add_argument('--impedence', type=float)
    ap.add_argument('--channel', type=int)
    ap.add_argument('--times', type=float, nargs='+')
    ap.add_argument('--fcut', type=float)
    ap.add_argument('--label', type=str)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', help='Ignored', action='store_true')
    options = ap.parse_args(sys.argv[1:])

    try: os.makedirs(options.output)
    except FileExistsError: pass

    if not options.times: times = [0, 150]
    else: times = options.times

    if not options.polarity: polarity = 'negative'
    else: polarity = options.polarity
    negative = polarity=='negative'

    if not options.impedence: impedence = 50
    else: impedence = options.impedence

    files = sorted(os.listdir(options.input))
    signalFiles = list()
    for file in files:
        m = re.match(r'C%d--Trace--(\d+).trc'%(options.channel), file)
        if m: signalFiles.append(file)
    totalEvents = len(signalFiles)
    signalFiles = np.array(signalFiles)

    averageSignal = None
    goodSignalNumber = 0
    for signalNumber,signalFile in tqdm(enumerate(signalFiles)):
        if options.verbose: print(signalNumber,signalFile)
        filePath = '%s/%s'%(options.input, signalFile)
        scopeSignal = scope.ScopeSignal(filePath, scopeImpedence=impedence)
        try: scopeSignal.ReadSignal(negative=negative)
        except:
            print('Skipping', scopeSignal)
            continue
        if averageSignal is None: averageSignal = scopeSignal
        else: averageSignal += scopeSignal
        goodSignalNumber += 1
    averageSignal /= goodSignalNumber

    if options.fcut: averageSignal = averageSignal.GetFilteredSignal(options.fcut)
    averageCharge = averageSignal.GetChargeBetween(times[0], times[1])
    averageBaseline = averageSignal.GetBaseline()
    averageTime = averageSignal.GetArrivalTimeBySigmoidFit(findStart=True)

    signalCanvas = rt.TCanvas('AverageSignalCanvas', '', 820, 600)
    signalCanvas.SetRightMargin(1.2)
    signalGraph = averageSignal.GetGraph()
    signalGraph.GetXaxis().SetRangeUser(times[0]-20, times[1])
    signalGraph.Draw('AL')
    rt.gStyle.SetOptFit(0)

    if options.label: root_style_ftm.labelRight(signalCanvas, options.label)
    root_style_ftm.labelFtm(signalCanvas)
    latex = rt.TLatex()
    latex.SetTextSize(.03)
    latex.DrawLatexNDC(.7, .88, 'Charge %1.2f pC'%(averageCharge))
    latex.DrawLatexNDC(.7, .83, 'Baseline %1.2f mV'%(averageBaseline))
    latex.DrawLatexNDC(.7, .78, 'Time %1.2f ns'%(averageTime))
    
    signalCanvas.SaveAs(options.output+'/AverageSignal.eps')

if __name__=='__main__': main()