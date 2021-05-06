#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

import ROOT as rt
rt.gROOT.SetBatch(rt.kTRUE)
from physlibs.root import root_style_ftm

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input', nargs='+')
    ap.add_argument('--names', nargs='+')
    ap.add_argument('--legend', nargs='+')
    ap.add_argument('--colors', nargs='+')
    ap.add_argument('--lines', nargs='+', type=int)
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--wide', action='store_true')
    ap.add_argument('--noparams', action='store_true')
    ap.add_argument('--refit', action='store_true')
    options = ap.parse_args(sys.argv[1:])

    try: os.makedirs(options.out)
    except FileExistsError: pass

    gainFiles = [ rt.TFile(f) for f in options.input ]
    gainGraphs = [ f.Get('Graph') for f in gainFiles ]

    gainMultiGraph = rt.TMultiGraph()
    for i,g in enumerate(gainGraphs):
        jp = 0
        while jp<g.GetN():
            if g.GetPointX(jp)<200: g.RemovePoint(jp)
            else: jp += 1
        g.SetName(options.names[i])

        if options.refit:
            # fit plot again; otherwise, use fit function in root file
            fit = rt.TF1('f', 'expo(0)', 300, 500)
            g.Fit('f', 'R')
        #g.GetFunction('f').Delete()
        gainMultiGraph.Add(g, 'p')

    if options.wide: gainCanvas = rt.TCanvas('GainCanvas', '', 1000, 600)
    else: gainCanvas = rt.TCanvas('GainCanvas', '', 800, 600)
    gainCanvas.SetLogy()
    gainCanvas.SetGrid()

    gainMultiGraph.SetTitle(';Amplification voltage (V);Effective gain')
    gainMultiGraph.Draw('a')
    #gainMultiGraph.GetXaxis().SetRange(200, 400)

    #gainCanvas.Update()
    #gainCanvas.Draw()
    if options.wide:
        legend = rt.TLegend(0.49, 0.13, 0.9, 0.31)
        legend.SetNColumns(2)
    else: legend = rt.TLegend(0.19, 0.73, 0.5, 0.91)
    legend.SetTextSize(0.032)
    legend.SetBorderSize(1)
    rtcolors = list()
    for color in options.colors: rtcolors.append(root_style_ftm.colorsDict[color])
    #titles = ['Laser - inverted field method', 'X-rays - anode+ground current method']
    #x1,y1,x2,y2 = 0.19,0.93,0.46,1.03
    x1,y1,x2,y2 = 0.65,0.03,0.93,0.14
    for i,g in enumerate(gainGraphs):
        legend.AddEntry(g, options.legend[i], 'pl')
        y1 += 0.13
        y2 += 0.13
        g.SetMarkerColor(rtcolors[i])
        g.SetLineColor(rtcolors[i])
        g.SetMarkerStyle(24)
        g.GetFunction('f').SetLineColor(rtcolors[i])
        if options.lines:
            g.GetFunction('f').SetLineStyle(options.lines[i])
            g.SetLineStyle(options.lines[i])
        fitBox = g.FindObject('stats')
        if options.noparams:
            fitBox.Delete()
            rt.gStyle.SetOptFit(0)
        else:
            fitBox.SetTextSize(0.025)
            fitBox.SetTextColor(rtcolors[i])
            fitBox.SetLineColor(rtcolors[i])
            fitBox.SetFillColor(0)
            fitBox.SetFillStyle(1001)
            fitBox.SetBorderSize(1)
            fitBox.SetX1NDC(x1)
            fitBox.SetY1NDC(y1)
            fitBox.SetX2NDC(x2)
            fitBox.SetY2NDC(y2)
    legend.Draw()

    root_style_ftm.labelFtm(gainCanvas)
    if options.label: root_style_ftm.labelRight(gainCanvas, options.label)

    gainMultiGraph.SaveAs('%s/GainComparison.root'%(options.out))
    gainCanvas.SaveAs('%s/GainComparison.eps'%(options.out))

    gainDiscrepancy = abs(gainGraphs[0].Eval(400)-gainGraphs[1].Eval(400))/gainGraphs[1].Eval(400)
    print(gainDiscrepancy, 'discrepancy at 400 V')

if __name__=='__main__': main()
