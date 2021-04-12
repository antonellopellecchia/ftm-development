#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import pandas as pd

import ROOT as rt
from physlibs.root import root_style_ftm

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input', nargs='+')
    ap.add_argument('--names', nargs='+')
    ap.add_argument('--legend', nargs='+')
    ap.add_argument('--colors', nargs='+')
    ap.add_argument('--out', type=str)
    ap.add_argument('--label', type=str)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('-b', action='store_true', help='ROOT batch mode')
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
        #g.GetFunction('f').Delete()
        gainMultiGraph.Add(g, 'p')

    gainCanvas = rt.TCanvas('GainCanvas', '', 800, 600)
    gainCanvas.SetLogy()
    gainCanvas.SetGrid()

    gainMultiGraph.SetTitle(';Amplification voltage (V);Effective gain')
    gainMultiGraph.Draw('a')
    #gainMultiGraph.GetXaxis().SetRange(200, 400)

    #gainCanvas.Update()
    #gainCanvas.Draw()
    legend = rt.TLegend(0.19, 0.73, 0.5, 0.91)
    legend.SetTextSize(0.032)
    legend.SetBorderSize(1)
    rtcolors = list()
    for color in options.colors: rtcolors.append(root_style_ftm.colorsDict[color])
    #titles = ['Laser - inverted field method', 'X-rays - anode+ground current method']
    #x1,y1,x2,y2 = 0.19,0.93,0.46,1.03
    x1,y1,x2,y2 = 0.65,0.03,0.93,0.14
    for i,g in enumerate(gainGraphs):
        legend.AddEntry(g, options.legend[i], 'p')
        y1 += 0.13
        y2 += 0.13
        g.SetMarkerColor(rtcolors[i])
        g.SetMarkerStyle(24)
        g.SetLineColor(rt.kBlack)
        g.GetFunction('f').SetLineColor(rtcolors[i])
        fitBox = g.FindObject('stats')
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
