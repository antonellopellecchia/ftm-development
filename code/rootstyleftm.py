import ROOT as rt

rt.gROOT.SetBatch(True)

colorsDict = {
    'red': rt.kRed+2,
    'green': rt.kGreen+3,
    'blue': rt.kBlue+2,
    'orange': rt.kOrange+2,
    'magenta': rt.kMagenta+2,
}

def getColor(c): return colorsDict[c]

def labelLeft(canvas, label, x=.15, y=.97, divide=-1):
    if divide>0: canvas.cd(divide)
    else: canvas.cd()
    latex = rt.TLatex()
    latex.SetTextSize(.035)
    latex.SetTextAlign(12)
    latex.DrawLatexNDC(x, y, label)

def labelRight(canvas, label, x=.95, y=.97, divide=-1):
    if divide>0: canvas.cd(divide)
    else: canvas.cd()
    latex = rt.TLatex()
    latex.SetTextSize(.035)
    latex.SetTextAlign(32)
    latex.DrawLatexNDC(x, y, label)

def labelFtm(canvas, x=.15, y=.97, divide=-1):
    labelLeft(canvas, '#bf{FTM-Next}', x, y, divide)

def getStyle():
    style = rt.TStyle("FTM", "FTM Style")
    style.Reset()
    style.SetFillColor(0)
    style.SetFillStyle(1001)
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(0)
    style.SetCanvasPreferGL(rt.kTRUE)
    style.SetCanvasDefH(600)
    style.SetCanvasDefW(600)
    style.SetPadBorderMode(0)
    style.SetPadColor(0)
    style.SetPadLeftMargin(0.15)
    style.SetPadBottomMargin(0.1)
    style.SetPadRightMargin(0.05)
    style.SetPadTopMargin(0.06)
    style.SetPadTickX(0)
    style.SetPadTickY(0)
    style.SetFrameFillColor(0)
    style.SetFrameBorderMode(0)
    style.SetDrawBorder(0)
    style.SetLegendBorderSize(0)

    style.SetGridColor(rt.kGray)
    style.SetGridStyle(3)
    style.SetGridWidth(1)
    style.SetPadGridX(rt.kFALSE)
    style.SetPadGridY(rt.kFALSE)
    
    font = 42
    tsize = 0.04
    style.SetTextFont(font)
    style.SetTextSize(tsize)
    style.SetTitleStyle(0)
    style.SetTitleBorderSize(0)
    style.SetTitleColor(1, "xyz")
    style.SetTitleColor(1, "t")
    style.SetTitleFillColor(0)
    style.SetTitleFont(font, "xyz")
    style.SetTitleFont(font, "t")
    style.SetTitleOffset(1.2, "xyz")
    style.SetTitleOffset(1.3, "y")
    style.SetTitleSize(tsize, "xyz")
    style.SetTitleSize(tsize, "t")

    style.SetLegendFont(font)
    style.SetStatStyle(0)
    style.SetStatBorderSize(0)
    style.SetStatColor(0)
    style.SetStatFont(font)
    style.SetStatFontSize(tsize)
    style.SetStatX(0.58)
    style.SetStatY(0.88)
    style.SetStatW(0.2)
    style.SetStatH(0.1)
    style.SetOptStat(111110)
    style.SetOptFit(1)
    style.SetStatFormat("6.3g")
    style.SetLabelFont(font, "xyz")
    style.SetLabelSize(tsize, "xyz")
    style.SetLabelOffset(0.01, "xyz")
    style.SetOptTitle(0)
    style.SetPaperSize(rt.TStyle.kA4)
    style.SetFuncWidth(2)
    style.SetHistLineColor(rt.kRed - 3)
    style.SetPalette(1)
    style.SetAxisColor(rt.kBlack, "X")
    style.SetAxisColor(rt.kBlack, "Y")
    style.SetAxisColor(rt.kBlack, "Z")
    style.SetNdivisions(505, "x")
    style.SetNdivisions(510, "y")
    lw = 2
    style.SetLineWidth(lw)
    style.SetLineStyleString(2, "[12 12]")
    style.SetFrameLineWidth(lw)
    style.SetHistLineWidth(lw)
    style.SetFuncWidth(lw)
    style.SetFuncColor(rt.kRed-2)
    style.SetGridWidth(lw) 
    style.SetMarkerSize(1.1)
    style.SetMarkerStyle(20)
    style.SetMarkerColor(rt.kBlack)
    style.cd()
    
    return style

style = rt.TStyle("FTM", "FTMStyle")
style = getStyle()
style.cd()
rt.gROOT.SetStyle("FTM")
rt.gROOT.ForceStyle()