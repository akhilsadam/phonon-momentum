import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.optimize as opt
import scipy.ndimage as ndi
import statistics as stat
import sys
import tkinter as tk
from tkinter import filedialog
#----------------------------------------------------------
pi = 3.14159265358979323846
wStart = 0 #automatically set
wEnd = 12000
peak_min_width = 30 #required peak width to count as an absorptive peak
STD = 1.4      #required sigma to count as an absorptive peak
PTP = 2        #required peak-to-peak distance (in sigma) to count as an elastic peak
WTH = 0.09       #width multiplier for CI (confidence interval) to completely include peak
σS = 30       #smoothing for initial peak finding def:30
θErr =1.5    #angle offset (max absolute value in radians)
AMax = 1000       #max |amplitude| 
#----------------------------------------------------------
#root = tk.Tk()
#root.withdraw()

#file_path = filedialog.askopenfilename()
#file = open(file_path)
def fileOpen():
    try:
        file = open(input("Filename (no extensions!)")+'.txt')
        print("----------------------------------------------------------")
        return file
    except:
        print("ERROR: Incorrect Filename")
        file = fileOpen()
    if 'file' in locals():
        return file
    else:
        fileOpen()
file = fileOpen()
#----------------------------------------------------------
#DATA STRIP
def dataStrip():
    n = 0
    wStrt,wEn=0,0
    lines = []
    for x in file:  
        ### Ignore first two lines ##
        if n==0:
            ws = x.strip().split("kHz")
            wStrt =int(1000*float(ws[0]))
            wEn =int(1000*float(ws[1][(15):-1]))
            n +=1
        elif n==1:
            n +=1
        else:
            ### Remove '\n' from elements ###
            t = x.strip()
            temp = t.split('\t')
            lines.append(temp)
    #----------------------------------------------------------
    #DATA ARRAY
    wL,rL,xL,yL = [],[],[],[]
    for line in lines:
        wL.append(float(line[0]))
        rL.append(float(line[1]))
        xL.append(float(line[2]))
        yL.append(float(line[3]))
    w,r,x,y=np.array(wL),np.array(rL),np.array(xL),np.array(yL)
    return w,r,x,y,wStrt,wEn
#----------------------------------------------------------
w,r,x,y,wStart,wEnd = dataStrip() #//for some reason this doesn't work in output (wStart and end need to be global?)
#----------------------------------------------------------
#FIT BKG (LIN) for x,y
def bkg(w_in,m,b):
    return (w_in*m + b)
def smth(x,y):
    x_smth = ndi.gaussian_filter1d(x, σS)
    y_smth = ndi.gaussian_filter1d(y, σS)
    return x_smth,y_smth
#----------------------------------------------------------
#SHIFT x,y (first)
def firstShift(w,x,y):
    x_bkg_p, x_bkg_COV = opt.curve_fit(bkg, w, x)
    x_bkg = bkg(w,x_bkg_p[0],x_bkg_p[1])
    x_shift = x.copy()-x_bkg
    
    y_bkg_p, y_bkg_COV = opt.curve_fit(bkg, w, y)
    y_bkg = bkg(w,y_bkg_p[0],y_bkg_p[1])
    y_shift = y.copy()-y_bkg

    x_smth, y_smth = smth(x_shift,y_shift)
    return x_shift,y_shift,x_smth,y_smth
#----------------------------------------------------------
#PEAK FINDING (assumes y is arp-like and x is ela-like)
def arp_peaks(w,_smth):
    σ = stat.stdev(_smth)
    TRHLD = STD*σ
    _clean = [0]*len(w)
    for f in range(0,len(w)):
        if((_smth[f]>TRHLD) or (_smth[f]<-TRHLD)):
            _clean[f]=_smth[f]
    _peak_data = []#[[peak_data, final w],[peak,w]..]
    peak_data = []
    collect = False
    for f in range(0,len(w)-1):
        if(collect):
            peak_data.append(_clean[f])
        if(_clean[f]==0 and _clean[f+1]!=0):
            peak_data.clear()
            collect = True
        if(_clean[f]!=0 and _clean[f+1]==0):
            _peak_data.append([f+wStart-len(np.array(peak_data)),f+wStart])
            collect = False
    for peak in _peak_data:
        if(peak[1]-peak[0]<peak_min_width):
            _peak_data.remove(peak)
    return _peak_data
def ela_peaks(_smth,peaks):
    σ = stat.stdev(_smth)
    TRHLD = PTP*σ
    for arppeak in peaks:
        ela_like = _smth[(arppeak[0]-wStart):(arppeak[1]-wStart)]
        _max = max(ela_like)
        _min = min(ela_like)
        if(_max-_min<TRHLD):
            peaks.remove(arppeak)
    return peaks
def peak_filter(w_in,x_smth_in,y_smth_in):
    peaks = []
    for peak in ela_peaks(x_smth_in,arp_peaks(w_in,y_smth_in)):
        d = peak[1]-peak[0]
        sig = (WTH-1)/2
        peaks.append([int(peak[0]-sig*d),int(peak[1]+sig*d)])
    return peaks
def raw_peakless(w_in,x,y,peaks,wStart):
    w_pl = w_in.copy()
    x_pl = x.copy()
    y_pl = y.copy()
    for peak in peaks:
        srt = int(peak[0]-wStart)
        end = int(peak[1]-wStart)
        w1 = w_pl[0:srt]
        w2 = w_pl[end:len(w_pl)]
        x1 = x_pl[0:srt]
        x2 = x_pl[end:len(w_pl)]
        y1 = y_pl[0:srt]
        y2 = y_pl[end:len(w_pl)]
        w_pl = np.hstack([w1,w2])
        x_pl = np.hstack([x1,x2])
        y_pl = np.hstack([y1,y2])
    return w_pl,x_pl,y_pl
#----------------------------------------------------------
#SHIFT x,y (second)
def secondShift(w_in,x,y,peaks,wStart):
    w_pl,x_pl,y_pl = raw_peakless(w_in,x,y,peaks,wStart)
    x_bkg2_p, x_bkg2_COV = opt.curve_fit(bkg, w_pl, x_pl)
    x_bkg2 = bkg(w_in,x_bkg2_p[0],x_bkg2_p[1])
    x_shift2 = x.copy()-x_bkg2
    
    y_bkg2_p, y_bkg2_COV = opt.curve_fit(bkg, w_pl, y_pl)
    y_bkg2 = bkg(w_in,y_bkg2_p[0],y_bkg2_p[1])
    y_shift2 = y.copy()-y_bkg2
    return x_shift2,y_shift2
#----------------------------------------------------------
#DEF FITS x,y
# λ is damping parameter, w0 is resonance frequency, θ is phase shift
# y = A cos θ - B sin θ
# x = A sin θ + B cos θ
def ela(w_in, λ, w0, A, sgn):#elastic (note w is not really f! so lambda needs to be divided by 2pi!)
    return sgn*A*(pow(w0,2) - pow(w_in,2))/(pow((pow(w0,2) - pow(w_in,2)),2) + pow(λ*w_in,2))
def arp(w_in, λ, w0, A, sgn):#absorptive
    return sgn*A*λ*w_in/(pow((pow(w0,2) - pow(w_in,2)),2) + pow(λ*w_in,2)) 
def y_fit(w_in, λ, w0, θ, A, sgn):#arp-like
    return arp(w_in,λ,w0,A,sgn)*math.cos(θ)-ela(w_in,λ,w0,A,sgn)*math.sin(θ)
def x_fit(w_in, λ, w0, θ, A, sgn):#ela-like
    return arp(w_in,λ,w0,A,sgn)*math.sin(θ)+ela(w_in,λ,w0,A,sgn)*math.cos(θ)
#stack fits
def posfitStack(w_in, λ, w0, θ, A):
    return np.hstack([x_fit(w_in,λ,w0,θ,A,1),y_fit(w_in,λ,w0,θ,A,1)])
def negfitStack(w_in, λ, w0, θ, A):
    return np.hstack([x_fit(w_in,λ,w0,θ,A,-1),y_fit(w_in,λ,w0,θ,A,-1)])
#stack data
def dataStack(w_in,x_shift2,y_shift2,wStart):
    srt = int(w_in[0])
    end = int(w_in[len(w_in)-1]+1)
    return np.hstack([(x_shift2[srt-wStart:end-wStart]), (y_shift2[srt-wStart:end-wStart])])
def Rsquare(obs,exp):
    sig = stat.stdev(exp)
    return 1-(np.sum((obs-exp)**2)/sig)
#----------------------------------------------------------
#FIT elastic and absorptive  #,bounds=([0,0,-θErr,-AMax], [1000, wEnd, θErr, AMax]))
def fitPeak(w_in,x_shift2,y_shift2,in_peaks, i, wStart,y_smth):
    wPeak = w_in[(in_peaks[i][0]-wStart):(in_peaks[i][1]-wStart)]
    #Guesses
    λguess = 100
    θguess = 0.5
    wGuess = 2140#(w_in[(in_peaks[i][0]-wStart)]+w_in[(in_peaks[i][1]-wStart)])/2
    amp_guess = y_smth[(in_peaks[i][0]-wStart)]
    p0 =[λguess,wGuess,θguess,amp_guess]
    #bounds=([50,4000,0.0001,0], [200, 12000, 1,10])
    #-------
    if(amp_guess<0): #y_smth[(_peaks[i][0]-wStart)]
        fit_p, fit_COV = opt.curve_fit(posfitStack, wPeak, dataStack(wPeak,x_shift2,y_shift2,wStart),p0=p0)
        sgn=1
    else:
        fit_p, fit_COV = opt.curve_fit(negfitStack, wPeak, dataStack(wPeak,x_shift2,y_shift2,wStart),p0=p0)
        sgn=-1
    fit_p[1]=np.abs(fit_p[1])
    if(fit_p[0]<0):
        print("ERROR - DAMPING NEGATIVE: Fitting has failed to produce realistic results")
        #sys.exit()
    x_fitline = x_fit(w_in,fit_p[0],fit_p[1],fit_p[2],fit_p[3],sgn)
    y_fitline = y_fit(w_in,fit_p[0],fit_p[1],fit_p[2],fit_p[3],sgn)
    Rsquared = Rsquare(dataStack(wPeak,x_fitline,y_fitline,wStart),dataStack(wPeak,x_shift2,y_shift2,wStart))
    return x_fitline,y_fitline,fit_p,fit_COV,Rsquared
#----------------------------------------------------------
def outVal(fit_p,n,t,rsquared):
    print(" Peak",(n+1),"of",t,"at",fit_p[1],"Hz | FIT Goodness: RSquared = ",rsquared)
    print("----------------------------------------------------------")
    print("   ABSORPTION AMPLITUDE: ",arp(fit_p[1],fit_p[0],fit_p[1],fit_p[2],fit_p[3]),"V")
    print("   ELASTIC AMPLITUDE:    ",arp(fit_p[1],fit_p[0],fit_p[1],fit_p[2],fit_p[3]),"V")
    #print("   -------------------")
    print("   RESONANCE AMPLITUDE:  ",fit_p[3],"V")
    print("   RESONANCE QUALITY:    ",(fit_p[1]*pi/fit_p[0]))
    print("   -------------------")
    print("   phase shift:          ",(fit_p[2]), "radians")
    #print("   --------------------")
    return
def plot(w_in,x_fitline,y_fitline,x_shifted,y_shifted):
    fig, axs = plt.subplots(2)
    axs[0].plot(w_in,x_shifted)
    axs[0].plot(w_in,x_fitline)
    axs[0].set_xlabel("Frequency (Hz)")
    axs[0].set_ylabel("X Amplitude (V)")
    axs[1].plot(w_in,y_shifted)
    axs[1].plot(w_in,y_fitline)
    axs[1].set_xlabel("Frequency (Hz)")
    axs[1].set_ylabel("Y Amplitude (V)")
    plt.show()
    return
#----------------------------------------------------------
def outputResonances(w_in,x,y,x_shift,y_shift,x_shift2,y_shift2,_peaks, iteration,wStart,y_smth):
    t = len(_peaks)
    for n in range(0,t):
        rsquared,rsquared2 = 0,0;
        try:
            x_fitline,y_fitline,fit_p,COV_MX,rsquared = fitPeak(w_in,x_shift,y_shift,_peaks,n,wStart,y_smth)
        except:
            print("ERROR: Fitting with average linear background has failed to converge")
        try:
            x_fitline2,y_fitline2,fit_p2,COV_MX2,rsquared2 = fitPeak(w_in,x_shift2,y_shift2,_peaks,n,wStart,y_smth)
        except:
            print("ERROR: Fitting with peakless linear background has failed to converge")
        if(rsquared==0 and rsquared2==0):
            print("ERROR: Both fits failed to converge.")
            if(iteration<2):
                retry = input("Do you wish to try with a smoother version? Y/N:")
                if(retry==("Y")):
                    print("Retrying...")
                    print("UNFIT data:")
                    x_shift_sm,y_shift_sm = smth(x_shift,y_shift)
                    x_shift_sm2,y_shift_sm2 = smth(x_shift2,y_shift2)
                    plt.plot(w,x_shift_sm2)
                    plt.plot(w,y_shift_sm2)
                    plt.show()
                    outputResonances(w_in,x,y,x_shift_sm,y_shift_sm,x_shift_sm2,y_shift_sm2,_peaks, (iteration+1),wStart,y_shift_sm)
                return
        elif(rsquared==0):
            print("\n/---------------------------FIT 2---peakless---------------------------------/")
            outVal(fit_p2,n,t,rsquared2)
            plot(w_in,x_fitline2,y_fitline2,x_shift2,y_shift2)
            print(COV_MX2)
        elif(rsquared2==0):
            print("\n/---------------------------FIT 1---average----------------------------------/")
            outVal(fit_p,n,t,rsquared)
            plot(w_in,x_fitline,y_fitline,x_shift,y_shift)
            print(COV_MX)
        else:
            print("\n/---------------------------FIT 1---average----------------------------------/")
            outVal(fit_p,n,t,rsquared)
            plot(w_in,x_fitline,y_fitline,x_shift,y_shift)
            print(COV_MX)
            print("\n/---------------------------FIT 2---peakless---------------------------------/")
            outVal(fit_p2,n,t,rsquared2)
            plot(w_in,x_fitline2,y_fitline2,x_shift2,y_shift2)
            print(COV_MX2)        
        print("----------------------------------------------------------------------------//")
    return
def output(w,x,y,wStart):
    x_shift,y_shift,x_smth,y_smth = firstShift(w,x,y)
    _peaks = peak_filter(w,x_smth,y_smth)
    x_shift2,y_shift2 = secondShift(w,x,y,_peaks,wStart)
    print("Resonances Found: ",len(_peaks))
    print("UNFIT data:")
    plt.plot(w,x_shift2)
    plt.plot(w,y_shift2)
    for peak in _peaks:
        plt.axvspan(peak[0], peak[1], color='red', alpha=0.1)
    plt.show()
    print("----------------------------------------------------------------------------//")
    if(len(_peaks)<1):
        retry = input("Do you wish to try with a smoother version? Y/N:")
        if(retry==("Y")):
            print("Retrying...")
            x_s,y_s = smth(x,y)
            output(w,x_s,y_s,wStart)
        return
    outputResonances(w,x,y,x_shift,y_shift,x_shift2,y_shift2,_peaks,0,wStart,y_smth)
    return
#----------------------------------------------------------
output(w,x,y,wStart)