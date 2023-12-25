
# 半高宽拟合
def fwhmfit(erg):
    a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000        
    return (a+b*pow(erg+c*erg**2,0.5))*1000

# 高斯函数拟合
def gauss(x,x0,sigma,h):
    y=h*np.exp(-(x-x0)**2/(2*sigma**2))
    return y
def baseline(x,a,b):
    y=a*x+b
    return y
def spectrumcompose(x,x0,sigma,h,a,b):
    y=h*np.exp(-(x-x0)**2/(2*sigma**2))+a*x+b
    return y

# 本底+高斯峰拟合
def spectrumfit(list_erg,list_eff,peak,width=20):
    list_erg=list_erg[peak-width:peak+width+1]
    list_eff=list_eff[peak-width:peak+width+1]
    peak=width

    # xleft=max(peak-width,0)
    # xright=min(peak+width,len(list_eff))
    # yleft=list_eff[xleft]
    # yright=list_eff[xright]
    # a0=(yright-yleft)/(xright-xleft);b0=yleft-xleft*a0
    # erg0=list_erg[peak]
    # sigma0=abs(fwhmfit(erg0)/2.35)
    # amplitude0=list_eff[peak]-b0 
    # list_erg=list_erg[xleft:xright]
    # list_eff=list_eff[xleft:xright]

    try:
        erg,sigma,amplitude,a,b=opt.curve_fit(spectrumcompose,list_erg,list_eff,p0=[erg0,sigma0,amplitude0,a0,b0])[0]
    except:
        # erg,sigma,amplitude,a,b=erg0,sigma0,amplitude0,a0,b0
        return 'Unable to fit the peak, peak='
    list_fit=[gauss(x,erg,sigma,amplitude)+base(x,a,b) for x in list_erg]
    # plt.plot(list_erg,list_eff)
    # plt.savefig(foldpath+'\\spectrum_fit\\Unfitted_peak=%s.jpg'%(peak))
    # plt.plot(list_erg,list_fit)
    # plt.savefig(foldpath+'\\spectrum_fit\\fitted_peak=%s.jpg'%(peak))
    # plt.close()
    dict={'erg':erg,'sigma':abs(sigma),'amplitude':amplitude,'a':a,'b':b}
    return dict

def areaFit(list_erg,list_eff,peak,width=20):
    dict=spectrumfit(list_erg,list_eff,peak,width=20)
    area=pow(2*np.pi,0.5)*dict['amplitude']*dict['sigma']/0.29317
    return area