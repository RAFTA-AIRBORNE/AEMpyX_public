#!/usr/bin/env python3

"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by vrath  - June 11, 2022

"""
import time
import sys
from sys import exit as error
import os
import warnings
from time import process_time
from datetime import datetime



import numpy
from cycler import cycler
import matplotlib
import matplotlib.pyplot


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]


for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)


from version import versionstrg


import util
import core1d
import inverse
import aesys
import viz


warnings.simplefilter(action="ignore", category=FutureWarning)
rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=inspect.getfile(inspect.currentframe()), out=False)
print(titstrng+"\n\n")

OutInfo = False
now = datetime.now()


"""
System related settings.
Parameter can be transformed according to:
    ParaTrans   = 0           no change
                = 1           natural logarithm
Data transformation is now allowed with three possible options:
    DataTrans   = 0           raw data
               = 1           natural log of data
                = 2           asinh transformation
An error model is applied for the raw data, which is
mixed additive/multiplicative. in case of data transformation,
errors are also transformed.
"""
AEM_system = "aem05"
AEM_system = "genesis"

print("AEM system: " + AEM_system + "\n \n")

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 75.
    DatErr_mult = 0.03
    alt = 60.
    Data_Active = numpy.ones((1,NN[2]))

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    alt = 90.
    Data_Active = numpy.ones((1,NN[2]))

nD = NN[0]




UseSigma = True
condflag = "rho" #"sig"
if UseSigma:
    condflag = "sig" #"sig"

lay1rho = 10.
IdString = AEM_system.upper()+"_3Layer_res"+str(int(lay1rho))+"_"+condflag
PlotTitleSens = "sensitivity "+condflag+": resitivity = "+str(int(lay1rho))+" Ohm m , depth = 30 m"
PlotTitleModl = "model "+condflag+": resitivity = "+str(int(lay1rho))+" Ohm m , depth = 30 m"


OutDir = AEMPYX_DATA+"/Projects/SynthData/Jacobian/"
OutName =  OutDir+"/SYNTH_"


if not os.path.isdir(OutDir):
    print("File: %s does not exist, but will be created" % OutDir)
    os.mkdir(OutDir)



"""
Define model on mesh Note: dz should be small enough to resolve layer, or
should be adapted.
"""
nlyr = 60
er =   100.         # resistivity
mu = 1.0            # relative magnetic permeability
ep = 1.0            # relative dielectric constant
im = 1.e-14         # chargeability
it = 1.e-14         # time constant
ic = 1.0            # freq c constant
dz = 5.0



Model_active, Model_base, model_var, _, _ = inverse.init_1dmod(nlyr)


# set resistivity parameters active
Model_active[0 * nlyr:1 * nlyr] = 1


# fill model vector
Model_base = numpy.zeros(7 * nlyr) #.reshape(1,-1)
Model_base[0 * nlyr:1 * nlyr] = er
Model_base[1 * nlyr:2 * nlyr] = mu
Model_base[2 * nlyr:3 * nlyr] = ep
Model_base[3 * nlyr:4 * nlyr] = im
Model_base[4 * nlyr:5 * nlyr] = it
Model_base[5 * nlyr:6 * nlyr] = ic
Model_base[6 * nlyr:7 * nlyr - 1] = dz


dz = Model_base[6 * nlyr:7 * nlyr - 1]
zm = inverse.set_zcenters(dz)
zn = inverse.set_znodes(dz)
laythk = numpy.append(dz, 999.0)


PlotSens = True
PlotModl = False
if PlotSens or PlotModl:
    
    
    PlotFormat = [".pdf", ".png", ".svg"]
    PdfCatalog = False
    if ".pdf" in PlotFormat:
        PdfCName = OutName+IdString+"_Catalog.pdf"
    else:
        error(" No pdfs generated. No catalog possible!")
        PdfCatalog = False
    if PdfCatalog:
        pdf_list = []
    
    """
    Determine graphical parameter.
    => print(matplotlib.pyplot.style.available)
    """
    """
    For just plot.ing to files, choose the cairo backend (eps, pdf, png, jpg...).
    If you need to see the plots directly (plot window, or jupyter), simply
    comment OutInfo the following line. In this case matplotlib may run into
    memory problems after a few hundreds of high-resolution plot..
    """
    FilesOnly = False
    if FilesOnly:
        matplotlib.use("cairo")
        
    matplotlib.pyplot.style.use("seaborn-paper")
    matplotlib.rcParams["figure.dpi"] = 400
    matplotlib.rcParams["axes.linewidth"] = 0.5
    matplotlib.rcParams["savefig.facecolor"] = "none"
    matplotlib.rcParams["savefig.transparent"] = True
    
    Fontsize = 8
    Labelsize = Fontsize
    Titlesize = 8
    Fontsizes = [Fontsize, Labelsize, Titlesize]
    
    # Markersize = 4
    cm = 1/2.54
    FigWidth = 8.5*cm
    
    ncols = 256
    Colors = ["k", "r", "g", "b","m"]
    ColorMap="hot_r"
    
    Grey = 0.7
    
    Lines = (  cycler("linestyle", ["-", "--", ":", "-."])
             * cycler("color", ["r", "g", "b", "m"])
             * cycler("linewidth", [1.,]))




"""
Define models:
These are loops over different parameters, in this case for a 3-Layer case.
Should be adapted according to your needs.
"""

# IdString = AEM_system.upper()+"_5Layer"
# Thick0 = [25.]
# Thick1 = [10., 20.]
# Thick2 = [10., 20.]
# Thick3 = [10., 20.]
# Rho0 = [100.]
# Rho1 = [10., 1000]
# Rho2 = [100.]
# Rho3 = [10., 1000]
# Rhob = [100.]

# IdString = AEM_system.upper()+"_3Layer60"
Thick0 = [30.]
Thick1 = [30.]
Thick2 = [10.]
Thick3 = [10.]
Rho0 = [ 100.]
Rho1 = [lay1rho]   #[ 10.]
Rho2 = [100.]
Rho3 = [100.]
Rhob = [100.]


#Alt = [60., 120.]
Alt = [90.,]    #Alt = [90.,120., 180., 240.]
# Delta = [1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8]
# Delta = [-1.0]
Delta = [1.e-6]


"""
Generate Data

"""
Nsamples = 0
Perturb = False
SplitData= True
mod_num = -1
for alt in Alt:
    for thk0 in Thick0:
        for thk1 in Thick1:
            for thk2 in Thick2:
                for thk3 in Thick3:
                   for rho0 in Rho0:
                        for rho1 in Rho1:
                            for rho2 in Rho2:
                                for rho3 in Rho3:
                                    for rhob in Rhob:
                                        for delt in Delta:

                                            mod_num += 1
    
                                            p_i = [mod_num, alt,
                                                   thk0, thk1, thk2, thk3,
                                                   rho0, rho1,rho2, rho3, rhob,
                                                   DataTrans, DatErr_add, DatErr_mult]
    
                                            description =    "{0:2d} ".format(p_i[0])\
                                                            +"{0:.0f} ".format(p_i[1])\
                                                            +"{0:.0f} ".format(p_i[2])\
                                                            +"{0:.0f} ".format(p_i[3])\
                                                            +"{0:.0f} ".format(p_i[4])\
                                                            +"{0:.0f} ".format(p_i[5])\
                                                            +"{0:.0f} ".format(p_i[6])\
                                                            +"{0:.0f} ".format(p_i[7])\
                                                            +"{0:.0f} ".format(p_i[8])\
                                                            +"{0:.0f} ".format(p_i[9])\
                                                            +"{0:.0f} ".format(p_i[10])\
                                                            +"{0:2d} ".format(p_i[11])\
                                                            +"{0:.0f}/{1:.2f} ".format(p_i[12],p_i[13])
    
                                            print("\n\n\nModel description:")
                                            print(description)
    
    
                                            m_i = Model_base.copy()
                                                                                       
                                            d_lay = numpy.flipud(numpy.cumsum([thk0, thk1, thk2, thk3]))
                                            r_lay = numpy.flipud([rho0, rho1, rho2, rho3])
                                            
                                            # print(d_lay)
                                            # print(r_lay)
                                                
                                            mod_rho = rhob*numpy.ones_like(m_i[0 * nlyr:1 * nlyr])
                                            
                                            for ilay in numpy.arange(len(r_lay)):
                                                zl = d_lay[ilay]
                                                # print(zl, r_lay[ilay])
                                                # print(numpy.where(zm<zl))
                                                # print(numpy.shape(mod_rho))
                                                mod_rho[numpy.where(zm<zl)] = r_lay[ilay]
                                                
                                            m_i[0 * nlyr:1 * nlyr] =  mod_rho
                                            

    
                                            d_state = 0
                                            m_state = 0
                                            
                                            m_current, m_state = inverse.transform_parameter(m_vec=m_i, m_trn=ParaTrans, m_state=m_state, mode="f")               
                                            d_current, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=alt,
                                                                              m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                                                              d_trn=0, d_state=d_state, d_act = Data_Active )
                                            e_current, d_current = inverse.set_errors(d_current, DatErr_add, DatErr_mult, perturb=Perturb)
                            
                            
                                            jac = inverse.calc_jac(fwdcall=FwdCall, alt=alt,
                                                                             m_vec=m_current, m_act=Model_active, m_trn=ParaTrans, m_state=m_state,
                                                                             d_vec=d_current, d_trn=DataTrans, 
                                                                             delta=delt, out= False)
                                            jac = inverse.extract_jac(J=jac, d_act=Data_Active,  m_act=Model_active)

                                            # sensitivities
                                            Sens = []
                                            

                                            s0_current = inverse.calc_sensitivity(Jac=jac, UseSigma=UseSigma, Type = "raw") #[:-1]
                                            s0_current = inverse.transform_sensitivity(S=s0_current, V=laythk, 
                                                                                  Transform=[" val","max"])
                                            Sens.append(numpy.abs(s0_current))
                                            
                                            s1_current = inverse.calc_sensitivity(Jac=jac, UseSigma=UseSigma, Type = "cov") #[:-1]            
                                            s1_current = inverse.transform_sensitivity(S=s1_current, V=laythk, 
                                                                                  Transform=[" val","max"])
                                            Sens.append(numpy.abs(s1_current))
                                            
                                            s2_current = inverse.calc_sensitivity(Jac=jac, UseSigma=UseSigma, Type = "euc") #[:-1]   
                                            s2_current = inverse.transform_sensitivity(S=s2_current, V=laythk,  
                                                                                  Transform=[" val","max", "sqr"]) 
                                            Sens.append(numpy.abs(s2_current))
                                            

                                            d = numpy.diag(d_current)
                                            s3_current = inverse.calc_sensitivity(Jac=d@jac, UseSigma=UseSigma, Type = "cum") 
                                            s3_current = inverse.transform_sensitivity(S=s3_current, V=laythk, 
                                                                                 Transform=[" val","max"]) 
                                            Sens.append(numpy.abs(s3_current))  
                                            
                                            
    
                                            # print("\nData from this model:")
                                            # print(d_ref)
                                            
                                            
                                            if mod_num==0:
                                                Count = mod_num
                                                Model = m_i
                                                Sens0 = s0_current 
                                                Sens1 = s1_current
                                                Sens2 = s2_current 
                                                Sens3 = s3_current 
                                                Info = [description]
                                                Data = numpy.insert(d_current,0,[mod_num, -1, alt])
                                            else:
                                                Count = mod_num
                                                Model = numpy.vstack((Model, m_i))
                                                Sens0 = numpy.vstack((Sens0, s0_current)) 
                                                Sens1 = numpy.vstack((Sens1, s1_current)) 
                                                Sens2 = numpy.vstack((Sens2, s2_current)) 
                                                Sens3 = numpy.vstack((Sens3, s3_current)) 
                                                Info.append(description)
                                                Data =  numpy.vstack((Data, numpy.insert(d_current,0,[mod_num, -1, alt])))

                            
                                            Sens.pop(0)
                                                                                        
                                            if PlotModl:
                                                
                                                if delt>0: 
                                                    PlotFile = OutName+IdString+"_Delta"+str(int(numpy.log10(delt)))+"_model"
                                                else:
                                                    PlotFile = OutName+IdString+"_DeltaAdaptive"+"_model"

                                                viz.plot_model(
                                                                PlotFile = PlotFile,
                                                                PlotTitle = PlotTitleModl+", delta ="+str(delt),
                                                                PlotFormat = PlotFormat,                    
                                                                Depth = zn,
                                                                Model = [],
                                                                Error = [],
                                                                Labels=[],
                                                                Colors = Colors,
                                                                Lines = Lines,
                                                                Fontsizes = Fontsizes,
                                                                SLimits = [],
                                                                DLimits = [],
                                                                SUnits = "",
                                                                PlotStrng="", #Formula, #"", #"Error: mult="+str(DatErr_mult)+" add="+str(DatErr_add),
                                                                StrngPos=[0.05,0.05])
  

                                                if PdfCatalog:
                                                         pdf_list.append(PlotFile+".pdf")
                                                        
                                            if PlotSens:
                                                
                                                
                                                if delt>0: 
                                                    PlotFile = OutName+IdString+"_Delta"+str(int(numpy.log10(delt)))+"_sens"
                                                else:
                                                    PlotFile = OutName+IdString+"_DeltaAdaptive"+"_sens"        
                                                    
                                                
                                                viz.plot_sens(
                                                                PlotFile = PlotFile,
                                                                PlotTitle = PlotTitleSens+", delta ="+str(delt),
                                                                PlotFormat = PlotFormat,
                                                                Depth = zn,
                                                                Sens = Sens,
                                                                Labels = ["coverage", "euclidean","cumulative"],    #  "cummulative"
                                                                Colors = Colors,
                                                                Lines = Lines,
                                                                Fontsizes = Fontsizes,
                                                                SLimits = [],
                                                                DLimits = [],
                                                                SUnits = "",
                                                                PlotStrng="", #Formula, #"", #"Error: mult="+str(DatErr_mult)+" add="+str(DatErr_add),
                                                                StrngPos=[0.05,0.05])
                                                                                                        
                                                if PdfCatalog:
                                                        pdf_list.append(PlotFile+".pdf")
                            
 

                                            if delt>0: 
                                                NPZFile = OutName+IdString+"_Delta"+str(int(numpy.log10(delt)))+".npz"
                                            else:
                                                NPZFile = OutName+IdString+"_DeltaAdaptive"+".npz"

                                            
                                            
                                            print("Results written to "+NPZFile)
                                            numpy.savez_compressed(NPZFile, Model_active=Model_active, Model_base=Model_base, info=Info,
                                                                   model=Model, data=Data, sens0 = Sens0, sens1 = Sens1, sens2 = Sens2, Sens3 = Sens3)



if PdfCatalog:
    print("PDF catalogue: "+PdfCName)
    viz.make_pdf_catalog(PdfList=pdf_list, FileName=PdfCName)
