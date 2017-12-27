import ROOT as r
import math
import os,sys,getopt
import rootUtils as ut
import shipunit as u
inputFile  = None
geoFile    = None
dy         = None
import shipVeto
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
import shipRoot_conf
shipRoot_conf.configure()
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3d




try:
        opts, args = getopt.getopt(sys.argv[1:], "n:f:g:A:Y:i", ["nEvents=","geoFile="])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter file name'
        sys.exit()
for o, a in opts:
        if o in ("-f",):
            inputFile = a
        if o in ("-g", "--geoFile",):
            geoFile = a
        if o in ("-Y",):
            dy = float(a)


countSBT=0
countSVT=0
countUVT=0
countRPC=0

FidCutEv=0
DocaCutEv=0
hdocaIPCut=0
IPCutEv=0#define histograms
SBTCutEv=0
selectedWtihoutSBT=0
selectedEv=0


hNsegmentsDigi=r.TH1D('hNsegmentsDigi','Number of segments per event',300,0.,300)
hNsegments=r.TH1D('hNsegments','Number of segments per event',300,0.,300)
Nhits_vs_Nparticle=r.TH2D('Nhits_vs_Nparticle','Nhits vs Nparticles', 300, 0.,300, 300,0.,300)
htry=r.TH3D('htry','Digitized Hits ',75,-3000,3000 ,24,-600,600,24,-600,600)
htryNoCut=r.TH3D('htryNoCut','Digitized Hits ',75,-3000,3000 ,24,-600,600,24,-600,600)

hcounmu=r.TH1D('hcountmu','Number of muons hitting spectrometer from 0 track', 100,0,100)
hcountSBT=r.TH1D('hcountSBT','Number of events which are rejected as a background using the SBT', 20000,0,20000)
hcountUVT=r.TH1D('hcountUVT','Number of events which are rejected as a background using the UVT', 20000,0,20000)
hcountRPC=r.TH1D('hcountRPC','Number of events which are rejected as a background using the RPC', 20000,0,20000)
hcountSVT=r.TH1D('hcountSVT','Number of events which are rejected as a background using the SVT', 20000,0,20000)

hselectedHNL=r.TH1D('hselectedHNL', 'The number of events selected as HNL ', 5000,0.,5000)
intp_xy=r.TH2D('intp_xy','Interaction point ;X [cm];Y[cm]',800,-2000.,2000.,800,-2000,2000)

intp_xz=r.TH2D('intp_xz ','Interaction point ;X [cm];Z [cm]',800,-2000,2000,4000,-4000.,4000.)

deltaT_vs_deltaZ_muonIP=r.TH2D('deltaT_vs_deltaZ_muonIP',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)

id_DIS=r.TH1D('id_DIS ','The Id of the particle produced in the DIS ;Id ',5000,-5000,5000)

deltaT_vs_deltaZ_allparticles=r.TH2D('deltaT_vs_deltaZ_allparticles',';Distance to start of decay volume [m];distance from beam axis[m]',8000,0.,200.,800,0.,20.)

deltaT_vs_deltaZ_211=r.TH2D('deltaT_vs_deltaZ_211',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)
deltaT_vs_deltaZ_321=r.TH2D('deltaT_vs_deltaZ_321',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)

deltaT_vs_deltaZ_130=r.TH2D('deltaT_vs_deltaZ_130',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

deltaT_vs_deltaZ_310=r.TH2D('deltaT_vs_deltaZ_310',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

deltaT_vs_deltaZ_3122=r.TH2D('deltaT_vs_deltaZ_3122',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

deltaT_vs_deltaZ_3112=r.TH2D('deltaT_vs_deltaZ_3112',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

deltaT_vs_deltaZ_3312=r.TH2D('deltaT_vs_deltaZ_3312',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

deltaT_vs_deltaZ_13=r.TH2D('deltaT_vs_deltaZ_13',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)
deltaT_vs_deltaZ_candidate=r.TH2D('deltaT_vs_deltaZ_candidate',';Distance to start of decay volume [m];distance from beam axis[m]',4000,0.,100.,800,0.,20.)

hangleMuonYZ=r.TH1D('hangleMuonYZ','Angle beetwen the moemntum in Z',100,-1,1)
hangleMuonXZ=r.TH1D('hangleMuonXZ','Angle beetwen the moemntum in Z',100,-1,1)


hvtxrec=r.TH2D('hvtxrec ','The reconstructed vertex of all HNL candidates ;Z [cm];X[cm]',800,-2000.,2000.,12000,-30000,30000)
hdoca=r.TH1D('hdoca','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hip=r.TH1D('hip','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)


hmomentumCandidate=r.TH1D('hmomentumCandidate', 'The momentum of the HNL candidate;P[GeV]', 8000,0.,400.)
hHNLevent_2track=r.TH2D('hHNLevent_2track', 'The number of NHL candidates for event vs number of tracks per event',300,0.,300.,20,0.,20.) 

hipFidCut=r.TH1D('hipFidCut','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)

hdocaFidCut=r.TH1D('hdocaFidCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)


hipFidCut5cm=r.TH1D('hipFidCut5cm','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)


hipDocaCut=r.TH1D('hipDocaCut','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)


hdocaIPCut=r.TH1D('hdocaIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaIP250Cut=r.TH1D('hdocaIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)


hdocaUVTcut=r.TH1D('hdocaUVTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaSVTcut=r.TH1D('hdocaSVTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaSBTcut=r.TH1D('hdocaSBTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)


hmass_rec=r.TH1D('hmass_rec','The invarian mass of the HNL candidate ;M[GeV]; Nentries',100,0.,10)


mass_ip_Noveto=r.TH2D('mass_ip_Noveto',";mass[GeV/c];IP[cm]",100,0,10,500,0.,4000);


mass_ip=r.TH2D('mass_ip',";mass[GeV/c];IP[cm]",100,0,10,500,0.,4000);

vetoDets = {"UVT":(False,0),"SVT":(False,0),"SBT":(False,0),"RPC":(False,0),"TRA":(False,0)}

hdoca2cut=r.TH1D('hdoca2cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hIP5cm=r.TH1D('hIP5cm', 'IP of HNL candidate ;IP[cm];Entries ',200,0.,1000)
hdoca3cut=r.TH1D('hdoca3cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaIP10=r.TH1D('hdocaIP10', 'The doca of the HNL event;DOCA[cm];Nentries',200,0.,1000)
docaNo5cmCut=r.TH1D('docaNo5cmCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
ipNo5cmCut=r.TH1D('ipNo5cmCut', 'IP of HNL candidate ;IP[cm];Entries ',200,0.,1000)
hdisToWall=r.TH1D('hdisToWall','Distance to the closesest boundary', 1000, 0.,100) 
hdisToWall3=r.TH1D('hdisToWall3','Distance to the closesest boundary', 10000, 0.,1000) 

hPDGMotherHNLcandidate=r.TH1D('hPDGMotherHNLcandidate', "the pdg code of the mother of the HNL canidate in case when both dauthers come from the same particle",5000,-5000,5000)
hPDGMotherHNLDOCA1=r.TH1D('hPDGMotherHNLDOCA1',"the pdg code of the mother of the HNL candidate in case when both are the same and the doca beetwen the vertex is less than 1cm",5000,-5000,5000)
hmum1mum2=r.TH2D('hmum1mum2', "T2D histigram with the mother id of the tracks", 5000,-5000,5000,5000,-5000,5000)
htrackID=r.TH2D('htrackID', "T2D histigram with the mother id of the tracks", 500,0,500,500,0,500)



import shipVeto

hPDGMotherTrackcandidate=r.TH1D('hPDGMotherTrackcandidate',"the pdg code of the mother of the HNL candidate in case when both are the same and the doca beetwen the vertex is less than 1cm",5000,-5000,5000)
#construct geometry and load it 
deltaT_vs_deltaZ_130Spectrom=r.TH2D('deltaT_vs_deltaZ_130Spectrom',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)
deltaT_vs_deltaZ_310Spectrom=r.TH2D('deltaT_vs_deltaZ_310Spectrom',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)
pdgD130=r.TH1D('pdgD130',"the pdg code of the mother of the HNL candidate in case when both are the same and the doca beetwen the vertex is less than 1cm",5000,-5000,5000)
pdgD310=r.TH1D('pdgD310',"the pdg code of the mother of the HNL candidate in case when both are the same and the doca beetwen the vertex is less than 1cm",5000,-5000,5000)



deltaT_vs_deltaZ_allSpecto=r.TH2D('deltaT_vs_deltaZ_allSpecto',';Distance to start of decay volume [m];distance from beam axis[m]',8000,-100.,100.,800,0.,20.)

hcountTrack130=r.TH1I('hcountTrack130','Number of tracks in spectrometer from K',10,0.,10)
hcountTrack310=r.TH1I('hcountTrack310','Number of tracks in spectrometer from K',10,0.,10)
hcountTrackall=r.TH1I('hcountTrackall','Number of tracks in spectrometer from all',10,0.,10)
XYallstraw=r.TH2D('XYallstraw','XY 130',100,-200,200,100,-600,600)
XY130straw=r.TH2D('XY130straw','XY310',100,-200,200,100,-600,600)
XY310straw=r.TH2D('XY310straw','XY all',100,-200,200,100,-600,600)
hcountmuhit=r.TH1D('hcountmuhit', "number of muons hititin ", 10,0.,10)
hdistwall1=r.TH1D('hdistwall1', "distance o the wall", 1000, 0.,1000)
hdistwall2=r.TH1D('hdistwall2', "distance o the wall", 1000, 0.,1000)
hdis2Wall=r.TH1D('hdis2Wall', "distance o the wall", 10000, 0.,10000)
hNumberTrack=r.TH1D('hNumberTrack', "the pdg code of the mother of the HNL canidate in case when both dauthers come from the same particle",500,0,500)
pdgSBTHit=r.TH1D('pdgSBTHit', 'The pdg code of all particles hitting the SBT', 10000,-5000,5000)
ElossVetoev=r.TH1D('ElossVetoev','The energy loss of fired segmests for not vetoed events;"Eloss[MeV]',500,0 ,1000)
hmuonHitsSbt=r.TH2D('hmuonHitsSbt', 'the track id of the muon vs the number of hitted segments;Number of hitted segment;Nentries',200,0.,200,200,0.,200) 

xy_digi_hit_45Mcopy= r.TH2D('xy_digi_hit_45Mcopy',"The XY position of the  digitized hit(45MeV);X[cm];Y[cm]",24,-600,600,24,-600,600)
yz_digi_hit_45Mcopy=r.TH2D('yz_digi_hit_45Mcopy',"The YZ position of the  digitized hit(45MeV);Z[cm];Y[cm]",75,-3000,3000,24,-600,600)
xz_digi_hit_45Mcopy=r.TH2D('xz_digi_hit_45Mcopy',"The XZ position of the  digitized hit(45MeV);Z[cm];X[cm]",75,-3000,3000,24,-600,600)

Esegment=r.TH1D('Esegment', "; Eloss[MeV]",500,0 ,1000)
Esegment_45=r.TH1D('Esegment_45', "; Eloss[MeV]",500,0 ,1000)
hnumFiredSegment=r.TH1D('hnumFiredSegment', 'Number of fired segment for not veto events', 100,0.,100)
elossMuon=r.TH1D('elossMuon','The energy loss of all muons hitting the SBT ;"Eloss[MeV]',500,0 ,1000)
elossAllParticles=r.TH1D('elossAllParticles','The energy loss of all particles hitting the SBT ;"Eloss[MeV]',500,0 ,1000)

xy_digi_hit_45M= r.TH2D('xy_digi_hit_45M',"The XY position of the  digitized hit(45MeV);X[cm];Y[cm]",24,-600,600,24,-600,600)
yz_digi_hit_45M=r.TH2D('yz_digi_hit_45M',"The YZ position of the  digitized hit(45MeV);Z[cm];Y[cm]",75,-3000,3000,24,-600,600)
xz_digi_hit_45M=r.TH2D('xz_digi_hit_45M',"The XZ position of the  digitized hit(45MeV);Z[cm];X[cm]",75,-3000,3000,24,-600,600)
pdgSBThit=r.TH1I('pdgSBThit',"The pdg Code of each particle hitting the SBT", 10000, -5000,5000)

pdgSBTEcut=r.TH1I('pdgSBTEcut',"The pdg Code of each particle hitting the SBT", 10000, -5000,5000)
timeOFSegment=r.TH1D('timeOFSegment', "the time distribution of the fired segments", 1000, 0., 1000)
xy_sbthit_1peak= r.TH2D('xy_sbthit_1peak',"The XY position of the  digitized hit(45MeV);X[cm];Y[cm]",24,-600,600,24,-600,600)
xz_sbthit_1peak=r.TH2D('xz_sbthit_1peak',"The YZ position of the  digitized hit(45MeV);Z[cm];Y[cm]",24,-600,600,75,-3000,3000)
xy_sbthit_2peak= r.TH2D('xy_sbthit_2peak',"The XY position of the  digitized hit(45MeV);X[cm];Y[cm]",24,-600,600,24,-600,600)
xz_sbthit_2peak=r.TH2D('xz_sbthit_2peak',"The XZ position of the  digitized hit(45MeV);Z[cm];Y[cm]",24,-600,600,75,-3000,3000)

muonMomentum1peak=r.TH1D('muonMomentum1peak', 'The momentum of the HNL candidate;P[GeV]', 8000,0.,400.)

muonMomentum2peak=r.TH1D('muonMomentum2peak', 'The momentum of the HNL candidate;P[GeV]', 8000,0.,400.)
e_vs_pmuon=r.TH2D('e_vs_pmuon', "Eloss vs P" , 100, 0.,1, 300,0 ,600)
e_vs_pmuon1=r.TH2D('e_vs_pmuon1', "Eloss vs P" , 8000,0.,400., 300,0 ,600)

nDigiHit=r.TH1D('nDigiHit','Number of digitized hits without appying Energy Cut ;Digitized hits in the SBT', 350,0.,350)
nDigiHit45=r.TH1D('nDigiHit45','Number of digitized hits  appying 45MeV Energy Cut ;Digitized hits in the SBT', 350,0.,350)
nHitRatePerSegm=r.TH1D('nHitRatePerSegm', ';HitRate[kHz];Frequency' ,5000, 0.,1000)
nHitRatePerSegmCut=r.TH1D('nHitRatePerSegmCut', ';HitRate[kHz];Frequency' ,750, 0.,150)
hNintSimulated=r.TH1D('hNintSimulated','The number of simulated DIS interaction',5000,0.,5000)
n=r.TObjArray(300)
volume=[]
t = r.TTree( 'mytree', 'My Tree' )
v = r.std.vector( r.std.string )()
t._v = v
t.Branch( 'mystrings', v )
f = r.TFile(inputFile)
sTree = f.Get("cbmsim")
fgeo=r.TFile(geoFile)
fGeo=fgeo.FAIRGeom
ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py",Yheight = 10, tankDesign = 5 , muShieldDesign = 7, nuTauTargetDesign=1)
startDecayVol=ShipGeo.vetoStation.z+20.*u.cm
trackST1=ShipGeo.TrackStation1.z-20*u.cm

print "Start of the decay Vessel=", startDecayVol

	
#define all functions that you need for the analysis
from array import array
def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
  node = fGeo.FindNode(X,Y,Z)
  #print "The current node is=", node
  #print "the position=",X,Y,Z
  if ShipGeo.tankDesign < 5:
     if not 'cave' in node.GetName():return dist  # TP 
  else:
     if not 'decayVol' in node.GetName(): return dist
  start = array('d',[X,Y,Z])
  nsteps = 8
  dalpha = 2*r.TMath.Pi()/nsteps
  rsq = X**2+Y**2
  minDistance = 100 *u.m
  for n in range(nsteps):
    alpha = n * dalpha
    sdir  = array('d',[r.TMath.Sin(alpha),r.TMath.Cos(alpha),0.])
    node = fGeo.InitTrack(start, sdir)
    nxt = fGeo.FindNextBoundary()
    #print "The node is=", nxt
    if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0: 
    	return 0
	print "It is 0"
    distance = fGeo.GetStep()
    #print "I am estimating the  distance", distance
    #print "====================================================================="
    if distance < minDistance  : minDistance = distance
  return minDistance

def FidCut_Z(Z):
   if Z > trackST1 : return False
   if Z < startDecayVol: return False
   return True

def isInFiducial(Z,x,y,z):
   if Z > trackST1 : return False
   if Z < startDecayVol: return False
   if dist2InnerWall(x,y,z)<5*u.cm: return False
   return True 


def WeightsAnalysis(P):

        sigmaTot=0	
	fcross = r.TFile.Open('/afs/cern.ch/work/p/pvenkova/12collaborationMeeting/muDIScrossSec.root')
	for x in fcross.GetListOfKeys():
	       sigma = x.ReadObj()
	       if not sigma.GetTitle() == 'PythiaCross' : break  #'AnalyticCross'
	       sig=sigma.Eval(P)
	       #sigmaTot=sigmaTot+sig
	       #print "the moemntum=", p, "sig=",sig, "sigmaTot=", sigmaTot

	return sig








def ImpactParameter(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'): 
  	P = tMom.P()

  else:                 
  	P = tMom.Mag()
  for i in range(3):   
  	t += tMom(i)/P*(point(i)-tPos(i))
  dist = 0
  for i in range(3):  
  	dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = r.TMath.Sqrt(dist)
  return dist

def EstimateWeight(theta,phi):
	
	w=0
	if theta<0.06 and phi< -1.57:
		#w=4.9e-05
		w=0.157895
	if theta>0.06 and phi< -1.57:
		w=0
	if theta<0.06 and (-1.57<phi<0):
	 	w=8.476e-05
	if theta> 0.06 and (-1.57<phi<0):
		w=0.000158
	if theta <0.06 and (0<phi<1.57):
		w=1.6711e-05
	if theta >0.06 and (0<phi<1.57):
		w=0.000643
	if theta<0.06 and (1.57<phi<3.14):
		w=6.331e-05
	
	if theta>0.06 and (1.57<phi<3.14):
		w=0.000173
	return w


'''def match2MC(n, candidate):
    d1= candidate.GetDaughter1()
    d2=candidate.GetDaughter2()
    #print d1
    #print d2
    # the fitTrack2MC leave is an array of the same length as FitTracks.
    # it contains the index of the corresponding MCTrack.
    d1_mc = sTree.fitTrack2MC[int(d1)]
    d2_mc=  sTree.fitTrack2MC[int(d2)]
    mum1=sTree.MCTrack[d1_mc].GetMotherId()
    mum2=sTree.MCTrack[d2_mc].GetMotherId()
    pdgCode1= sTree.MCTrack[mum1].GetPdgCode()
    pdgCode2= sTree.MCTrack[mum2].GetPdgCode() 
    muW=sTree.MCTrack[1].GetWeight()
    xsec=sTree.MCTrack[0].GetWeight()
    rhoL=sTree.MCTrack[2].GetWeight()

    return pdgCode1, pdgCode2,xsec, muW, rhoL,mum1,mum2,d1_mc,d2_mc'''

# using Thomas file 1e-30
def prob2int(cross,weight,muWeight):
	prob2int=(cross*1e-27*6.022e+23*weight)
	return prob2int

def scaleFactor(muWeight):
	scale=(muWeight/1000)
	return scale
V0dict={130:'KL',310:'KS',3122:'Lambda',321:'K+',22:'gamma'}




veto=shipVeto.Task(sTree)
ev_dic={}
ev=0
ev_dic1=[]
ev_det={}
detDigi={}
numHitSeg=1
detDigiNocut={}
numHitSegNoCut=1
NintSimulated=0


for nb in range (0,sTree.GetEntries()):
    if nb!=9014:continue
    sTree.GetEntry(nb)
    print nb
    ev=ev+1
    trackSBT={}
    countTrack130=0
    countTrack310=0
    countTrackall=0
    ev=ev+1
    countmu=0
    trackStrawTubes=[]
    muWeight=sTree.MCTrack[1].GetWeight()
    weight=sTree.MCTrack[3].GetWeight()
    
    #cross=sTree.MCTrack[0].GetWeight()
    P=sTree.MCTrack[0].GetP()
    #cross=WeightsAnalysis(P)
    cross=sTree.MCTrack[0].GetWeight()
    scale=scaleFactor(muWeight)
    prob=prob2int(cross,weight,muWeight)

    '''print "============================================================================================================"

    #print "The event is=", 

    print "the momentum is =", P
    print "The probability is =", prob
    print "The weight is =", weight
    print "the cross section is =", cross
    print "the scale is=", scale
    print "the muWeight is=",muWeight
    
    print "prob*scale=", prob*scale'''
    

    NintSimulated=NintSimulated+prob*scale
    print "the prob*scale=", prob*scale
    print "the moemntum=",sTree.MCTrack[0].GetP()  
    print "N=", NintSimulated
    P=sTree.MCTrack[0].GetP()
    pz=sTree.MCTrack[0].GetPz() 
    px=sTree.MCTrack[0].GetPx()
    py=sTree.MCTrack[0].GetPy()
    theta=r.TMath.ACos(pz / P)
    phi=r.TMath.ATan2(py, px)
    #weightNoMagnet=EstimateWeight(theta, phi)
    intp_xy.Fill(sTree.MCTrack[1].GetStartX(),sTree.MCTrack[1].GetStartY(),prob*scale)
    intp_xz.Fill(sTree.MCTrack[1].GetStartX(),sTree.MCTrack[1].GetStartZ(),prob*scale)
    angleMuonZ=(pz/(px**2+pz**2+py**2))
    angleMuonY=(py/(px**2+pz**2+py**2))
    angleMuonX=(px/(px**2+pz**2+py**2))
    hangleMuonXZ.Fill(r.TMath.ATan2(angleMuonX,angleMuonZ))
    hangleMuonYZ.Fill(r.TMath.ATan2(angleMuonY,angleMuonZ))
    deltaZ= sTree.MCTrack[1].GetStartZ()-startDecayVol
    deltaT=math.sqrt(sTree.MCTrack[1].GetStartX()*sTree.MCTrack[1].GetStartX()+sTree.MCTrack[1].GetStartY()*sTree.MCTrack[1].GetStartY())
    deltaT_vs_deltaZ_muonIP.Fill(deltaZ/100,deltaT/100,prob*scale)
    #if P<20:continue 


	

    for track in sTree.MCTrack:
	
		id_DIS.Fill(track.GetPdgCode(),prob*scale)

		if track.GetMotherId()!=-1:
			pdgCode=sTree.MCTrack[track.GetMotherId()].GetPdgCode()
			deltaZ=track.GetStartZ()-startDecayVol
			deltaT=math.sqrt(track.GetStartX()*track.GetStartX()+track.GetStartY()*track.GetStartY())
			deltaT_vs_deltaZ_allparticles.Fill(deltaZ/100,deltaT/100,prob*scale)
				
			if pdgCode ==130:
					
				deltaT_vs_deltaZ_130.Fill(deltaZ/100,deltaT/100,prob*scale) 
				
						
				pdgD130.Fill( track.GetPdgCode())
									
			
			if pdgCode ==310:
				deltaT_vs_deltaZ_310.Fill(deltaZ/100,deltaT/100,prob*scale) 	
									
				pdgD310.Fill( track.GetPdgCode())

				
			if  pdgCode==211:
				deltaT_vs_deltaZ_211.Fill(deltaZ/100,deltaT/100,prob*scale)

			if pdgCode==321:
				deltaT_vs_deltaZ_321.Fill(deltaZ/100,deltaT/100,prob*scale)
			if pdgCode==3122:deltaT_vs_deltaZ_3122.Fill(deltaZ/100,deltaT/100,prob*scale)


    

    numHit45=0




    for aDigi in sTree.Digi_SBTHits:
    
        	
     	detIDD    = aDigi.GetDetectorID()
        if not detIDD>100000: continue  # not a LiSc or plastic detector
	if not detIDD<999999:continue
	node= aDigi.GetNode()
	
	if detIDD not in detDigi:

		detDigiNocut[detIDD]=numHitSegNoCut*prob*scale
	else:
		numHitSegNoCut=numHitSegNoCut*prob*scale
		detDigiNocut[detIDD]+=numHitSegNoCut



	xy_digi_hit_45M.Fill(aDigi.GetY(),aDigi.GetX(),prob*scale) 
        yz_digi_hit_45M.Fill(aDigi.GetZ(),aDigi.GetY(),prob*scale)
        xz_digi_hit_45M.Fill(aDigi.GetZ(),aDigi.GetX(),prob*scale)
	htryNoCut.Fill(aDigi.GetZ(),aDigi.GetY(),aDigi.GetX(),prob*scale)
	elos=aDigi.GetEloss()
	Esegment.Fill(elos*1000, prob*scale)
	if aDigi.isValid():
		if detIDD not in detDigi:
			detDigi[detIDD]=numHitSeg*prob*scale
		else:
			numHitSeg=numHitSeg*prob*scale
			detDigi[detIDD]+=numHitSeg


			

        	numHit45=numHit45+1 

		#htry.Fill(aDigi.GetZ(),aDigi.GetY(),aDigi.GetX(),prob*scale)
		#xy_digi_hit_45Mcopy.Fill(aDigi.GetY(),aDigi.GetX(),prob*scale) 		
                #yz_digi_hit_45Mcopy.Fill(aDigi.GetZ(),aDigi.GetY(),prob*scale)
                #xz_digi_hit_45Mcopy.Fill(aDigi.GetZ(),aDigi.GetX(),prob*scale)

		#Esegment_45.Fill(elos*1000,prob*scale)


    
    Nhits_vs_Nparticle.Fill(sTree.digiSBT2MC.size(),sTree.Digi_SBTHits.GetEntriesFast()) 
    nDigiHit45.Fill(numHit45,prob*scale)
    nDigiHit.Fill(sTree.Digi_SBTHits.GetEntriesFast(),prob*scale)

    #hNsegmentsDigi.Fill(len(det))	
    if numHit45>60:
    			
    	for aDigi in sTree.Digi_SBTHits:
    
        	
     		detIDD    = aDigi.GetDetectorID()
        	if not detIDD>100000: continue  # not a LiSc or plastic detector
		if not detIDD<999999:continue
		if aDigi.isValid():

			htry.Fill(aDigi.GetZ(),aDigi.GetY(),aDigi.GetX(),prob*scale)
			xy_digi_hit_45Mcopy.Fill(aDigi.GetY(),aDigi.GetX(),prob*scale) 		
                	yz_digi_hit_45Mcopy.Fill(aDigi.GetZ(),aDigi.GetY(),prob*scale)
                	xz_digi_hit_45Mcopy.Fill(aDigi.GetZ(),aDigi.GetX(),prob*scale)















    '''detID=0
    trackID=0
    detectors=[]
    for sbtHit in sTree.vetoPoint:
    	detID=sbtHit.GetDetectorID()
	trackSBTID=sbtHit.GetTrackID()
	p=math.sqrt(sbtHit.GetPx()**2+sbtHit.GetPy()**2+sbtHit.GetPz()**2)
	E=sbtHit.GetEnergyLoss()*1000
        if not detID>100000: continue  # not a LiSc or plastic detector
	if not detID<999999:continue

	if detID not in detectors:
		detectors.append(detID)
	pdgSBThit.Fill(sbtHit.PdgCode())
	elossAllParticles.Fill(sbtHit.GetEnergyLoss()*1000)
	if sbtHit.GetEnergyLoss()*1000>150:pdgSBTEcut.Fill(sbtHit.PdgCode())
	if abs(sbtHit.PdgCode())!=13: continue
	elossMuon.Fill(sbtHit.GetEnergyLoss()*1000)
	e_vs_pmuon.Fill(E/(p*1000),p)
	e_vs_pmuon1.Fill(p,E)
	if sbtHit.GetEnergyLoss()*1000>125 and  sbtHit.GetEnergyLoss()*1000<145:
		xy_sbthit_2peak.Fill(sbtHit.GetX(),sbtHit.GetY())
		xz_sbthit_2peak.Fill(sbtHit.GetX(),sbtHit.GetZ())
		muonMomentum2peak.Fill(p)
	if sbtHit.GetEnergyLoss()*1000>110 and  sbtHit.GetEnergyLoss()*1000<125:
		xy_sbthit_1peak.Fill(sbtHit.GetX(),sbtHit.GetY())
		xz_sbthit_1peak.Fill(sbtHit.GetX(),sbtHit.GetZ())
		muonMomentum1peak.Fill(p)
       
        if not trackSBTID>0:continue
        if not trackSBT.has_key(trackSBTID):
		trackSBT[trackSBTID]=[]
	
	if not  (trackSBTID,detID) in trackSBT.items():
		trackSBT[trackSBTID].append(detID)
    #for m in trackSBT:

   	#print "The lenght is =",len(trackSBT[m])
	#if len(trackSBT[m])>100: print "TWFKDHJ"

    hNsegments.Fill(len(detectors))

    	
    for n in trackSBT:
    	hmuonHitsSbt.Fill(n,len(trackSBT[n]))

    #print  "The number of all Digitized HIts=", sTree.Digi_SBTHits.GetEntriesFast()



    nHNLevent=sTree.Particles.GetEntries()
    #hnumTrack.Fill(len(sTree.goodTracks))
    numTrack=len(sTree.goodTracks)
    hHNLevent_2track.Fill(nHNLevent,len(sTree.goodTracks))
    if numTrack!=2:continue
    #print "Number of HNL candidate per event=", nHNLevent
    if sTree.Particles.GetEntries()!=1:continue
    for candidate in sTree.Particles:

    	#print "The probability is=", prob
	#print "mu weight is =", muWeight
	#print "cross=", cross
	#print "the scale is =", scale
	#print "The weight is =", weight
        indexCandidate=sTree.Particles.index(candidate)
	mum1,mum2,xsec,muW,rhoL,d1_mc,d2_mc,trackID1,trackID2=match2MC(sTree,candidate)
	
	if mum1==mum2:
		hPDGMotherHNLcandidate.Fill(mum1)
		hPDGMotherTrackcandidate.Fill(d1_mc)
		hNumberTrack.Fill(trackID1)

	else:
		hmum1mum2.Fill(mum1,mum2)
		htrackID.Fill(trackID1,trackID2)

		
 	vtx = r.TVector3()
	momentumHNL=r.TLorentzVector()
	candidate.GetVertex(vtx)
	candidate.GetMomentum(momentumHNL)
	mass_rec=momentumHNL.Mag2()
	hmass_rec.Fill(mass_rec,prob*scale)
	hvtxrec.Fill(vtx.Z(),vtx.X(),prob*scale)
    	deltaZ1= vtx.Z()-startDecayVol
    	deltaT1=math.sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y())
    	deltaT_vs_deltaZ_candidate.Fill(deltaZ1/100,deltaT1/100)
	




	doca=candidate.GetDoca()

	hdoca.Fill(doca,prob*scale)
	vtarget=r.TVector3(0,0,ShipGeo.target.z0)
	ip=ImpactParameter(vtarget,vtx,momentumHNL)
	hip.Fill(ip,prob*scale)
	hmomentumCandidate.Fill(momentumHNL.P(),prob*scale)
	distToWall,node,a=veto.fiducialCheckSignal(indexCandidate)
	volume.append(node)
	#hdisToWall.Fill(distToWall)

	minDis=dist2InnerWall(vtx.X(),vtx.Y(),vtx.Z())
	#print "The MIN DIS=",minDis
    	vetoDets['SBT'] = veto.SBT_decision()
    	vetoDets['SVT'] = veto.SVT_decision()
    	vetoDets['UVT'] = veto.UVT_decision()
    	vetoDets['RPC'] = veto.RPC_decision() 
    	if vetoDets['SBT'][2] != 0:countSBT+=1
    	if vetoDets['SVT'][2] != 0:countSVT+=1
    	if vetoDets['UVT'][2] != 0:countUVT+=1
    	if vetoDets['RPC'][2] != 0:countRPC+=1



	#First Cut on the Fiducial Volume

	#if isInFiducial(vtx.Z())==True: 
		#docaNo5cmCut.Fill(doca)
	#	ipNo5cmCut.Fill(ip)
	#	if doca<1:
	#		docaNo5cmCut.Fill(doca)
	#if distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True:
			#hipFidCut5cm.Fill(prob*scale)
			#hipFidCut.Fill(ip,prob*scale)
			#FidCutEv=FidCutEv+1
			#hdocaFidCut.Fill(doca,prob*scale)
			
	if doca<1:
		hipDocaCut.Fill(ip,prob*scale)
		DocaCutEv=DocaCutEv+1

	
	if ip<10:
		hdocaIPCut.Fill(doca,prob*scale)
		IPCutEv=IPCutEv+1

	if ip<250:
		hdocaIP250Cut.Fill(doca,prob*scale)

	if FidCut_Z(vtx.Z())==True: 
		hdisToWall3.Fill(a)
		hdis2Wall.Fill(minDis)





	if isInFiducial(vtx.Z(),vtx.X(),vtx.Y(),vtx.Z())==True:

		hipFidCut5cm.Fill(prob*scale)
                hipFidCut.Fill(ip,prob*scale)
                FidCutEv=FidCutEv+1
                hdocaFidCut.Fill(doca,prob*scale)
	#see how many events had been veto by the systems
	if vetoDets['UVT'][0]!=False:hdocaUVTcut.Fill(doca,prob*scale)
	if vetoDets['SBT'][0]!=True:
		hdocaSBTcut.Fill(doca,prob*scale)
		SBTCutEv=SBTCutEv+1
		mass_ip_Noveto.Fill(mass_rec,ip,prob*scale)
		#ElossVetoev.Fill(Eloss*1000)
		#hnumFiredSegment.Fill(hitSeg)



#	else:

        if vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=0 or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False:
		mass_ip.Fill(mass_rec,ip,prob*scale)
	
	if vetoDets['SVT'][0]!=False:hdocaSVTcut.Fill(doca,prob*scale)
	#if vetoDets['RPC'][0]!=False:hdocaRPCcut.Fill(doca)
	

	if not (isInFiducial(vtx.Z(),vtx.X(),vtx.Y(),vtx.Z())==True):continue
	hIP5cm.Fill(ip,prob*scale)
	if doca<1:hPDGMotherHNLDOCA1.Fill(mum1,mum2)
	if doca>1:continue
	hdoca2cut.Fill(doca,prob*scale)



#if ip>10: continue
	#hdocaIP10.Fill(ip,prob*scale)
	# for three body decay 
	if ip>250:continue
	hdoca3cut.Fill(doca,prob*muWeight*0.001)
	
	selectedWtihoutSBT=selectedWtihoutSBT+1

	#if vetoDets['SBT'][0]==True:continue
        if vetoDets['SBT'][0] or vetoDets['UVT'][0] or vetoDets['SVT'][0] or vetoDets['RPC'][0]:continue
	selectedEv=selectedEv+1
	hselectedHNL.Fill(selectedEv)
	print "Found HNL=", selectedEv'''
for n in detDigi:
	nHitRatePerSegmCut.Fill((detDigi[n])/1.2)
for m in detDigiNocut:
	nHitRatePerSegm.Fill((detDigiNocut[m])/1.2)
print "num=", NintSimulated
hNintSimulated.Fill(NintSimulated)
#print trackSBT.items()
hcountSBT.Fill(countSBT)
hcountUVT.Fill(countUVT)
hcountRPC.Fill(countRPC)
hcountSVT.Fill(countSVT)
outfile=r.TFile("muonDISstudy.root", "RECREATE")
for s in volume:
   v.push_back( s )
   t.Fill()
   t.Write()
mass_ip.Write()
hdocaSBTcut.Write()
hdocaUVTcut.Write()
mass_ip_Noveto.Write()
hdocaSVTcut.Write()
intp_xy.Write()
intp_xz.Write()
deltaT_vs_deltaZ_muonIP.Write()
id_DIS.Write()
deltaT_vs_deltaZ_allparticles.Write()
deltaT_vs_deltaZ_211.Write()
deltaT_vs_deltaZ_321.Write()
deltaT_vs_deltaZ_130.Write()
deltaT_vs_deltaZ_310.Write()
deltaT_vs_deltaZ_3122.Write()
deltaT_vs_deltaZ_3112.Write()
deltaT_vs_deltaZ_3312.Write()
hmass_rec.Write()
hvtxrec.Write()
hdoca.Write()
hip.Write()
hmomentumCandidate.Write()
hHNLevent_2track.Write()
hipFidCut.Write()
hipFidCut5cm.Write()
hdocaFidCut.Write()
hdocaIPCut.Write()
hipDocaCut.Write()
hdoca2cut.Write()
hIP5cm.Write()
hdocaIP250Cut.Write()
hselectedHNL.Write()
hcountSVT.Write()
hcountSBT.Write()
hcountRPC.Write()
hcountUVT.Write()
hdoca3cut.Write()
hPDGMotherHNLcandidate.Write()
hPDGMotherHNLDOCA1.Write()
hdocaIP10.Write()
hPDGMotherTrackcandidate.Write()
hmum1mum2.Write()
deltaT_vs_deltaZ_310Spectrom.Write()
deltaT_vs_deltaZ_130Spectrom.Write()
deltaT_vs_deltaZ_allSpecto.Write()
hcountTrack130.Write()
hcountTrackall.Write()
hcountTrack310.Write()
hcountTrack130.Write()
hcountTrack310.Write()
XYallstraw.Write()
XY130straw.Write()
XY130straw.Write()
hcountmuhit.Write()
deltaT_vs_deltaZ_candidate.Write()
docaNo5cmCut.Write()
ipNo5cmCut.Write()
hdistwall1.Write()
hdistwall2.Write()
hdisToWall.Write()
hdisToWall.Write()
hdis2Wall.Write()
xy_digi_hit_45M.Write()
yz_digi_hit_45M.Write()
xz_digi_hit_45M.Write()
Esegment_45.Write()
hNumberTrack.Write()
hdisToWall3.Write()
pdgSBTHit.Write()
htrackID.Write()
ElossVetoev.Write()
hmuonHitsSbt.Write()
xy_digi_hit_45Mcopy.Write()
yz_digi_hit_45Mcopy.Write()
xz_digi_hit_45Mcopy.Write()
Esegment.Write()
elossMuon.Write()
hnumFiredSegment.Write()
pdgSBThit.Write()
elossAllParticles.Write()
timeOFSegment.Write()
pdgSBTEcut.Write()
xy_sbthit_1peak.Write()
xz_sbthit_1peak.Write()
xy_sbthit_2peak.Write()
xz_sbthit_2peak.Write()
muonMomentum1peak.Write()
e_vs_pmuon.Write()
e_vs_pmuon1.Write()
pdgD130.Write()
pdgD310.Write()
htry.Write()
htryNoCut.Write()
Nhits_vs_Nparticle.Write()
muonMomentum2peak.Write()
nDigiHit.Write()
nDigiHit45.Write()
hangleMuonYZ.Write()
hNsegments.Write()
hNsegmentsDigi.Write()
hangleMuonXZ.Write()
nHitRatePerSegm.Write()
nHitRatePerSegmCut.Write()
hNintSimulated.Write()
