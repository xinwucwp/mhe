import time
from utils import *

#pngDir = None
setupForSubset("f3d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 154,950,540; azimuth=240; elevation=25 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

pngDir = None
pngDir = "../../png/3d/f3d/"
gxfile = "gx"   #seismic volume
p2file = "p2"   #inline slope volume
p3file = "p3"   #crossline slope volume
epfile = "ep"   #planarity volume
mh1file = "mh1" #horizon one with local slopes and multi-grid correlations
mh2file = "mh2" #horizon two with local slopes and multi-grid correlations
sh1file = "sh1" #horizon one with local slopes
sh2file = "sh2" #horizon two with local slopes
plotOnly = True
plotOnly = False

def main(args):
  slopes()       
  goHorizonOne()
  goHorizonTwo()

#estimate local slopes
def slopes():
  gx = readImage(gxfile)
  if not plotOnly:
    sigma1,sigma2=8.0,2.0
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sigma1,sigma2,sigma2,5) 
    lsf.findSlopes(gx,p2,p3,ep);
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
  else:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
  plot3(gx,g=p2,cmin=-1, cmax=1,  cmap=jet,png="p2")
  plot3(gx,g=p3,cmin=-1, cmax=1,  cmap=jet,png="p3")
  plot3(gx,g=ep,cmin=0.1,cmax=0.8,cmap=jet,png="ep")

# pick horizon one with only one control point
def goHorizonOne():
  # one control point
  k1 = [90 ]
  k2 = [413]
  k3 = [312]
  gx  = readImage(gxfile)
  se = HorizonExtractor3()
  if not plotOnly:
    p2  = readImage(p2file)
    p3  = readImage(p3file)
    ep  = readImage(epfile)
    wp = pow(ep,8.0) 
    lmt = n1-1
    se.setWeights(0.0)
    se.setSmoothings(8.0,8.0)
    se.setCG(0.01,100)
    se.setExternalIterations(20)
    surf = se.surfaceInitialization(n1,n2,n3,k1,k2,k3)
    surf2 = copy(surf)
    start = time.time()
    #horizon extraction with slopes and multi-grid correlations
    surf1 = se.surfaceUpdateFromSlopesAndCorrelations(5,20,gx,wp,p2,p3,k2,k3,surf)
    print "Process time:"
    print (time.time()-start)
    start = time.time()
    #horizon extraction with only slopes
    surf2 = se.surfaceUpdateFromSlopes(wp,p2,p3,k2,k3,surf2)
    print "Process time:"
    print (time.time()-start)
    writeImage(mh1file,surf1) 
    writeImage(sh1file,surf2) 
  else:
    surf1 = readImage2(mh1file)
    surf2 = readImage2(sh1file)
  mp = ColorMap(-2.5,2.5,rwb)
  r1,g1,b1 = se.amplitudeRgb(mp,gx,surf1) 
  r2,g2,b2 = se.amplitudeRgb(mp,gx,surf2) 
  '''
  cmin = min(surf1)+5
  cmax = max(surf1)
  mp = ColorMap(cmin,cmax,jet)
  r1,g1,b1 = se.heightRgb(mp,surf1) 
  r2,g2,b2 = se.heightRgb(mp,surf2) 
  '''
  sf1 = [surf1,r1,g1,b1]
  sf2 = [surf2,r2,g2,b2]
  plot3(gx,hz=sf1,ks=[k1,k2,k3],cmap=rwb,png="surfm1")
  plot3(gx,hz=sf2,ks=[k1,k2,k3],cmap=rwb,png="surfs1")

# pick horizon two with two control point
def goHorizonTwo():
  # two control points
  k1 = [144,172]
  k2 = [667,340]
  k3 = [217,302]
  gx  = readImage(gxfile)
  se = HorizonExtractor3()
  if not plotOnly:
    p2  = readImage(p2file)
    p3  = readImage(p3file)
    ep  = readImage(epfile)
    wp = pow(ep,4.0) 
    lmt = n1-1
    se.setWeights(0.0)
    se.setSmoothings(8.0,8.0)
    se.setCG(0.01,100)
    se.setExternalIterations(20)
    surf = se.surfaceInitialization(n1,n2,n3,k1,k2,k3)
    surf2 = copy(surf)
    start = time.time()
    #horizon extraction with slopes and multi-grid correlations
    surf1 = se.surfaceUpdateFromSlopesAndCorrelations(5,20,gx,wp,p2,p3,k2,k3,surf)
    print "Process time:"
    print (time.time()-start)
    start = time.time()
    #horizon extraction with only slopes
    surf2 = se.surfaceUpdateFromSlopes(wp,p2,p3,k2,k3,surf2)
    print "Process time:"
    print (time.time()-start)
    writeImage(mh2file,surf1) 
    writeImage(sh2file,surf2) 
  else:
    surf1 = readImage2(mh2file)
    surf2 = readImage2(sh2file)
  mp = ColorMap(-2.5,2.5,rwb)
  r1,g1,b1 = se.amplitudeRgb(mp,gx,surf1) 
  r2,g2,b2 = se.amplitudeRgb(mp,gx,surf2) 
  '''
  cmin = min(surf1)+5
  cmax = max(surf1)
  mp = ColorMap(cmin,cmax,jet)
  r1,g1,b1 = se.heightRgb(mp,surf1) 
  r2,g2,b2 = se.heightRgb(mp,surf2) 
  '''
  sf1 = [surf1,r1,g1,b1]
  sf2 = [surf2,r2,g2,b2]
  plot3(gx,cmap=rwb,clab="Amplitude",png="seis")
  plot3(gx,hz=sf1,ks=[k1,k2,k3],cmap=rwb,png="surf2m")
  plot3(gx,hz=sf2,ks=[k1,k2,k3],cmap=rwb,png="surf2s")

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
bwr = ColorMap.BLUE_WHITE_RED
rwb = ColorMap.RED_WHITE_BLUE
def plot3(f,g=None,hz=None,ks=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-2.5,2.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.5,2.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(100)
  if hz:
    tg=TriangleGroup(True,s3,s2,hz[0],hz[1],hz[2],hz[3])
    sf.world.addChild(tg)
  if ks:
    pg = setPointGroup(ks[0],ks[1],ks[2],18)
    sf.world.addChild(pg)
  ipg.setSlices(n1-1,n2-1,0)
  if cbar:
    sf.setSize(967,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.08,0.00,0.05))
  ov.setTranslate(Vector3( 0.02,0.00,0.05))
  ov.setAzimuthAndElevation(130.0,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")
 

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def setPointGroup(k1,k2,k3,size):
  np  = len(k1)
  xyz = zerofloat(np*3)
  rgb = zerofloat(np*3)
  ki = 0
  for i in range(np):
    xyz[ki  ] = k3[i]
    xyz[ki+1] = k2[i]
    xyz[ki+2] = k1[i]
    rgb[ki  ]  = 0#1/225 
    rgb[ki+1]  = 1#225/225 
    rgb[ki+2]  = 0#1/225 
    ki = ki+3
  pg = PointGroup(size,xyz,rgb);
  states = StateSet();
  cs = ColorState();
  cs.setColor(Color.GREEN);
  states.add(cs);
  lms = LightModelState();
  lms.setTwoSide(True);
  states.add(lms);
  ms = MaterialState();
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
  ms.setShininess(100.0);
  states.add(ms);
  pg.setStates(states);
  return pg;

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
