package mhe;

import java.awt.*;
import java.util.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Extracts a seismic horizon surface with both local slopes 
 * and multi-grid correlations.
 * * <p>
 * We propose a novel method to automatically extract seismic horizons 
 * that follow consistent phase (e.g., peaks or troughs) and correctly 
 * track reflections across faults. In this method.
 * <p>
 * We estimate local reflection slopes by using structure tensors.
 * <p>
 * We compute multi-grid correlations/slopes by using dynamic time warping 
 * to correlate seismic traces at multiple laterally coarse grids. 
 * These coarse-grid correlations can correctly correlate reflections that 
 * may be significantly dislocated by faults or other dis- continuous structures. 
 * <p>
 * We compute a horizon by fitting, in the least-squares sense, the slopes 
 * of the horizon with both the precomputed local reflection slopes and 
 * multi-grid slopes. 
 * In this least-squares system, the local slopes on the fine grid and 
 * the multiple coarse-grid slopes will fit a consistent horizon in areas 
 * without lateral discontinuities.
 * Across laterally discontinuous areas where the local slopes fail to 
 * correctly correlate reflections and mislead the horizon extraction, the 
 * coarse-grid slopes will help to find the corresponding reflections and 
 * correct the horizon extraction. 
 * <p>
 * In addition, the multigrid correlations or slopes computed by dynamic 
 * warping can also assist in computing phase-consistent horizons.
 * <p>
 * Control points are incorporated in a preconditioner in the conjugate 
 * gradient method used to solve the linear system for horizon extracting.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class HorizonExtractor3 {
  /**
   * Sets half-widths for smoothings in 1st and 2nd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   */
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  /**
   * Sets outer interations for the horizon extraction.
   * @param exniter maximum number of outer iterations.
   */
  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }

  /**
   * Sets the relative weight of the curvature term for smoothing.
   * the default weight is set to be 0
   * @param w1 the weight.
   */
  public void setWeights(float w){
    _weight = w;
  }


  /**
   * Interpolates an initial horizon surface passing through control points.
   * One control point will initialize a horizontal surface
   * @param n1 1st dimension of the seismic volume.
   * @param n2 2nd dimension of the seismic volume.
   * @param n3 3rd dimension of the seismic volume.
   * @param k1 array of 1st coordinates of control points.
   * @param k2 array of 2nd coordinates of control points.
   * @param k3 array of 3rd coordinates of control points.
   */
  public float[][] surfaceInitialization(
    int n1, int n2, int n3, float[] k1, float[] k2, float[] k3) {
    //if only one control point, initialize a horizontal surface
    if (k1.length==1) {
      float[][] surf = zerofloat(n2,n3);
      add(surf,k1[0],surf);
      return surf; 
    } else {
      int nc = k1.length;
      MedianFinder mf = new MedianFinder(nc);
      float m1 = mf.findMedian(k1);
      float[] d1 = new float[nc];
      for (int ic=0; ic<nc; ++ic) {
        d1[ic] = k1[ic]-m1;
      }
      Sampling s2 = new Sampling(n2,1.0f,0.0f);
      Sampling s3 = new Sampling(n3,1.0f,0.0f);
      RadialInterpolator2.Biharmonic bs = new RadialInterpolator2.Biharmonic();
      RadialGridder2 rg = new RadialGridder2(bs,d1,k2,k3);    
      float[][] surf = add(m1,rg.grid(s2,s3));
      surfaceCorrect(surf,n1-1);
      checkControlPoints(k2, k3, surf); 
      return surf;
    }
  }

  /**
   * Extracts a horizon surface by using only local reflection slopes.
   * @param ep 3D array of planarity used to weight the horizon extraction.
   * @param p2 3D array of inline slopes.
   * @param p3 3D array of crossline slopes.
   * @param k2 1D array of 2nd coordinates of control points.
   * @param k3 1D array of 3rd coordinates of control points.
   * @param hs 2D array of an initial horizon surface (z coordinates) to be updated.
   */
  public float[][] surfaceUpdateFromSlopes(
    float[][][] ep, float[][][] p2 ,float[][][] p3,
    float[] k2, float[] k3, float[][] hs)
  {	
    int n3 = p2.length; 
    int n2 = p2[0].length; 
    int n1 = p2[0][0].length; 
    float[][] hst = copy(hs);
    float[][] b   = new float[n3][n2]; 
    float[][] p21 = new float[n3][n2]; 
    float[][] p31 = new float[n3][n2]; 
    float[][] ep1 = new float[n3][n2]; 
    checkControlPoints(k2, k3, hs); 
    VecArrayFloat2 vb    = new VecArrayFloat2(b);
    VecArrayFloat2 vh = new VecArrayFloat2(hs);
    int niter = 100;
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      updateSlopesAndWeights(p2,p3,ep,hs,p21,p31,ep1);
      A2 a2 = new A2(_weight,ep1,null);
      M2 m2 = new M2(_sigma1,_sigma2,ep1,k2,k3);
      if(n>5) {niter=_niter;}
      CgSolver cs = new CgSolver(_small,niter);
      vb.zero();
      makeRhs(ep1,p21,p31,b);
      cs.solve(a2,m2,vb,vh);
      checkControlPoints(k2, k3, hs); 
      hs = vh.getArray();
      surfaceCorrect(hs,n1-1);
      float ad = sum(abs(sub(hst,hs)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      //if (ad<0.02f) break;
      hst = copy(hs);
    }
    return vh.getArray();
  }


  /**
   * Extracts a horizon surface by using both local slopes 
   * and multi-grid correlations.
   * We compute pair-wise coolreation of seismic traces on a coarse 
   * grid by using dynamic time warping.
   * To improve the efficiency of the pair-wise correlations, we do 
   * not vertically correlate the traces in the entire sampling range 
   * of depth or time. Instead, we compute correlations of the traces 
   * only in a small depth window centered at the horizon.
   * @param um coarse-grid interval: compute correlation for every um traces.
   * @param wd defines the vertical window to be correlated.
   * @param ep 3D array of planarity used to weight the horizon extraction.
   * @param p2 3D array of inline slopes.
   * @param p3 3D array of crossline slopes.
   * @param k2 1D array of 2nd coordinates of control points.
   * @param k3 1D array of 3rd coordinates of control points.
   * @param hs 2D array of an initial horizon surface (z coordinates) to be updated.
   */
  public float[][] surfaceUpdateFromSlopesAndCorrelations(
    int um, int wd, float[][][] fx,
    float[][][] ep, float[][][] p2 ,float[][][] p3, 
    float[] k2, float[] k3, float[][] hs)
  {	
    int n3 = p2.length; 
    int n2 = p2[0].length; 
    int n1 = p2[0][0].length; 
    float[][] hst = copy(hs);
    float[][] b   = new float[n3][n2]; 
    float[][] p21 = new float[n3][n2]; 
    float[][] p31 = new float[n3][n2]; 
    float[][] ep1 = new float[n3][n2]; 
    checkControlPoints(k2, k3, hs); 
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vh = new VecArrayFloat2(hs);
    //GlobalCorrelationFinder gcf = new GlobalCorrelationFinder(-10,10);
    int dm = 10;
    int dc = 10;
    float[] u1 = null;
    float[][] pc = null;
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      dm = min(dm,100);
      if(n>1) {
        pc = getTraceIndexes(10,10,dm,dc,k2,k3,n2,n3,0.01f);
        int ns = pc.length;
        u1 = new float[ns];
      }
      updateSlopesAndWeights(um,wd,fx,p2,p3,ep,pc,hs,p21,p31,ep1,u1);
      A2 a2 = new A2(_weight,ep1,pc);
      M2 m2 = new M2(_sigma1,_sigma2,ep1,k2,k3);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(ep1,p21,p31,pc,u1,b);
      cs.solve(a2,m2,vb,vh);
      checkControlPoints(k2, k3, hs); 
      hs = vh.getArray();
      surfaceCorrect(hs,n1-1);
      float ad = sum(abs(sub(hst,hs)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      hst = copy(hs);
      dm += 20;
      dc += 20;
    }
    return vh.getArray();
  }


  //Colors the horizon surface by depth
  //For displaying only
  public float[][][] heightRgb(ColorMap mp, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = sf[i3][i2];
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  //Colors the horizon surface by amplitudes
  //For displaying only
  public float[][][] amplitudeRgb(ColorMap mp, float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  //Colors the horizon surface by amplitudes
  //For displaying only
  public float[][][] amplitudeRgb(
    int d1, ColorMap mp, float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] sa = new float[n3*n2];
    int k = -1;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) { 
    for (int i2=0; i2<n2; ++i2) {
      k++;
      for (int k1=-d1; k1<=d1; k1++) {
        sa[k] += si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]+k1);
      }
    }}
    sa = div(sa,d1*2f+1f);
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }


  ///////////////////////////////////////////////////////////////////////////
  //private funcitons


  //update local slopes and weights on the current horizon surface
  private static void updateSlopesAndWeights (
    float[][][] p, float[][][] q, float[][][] ep,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        float x1i = surf[i3][i2];
	      float wi = si.interpolate(n1,1.0,0.0,ep[i3][i2],x1i);
        wi1[i3][i2] = (wi>0.0005f)?wi:0.0005f;
        pi1[i3][i2] = si.interpolate(n1,1.0,0.0,p[i3][i2],x1i);
	      qi1[i3][i2] = si.interpolate(n1,1.0,0.0,q[i3][i2],x1i);
      }
    }
  }

  //update local slopes, multi-grid slops, and weights on 
  //the current horizon surface
  private static void updateSlopesAndWeights (
    int um, int dc, float[][][] fx,
    float[][][] p, float[][][] q, float[][][] ep, float[][] pc,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1, float[] ui1)
  {
    final int n3 = p.length;
    final int n2 = p[0].length;
    final int n1 = p[0][0].length;
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++){
        float x1i = surf[i3][i2];
	      float wi = si.interpolate(n1,1.0,0.0,ep[i3][i2],x1i);
        wi1[i3][i2] = (wi>0.0005f)?wi:0.0005f;
        pi1[i3][i2] = si.interpolate(n1,1.0,0.0,p[i3][i2],x1i);
	      qi1[i3][i2] = si.interpolate(n1,1.0,0.0,q[i3][i2],x1i);
      }
    }});
    if(pc!=null) updateCorrelations(um,dc,fx,pc,surf,ui1);
  }

  //updates multi-grid correlations of seismic traces within a vertical window
  //centered at the current horizon surface.
  private static void updateCorrelations(
    int um, int dc, float[][][] fx, float[][] pc, float[][] surf, float[] ui1) {
    final int ns = ui1.length;
    int n1 = fx[0][0].length;
    final SincInterpolator usi = new SincInterpolator();
    final DynamicWarping dw = new DynamicWarping(-um,um);
    dw.setStrainMax(0.25);
    dw.setErrorSmoothing(2);
    dw.setShiftSmoothing(4);
    final Random rd = new Random();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      int p2 = (int)pc[is][0];
      int p3 = (int)pc[is][1];
      int m2 = (int)pc[is][2];
      int m3 = (int)pc[is][3];
      int cp = round(surf[p3][p2]);
      int cm = round(surf[m3][m2]);
      int nc = dc*2+1;
      float[] fp = new float[nc];
      float[] fm = new float[nc];
      for (int i1 = cp-dc, p1=0; i1<=cp+dc; i1++, p1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fp[p1] = fx[p3][p2][k1];
        } else {
          fp[p1] = fx[p3][p2][i1];
        }
      }
      for (int i1 = cm-dc, m1=0; i1<=cm+dc; i1++, m1++) {
        if(i1<0||i1>=n1) {
          int k1 = rd.nextInt(n1);
          fm[m1] = fx[m3][m2][k1];
        } else {
          fm[m1] = fx[m3][m2][i1];
        }
      }
      float[] ui = dw.findShifts(fp,fm);
      int dm = min(cm,abs(n1-cm));
      int dp = min(cp,abs(n1-cp));
      if(min(dm,dp)<2) pc[is][4] = 0f;
      ui1[is] = usi.interpolate(nc,1.0,0.0,ui,surf[m3][m2]-cm+dc)+cm-cp;
    }});
  }



    /**
   * Finds 3D indexes of traces to be correlated.
   */
  private float[][] getTraceIndexes(
    int d2, int d3, int dm, int dc, 
    float[] k2, float[] k3, int n2, int n3, float scale) { 
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(dm);
    float[][] wm = new float[n3*2+1][n2*2+1]; 
    wm[n3][n2] = 1f;
    rgf.apply00(wm,wm);
    wm = div(wm,max(wm));
    ArrayList<int[]> ks = new ArrayList<int[]>();
    for (int i3=0; i3<n3; i3+=d3) {
    for (int i2=0; i2<n2; i2+=d2) {
      ks.add(new int[]{i2,i3});
    }}
    int ns = ks.size();
    ArrayList<float[]> kc = new ArrayList<float[]>();
    for (int is=0; is<ns; ++is) {
    for (int js=is+1; js<ns; ++js) {
      int p2 = ks.get(is)[0];
      int p3 = ks.get(is)[1];
      int m2 = ks.get(js)[0];
      int m3 = ks.get(js)[1];
      int b2 = abs(p2-m2);
      int b3 = abs(p3-m3);
      float ds = sqrt(b2*b2+b3*b3);
      if(ds<=dm)
        kc.add(new float[]{p2,p3,m2,m3,scale*wm[n3+b3][n2+b2]});
    }}

    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(1);
    RecursiveGaussianFilter rgfc = new RecursiveGaussianFilter(dc);
    float[][] wc = new float[n3*2+1][n2*2+1]; 
    wc[n3  ][n2  ] = 10f;
    rgf1.apply00(wc,wc);
    rgfc.apply00(wc,wc);
    wc = sub(wc,min(wc));
    wc = div(wc,max(wc));
    int np = k2.length;
    for (int ip=0; ip<np; ++ip) {
      int k2i = (int)k2[ip];
      int k3i = (int)k3[ip];
      int b2 = max(k2i-dc,0);
      int e2 = min(k2i+dc,n2-1);
      int b3 = max(k3i-dc,0);
      int e3 = min(k3i+dc,n3-1);
      for (int i3=b3; i3<=e3; i3+=d3) {
      for (int i2=b2; i2<=e2; i2+=d2) {
        int s2 = abs(i2-k2i);
        int s3 = abs(i3-k3i);
        kc.add(new float[]{k2i,k3i,i2,i3,5*scale*wc[n3+s3][n2+s2]});
      }}
    }
    return kc.toArray(new float[0][]);
  }



  //make sure the surface does not go beyond the top and bottom bounds
  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (surf[i2][i1]<0.f) surf[i2][i1]=0.f;
        if (surf[i2][i1]>lmt) surf[i2][i1]=lmt;
      }
    }
  }

  private static class A2 implements CgSolver.A{
    A2(float w1, float[][] wp, float[][] pc){
      _w1 = w1;
      _wp = wp;
      _pc = pc;
      //testSpd();
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if(_pc!=null) screenLhs(_pc,_wp,x,y);
      if (_w1>0.0f) {
        applyLhs(_wp,y,yy);
        v2y.add(1.f,v2yy,_w1);
      }
    }
    private float _w1;
    private float[][] _wp;
    private float[][] _pc;
    public void testSpd() {
    // symmetric: y'Ax = x'(A'y) = x'Ay
    // positive-semidefinite: x'Ax >= 0
      int n2 = _wp.length;
      int n1 = _wp[0].length;
      float[][] x = sub(randfloat(n1,n2),0.5f);
      float[][] y = sub(randfloat(n1,n2),0.5f);
      float[][] ax = zerofloat(n1,n2);
      float[][] ay = zerofloat(n1,n2);
      VecArrayFloat2 vx = new VecArrayFloat2(x);
      VecArrayFloat2 vy = new VecArrayFloat2(y);
      VecArrayFloat2 vax = new VecArrayFloat2(ax);
      VecArrayFloat2 vay = new VecArrayFloat2(ay);
      apply(vx,vax);
      apply(vy,vay);
      double yax = vy.dot(vax);
      double xay = vx.dot(vay);
      double xax = vx.dot(vax);
      double yay = vy.dot(vay);
      System.out.println("A3: yax="+yax+" xay="+xay);
      System.out.println("A3: xax="+xax+" yay="+yay);
    }
  }

  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, float[] k2, float[] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k2,_k3,y);
    }
    private float _sigma1,_sigma2;
    private float[][] _wp;
    private float[] _k2,_k3;
  }

  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
  }

  private static void checkControlPoints(
    float[] k1, float[] k2, float[] k3, float[][] f) 
  {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
        f[i3][i2] = k1[ip];
      }
    }
  }

  //Original constrain operation
  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

  //Modified constrain operation suggested by Dr. Yong Ma at Conocophillips
  //This modification is helpful when the updated horizon is vertically 
  //shifted far away from the control points
  private static void constrainM(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      float zshift = 0f;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        zshift -= x[i3][i2];
      }
      zshift /= np;
      int n3 = x.length;
      int n2 = x[0].length;
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        x[i3][i2] += zshift;
      }}
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }


  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    final LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[] xt = zerofloat(n1);
      float[] yt = zerofloat(n1);
      float[] st = zerofloat(n1);
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }});
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    final int n2 = x.length;
    final int n1 = x[0].length;
    final float c = 0.5f*sigma*sigma;
    final LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[] xt = zerofloat(n2);
      float[] yt = zerofloat(n2);
      float[] st = zerofloat(n2);
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }});
  }

  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    final int n2 = x.length;
    Parallel.loop(1,n2,2,new Parallel.LoopInt() { // i2 = 1, 3, 5, ...
    public void compute(int i2) {
      applyLhsSlice2(i2,wp,x,y);
    }});
    Parallel.loop(2,n2,2,new Parallel.LoopInt() { // i2 = 2, 4, 6, ...
    public void compute(int i2) {
      applyLhsSlice2(i2,wp,x,y);
    }});
  }

  private static void applyLhsSlice2(
    int i2, float[][] wp, float[][] x, float[][] y) {
    int n1 = x[0].length;
    for (int i1=1; i1<n1; ++i1) {
      float wpi = (wp!=null)?wp[i2][i1]:1.000f;
      if(wpi<0.05f) {wpi=0.05f;}
      float ws = wpi*wpi*0.25f;
      float xa = 0.0f;
      float xb = 0.0f;
      xa += x[i2  ][i1  ];
      xb -= x[i2  ][i1-1];
      xb += x[i2-1][i1  ];
      xa -= x[i2-1][i1-1];
      float x1 = xa+xb;
      float x2 = xa-xb;
      float y1 = ws*x1;
      float y2 = ws*x2;
      float ya = y1+y2;
      float yb = y1-y2;
      y[i2  ][i1  ] += ya;
      y[i2  ][i1-1] -= yb;
      y[i2-1][i1  ] += yb;
      y[i2-1][i1-1] -= ya;
    }
  }

  private static void screenLhs(
    float [][] ks, float[][] w, float[][] x, float[][] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p2 = (int)ks[is][0];
      int p3 = (int)ks[is][1];
      int m2 = (int)ks[is][2];
      int m3 = (int)ks[is][3];
      float scale = ks[is][4];
      float wi = (w[p3][p2]+w[m3][m2])*0.5f;
      float ws = wi*wi; 
      float dx = 0.0f;
      dx += x[p3][p2];
      dx -= x[m3][m2];
      dx *= scale*ws;
      y[m3][m2] -= dx;
      y[p3][p2] += dx;
    }
  }

  private static void makeRhs(
    float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    zero(y);
    int n2 = y.length;
    int n1 = y[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float ws = wpi*wpi*0.5f;
        float y1 = ws*p2i;
        float y2 = ws*p3i;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void makeRhs(
    float[][] wp, float[][] p2, float[][] p3, 
    float[][] pc, float[] us, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float wps = wpi*wpi*0.5f;
        float y1 = wps*p2i;
        float y2 = wps*p3i;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
    if(pc!=null) screenRhs(wp,pc,us,y);
  }


  private static void screenRhs(
    float[][] w, float[][] ks, float[] us, float[][] y) {
    int ns = ks.length;
    for (int is=0; is<ns; ++is) {
      int p2 = (int)ks[is][0];
      int p3 = (int)ks[is][1];
      int m2 = (int)ks[is][2];
      int m3 = (int)ks[is][3];
      float scale = ks[is][4];
      float wi = (w[p3][p2]+w[m3][m2])*0.5f;
      float ws = wi*wi; 
      float ui = -us[is];
      ui *= scale*ws;
      y[m3][m2] -= ui;
      y[p3][p2] += ui;
    }
  }


  private void checkConstraintForce(float[] k2, float[] k3, float[] cf) {
    int np = k2.length;
    div(cf,max(cf),cf);
    for (int ip=0; ip<np; ++ip) {
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      System.out.println(" i2="+i2+" i3="+i3+" ad="+cf[ip]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.0f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
