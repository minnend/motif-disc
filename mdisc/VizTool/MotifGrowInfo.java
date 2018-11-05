package mdisc.VizTool;

import java.util.*;

import kdm.data.*;
import kdm.models.*;
import kdm.util.*;

/**
 * Stores information pertinent to motif growing (extending a motif forward or backward in
 * time); this class also does the computation to determine if growing is a good idea and
 * can compare to possible growth opporunities to determine which is better.
 */
public class MotifGrowInfo
{
   protected int iMotif, maxFront, maxBack;
   protected ArrayList<WindowLocation> occs;
   protected ArrayList<Sequence> train;
   protected OatesModelUSamp oates;
   protected FeatureVec fvFront, fvBack;
   protected FeatureVec fvSDevFront, fvSDevBack;
   protected GlobalData gdata;

   public MotifGrowInfo(int _iMotif, ArrayList<WindowLocation> _occs, GlobalData _gdata, SpanList[] spanAvail)
   {
      iMotif = _iMotif;
      occs = _occs;
      maxFront = maxBack = 0;
      gdata = _gdata;
      int nDims = gdata.tseries.get(0).getNumDims();

      // build training set
      train = WindowLocation.getExamples(occs, gdata.tseries); 
         
      // calc max grow values
      maxFront = Integer.MAX_VALUE;
      maxBack = Integer.MAX_VALUE;
      for(WindowLocation wloc : occs)
      {
         Sequence seq = gdata.tseries.get(wloc.iSeries);
         SpanList span = spanAvail[wloc.iSeries];

         int ix = wloc.iStart - 1;
         while(ix >= 0 && span.contains(ix))
            ix--;
         int nFront = wloc.iStart - (ix + 1);
         if (nFront < maxFront) maxFront = nFront;

         ix = wloc.getLastIndex() + 1;
         while(ix < seq.length() && span.contains(ix))
            ix++;
         int nBack = (ix - 1) - wloc.getLastIndex();
         if (nBack < maxBack) maxBack = nBack;
      }

      // calc model variances
      GaussianDyn1D[] gmVar = null;
      if (maxFront > 0 || maxBack > 0)
      {
         // use oates model for alignment and mean/var estimation
         oates = new OatesModelUSamp(train, 3, gdata.initv, gdata.minv);

         // learn gaussian for variances per dim
         gmVar = new GaussianDyn1D[nDims];
         for(int d = 0; d < nDims; d++)
         {
            gmVar[d] = new GaussianDyn1D();
            gmVar[d].setReport(ProbFVModel.Report.prob);
            for(int i = 0; i < oates.length(); i++)
            {
               double x = oates.getModel(i, d).getVar();
               gmVar[d].add(x, false);
            }
            gmVar[d].update();
         }
      }

      if (maxFront > 0)
      {
         // Calc variance of prefix values
         GaussianDyn1D[] gmFront = new GaussianDyn1D[nDims];
         for(int d = 0; d < nDims; d++)
            gmFront[d] = new GaussianDyn1D();
         for(WindowLocation wloc : occs)
         {
            Sequence seq = gdata.tseries.get(wloc.iSeries);
            FeatureVec fv = seq.get(wloc.getFirstIndex() - 1);
            for(int d = 0; d < nDims; d++)
               gmFront[d].add(fv.get(d), false);
         }
         fvFront = new FeatureVec(nDims);
         fvSDevFront = new FeatureVec(nDims);
         for(int d = 0; d < nDims; d++)
         {
            gmFront[d].update();
            fvFront.set(d, gmFront[d].getVar());
            fvSDevFront.set(d, (gmFront[d].getVar()-gmVar[d].getMean())/gmVar[d].getSDev());
         }
      }

      if (maxBack > 0)
      {
         // Calc variance of suffix values
         GaussianDyn1D[] gmBack = new GaussianDyn1D[nDims];
         for(int d = 0; d < nDims; d++)
            gmBack[d] = new GaussianDyn1D();
         for(WindowLocation wloc : occs)
         {
            Sequence seq = gdata.tseries.get(wloc.iSeries);
            FeatureVec fv = seq.get(wloc.getLastIndex() + 1);
            for(int d = 0; d < nDims; d++)
               gmBack[d].add(fv.get(d), false);
         }
         fvBack = new FeatureVec(nDims);
         fvSDevBack = new FeatureVec(nDims);
         for(int d = 0; d < nDims; d++)
         {
            gmBack[d].update();
            fvBack.set(d, gmBack[d].getVar());
            fvSDevBack.set(d, (gmBack[d].getVar()-gmVar[d].getMean())/gmVar[d].getSDev());
         }
      }
   }

   /**
    * Compare vectors by projecting on to "1" and seeing which is longer (more positive)
    * @return -1 if a is a<b, 0 if a==b, or 1 if a>b; comparison based on which vector is "more positive"
    */
   public static int compareVecs(FeatureVec a, FeatureVec b)
   {
      int nd = a.getNumDims();
      assert (nd == b.getNumDims()) : String.format("feature vecs have different dimensionality (%d vs. %d)", nd, b.getNumDims());      
      FeatureVec fvOne = FeatureVec.ones(nd);
      double alen = a.projLen(fvOne);
      double blen = b.projLen(fvOne);
      if (alen < blen) return -1;
      if (alen > blen) return 1;
      return 0;
   }   

   public boolean canGrow()
   {
      return (canGrowFront() || canGrowBack());
   }

   public boolean canGrowFront()
   {
      return (maxFront > 0);
   }

   public boolean canGrowBack()
   {
      return (maxBack > 0);
   }
}
