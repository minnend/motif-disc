package mdisc.VizTool;

import java.util.*;

import kdm.data.*;
import kdm.mlpr.suffix_tree.*;
import kdm.models.*;
import kdm.tools.SupTest;
import mdisc.VizTool.Tanaka.TanakaView;

import javax.swing.*;

/**
 * Container for "global" data used throughout the VizTool
 */
public class GlobalData
{
   public int nSymbols = 0;
   public String[] classes;
   public int[] nExamples;
   public TreeMap<String, ArrayList<Sequence>> labData;
   public TreeMap<String, ArrayList<DiscreteSeq>> labQData;
   public ArrayList<Sequence> tseries;
   public ArrayList<DiscreteSeq> qseries;
   public ArrayList<String> sentences;
   public SuffixTree sufTree;
   public LabeledView labeledView;
   public TreeView treeView;
   public HashView hashView;
   public AllDataView allDataView;
   public DataOverview dataOverview;
   public BoundaryView boundaryView;
   public WordSpottingView wordSpottingView;   
   public DiscoverView discoverView;
   public ChiuView chiuView;
   public TanakaView tanakaView;
   public SubDimView subdimView;   
   protected FeatureVec initv, minv, gmean, gvar, gmin, gmax;
   public int nTotalSeriesLength, nAvgSeriesLength;
   public int[] seqLens;
   public int vizgrp[][];   
   public JFrame frame;

   /** @return the i-th sentence or "(Unknown)" if it doesn't exist */
   public String getSentence(int i)
   {
      if (sentences==null || i>=sentences.size()) return "(Unknown)";
      else return sentences.get(i);
   }
   
   /** @return dimensionality of raw data */
   public final int getNumDims(){ return (tseries==null || tseries.isEmpty() ? 0 : tseries.get(0).getNumDims()); }
   
   /** @return number of classes */
   public final int getNumClasses(){ return classes==null ? 0 : classes.length; }
   
   /** @return number of sequences */
   public final int getNumSeqs(){ return tseries==null ? 0 : tseries.size(); }
   
   public FeatureVec getGlobalMean(){ return new FeatureVec(gmean); }
   public FeatureVec getGlobalMin(){ return new FeatureVec(gmin); }
   public FeatureVec getGlobalMax(){ return new FeatureVec(gmax); }
   public FeatureVec getGlobalVar(){ return new FeatureVec(gvar); }
   public FeatureVec getInitVar(){ return new FeatureVec(initv); }
   public FeatureVec getMinVar(){ return new FeatureVec(minv); }
   
   public void setGlobalMean(FeatureVec x){ gmean = new FeatureVec(x); }
   public void setGlobalMin(FeatureVec x){ gmin = new FeatureVec(x); }
   public void setGlobalMax(FeatureVec x){ gmax = new FeatureVec(x); }
   public void setGlobalVar(FeatureVec x){ gvar= new FeatureVec(x); }
   public void setInitVar(FeatureVec x){ initv = new FeatureVec(x); }
   public void setMinVar(FeatureVec x){ minv = new FeatureVec(x); }
   
   /**
    * use the labData field to initialize the class string array
    */
   public void updateClassStrings()
   {
      if (labData==null) return;
      classes = new String[labData.size()];
      nExamples = new int[labData.size()];
      Iterator<String> it = labData.keySet().iterator();
      for(int i=0; i<classes.length; i++)
      {
         classes[i] = it.next();
         nExamples[i] = labData.get(classes[i]).size();
      }
   }
   
   /** calculate global stats (min, max, mean, var, etc.) */
   public void calcGlobalStats()
   {
      Gaussian1D[] gmSeries = SupTest.calcGauss(tseries.toArray(new Sequence[0]));
      setInitVar(SupTest.initv);
      setMinVar(SupTest.minv);
      int nDims = getNumDims();
      FeatureVec gmean = new FeatureVec(nDims);
      FeatureVec gvar = new FeatureVec(nDims);
      for(int i = 0; i < nDims; i++){
         gmean.set(i, gmSeries[i].getMean());
         gvar.set(i, gmSeries[i].getVar());
      }
      setGlobalMean(gmean);
      setGlobalVar(gvar);
      FeatureVec gmin = tseries.get(0).getMin();
      FeatureVec gmax = tseries.get(0).getMax();
      int nSeqs = getNumSeqs();
      for(int i = 1; i < nSeqs; i++){
         gmin._min(tseries.get(i).getMin());
         gmax._max(tseries.get(i).getMax());
      }
      setGlobalMin(gmin);
      setGlobalMax(gmax);
   }
   
   /**
    * Calculate the average length of the given sequences and store it in avgSeriesLength
    * 
    * @param tseries sequences over which to calc average length
    */
   public void calcLengthStats(ArrayList<Sequence> tseries)
   {
      int nSeqs = tseries.size();
      seqLens = new int[nSeqs];
      nTotalSeriesLength = 0;
      for(int i=0; i<nSeqs; i++)
      {
         Sequence seq = tseries.get(i);
         seqLens[i] = seq.length();
         nTotalSeriesLength += seqLens[i];
      }
      nAvgSeriesLength = (int)Math.round((double)nTotalSeriesLength / tseries.size());      
   }

   /**
    * Setup the visual input grouping from a control string
    * 
    * @param s control string -- letters refer to dimensions (A-Z=1-26, a-z=27-52), commas
    *           separate different graphs
    * @return true if the visual group was setup properly from the control string
    */
   public boolean setupVizGrp(String s)
   {
      StringTokenizer st = new StringTokenizer(s, ", \r\n\t");
      int nViz = st.countTokens();
      vizgrp = new int[nViz][];
      for(int iViz = 0; iViz < nViz; iViz++)
      {
         String sViz = st.nextToken();
         vizgrp[iViz] = new int[sViz.length()];
         for(int i = 0; i < sViz.length(); i++)
         {
            char c = sViz.charAt(i);
            if (c >= 'A' && c <= 'Z') vizgrp[iViz][i] = (int)(c - 'A');
            else if (c >= 'a' && c <= 'z') vizgrp[iViz][i] = 26 + (int)(c - 'a');
            else{
               System.err.printf("Error: invalid visual grouping value (%c)\n", c);
               return false;
            }
         }
      }
      return true;
   }
}
