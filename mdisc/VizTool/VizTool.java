package mdisc.VizTool;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.io.*;
import kdm.util.*;
import kdm.gui.*;
import kdm.models.*;
import kdm.mlpr.suffix_tree.*;
import kdm.tools.*;
import mdisc.VizTool.Tanaka.TanakaView;

import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.border.*;
import gnu.getopt.*;

/**
 * Tool for exploring time series data with different views and operations.
 */
public class VizTool
{
   protected GlobalData gdata;

   public static void usage()
   {
      System.err.println();
      System.err.println("Usage: java kdm.tools.VizTool.VizTool [options] <data def file>");
      System.err.println(" Options: ");
      System.err.println("  -help/?                 Usage information");
      System.err.println("  -quant <def file>       Specify a quantized data def file");
      System.err
            .println("  -vizgrp <info>          Specify a visual grouping -- info: letters separated by commas");
      System.err.println("                           e.g. \"AC,BD\" to group dims 1+3 and then 2+4");
      System.err.println();
      System.exit(1);
   }

   /**
    * Report any NaN or inf values in the data set
    */
   protected static void checkData(ArrayList<Sequence> tseries)
   {
      for(int i = 0; i < tseries.size(); i++){
         Sequence seq = tseries.get(i);
         int T = seq.length();
         int D = seq.getNumDims();
         for(int t = 0; t < T; t++){
            FeatureVec fv = seq.get(t);
            for(int d = 0; d < D; d++){
               double x = fv.get(d);
               if (Double.isNaN(x) || Double.isInfinite(x))
                  System.err.printf("Error: invalid data -- seq=%s (%d)  t=%d  d=%d  v=%f\n", seq.getName(),
                        i, t, d, x);
            }
         }
      }
   }

   /**
    * Load the data and setup the various views
    * 
    * @param sDefFile data def file for continuous data
    * @param sQuantDefFile data def file for quantized data
    * @return true if successful
    */
   public boolean setup(String sDefFile, String sQuantDefFile, String sVizGrp)
   {
      TimerMS timer = new TimerMS();
      gdata = new GlobalData();
      if (sVizGrp != null && !gdata.setupVizGrp(sVizGrp)){
         System.err.println("VizTool) Failed to setup visualization groups");
         return false;
      }

      // load data
      TreeMap<String, ArrayList<Sequence>> data = LabeledDataLoader.load(new File(sDefFile));
      if (data == null) System.err.printf("VizTool) Failed to load labeled data\n (%s)\n", sDefFile);

      gdata.labData = data;
      gdata.updateClassStrings();
      ArrayList<Sequence> tseries = LabeledDataLoader.tseries;
      gdata.calcLengthStats(tseries);
      System.err.printf("Loaded %d time series (len: total=%d, avg=%d).\n", tseries.size(),
            gdata.nTotalSeriesLength, gdata.nAvgSeriesLength);
      checkData(tseries);
      gdata.tseries = tseries;
      gdata.sentences = LabeledDataLoader.sentences;
      LabeledDataLoader.dumpDataSummary(data);

      Gaussian1D[] gmSeries = SupTest.calcGauss(tseries.toArray(new Sequence[0]));
      gdata.initv = SupTest.initv;
      gdata.minv = SupTest.minv;
      int nDims = gdata.getNumDims();
      gdata.gmean = new FeatureVec(nDims);
      gdata.gvar = new FeatureVec(nDims);
      for(int i=0; i<nDims; i++){
         gdata.gmean.set(i, gmSeries[i].getMean());
         gdata.gvar.set(i, gmSeries[i].getVar());
      }

      // load quantized data
      TreeMap<String, ArrayList<DiscreteSeq>> qdata = null;
      ArrayList<DiscreteSeq> qseries = null;
      if (sQuantDefFile != null){
         qdata = LabeledDataLoader.loadDiscrete(new File(sQuantDefFile), data);
         qseries = LabeledDataLoader.qseries;

         if (qseries != null){
            // figure out how many symbols there are
            timer.reset();
            HashSet<Integer> uniq = new HashSet<Integer>();
            int vmax = -1;
            for(DiscreteSeq dseq : qseries){
               int T = dseq.length();
               for(int t = 0; t < T; t++){
                  int x = dseq.geti(t);
                  uniq.add(x);
                  if (x > vmax) vmax = x;
               }
            }
            if (vmax + 1 != uniq.size()){
               HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
               Iterator<Integer> it = uniq.iterator();
               int index = 0;
               while(it.hasNext())
                  map.put(it.next(), index++);
               for(DiscreteSeq dseq : qseries)
                  dseq.remapSymbols(map);
            }
            gdata.nSymbols = uniq.size();
            System.err.printf("Number of symbols in quantized data: %d (%dms)\n", gdata.nSymbols, timer
                  .time());
            gdata.qseries = qseries;
            gdata.labQData = qdata;
         }
      }

      if (qseries == null){
         System.err.printf("Error: no quantized data (use -quant)\n");
         System.exit(1);
      }

      // make sure the quant and cont data matches
      if (qseries.size() != tseries.size()){
         System.err.printf(
               "Warning with data: # quant series data doesn't match # cont series (%d vs %d)\n", qseries
                     .size(), tseries.size());
      }

      if (qdata != null && data != null && qdata.size() != data.size() && !qdata.isEmpty()){
         System.err.printf(
               "Warning with labeled data: quant data doesn't match cont data (%d vs %d labeled seqs)\n",
               qdata.size(), data.size());
      }

      // create suffix tree
      System.err.print("Constructing suffix tree... ");

      // build alphabet
      timer.reset();
      StringBuffer sbAlpha = new StringBuffer();
      for(int i = 0; i < gdata.nSymbols; i++)
         sbAlpha.append((char)('a' + i));
      String alphabet = sbAlpha.toString();
      SuffixTree sufTree = new SuffixTree();

      gdata.sufTree = sufTree;
      TreeBuilder builder = new TreeBuilder(sufTree);
      for(DiscreteSeq dseq : qseries){
         String s = dseq.getAsString(alphabet);
         builder.addToken(s + (char)10000);
      }
      System.err.println("done.");
      System.err.printf("  # nodes = %d,  total length = %d,  # tokens = %d (%dms)\n",
            sufTree.getNumNodes(), sufTree.getTotalLength(), sufTree.getNumTokens(), timer.time());

      // create the window
      System.err.print("Creating visualization... ");
      timer.reset();
      JFrame frame = new JFrame(String.format("VizTools - %s (%s)", Library.getTitle(sDefFile), Library
            .getTitle(sQuantDefFile)));
      gdata.frame = frame;
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      Library.centerSizeWin(frame, null, 0.7);
      frame.setLocation(frame.getLocation().x, 0);
      JPanel p, mainp = new JPanel(new BorderLayout());
      frame.setContentPane(mainp);
      JTabbedPane tabs = new JTabbedPane(JTabbedPane.TOP);
      mainp.add(tabs, BorderLayout.CENTER);

      tabs.add("Raw Data", new AllDataView(gdata));
      // tabs.add("Data Overview", new DataOverview(gdata));
      // tabs.add("Boundary Info", new BoundaryView(gdata));
      tabs.add("Labeled Data", new LabeledView(gdata));
      if (gdata.labData != null) tabs.add("Word Spotting", new WordSpottingView(gdata));
      tabs.add("Hash View", new HashView(gdata));
      // tabs.add("Subseq Tree", new TreeView(gdata));
      tabs.add("Discover", new DiscoverView(gdata));
      tabs.add("Chiu", new ChiuView(gdata));
      tabs.add("Tanaka", new TanakaView(gdata));
      // tabs.add("SubDim", new SubDimView(gdata));

      tabs.setSelectedComponent(gdata.chiuView); // TODO: for debug ease

      System.err.printf("done (%dms).\n", timer.time());
      frame.setVisible(true);
      return true;
   }

   public static void main(String args[])
   {
      Runtime rt = Runtime.getRuntime();
      System.err.printf("Java version: %s  (%s)\n", System.getProperty("java.version"), System
            .getProperty("java.vendor"));
      System.err.printf("VM: %s / %s\n", System.getProperty("java.vm.version"), System
            .getProperty("java.vm.name"));
      System.err.printf("Memory (free/total/max): %dk / %dk / %dk\n", rt.freeMemory() / 1024, rt
            .totalMemory() / 1024, rt.maxMemory() / 1024);

      String sQuantDefFile = null;
      String sVizGrp = null;
      int c;

      LongOpt[] longopts = new LongOpt[] { new LongOpt("help", LongOpt.NO_ARGUMENT, null, 1001),
            new LongOpt("quant", LongOpt.REQUIRED_ARGUMENT, null, 1002),
            new LongOpt("vizgrp", LongOpt.REQUIRED_ARGUMENT, null, 1003) };

      Getopt g = new Getopt("VizTool", args, "?", longopts, true);
      while((c = g.getopt()) != -1){
         switch(c){
         case '?':
         case 1001: // help
            usage();
            System.exit(1);
            break;
         case 1002: // quant
            sQuantDefFile = g.getOptarg();
            break;
         case 1003: // vizgrp
            sVizGrp = g.getOptarg();
            break;
         }
      }

      if (g.getOptind() == args.length){
         System.err.println("\nError: missing data def file");
         System.exit(1);
      }
      String sDefFile = args[g.getOptind()];

      ToolTipManager.sharedInstance().setInitialDelay(500);
      VizTool vsd = new VizTool();
      if (!vsd.setup(sDefFile, sQuantDefFile, sVizGrp)){
         System.err.println("VizTool setup failed.");
         System.exit(1);
      }
   }
}
