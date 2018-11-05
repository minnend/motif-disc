package mdisc.kdd2007;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.gui.*;
import kdm.io.*;
import kdm.metrics.*;
import kdm.mlpr.dataTree.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.tools.*;
import kdm.util.*;
import mdisc.VizTool.*;
import mdisc.aaai2007.*;
import kdm.mlpr.htk.*;
import java.util.regex.*;
import org.apache.commons.math.stat.*;

// TODO improve HTK interface (efficiency)
// - only write support files once
// TODO grow each motif after each iteration? slow?
// - or will better neighbors be enough?
// TODO explicitly remove silence

/** visualization for KNN-based density mode estimation */
public class KDD2007 implements ActionListener, ChangeListener, ScoreGraphListener
{
   public static final Runtime rt = Runtime.getRuntime();
   public static final String sDiscoverView = "discover";
   public static final String sReviewView = "review";

   public static final boolean bNormalizeData = false;
   public static final boolean bForcePrune = true;
   public static final double THRESH_GROW_VAR_RATIO = 1.21;
   public static final double THRESH_PRUNE_OVERLAP = 0.8;
   public static final double DEN_MERGE_FACTOR = 1.05;
   public static final int DEN_MERGE_SPEED = 0;
   public static final double HMM_GLOBAL_SDEV_BLEND = 0.98; // fraction of estimated stdev to use
   public static final double HMM_OVERLAP_PERCENT = 0.5;
   public static final double THRESH_OVERLAP = 0.5;
   public static final double THRESH_NO_OVERLAP = 0.05;
   public static final double CLEAN_FRAC = 0.25;
   public static final int NMIX = 3;
   public static final double bgVarScale = 5.0;
   public static final String BG = "bg";
   public static final Pattern rexIndex = Pattern.compile("^.*(\\d+).*$");
   public static final double rBand = 0.15;
   public static final MetricSeq metseq = new DTWSimple(rBand, MetricSeq.LengthPrep.extend);
   public static final MetricFV metric = new EuclideanFV(true);

   protected GlobalData gdata;
   protected JFrame frame;
   protected ColoredRangeBar crbGT[];
   protected Color[] colors;
   protected VertChooseContainer vcc;
   protected ScoreGraphComplex sgc;

   protected JSpinner spinLenMin, spinLenMax;
   protected JLabel lbWinLens;
   protected LabeledSlider lsSeed, lsNumSegs, lsStep, lsKNN, lsMotif, lsNumFind;
   protected JButton btFullKNN, btFastKNN, btSave, btLoad, btCheckOverlap, btLocalOpt;
   protected JButton btFindNext, btScoreMotifs, btReEst, btReEstBG, btReset;
   protected JButton btReviewView, btDiscoverView, btMarkSilence, btSilenceGraph;
   protected JCheckBox btcGrow, btcExact;
   protected JRadioButton btrOptAcc, btrOptRecall;
   protected MotifViz mviz;
   protected CardLayout clView;
   protected JPanel mainp;

   protected int step;
   protected int wlens[];
   protected NNInfo[][] knn;
   protected WindowLocation[] i2loc;
   protected Sequence[] ss;
   protected FeatureVec[] descriptor;
   protected ArrayList<AbstractHMM> allHMMs;
   protected AbstractHMM hmmBG;
   protected ArrayList<AbstractHMM> motifHMMs = new ArrayList<AbstractHMM>();
   protected ArrayList<MarkupSet> motifLabels = new ArrayList<MarkupSet>();
   protected ContRecRet crrBase = null;
   protected ArrayList<ArrayList<WindowLocation>> seeds; // TODO linked list for fast removal?
   protected MyDoubleList seedRadius;

   /** real valued dist between successive window lengths */
   protected double dWinLen = 0;

   protected MyDoubleList listInfoGain = new MyDoubleList();
   protected HtkSetupInfo hsi = null;

   /**
    * Load the data and setup the view
    * 
    * @param sDefFile data def file for continuous data
    * @return true if successful
    */
   public boolean setup(String sDefFile)
   {
      TimerMS timer = new TimerMS();
      gdata = new GlobalData();

      // load data
      TreeMap<String, ArrayList<Sequence>> data = LabeledDataLoader.load(new File(sDefFile));
      if (data == null){
         System.err.printf("Failed to load labeled data\n (%s)\n", sDefFile);
         return false;
      }

      // setup the global data structure
      gdata.labData = data;
      gdata.updateClassStrings();
      gdata.tseries = LabeledDataLoader.tseries;
      gdata.calcLengthStats(gdata.tseries);
      System.err.printf("Loaded %d time series (len: total=%d, avg=%d).\n", gdata.tseries.size(),
            gdata.nTotalSeriesLength, gdata.nAvgSeriesLength);
      gdata.sentences = LabeledDataLoader.sentences;

      // normalize the data
      gdata.calcGlobalStats();
      if (bNormalizeData){
         FeatureVec vMin = gdata.getGlobalMin();
         FeatureVec vMax = gdata.getGlobalMax();
         FeatureVec vRange = vMax.sub(vMin);
         for(Sequence seq : gdata.tseries){
            int T = seq.length();
            for(int t = 0; t < T; t++){
               FeatureVec x = seq.get(t);
               x._sub(vMin)._div(vRange);
            }
         }
         gdata.calcGlobalStats();
      }

      LabeledDataLoader.dumpDataSummary(gdata.labData);

      // create the window
      System.err.print("Creating visualization... ");
      timer.reset();
      frame = new JFrame(String.format("KDD2007 - %s", Library.getTitle(sDefFile)));
      gdata.frame = frame;
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setSize(1200, 700);
      Library.centerWin(frame, null);
      frame.setLocation(frame.getLocation().x, 0);
      clView = new CardLayout();
      mainp = new JPanel(clView);
      frame.setContentPane(mainp);

      mainp.add(sDiscoverView, buildDiscoverGUI());
      mainp.add(sReviewView, buildReviewGUI());

      System.err.printf("done (%dms).\n", timer.time());
      frame.setVisible(true);
      return true;
   }

   protected void createSimpleBgModel()
   {
      System.err.printf("Learning simple background model... ");
      TimerMS timer = new TimerMS();
      int nDims = gdata.getNumDims();
      hmmBG = new HmmLR(1, nDims);
      ((HmmLR)hmmBG).setMinVar(gdata.getMinVar());
      hmmBG.setName(BG);
      hmmBG.setPiLeave(0, Math.log(0.2));
      GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      dg.setMean(gdata.getGlobalMean());
      dg.setVar(gdata.getGlobalVar().mul(bgVarScale));
      System.err.printf("done (%dms).\n", timer.time());
   }

   protected void createMixtureBgModel(int nMix)
   {
      System.err.printf("Learning background mixture model (%d)... ", nMix);
      TimerMS timer = new TimerMS();
      int nDims = gdata.getNumDims();
      hmmBG = new HmmLR(1, nDims);
      ((HmmLR)hmmBG).setMinVar(gdata.getMinVar());
      hmmBG.setName(BG);
      hmmBG.setPiLeave(0, Math.log(0.2));
      Sequence all = new Sequence();
      for(Sequence seq : gdata.tseries)
         all.append(seq);
      GMM gmm = null;
      double llWorst = Library.INF;
      double llBest = Library.NEGINF;
      for(int i = 0; i < 3; i++){
         GMM g = learnMixtureModel(nMix, all);
         double loglik = g.eval(all);
         if (loglik > llBest){
            llBest = loglik;
            gmm = g;
         }
         if (loglik < llWorst) llWorst = loglik;
      }
      hmmBG.setState(0, gmm);
      mviz.setHMM(hmmBG);
      System.err.printf("done (%dms).\n", timer.time());
      System.err.printf("loglik = %f  (dworst=%f)\n", gmm.eval(all), llBest - llWorst);
   }

   /** @return mixture model fit to the given data */
   protected GMM learnMixtureModel(int nMix, Sequence seq)
   {
      int nDims = gdata.getNumDims();
      GMM gmm = new GMM(nDims, nMix);
      for(int i = 0; i < nMix; i++){
         GaussianDiagonal g = new GaussianDiagonal(nDims);
         FeatureVec x = seq.get(Library.random(seq.length()));
         // FeatureVec x = new FeatureVec(gdata.getGlobalMean());
         FeatureVec noise = gdata.getGlobalVar().div(100);
         x._add(noise.mul(Library.random())._sub(noise.div(2)));
         g.setMean(x);
         g.setVar(gdata.getGlobalVar().mul(3)); // large to start to allow early adaptation
         gmm.setComp(i, 1.0 / nMix, g);
      }
      gmm.learn(seq);
      for(int i = 0; i < nMix; i++){
         FeatureVec var = ((GaussianDiagonal)gmm.getComp(i)).getVar();
         // enforce global minimum and BG variance scaling
         var._max(gdata.getMinVar())._mul(bgVarScale);
         ((GaussianDiagonal)gmm.getComp(i)).setVar(var);
      }
      return gmm;
   }

   protected JPanel buildDiscoverGUI()
   {
      JPanel mainp = new JPanel(new BorderLayout());
      JPanel leftp = new JPanel(new BorderLayout());

      // create the center view
      colors = Library.generateColors(gdata.getNumClasses());
      int nSeqs = gdata.getNumSeqs();
      crbGT = new ColoredRangeBar[nSeqs];
      JPanel seqsp = new JPanel(new VerticalLayout());
      seqsp.setBackground(Color.darkGray.darker());
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = gdata.tseries.get(iSeq);
         // create the labeled bars
         crbGT[iSeq] = new ColoredRangeBar(new Range(0, seq.length() - 1));
         crbGT[iSeq].setBackground(Color.black);
         crbGT[iSeq].setPreferredHeight(7);
         crbGT[iSeq].setHighlightLoc(ColoredRangeBar.Highlight.Top);
         crbGT[iSeq].setHighlightColor(Color.yellow);
         crbGT[iSeq].setHighlightHeight(2);
         crbGT[iSeq].setBorder(BorderFactory.createMatteBorder(1, 0, 2, 0, crbGT[iSeq].getBackground()));
         seqsp.add(crbGT[iSeq]);
         if (iSeq + 1 < gdata.tseries.size()) seqsp.add(Box.createVerticalStrut(1));
      }

      // add labeled ranges to the bar
      Iterator<String> itClass = gdata.labData.keySet().iterator();
      for(int iClass = 0; iClass < gdata.labData.size(); iClass++){
         String sClass = itClass.next();
         ArrayList<Sequence> subs = gdata.labData.get(sClass);
         for(Sequence sub : subs){
            int iSeq = sub.getParentIndex();
            Range r = new Range(sub.getParentOffset(), sub.getParentOffset() + sub.length() - 1);
            r.payload = new Pair(sClass, colors[iClass].darker());
            crbGT[iSeq].add(r, false, false, false);
         }
      }
      for(int i = 0; i < nSeqs; i++){
         crbGT[i].sortBars();
         crbGT[i].updateToolTip();
      }
      JScrollPane sp = new JScrollPane(seqsp);
      sp.setBorder(null);
      leftp.add(sp, BorderLayout.CENTER);

      // create the right control panel
      vcc = new VertChooseContainer();
      JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftp, vcc);
      split.setDividerLocation(Math.max(gdata.frame.getWidth() / 2, gdata.frame.getWidth() - 280));
      split.setResizeWeight(1.0);
      split.setDividerSize(5);
      split.setBorder(null);
      mainp.add(split, BorderLayout.CENTER);

      // create the control pane for knn
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      sp = new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      sp.setBorder(null);
      vcc.addPane("KNN", sp);

      Box p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      spinLenMin = new JSpinner(new SpinnerNumberModel(16, 4, 16, 1));
      spinLenMin.addChangeListener(this);
      spinLenMin.requestFocusInWindow();
      p.add(new JLabel("Window: "));
      p.add(spinLenMin);
      p.add(new JLabel(" to "));
      spinLenMax = new JSpinner(new SpinnerNumberModel(16, 16, 250, 1));
      spinLenMax.addChangeListener(this);
      p.add(spinLenMax);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);
      lbWinLens = new JLabel("", JLabel.CENTER);
      rightp.add(lbWinLens);
      lsStep = new LabeledSlider("Step Size: %d", 1, 125);
      rightp.add(lsStep);
      lsKNN = new LabeledSlider("K: %d", 1, 20);
      rightp.add(lsKNN);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btFastKNN = new JButton("  Fast KNNs  ");
      btFastKNN.addActionListener(this);
      p.add(btFastKNN);
      p.add(Box.createHorizontalStrut(8));
      btFullKNN = new JButton("Full KNNs");
      btFullKNN.addActionListener(this);
      p.add(btFullKNN);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btLoad = new JButton("Load KNN");
      btLoad.addActionListener(this);
      p.add(btLoad);
      p.add(Box.createHorizontalStrut(8));
      btSave = new JButton("Save KNN");
      btSave.addActionListener(this);
      btSave.setEnabled(false);
      p.add(btSave);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);
      btCheckOverlap = new JButton("Check for Overlapping Windows");
      btCheckOverlap.addActionListener(this);
      btCheckOverlap.setEnabled(false);
      rightp.add(btCheckOverlap);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btMarkSilence= new JButton("Mark Silence");
      btMarkSilence.addActionListener(this);
      p.add(btMarkSilence);
      p.add(Box.createHorizontalStrut(8));
      btSilenceGraph = new JButton("Silence Graph");
      btSilenceGraph.addActionListener(this);
      p.add(btSilenceGraph);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      // create the control pane for motif discovery
      rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      sp = new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      sp.setBorder(null);
      vcc.addPane("Motif Discovery", sp);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      lsNumSegs = new LabeledSlider("# segments: %d", 2, 50);
      p.add(lsNumSegs);
      btReset = new JButton("Reset");
      btReset.addActionListener(this);
      p.add(btReset);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btLocalOpt = new JButton("Init Seeds");
      btLocalOpt.addActionListener(this);
      btLocalOpt.setEnabled(false);
      p.add(btLocalOpt);
      p.add(Box.createHorizontalStrut(8));
      btcGrow = new JCheckBox("Grow Seeds", true);
      btcGrow.setEnabled(false);
      p.add(btcGrow);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      lsSeed = new LabeledSlider("Seed: %d\n", 0, 0);
      lsSeed.addChangeListener(this);
      lsSeed.setEnabled(false);
      rightp.add(lsSeed);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      lsNumFind = new LabeledSlider("# Motifs to Find: %d", 1, 31);
      lsNumFind.addChangeListener(this);
      lsNumFind.setEnabled(false);
      p.add(lsNumFind);
      p.add(Box.createHorizontalStrut(8));
      btcExact = new JCheckBox("Exact", false);
      btcExact.setEnabled(false);
      p.add(btcExact);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      btFindNext = new JButton("Find Next Motif");
      btFindNext.addActionListener(this);
      btFindNext.setEnabled(false);
      rightp.add(btFindNext);

      lsMotif = new LabeledSlider("No Motifs", 0, 0);
      lsMotif.addChangeListener(this);
      lsMotif.setEnabled(false);
      rightp.add(lsMotif);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btReEst = new JButton("ReEst Motifs");
      btReEst.addActionListener(this);
      btReEst.setEnabled(false);
      p.add(btReEst);
      p.add(Box.createHorizontalStrut(8));
      btReEstBG = new JButton("ReEst BG");
      btReEstBG.addActionListener(this);
      btReEstBG.setEnabled(false);
      p.add(btReEstBG);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      btScoreMotifs = new JButton("Score Motifs");
      btScoreMotifs.addActionListener(this);
      btScoreMotifs.setEnabled(false);
      p.add(btScoreMotifs);
      p.add(Box.createHorizontalStrut(8));
      btReviewView = new JButton("Review");
      btReviewView.addActionListener(this);
      // btReviewView.setEnabled(false);
      p.add(btReviewView);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      p = Box.createHorizontalBox();
      p.add(Box.createHorizontalGlue());
      p.add(new JLabel("Optimize: "));
      btrOptAcc = new JRadioButton("Accuracy", true);
      btrOptRecall = new JRadioButton("Recall", false);
      btrOptAcc.setEnabled(false);
      btrOptRecall.setEnabled(false);
      ButtonGroup bg = new ButtonGroup();
      bg.add(btrOptAcc);
      bg.add(btrOptRecall);
      p.add(btrOptAcc);
      p.add(btrOptRecall);
      p.add(Box.createHorizontalGlue());
      rightp.add(p);

      int nDims = gdata.getNumDims();
      // mviz = new MotifViz(null, FeatureVec.zeros(nDims), FeatureVec.ones(nDims)._mul(255));
      mviz = new MotifViz(null, gdata.getGlobalMin(), gdata.getGlobalMax());
      mviz.setBorder(BorderFactory.createLoweredBevelBorder());
      rightp.add(new Constrainer(mviz, 160, 160));

      // set initial values for various components
      vcc.show(0);
      lsNumSegs.setValue(8);
      spinLenMin.setValue(24);
      spinLenMax.setValue(40);
      lsStep.setValue(2);
      lsKNN.setValue(5);
      lsNumFind.setValue(1);
      calcWinLens((Integer)spinLenMin.getValue(), (Integer)spinLenMax.getValue());
      return mainp;
   }

   protected JPanel buildReviewGUI()
   {
      JPanel mainp = new JPanel(new BorderLayout());
      JPanel leftp = new JPanel(new BorderLayout());

      // create bottom view
      sgc = new ScoreGraphComplex();
      sgc.addScoreGraphListener(this);
      sgc.addCurve("Info Gain", listInfoGain);
      leftp.add(sgc, BorderLayout.SOUTH);

      // create right view
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      JScrollPane sp = new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      sp.setBorder(null);

      JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftp, sp);
      split.setDividerLocation(Math.max(gdata.frame.getWidth() / 2, gdata.frame.getWidth() - 280));
      split.setResizeWeight(1.0);
      split.setDividerSize(5);
      split.setBorder(null);
      mainp.add(split, BorderLayout.CENTER);

      btDiscoverView = new JButton("Discover View...");
      btDiscoverView.addActionListener(this);
      // btDiscoverView.setEnabled(false);
      rightp.add(btDiscoverView);

      return mainp;
   }

   protected void calcWinLens(int wmin, int wmax)
   {
      if (wmin >= wmax){
         wlens = new int[1];
         wlens[0] = wmin;
      }
      else{
         int range = wmax - wmin;
         int n = (int)Math.round(Math.max(Math.min((double)(range + 1) / 2, 5), 1));
         wlens = new int[n];
         for(int i = 0; i < n; i++)
            wlens[i] = wmin + (int)Math.round((double)i * range / (n - 1));
      }

      // update the label
      StringBuffer sb = new StringBuffer("Lengths: [");
      for(int i = 0; i < wlens.length; i++){
         if (i > 0) sb.append(", ");
         sb.append(wlens[i]);
      }
      sb.append(']');
      lbWinLens.setText(sb.toString());
   }

   /**
    * search for pairs of clusters (point + knn) that are too close (lambda*d(knn)); prune the lower density
    * cluster
    * 
    * @param a list of indices of cluster centers; -1 implies it's been pruned
    * @param K number of nearest neighbors
    * @param lambda multiplier for distance comparison (d(i,j) < d(i, knn(i))*lambda => prune)
    * @param speed 3=all (slowest), 2=no knn-to-knn, 1=no window-to-knn (fastest), 0=no prune
    */
   protected int pruneClusterPairs(int[] a, int K, double lambda, int speed, boolean bForcePrune)
   {
      if (speed == 0) return 0; // don't prune at all

      int nPruned = 0;
      // check inter-cluster dist; if close to knn dist, prune lower density cluster
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue; // already removed this index
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // calc min dist from knn(a[i]) to knn(a[j])
            double dmin = Library.INF;
            if (speed > 2){
               for(int ik = 0; ik < K; ik++)
                  for(int jk = 0; jk < K; jk++){
                     double d = calcSubseqDist(ss[knn[a[i]][ik].index], ss[knn[a[j]][jk].index]);
                     if (d < dmin) dmin = d;
                  }
            }

            if (speed > 1){
               // min over a[i] to knn(a[j]) and a[j] to knn(a[i])
               for(int ik = 0; ik < K; ik++){
                  double d = calcSubseqDist(ss[a[j]], ss[knn[a[i]][ik].index]);
                  if (d < dmin) dmin = d;
                  d = calcSubseqDist(ss[a[i]], ss[knn[a[j]][ik].index]);
                  if (d < dmin) dmin = d;
               }
            }

            // min over a[i] to a[j]
            double d = calcSubseqDist(ss[a[i]], ss[a[j]]);
            if (d < dmin) dmin = d;

            dmin /= lambda;
            if (dmin < knn[a[i]][K - 1].dist || dmin < knn[a[j]][K - 1].dist){
               int iPrune = pruneOne(a, i, j, bForcePrune);
               if (iPrune != -1){
                  nPruned++;
                  if (iPrune == i) break;
               }
            }
         }
      }
      return nPruned;
   }

   /**
    * search for connected components of clusters (point + knn) that are too close (lambda*d(knn)) and merge
    * them
    * 
    * @param a list of indices of cluster centers; -1 implies it's been pruned
    * @param K number of nearest neighbors
    * @param lambda multiplier for distance comparison (d(i,j) < d(i, knn(i))*lambda => merge)
    * @param speed 3=all (slowest), 2=no knn-to-knn, 1=no window-to-knn (fastest)
    * @return list of sets of window indices (from 'a', that is, a[i] rather than 'i')
    */
   protected ArrayList<HashSet<Integer>> mergeConnectedComponents(int[] a, int K, double lambda, int speed)
   {
      AList2<Integer> adj = new AList2<Integer>();

      // check inter-cluster dist; if close to knn dist, prune lower density cluster
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue;
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // calc min dist from a[i] to a[j]
            double dmin = Library.INF;
            if (speed > 2){
               for(int ik = 0; ik < K; ik++)
                  for(int jk = 0; jk < K; jk++){
                     double d = calcSubseqDist(ss[knn[a[i]][ik].index], ss[knn[a[j]][jk].index]);
                     if (d < dmin) dmin = d;
                  }
            }

            if (speed > 1){
               for(int ik = 0; ik < K; ik++){
                  double d = calcSubseqDist(ss[a[j]], ss[knn[a[i]][ik].index]);
                  if (d < dmin) dmin = d;
                  d = calcSubseqDist(ss[a[i]], ss[knn[a[j]][ik].index]);
                  if (d < dmin) dmin = d;
               }
            }

            double d = calcSubseqDist(ss[a[i]], ss[a[j]]);
            if (d < dmin) dmin = d;

            // are they close enough?
            dmin /= lambda;
            if (dmin < knn[a[i]][K - 1].dist || dmin < knn[a[j]][K - 1].dist) adj.add(i, j);
         }
      }

      // now search for connected components
      ArrayList<HashSet<Integer>> ret = new ArrayList<HashSet<Integer>>();
      boolean[] used = new boolean[a.length];
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0 || used[i]) continue;

         HashSet<Integer> members = new HashSet<Integer>();
         findConnComps(a, i, adj, used, members);
         ret.add(members);
      }
      return ret;
   }

   /** helper function for mergeConnectedComponents */
   protected void findConnComps(int[] a, int ix, AList2<Integer> adj, boolean[] used, Set<Integer> members)
   {
      if (used[ix]) return;
      used[ix] = true;
      members.add(a[ix]);
      int n = adj.size(ix);
      for(int i = 0; i < n; i++)
         findConnComps(a, adj.get(ix, i), adj, used, members);
   }

   /** if a[i] overlaps a[j]+knn or a[j] overlaps a[i]+knn, prune the one with lower density */
   protected int pruneOverlapKNN(int[] a, int K, double minPO, boolean bForcePrune)
   {
      int nPruned = 0;
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue; // already removed this index
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // check overlap: a[i] <-> a[j]
            int n = i2loc[a[i]].getNumOverlap(i2loc[a[j]]);
            double po = (double)n / Math.min(i2loc[a[i]].length(), i2loc[a[j]].length());
            boolean bOverlap = (po >= minPO);

            // check overlap: a[i] <-> a[j].knn
            if (!bOverlap) for(int k = 0; k < K; k++){
               n = i2loc[a[i]].getNumOverlap(i2loc[knn[a[j]][k].index]);
               po = (double)n / Math.min(i2loc[a[i]].length(), i2loc[knn[a[j]][k].index].length());
               if (po >= minPO){
                  bOverlap = true;
                  break;
               }
            }

            // check overlap: a[j] <-> a[i].knn
            if (!bOverlap) for(int k = 0; k < K; k++){
               n = i2loc[a[j]].getNumOverlap(i2loc[knn[a[i]][k].index]);
               po = (double)n / Math.min(i2loc[a[j]].length(), i2loc[knn[a[i]][k].index].length());
               if (po >= minPO){
                  bOverlap = true;
                  break;
               }
            }

            if (bOverlap){
               int iPrune = pruneOne(a, i, j, bForcePrune);
               if (iPrune != -1){
                  nPruned++;
                  if (iPrune == i) break;
               }
            }
         }
      }
      return nPruned;
   }

   /** if any pair of window locations from two seeds overlap, prune one seed */
   protected int pruneOverlapSeeds(double minPO, boolean bForcePrune)
   {
      // TODO prune or merge? which seed? should calc seed radius
      int nPruned = 0;
      for(int iSeed = 0; iSeed < seeds.size(); iSeed++){
         ArrayList<WindowLocation> seedsi = seeds.get(iSeed);
         for(int jSeed = iSeed + 1; jSeed < seeds.size();){
            ArrayList<WindowLocation> seedsj = seeds.get(jSeed);

            boolean bOverlap = false;
            search: for(WindowLocation loci : seedsi)
               for(WindowLocation locj : seedsj){
                  int n = loci.getNumOverlap(locj);
                  double p = (double)n / Math.min(loci.length(), locj.length());
                  if (p >= minPO){
                     bOverlap = true;
                     break search;
                  }
               }
            if (bOverlap){
               if (seedRadius.get(iSeed) <= seedRadius.get(jSeed)){
                  seeds.remove(jSeed);
                  seedRadius.removeElementAt(jSeed);
               }
               else{
                  seeds.remove(iSeed);
                  seedRadius.removeElementAt(iSeed);
                  iSeed--;
                  break;
               }
            }
            else jSeed++;
         }
      }
      return nPruned;
   }

   protected int pruneSimilarSeeds(double mergeFactor, int mergeSpeed, boolean bForcePrune)
   {
      // TODO
      return 0;
   }

   /**
    * Select index (i or j) to prune and set a[iPruned] = -1
    * 
    * @return index (i or j) that was pruned, -1 if no clear decision
    */
   protected int pruneOne(int[] a, int i, int j, boolean bForced)
   {
      int K = knn[i].length;

      int lenI = i2loc[a[i]].length();
      int lenJ = i2loc[a[j]].length();
      double distI = knn[a[i]][K - 1].dist;
      double distJ = knn[a[j]][K - 1].dist;

      // is there a trivial (>= length, < dist) winner?
      if (lenI >= lenJ && distI <= distJ){
         a[j] = -1;
         return j;
      }
      if (lenJ >= lenI && distJ <= distI){
         a[i] = -1;
         return i;
      }

      // no trivial winner, but maybe not a tough decision...

      // close in density, so decide based on length
      if (Math.abs(distI - distJ) / Math.max(distI, distJ) < 0.1){
         assert (lenI != lenJ) : String.format("lenI=%d  lenJ=%d  distI=%f  distJ=%f", lenI, lenJ, distI,
               distJ);
         if (lenI > lenJ){
            a[j] = -1;
            return j;
         }
         else{
            a[i] = -1;
            return i;
         }
      }

      // close in length, so decide based on density
      if (Math.abs(lenI - lenJ) <= (int)Math.ceil(dWinLen)){
         if (distI < distJ){ // smaller dist = higher density
            a[j] = -1;
            return j;
         }
         else{
            a[i] = -1;
            return i;
         }
      }

      // maybe one is much better, so ignore length
      if (Math.abs(distI - distJ) / Math.min(distI, distJ) > 0.5){
         if (distI < distJ){ // smaller dist = higher density
            a[j] = -1;
            return j;
         }
         else{
            a[i] = -1;
            return i;
         }
      }

      // no easy choice...
      if (bForced){
         if (distI < distJ){ // TODO forced choice based on length or dist?
            a[j] = -1;
            return j;
         }
         else{
            a[i] = -1;
            return i;
         }
      }

      // System.err.printf("Can't decide: %d, %d %.2f, %.2f\n", lenI, lenJ, distI, distJ);
      return -1;
   }

   /** find local optima in density, prune, and learn HMMs */
   protected void calcLocalOptsHmms()
   {
      boolean bGrow = btcGrow.isSelected();
      int N = knn.length;
      int K = knn[0].length;
      SpanList spanLocalOpt = new SpanList(0, N - 1, false);
      TimerMS timer = new TimerMS();

      // find local optima
      timer.reset();
      for(int i = 0; i < N; i++){
         assert (knn[i][K - 1].dist >= knn[i][0].dist);
         double deni = knn[i][K - 1].dist;
         boolean bLocalOpt = true;
         for(int k = 0; k < K; k++){
            int j = knn[i][k].index;
            double d = knn[j][K - 1].dist;
            if (d < deni){
               bLocalOpt = false;
               break;
            }
         }
         if (bLocalOpt) spanLocalOpt.add(i);
      }
      int nLocalOpt = spanLocalOpt.size();
      System.err.printf("Found %d local optima in %d data points (%.1f%%, %dms).\n", nLocalOpt, N, 100.0
            * nLocalOpt / N, timer.time());

      // prune (temporally) overlapping windows; keep highest density
      int[] a = spanLocalOpt.toIndexArray();
      timer.reset();
      int nPruned = pruneOverlapKNN(a, K, THRESH_PRUNE_OVERLAP, bForcePrune);
      nLocalOpt -= nPruned;
      System.err.printf("Pruned %d local optima due to overlapping knn (%d left, %.2f%%, %dms).\n", nPruned,
            nLocalOpt, 100.0 * nLocalOpt / N, timer.time());

      // try to reduce the number of local optima by detecting redundancy
      timer.reset();
      nPruned = pruneClusterPairs(a, K, DEN_MERGE_FACTOR, DEN_MERGE_SPEED, bForcePrune);
      nLocalOpt -= nPruned;
      ArrayList<HashSet<Integer>> seedIndices = new ArrayList<HashSet<Integer>>();
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue;
         HashSet<Integer> set = new HashSet<Integer>();
         set.add(a[i]);
         seedIndices.add(set);
      }
      assert (nLocalOpt == seedIndices.size());
      System.err.printf("Pruned %d local optima due to density clustering (%d left, %.2f%%, %dms).\n",
            nPruned, nLocalOpt, 100.0 * nLocalOpt / N, timer.time());

      /*
       * ArrayList<HashSet<Integer>> seedIndices = mergeConnectedComponents(a,K,DEN_MERGE_FACTOR,
       * DEN_MERGE_SPEED); nLocalOpt = seedIndices.size(); System.err.printf("Connected components: %d
       * (%.2f%%, %dms)\n", seedIndices.size(), 100.0*nLocalOpt/N,timer.time());
       */

      // convert the seed indices into window locations
      seeds = new ArrayList<ArrayList<WindowLocation>>();
      for(HashSet<Integer> group : seedIndices){
         ArrayList<WindowLocation> wlocs = new ArrayList<WindowLocation>();
         Iterator<Integer> it = group.iterator();
         while(it.hasNext()){
            int ix = it.next();
            wlocs.add(new WindowLocation(i2loc[ix]));
            for(int k = 0; k < knn[ix].length; k++)
               wlocs.add(new WindowLocation(i2loc[knn[ix][k].index]));
         }
         seeds.add(wlocs);
      }

      // grow seeds (temporal extension) if requested
      if (bGrow){
         System.err.printf("Growing %d seeds... ", seeds.size());
         timer.reset();
         int iSeed = 0;
         for(ArrayList<WindowLocation> wlocs : seeds)
            grow(iSeed++, wlocs, THRESH_GROW_VAR_RATIO);
         System.err.printf("done (%dms).\n", timer.time());

         // compute seed radii
         seedRadius = new MyDoubleList();
         System.err.printf("Computing seed radii (%d)... ", seeds.size());
         timer.reset();
         for(ArrayList<WindowLocation> wlocs : seeds){
            ArrayList<Sequence> seqs = Sequence.extractPats(gdata.tseries, wlocs);
            int nLocs = wlocs.size();
            int iBase = Library.random(nLocs);
            int iFar = -1;
            double vFar = 0;
            for(int i = 0; i < nLocs; i++){
               if (i == iBase) continue;
               double v = calcSubseqDist(seqs.get(i), seqs.get(iBase));
               if (v > vFar){
                  vFar = v;
                  iFar = i;
               }
            }
            iBase = iFar;
            iFar = -1;
            vFar = 0;
            for(int i = 0; i < nLocs; i++){
               if (i == iBase) continue;
               double v = calcSubseqDist(seqs.get(i), seqs.get(iBase));
               if (v > vFar){
                  vFar = v;
                  iFar = i;
               }
            }
            seedRadius.add(vFar / 2);
         }
         assert (seeds.size() == seedRadius.size());
         System.err.printf("done (%dms).\n", timer.time());

         System.err.printf("Pruning %d extended seeds... ", seeds.size());
         timer.reset();
         pruneOverlapSeeds(THRESH_PRUNE_OVERLAP, bForcePrune);
         assert (seeds.size() == seedRadius.size());
         System.err.printf("(After Overlap: %d)... ", seeds.size());
         // pruneSimilarSeeds(DEN_MERGE_FACTOR, DEN_MERGE_SPEED, bForcePrune);
         nLocalOpt = seeds.size();
         System.err.printf("done (%d left, %dms).\n", nLocalOpt, timer.time());
      }

      // learn HMMs for the remaining seeds
      System.err.printf("Learning HMMs (%d)... ", nLocalOpt);
      timer.reset();
      allHMMs = new ArrayList<AbstractHMM>(nLocalOpt);
      for(int i = 0; i < nLocalOpt; i++)
         allHMMs.add(genHMM(seeds.get(i), String.format("HMM%05d", i)));
      System.err.printf("done (%dms).\n", timer.time());
   }

   /** temporally extend the window locations */
   protected void grow(int iSeed, ArrayList<WindowLocation> wlocs, double maxVarRatio)
   {
      int nDims = gdata.getNumDims();
      int nLocs = wlocs.size();

      // find the shortest sequence
      int minLen = wlocs.get(0).length();
      for(int i = 1; i < nLocs; i++){
         int len = wlocs.get(i).length();
         if (len < minLen) minLen = len;
      }
      int nStates = minLen;

      // learn a left-right HMM for alignment
      // TODO proper topology?
      HmmLR hmm = new HmmLR(nStates, 0, nDims);
      ArrayList<Sequence> examples = WindowLocation.getExamples(wlocs, gdata.tseries);
      hmm.init_segk(examples);

      // find max variance in HMM
      MyDoubleList proj = new MyDoubleList();
      MyDoubleList projOrig = new MyDoubleList();
      FeatureVec maxVar = FeatureVec.zeros(nDims);
      FeatureVec sumVar = FeatureVec.zeros(nDims);
      int nGrow = 0;
      for(int i = 0; i < nStates; i++){
         FeatureVec fvVar = ((GaussianDiagonal)hmm.getState(i)).getVar();
         maxVar._max(fvVar);
         sumVar._add(fvVar);
         proj.add(fvVar.abs().sum());
         projOrig.add(fvVar.abs().sum());
      }

      FeatureVec meanVar = sumVar.div(nStates);
      double maxNorm = maxVar.abs().sum() * maxVarRatio;
      double meanNorm = meanVar.abs().sum() * maxVarRatio;

      int iLeft = 0, iRight = 0;
      int nPruneLeft = 0, nPruneRight = 0;
      // scan left of each loc and see if the variance is small enough
      {
         int nMaxLeft = Integer.MAX_VALUE;
         for(WindowLocation wloc : wlocs){
            int ix = wloc.getFirstIndex();
            if (ix < nMaxLeft) nMaxLeft = ix;
         }
         while(iLeft < nMaxLeft){
            Sequence data = new Sequence();
            for(WindowLocation wloc : wlocs)
               data.add(gdata.tseries.get(wloc.iSeries).get(wloc.getFirstIndex() - iLeft - 1));
            FeatureVec fvVar = data.getVar();
            sumVar._add(fvVar);
            FeatureVec fvMeanGrow = sumVar.div(nStates + iLeft);
            double v = fvMeanGrow.abs().sum();
            if (v > meanNorm) break;
            else{
               iLeft++;
               projOrig.add(0, 0);
               proj.add(0, fvVar.abs().sum());
            }
         }
         assert (iLeft <= nMaxLeft) : String.format("iLeft=%d  max=%d", iLeft, nMaxLeft);

         // we might need to trim
         while(iLeft - nPruneLeft > -minLen / 2 + 1){
            double v = proj.get(nPruneLeft);
            if (v <= meanNorm) break;
            nPruneLeft++;
         }

         if (iLeft - nPruneLeft != 0){
            nGrow += (iLeft - nPruneLeft);
            for(WindowLocation wloc : wlocs){
               // System.err.printf("grow seed %d left (%d): %d (%d) -> %d (%d)\n", iSeed + 1, iLeft,
               // wloc.iStart, wloc.nLength, wloc.iStart - iLeft, wloc.nLength + iLeft);
               wloc.iStart -= (iLeft - nPruneLeft);
               wloc.nLength += (iLeft - nPruneLeft);
               assert (wloc.iStart >= 0);
            }
         }
      }

      // scan right of each loc and see if the variance is small enough
      {
         int nMaxRight = Integer.MAX_VALUE;
         for(WindowLocation wloc : wlocs){
            int n = gdata.seqLens[wloc.iSeries] - wloc.getLastIndex() - 1;
            if (n < nMaxRight) nMaxRight = n;
         }
         while(iRight < nMaxRight){
            Sequence data = new Sequence();
            for(WindowLocation wloc : wlocs)
               data.add(gdata.tseries.get(wloc.iSeries).get(wloc.getLastIndex() + iRight));
            FeatureVec fvVar = data.getVar();
            sumVar._add(fvVar);
            FeatureVec fvMeanGrow = sumVar.div(nStates + iRight);
            double v = fvMeanGrow.abs().sum();
            if (v > meanNorm) break;
            else{
               iRight++;
               projOrig.add(0);
               proj.add(fvVar.abs().sum());
            }
         }
         assert (iRight <= nMaxRight) : String.format("iRight=%d  max=%d", iRight, nMaxRight);

         // we might need to trim
         while(iRight - nPruneRight > -minLen / 2 + 1){
            double v = proj.get(proj.size() - 1 - nPruneRight);
            if (v <= meanNorm) break;
            nPruneRight++;
         }

         if (iRight - nPruneRight != 0){
            nGrow += (iRight - nPruneRight);
            for(WindowLocation wloc : wlocs){
               // System.err.printf("grow seed %d right (%d): %d (%d) -> %d (%d)\n", iSeed + 1, iRight,
               // wloc.iStart, wloc.nLength, wloc.iStart, wloc.nLength + iRight);
               wloc.nLength += iRight - nPruneRight;
               assert (wloc.getLastIndex() < gdata.seqLens[wloc.iSeries]);
            }
         }
      }

      // TODO debug
      /*
       * if (nGrow > 4){ System.err.printf("seed=%d len=%d iLeft=%d iRight=%d\n", iSeed + 1, proj.size(),
       * iLeft, iRight); ScatterPlot plot = new ScatterPlot(); Sequence seq = new Sequence("Orig",
       * projOrig.toArray()); plot.addData(seq, Color.green, ScatterPlot.DEF_STROKE, 'o'); seq = new
       * Sequence("Extended", proj.toArray()); plot.addData(seq, Color.blue, ScatterPlot.DEF_STROKE, 'x');
       * double max = seq.getMax().max(); plot.addLine(nPruneLeft, -max / 100, nPruneLeft, max,
       * Color.cyan).setStroke(new BasicStroke(2)); plot.addLine(iLeft, -max / 100, iLeft, max, Color.orange);
       * plot.addLine(seq.length() - 1 - nPruneRight, -max / 100, seq.length() - 1 - nPruneRight, max,
       * Color.cyan).setStroke(new BasicStroke(2)); plot.addLine(seq.length() - iRight - 1, -max / 100,
       * seq.length() - iRight - 1, max, Color.orange); plot.addLine(0, meanNorm, seq.length() - 1, meanNorm,
       * Color.magenta).setStroke(new BasicStroke(2)); Library.popup(String.format("Grow %d - #left=%d
       * #right=%d", iSeed + 1, iLeft, iRight), plot, 600, 400); }
       */
   }

   /*
    * find the next best motif using the current search method
    */
   protected boolean findNext(int nFindMax)
   {
      TimerMS timerTotal = new TimerMS();
      int iMotif = motifHMMs.size() + 1;

      // score the baseline (no additional motifs)
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      if (motifHMMs.isEmpty()){
         TimerMS timer = new TimerMS();
         if (hsi == null) hsi = HTK.setupContRec(null, hmms, gdata.tseries);
         crrBase = HTK.contRec(hsi, hmms, gdata.tseries);
         if (crrBase == null){
            System.err.println("Error: HTK / HVite failed");
            return false;
         }
         motifLabels = crrBase.markupSets;
         listInfoGain.add(0);
      }
      else hmms.addAll(motifHMMs);

      // score each window (seed) location and explain data
      String sMotifName = String.format("hmm%03d.txt", iMotif);
      int nSeeds = allHMMs.size();
      boolean bSetupHSI = false;
      if (nSeeds < 1) return false;
      ContRecRet[] crrs = new ContRecRet[nSeeds];
      TimerMS timer = new TimerMS();
      for(int iSeed = 0; iSeed < nSeeds; iSeed++){
         System.err.printf(" Testing candidate HMM %d / %d ... ", iSeed + 1, nSeeds);
         timer.reset();

         // run the continuous recognizer and save results
         AbstractHMM hmm = allHMMs.get(iSeed);
         hmm.setName(sMotifName);
         if (!bSetupHSI){
            hmms.add(hmm);
            hsi = HTK.setupContRec(hsi, hmms, gdata.tseries);
            bSetupHSI = true;
         }
         else{
            // need to replace the last HMM with a new candidate
            hmms.set(hmms.size() - 1, hmm);
            File fHmm = new File(HTK.fTmpDir, hmm.getName());
            assert (fHmm.exists()) : String.format("HMM file should exist (%s)", fHmm.getAbsolutePath());
            if (!HtkHmm.save(hmm, fHmm)) return false;
         }
         assert (hsi != null);
         crrs[iSeed] = HTK.contRecAfterSetup(hsi);
         if (crrs[iSeed] == null){
            System.err.printf("done (%dms) !! ERROR !!\n", timer.time());
            assert (false);
         }
         else System.err.printf("done (%dms)  (ll=%.1f)\n", timer.time(), crrs[iSeed].llTotal);
      }

      int nFound = 0;
      // greedily select best motif from remaining set
      while(nFound < nFindMax || nFindMax < 0){
         // find best motif (null => removed,
         int iBest = -1;
         for(int i = 0; i < nSeeds; i++){
            if (crrs[i] == null) continue;
            if (iBest < 0 || crrs[i].llTotal > crrs[iBest].llTotal) iBest = i;
         }
         if (iBest < 0) break; // no next best found
         System.err.printf("Adding seed: %d (gain: %f)\n", iBest, crrs[iBest].llTotal - crrBase.llTotal);
         nFound++;

         // initialize new motif from best window location
         allHMMs.get(iBest).setName(String.format("hmm%03d.txt", motifHMMs.size() + 1));
         motifHMMs.add(allHMMs.get(iBest));
         listInfoGain.add(crrs[iBest].llTotal - crrBase.llTotal);

         // calculate overlap percent to/from every other seed
         SpanList[] spansI = buildSpanLists(sMotifName, crrs[iBest].markupSets);
         double[] overIJ = new double[nSeeds];
         double[] overJI = new double[nSeeds];
         overIJ[iBest] = overJI[iBest] = 1.0; // perfect overlap with self
         for(int i = 0; i < nSeeds; i++){
            if (crrs[i] == null){
               overIJ[i] = overJI[i] = Double.NaN;
               continue;
            }
            if (i == iBest) continue;
            SpanList[] spansJ = buildSpanLists(sMotifName, crrs[i].markupSets);
            overIJ[i] = calcOverlap(spansI, spansJ);
            overJI[i] = calcOverlap(spansJ, spansI);
            // System.err.printf("%d vs. %d: %.2f, %.2f\n", iBest, i, overIJ[i], overJI[i]);
         }

         // remove all overlapping (>thresh1) seeds permanently
         for(int i = 0; i < nSeeds; i++){
            if (crrs[i] == null) continue;
            if (overIJ[i] >= THRESH_OVERLAP && overJI[i] >= THRESH_OVERLAP){
               if (i != iBest)
                  System.err.printf("Overlap: %d & %d  (%f, %f)\n", iBest, i, overIJ[i], overJI[i]);
               crrs[i] = null;
               allHMMs.set(i, null);
            }
         }
         assert (allHMMs.get(iBest) == null);
         assert (crrs[iBest] == null);

         // removel all non-unique (>thresh2) seeds from consideration in this round
         if (listInfoGain.size() > 0){
            for(int i = 0; i < nSeeds; i++){
               if (crrs[i] == null) continue;
               if (overIJ[i] <= THRESH_NO_OVERLAP && overJI[i] <= THRESH_NO_OVERLAP){
                  double baseGain = listInfoGain.last();
                  double igain = crrs[i].llTotal - crrBase.llTotal;
                  double igRatio = igain / baseGain;
                  if (igRatio < 0.8) crrs[i] = null;
                  else System.err.printf("No overlap: %d & %d igain=%f (d=%f, %.2f%%)\n", iBest, i, igain,
                        igain - baseGain, 100.0 * igRatio);
               }
               else crrs[i] = null;
            }
         }
         else break;
      }

      ContRecRet crr = runContRec(-1, true);
      /*
       * for(TimeMarker tm : crr.markupSets.get(0).getList()){ if
       * (tm.getTag().equals(motifHMMs.get(motifHMMs.size() - 1).getName())) System.err.printf("%d %d\n",
       * tm.getStartIndex(), tm.length()); }
       */
      sgc.updateGraph();
      return true;
   }

   /** @return percent of spansB covered by spansA */
   protected double calcOverlap(SpanList[] spansA, SpanList[] spansB)
   {
      int nTotal = 0, nOver = 0;
      assert (spansA.length == spansB.length) : String.format(
            "Error: span array lengths don't match (%d vs. %d)\n", spansA.length, spansB.length);
      for(int i = 0; i < spansA.length; i++){
         nTotal += spansB[i].size();
         SpanList span = spansA[i].intersect(spansB[i]);
         nOver += span.size();
      }
      return (double)nOver / nTotal;
   }

   /** @return span list for each markup set in the given list */
   protected SpanList[] buildSpanLists(String tag, ArrayList<MarkupSet> marks)
   {
      int n = marks.size();
      SpanList[] spans = new SpanList[n];
      for(int i = 0; i < n; i++){
         spans[i] = new SpanList(0, gdata.seqLens[i] - 1, false);
         for(TimeMarker tm : marks.get(i).getList()){
            if (tag.equals(tm.getTag())) spans[i].add(tm.getStartIndex(), tm.getStopIndex() - 1);
         }
      }
      return spans;
   }

   protected int setupWindows()
   {
      int nSeqs = gdata.getNumSeqs();
      int nWins = 0;
      for(int i = 0; i < nSeqs; i++)
         for(int j = 0; j < wlens.length; j++)
            nWins += Library.getNumSlidingWindowSites(gdata.seqLens[i], wlens[j], step);
      i2loc = new WindowLocation[nWins];
      ss = new Sequence[nWins];

      // extract windows
      int ix = 0;
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         for(int iLen = 0; iLen < wlens.length; iLen++){
            Sequence seq = gdata.tseries.get(iSeq);
            int T = gdata.seqLens[iSeq] - wlens[iLen];
            for(int t = 0; t <= T; t += step){
               i2loc[ix] = new WindowLocation(iSeq, t, wlens[iLen]);
               ss[ix] = seq.subseq(i2loc[ix]);
               ix++;
            }
         }
      }
      assert (ix == nWins);
      return nWins;
   }

   /** notify gui that we now have data */
   protected void updateGUIforData()
   {
      btSave.setEnabled(true);
      btCheckOverlap.setEnabled(true);
      btLocalOpt.setEnabled(true);
      btcGrow.setEnabled(true);
      lsNumSegs.setEnabled(true);
      lsSeed.setValues(0, 0, 0);
      lsSeed.setEnabled(false);
      lsMotif.setValues(0, 0, 0);
      lsMotif.setEnabled(false);
      lsNumFind.setEnabled(false);
      btcExact.setEnabled(false);
      btFindNext.setEnabled(false);
      btReEst.setEnabled(false);
      btReEstBG.setEnabled(false);
      btScoreMotifs.setText("Score Motifs");
      btScoreMotifs.setEnabled(false);
      mviz.setHMM(null);
      listInfoGain.clear();
      motifHMMs.clear();
      motifLabels.clear();
      crrBase = null;
      new Thread() {
         public void run()
         {
            createMixtureBgModel(NMIX); // TODO not exactly thread safe
         }
      }.start();
      highlightMotif(-1, null);
      sgc.updateGraph();
   }

   protected void allKNN(FeatureVec query, ArrayList<FeatureVec> data, NBestList<NNInfo> nn)
   {
      WindowLocation wQuery = i2loc[(Integer)query.getMeta("Index")];
      for(FeatureVec x : data){
         int j = (Integer)x.getMeta("Index");
         if (wQuery.overlaps(i2loc[j])) continue;
         double d = metric.dist(query, x);
         nn.add(new NNInfo(j, d, i2loc[j]));
      }
   }

   protected void fastKNN(FeatureVec query, VectorTreeNode node, NBestList<NNInfo> nn)
   {
      if (node.isLeaf()) allKNN(query, node.getData(), nn);
      else{
         VectorTreeNode left = node.getKid(0);
         VectorTreeNode right = node.getKid(1);
         FeatureVec mean1 = (FeatureVec)left.meta.get(MetricTree.KeyMean);
         FeatureVec mean2 = (FeatureVec)right.meta.get(MetricTree.KeyMean);
         double d1 = metric.dist(query, mean1);
         double d2 = metric.dist(query, mean2);
         VectorTreeNode good, bad;
         if (d1 < d2){
            good = left;
            bad = right;
         }
         else{
            good = right;
            bad = left;
         }
         fastKNN(query, good, nn);
         if (nn.size() < nn.getN()) fastKNN(query, bad, nn);
      }
   }

   /** @return feature vector that describes the given sequence */
   protected FeatureVec calcDescriptor(Sequence seq, int len)
   {
      int seqD = seq.getNumDims();
      int xD = len * seqD;
      FeatureVec x = new FeatureVec(xD);

      // procedure taken from segk_overlap for left-right HMM initialization
      double denom = len - (len - 1) * HMM_OVERLAP_PERCENT;
      double fpd = (double)seq.length() / denom;
      double foverlap = fpd * HMM_OVERLAP_PERCENT;
      int ix = 0;
      for(int iState = 0; iState < len; iState++){
         int a = (int)Math.round(iState * (fpd - foverlap));
         int b = (int)Math.round(iState * (fpd - foverlap) + fpd);
         b = Math.max(b, a + 1);
         int n = b - a;
         for(int j = 0; j < seqD; j++){
            if (n == 1) x.set(ix++, seq.get(a, j));
            else{
               double sum = 0;
               for(int k = a; k < b; k++)
                  sum += seq.get(k, j);
               x.set(ix++, sum / n);
            }
         }
      }
      assert (ix == xD);
      return x;
   }

   /** @return normalized distance between given subsequences */
   protected double calcSubseqDist(Sequence seq1, Sequence seq2)
   {
      int maxlen = Library.max(seq1.length(), seq2.length());
      return metseq.dist(seq1, seq2) / maxlen;
   }

   /** compute the KNN (non-overlapping) for every subsequence by computing all pair-wise distances */
   protected void calcFullKNN(int step, int K)
   {
      this.step = step;
      System.err.printf("Run Full KNN: #wlens=%d  k=%d\n", wlens.length, K);
      TimerMS timer = new TimerMS();
      TimerMS timer2 = new TimerMS();
      int nSeqs = gdata.getNumSeqs();

      System.err.println("Generating windows... ");
      timer2.reset();
      int nWins = setupWindows();
      System.err.printf("Generated %d windows in %dms.\n", nWins, timer2.time());
      knn = new NNInfo[nWins][K];

      System.err.printf("Collecting %d windows... ", nWins);
      timer2.reset();
      for(int i = 0; i < ss.length; i++)
         ss[i].setMeta("Index", i);
      System.err.printf("done (%dms).\n", timer2.time());
try{
      PrintWriter out = new PrintWriter(new FileWriter("d:/tidigits-kdd-full.txt"));
      
      // calc k smallest distances
      NBestList<NNInfo> nnSeq = new NBestList<NNInfo>(K);
      nnSeq.setBiggerIsBetter(false);
      nnSeq.setAvoidConflict(true);
      for(int iWin = 0; iWin < nWins; iWin++){
         if (iWin % 10 == 9) System.err.printf("Calculating KNN for window %d / %d\n", iWin + 1, nWins);
         WindowLocation wloc = i2loc[iWin];

         // search the nearest descriptors for the nearest sequences
         nnSeq.clear();
         for(int i = 0; i < ss.length; i++){
            WindowLocation wloc2 = i2loc[i];
            if (wloc.overlaps(wloc2)) continue;
            double dseq = calcSubseqDist(ss[iWin], ss[i]);
            nnSeq.add(new NNInfo(i, dseq, wloc2));
         }

         // save the real knn
         for(int j = 0; j < K; j++){
            knn[iWin][j] = nnSeq.get(j);
            out.printf("%d ", knn[iWin][j].index);
         }
         out.println();
         out.flush();
      }
      out.close();
}catch(Exception e){e.printStackTrace();}
      updateGUIforData();
      System.err.printf("Time for KNN computation: %s\n", Library.formatDuration(timer.time(), 2));
   }

   /** compute the KNN (non-overlapping) for every subsequence */
   protected void calcFastKNN(int step, int K)
   {
      this.step = step;
      System.err.printf("Run Fast KNN: #wlens=%d  k=%d\n", wlens.length, K);
      TimerMS timer = new TimerMS();
      TimerMS timer2 = new TimerMS();
      int nSeqs = gdata.getNumSeqs();

      System.err.printf("Pre-KNN) Mem (free/total/max): %dk / %dk / %dk\n", rt.freeMemory() / 1024, rt
            .totalMemory() / 1024, rt.maxMemory() / 1024);

      System.err.println("Generating windows... ");
      timer2.reset();
      int nWins = setupWindows();
      System.err.printf("Generated %d windows in %dms.\n", nWins, timer2.time());
      knn = new NNInfo[nWins][K];

      System.err.printf("Collecting %d windows... ", nWins);
      timer2.reset();
      for(int i = 0; i < ss.length; i++)
         ss[i].setMeta("Index", i);
      System.err.printf("done (%dms).\n", timer2.time());

      int descLen = calcHmmLen(wlens[0]);
      System.err.printf("Extracting sequence descriptors (%d, D=%d)... ", nWins, descLen);
      timer2.reset();
      ArrayList<FeatureVec> data = new ArrayList<FeatureVec>(nWins);
      descriptor = new FeatureVec[nWins];
      for(int i = 0; i < nWins; i++){
         descriptor[i] = calcDescriptor(ss[i], descLen);
         descriptor[i].setMeta("Index", i);
         data.add(descriptor[i]);
      }
      System.err.printf("done (%dms).\n", timer2.time());

      int KVec = K * 4;
      int nMaxMembers = Math.max((KVec + 1) * 10, (int)(3 * Math.sqrt(ss.length)));
      // int nMaxMembers = (K+1)*2;
      System.err.printf("Building tree (max members=%d, K=%d, KVec=%d)... ", nMaxMembers, K, KVec);
      timer2.reset();
      MetricTree mtree = new MetricTree();
      mtree.constructNaive(data, metric, nMaxMembers);
      System.err.printf("done (%dms, %d nodes).\n", timer2.time(), mtree.getNumNodes());

      // calc k smallest distances
      NBestList<NNInfo> nnVec = new NBestList<NNInfo>(KVec);
      nnVec.setBiggerIsBetter(false);
      nnVec.setAvoidConflict(true);
      NBestList<NNInfo> nnSeq = new NBestList<NNInfo>(K);
      nnSeq.setBiggerIsBetter(false);
      nnSeq.setAvoidConflict(true);
      for(int iWin = 0; iWin < nWins; iWin++){
         if (iWin % 200 == 199) System.err.printf("Calculating KNN for window %d / %d\n", iWin + 1, nWins);

         // find the KVec nearest descriptors
         nnVec.clear();
         fastKNN(descriptor[iWin], mtree.getRoot(), nnVec);
         assert (nnVec.size() == KVec) : String.format("nn.size=%d  K=%d\n", nnVec.size(), KVec);

         // search the nearest descriptors for the nearest sequences
         nnSeq.clear();
         for(int i = 0; i < KVec; i++){
            NNInfo nni = nnVec.get(i);
            double dseq = calcSubseqDist(ss[iWin], ss[nni.index]);
            nnSeq.add(new NNInfo(nni.index, dseq, i2loc[nni.index]));
         }

         // save the real knn
         for(int j = 0; j < K; j++)
            knn[iWin][j] = nnSeq.get(j);
      }
      updateGUIforData();
      System.err.printf("Time for KNN computation: %s\n", Library.formatDuration(timer.time(), 2));
   }

   /** load knn data from the given file */
   protected boolean load(File f)
   {
      try{
         BufferedReader in = new BufferedReader(new FileReader(f));
         String line = in.readLine();
         StringTokenizer st = new StringTokenizer(line, " \t\r\n");
         int N = Integer.parseInt(st.nextToken());
         int K = Integer.parseInt(st.nextToken());
         int wlenMin = Integer.parseInt(st.nextToken());
         int wlenMax = Integer.parseInt(st.nextToken());
         step = Integer.parseInt(st.nextToken());
         calcWinLens(wlenMin, wlenMax);
         setupWindows();
         knn = new NNInfo[N][K];

         for(int i = 0; i < N; i++){
            st = new StringTokenizer(in.readLine(), " \t\r\n");
            for(int j = 0; j < K; j++){
               int ix = Integer.parseInt(st.nextToken());
               double d = Double.parseDouble(st.nextToken());
               knn[i][j] = new NNInfo(ix, d, i2loc[ix]);
            }
         }
         in.close();

         // update the gui
         spinLenMin.setValue(wlenMin);
         spinLenMax.setValue(wlenMax);
         lsStep.setValue(step);
         lsKNN.setValue(K);
         updateGUIforData();
         System.err.printf("Loaded KNNs: N=%d  K=%d  wlen=%d->%d  step=%d\n", N, K, wlenMin, wlenMax, step);
      } catch (Exception e){
         System.err.println(e);
         return false;
      }
      return true;
   }

   /** save knn info in the given file */
   protected boolean save(File f)
   {
      try{
         PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f)));
         out.printf("%d %d %d %d %d\n", knn.length, knn[0].length, wlens[0], wlens[wlens.length - 1], step);
         for(int i = 0; i < knn.length; i++){
            for(int j = 0; j < knn[i].length; j++)
               out.printf("%d %f ", knn[i][j].index, knn[i][j].dist);
            out.println();
         }
         out.close();
      } catch (Exception e){
         return false;
      }
      return true;
   }

   /** @return number of states for a HMM for sequences of the given length */
   protected int calcHmmLen(int wlen)
   {
      // int nStates = Math.max(3, Math.round(wlen * 0.75f));
      // int nStates = Math.max(3, wlen / 2);
      // int nStates = (int)Math.round(2.0/rBand-1);
      // if (nStates % 2 == 0) nStates++; // ensure # states is odd

      int nSegs = lsNumSegs.getValue();
      int nStates = nSegs * 2 - 1;
      nStates = Math.min(wlen * 2 - 1, nStates);
      return nStates;
   }

   /** @return left-right HMM initialized from the given window index + knn */
   protected HmmLR genHMM(int ix, String name)
   {
      ArrayList<WindowLocation> group = new ArrayList<WindowLocation>();
      group.add(i2loc[ix]);
      for(int i = 0; i < knn[ix].length; i++)
         group.add(i2loc[knn[ix][i].index]);
      return genHMM(group, name);
   }

   /** @return left-right HMM initialized from the given window index + knn */
   protected HmmLR genHMM(ArrayList<WindowLocation> group, String name)
   {
      // build initial HMM
      int nStates = calcHmmLen(wlens[0]);
      int nDims = gdata.getNumDims();
      int nSkip = 1;
      HmmLR hmm = new HmmLR(nStates, nSkip, nDims);
      hmm.setMinVar(gdata.getMinVar());
      double edur = (double)wlens[wlens.length - 1] / (nStates - 1);
      double pself = Math.max(1.0 - 1.0 / edur, 0.1);
      hmm.initTran(pself, nSkip);
      hmm.setName(name);

      // build the training set
      ArrayList<Sequence> vtrain = new ArrayList<Sequence>();
      Iterator<WindowLocation> it = group.iterator();
      while(it.hasNext()){
         WindowLocation wloc = it.next();
         vtrain.add(wloc.getSeq(gdata.tseries));
      }

      // learn HMM params
      hmm.init_segk_overlap(vtrain, HMM_OVERLAP_PERCENT);

      // double[][] tranOrig = hmm.saveTran();
      // hmm.train_bw(vtrain); // TODO BW or Viterbi or none?
      // hmm.blendTran(tranOrig, 0.5);

      if (HMM_GLOBAL_SDEV_BLEND < 1.0){
         for(int i = 0; i < nStates; i++){
            GaussianDiagonal dg = (GaussianDiagonal)hmm.getState(i);
            FeatureVec fvCurSDev = gdata.getGlobalVar().sqrt();
            FeatureVec fv = fvCurSDev._mul(1.0 - HMM_GLOBAL_SDEV_BLEND);
            FeatureVec fvNewSDev = dg.getVar().sqrt()._mul(HMM_GLOBAL_SDEV_BLEND)._add(fv)._sqr();
            fvNewSDev._max(fvCurSDev); // should not decrease estimated variance
            dg.setVar(fvNewSDev);
         }
      }

      return hmm;
   }

   protected boolean reEstBG(ArrayList<MarkupSet> labels)
   {
      int nSeqs = gdata.getNumSeqs();
      int nDims = gdata.getNumDims();
      // GaussianDynDiag gdd = new GaussianDynDiag(nDims);

      // collect all occurrences of background
      Sequence all = new Sequence();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         MarkupSet marks = labels.get(iSeq);
         Sequence seq = gdata.tseries.get(iSeq);
         for(TimeMarker tm : marks.getList()){
            if (!BG.equals(tm.getTag())) continue;
            int a = tm.getStartIndex();
            int b = tm.getStopIndex();
            all.append(seq.subseq(a, b));
            // for(int i = a; i < b; i++)
            // gdd.add(seq.get(i), false);
         }
      }

      // GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      // gdd.update();
      // dg.setMean(gdd.getMean());
      // dg.setVar(gdd.getVar());
      hmmBG.setState(0, learnMixtureModel(NMIX, all));

      // now that we've reestimated, rerun the cont rec system
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      hmms.addAll(motifHMMs);
      crrBase = HTK.contRec(hsi, hmms, gdata.tseries);
      if (crrBase == null){
         System.err.println("Error: HTK / HVite failed after re-estimation.");
         return false;
      }
      motifLabels = crrBase.markupSets;
      highlightMotif(lsMotif.getValue() - 1, motifLabels);

      // display the updated HMM
      if (lsMotif.getValue() == 0) mviz.setHMM(hmmBG);

      return true;
   }

   protected boolean reEstMotifs(ArrayList<MarkupSet> labels)
   {
      int nMotifs = motifHMMs.size();
      int nSeqs = gdata.getNumSeqs();
      System.err.printf("#motifs=%d  #seq=%d\n", nMotifs, nSeqs);
      for(int iMotif = 0; iMotif < nMotifs; iMotif++){
         AbstractHMM hmm = motifHMMs.get(iMotif);
         String name = hmm.getName();
         ArrayList<Sequence> vtrain = new ArrayList<Sequence>();

         // collect all occurrences of this motif
         for(int iSeq = 0; iSeq < nSeqs; iSeq++){
            MarkupSet marks = labels.get(iSeq);
            Sequence seq = gdata.tseries.get(iSeq);
            for(TimeMarker tm : marks.getList()){
               if (name.equals(tm.getTag()))
                  vtrain.add(seq.subseq(tm.getStartIndex(), tm.getStopIndex(), iSeq));
            }
         }

         if (vtrain.isEmpty()) continue;
         System.err.printf("Training %s: %d occurrences... ", name, vtrain.size());
         TimerMS timer = new TimerMS();
         hmm.train_bw(vtrain);
         ((HmmLR)hmm).removeUnusedStates();
         System.err.printf("(%dms).\n", timer.time());
      }

      // display the updated HMM
      int iViz = lsMotif.getValue() - 1;
      if (iViz >= 0) mviz.setHMM(motifHMMs.get(iViz));

      // now that we've reestimated, rerun the cont rec system
      return runContRec(-1, true) != null;
   }

   /**
    * Run the continuous recognition system
    * 
    * @param nMotifs number of motifs to use (-1 for all)
    * @param bSaveResults true => copy results into global variables
    * @return continuous recognition results or null on error
    */
   protected ContRecRet runContRec(int nMotifs, boolean bSaveResults)
   {
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      if (nMotifs < 0 || nMotifs >= motifHMMs.size()) hmms.addAll(motifHMMs);
      else{
         for(int i = 0; i < nMotifs; i++)
            hmms.add(motifHMMs.get(i));
      }
      ContRecRet crr = HTK.contRec(hsi, hmms, gdata.tseries);
      if (crr == null){
         System.err.println("Error: HTK / HVite failed after re-estimation.");
         return null;
      }

      if (bSaveResults){
         crrBase = crr;
         motifLabels = crrBase.markupSets;
         highlightMotif(lsMotif.getValue() - 1, motifLabels);
         HashMap<String, ContRecRet.MotifInfo> map = crrBase.getMotifInfo();
         listInfoGain.clear();
         for(AbstractHMM hmm : motifHMMs){
            ContRecRet.MotifInfo mi = map.get(hmm.getName());
            // listInfoGain.add(-mi.loglik / mi.getNumLocs());
            listInfoGain.add(calcDensity(mi));
         }
      }
      return crr;
   }

   protected double calcDensity(ContRecRet.MotifInfo mi)
   {
      if (mi == null) return 0;
      MetricSeq.LengthPrep lprep = metseq.getLengthPrep();
      metseq.setLengthPrep(MetricSeq.LengthPrep.extend);
      if (mi.getNumLocs() < 2) return 0; // can't calc density with 1 point
      Sequence[] seqs = new Sequence[mi.getNumLocs()];
      for(int i = 0; i < seqs.length; i++)
         seqs[i] = mi.wlocs.get(i).getSeq(gdata.tseries);
      double radius = 0;
      for(int iter = 0; iter < 10; iter++){
         int iBase = Library.random(seqs.length);

         // find farthest from base
         int iFar = -1;
         double far = 0;
         for(int i = 0; i < seqs.length; i++){
            if (i == iBase) continue;
            double r = calcSubseqDist(seqs[iBase], seqs[i]);
            if (r > far){
               far = r;
               iFar = i;
            }
         }

         // update base and find farthest
         iBase = iFar;
         iFar = -1;
         far = 0;
         for(int i = 0; i < seqs.length; i++){
            if (i == iBase) continue;
            double r = calcSubseqDist(seqs[iBase], seqs[i]);
            if (r > far){
               far = r;
               iFar = i;
            }
         }

         if (far > radius) radius = far;
      }
      metseq.setLengthPrep(lprep);
      return (double)seqs.length / radius;
   }

   protected void highlightSeed(int iSeed)
   {
      int nSeqs = gdata.getNumSeqs();

      // clear the highlights
      for(int i = 0; i < nSeqs; i++)
         crbGT[i].clearHighlight();

      // highlight a real local optima
      if (iSeed >= 0){
         ArrayList<WindowLocation> group = seeds.get(iSeed);
         for(WindowLocation wloc : group){
            Range r = wloc.getRange();
            crbGT[wloc.iSeries].addHighlight(r);
         }
      }
   }

   protected void highlightMotif(int iMotif, ArrayList<MarkupSet> labels)
   {
      int nSeqs = gdata.getNumSeqs();

      // clear the highlights
      for(int i = 0; i < nSeqs; i++)
         crbGT[i].clearHighlight();

      // highlight a real motif
      if (motifHMMs.isEmpty()) return;
      String sMotifName = (iMotif >= 0) ? motifHMMs.get(iMotif).getName() : hmmBG.getName();
      assert (nSeqs == labels.size()) : String.format("mismatch: nSeqs=%d  labels.size=%d", nSeqs, labels
            .size());
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         for(TimeMarker tm : labels.get(iSeq).getList()){
            if (tm.getTag() == null) continue;
            if (!tm.getTag().equals(sMotifName)) continue;
            Range r = new Range(tm.getStartIndex(), tm.getStopIndex() - 1);
            crbGT[iSeq].addHighlight(r);
         }
      }
   }

   /** @return index of motif given a name */
   protected int getMotifIndex(String name)
   {
      if (name.equals(BG)) return -1;

      // try searching for the name directly
      // TODO hash table
      int nMotifs = motifHMMs.size();
      for(int i = 0; i < nMotifs; i++)
         if (motifHMMs.get(i).getName().equals(name)) return i;

      // if we can't find an HMM, try to extract a number
      try{
         Matcher m = rexIndex.matcher(name);
         if (m.matches()){ return Integer.parseInt(m.group(1)); }
      } catch (NumberFormatException e){}
      assert false : String.format("Can't find motif for given name (%s)", name);
      return -1;
   }

   protected void scoreMotifs(ArrayList<MarkupSet> labels)
   {
      // setup the cont rec info
      ContRecInfo cri = new ContRecInfo(gdata.getNumClasses());
      Iterator<String> itClass = gdata.labData.keySet().iterator();
      for(int iClass = 0; iClass < gdata.labData.size(); iClass++){
         ArrayList<Sequence> subs = gdata.labData.get(itClass.next());
         cri.nLabeledWords += subs.size();
         for(Sequence sub : subs)
            cri.nLabeledFrames += sub.length();
      }

      // add locations to cri
      int nspots = 0;
      int nSeqs = gdata.getNumSeqs();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         for(TimeMarker tm : labels.get(iSeq).getList()){
            if (tm.getTag() == null || tm.getTag().equals(BG)) continue;
            double score = (tm.getMeta() == null ? 0 : (Double)tm.getMeta());
            int iClass = getMotifIndex(tm.getTag());
            assert (iClass >= 0) : String.format("iClass=%d  name=%s", iClass, tm.getTag());
            cri.add(new WordSpot(iSeq, tm.getStartIndex(), (int)tm.length(), score, iClass));
            nspots++;
         }
      }

      // evaluate the results
      // int[][] confm = cri.scoreWordSpotWordsFast(gdata.labData);
      int[][] omap = cri.buildOverlapMap(gdata.labData, gdata.getNumSeqs());

      /*
       * System.err.printf("Overlap map (%d x %d) [RAW]:\n", omap.length, omap[0].length); for(int i = 0; i <
       * omap.length; i++){ for(int j = 0; j < omap[i].length; j++) System.err.printf("%3d ", omap[i][j]);
       * System.err.println(); }
       */

      cri.cleanOMap(omap, CLEAN_FRAC);
      int[] nFoundLabels = new int[omap.length];
      for(WordSpot spot : cri.getSpots())
         nFoundLabels[spot.iClass]++;

      /*
       * System.err.printf("Overlap map (%d x %d) [CLEAN]:\n", omap.length, omap[0].length); for(int i = 0; i <
       * omap.length; i++){ for(int j = 0; j < omap[i].length; j++) System.err.printf("%3d ", omap[i][j]);
       * System.err.println(); }
       */

      // calc best class permutation and remap spots
      System.err.print("Calculating best cluster <-> class mapping... ");

      /*
       * System.err.print("\nFound: "); for(int i = 0; i < nFoundLabels.length; i++) System.err.printf("%3d ",
       * nFoundLabels[i]); System.err.println(); System.err.print(" True: "); for(int i = 0; i <
       * gdata.nExamples.length; i++) System.err.printf("%3d ", gdata.nExamples[i]); System.err.println();
       */

      // short[] map = SupTest.calcBestMappingFromConfMatrix(confm, bOptAcc);
      boolean bOptAcc = btrOptAcc.isSelected();
      short[] map = ContRecInfo.calcBestMappingFromOMap(omap, nFoundLabels, gdata.nExamples, bOptAcc);
      cri.remapSpots(map, gdata.labData.size());
      System.err.println("done.");

      /*
       * for(int i = 0; i < map.length; i++){ if (map[i] < 0) continue; System.err.printf("%s -> %s\n",
       * motifHMMs.get(i).getName(), gdata.classes[map[i]]); }
       */

      // recalc stats with permuted classes
      System.err.print("Recomputing performance stats... ");
      TimerMS timer = new TimerMS();
      int[][] confm = cri.scoreWordSpotWordsFull(gdata.labData);
      cri.scoreWordSpotFrames(gdata.tseries, gdata.labData);
      System.err.printf("done (%dms).\n", timer.time());

      SupTest.dumpConfMatrix(confm, "Confusion Matrix", gdata.classes, gdata.nExamples, false);

      // cri.dump();

      assert (cri.isValid(true)) : "invalid continuous recognition stats";
      // System.err.println(cri.getResultsString(gdata.classes));
   }

   protected void showSilenceGraph(int wlen)
   {
      // extract windows
      int nSeqs = gdata.getNumSeqs();
      int nWins = 0;
      for(int i = 0; i < nSeqs; i++)
         nWins += Library.getNumSlidingWindowSites(gdata.seqLens[i], wlen, 1);
      Sequence var = new Sequence("Window Variance");
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = gdata.tseries.get(iSeq);
         int T = gdata.seqLens[iSeq] - wlen;
         for(int t = 0; t <= T; t++)            
            var.add(seq.getVar(t, t+wlen));
      }
      assert (var.length() == nWins);

      // display results
      int nDims = gdata.getNumDims();
      //int nCols = (int)Math.floor(Math.sqrt(nDims));
      //int nRows = (int)Math.ceil((double)nDims / nCols);
      JPanel p = new JPanel(new GridFillLayout(nDims, 2, 2, 2));
      for(int d=0; d<nDims; d++){
         double[] x = var.extractDim(d);
         double xmax = StatUtils.max(x);
         MyDoubleList list = new MyDoubleList();
         //for(int i=0; i<x.length; i++)
         //   if (x[i] < xmax/5) list.add(x[i]);
         //x = list.toArray();
         SeqStatsView ssv = new SeqStatsView(x, 100);
         ssv.setTitle(String.format("Variance (d=%d)", d+1));
         p.add(ssv);
         
         ScatterPlot plot = new ScatterPlot();         
         Arrays.sort(x);
         plot.addData(new Sequence("Sorted Variance", x), Color.blue, ScatterPlot.DEF_STROKE);
         p.add(plot);
      }
      Library.popup("Window Variance", p, 800,600);
   }
   
   protected void markSilence(int wlen)
   {
      // extract windows
      int nSeqs = gdata.getNumSeqs();
      int nWins = 0;
      for(int i = 0; i < nSeqs; i++)
         nWins += Library.getNumSlidingWindowSites(gdata.seqLens[i], wlen, 1);
      Sequence var = new Sequence("Window Variance");
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = gdata.tseries.get(iSeq);
         int T = gdata.seqLens[iSeq] - wlen;
         for(int t = 0; t <= T; t++)            
            var.add(seq.getVar(t, t+wlen));
      }
      assert (var.length() == nWins);
      
      FeatureVec varMin = var.getMin();
      FeatureVec varMax = var.getMax();
      FeatureVec fvNoise = varMax.sub(varMin)._div(4.0)._add(varMin);
      
      // clear the highlights
      for(int i = 0; i < nSeqs; i++)
         crbGT[i].clearHighlight();

      // highlight a real motif
      int ix = 0;
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         int T = gdata.seqLens[iSeq] - wlen;
         for(int t = 0; t <= T; t++){
            if (var.get(ix++).leqThan(fvNoise))
               crbGT[iSeq].addHighlight(new Range(t, t));
         }
      }
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();
      if (src == btFastKNN){
         int step = lsStep.getValue();
         int k = lsKNN.getValue();
         calcFastKNN(step, k);
      }
      else if (src == btFullKNN){
         int step = lsStep.getValue();
         int k = lsKNN.getValue();
         calcFullKNN(step, k);
      }
      else if (src == btLoad){
         JFileChooser fc = Library.buildFileChooser("txt", "KNN Data (*.txt)");
         if (fc.showOpenDialog(frame) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            if (!load(f)){
               System.err.println("Error: failed to load KNN data in file:");
               System.err.println(" " + Library.getCanonical(f));
               JOptionPane.showMessageDialog(frame, "Failed to load KNN data", "Load Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
         }
      }
      else if (src == btSave){
         JFileChooser fc = Library.buildFileChooser("txt", "KNN Data (*.txt)");
         if (fc.showSaveDialog(frame) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            if (!f.getAbsolutePath().endsWith(".txt")) f = new File(f.getAbsolutePath() + ".txt");
            if (!save(f)){
               System.err.println("Error: failed to save KNN data in file:");
               System.err.println(" " + Library.getCanonical(f));
               JOptionPane.showMessageDialog(frame, "Failed to save KNN data", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
            else JOptionPane.showMessageDialog(frame, "Saved KNN data", "Save Successful",
                  JOptionPane.INFORMATION_MESSAGE);
         }
      }
      else if (src == btCheckOverlap){
         for(int x = 0; x < knn.length; x++){
            NNInfo[] a = knn[x];
            for(int i = 0; i < a.length; i++){
               if (i2loc[x].overlaps(a[i].wloc)){
                  System.err.printf("Error: query point %d overlaps neighbor %d (%d)! %s  vs.  %s\n", x + 1,
                        i + 1, a[i].index, i2loc[x], a[i].wloc);
                  return;
               }
               for(int j = i + 1; j < a.length; j++)
                  if (a[i].hasConflict(a[j])){
                     System.err.printf("Error: neighbors in conflict! (%d,%d)\n %s  vs.  %s\n", a[i].index,
                           a[j].index, a[i].wloc, a[j].wloc);
                     return;
                  }
            }
         }
         System.err.println("Nearest neighbors are overlap free!");
      }
      else if (src == btSilenceGraph){
         showSilenceGraph((Integer)spinLenMin.getValue());
      }
      else if (src == btMarkSilence){
         markSilence((Integer)spinLenMin.getValue());
      }
      else if (src == btLocalOpt){
         btLocalOpt.setEnabled(false);
         btcGrow.setEnabled(false);
         calcLocalOptsHmms();
         lsNumSegs.setEnabled(false);
         btFindNext.setEnabled(true);
         lsNumFind.setEnabled(true);
         btcExact.setEnabled(true);
         lsSeed.setEnabled(true);
         lsSeed.setValues(0, 0, seeds.size());
      }
      else if (src == btFindNext){
         TimerMS timer = new TimerMS();
         boolean bExact = btcExact.isSelected();
         int nFind = lsNumFind.getValue();
         if (nFind == lsNumFind.getMaximum()) nFind = -1;
         int nPre = motifHMMs.size();
         while(true){
            int nFound = motifHMMs.size() - nPre;
            if (nFound >= nFind) break;
            int nFindLeft = (bExact ? nFind - nFound : -1);
            if (!findNext(nFindLeft)){
               System.err.printf("Error: failed to find next motif\n");
               break;
            }

            // rebuild seed list
            System.err.printf("# seeds before: %d\n", allHMMs.size());
            for(int i = 0; i < allHMMs.size();){
               if (allHMMs.get(i) == null) allHMMs.remove(i);
               else i++;
            }
            System.err.printf("# seeds after: %d\n", allHMMs.size());

            lsMotif.setText(String.format("Motif: %%d / %d", motifHMMs.size()));
            lsMotif.setValues(motifHMMs.size(), 0, motifHMMs.size());
            lsMotif.setEnabled(true);
            btScoreMotifs.setEnabled(true);
            btrOptAcc.setEnabled(true);
            btrOptRecall.setEnabled(true);
            btReEst.setEnabled(true);
            btReEstBG.setEnabled(true);
            if (allHMMs.size() < 1) break;
         }
         System.err.printf("Total time for find-next: %s\n", Library.formatDuration(timer.time()));
      }
      else if (src == btScoreMotifs){
         int nMotifs = motifHMMs.size();
         int nScore = sgc.getDisplayIndex() + 1;
         if (nMotifs == nScore) scoreMotifs(motifLabels);
         else{
            ContRecRet crr = runContRec(nScore, false);
            if (crr == null) return;
            scoreMotifs(crr.markupSets);
         }
      }
      else if (src == btReEst){
         reEstMotifs(motifLabels);
      }
      else if (src == btReEstBG){
         reEstBG(motifLabels);
      }
      else if (src == btReset){
         updateGUIforData();
      }
      else if (src == btDiscoverView){
         clView.show(mainp, sDiscoverView);
      }
      else if (src == btReviewView){
         clView.show(mainp, sReviewView);
      }
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();
      if (src == lsMotif){
         int iMotif = lsMotif.getValue() - 1;
         highlightMotif(iMotif, motifLabels);
         if (iMotif >= 0) mviz.setHMM(motifHMMs.get(iMotif));
         else mviz.setHMM(hmmBG);
      }
      else if (src == lsNumFind){
         int nf = lsNumFind.getValue();
         if (nf == lsNumFind.getMaximum()) btFindNext.setText("Find All Motifs");
         else if (nf == 1) btFindNext.setText("Find Next Motif");
         else btFindNext.setText(String.format("Find Next %d Motifs", nf));
      }
      else if (src == spinLenMin){
         int wmin = (Integer)spinLenMin.getValue();
         int wmax = (Integer)spinLenMax.getValue();
         ((SpinnerNumberModel)spinLenMax.getModel()).setMinimum(wmin);
         calcWinLens(wmin, wmax);
         lsStep.setValues(lsStep.getValue(), 1, wmin);
      }
      else if (src == spinLenMax){
         int wmin = (Integer)spinLenMin.getValue();
         int wmax = (Integer)spinLenMax.getValue();
         ((SpinnerNumberModel)spinLenMin.getModel()).setMaximum(wmax);
         calcWinLens(wmin, wmax);
         lsStep.setValues(lsStep.getValue(), 1, wmin);
      }
      else if (src == lsSeed){
         int iSeed = lsSeed.getValue() - 1;
         highlightSeed(iSeed);
      }
   }

   /** @return number of motifs to include based on analyzing given number of motifs */
   protected int calcStop(int nDisp)
   {
      if (nDisp == 1) return 1;

      // create sequence from density
      Sequence data = new Sequence();
      data.add(new FeatureVec(1, listInfoGain.get(1)));
      for(int i = 1; i < nDisp; i++)
         data.add(new FeatureVec(1, listInfoGain.get(i)));

      // smooth values
      TransformSmooth smoother = new TransformSmooth(5);
      Sequence smooth = smoother.transform(data);
      for(int i = 0; i < 4; i++)
         smooth = smoother.transform(smooth);

      // calc deriv
      double[] orig = data.extractDim(0);
      double[] a = smooth.extractDim(0);
      double[] da = new double[nDisp];
      da[0] = 0;
      for(int i = 1; i < nDisp; i++)
         da[i] = a[i] - a[i - 1];

      // for(int i = 0; i < nDisp; i++)
      // System.err.printf(" %f %f %f\n", orig[i], a[i], da[i]);

      // TODO find max neg region, even if there's positive

      // find max negative stretch
      double[] dsum = new double[nDisp];
      int iMin = nDisp - 1;
      double vMin = 0;
      for(int i = 1; i < nDisp; i++){
         if (da[i] < 0) dsum[i] = dsum[i - 1] + da[i];
         else dsum[i] = 0;
         if (dsum[i] < vMin){
            vMin = dsum[i];
            iMin = i;
         }
      }
      // find start of this stretch
      int jMin = iMin - 1;
      while(jMin > 0){
         if (dsum[jMin] >= 0){
            jMin++;
            break;
         }
         jMin--;
      }
      jMin = Math.max(jMin, 0);
      if (jMin < iMin && dsum[jMin] >= 0) jMin++;
      jMin = Math.min(jMin, iMin);

      // find weighted mean of jMin -> iMin range
      double x = 0, y = 0;
      for(int i = jMin; i <= iMin; i++){
         assert (i == nDisp - 1 || dsum[i] < 0) : String.format("nDisp=%d  i=%d  dsum[i]=%f", nDisp, i,
               dsum[i]);
         x += da[i] * i;
         y += da[i];
      }
      double f = (y == 0 ? nDisp - 1 : x / y);
      int fi = (int)Math.ceil(f);

      // int iFinal = fi + (iMin-fi+1)/4;
      int iFinal = Math.min(fi, nDisp - 1);
      // System.err.printf("nDisp=%d iMin=%d jMin=%d f=%.2f (%d) iFinal=%d (=> %d motifs)\n", nDisp, iMin,
      // jMin, f, fi, iFinal, iFinal+1);
      return iFinal + 1;
   }

   public void vizChanged(ScoreGraphComplex src, int nDisp)
   {
      if (src == sgc){
         nDisp++;
         btScoreMotifs.setEnabled(nDisp >= 1);
         btrOptAcc.setEnabled(nDisp >= 1);
         btrOptRecall.setEnabled(nDisp >= 1);
         btScoreMotifs.setText(String.format("Score Motifs (%d)", nDisp));
         System.err.printf("Auto-stop: %d motifs (nDisp=%d)\n", calcStop(nDisp), nDisp);
      }
   }

   public static void main(String args[]) throws Exception
   {
      System.err.printf("Java Version: %s   VM: %s (%s)\n", System.getProperty("java.version"), System
            .getProperty("java.vm.version"), System.getProperty("java.vm.name"));
      System.err.printf("Memory (free/total/max): %dk / %dk / %dk\n", rt.freeMemory() / 1024, rt
            .totalMemory() / 1024, rt.maxMemory() / 1024);

      final String sDefFile = "d:/classes/research/kdm/data/ae-dim13/ae.def";
      
      //final String sDefFile = "/home/dminn/research/kdm/data/ah-data/exercise/sub8/accgyr.def";
      //final String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/dim13/ae.def";
      // final String sDefFile =
      // "/home/dminn/research/kdm/data/synth/chars/easy/easy-zgwxfpmik-serif-noshot/easy.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/ae.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/asl-joseph/asl16D.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/keogh/shuttle.def";

      final KDD2007 kdd = new KDD2007();
      SwingUtilities.invokeAndWait(new Runnable() {
         public void run()
         {
            if (!kdd.setup(sDefFile)){
               System.err.println("KDD setup failed.");
               System.exit(1);
            }
         }
      });

      final BlockDlg dlg = new BlockDlg(kdd.frame, true);
      dlg.display("Loading default data", new Thread() {
         public void run()
         {
            /*if (!kdd.load(new File("/home/dminn/exercise-kdd-12-20-norm.txt"))){
               // if (!kdd.load(new File("/home/dminn/tidigits-kdd-21-35-norm.txt"))){
               // if (!kdd.load(new File("/home/dminn/easy9-serif-noshot-8-26-norm.txt"))){
               System.err.println("Error: failed to load default data!");
               System.exit(1);
            }*/
            dlg.close();
         }
      });
   }
}
