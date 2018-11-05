package mdisc.aaai2007;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import kdm.data.*;
import kdm.data.transform.TransformSmooth;
import kdm.gui.*;
import kdm.io.*;
import kdm.metrics.*;
import kdm.mlpr.dataTree.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.tools.*;
import kdm.util.*;
import mdisc.VizTool.*;
import mdisc.kdd2007.*;
import kdm.mlpr.htk.*;
import java.util.regex.*;

/** visualization for KNN-based density mode estimation */
public class ExploreKNN implements ActionListener, ChangeListener, ScoreGraphListener
{
   public static final double HMM_OVERLAP_PERCENT = 0.5;
   public static final double CLEAN_FRAC = 0.25;
   public static final boolean bIncBg = false;
   public static final boolean bIncFirstMotif = true;
   public static final String sTmpDir = System.getProperty("java.io.tmpdir");
   public static final String BG = "bg";
   public static final Pattern rexIndex = Pattern.compile("^.*(\\d+).*$");
   public static final MetricSeq metseq = new DTWSimple(0.1, MetricSeq.LengthPrep.none);
   public static final double MODSEL_lambda = 1.0;
   public static final double DEN_MERGE_FACTOR = 1.1; // TODO right value for this constant?

   public PrintWriter out;

   protected GlobalData gdata;
   protected JFrame frame;
   protected ColoredRangeBar crbGT[];
   protected Color[] colors;
   protected VertChooseContainer vcc;
   protected ScoreGraphComplex sgc;

   protected LabeledSlider lsWinLen, lsStep, lsKNN, lsLocalOpt, lsMotif, lsNumFind;
   protected JButton btKNN, btFastKNN, btSave, btLoad, btCheckOverlap, btLocalOpt;
   protected JButton btFindNext, btScoreMotifs, btReEst, btReEstBG, btNbrGraph;
   protected JButton btSaveHMMs;
   protected JRadioButton btrOptAcc, btrOptRecall;

   protected int wlen, step;
   protected NNInfo[][] knn;
   protected WindowLocation[] i2loc;
   protected Sequence[] ss;
   protected int[] seeds;
   protected AbstractHMM[] allHMMs;
   protected AbstractHMM hmmBG;
   protected ArrayList<AbstractHMM> motifHMMs = new ArrayList<AbstractHMM>();
   protected ArrayList<MarkupSet> motifLabels = new ArrayList<MarkupSet>();
   protected ContRecRet crrBase = null;
   protected MotifViz mviz;

   protected MyDoubleList listInfoGain = new MyDoubleList();
   protected MyDoubleList listDataLL = new MyDoubleList();
   protected MyDoubleList listBIC = new MyDoubleList();
   protected MyDoubleList listThreshMDL = new MyDoubleList();

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
      if (data == null) System.err.printf("Failed to load labeled data\n (%s)\n", sDefFile);

      // setup the global data structure
      gdata.labData = data;
      gdata.updateClassStrings();
      gdata.tseries = LabeledDataLoader.tseries;
      gdata.calcLengthStats(gdata.tseries);
      System.err.printf("Loaded %d time series (len: total=%d, avg=%d).\n", gdata.tseries.size(),
            gdata.nTotalSeriesLength, gdata.nAvgSeriesLength);
      gdata.sentences = LabeledDataLoader.sentences;
      Gaussian1D[] gmSeries = SupTest.calcGauss(gdata.tseries.toArray(new Sequence[0]));
      gdata.setInitVar(SupTest.initv);
      gdata.setMinVar(SupTest.minv);
      int nDims = gdata.getNumDims();
      FeatureVec gmean = new FeatureVec(nDims);      
      FeatureVec gvar = new FeatureVec(nDims);
      for(int i = 0; i < nDims; i++){
         gmean.set(i, gmSeries[i].getMean());
         gvar.set(i, gmSeries[i].getVar());
      }
      gdata.setGlobalMean(gmean);      
      gdata.setGlobalVar(gvar);
      FeatureVec gmin = gdata.tseries.get(0).getMin();
      FeatureVec gmax = gdata.tseries.get(0).getMax();
      int nSeqs = gdata.getNumSeqs();
      for(int i=1; i<nSeqs; i++){
         gmin._min(gdata.tseries.get(i).getMin());
         gmax._max(gdata.tseries.get(i).getMax());
      }
      System.err.printf("min: %s\nmax: %s\n", gmin, gmax);
      gdata.setGlobalMin(gmin);
      gdata.setGlobalMax(gmax);
      LabeledDataLoader.dumpDataSummary(gdata.labData);

      // build the background model
      createBgModel();

      // create the window
      System.err.print("Creating visualization... ");
      timer.reset();
      frame = new JFrame(String.format("ExploreKNN - %s", Library.getTitle(sDefFile)));
      gdata.frame = frame;
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setSize(1200, 700);
      Library.centerWin(frame, null);
      frame.setLocation(frame.getLocation().x, 0);
      JPanel mainp = new JPanel(new BorderLayout());
      frame.setContentPane(mainp);

      mainp.add(buildGUI(), BorderLayout.CENTER);

      System.err.printf("done (%dms).\n", timer.time());
      frame.setVisible(true);
      return true;
   }

   protected void createBgModel()
   {
      int nDims = gdata.getNumDims();
      hmmBG = new HmmLR(1, nDims);
      ((HmmLR)hmmBG).setMinVar(gdata.getMinVar());
      hmmBG.setName(BG);
      hmmBG.setPiLeave(0, Math.log(0.2));
      GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      dg.setMean(gdata.getGlobalMean());
      dg.setVar(gdata.getGlobalVar());
   }

   protected JPanel buildGUI()
   {
      JPanel mainp = new JPanel(new BorderLayout());
      JPanel leftp = new JPanel(new BorderLayout());

      // create the center view
      colors = Library.generateColors(gdata.getNumClasses());
      int nSeqs = gdata.getNumSeqs();
      crbGT = new ColoredRangeBar[nSeqs];
      JPanel seqsp = new JPanel(new VerticalLayout());
      seqsp.setBackground(Color.darkGray);
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

         // add labeled ranges to the bar
         Iterator<String> itClass = gdata.labData.keySet().iterator();
         for(int iClass = 0; iClass < gdata.labData.size(); iClass++){
            String sClass = itClass.next();
            ArrayList<Sequence> subs = gdata.labData.get(sClass);
            for(Sequence sub : subs){
               if (sub.getParentIndex() != iSeq) continue;
               Range r = new Range(sub.getParentOffset(), sub.getParentOffset() + sub.length() - 1);
               r.payload = new Pair(sClass, colors[iClass].darker());
               crbGT[iSeq].add(r);
            }
         }
      }
      leftp.add(new JScrollPane(seqsp), BorderLayout.CENTER);

      // create bottom view
      sgc = new ScoreGraphComplex();
      sgc.addScoreGraphListener(this);
      sgc.addCurve("Info Gain", listInfoGain);
      // sgc.addCurve("Log-likelihood", listDataLL);
      // sgc.addCurve("BIC", listBIC);
      // sgc.addCurve("MDL", listThreshMDL);
      leftp.add(sgc, BorderLayout.SOUTH);

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
      vcc.addPane("KNN", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      lsWinLen = new LabeledSlider("Window Length: %d", 4, 250);
      rightp.add(lsWinLen);
      lsStep = new LabeledSlider("Step Size: %d", 1, 125);
      rightp.add(lsStep);
      lsKNN = new LabeledSlider("K: %d", 1, 20);
      rightp.add(lsKNN);
      btKNN = new JButton("Calc KNNs");
      btKNN.addActionListener(this);
      rightp.add(btKNN);
      btFastKNN = new JButton("Calc Fast KNNs");
      btFastKNN.addActionListener(this);
      rightp.add(btFastKNN);
      JPanel p = new JPanel();
      btLoad = new JButton("Load KNN");
      btLoad.addActionListener(this);
      p.add(btLoad);
      btSave = new JButton("Save KNN");
      btSave.addActionListener(this);
      btSave.setEnabled(false);
      p.add(btSave);
      rightp.add(p);
      btCheckOverlap = new JButton("Check for Overlapping Windows");
      btCheckOverlap.addActionListener(this);
      btCheckOverlap.setEnabled(false);
      rightp.add(btCheckOverlap);

      // create the control pane for motif discovery
      rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      vcc.addPane("Motif Discovery", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      btLocalOpt = new JButton("Detect Local Optima");
      btLocalOpt.addActionListener(this);
      btLocalOpt.setEnabled(false);
      rightp.add(btLocalOpt);
      lsLocalOpt = new LabeledSlider("Local Optima", 0, 0);
      lsLocalOpt.addChangeListener(this);
      lsLocalOpt.setEnabled(false);
      rightp.add(lsLocalOpt);

      btNbrGraph = new JButton("Neighborhood");
      btNbrGraph.addActionListener(this);
      btNbrGraph.setEnabled(false);
      rightp.add(btNbrGraph);

      lsNumFind = new LabeledSlider("# Motifs to Find: %d", 1, 31);
      lsNumFind.addChangeListener(this);
      lsNumFind.setEnabled(false);
      rightp.add(lsNumFind);

      btFindNext = new JButton("Find Next Motif");
      btFindNext.addActionListener(this);
      btFindNext.setEnabled(false);
      rightp.add(btFindNext);

      lsMotif = new LabeledSlider("No Motifs", 0, 0);
      lsMotif.addChangeListener(this);
      lsMotif.setEnabled(false);
      rightp.add(lsMotif);

      p = new JPanel();
      btReEst = new JButton("ReEst Motifs");
      btReEst.addActionListener(this);
      btReEst.setEnabled(false);
      p.add(btReEst);
      btReEstBG = new JButton("ReEst BG");
      btReEstBG.addActionListener(this);
      btReEstBG.setEnabled(false);
      p.add(btReEstBG);
      rightp.add(p);

      btScoreMotifs = new JButton("Score Motifs");
      btScoreMotifs.addActionListener(this);
      btScoreMotifs.setEnabled(false);
      rightp.add(btScoreMotifs);

      p = new JPanel();
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
      rightp.add(p);

      btSaveHMMs = new JButton("Save HMMs");
      btSaveHMMs.addActionListener(this);
      btSaveHMMs.setEnabled(false);
      rightp.add(btSaveHMMs);
      
      mviz = new MotifViz(null, gdata.getGlobalMin(), gdata.getGlobalMax());
      mviz.setBorder(BorderFactory.createLoweredBevelBorder());
      rightp.add(new Constrainer(mviz, 160, 160));

      // set initial values for various components
      vcc.show(1);
      lsWinLen.setValue(15);
      lsStep.setValue(2);
      lsKNN.setValue(5);
      lsNumFind.setValue(1);

      return mainp;
   }

   /**
    * search for pairs of clusters (point + knn) that are too close (lambda*d(knn)); prune the lower density
    * cluster
    * 
    * @param a list of indices of cluster centers; -1 implies it's been pruned
    * @param K number of nearest neighbors
    * @param lambda multiplier for distance comparison (d(i,j) < d(i, knn(i))*lambda => prune)
    */
   protected void pruneClusterPairs(int[] a, int K, double lambda)
   {
      // check inter-cluster dist; if close to knn dist, prune lower density cluster
      WindowLocation wloc = new WindowLocation(0, wlen);
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue; // already removed this index
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // calc min dist from knn(a[i]) to knn(a[j])
            double dmin = Library.INF;
            for(int ik = 0; ik < K; ik++)
               for(int jk = 0; jk < K; jk++){
                  double d = metseq.dist(ss[knn[a[i]][ik].index], wloc, ss[knn[a[j]][jk].index], wloc);
                  if (d < dmin) dmin = d;
               }

            // min over a[i] to knn(a[j]) and a[j] to knn(a[i])
            for(int ik = 0; ik < K; ik++){
               double d = metseq.dist(ss[a[j]], wloc, ss[knn[a[i]][ik].index], wloc);
               if (d < dmin) dmin = d;
               d = metseq.dist(ss[a[i]], wloc, ss[knn[a[j]][ik].index], wloc);
               if (d < dmin) dmin = d;
            }

            // min over a[i] to a[j]
            double d = metseq.dist(ss[a[i]], wloc, ss[a[j]], wloc);
            if (d < dmin) dmin = d;

            dmin /= lambda;
            if (dmin < knn[a[i]][K - 1].dist || dmin < knn[a[j]][K - 1].dist){
               // prune lower density set
               if (knn[a[i]][K - 1].dist < knn[a[j]][K - 1].dist) a[j] = -1;
               else{
                  a[i] = -1;
                  break;
               }
            }
         }
      }
   }

   /**
    * search for connected components of clusters (point + knn) that are too close (lambda*d(knn)); prune all
    * but the highest density cluster
    * 
    * @param a list of indices of cluster centers; -1 implies it's been pruned
    * @param K number of nearest neighbors
    * @param lambda multiplier for distance comparison (d(i,j) < d(i, knn(i))*lambda => prune)
    */
   protected void pruneConnectedComponents(int[] a, int K, double lambda)
   {
      AList2<Integer> adj = new AList2<Integer>();

      // check inter-cluster dist; if close to knn dist, prune lower density cluster
      WindowLocation wloc = new WindowLocation(0, wlen);
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue;
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // calc min dist from a[i] to a[j]
            double dmin = Library.INF;
            for(int ik = 0; ik < K; ik++)
               for(int jk = 0; jk < K; jk++){
                  double d = metseq.dist(ss[knn[a[i]][ik].index], wloc, ss[knn[a[j]][jk].index], wloc);
                  if (d < dmin) dmin = d;
               }

            for(int ik = 0; ik < K; ik++){
               double d = metseq.dist(ss[a[j]], wloc, ss[knn[a[i]][ik].index], wloc);
               if (d < dmin) dmin = d;
               d = metseq.dist(ss[a[i]], wloc, ss[knn[a[j]][ik].index], wloc);
               if (d < dmin) dmin = d;
            }

            double d = metseq.dist(ss[a[i]], wloc, ss[a[j]], wloc);
            if (d < dmin) dmin = d;

            // are they close enough?
            dmin /= lambda;
            if (dmin < knn[a[i]][K - 1].dist || dmin < knn[a[j]][K - 1].dist) adj.add(i, j);
         }
      }

      // now search for connected components
      boolean[] used = new boolean[a.length];
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0 || used[i]) continue;

         HashSet<Integer> members = new HashSet<Integer>();
         findConnComps(i, adj, used, members);
         if (members.size() < 2) continue;

         // System.err.printf("ConnComp: base=%d #members=%d\n", i, members.size());

         // find the max density member
         int iBest = -1;
         double mink = Library.INF;
         Iterator<Integer> it = members.iterator();
         while(it.hasNext()){
            int ix = it.next();
            if (knn[a[ix]][K - 1].dist < mink){
               mink = knn[a[ix]][K - 1].dist;
               iBest = ix;
            }
         }

         // prune all but max density member
         it = members.iterator();
         while(it.hasNext()){
            int ix = it.next();
            if (ix != iBest) a[ix] = -1;
         }
      }
   }

   protected void findConnComps(int ix, AList2<Integer> adj, boolean[] used, Set<Integer> members)
   {
      if (used[ix]) return;
      used[ix] = true;
      members.add(ix);
      int n = adj.size(ix);
      for(int i = 0; i < n; i++)
         findConnComps(adj.get(ix, i), adj, used, members);
   }

   /** if a[i] and a[j] temporally overlap, prune the one with lower density */
   protected void pruneOverlapWindows(int[] a, int K, int nMinOverlap)
   {
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue; // already removed this index
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;
            if (i2loc[a[i]].getNumOverlap(i2loc[a[j]]) >= nMinOverlap){
               if (knn[a[i]][K - 1].dist < knn[a[j]][K - 1].dist) a[j] = -1;
               else{
                  a[i] = -1;
                  break;
               }
            }
         }
      }
   }

   /** if a[i] overlaps a[j]+knn or a[j] overlaps a[i]+knn, prune the one with lower density */
   protected void pruneOverlapKNN(int[] a, int K, int nMinOverlap)
   {
      for(int i = 0; i < a.length; i++){
         if (a[i] < 0) continue; // already removed this index
         for(int j = i + 1; j < a.length; j++){
            if (a[j] < 0) continue;

            // check for overlap
            boolean bOverlap = i2loc[a[i]].overlaps(i2loc[a[j]]);
            if (!bOverlap) for(int k = 0; k < K; k++){
               int n = i2loc[a[i]].getNumOverlap(i2loc[knn[a[j]][k].index]);
               if (n >= nMinOverlap){
                  bOverlap = true;
                  break;
               }
            }
            if (!bOverlap) for(int k = 0; k < K; k++){
               int n = i2loc[a[j]].getNumOverlap(i2loc[knn[a[i]][k].index]);
               if (n >= nMinOverlap){
                  bOverlap = true;
                  break;
               }
            }

            if (bOverlap){
               if (knn[a[i]][K - 1].dist < knn[a[j]][K - 1].dist) a[j] = -1;
               else{
                  a[i] = -1;
                  break;
               }
            }
         }
      }
   }

   protected void calcLocalOptima()
   {
      int N = knn.length;
      int K = knn[0].length;
      SpanList spanLocalOpt = new SpanList(0, N - 1, false);

      // find local optima
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
      System.err.printf("Found %d local optima in %d data points (%.1f%%).\n", nLocalOpt, N, 100.0
            * nLocalOpt / N);

      // prune (temporally) overlapping windows; keep highest density
      int[] a = spanLocalOpt.toIndexArray();
      pruneOverlapWindows(a, K, Math.max(1, wlen / 3));
      // pruneOverlapWindows(a, K, wlen/2);
      // pruneOverlapKNN(a, K, Math.round(wlen*0.7f));

      // try to reduce the number of local optima by detecting redundancy
      // TODO could merge instead of prune
      pruneClusterPairs(a, K, DEN_MERGE_FACTOR);
      // pruneConnectedComponents(a, K, DEN_MERGE_FACTOR);

      // sort the windows by density (descending density = ascending distance to knn)
      ArrayList<Pair<Double, Integer>> ap = new ArrayList<Pair<Double, Integer>>();
      for(int i = 0; i < a.length; i++)
         if (a[i] >= 0) ap.add(new Pair(knn[a[i]][K - 1].dist, a[i], true));
      Collections.sort(ap);

      // build seed list from array
      nLocalOpt = ap.size();
      seeds = new int[nLocalOpt];
      for(int i = 0; i < nLocalOpt; i++)
         seeds[i] = ap.get(i).second;
      System.err.printf("Found %d local optima after overlap and equi-density pruning (%.1f%%).\n",
            nLocalOpt, 100.0 * nLocalOpt / N);

      // update gui for newly found local optima
      lsLocalOpt.setText(String.format("Local Optima: %%d / %d", nLocalOpt));
      lsLocalOpt.setValues(0, 0, nLocalOpt);
      lsLocalOpt.setEnabled(true);
   }

   /*
    * find the next best motif using the current search method
    */
   protected boolean findNext()
   {
      TimerMS timerTotal = new TimerMS();
      int iMotif = motifHMMs.size() + 1;

      // score the baseline (no additional motifs)
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      if (motifHMMs.isEmpty()){
         TimerMS timer = new TimerMS();
         crrBase = HTK.contRec(hmms, gdata.tseries);
         if (crrBase == null){
            System.err.println("Error: HTK / HVite failed");
            return false;
         }
         motifLabels = crrBase.markupSets;
         if (bIncFirstMotif) listInfoGain.add(0);
         listDataLL.add(-crrBase.llTotal / MODSEL_lambda);
         listBIC.add(calcBIC(crrBase.llTotal, hmms));
         listThreshMDL.add(calcMDL(crrBase.llTotal, hmms));
      }
      else hmms.addAll(motifHMMs);

      // score each window (seed) location and explain data
      String sMotifName = String.format("hmm%03d.txt", iMotif);
      int nSeeds = seeds.length;
      if (nSeeds < 1) return false;
      HtkSetupInfo hsi = null;
      ContRecRet[] crrs = new ContRecRet[nSeeds];
      TimerMS timer = new TimerMS();
      for(int iSeed = 0; iSeed < nSeeds; iSeed++){
         System.err.printf(" Testing candidate HMM %d / %d ... ", iSeed + 1, nSeeds);
         timer.reset();

         // run the continuous recognizer and save results
         AbstractHMM hmm = allHMMs[seeds[iSeed]];
         hmm.setName(sMotifName);
         if (hsi == null){
            hmms.add(hmm);
            hsi = HTK.setupContRec(null, hmms, gdata.tseries);
         }
         else{
            // need to replace the last HMM with a new candidate
            hmms.set(hmms.size() - 1, hmm);
            File fHmm = new File(sTmpDir, hmm.getName());
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

      // TODO display graph of (sorted) likelihoods

      // greedily select best motif from remaining set
      final double THRESH_OVERLAP = 0.7;
      final double THRESH_NO_OVERLAP = 0.05;
      while(true){
         // find best motif (null => removed,
         int iBest = -1;
         for(int i = 0; i < nSeeds; i++){
            if (crrs[i] == null) continue;
            if (iBest < 0 || crrs[i].llTotal > crrs[iBest].llTotal) iBest = i;
         }
         if (iBest < 0) break; // no next best found
         System.err.printf("Adding seed: %d (%d) (gain: %f)\n", iBest, seeds[iBest], crrs[iBest].llTotal
               - crrBase.llTotal);
         // out.printf("%d\n", seeds[iBest]); out.flush();

         // initialize new motif from best window location
         allHMMs[seeds[iBest]].setName(String.format("hmm%03d.txt", motifHMMs.size() + 1));
         motifHMMs.add(allHMMs[seeds[iBest]]);
         hmms.clear();
         hmms.add(hmmBG);
         hmms.addAll(motifHMMs);
         if (bIncFirstMotif || motifHMMs.size() > 1)
            listInfoGain.add(crrs[iBest].llTotal - crrBase.llTotal);
         listDataLL.add(-crrs[iBest].llTotal / MODSEL_lambda);
         listBIC.add(calcBIC(crrs[iBest].llTotal, hmms));
         listThreshMDL.add(calcMDL(crrs[iBest].llTotal, hmms));

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
            if (overIJ[i] >= 0.5 && overJI[i] >= 0.5){ // .5, .3?
               // System.err.printf("Overlap: %d & %d\n", iBest, i);
               crrs[i] = null;
               seeds[i] = -1;
            }
         }

         // removel all non-unique (>thresh2) seeds from consideration in this round
         for(int i = 0; i < nSeeds; i++){
            if (crrs[i] == null) continue;
            if (overIJ[i] <= THRESH_NO_OVERLAP && overJI[i] <= THRESH_NO_OVERLAP){
               /*
                * if (listInfoGain.size() > 0){ double baseGain = listInfoGain.last(); double igain =
                * crrs[i].llTotal-crrBase.llTotal; System.err.printf("No overlap: %d & %d igain=%f (d=%f,
                * %%=%f)\n", iBest, i, igain, igain-baseGain, igain/baseGain); }
                */
            }
            else crrs[i] = null;
         }

         assert (seeds[iBest] < 0);
         assert (crrs[iBest] == null);

         break; // TODO if we do this, then we don't really need the loop
      }

      ContRecRet crr = runContRec(-1, true);
      for(TimeMarker tm : crr.markupSets.get(0).getList()){
         if (tm.getTag().equals(motifHMMs.get(motifHMMs.size() - 1).getName()))
            System.err.printf("%d %d\n", tm.getStartIndex(), tm.length());
      }
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
         nWins += Library.getNumSlidingWindowSites(gdata.seqLens[i], wlen, step);
      System.err.printf("nWins=%d (wlen=%d, step=%d)\n", nWins, wlen, step);
      i2loc = new WindowLocation[nWins];
      ss = new Sequence[nWins];
      allHMMs = new AbstractHMM[nWins];

      // extract windows
      int ix = 0;
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = gdata.tseries.get(iSeq);
         int T = gdata.seqLens[iSeq] - wlen;
         for(int t = 0; t <= T; t += step){
            i2loc[ix] = new WindowLocation(iSeq, t, wlen);
            ss[ix] = seq.subseq(i2loc[ix]);
            ix++;
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
      listBIC.clear();
      listDataLL.clear();
      listInfoGain.clear();
      listThreshMDL.clear();
      motifHMMs.clear();
      crrBase = null;
      sgc.updateGraph();
   }

   protected void allKNN(Sequence query, ArrayList<Sequence> data, NBestList<NNInfo> nn, MetricSeq metseq)
   {
      WindowLocation wQuery = i2loc[(Integer)query.getMeta("Index")];
      WindowLocation wloc = new WindowLocation(0, wlen);
      for(Sequence seq : data){
         int j = (Integer)seq.getMeta("Index");
         if (wQuery.overlaps(i2loc[j])) continue;
         double d = metseq.dist(query, wloc, seq, wloc);
         nn.add(new NNInfo(j, d, i2loc[j]));
      }
   }

   protected void fastKNN(Sequence query, SeqTreeNode node, NBestList<NNInfo> nn, MetricSeq metseq)
   {
      if (node.isLeaf()) allKNN(query, node.getData(), nn, metseq);
      else{
         SeqTreeNode left = node.getLeftKid();
         SeqTreeNode right = node.getRightKid();
         Sequence seq1 = (Sequence)left.meta.get("attractor");
         Sequence seq2 = (Sequence)right.meta.get("attractor");
         double d1 = metseq.dist(query, seq1);
         double d2 = metseq.dist(query, seq2);
         SeqTreeNode good, bad;
         if (d1 < d2){
            good = left;
            bad = right;
         }
         else{
            good = right;
            bad = left;
         }
         fastKNN(query, good, nn, metseq);
         if (nn.size() < nn.getN()) fastKNN(query, bad, nn, metseq);
      }
   }

   /** compute the KNN (non-overlapping) for every subsequence */
   protected void calcFastKNN(int wlen, int step, int K)
   {
      this.wlen = wlen;
      this.step = step;
      System.err.printf("Run Fast KNN: wlen=%d  k=%d\n", wlen, K);
      TimerMS timer = new TimerMS();
      TimerMS timer2 = new TimerMS();
      int nSeqs = gdata.getNumSeqs();

      System.err.print("Generating windows... ");
      timer2.reset();
      int nWins = setupWindows();
      System.err.printf("done (%dms).\n", timer2.time());
      knn = new NNInfo[nWins][K];

      // build tree
      System.err.print("Collecting windows... ");
      timer2.reset();
      ArrayList<Sequence> windows = new ArrayList<Sequence>(ss.length);
      for(int i = 0; i < ss.length; i++){
         ss[i].setMeta("Index", i);
         windows.add(ss[i]);
      }
      System.err.printf("done (%dms).\n", timer2.time());

      System.err.print("Building tree... ");
      timer2.reset();
      SeqTree tree = new SeqTree();
      int nSearch = Math.max((K + 1) * 10, (int)Math.sqrt(ss.length));
      tree.build(windows, SeqTree.SplitOrder.MaxMembers, metseq, SeqTree.Stop.MaxMembers, nSearch);
      System.err.printf("done (%s).\n", Library.formatDuration(timer2.time(), 2));

      // calc k smallest distances
      WindowLocation wlocLB = new WindowLocation(0, wlen);
      NBestList<NNInfo> nn = new NBestList<NNInfo>(K);
      nn.setBiggerIsBetter(false);
      nn.setAvoidConflict(true);
      for(int i = 0; i < nWins; i++){
         if (i % 50 == 49) System.err.printf("Calculating KNN for window %d / %d\n", i + 1, nWins);
         nn.clear();
         LBInfo lbi = metseq.calcLBInfo(ss[i]);

         fastKNN(ss[i], tree.getRoot(), nn, metseq);
         // System.err.printf("comp: %d / %d (%.1f%%)\n", nComp, ss.length, 100.0*nComp / ss.length);

         // save the knn
         assert (nn.size() == K) : String.format("nn.size=%d  K=%d\n", nn.size(), K);
         for(int j = 0; j < K; j++)
            knn[i][j] = nn.get(j);
      }
      updateGUIforData();
      System.err.printf("Time for KNN computation: %s\n", Library.formatDuration(timer.time(), 2));
   }

   /** compute the KNN (non-overlapping) for every subsequence */
   protected void calcKNN(int wlen, int step, int K)
   {
      this.wlen = wlen;
      this.step = step;
      System.err.printf("Run KNN: wlen=%d  k=%d\n", wlen, K);
      TimerMS timer = new TimerMS();
      int nSeqs = gdata.getNumSeqs();
      int nWins = setupWindows();
      knn = new NNInfo[nWins][K];

      // calc k smallest distances
      WindowLocation wlocLB = new WindowLocation(0, wlen);
      NBestList<NNInfo> nn = new NBestList<NNInfo>(K);
      nn.setBiggerIsBetter(false);
      nn.setAvoidConflict(true);
      for(int i = 0; i < nWins; i++){
         if (i % 10 == 0) System.err.printf("Calculating KNN for window %d / %d\n", i + 1, nWins);
         nn.clear();
         LBInfo lbi = metseq.calcLBInfo(ss[i]);
         for(int j = 0; j < nWins; j++){
            // we don't want overlapping subsequences
            if (i2loc[i].overlaps(i2loc[j])) continue;

            // only calc full if lower bound is too high or not enough neighbors yet
            if ((nn.size() < K) || (metseq.lowerBound(lbi, ss[j], wlocLB) < nn.largest().dist))
               nn.add(new NNInfo(j, metseq.dist(ss[i], ss[j]), i2loc[j]));

            assert (nn.largest().dist + 1e-9 > nn.smallest().dist);
         }

         // save the knn
         assert (nn.size() == K);
         for(int j = 0; j < K; j++)
            knn[i][j] = nn.get(j);
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
         wlen = Integer.parseInt(st.nextToken());
         step = Integer.parseInt(st.nextToken());
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
         lsWinLen.setValue(wlen);
         lsStep.setValue(step);
         lsKNN.setValue(K);
         updateGUIforData();
         System.err.printf("Loaded KNNs: N=%d  K=%d  wlen=%d  step=%d\n", N, K, wlen, step);
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
         out.printf("%d %d %d %d\n", knn.length, knn[0].length, wlen, step);
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

   /** save hmms in given path, append hmm index and add .hmm extension */
   protected boolean saveHMMs(File fBase)
   {
      String s = Library.getCanonical(fBase.getAbsolutePath());
      String sPath = Library.getPath(s);
      String sBase = Library.getFileName(s);

      System.err.printf("path: %s\n", sPath);
      System.err.printf("title base: %s\n", sBase);

      for(int i = 0; i < motifHMMs.size(); i++){
         AbstractHMM hmm = motifHMMs.get(i);
         assert false : "not yet implemented";
      }

      return true;
   }

   /** @return left-right HMM initialized from the given window index + knn */
   protected HmmLR genHMM(int ix, String name)
   {
      // build initial HMM
      int nStates = Math.max(3, wlen / 2);
      if (nStates % 2 == 0) nStates++; // ensure # states is odd
      nStates = Math.min(wlen * 2 - 1, nStates);
      int nDims = gdata.getNumDims();
      int nSkip = 1;
      HmmLR hmm = new HmmLR(nStates, nSkip, nDims);
      hmm.setMinVar(gdata.getMinVar());
      double edur = (double)wlen / (nStates - 1);
      double pself = Math.max(1.0 - 1.0 / edur, 0.1);
      hmm.initTran(pself, nSkip);
      hmm.setName(name);

      // build the training set
      ArrayList<Sequence> vtrain = new ArrayList<Sequence>();
      vtrain.add(ss[ix]);
      for(int i = 0; i < knn[ix].length; i++)
         vtrain.add(ss[knn[ix][i].index]);

      // learn HMM params
      hmm.init_segk_overlap(vtrain, HMM_OVERLAP_PERCENT);

      // double[][] tranOrig = hmm.saveTran();
      // hmm.train_bw(vtrain); // TODO BW or Viterbi or none!?
      // hmm.blendTran(tranOrig, 0.5);

      // TODO inflate obs variances (blend with global var)?
      /*
       * double alpha = 0.98; for(int i = 0; i < nStates; i++){ GaussianDiagonal dg =
       * (GaussianDiagonal)hmm.getState(i); FeatureVec fv = gdata.getGlobalVar().sqrt()._mul(1.0 - alpha);
       * dg.setVar(dg.getVar().sqrt()._mul(alpha)._add(fv)._sqr()); }
       */

      return hmm;
   }

   protected int getNumParams(AbstractHMM hmm)
   {
      int nStates = hmm.getNumStates();
      int nDims = gdata.getNumDims();
      int nObs = 2 * nStates * nDims; // mean & variance per state per dimension
      int nTran = 2 * nStates - 3; // self, next, skip (but skip = 1-pself-pnext)
      return nObs + nTran;
   }

   protected int getNumParams(ArrayList<AbstractHMM> hmms)
   {
      int n = 0;
      for(AbstractHMM hmm : hmms)
         n += getNumParams(hmm);
      return n;
   }

   protected double calcBIC(double loglik, ArrayList<AbstractHMM> hmms)
   {
      return -loglik / MODSEL_lambda + (hmms.size() + getNumParams(hmms))
            * Math.log(gdata.nTotalSeriesLength) / 2.0;
   }

   protected double calcMDL(double loglik, ArrayList<AbstractHMM> hmms)
   {
      int nDims = gdata.getNumDims();
      double modelBits = -Math.log(hmms.size());
      for(AbstractHMM hmm : hmms){
         int nStates = hmm.getNumStates();
         int nTran = 2 * nStates - 3;
         int nParam = nDims * nStates * 2;
         // int nParam = nStates * 2; // ignore dimensionality?
         modelBits += 18 * (nTran + nParam); // ~18 compressed bits per 32-bit float
         // TODO encode at certain precision and calc entropy
      }

      return -loglik + modelBits;
   }

   protected void learnHMMs()
   {
      System.err.printf("Learning HMMs (%d, wlen=%d)... ", seeds.length, wlen);
      TimerMS timer = new TimerMS();
      for(int i = 0; i < seeds.length; i++){
         int ix = seeds[i];
         allHMMs[ix] = genHMM(ix, String.format("HMM%05d", ix));
      }
      System.err.printf("done (%dms).\n", timer.time());
   }

   protected boolean reEstBG(ArrayList<MarkupSet> labels)
   {
      int nSeqs = gdata.getNumSeqs();
      int nDims = gdata.getNumDims();
      GaussianDynDiag gdd = new GaussianDynDiag(nDims);

      // collect all occurrences of background
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         MarkupSet marks = labels.get(iSeq);
         Sequence seq = gdata.tseries.get(iSeq);
         for(TimeMarker tm : marks.getList()){
            if (!BG.equals(tm.getTag())) continue;
            int a = tm.getStartIndex();
            int b = tm.getStopIndex();
            for(int i = a; i < b; i++)
               gdd.add(seq.get(i), false);
         }
      }

      GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      gdd.update();
      dg.setMean(gdd.getMean());
      dg.setVar(gdd.getVar());

      // now that we've reestimated, rerun the cont rec system
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      hmms.addAll(motifHMMs);
      crrBase = HTK.contRec(hmms, gdata.tseries);
      if (crrBase == null){
         System.err.println("Error: HTK / HVite failed after re-estimation.");
         return false;
      }
      motifLabels = crrBase.markupSets;
      if (lsMotif.getValue() > 0) highlightMotif(lsMotif.getValue() - 1, motifLabels);
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
         System.err.printf("(%dms).\n", timer.time());
      }

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
      ContRecRet crr = HTK.contRec(hmms, gdata.tseries);
      if (crr == null){
         System.err.println("Error: HTK / HVite failed after re-estimation.");
         return null;
      }

      if (bSaveResults){
         crrBase = crr;
         motifLabels = crrBase.markupSets;
         if (lsMotif.getValue() > 0) highlightMotif(lsMotif.getValue() - 1, motifLabels);
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
            double r = metseq.dist(seqs[iBase], seqs[i]);
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
            double r = metseq.dist(seqs[iBase], seqs[i]);
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

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();
      if (src == btKNN){
         int wlen = lsWinLen.getValue();
         int step = lsStep.getValue();
         int k = lsKNN.getValue();
         calcKNN(wlen, step, k);
      }
      else if (src == btFastKNN){
         int wlen = lsWinLen.getValue();
         int step = lsStep.getValue();
         int k = lsKNN.getValue();
         calcFastKNN(wlen, step, k);
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
      else if (src == btLocalOpt){
         btLocalOpt.setEnabled(false);
         calcLocalOptima();
         learnHMMs();
         btFindNext.setEnabled(true);
         lsNumFind.setEnabled(true);
      }
      else if (src == btFindNext){
         TimerMS timer = new TimerMS();
         int nFind = lsNumFind.getValue();
         if (nFind == lsNumFind.getMaximum()) nFind = Integer.MAX_VALUE;
         for(int iter = 0; iter < nFind; iter++){
            if (!findNext()){
               System.err.printf("Error: failed to find next motif\n");
               break;
            }

            // rebuild seed list
            System.err.printf("# seeds before: %d\n", seeds.length);
            MyIntList a = new MyIntList();
            for(int i = 0; i < seeds.length; i++)
               if (seeds[i] >= 0) a.add(seeds[i]);
            seeds = a.toArray();
            System.err.printf("# seeds after: %d\n", seeds.length);

            lsMotif.setText(String.format("Motif: %%d / %d", motifHMMs.size()));
            lsMotif.setValues(motifHMMs.size(), 0, motifHMMs.size());
            lsMotif.setEnabled(true);
            btScoreMotifs.setEnabled(true);
            btSaveHMMs.setEnabled(true);
            btrOptAcc.setEnabled(true);
            btrOptRecall.setEnabled(true);
            btReEst.setEnabled(true);
            btReEstBG.setEnabled(true);
            if (seeds.length < 1) break;
         }
         System.err.printf("Total time for find-next: %s\n", Library.formatDuration(timer.time()));
      }
      else if (src == btScoreMotifs){
         int nMotifs = motifHMMs.size();
         int nScore = sgc.getDisplayIndex();
         if (!bIncBg) nScore++;
         if (!bIncFirstMotif) nScore += (bIncBg ? 2 : 1);
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
      else if (src == btNbrGraph){
         int ix = lsLocalOpt.getValue() - 1;
         ix = seeds[ix];
         displayNeighborhood(ix);
      }
      else if (src == btSaveHMMs){
         JFileChooser fc = Library.buildFileChooser("hmm", "HMMs (*.hmm)");
         if (fc.showSaveDialog(frame) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            if (!saveHMMs(f)){
               System.err.println("Error: failed to save HMMs (base file:");
               System.err.println(" " + Library.getCanonical(f) + ")");
               JOptionPane.showMessageDialog(frame, "Failed to save HMMs", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
            else JOptionPane.showMessageDialog(frame, "Saved HMMs", "Save Successful",
                  JOptionPane.INFORMATION_MESSAGE);
         }
      }
   }

   protected void displayNeighborhood(int ix)
   {
      int nWins = ss.length;

      // calc dist to all other subsequences
      ArrayList<Pair<Double, Integer>> dlist = new ArrayList<Pair<Double, Integer>>();
      WindowLocation sswloc = new WindowLocation(0, wlen);
      boolean[] bAvail = new boolean[nWins];
      Arrays.fill(bAvail, true);
      bAvail[ix] = false;
      for(int i = 0; i < nWins; i++){
         if (i2loc[ix].overlaps(i2loc[i])){
            bAvail[i] = false;
            continue;
         }
         dlist.add(new Pair(metseq.dist(ss[ix], sswloc, ss[i], sswloc), i, true));
      }
      System.err.printf("dlist size: %d\n", dlist.size());

      // sort by distance and calc non-overlapping set
      Collections.sort(dlist);
      int nNbrs = dlist.size();
      ArrayList<Pair<Double, Integer>> sslist = new ArrayList<Pair<Double, Integer>>();
      for(int iNbr = 0; iNbr < nNbrs; iNbr++){
         Pair<Double, Integer> nbr = dlist.get(iNbr);
         if (!bAvail[nbr.second]) continue;
         sslist.add(nbr);
         bAvail[nbr.second] = false;
         for(int i = nbr.second - 1; i >= 0; i--){
            if (!i2loc[nbr.second].overlaps(i2loc[i])) break;
            bAvail[i] = false;
         }
         for(int i = nbr.second + 1; i < nWins; i++){
            if (!i2loc[nbr.second].overlaps(i2loc[i])) break;
            bAvail[i] = false;
         }
      }
      System.err.printf("sslist size: %d\n", sslist.size());

      // make sure we don't have any overlapping neighbors
      nNbrs = sslist.size();
      for(int i = 0; i < nNbrs; i++)
         for(int j = i + 1; j < nNbrs; j++){
            assert (!i2loc[sslist.get(i).second].overlaps(i2loc[sslist.get(j).second])) : String.format(
                  "Error: %d and %d overlap!  %s <-> %s", i, j, i2loc[sslist.get(i).second], i2loc[sslist
                        .get(j).second]);
         }
      System.err.printf("No overlaps!\n");

      for(int i = 0; i < bAvail.length; i++){
         assert (!bAvail[i]) : String.format("left-over spots! %d  %s", i, i2loc[dlist.get(i).second]);
      }
      for(int i = 1; i < sslist.size(); i++){
         assert (sslist.get(i - 1).first - 1e-9 < sslist.get(i).first);
      }

      Sequence data = new Sequence("Distance & Density");
      nNbrs = sslist.size();
      for(int i = 0; i < nNbrs; i++){
         double den = (double)(i + 2) / sslist.get(i).first;
         // data.add(new FeatureVec(2, sslist.get(i).first, den));
         data.add(new FeatureVec(1, den));
      }

      JFrame win = new JFrame(String.format("Neighborhood for Local Optima #%d %s", ix, i2loc[ix]));
      win.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      win.setSize(600, 400);
      Library.centerWin(win, frame);
      JPanel mainp = new JPanel(new BorderLayout());
      win.setContentPane(mainp);

      LineGraph lg = new LineGraph(data);
      lg.setShowAll(true);
      mainp.add(new GraphComplexLite(lg), BorderLayout.CENTER);

      win.setVisible(true);
   }

   protected void highlightWindow(int iWin)
   {
      int nSeqs = gdata.getNumSeqs();

      // clear the highlights
      for(int i = 0; i < nSeqs; i++)
         crbGT[i].clearHighlight();

      // highlight a real local optima
      if (iWin >= 0){
         WindowLocation wloc = i2loc[iWin];
         Range r = wloc.getRange();
         crbGT[wloc.iSeries].addHighlight(r);
         for(int i = 0; i < knn[iWin].length; i++){
            wloc = i2loc[knn[iWin][i].index];
            r = wloc.getRange();
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
      if (iMotif >= 0){
         String sMotifName = motifHMMs.get(iMotif).getName();
         assert (nSeqs == labels.size()) : String.format("mismatch: nseqs=%d  labels.size=%d", nSeqs, labels
               .size());
         for(int iSeq = 0; iSeq < nSeqs; iSeq++){
            for(TimeMarker tm : labels.get(iSeq).getList()){
               if (tm.getTag() == null || tm.getTag().equals(BG)) continue;
               if (!tm.getTag().equals(sMotifName)) continue;
               Range r = new Range(tm.getStartIndex(), tm.getStopIndex() - 1);
               crbGT[iSeq].addHighlight(r);
            }
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
      System.err.printf("Overlap map (%d x %d) [RAW]:\n", omap.length, omap[0].length);
      for(int i = 0; i < omap.length; i++){
         for(int j = 0; j < omap[i].length; j++){
            System.err.printf("%3d  ", omap[i][j]);
         }
         System.err.println();
      }
      cri.cleanOMap(omap, CLEAN_FRAC);
      int[] nFoundLabels = new int[omap.length];
      for(WordSpot spot : cri.getSpots())
         nFoundLabels[spot.iClass]++;

      System.err.printf("Overlap map (%d x %d) [CLEAN]:\n", omap.length, omap[0].length);
      for(int i = 0; i < omap.length; i++){
         for(int j = 0; j < omap[i].length; j++){
            System.err.printf("%3d  ", omap[i][j]);
         }
         System.err.println();
      }

      // calc best class permutation and remap spots
      System.err.print("Calculating best cluster <-> class mapping... ");
      System.err.print("\nFound: ");
      for(int i = 0; i < nFoundLabels.length; i++)
         System.err.printf("%3d  ", nFoundLabels[i]);
      System.err.println();
      System.err.print(" True: ");
      for(int i = 0; i < gdata.nExamples.length; i++)
         System.err.printf("%3d  ", gdata.nExamples[i]);
      System.err.println();

      // short[] map = SupTest.calcBestMappingFromConfMatrix(confm, bOptAcc);
      boolean bOptAcc = btrOptAcc.isSelected();
      short[] map = ContRecInfo.calcBestMappingFromOMap(omap, nFoundLabels, gdata.nExamples, bOptAcc);
      cri.remapSpots(map, gdata.labData.size());
      System.err.println("done.");

      for(int i = 0; i < map.length; i++){
         if (map[i] < 0) continue;
         System.err.printf("%s -> %s\n", motifHMMs.get(i).getName(), gdata.classes[map[i]]);
      }

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

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();

      if (src == lsLocalOpt){
         int ilo = lsLocalOpt.getValue();
         if (ilo == 0){
            highlightWindow(-1);
            btNbrGraph.setEnabled(false);
         }
         else{
            int ix = seeds[ilo - 1];
            // System.err.printf("Local Optima #%d / %d = Window #%d / %d (%.2f, %s)\n", ilo, spanLocalOpt
            // .size(), ix, knn.length, knn[ix][knn[ix].length - 1].dist, i2loc[ix]);
            highlightWindow(ix);
            btNbrGraph.setEnabled(true);
         }
      }
      else if (src == lsMotif){
         int iMotif = lsMotif.getValue() - 1;
         highlightMotif(iMotif, motifLabels);
         if (iMotif>=0) mviz.setHMM(motifHMMs.get(iMotif));
         else mviz.setHMM(null);
      }
      else if (src == lsNumFind){
         int nf = lsNumFind.getValue();
         if (nf == lsNumFind.getMaximum()) btFindNext.setText("Find All Motifs");
         else if (nf == 1) btFindNext.setText("Find Next Motif");
         else btFindNext.setText(String.format("Find Next %d Motifs", nf));
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
         if (!bIncBg) nDisp++;
         if (!bIncFirstMotif) nDisp += (bIncBg ? 2 : 1);
         btScoreMotifs.setEnabled(nDisp >= 1);
         btrOptAcc.setEnabled(nDisp >= 1);
         btrOptRecall.setEnabled(nDisp >= 1);
         btScoreMotifs.setText(String.format("Score Motifs (%d)", nDisp));
         System.err.printf("Auto-stop (No Pos): %d motifs (nDisp=%d)\n", calcStop(nDisp), nDisp);
      }
   }

   protected void run(String sFile, int wlen, int step) throws Exception
   {
      this.wlen = wlen;
      this.step = step;
      setupWindows();

      MyIntList seeds = new MyIntList();
      BufferedReader in = new BufferedReader(new FileReader(sFile));
      while(true){
         String line = in.readLine();
         if (line == null) break;
         seeds.add(Integer.parseInt(line));
      }

      System.err.printf("Found %d seeds.\n", seeds.size());

      // for(int nMotifs=1; nMotifs<=seeds.size(); nMotifs++){
      int nMotifs = 15;
      // System.err.printf("Running with %d / %d motifs...\n", nMotifs, seeds.size());
      createBgModel();
      motifHMMs = new ArrayList<AbstractHMM>();
      for(int i = 0; i < nMotifs; i++){
         AbstractHMM hmm = genHMM(seeds.get(i), String.format("HMM%05d", i + 1));
         motifHMMs.add(hmm);
      }
      for(int i = 0; i < 3; i++){
         ContRecRet crr = runContRec(-1, true);
         reEstMotifs(crr.markupSets);
      }
      ContRecRet crr = runContRec(-1, true);
      scoreMotifs(crr.markupSets);
      // System.err.printf("Auto-stop: %d motifs (total #motifs=%d)\n", calcStop(nMotifs), nMotifs);

      MarkupSet marks = crr.markupSets.get(6);
      for(TimeMarker tm : marks.getList()){
         System.err.println(tm);
      }

      // }

      System.exit(1);
   }

   public static void main(String args[]) throws Exception
   {
      Runtime rt = Runtime.getRuntime();
      System.err.printf("Java Version: %s   VM: %s (%s)\n", System.getProperty("java.version"), System
            .getProperty("java.vm.version"), System.getProperty("java.vm.name"));
      System.err.printf("Memory (free/total/max): %dk / %dk / %dk\n", rt.freeMemory() / 1024, rt
            .totalMemory() / 1024, rt.maxMemory() / 1024);

      final String sDefFile = "/home/dminn/research/kdm/data//synth/chars/easy/easy.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/ah-data/exercise/sub8/accgyr.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/dim13/ae.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/ae.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/asl-joseph/asl16D.def";
      // final String sDefFile = "/home/dminn/research/kdm/data/keogh/shuttle.def";

      SwingUtilities.invokeLater(new Runnable() {
         public void run()
         {
            final ExploreKNN ek = new ExploreKNN();
            if (!ek.setup(sDefFile)){
               System.err.println("Explore-KNN setup failed.");
               System.exit(1);
            }
      
        // ek.out = new PrintWriter(new FileWriter("/tmp/seeds.txt"));
      
            final BlockDlg dlg = new BlockDlg(ek.frame, true); 
            dlg.display("Loading default data", new Thread() {
               public void run() {
                 if(!ek.load(new File("/home/dminn/knn-easy.txt"))){
                     //if(!ek.load(new File("/home/dminn/exercise-knn-15-2-5-fast.txt"))){
                     //if (!ek.load(new File("/home/dminn/tidigits-knn-29-2-5-fast.txt"))){ 
                     System.err.println("Error: failed to load default data!"); 
                     System.exit(1);
                     }
                  dlg.close();
               }});
         }
      });
      
      // ek.run("/home/dminn/research/pubs/2007/aaai/experiments/seeds-speech-29-fast.txt", 29, 2);
      // ek.run("/home/dminn/research/pubs/2007/aaai/experiments/seeds-exercise-15-fast.txt", 15, 2);
   }

}
