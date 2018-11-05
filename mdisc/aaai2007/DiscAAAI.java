package mdisc.aaai2007;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.io.*;
import kdm.util.*;
import kdm.gui.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.mlpr.*;
import kdm.mlpr.dataTree.*;
import kdm.mlpr.htk.*;
import kdm.metrics.*;
import kdm.tools.*;

import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import gnu.getopt.*;
import mdisc.VizTool.*;

/**
 * Tool for exploring time series data with different views and operations.
 */
public class DiscAAAI implements ActionListener, ChangeListener
{
   public static final double lambda = 1.0;
   public static final String BG = "bg";
   
   public static final int METHOD_RANDOM = 0;
   public static final int METHOD_UNLIKELY = 1;
   public static final int METHOD_CLUSTER = 2;
   public static final int METHOD_IJCAI = 3;
   public static final int METHOD_KNN_FULL = 4;
   public static final int METHOD_EXHAUSTIVE = 5;
   public static final int METHOD_COMBINE = 6;
   public static final int METHOD_NONE = 7;

   protected GlobalData gdata;
   protected ColoredRangeBar crbGT[];
   protected Color[] colors;
   protected VertChooseContainer vcc;
   protected JTextField tfMotif, tfRadius, tfEstRProp;
   protected LabeledSlider lsWinLen, lsRandomN, lsUnlikelyN, lsExhaustStep, lsTree, lsMotif, lsFindNMotifs;
   protected LabeledSlider lsSegExpLength, lsSegMixSize, lsSegCompLen, lsSegMotif;
   protected LabeledSlider lsPAA, lsSAX, lsKnnFull;
   protected JButton btNextMotif, btReset, btRemoveLast, btScore, btSegRun, btRealRadius, btAddIJCAI;
   protected JRadioButton btrFixedLeaf, btrNumLeaves, btrSplitSize, btrSplitWidth;
   protected JRadioButton btrOneLen, btrMultiLen;
   protected JRadioButton btrOptAcc, btrOptRecall;
   protected JRadioButton btrIJCAIseeds, btrIJCAIhmms;
   protected JCheckBox btcNodeRandom, btcNodeProto, btcEstR;
   protected JCheckBox btcRandom, btcUnlikely, btcCluster, btcExhaustive, btcIJCAI, btcKnnFull;
   protected ComboPanel cpSearch;
   protected JComboBox cbWinLen, cbScale;
   protected ScoreGraphComplex sgc;
   protected ButtonGroup bgTreeStop, bgTreeSplit, bgLen, bgOpt, btIJCAI;
   protected MetricSeq metseq = new DTWSimple(0.1);
   protected int nMaxLeafSize = 10;
   protected int nTotalLeaves = 10;
   protected int iHighlightMotif = -1;
   protected MyDoubleList listDataLL = new MyDoubleList();
   protected MyDoubleList listAIC = new MyDoubleList();
   protected MyDoubleList listBIC = new MyDoubleList();
   protected MyDoubleList listThreshMDL = new MyDoubleList();
   protected AbstractHMM hmmBG;
   protected ContRecRet curCRR;

   /** these HMMs are added to see list before greedy search for next best model */
   protected ArrayList<AbstractHMM> globalHmmSeeds = new ArrayList<AbstractHMM>();

   /** HMMs learned from IJCAI seeds */
   protected ArrayList<AbstractHMM> hmmsIJCAI;

   protected ArrayList<AbstractHMM> motifHMMs = new ArrayList<AbstractHMM>();
   protected ArrayList<MarkupSet> motifLabels = new ArrayList<MarkupSet>();
   protected ArrayList<MarkupSet> segLabels = new ArrayList<MarkupSet>();
   protected SeqTree treeCluster;

   public static void usage()
   {
      System.err.println();
      System.err.println("Usage: java ~.DiscApp [options] <data def file>");
      System.err.println(" Options: ");
      System.err.println("  -help/?                 Usage information");
      System.err.println();
      System.exit(1);
   }

   protected void reset()
   {
      globalHmmSeeds.clear();
      motifHMMs.clear();
      motifLabels.clear();
      curCRR = null;
      lsMotif.setText("No Motifs");
      lsMotif.setValues(0, 0, 0);
      lsMotif.setEnabled(false);
      btRemoveLast.setEnabled(false);
      btReset.setEnabled(false);
      btrOptAcc.setEnabled(false);
      btrOptRecall.setEnabled(false);
      btScore.setEnabled(false);
      tfMotif.setText("");
      listDataLL.clear();
      listAIC.clear();
      listBIC.clear();
      listThreshMDL.clear();
      sgc.updateGraph();
      gdata.chiuView.occList = null;
      gdata.chiuView.seeds = null;
      hmmsIJCAI = null;
      enable("ijcai", true);
      enable("cluster", true);
   }

   protected boolean enable(String sMethod, boolean bEnable)
   {
      if (sMethod.equals("ijcai")){
         lsPAA.setEnabled(bEnable);
         lsSAX.setEnabled(bEnable);
         btcEstR.setEnabled(bEnable);
         btrIJCAIseeds.setEnabled(bEnable);
         btrIJCAIhmms.setEnabled(bEnable);
         btAddIJCAI.setEnabled(bEnable);
         if (bEnable){
            boolean bEstR = btcEstR.isSelected();
            tfEstRProp.setEnabled(bEstR);
            tfRadius.setEnabled(!bEstR);
            btRealRadius.setEnabled(!bEstR);
         }
         else{
            tfEstRProp.setEnabled(false);
            tfRadius.setEnabled(false);
            btRealRadius.setEnabled(false);
         }
      }
      else if (sMethod.equals("cluster")){
         btrNumLeaves.setEnabled(bEnable);
         btrFixedLeaf.setEnabled(bEnable);
         btrSplitSize.setEnabled(bEnable);
         btrSplitWidth.setEnabled(bEnable);
         btcNodeProto.setEnabled(bEnable);
         btcNodeRandom.setEnabled(bEnable);
         lsTree.setEnabled(bEnable);
      }
      else{
         assert false : String.format("unrecognized enable/disable name (%s)", sMethod);
         return false;
      }
      return true;
   }

   protected void createBgModel()
   {
      int nDims = gdata.getNumDims();
      hmmBG = new HmmLR(1, nDims);
      hmmBG.setName(BG);
      hmmBG.setPiLeave(0, Math.log(0.2));
      GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      dg.setMean(gdata.getGlobalMean());
      dg.setVar(gdata.getGlobalVar());
   }

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
      LabeledDataLoader.dumpDataSummary(gdata.labData);
      gdata.chiuView = new ChiuView(gdata);

      // build the background model
      createBgModel();

      // create the window
      System.err.print("Creating visualization... ");
      timer.reset();
      JFrame frame = new JFrame(String.format("Discover (AAAI) - %s", Library.getTitle(sDefFile)));
      gdata.frame = frame;
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setSize(1200, 700);
      Library.centerWin(frame, null);
      frame.setLocation(frame.getLocation().x, 0);
      JPanel mainp = new JPanel(new BorderLayout());
      frame.setContentPane(mainp);

      mainp.add(buildHtkPanel(), BorderLayout.CENTER);

      System.err.printf("done (%dms).\n", timer.time());
      frame.setVisible(true);
      return true;
   }

   protected JPanel buildHtkPanel()
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
      sgc.addCurve("Log-likelihood", listDataLL);
      sgc.addCurve("AIC", listAIC);
      sgc.addCurve("BIC", listBIC);
      sgc.addCurve("MDL", listThreshMDL);
      leftp.add(sgc, BorderLayout.SOUTH);

      // create the right control panel
      vcc = new VertChooseContainer();
      JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftp, vcc);
      split.setDividerLocation(Math.max(gdata.frame.getWidth() / 2, gdata.frame.getWidth() - 280));
      split.setResizeWeight(1.0);
      split.setDividerSize(7);
      split.setBorder(null);
      mainp.add(split, BorderLayout.CENTER);

      // create the control pane for motif discovery
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      vcc.addPane("Discover Motifs", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      lsWinLen = new LabeledSlider("Window Length: %d", 4, 50);
      lsWinLen.addChangeListener(this);
      rightp.add(lsWinLen);

      JPanel p = new JPanel();
      cbWinLen = new JComboBox(new String[] { "Exercise (15)", "TIDIGITS (22)", "Custom"});
      cbWinLen.addActionListener(this);
      p.add(cbWinLen);
      cbScale = new JComboBox(new String[] { " Window Length", " +/-25%", " +/-50%", " +/-25% & +/-50% " });
      p.add(cbScale);
      rightp.add(p);

      btrOneLen = new JRadioButton("Select");
      btrOneLen.setSelected(true);
      btrMultiLen = new JRadioButton("Multi");
      bgLen = new ButtonGroup();
      bgLen.add(btrOneLen);
      bgLen.add(btrMultiLen);
      p = new JPanel();
      p.add(new JLabel("Length:"));
      p.add(btrOneLen);
      p.add(btrMultiLen);
      rightp.add(p);

      rightp.add(new Divider());

      cpSearch = new ComboPanel("Search: ");
      rightp.add(cpSearch);

      // random N search pane
      p = new JPanel(new VerticalLayout(-1, 2));
      lsRandomN = new LabeledSlider("Random N: %d", 1, 100);
      p.add(lsRandomN);
      cpSearch.addPanel("Random N", p);

      // least likely N search pane
      p = new JPanel(new VerticalLayout(-1, 2));
      lsUnlikelyN = new LabeledSlider("Unlikely N: %d", 1, 100);
      p.add(lsUnlikelyN);
      cpSearch.addPanel("Least Likely N", p);

      // cluster search pane
      p = new JPanel(new VerticalLayout(-1, 0));
      btrFixedLeaf = new JRadioButton("Leaf Size", false);
      btrFixedLeaf.setToolTipText("Don't expand node if fewer than given number of members");
      btrFixedLeaf.addActionListener(this);
      btrNumLeaves = new JRadioButton("N Leaves", true);
      btrNumLeaves.setToolTipText("Stop expanding tree at fixed number of leaves");
      btrNumLeaves.addActionListener(this);
      JPanel q = new JPanel();
      q.add(new JLabel("Stop:"));
      q.add(btrNumLeaves);
      q.add(btrFixedLeaf);
      p.add(q);
      bgTreeStop = new ButtonGroup();
      bgTreeStop.add(btrFixedLeaf);
      bgTreeStop.add(btrNumLeaves);

      btrSplitSize = new JRadioButton("Size", false);
      btrSplitSize.setToolTipText("Split node with most members");
      btrSplitWidth = new JRadioButton("Width", true);
      btrSplitWidth.setToolTipText("Split node with largest width");
      q = new JPanel();
      q.add(new JLabel("Split:"));
      q.add(btrSplitSize);
      q.add(btrSplitWidth);
      p.add(q);
      bgTreeSplit = new ButtonGroup();
      bgTreeSplit.add(btrSplitSize);
      bgTreeSplit.add(btrSplitWidth);

      btcNodeRandom = new JCheckBox("Random", true);
      btcNodeRandom.addActionListener(this);
      btcNodeProto = new JCheckBox("Protoype", true);
      btcNodeProto.addActionListener(this);
      q = new JPanel();
      q.add(btcNodeProto);
      q.add(btcNodeRandom);
      p.add(q);

      lsTree = new LabeledSlider("Num Leaves: %d", 2, 100);
      p.add(lsTree);
      cpSearch.addPanel("Cluster", p);

      // IJCAI seach pane
      p = new JPanel(new VerticalLayout(-1, 4));
      lsPAA = new LabeledSlider("# PAA Segments: %d", 3, 10);
      p.add(lsPAA);
      lsSAX = new LabeledSlider("# SAX Symbols: %d", 3, 5);
      p.add(lsSAX);
      q = new JPanel();
      btcEstR = new JCheckBox("Estimate R", true);
      btcEstR.addActionListener(this);
      tfEstRProp = new JTextField(5);
      q.add(btcEstR);
      q.add(Box.createHorizontalStrut(4));
      q.add(new JLabel("Prop:"));
      q.add(tfEstRProp);
      p.add(q);
      q = new JPanel();
      tfRadius = new JTextField(5);
      tfRadius.setEnabled(false);
      q.add(new JLabel(" R:"));
      q.add(tfRadius);
      btRealRadius = new JButton("Real");
      btRealRadius.addActionListener(this);
      btRealRadius.setEnabled(false);
      q.add(btRealRadius);
      p.add(q);
      q = new JPanel();
      q.add(new JLabel("Return:"));
      btrIJCAIseeds = new JRadioButton("Seed Locs");
      q.add(btrIJCAIseeds);
      btrIJCAIhmms = new JRadioButton("HMMs", true);
      q.add(btrIJCAIhmms);
      p.add(q);
      btIJCAI = new ButtonGroup();
      btIJCAI.add(btrIJCAIseeds);
      btIJCAI.add(btrIJCAIhmms);
      btAddIJCAI = new JButton("Add all HMMs");
      btAddIJCAI.addActionListener(this);
      p.add(btAddIJCAI);
      cpSearch.addPanel("IJCAI", p);

      // knn-full search pane
      p = new JPanel(new VerticalLayout(-1, 4));
      lsKnnFull = new LabeledSlider("K: %d", 1, 10);
      p.add(lsKnnFull);
      cpSearch.addPanel("KNN (Full)", p);
      
      // exhaustive search pane
      p = new JPanel(new VerticalLayout(-1, 4));
      lsExhaustStep = new LabeledSlider("Step: %d", 1, 50);
      p.add(lsExhaustStep);
      cpSearch.addPanel("Exhaustive", p);

      // combination search pane
      p = new JPanel(new VerticalLayout(-1, 1));
      btcRandom = new JCheckBox("Random N");
      btcRandom.setSelected(true);
      p.add(btcRandom);
      btcUnlikely = new JCheckBox("Unlikely N");
      btcUnlikely.setSelected(true);
      p.add(btcUnlikely);
      btcCluster = new JCheckBox("Cluster");
      btcCluster.setSelected(true);
      p.add(btcCluster);
      btcIJCAI = new JCheckBox("IJCAI");
      btcIJCAI.setSelected(false);
      p.add(btcIJCAI);
      btcKnnFull = new JCheckBox("KNN (Full)");
      btcKnnFull.setSelected(false);
      p.add(btcKnnFull);
      btcExhaustive = new JCheckBox("Exhaustive");
      btcExhaustive.setSelected(false);
      p.add(btcExhaustive);
      cpSearch.addPanel("Combine", p);
      
      // no search pane
      p = new JPanel();
      p.add(new JLabel("--- No Search ---"));
      cpSearch.addPanel("None", p);

      rightp.add(new Divider());

      lsFindNMotifs = new LabeledSlider("Find next %d motif(s)", 1, 50);
      rightp.add(lsFindNMotifs);

      btNextMotif = new JButton("Find Next Motifs");
      btNextMotif.addActionListener(this);
      rightp.add(btNextMotif);

      lsMotif = new LabeledSlider("No Motifs", 0, 0);
      lsMotif.addChangeListener(this);
      lsMotif.setEnabled(false);
      rightp.add(lsMotif);

      tfMotif = new JTextField();
      tfMotif.setEditable(false);
      rightp.add(tfMotif);

      q = new JPanel();
      btRemoveLast = new JButton("Remove Last Motif");
      btRemoveLast.addActionListener(this);
      btRemoveLast.setEnabled(false);
      q.add(btRemoveLast);

      btReset = new JButton("Reset");
      btReset.addActionListener(this);
      btReset.setEnabled(false);
      q.add(btReset);
      rightp.add(q);

      btrOptAcc = new JRadioButton("Accuracy", true);
      btrOptAcc.setEnabled(false);
      btrOptRecall = new JRadioButton("Recall", false);
      btrOptRecall.setEnabled(false);
      bgOpt = new ButtonGroup();
      bgOpt.add(btrOptAcc);
      bgOpt.add(btrOptRecall);

      q = new JPanel();
      q.add(new JLabel("Optimize: "));
      q.add(btrOptAcc);
      q.add(btrOptRecall);
      rightp.add(q);

      btScore = new JButton("Score Motifs");
      btScore.addActionListener(this);
      btScore.setEnabled(false);
      rightp.add(btScore);

      // setup initial values
      cpSearch.setSelectedIndex(4);
      lsWinLen.setValue(15); // TODO 22
      lsExhaustStep.setValue(7);
      lsRandomN.setValue(10);
      lsUnlikelyN.setValue(10);
      lsMotif.setValue(0);
      lsTree.setValue(nMaxLeafSize);
      lsPAA.setValue(4);
      lsSAX.setValue(3);
      lsKnnFull.setValue(5);
      tfRadius.setText("10.0");
      tfEstRProp.setText("5.0");      

      // /////////////////////////////////////////////////////////
      // create the control pane for EM segmentation
      rightp = new VerticalScrollPanel(new VerticalLayout(-1, 8));
      rightp.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
      vcc.addPane("Segment EM", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      lsSegExpLength = new LabeledSlider("Expected Length: %d", 8, 100);
      rightp.add(lsSegExpLength);
      lsSegMixSize = new LabeledSlider("Mixture Size: %d", 2, 20);
      rightp.add(lsSegMixSize);
      lsSegCompLen = new LabeledSlider("Component Length: %d", 3, 20);
      rightp.add(lsSegCompLen);
      btSegRun = new JButton("Train Mixture Model");
      btSegRun.addActionListener(this);
      rightp.add(btSegRun);
      lsSegMotif = new LabeledSlider("No Motifs", 0, 0);
      lsSegMotif.addChangeListener(this);
      lsSegMotif.setEnabled(false);
      rightp.add(lsSegMotif);

      lsSegExpLength.setValue(15);
      lsSegMixSize.setValue(8);
      lsSegCompLen.setValue(10);

      return mainp;
   }

   /** @return generic step size for a window of the given length */
   public int getStepForLen(int wlen)
   {
      return Math.max(2, wlen / 10);
   }

   /**
    * Choose list of candidate locations (seed motifs) to search
    * 
    * @param iMotif index of the motif we're searching for
    * @param kSearch search method
    * @param labels list of markup sets for the motifs so far (used for N Least Likely method)
    * @return list of candidate locations to search
    */
   protected ArrayList<WindowLocation> getCandidateLocs(int iMotif, int kSearch, ArrayList<MarkupSet> labels)
   {
      final ArrayList<WindowLocation> clocs = new ArrayList<WindowLocation>();
      int wlen = lsWinLen.getValue();
      int nSeqs = gdata.getNumSeqs();

      if (kSearch == METHOD_RANDOM){ // N random locations ////////////////////
         int N = lsRandomN.getValue();
         int nSpots = gdata.nTotalSeriesLength - gdata.getNumSeqs() * (wlen - 1);
         if (N * 4 >= nSpots){
            System.err.printf("Warning: searching for %d random locations from only %d spots.", N, nSpots);
            return null;
         }

         HashSet<Integer> uniq = new HashSet<Integer>();
         int[] ret = null;
         while(clocs.size() < N){
            int ix = Library.random(gdata.nTotalSeriesLength);
            if (uniq.contains(ix)) continue;
            uniq.add(ix);
            uniq.add(ix - 1);
            uniq.add(ix + 1);
            ret = Library.getArrayOffset(ix, gdata.seqLens, ret);
            if (ret[1] + wlen <= gdata.seqLens[ret[0]]) clocs.add(new WindowLocation(ret[0], ret[1], wlen));
         }
      }
      else if (kSearch == METHOD_UNLIKELY){ // unlikely N ////////////////////
         // find the N matches that fit the worst
         int N = lsUnlikelyN.getValue();
         Comparator comp = new Comparator() {
            public int compare(Object o1, Object o2)
            {
               ScoredWindow swin1 = (ScoredWindow)o1;
               ScoredWindow swin2 = (ScoredWindow)o2;

               double meanLL1 = swin1.score;
               double meanLL2 = swin2.score;

               if (meanLL1 < meanLL2) return 1;
               if (meanLL1 > meanLL2) return -1;
               return 0;
            }
         };
         NBestList<ScoredWindow> nbl = new NBestList<ScoredWindow>(N, comp);
         for(int iSeries = 0; iSeries < nSeqs; iSeries++){
            Sequence seq = gdata.tseries.get(iSeries);
            MarkupSet marks = labels.get(iSeries);
            for(TimeMarker tm : marks.getList()){
               double score = (Double)tm.getMeta();
               if (tm.length() > 1){
                  // get the viterbi parse and check each frame
                  AbstractHMM hmm = getMotifFromName(tm.getTag());
                  Sequence subseq = seq.subseq(tm.getStartIndex(), tm.getStopIndex(), iSeries);
                  hmm.viterbi(subseq);
                  int[] path = hmm.getPath();
                  int step = getStepForLen(wlen);

                  int iWorst = -1;
                  double vWorst = Library.INF;
                  double[] vlist = new double[path.length];
                  for(int i = 0; i < path.length; i++){
                     vlist[i] = hmm.getState(path[i]).eval(subseq.get(i));
                     if (vlist[i] < vWorst){
                        vWorst = vlist[i];
                        iWorst = i;
                     }
                  }

                  for(int i = iWorst; i < vlist.length; i += step)
                     nbl.add(new ScoredWindow(iSeries, tm.getStartIndex() + i, 1, vlist[i]));
                  for(int i = iWorst - step; i >= 0; i -= step)
                     nbl.add(new ScoredWindow(iSeries, tm.getStartIndex() + i, 1, vlist[i]));
               }
               else nbl.add(new ScoredWindow(iSeries, tm.getStartIndex(), 1, score));
            }
         }

         // now center a window on each of them
         for(int i = 0; i < nbl.size(); i++){
            ScoredWindow swin = nbl.get(i);
            int iStart = swin.iStart + (swin.nLength - wlen) / 2;
            iStart = Math.max(iStart, 0);
            iStart = Math.min(iStart, gdata.seqLens[swin.iSeries] - wlen);
            WindowLocation wloc = new WindowLocation(swin.iSeries, iStart, wlen);
            clocs.add(wloc);
         }
      }
      else if (kSearch == METHOD_CLUSTER){ // cluster ////////////////////
         System.err.print("Extracting windows and building tree... ");
         TimerMS timer = new TimerMS();

         if (treeCluster == null){
            enable("cluster", false);

            // extract windows
            int step = getStepForLen(wlen);
            final ArrayList<Sequence> windows = new ArrayList<Sequence>();
            for(int iSeq = 0; iSeq < nSeqs; iSeq++){
               Sequence seq = gdata.tseries.get(iSeq);
               int T = seq.length() - wlen;
               for(int i = 0; i <= T; i += step){
                  Sequence ss = seq.subseq(i, i + wlen, iSeq);
                  windows.add(ss);
               }
            }

            // build tree
            treeCluster = new SeqTree();
            SeqTree.Stop stop = (btrFixedLeaf.isSelected() ? SeqTree.Stop.MaxMembers : SeqTree.Stop.NLeaves);
            int stopParam = lsTree.getValue();
            if (stop == SeqTree.Stop.NLeaves) stopParam = stopParam * 2 - 1;
            SeqTree.SplitOrder split = (btrSplitSize.isSelected() ? SeqTree.SplitOrder.MaxMembers
                  : SeqTree.SplitOrder.Widest);
            if (!treeCluster.build(windows, split, metseq, stop, stopParam)){
               System.err.printf("Error: failed to build sequence tree.\n");
               return null;
            }
            System.err.printf("done (%dms| %d windows, %d nodes).\n", timer.time(), windows.size(),
                  treeCluster.getNumNodes());
         }

         // search leaves for seeds
         final boolean bRandom = btcNodeRandom.isSelected();
         final boolean bPrototype = btcNodeProto.isSelected();
         if (!bRandom && !bPrototype){
            System.err.printf("Error: neither random nor prototype seed selection is checked.");
            return null;
         }
         System.err.print("Searching for seeds in tree leaves... ");
         timer.reset();
         DataTreeApply op = new DataTreeApply() {
            public void apply(DataTreeNode dnode, Object param)
            {
               SeqTreeNode node = (SeqTreeNode)dnode;
               if (node.isInternal()) return;
               int M = node.getNumMembers();

               if (bRandom) clocs.add(node.getMember(Library.random(M)).getWindowLoc());

               if (bPrototype){
                  MyIntList bases = new MyIntList();
                  int N = 2;
                  int ix;
                  if (M > 2 * N){
                     // we have many data points, so approximate the prototype
                     ix = node.findFarthest(Library.random(node.getNumMembers()), metseq);
                     bases.add(ix);
                     for(int i = 0; i < N; i++){
                        ix = node.findFarthest(bases, metseq);
                        bases.add(ix);
                     }
                     ix = node.findClosest(bases, metseq);
                  }
                  else{
                     // we have few data points, so compute the exact prototype
                     ix = -1;
                     double vBest = Library.INF;
                     for(int i = 0; i < M; i++){
                        double v = node.getDistSum(node.getMember(i), metseq);
                        if (v < vBest){
                           vBest = v;
                           ix = i;
                        }
                     }
                  }
                  clocs.add(node.getMember(ix).getWindowLoc());
               }
            }
         };
         treeCluster.apply(op, null);
         System.err.printf("done (%dms).\n", timer.time());
      }
      else if (kSearch == METHOD_IJCAI){ // IJCAI ////////////////////
         if (gdata.chiuView.seeds == null) if (!runIJCAI()) return null;

         enable("ijcai", false);
         ArrayList<ArrayList<WindowLocation>> occList = gdata.chiuView.occList;
         ArrayList<Pair<Sequence, Sequence>> seeds = gdata.chiuView.seeds;

         if (btrIJCAIseeds.isSelected()){
            // convert the IJCAI seeds to candidate locations
            for(Pair<Sequence, Sequence> seed : seeds){
               clocs.add(seed.first.getWindowLoc());
               clocs.add(seed.second.getWindowLoc());
            }
         }
         else{
            String name = getMotifName(iMotif);
            if (hmmsIJCAI == null){
               TimerMS timer = new TimerMS();
               hmmsIJCAI = new ArrayList<AbstractHMM>();
               System.err.printf("Learning %d hmms...", occList.size());
               for(ArrayList<WindowLocation> occs : occList)
                  hmmsIJCAI.add(genHMM(occs, name));
               System.err.printf(" done (%dms).\n", timer.time());
            }
            else{
               // update the name of each HMM to reflect the current motif index
               for(AbstractHMM hmm : hmmsIJCAI)
                  hmm.setName(name);
            }
            globalHmmSeeds.addAll(hmmsIJCAI);
         }
      }
      else if (kSearch == METHOD_KNN_FULL){ // KNN (Full) ////////////////////
      
      }
      else if (kSearch == METHOD_EXHAUSTIVE){ // exhaustive ////////////////////
         int step = lsExhaustStep.getValue();
         for(int iSeq = 0; iSeq < nSeqs; iSeq++){
            int xMax = gdata.seqLens[iSeq] - wlen;
            for(int x = (step + 2) / 4; x <= xMax; x += step)
               clocs.add(new WindowLocation(iSeq, x, wlen));
         }
      }
      else if (kSearch == METHOD_COMBINE){ // mixture ////////////////////
         if (btcRandom.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_RANDOM, labels);
            if (list != null) clocs.addAll(list);
         }
         if (btcUnlikely.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_UNLIKELY, labels);
            if (list != null) clocs.addAll(list);
         }
         if (btcCluster.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_CLUSTER, labels);
            if (list != null) clocs.addAll(list);
         }
         if (btcIJCAI.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_IJCAI, labels);
            if (list != null) clocs.addAll(list);
         }
         if (btcKnnFull.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_KNN_FULL, labels);
            if (list != null) clocs.addAll(list);
         }
         if (btcExhaustive.isSelected()){
            ArrayList<WindowLocation> list = getCandidateLocs(iMotif, METHOD_EXHAUSTIVE, labels);
            if (list != null) clocs.addAll(list);
         }
      }
      return clocs;
   }

   protected boolean runIJCAI()
   {
      int wlen = lsWinLen.getValue();
      int nPAA = lsPAA.getValue();
      int nSAX = lsSAX.getValue();
      boolean bEstR = btcEstR.isSelected();
      double R = 0;
      double fProp = 0;
      if (bEstR){
         try{
            fProp = Double.parseDouble(tfEstRProp.getText()) / 100.0;
         } catch (NumberFormatException e){
            System.err.printf("Error: invalid proportion size (\"%s\")\n", tfEstRProp.getText());
            return false;
         }
      }
      else{
         try{
            R = Double.parseDouble(tfRadius.getText());
         } catch (NumberFormatException e){
            System.err.printf("Error: invalid neighborhood radius (\"%s\")\n", tfRadius.getText());
            return false;
         }
      }

      gdata.chiuView.bVerbose = false;
      TimerMS timer = new TimerMS();
      System.err.print("Running IJCAI algorithm to locate motif seeds... ");
      gdata.chiuView.run(wlen, nPAA, nSAX, bEstR, false, R, fProp, 0, 3);
      System.err.printf("done (%d seeds, %s).\n", gdata.chiuView.seeds.size(), Library.formatDuration(timer
            .time()));

      return true;
   }

   /** @return augment a single location; may contain additional sizes centered on existing locations */
   protected ArrayList<WindowLocation> augmentLocSize(WindowLocation wloc, int iAugSize)
   {
      ArrayList<WindowLocation> wlocsAug = new ArrayList<WindowLocation>();
      wlocsAug.add(wloc);
      if (iAugSize == 0) return wlocsAug; // no augmentation
      if (iAugSize == 2 || iAugSize == 3){ // half + double
         int wlen = Library.getOddHalf(wloc.nLength);
         int dlen = Math.round(wloc.nLength / 4.0f);
         if (wlen > 2){
            int iStart = wloc.iStart + dlen;
            wlocsAug.add(new WindowLocation(wloc.iSeries, iStart, wlen));
         }

         wlen = (int)(wloc.nLength * 1.5f);
         int iStart = wloc.iStart - dlen;
         if (iStart >= 0 && iStart + wlen <= gdata.seqLens[wloc.iSeries])
            wlocsAug.add(new WindowLocation(wloc.iSeries, iStart, wlen));
      }
      if (iAugSize == 1 || iAugSize == 3){ // 3/4 + 5/4
         int wlen = Math.round(wloc.nLength * 0.75f);
         int dlen = Math.round(wloc.nLength / 8.0f);
         if (wlen > 2){
            int iStart = wloc.iStart + dlen;
            wlocsAug.add(new WindowLocation(wloc.iSeries, iStart, wlen));
         }

         wlen = (int)(wloc.nLength * 1.25f);
         int iStart = wloc.iStart - dlen;
         if (iStart >= 0 && iStart + wlen <= gdata.seqLens[wloc.iSeries])
            wlocsAug.add(new WindowLocation(wloc.iSeries, iStart, wlen));
      }
      return wlocsAug;
   }

   /** @return augment list of locations; may contain additional sizes centered on existing locations */
   protected ArrayList<WindowLocation> augmentLocSizes(ArrayList<WindowLocation> clocs, int iAugSize)
   {
      if (iAugSize == 0) return clocs; // no augmentation
      ArrayList<WindowLocation> clocsAug = new ArrayList<WindowLocation>();
      for(WindowLocation wloc : clocs){
         ArrayList<WindowLocation> aug = augmentLocSize(wloc, iAugSize);
         clocsAug.addAll(aug);
      }
      return clocsAug;
   }

   /** @return left-right HMM initialized from the given window locations */
   protected HmmLR genHMM(ArrayList<WindowLocation> occs, String name)
   {
      int nOccs = occs.size();

      // calc length stats for training set
      int minLen = Integer.MAX_VALUE;
      int meanLen = 0;
      for(WindowLocation wloc : occs){
         if (wloc.length() < minLen) minLen = wloc.length();
         meanLen += wloc.length();
      }
      meanLen /= nOccs;
      assert (minLen > 0) : String.format("invalid minimum occurrence length (%d)", minLen);

      // build initial HMM
      int nStatesMax = minLen * 2 - 1;
      int nStates = Math.min(nStatesMax, meanLen);
      if (nStates % 2 == 0) nStates++; // ensure # states is odd
      int nDims = gdata.getNumDims();
      int nSkip = 1;
      HmmLR hmm = new HmmLR(nStates, 1, nDims);
      double edur = (double)meanLen / (nStates - 1);
      double pself = Math.max(1.0 - 1.0 / edur, 0.1);
      hmm.initTran(pself, nSkip);
      hmm.setName(name);

      // build the training set
      ArrayList<Sequence> vtrain = new ArrayList<Sequence>();
      for(WindowLocation wloc : occs)
         vtrain.add(wloc.getSeq(gdata.tseries));

      // learn HMM params
      hmm.init_segk_overlap(vtrain, 0.25);
      double[][] tranOrig = hmm.saveTran();
      hmm.train_bw(vtrain); // TODO BW or Viterbi? 
      hmm.blendTran(tranOrig, 0.5);
      return hmm;
   }

   /** @return left-right HMM initialized from the given window location */
   protected HmmLR genHMM(WindowLocation wloc, String name)
   {
      int nStates = Math.max(3, Math.round(wloc.nLength * 0.75f));
      if (nStates % 2 == 0) nStates++; // ensure # states is odd
      int nDims = gdata.getNumDims();
      int nSkip = 1;
      HmmLR hmm = new HmmLR(nStates, nSkip, nDims);
      double edur = (double)wloc.length() / (nStates - 1);
      double pself = Math.max(1.0 - 1.0 / edur, 0.1);
      hmm.initTran(pself, nSkip);
      hmm.setName(name);
      Sequence seqTrain = wloc.getSeq(gdata.tseries);
      hmm.setVar(gdata.getInitVar());
      hmm.setUpdateVar(false);
      hmm.init_segk_overlap(seqTrain, 0.25);
      // System.err.printf("segk: %s\n", hmm);
      double[][] tranOrig = hmm.saveTran();
      hmm.train_bw(seqTrain);
      // System.err.printf("bw: %s\n", hmm);
      hmm.blendTran(tranOrig, 0.5);
      // System.err.printf("tblend: %s\n", hmm);
      hmm.setUpdateVar(true);
      return hmm;
   }

   /** @return parallel (mixture) HMM initialized from the seed location and size augmentation parameter */
   protected AbstractHMM genMixture(WindowLocation seedLoc, String name, int iAugSize)
   {
      ArrayList<WindowLocation> wlocs = augmentLocSize(seedLoc, iAugSize);
      if (wlocs.size() == 1) return genHMM(seedLoc, name);

      // build each path
      ArrayList<HmmLR> hmmPaths = new ArrayList<HmmLR>();
      int nStates = 0;
      for(WindowLocation wloc : wlocs){
         HmmLR hmm = genHMM(wloc, name);
         nStates += hmm.getNumStates();
         hmmPaths.add(hmm);
      }

      // put the paths together into a single HMM
      HMM hmm = new HMM(nStates, gdata.getNumDims());
      hmm.setName(name);
      Arrays.fill(hmm.getPiStart(), Library.LOG_ZERO);
      Arrays.fill(hmm.getPiEnd(), Library.LOG_ZERO);
      Arrays.fill(hmm.getPiLeave(), Library.LOG_ZERO);
      double[][] tran = hmm.getFullTransMatrix();
      for(int i = 0; i < tran.length; i++)
         Arrays.fill(tran[i], Library.LOG_ZERO);
      int nPaths = hmmPaths.size();
      double piStart = Math.log(1.0 / nPaths);
      int iStateBase = 0;
      for(HmmLR hpath : hmmPaths){
         hmm.setPiStart(iStateBase, piStart);
         double[][] tranLR = hpath.getTranLR();
         int nPathStates = hpath.getNumStates();
         for(int i = 0; i < nPathStates; i++){
            hmm.setState(i + iStateBase, hpath.getState(i));
            for(int j = 0; j < tranLR[i].length; j++)
               tran[iStateBase + i][iStateBase + i + j] = tranLR[i][j];
         }
         iStateBase += nPathStates;
         hmm.setPiEnd(iStateBase - 1, Library.LOG_ONE);
         hmm.setPiLeave(iStateBase - 1, hpath.getPiLeave(nPathStates - 1));
      }

      return hmm;
   }

   /**
    * @return parallel (mixture) HMM with the given HMMs along each path and the background model between them
    */
   protected AbstractHMM genParseMixture(String name, ArrayList<AbstractHMM> hmmPaths, AbstractHMM hmmBG)
   {
      int nDims = gdata.getNumDims();
      int nPaths = hmmPaths.size();
      int nStatesBG = hmmBG.getNumStates();
      int nStates = nStatesBG;
      for(AbstractHMM hmmPath : hmmPaths)
         nStates += hmmPath.getNumStates();

      // create the mixture HMM
      HMM hmm = new HMM(nStates, nDims);

      // setup p(start), p(end), p(leave) from BG model
      double[] pi = hmm.getPiStart();
      Arrays.fill(pi, Library.LOG_ZERO);
      Library.copy(hmmBG.getPiStart(), pi);

      pi = hmm.getPiEnd();
      Arrays.fill(pi, Library.LOG_ZERO);
      Library.copy(hmmBG.getPiEnd(), pi);

      pi = hmm.getPiLeave();
      Arrays.fill(pi, Library.LOG_ZERO);
      Library.copy(hmmBG.getPiLeave(), pi);

      // clear the HMM transition matrix
      double[][] tran = hmm.getFullTransMatrix();
      for(int i = 0; i < nStates; i++)
         Arrays.fill(tran[i], Library.LOG_ZERO);

      // how many states in the BG model are start states?
      int nStartBG = 0;
      for(int i = 0; i < nStatesBG; i++)
         if (hmmBG.getPiStart(i) > Library.LOG_ZERO) nStartBG++;

      // copy the bg transitions and account for leaving
      double[][] tranBG = hmmBG.getFullTransMatrix();
      for(int i = 0; i < nStatesBG; i++){
         if (hmmBG.getPiEnd(i) > Library.LOG_ZERO){
            double pstay = Math.log(1.0 - Math.exp(hmmBG.getPiLeave(i)));
            for(int j = 0; j < nStatesBG; j++)
               tran[i][j] = tranBG[i][j] + pstay;
         }
         else Library.copy(tranBG[i], tran[i]);
      }

      // copy the bg obs dists
      for(int i = 0; i < nStatesBG; i++)
         hmm.getState(i).copyFrom(hmmBG.getState(i));

      // build the paths
      int iBase = nStatesBG; // BG model comes first
      for(int iPath = 0; iPath < nPaths; iPath++){
         AbstractHMM hmmPath = hmmPaths.get(iPath);
         int nStatesPath = hmmPath.getNumStates();

         // we can get from the BG model to this path
         for(int i = 0; i < nStatesBG; i++){
            if (hmmBG.getPiEnd(i) == Library.LOG_ZERO) continue;
            tran[i][iBase] = hmmBG.getPiLeave(i) - Math.log(nPaths);
         }

         // fill in the path transitions
         double[][] tranPath = hmmPath.getFullTransMatrix();
         for(int i = 0; i < nStatesPath; i++)
            Library.copy(tranPath[i], tran[iBase + i], 0, iBase, nStatesPath);

         // fill in the obs dists
         for(int i = 0; i < nStatesPath; i++)
            hmm.getState(iBase + i).copyFrom(hmmPath.getState(i));

         // we can get from this model back to the bg model
         for(int i = 0; i < nStatesPath; i++){
            if (hmmPath.getPiEnd(i) == Library.LOG_ZERO) continue;
            for(int j = 0; j < nStatesBG; j++){
               if (hmmBG.getPiStart(j) == Library.LOG_ZERO) continue;
               tran[iBase + i][j] = hmmPath.getPiLeave(i) + hmmBG.getPiStart(j);
            }
         }

         // normalize intra-path transitions for leave transitions
         for(int i = 0; i < nStatesPath; i++){
            if (hmmPath.getPiEnd(i) == Library.LOG_ZERO) continue;
            double pstay = Math.log(1.0 - Math.exp(hmmPath.getPiLeave(i)));
            for(int j = 0; j < nStatesPath; j++){
               tran[iBase + i][iBase + j] += pstay;
            }
         }

         iBase += nStatesPath;
      }

      return hmm;
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
      hmms.addAll(motifHMMs);
      ContRecRet crrBest = curCRR;
      if (crrBest == null){
         TimerMS timer = new TimerMS();
         crrBest = HTK.contRec(hmms, gdata.tseries);
         if (crrBest == null){
            System.err.println("Error: HTK / HVite failed");
            return false;
         }
         motifLabels = crrBest.markupSets;
         listDataLL.add(-crrBest.llTotal / lambda);
         listAIC.add(calcAIC(crrBest.llTotal, hmms));
         listBIC.add(calcBIC(crrBest.llTotal, hmms));
         listThreshMDL.add(calcMDL(crrBest.llTotal, hmms));
      }

      // get the base candidate locations
      ArrayList<AbstractHMM> hmmSeeds = new ArrayList<AbstractHMM>();
      ArrayList<WindowLocation> clocs = getCandidateLocs(iMotif, cpSearch.getSelectedIndex(),
            crrBest.markupSets);
      if (clocs == null) return false;
      int nSeeds = clocs.size() + hmmSeeds.size();
      if (nSeeds > 0){
         System.err.printf("Searching %d candidate locations + %d HMMs (%d total)  ", clocs.size(), hmmSeeds
               .size(), nSeeds);

         // build the models
         if (!clocs.isEmpty()){
            boolean bSelectSize = btrOneLen.isSelected();
            int iAugSize = cbScale.getSelectedIndex();
            if (bSelectSize){
               // augment base locations with larger/smaller locatiosn if requested
               clocs = augmentLocSizes(clocs, iAugSize);
               for(WindowLocation wloc : clocs)
                  hmmSeeds.add(genHMM(wloc, getMotifName(iMotif)));

               System.err.printf("(%d total models).\n", hmmSeeds.size());
            }
            else{
               for(WindowLocation wloc : clocs)
                  hmmSeeds.add(genMixture(wloc, getMotifName(iMotif), iAugSize));
               System.err.printf("(%d total mixture models).\n", hmmSeeds.size());
            }
         }
         else System.err.println();
         nSeeds = hmmSeeds.size();

         // score each window location
         String sTmpDir = System.getProperty("java.io.tmpdir");
         HtkSetupInfo hsi = null;
         int iBest = -1;
         TimerMS timer = new TimerMS();
         for(int iSeed = 0; iSeed < nSeeds; iSeed++){
            System.err.printf(" Testing candidate HMM %d / %d ... ", iSeed + 1, nSeeds);
            timer.reset();

            // run the continuous recognizer and save positive results
            AbstractHMM hmm = hmmSeeds.get(iSeed);
            if (hsi == null){
               hmms.add(hmm);
               hsi = HTK.setupContRec(null, hmms, gdata.tseries);
               assert (hsi != null);
            }
            else{
               // need to replace the last HMM with a new candidate
               hmms.set(hmms.size() - 1, hmm);
               File fHmm = new File(sTmpDir, hmm.getName());
               assert (fHmm.exists()) : String.format("HMM file should exist (%s)", fHmm.getAbsolutePath());
               if (!HtkHmm.save(hmm, fHmm)) return false;
            }
            ContRecRet crr = HTK.contRecAfterSetup(hsi);
            if (crr == null){
               System.err.printf("done (%dms) ERROR\n", timer.time());
               assert (false);
            }
            else{
               if (crrBest == null || crr.llTotal > crrBest.llTotal){
                  crrBest = crr;
                  iBest = iSeed;
               }
               System.err.printf("done (%dms)  (ll=%.1f)\n", timer.time(), crr.llTotal);
            }
         }
         if (iBest < 0) return false;

         // initialize new motif from best window location
         AbstractHMM hmmBest = hmmSeeds.get(iBest);
         motifHMMs.add(hmmBest);
         motifLabels = crrBest.markupSets;
         curCRR = crrBest;

         // remove this HMM from global repositories
         if (hmmsIJCAI != null) hmmsIJCAI.remove(hmmBest);
      }
      else System.err.println("Warning: no more candidate seeds.");

      listDataLL.add(-crrBest.llTotal / lambda);
      listAIC.add(calcAIC(crrBest.llTotal, hmms));
      listBIC.add(calcBIC(crrBest.llTotal, hmms));
      listThreshMDL.add(calcMDL(crrBest.llTotal, hmms));
      sgc.updateGraph();
      lsMotif.setText(String.format("Motif: %%d / %d", motifHMMs.size()));
      lsMotif.setValues(motifHMMs.size(), 0, motifHMMs.size());
      lsMotif.setEnabled(true);
      btRemoveLast.setEnabled(true);
      btReset.setEnabled(true);
      btrOptAcc.setEnabled(true);
      btrOptRecall.setEnabled(true);
      btScore.setEnabled(true);

      return true;
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

   protected double calcAIC(double loglik, ArrayList<AbstractHMM> hmms)
   {
      int d = getNumParams(hmms);
      int n = gdata.nTotalSeriesLength;
      return -loglik / lambda + d + (double)d * (d + 1) / (n - d - 1);
   }

   protected double calcBIC(double loglik, ArrayList<AbstractHMM> hmms)
   {
      return -loglik / lambda + (hmms.size() + getNumParams(hmms)) * Math.log(gdata.nTotalSeriesLength)
            / 2.0;
   }

   protected double calcMDL(double loglik, ArrayList<AbstractHMM> hmms)
   {
      int nDims = gdata.getNumDims();
      double modelBits = -Math.log(hmms.size());
      for(AbstractHMM hmm : hmms){
         int nStates = hmm.getNumStates();
         int nTran = 2 * nStates - 3;
         int nParam = nDims * nStates;// * 2; // TODO with fixed variance, we don't need to send var param
         modelBits += 18 * (nTran + nParam); // ~18 compressed bits per 32-bit float
         // TODO encode at certain precision and calc entropy
      }

      return -loglik + modelBits;
   }

   protected void highlight(int iMotif, ArrayList<MarkupSet> labels)
   {
      if (iMotif == iHighlightMotif) return;
      iHighlightMotif = iMotif;
      int nSeqs = gdata.getNumSeqs();

      // clear the highlights
      for(int i = 0; i < nSeqs; i++)
         crbGT[i].clearHighlight();

      // highlight a real motif
      if (iMotif >= 0){
         String sMotifName = getMotifName(iMotif + 1);
         assert(nSeqs == labels.size()) : String.format("mismatch: nseqs=%d  labels.size=%d", nSeqs, labels.size());
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

   /** run the continuous recognizer with the current motif set and save the results */
   protected boolean runContRec()
   {
      // build list of all bg, current motifs, and this candidate HMM
      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      hmms.add(hmmBG);
      hmms.addAll(motifHMMs);
      curCRR = HTK.contRec(hmms, gdata.tseries);
      if (curCRR == null) return false;
      motifLabels = curCRR.markupSets;
      return true;
   }

   /** @return name of motif given an index */
   protected String getMotifName(int ix)
   {
      return String.format("motif%02d", ix);
   }

   /** @return index of motif given a name */
   protected int getMotifIndex(String name)
   {
      if (name.equals(BG)) return -1;
      int n = motifHMMs.size();
      for(int i = 0; i < n; i++)
         if (motifHMMs.get(i).getName().equals(name)) return i;

      // if we can't find an HMM, try to extract a number
      try{
         return Integer.parseInt(name.substring(name.length() - 2));
      } catch (NumberFormatException e){
         return -1;
      }
   }

   protected AbstractHMM getMotifFromName(String name)
   {
      int ix = getMotifIndex(name);
      if (ix < 0) return hmmBG;
      return motifHMMs.get(ix);
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
      System.err.printf("#spots: %d\n", nspots);

      // evaluate the results
      int[][] confm = cri.scoreWordSpotWordsFast(gdata.labData);

      // calc best class permutation and remap spots
      System.err.print("Calculating best cluster <-> class mapping... ");
      boolean bOptAcc = btrOptAcc.isSelected();
      short[] map = SupTest.calcBestMappingFromConfMatrix(confm, bOptAcc);
      cri.remapSpotsOld(map);
      System.err.println("done.");

      // recalc stats with permuted classes
      System.err.print("Recomputing performance stats... ");
      TimerMS timer = new TimerMS();
      confm = cri.scoreWordSpotWordsFull(gdata.labData);
      cri.scoreWordSpotFrames(gdata.tseries, gdata.labData);
      System.err.printf("done (%dms).\n", timer.time());

      SupTest.dumpConfMatrix(confm, "Confusion Matrix", gdata.classes, gdata.nExamples, false);

      // cri.dump();

      assert (cri.isValid(true)) : "invalid continuous recognition stats";
      // System.err.println(cri.getResultsString(gdata.classes));
   }

   /**
    * build a composite HMM with multiple, parallel paths and then run EM to fit each path to a "segment" of
    * the data
    */
   public void segmentEM()
   {
      System.err.print("Building the segmentation HMM... ");
      TimerMS timer = new TimerMS();

      // build the mixture model
      int eLength = lsSegExpLength.getValue();
      int nPaths = lsSegMixSize.getValue();
      int nPathStates = lsSegCompLen.getValue();
      int nDims = gdata.getNumDims();
      int nTotalStates = nPaths * nPathStates + 2;
      HMM hmm = new HMM(nTotalStates, nDims);

      // setup start, end, leave priors
      double[] pi = hmm.getPiStart();
      Arrays.fill(pi, Library.LOG_ZERO);
      pi[0] = Library.LOG_ONE;

      pi = hmm.getPiEnd();
      Arrays.fill(pi, Library.LOG_ZERO);
      pi[nTotalStates - 1] = Library.LOG_ONE;

      pi = hmm.getPiLeave();
      Arrays.fill(pi, Library.LOG_ZERO);
      pi[nTotalStates - 1] = Math.log(0.1);

      // build the transition matrix
      double edur = (double)eLength / (nPathStates - 1);
      double pDummySelf = 0.75;
      double pself = Math.max(1.0 - 1.0 / edur, 0.1);
      double llself = Math.log(pself);
      double llnext1 = Math.log(1.0 - pself);
      double llnext12 = Math.log((1.0 - pself) * 2 / 3);
      double llnext22 = Math.log((1.0 - pself) / 3);
      double[][] tran = hmm.getFullTransMatrix();
      for(int i = 0; i < tran.length; i++)
         Arrays.fill(tran[i], Library.LOG_ZERO);
      int iBase = 1;
      for(int iPath = 0; iPath < nPaths; iPath++){
         tran[0][iBase] = Math.log((1.0 - pDummySelf) / nPaths); // start node to this path
         for(int i = 0; i < nPathStates; i++){
            tran[iBase + i][iBase + i] = llself;
            if (i + 1 < nPathStates){
               if (i + 2 < nPathStates){
                  tran[iBase + i][iBase + i + 1] = llnext12;
                  tran[iBase + i][iBase + i + 2] = llnext22;
               }
               else tran[iBase + i][iBase + i + 1] = llnext1;
            }
         }
         iBase += nPathStates;
         tran[iBase - 1][nTotalStates - 1] = llnext1; // this path to end node
      }
      tran[0][0] = Math.log(pDummySelf);
      tran[nTotalStates - 1][0] = Library.LOG_ONE;

      // fill in the observation distributions
      GaussianDiagonal gd = (GaussianDiagonal)hmm.getState(0);
      gd.setMean(gdata.getGlobalMean());
      gd.setVar(gdata.getGlobalVar());
      gd = (GaussianDiagonal)hmm.getState(nTotalStates - 1);
      gd.setMean(gdata.getGlobalMean());
      gd.setVar(gdata.getGlobalVar());

      GaussianDiagonal dgNoise = new GaussianDiagonal(nDims);
      dgNoise.setMean(FeatureVec.zeros(nDims));
      dgNoise.setVar(gdata.getInitVar());

      iBase = 1;
      for(int iPath = 0; iPath < nPaths; iPath++){
         for(int i = 0; i < nPathStates; i++){
            gd = (GaussianDiagonal)hmm.getState(iBase + i);
            gd.setMean(gdata.getGlobalMean()._add(dgNoise.sample()));
            gd.setVar(gdata.getGlobalVar());
         }
         iBase += nPathStates;
      }
      System.err.printf("done (%dms).\n", timer.time());

      // train the composite model
      System.err.print("Training composite model via EM... ");
      timer.reset();
      hmm.train_bw(gdata.tseries);
      timer.mark();
      System.err.printf("done (%s) (%dms).\n", Library.formatDuration(timer.time()), timer.time());

      // build state map to convert composite states to component index
      HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
      map.put(0, -1);
      map.put(nTotalStates - 1, -1);
      iBase = 1;
      for(int i = 0; i < nPaths; i++){
         for(int j = 0; j < nPathStates; j++)
            map.put(iBase + j, i);
         iBase += nPathStates;
      }

      // fit the model to the data
      segLabels = new ArrayList<MarkupSet>();
      int nSeqs = gdata.getNumSeqs();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         MarkupSet marks = new MarkupSet();
         segLabels.add(marks);
         double loglik = hmm.viterbi(gdata.tseries.get(iSeq));
         int[] path = hmm.getPath();

         // map the viterbi path to components of the mixture model
         int iCurComp = -1;
         int iStart = -1;
         for(int i = 0; i < path.length; i++){
            int iComp = map.get(path[i]);
            if (iComp < 0) continue; // background
            if (iComp != iCurComp){
               // save the previous run
               if (iCurComp >= 0)
                  marks.add(new TimeMarker(getMotifName(iCurComp + 1), TimeMarker.Units.Index, iStart, i));

               // start the new run
               iStart = i;
               iCurComp = iComp;
            }
         }
         // don't forget the last segment
         if (iCurComp >= 0)
            marks.add(new TimeMarker(getMotifName(iCurComp), TimeMarker.Units.Index, iStart, path.length));
      }

      scoreMotifs(segLabels);

      lsSegMotif.setText(String.format("Motif: %%d / %d", nPaths));
      lsSegMotif.setValues(1, 0, nPaths);
      lsSegMotif.setEnabled(true);
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();

      if (src == btNextMotif){
         TimerMS timer = new TimerMS();
         int N = lsFindNMotifs.getValue();
         for(int i = 0; i < N; i++)
            if (!findNext()) break;
         timer.mark();
         System.err.printf("Total time for discovery: %s (%dms)\n", Library.formatDuration(timer.time()),
               timer.time());
      }
      else if (src == cbWinLen){
         switch(cbWinLen.getSelectedIndex()){
            case 0: // exercise
               lsWinLen.setValue(15);
               break;
            case 1: // tidigits
               lsWinLen.setValue(22);
               break;
            case 2: // custom
               // do nothing
               break;
         }
      }
      else if (src == btReset) reset();
      else if (src == btRemoveLast){
         int iRemove = motifHMMs.size() - 1;
         if (iRemove == 0) reset();
         else{
            motifHMMs.remove(iRemove);
            if (motifHMMs.isEmpty()) reset();
            else{
               runContRec();
               iHighlightMotif = -1; // make sure we redraw the highlights
               lsMotif.setValues(motifHMMs.size(), 0, motifHMMs.size());
               lsMotif.setText(String.format("Motif: %%d / %d", motifHMMs.size()));
            }
            listDataLL.removeLast();
            listAIC.removeLast();
            listBIC.removeLast();
            listThreshMDL.removeLast();
            sgc.updateGraph();
         }
      }
      else if (src == btScore) scoreMotifs(motifLabels);
      else if (src == btrFixedLeaf){
         nTotalLeaves = lsTree.getValue();
         lsTree.setValue(nMaxLeafSize);
         lsTree.setText("Max # Leaf Members: %d");
      }
      else if (src == btrNumLeaves){
         nMaxLeafSize = lsTree.getValue();
         lsTree.setValue(nTotalLeaves);
         lsTree.setText("Num Leaves: %d");
      }
      else if (src == btSegRun) segmentEM();
      else if (src == btcEstR){
         boolean bEstR = btcEstR.isSelected();
         tfRadius.setEnabled(!bEstR);
         btRealRadius.setEnabled(!bEstR);
         tfEstRProp.setEnabled(bEstR);
      }
      else if (src == btRealRadius) gdata.chiuView.dumpRealR();
      else if (src == btAddIJCAI){
         if (!runIJCAI()) return;
         enable("ijcai", false);
         ArrayList<ArrayList<WindowLocation>> occList = gdata.chiuView.occList;
         int nSeeds = occList.size();
         TimerMS timer = new TimerMS();
         System.err.printf("Learning %d hmms...", nSeeds);
         for(int iSeed = 0; iSeed < nSeeds; iSeed++){
            String name = getMotifName(iSeed + 1);
            ArrayList<WindowLocation> occs = occList.get(iSeed);
            motifHMMs.add(genHMM(occs, name));
         }
         System.err.printf(" done (%dms).\n", timer.time());
         cpSearch.setSelectedIndex(METHOD_NONE);
      }
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();

      if (src == lsMotif){
         int iMotif = lsMotif.getValue();
         highlight(iMotif - 1, motifLabels);
         if (iMotif == 0) tfMotif.setText("");
         else{
            int nOccs = 0;
            String sMotifName = getMotifName(iMotif);
            for(MarkupSet marks : motifLabels){
               for(TimeMarker tm : marks.getList()){
                  if (tm.getTag() == null || tm.getTag().equals(BG)) continue;
                  if (!tm.getTag().equals(sMotifName)) continue;
                  nOccs++;
               }
            }
            tfMotif.setText(String.format("Motif) %d occs", nOccs));
         }
      }
      else if (src == lsSegMotif){
         int iMotif = lsSegMotif.getValue();
         highlight(iMotif - 1, segLabels);
      }
      else if (src == lsWinLen){
         cbWinLen.setSelectedIndex(2); // select "custom"
      }
   }

   public static void main(String args[])
   {
      Runtime rt = Runtime.getRuntime();
      System.err.printf("Java Version: %s   VM: %s (%s)\n", System.getProperty("java.version"), System
            .getProperty("java.vm.version"), System.getProperty("java.vm.name"));
      System.err.printf("Memory (free/total/max): %dk / %dk / %dk\n", rt.freeMemory() / 1024, rt
            .totalMemory() / 1024, rt.maxMemory() / 1024);

      String sDefFile = "/home/dminn/research/kdm/data/ah-data/exercise/sub8/accgyr.def";
      // String sDefFile = "E:\\research\\gatech-discovery\\exercise-data\\sub8\\accgyr.def";
      //String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/dim13/ae.def";
      //String sDefFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/ae.def";
      //String sDefFile = "/home/dminn/research/kdm/data/asl-joseph/asl16D.def";

      int c;
      LongOpt[] longopts = new LongOpt[] { new LongOpt("help", LongOpt.NO_ARGUMENT, null, 1001) };

      Getopt g = new Getopt("DiscApp", args, "?", longopts, true);
      while((c = g.getopt()) != -1){
         switch(c){
         case '?':
         case 1001: // help
            usage();
            System.exit(1);
            break;
         }
      }

      if (g.getOptind() < args.length) sDefFile = args[g.getOptind()];

      ToolTipManager.sharedInstance().setInitialDelay(500);
      DiscAAAI da = new DiscAAAI();
      if (!da.setup(sDefFile)){
         System.err.println("DiscAAAI setup failed.");
         System.exit(1);
      }
   }

}
