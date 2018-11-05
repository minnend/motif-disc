package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

import javax.swing.filechooser.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.io.*;
import mdisc.io.*;
import kdm.mlpr.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.util.*;
import kdm.tools.*;
import kdm.tools.SupTest.*;
import mdisc.data.*;

/**
 * View of the results of word spotting with models built from labeled data
 */
public class WordSpottingView extends JPanel implements ActionListener, ChangeListener,
      ScoreGraphListener, MouseListener
{
   public static final int nModelTypes = 3;

   protected GlobalData gdata;
   protected ArrayList<Sequence> tseries;
   protected ArrayList<DiscreteSeq> qseries;
   protected int nSymbols;
   protected JPanel pModel;
   protected TreeMap<String, ArrayList<Sequence>> curData, labData;
   protected TreeMap<String, ArrayList<DiscreteSeq>> curQData, labQData;
   protected HashMap<String, TreeMap<String, ArrayList<Sequence>>> labLoadData;
   protected HashMap<String, TreeMap<String, ArrayList<DiscreteSeq>>> labLoadQData;
   protected ArrayList<MotifInfo> seedMotifs;
   protected ColoredRangeBar[] crbSpot;
   protected ScoreGraphComplex sgc;
   protected Color[] colors;
   protected ComboPanel cpModel;
   protected JComboBox cbScore, cbHmmTrain, cbHmmEval, cbVizModel;
   protected JComboBox cbSeeds, cbModels;
   protected JCheckBox btNorm;
   protected JTextField tfThresh;
   protected JTextArea taResults;
   protected JSpinner spinMinLen, spinMaxLen, spinMinSpots, spinMaxSpots;
   protected JSpinner spinHmmStates, spinHmmSkip;
   protected JButton btSpot, btLoadCRI, btSaveCRI, btSaveLabels, btLoadTrain, btShowTrain, btVizTrain;
   protected JButton btTrain, btTrainGarbage, btLoadModels, btSaveModels, btSaveHtk;
   protected VertChooseContainer vcc;
   protected ProbSeqModel models[][]; // [iModel][iClass]
   protected GaussianDyn1D pScore[][]; // [iModel][iClass]
   protected ContRecInfo cri;
   protected int nDisp = -1;
   protected HashMap<String, TrainInfo> modelCache;
   protected MatrixViz mviz;
   protected StateTopoViz topoViz;
   protected ParallelCoordinateView pcviz;
   protected File fDataSrc;

   public WordSpottingView(GlobalData gdata)
   {
      super(new BorderLayout());
      gdata.wordSpottingView = this;
      setOpaque(true);

      this.gdata = gdata;
      tseries = gdata.tseries;
      qseries = gdata.qseries;
      nSymbols = gdata.nSymbols;
      curData = labData = gdata.labData;
      curQData = labQData = gdata.labQData;
      modelCache = new HashMap<String, TrainInfo>();
      labLoadData = new HashMap<String, TreeMap<String, ArrayList<Sequence>>>();
      labLoadQData = new HashMap<String, TreeMap<String, ArrayList<DiscreteSeq>>>();

      int nClasses = labData.size();
      colors = Library.generateColors(nClasses);
      selectLabData();

      JPanel leftp = new JPanel(new BorderLayout());

      // create the center view
      crbSpot = new ColoredRangeBar[tseries.size()];
      JPanel seqsp = new JPanel(new VerticalLayout());
      seqsp.setBackground(Color.gray);
      for(int iSeq = 0; iSeq < tseries.size(); iSeq++){
         Sequence seq = tseries.get(iSeq);
         JPanel p = new JPanel(new GridFillLayout(2, 1, 0, 0));
         p.setBackground(Color.darkGray);
         p.setBorder(BorderFactory.createMatteBorder(1, 0, 1, 0, Color.darkGray));

         // create the labeled bars
         ColoredRangeBar crb = new ColoredRangeBar(new Range(0, seq.length() - 1));
         crb.setBorder(BorderFactory.createMatteBorder(1, 0, 1, 0, crb.getBackground()));
         crb.setToolTipText(String.format("Series %d: %s", iSeq + 1, gdata.getSentence(iSeq)));
         p.add(crb);

         // create the spotted bars
         crbSpot[iSeq] = new ColoredRangeBar(new Range(0, seq.length() - 1));
         crbSpot[iSeq].setBorder(BorderFactory.createMatteBorder(1, 0, 1, 0, crb.getBackground()));
         p.add(crbSpot[iSeq]);
         seqsp.add(p);
         if (iSeq + 1 < tseries.size()) seqsp.add(Box.createVerticalStrut(1));

         // add labeled ranges to the bar
         Iterator<String> itClass = labData.keySet().iterator();
         for(int iClass = 0; iClass < labData.size(); iClass++){
            String sClass = itClass.next();
            ArrayList<Sequence> subs = labData.get(sClass);
            for(Sequence sub : subs){
               if (sub.getParentIndex() != iSeq) continue;
               Range r = new Range(sub.getParentOffset(), sub.getParentOffset() + sub.length() - 1);
               r.payload = new Pair(sClass, colors[iClass].darker().darker());
               crb.add(r);
            }
         }
      }
      leftp.add(new JScrollPane(seqsp), BorderLayout.CENTER);

      // create the bottom view (word spot)
      sgc = new ScoreGraphComplex();
      sgc.addScoreGraphListener(this);
      leftp.add(sgc, BorderLayout.SOUTH);

      // create the right control panel
      vcc = new VertChooseContainer();

      JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftp, vcc);
      split.setDividerLocation(Math.max(gdata.frame.getWidth() / 2, gdata.frame.getWidth() - 260));
      split.setResizeWeight(1.0);
      split.setDividerSize(7);
      split.setBorder(null);
      this.add(split, BorderLayout.CENTER);

      cbHmmTrain = new JComboBox();
      cbHmmTrain.addItem("Baum-Welch");
      cbHmmTrain.addItem("Viterbi");
      cbHmmTrain.addActionListener(this);

      cbHmmEval = new JComboBox();
      cbHmmEval.addItem("Forward");
      cbHmmEval.addItem("Viterbi");
      cbHmmEval.addActionListener(this);

      cbScore = new JComboBox();
      cbScore.addItem(" None (raw)");
      cbScore.addItem(" Sigmoid / Var");
      cbScore.addActionListener(this);

      spinHmmStates = new JSpinner(new SpinnerNumberModel(8, 1, 25, 1));
      spinHmmStates.addChangeListener(this);

      spinHmmSkip = new JSpinner(new SpinnerNumberModel(1, 0, 4, 1));
      spinHmmSkip.addChangeListener(this);
      spinHmmSkip.setEnabled(false); // not necessary with new HMM init method

      // create storage for the models and normalization data
      spinMinLen = new JSpinner(new SpinnerNumberModel(4, 0, 12, 2));
      spinMinLen.addChangeListener(this);
      spinMaxLen = new JSpinner(new SpinnerNumberModel(32, 4, 128, 4));
      spinMaxLen.addChangeListener(this);
      spinMinSpots = new JSpinner(new SpinnerNumberModel(50, 0, 2000, 10));
      spinMinSpots.addChangeListener(this);
      spinMaxSpots = new JSpinner(new SpinnerNumberModel(1200, 0, 5000, 10));
      spinMaxSpots.addChangeListener(this);

      btNorm = new JCheckBox("Normalize", true);
      btNorm.addActionListener(this);

      cbModels = new JComboBox();
      cbModels.addItem("Train Models...");
      cbModels.addActionListener(this);

      cbSeeds = new JComboBox();
      cbSeeds.addItem("Labeled Data");
      cbSeeds.addActionListener(this);

      btLoadModels = new JButton("Load");
      btLoadModels.setToolTipText("Load models");
      btLoadModels.addActionListener(this);
      btSaveModels = new JButton("Save");
      btSaveModels.setToolTipText("Save models (train if necessary)");
      btSaveModels.addActionListener(this);
      btSaveHtk = new JButton("HTK");
      btSaveHtk.setToolTipText("Save HMMs in HTK format");
      btSaveHtk.addActionListener(this);

      tfThresh = new JTextField(String.format("%.3f", SupTest.spotThresh), 6);

      btSpot = new JButton("  Spot Words!  ");
      btSpot.setToolTipText("Perform word spotting");
      btSpot.addActionListener(this);
      btSaveCRI = new JButton("Save...");
      btSaveCRI.setToolTipText("Save Continuous Rec Info");
      btSaveCRI.setEnabled(false);
      btSaveCRI.addActionListener(this);
      btLoadCRI = new JButton("Load...");
      btLoadCRI.setToolTipText("Load Continuous Rec Info");
      btLoadCRI.addActionListener(this);
      btSaveLabels = new JButton("Save...");
      btSaveLabels.setToolTipText("Save Spots in Label Files");
      btSaveLabels.setEnabled(false);
      btSaveLabels.addActionListener(this);
      btLoadTrain = new JButton("Load");
      btLoadTrain.setToolTipText("Load Training Data (Seeds)");
      btLoadTrain.addActionListener(this);
      btShowTrain = new JButton("Show");
      btShowTrain.setToolTipText("Display Training Data Locations");
      btShowTrain.addActionListener(this);
      btTrain = new JButton("Train");
      btTrain.setToolTipText("Train models using current settings");
      btTrain.addActionListener(this);
      btTrainGarbage = new JButton("Garbage Model");
      btTrainGarbage.setToolTipText("Train garbage model using current settings");
      btTrainGarbage.addActionListener(this);
      btVizTrain = new JButton("Viz");
      btVizTrain.setToolTipText("Visualize training data and models");
      btVizTrain.addActionListener(this);

      JPanel p, q;

      // setup operations tab
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout());
      rightp.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 4));
      vcc.addPane("Operations", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      rightp.add(Box.createVerticalStrut(4));
      rightp.add(new JLabel("Training Data...", JLabel.CENTER));
      q = new JPanel();
      q.add(btLoadTrain);
      q.add(btShowTrain);
      q.add(btVizTrain);
      rightp.add(q);
      rightp.add(Box.createVerticalStrut(6));
      q = new JPanel(new BorderLayout());
      q.add(new JLabel("Seeds: "), BorderLayout.WEST);
      q.add(cbSeeds, BorderLayout.CENTER);
      rightp.add(q);
      rightp.add(Box.createVerticalStrut(8));
      q = new JPanel(new BorderLayout());
      q.add(new JLabel("Models: "), BorderLayout.WEST);
      q.add(cbModels, BorderLayout.CENTER);
      rightp.add(q);
      rightp.add(Box.createVerticalStrut(8));
      p = new JPanel();
      p.add(btLoadModels);
      p.add(btSaveModels);
      p.add(btSaveHtk);
      rightp.add(p);
      q = new JPanel();
      q.add(btTrain);
      q.add(btTrainGarbage);
      rightp.add(q);
      rightp.add(btSpot);
      rightp.add(Box.createVerticalStrut(8));

      taResults = new JTextArea(18, 10);
      taResults.setEditable(false);
      rightp.add(new JScrollPane(taResults));

      rightp.add(Box.createVerticalStrut(8));

      GridFlexLayout gfl = new GridFlexLayout(1, 3, 4, 0);
      gfl.setColumn(0, GridFlexLayout.Style.pref);
      gfl.setRow(0, GridFlexLayout.Style.pref);
      p = new JPanel(gfl);
      p.add(new JLabel("CRI:"));
      p.add(btLoadCRI);
      p.add(btSaveCRI);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Labels: "));
      p.add(btSaveLabels);
      rightp.add(p);

      rightp.add(Box.createVerticalStrut(12));

      // create a class-color legend
      if (!labData.isEmpty()){
         rightp.add(new Divider("Legend"));
         rightp.add(Box.createVerticalStrut(8));
         gfl = new GridFlexLayout(labData.size(), 3, 4, 2);
         gfl.setColumn(0, GridFlexLayout.Style.fixed, 10);
         gfl.setColumn(1, GridFlexLayout.Style.pref);
         gfl.setColumn(2, GridFlexLayout.Style.fixed, 60);
         gfl.setRows(GridFlexLayout.Style.pref);
         JPanel legend = new JPanel(gfl);
         rightp.add(legend);
         Iterator<String> itClass = labData.keySet().iterator();
         for(int iClass = 0; iClass < labData.size(); iClass++){
            legend.add(new EmptyComponent());
            legend.add(new JLabel(itClass.next() + ":", JLabel.RIGHT));
            JComponent comp = new JPanel();
            comp.setBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.darkGray));
            comp.setBackground(colors[iClass]);
            legend.add(comp);
         }
         rightp.add(Box.createVerticalStrut(8));
      }

      // create model param pane
      rightp = new VerticalScrollPanel(new VerticalLayout(250));
      rightp.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 4));
      vcc.addPane("Model Settings", new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));

      pModel = new JPanel(new CardLayout());
      rightp.add(pModel);
      rightp.add(Box.createVerticalStrut(12));

      JPanel pHmm = new JPanel(new VerticalLayout());
      q = new JPanel();
      q.add(new JLabel("# States: "));
      q.add(spinHmmStates);
      pHmm.add(q);
      q = new JPanel();
      q.add(new JLabel("# Skips: "));
      q.add(spinHmmSkip);
      pHmm.add(q);
      q = new JPanel();
      q.add(new JLabel("Train: "));
      q.add(cbHmmTrain);
      pHmm.add(q);
      q = new JPanel();
      q.add(new JLabel("Eval: "));
      q.add(cbHmmEval);
      pHmm.add(q);

      cpModel = new ComboPanel("Model: ");
      cpModel.addPanel("Oates", null); // no parameters for Oates Model
      cpModel.addPanel("HMM (LR,diag Gauss)", pHmm);
      cpModel.addPanel("HMM (LR,discrete)", pHmm);
      cpModel.addActionListener(this);
      rightp.add(cpModel);
      rightp.add(Box.createVerticalStrut(16));

      p = new JPanel();
      p.add(new JLabel("Score Method: "));
      p.add(cbScore);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Min Word Len: "));
      p.add(spinMinLen);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Max Word Len: "));
      p.add(spinMaxLen);
      rightp.add(p);

      p = new JPanel();
      p.add(btNorm);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Threshold: "));
      p.add(tfThresh);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Min Number: "));
      p.add(spinMinSpots);
      rightp.add(p);

      p = new JPanel();
      p.add(new JLabel("Max Number: "));
      p.add(spinMaxSpots);
      rightp.add(p);

      // change default view
      cpModel.setSelectedIndex(0);
      cbScore.setSelectedIndex(1);
      cbHmmTrain.setSelectedIndex(1);
      cbHmmEval.setSelectedIndex(1);
      btNorm.setSelected(false);
   }

   public int getCurModelType()
   {
      return cpModel.getSelectedIndex();
   }

   /** @return true if the given model handles discrete data */
   public static boolean isDiscrete(int iModel)
   {
      return (iModel == 2);
   }

   /** @return true if the given model is a kind of HMM */
   public static boolean isHMM(int iModel)
   {
      return (iModel == 1 || iModel == 2);
   }

   /**
    * Create the pscore gaussian for the given model type
    * 
    * @param iModel model type to use
    * @return true if successful
    */
   protected boolean calcPScores(int iModel)
   {
      // build distribution over training data scores
      Iterator<String> itClass = curData.keySet().iterator();
      for(int iClass = 0; iClass < curData.size(); iClass++){
         pScore[iModel][iClass] = new GaussianDyn1D();
         if (isDiscrete(iModel)){
            for(DiscreteSeq dseq : labQData.get(itClass.next())){
               double x = models[iModel][iClass].eval(dseq);
               if (btNorm.isSelected()) x /= dseq.length();
               pScore[iModel][iClass].add(x, false);
            }
         }
         else{
            for(Sequence seq : curData.get(itClass.next())){
               double x = models[iModel][iClass].eval(seq);
               if (btNorm.isSelected()) x /= seq.length();
               pScore[iModel][iClass].add(x, false);
            }
         }
         pScore[iModel][iClass].update();
      }
      return true;
   }

   /**
    * Visualize continuous rec info
    * 
    * @param cri the results to visualize
    * @param bUpdateWSView update the word spot view?
    */
   protected void displayResults(ContRecInfo cri, boolean bUpdateWSView)
   {
      if (bUpdateWSView){
         MyDoubleList scores = new MyDoubleList();
         for(WordSpot spot : cri.getSpots()) scores.add(spot.score);
         sgc.addCurve(String.format("CRI-%03d", sgc.getNumCurves()+1), scores);
      }
      for(int i = 0; i < crbSpot.length; i++)
         crbSpot[i].clear();
      int nSpots = Math.min(cri.getNumSpots(), nDisp < 0 ? Integer.MAX_VALUE : nDisp);
      for(int iSpot = 0; iSpot < nSpots; iSpot++){
         WordSpot wspot = cri.getSpot(iSpot);
         Range r = wspot.getRange();
         r.payload = new Pair(gdata.classes[wspot.iClass], colors[wspot.iClass]);
         crbSpot[wspot.iSeries].add(r);
      }
      if (cri.isValid(false)){
         taResults.setText(cri.getResultsString(gdata.classes) + "\n" + cri.getDumpString());
         taResults.setCaretPosition(0);
      }
   }

   /**
    * Setup SupTest parameters based on cached models
    * 
    * @param sName name of model to use (key for hash look-up)
    * @return true if model found and successfully "loaded"
    */
   protected boolean setupFromCache(String sName)
   {
      TrainInfo ti = modelCache.get(sName);
      if (ti == null){
         System.err.printf("Error: unable to find models in cache (%s)\n", sName);
         if (modelCache.isEmpty()) System.err.printf(" No available models!\n");
         else{
            System.err.printf(" Available models:\n");
            Iterator<String> it = modelCache.keySet().iterator();
            while(it.hasNext())
               System.err.printf("  %s\n", it.next());
         }
         return false;
      }

      System.err.printf("Setting up model from cache: %s  (%d classes)\n", sName, ti.getNumClasses());

      int iModel = SupTest.getIndexFromModel(ti.model);
      int iScore = SupTest.getIndexFromScore(ti.scoreMethod);
      cpModel.setSelectedIndex(iModel);
      cbScore.setSelectedIndex(iScore);
      spinHmmStates.setValue(ti.nHmmStates);
      spinHmmSkip.setValue(ti.nHmmSkip);
      ti.updateSupTest();
      Library.copy(ti.models, models[iModel]);
      Library.copy(ti.pScores, pScore[iModel]);

      return true;
   }

   /**
    * Train models based on the current settings
    * 
    * @param bForce force training and avoid model cache
    */
   public boolean train(boolean bForce)
   {
      int wMinLen = ((Integer)spinMinLen.getValue()).intValue();
      int wMaxLen = ((Integer)spinMaxLen.getValue()).intValue();
      if (wMinLen > wMaxLen){
         System.err.println("Error: max word length is less than min length!");
         return false;
      }
      SupTest.rWordLen = new Range(wMinLen, wMaxLen);

      SupTest.per = SupTest.Per.cls;
      SupTest.nSpotsMin = ((Integer)spinMinSpots.getValue()).intValue();
      SupTest.nSpotsMax = ((Integer)spinMaxSpots.getValue()).intValue();
      try{
         SupTest.spotThresh = Double.parseDouble(tfThresh.getText());
      } catch (NumberFormatException nfe){
         SupTest.spotThresh = Double.NaN;
      }
      SupTest.bNorm = btNorm.isSelected();

      int iModel = cpModel.getSelectedIndex();
      SupTest.model = SupTest.getModelFromIndex(iModel);
      if (isHMM(iModel)){
         SupTest.hmmEval = cbHmmEval.getSelectedIndex() == 0 ? SupTest.HmmEval.forward
               : SupTest.HmmEval.viterbi;
         SupTest.hmmTrain = cbHmmTrain.getSelectedIndex() == 0 ? SupTest.HmmTrain.bw
               : SupTest.HmmTrain.viterbi;
         SupTest.nHmmStates = ((Integer)spinHmmStates.getValue()).intValue();
         SupTest.nHmmSkip = ((Integer)spinHmmSkip.getValue()).intValue();
      }
      SupTest.scoreMethod = SupTest.getScoreFromIndex(cbScore.getSelectedIndex());

      System.err.printf("Train Models (%d = %s):", iModel, cpModel.getSelectedItem());
      System.err.printf(" Word length range: %d -> %d\n", wMinLen, wMaxLen);
      if (isHMM(iModel)){
         System.err.printf("  # states=%d   # skip=%d\n", SupTest.nHmmStates, SupTest.nHmmSkip);
         System.err.printf("  Train: %s    Eval: %s\n", SupTest.hmmTrain, SupTest.hmmEval);
      }
      System.err.printf(" Scoring Method: %s\n", cbScore.getSelectedItem());
      System.err.printf(" Stop Thresh: %.4f    Min Num: %d    Max Num: %d\n", SupTest.spotThresh,
            SupTest.nSpotsMin, SupTest.nSpotsMax);
      System.err.printf(" Normalize: %s\n", btNorm.isSelected() ? "Yes" : "No");

      // has the user selected a particular model?
      if (!bForce && cbModels.getSelectedIndex() > 0)
         return setupFromCache((String)cbModels.getSelectedItem());

      // we may have more discovered motifs than known motifs
      models = new ProbSeqModel[nModelTypes][curData.size()];
      pScore = new GaussianDyn1D[nModelTypes][curData.size()];

      TimerMS timer = new TimerMS();
      System.err.printf("Training models (%d, %s)... ", curData.size(), cpModel.getItemAt(iModel));
      Iterator<String> itClass = curData.keySet().iterator();
      for(int iClass = 0; iClass < curData.size(); iClass++){
         String sClass = itClass.next();
         ArrayList<Sequence> examples = curData.get(sClass);
         ArrayList<? extends Sequence> dexamples = curQData.get(sClass);

         // find length of shortest example
         int minLen = examples.get(0).length();
         for(Sequence seq : examples)
            if (seq.length() < minLen) minLen = seq.length();
         System.err.printf(" %s) min.len: %d  #examples: %d\n", sClass, minLen, examples.size());

         switch(iModel){
         case 0: // oates
            models[iModel][iClass] = new OatesModelUSamp(examples, 0, gdata.initv, gdata.minv);
            break;
         case 1: // hmm
            HmmLRa hmm = new HmmLRa(SupTest.nHmmStates, Math.max(minLen - 1, 3), examples.get(0)
                  .getNumDims());
            // HmmLR hmm = new HmmLR(SupTest.nHmmStates, SupTest.nHmmSkip, examples.get(0).getNumDims());
            System.err.printf(" Training HMM (LR-aug) %d (%s): #states=%d  min_len=%d (%d)\n", iClass,
                  sClass, SupTest.nHmmStates, HmmUtils.getMinPathLength(hmm), minLen);
            hmm.init_segk(examples);
            if (SupTest.hmmTrain == HmmTrain.viterbi) hmm.train_viterbi(examples);
            else hmm.train_bw(examples);
            models[iModel][iClass] = hmm;
            break;
         case 2: // discrete hmm
            HmmLRD hmmd = new HmmLRD(SupTest.nHmmStates, SupTest.nHmmSkip, nSymbols);
            hmmd.init_segk(dexamples);
            if (SupTest.hmmTrain == HmmTrain.viterbi) hmmd.train_viterbi(dexamples);
            else hmmd.train_bw(dexamples);
            models[iModel][iClass] = hmmd;
            break;
         }
      }
      if (!calcPScores(iModel)) return false;
      System.err.printf("done (%dms).\n", timer.time());

      // add model info to cache and combo box
      String sName = String
            .format("%s - %s", (String)cbSeeds.getSelectedItem(), SupTest.getShortModelDesc());
      TrainInfo curTI = new TrainInfo(sName, SupTest.scoreMethod);
      curTI.setDetails(models[iModel], pScore[iModel]);
      modelCache.put(curTI.sDataSrc, curTI);
      cbModels.addItem(sName);
      cbModels.setSelectedItem(sName);

      return true;
   }

   /**
    * Perform the word spotting
    */
   protected void spot()
   {
      // setup SupTest and tran the models
      train(false);

      // setup the cont rec info
      cri = new ContRecInfo(gdata.getNumClasses());
      Iterator<String> itClass = labData.keySet().iterator();
      for(int iClass = 0; iClass < labData.size(); iClass++){
         ArrayList<Sequence> subs = labData.get(itClass.next());
         cri.nLabeledWords += subs.size();
         for(Sequence sub : subs)
            cri.nLabeledFrames += sub.length();
      }

      // do the word spotting
      int[][] confm = null;
      int iModel = cpModel.getSelectedIndex();
      if (isDiscrete(iModel)) confm = SupTest.runWordSpotDiscrete(models[iModel], pScore[iModel], qseries,
            labQData, cri);
      else confm = SupTest.runWordSpot(models[iModel], pScore[iModel], tseries, labData, cri);
      btSaveCRI.setEnabled(true);
      btSaveLabels.setEnabled(true);

      // calc best class permutation and remap spots
      System.err.print("Calculating best cluster <-> class mapping... ");
      short[] map = SupTest.calcBestMappingFromConfMatrix(confm, false);
      cri.remapSpotsOld(map);
      System.err.println("done.");

      // recalc stats with permuted classes
      System.err.print("Recomputing performance stats... ");
      confm = cri.scoreWordSpotWordsFast(labData);
      cri.scoreWordSpotFrames(tseries, labData); // TODO: should just be able to move previous numbers
      // around...
      System.err.println("done.");

      // display results
      colors = Library.generateColors(curData.size());
      displayResults(cri, true);
   }

   /**
    * Open a new window that displays various visualizations of the training data.
    */
   protected void vizTrainingData()
   {
      JFrame frame = new JFrame("Training Data Visualization");
      frame.setSize(800, 700);
      Library.centerWin(frame, null);
      frame.setVisible(true);
      train(false);

      JPanel pMain = new JPanel(new BorderLayout());
      pMain.setBackground(Color.gray);
      frame.setContentPane(pMain);

      int iModel = cpModel.getSelectedIndex();
      int nClasses = models[iModel].length;
      if (curData != null){
         GridFlexLayout gfl = new GridFlexLayout(nClasses + 1, 4, 1, 1);
         gfl.setRow(0, GridFlexLayout.Style.pref);
         gfl.setColumn(0, GridFlexLayout.Style.pref);
         JPanel pHist = new JPanel(gfl);
         pMain.add(pHist, BorderLayout.CENTER);
         pHist.add(new EmptyComponent());
         pHist.add(new JLabel("Likelihoods", JLabel.CENTER));
         pHist.add(new JLabel("Scores", JLabel.CENTER));
         pHist.add(new JLabel("Lengths", JLabel.CENTER));
         Iterator<String> it = curData.keySet().iterator();
         for(int iClass = 0; iClass < nClasses; iClass++){
            String sClass = it.next();
            pHist.add(new VerticalLabel(sClass, VerticalLabel.ROTATE_LEFT));
            ArrayList<Sequence> examples = curData.get(sClass);
            MyDoubleList scores = new MyDoubleList();
            MyDoubleList probs = new MyDoubleList();
            MyIntList lengths = new MyIntList();
            for(Sequence seq : examples){
               double x = models[iModel][iClass].eval(seq) / (btNorm.isSelected() ? seq.length() : 1.0);
               double y = SupTest.adjustScore(x, pScore[iModel][iClass]);
               probs.add(x);
               scores.add(y);
               lengths.add(seq.length());
            }

            SeqStatsView ssv = new SeqStatsView(probs.toArray());
            pHist.add(ssv);
            ssv = new SeqStatsView(scores.toArray());
            ssv.setHistColor(Color.blue);
            pHist.add(ssv);
            ssv = new SeqStatsView(lengths.toArray());
            ssv.setHistColor(Color.magenta);
            pHist.add(ssv);
         }
      }

      if (isHMM(iModel)){
         GridFlexLayout gflh = new GridFlexLayout(7, 1, 0, 10);
         gflh.setRow(0, GridFlexLayout.Style.fixed, 10);
         gflh.setRow(1, GridFlexLayout.Style.pref);
         gflh.setRow(2, GridFlexLayout.Style.fixed, 10);
         gflh.setRow(3, GridFlexLayout.Style.fixed, 200);
         gflh.setRow(4, GridFlexLayout.Style.fixed, 40);
         gflh.setRow(5, GridFlexLayout.Style.fixed, 10);
         gflh.setRow(6, GridFlexLayout.Style.fixed, 200);

         JPanel pModel = new JPanel(gflh);
         pModel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
         pMain.add(new Constrainer(pModel, 250, -1), BorderLayout.EAST);

         // spacer
         pModel.add(new EmptyComponent());

         // top: select model
         JPanel p = new JPanel();
         p.add(new JLabel("Model:"));
         cbVizModel = new JComboBox(gdata.classes);
         cbVizModel.addActionListener(this);
         p.add(cbVizModel);
         pModel.add(p);

         pModel.add(new EmptyComponent()); // spacer

         // matrix view
         double[][] m = ((AbstractHMM)models[iModel][0]).getFullTransMatrix();
         mviz = new MatrixViz(m);
         pModel.add(mviz);

         // topological view
         topoViz = new StateTopoViz(m);
         topoViz.addMouseListener(this);
         pModel.add(topoViz);

         pModel.add(new EmptyComponent()); // spacer

         // obs dist view
         pcviz = new ParallelCoordinateView();
         pcviz.setBorder(BorderFactory.createLineBorder(Color.darkGray));
         pcviz.setFixedRange(collectObsDists(iModel, -1));
         pcviz.setData(collectObsDists(iModel, 0));
         pModel.add(pcviz);
      }
   }

   /**
    * Collet a list of observation distribution means from the given model
    * 
    * @param iModel model from which to collect data
    * @param iClass index of class to collect from, -1 for all of them
    * @return list of means of obs dists
    */
   protected ArrayList<FeatureVec> collectObsDists(int iModel, int iClass)
   {
      ArrayList<FeatureVec> data = new ArrayList<FeatureVec>();
      for(int ic = 0; ic < models[iModel].length; ic++){
         if (iClass >= 0 && ic != iClass) continue;
         ProbFVModel[] odist = ((AbstractHMM)models[iModel][ic]).getStates();
         if (odist != null && odist.length > 0){
            // collect data
            if (odist[0] instanceof GaussianDiagonal){
               for(int i = 0; i < odist.length; i++)
                  data.add(((GaussianDiagonal)odist[i]).getMean());
            }
            else{
               // TODO: other obs dists?
            }
         }
      }

      return (data.isEmpty() ? null : data);
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();

      if (src == btSpot) spot();
      else if (src == cbVizModel){
         int iModel = cpModel.getSelectedIndex();
         int iClass = cbVizModel.getSelectedIndex();
         double[][] m = ((AbstractHMM)models[iModel][iClass]).getFullTransMatrix();
         mviz.setTransMatrix(m);
         topoViz.setTransMatrix(m);
         pcviz.setData(collectObsDists(iModel, iClass));
      }
      else if (src == cpModel.getComboBox()){
         btSaveHtk.setEnabled(isHMM(cpModel.getSelectedIndex()));
      }
      else if (src == cbSeeds){
         int iSelected = cbSeeds.getSelectedIndex();
         if (iSelected == 0) selectLabData();
         else selectLoadData((String)cbSeeds.getSelectedItem());
         cbSeeds.setToolTipText((String)cbSeeds.getSelectedItem());
      }
      else if (src == cbModels){
         cbModels.setToolTipText((String)cbModels.getSelectedItem());
      }
      else if (src == btLoadCRI){
         JFileChooser fc = Library.buildFileChooser("cri", "Continuous Rec Info Files (*.cri)");
         if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            if ((cri = ContRecInfo.load(f)) == null){
               System.err.println("Error: failed to load Continuous Rec Info in file:");
               System.err.println(" " + Library.getCanonical(f));
               JOptionPane.showMessageDialog(this, "Failed to load Continuous Rec Info", "Load Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
            else{
               colors = Library.generateColors(curData.size());
               displayResults(cri, true);
               btSaveCRI.setEnabled(true);
               btSaveLabels.setEnabled(true);
            }
         }
      }
      else if (src == btLoadTrain){
         JFileChooser fc = Library.buildFileChooser("seeds", "Discovered Seed Locations (*.seeds)");
         if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            ArrayList<Motif> mis = MotifIO.loadSeeds(f, new MotifInfo());
            if (mis == null) return;
            seedMotifs = new ArrayList<MotifInfo>();
            for(Motif m : mis)
               seedMotifs.add((MotifInfo)m);

            TreeMap<String, ArrayList<Sequence>> data = new TreeMap<String, ArrayList<Sequence>>();
            TreeMap<String, ArrayList<DiscreteSeq>> qdata = new TreeMap<String, ArrayList<DiscreteSeq>>();
            for(int i = 0; i < seedMotifs.size(); i++){
               MotifInfo mi = seedMotifs.get(i);
               ArrayList<Sequence> examples = WindowLocation.getExamples(mi.getOccs(), tseries);
               data.put(String.format("Motif %d", i + 1), examples);
               ArrayList<DiscreteSeq> dexamples = WindowLocation.getExamplesD(mi.getOccs(), qseries);
               qdata.put(String.format("Motif %d", i + 1), dexamples);
            }

            String sFile;
            try{
               sFile = f.getCanonicalPath();
            } catch (IOException ex){
               ex.printStackTrace();
               sFile = f.getAbsolutePath();
            }
            String sName = String.format("%s (%s)", Library.getFileName(sFile), Library.getPath(sFile));
            labLoadData.put(sName, data);
            labLoadQData.put(sName, qdata);
            cbSeeds.addItem(sName);
            cbSeeds.setSelectedItem(sName);
            selectLoadData(sName);
         }
      }
      else if (src == btShowTrain){
         int nMotifs = seedMotifs.size();
         cri = new ContRecInfo(gdata.getNumClasses());
         colors = Library.generateColors(nMotifs);
         for(int iMotif = 0; iMotif < nMotifs; iMotif++){
            MotifInfo mi = seedMotifs.get(iMotif);
            ArrayList<WindowLocation> occs = mi.getOccs();
            for(WindowLocation occ : occs){
               WordSpot spot = new WordSpot(occ.iSeries, occ.iStart, occ.nLength, 0.0, iMotif);
               cri.add(spot);
            }
         }
         displayResults(cri, false);
      }
      else if (src == btVizTrain){
         vizTrainingData();
      }
      else if (src == btTrain){
         train(true);
      }
      else if (src == btSaveCRI){
         JFileChooser fc = Library.buildFileChooser("cri", "Continuous Rec Info Files (*.cri)");
         if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            if (!f.getAbsolutePath().endsWith(".cri")) f = new File(f.getAbsolutePath() + ".cri");
            if (!cri.save(f)){
               System.err.println("Error: failed to save Continuous Rec Info in file:");
               System.err.println(" " + Library.getCanonical(f));
               JOptionPane.showMessageDialog(this, "Failed to save Continuous Rec Info", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
            else JOptionPane.showMessageDialog(this, "Saved Continuous Rec Info", "Save Successful",
                  JOptionPane.INFORMATION_MESSAGE);
         }
      }
      else if (src == btSaveLabels){
         // TODO use same code as in AbstractDiscView
         JFileChooser fc = Library.buildFileChooser("labels", "Label Files (*.labels)");
         if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            String sPath = Library.getPath(f.getAbsolutePath()); 
            String sTitle = Library.getTitle(f.getAbsolutePath());
            String sExt = Library.getExt(f.getAbsolutePath());
            if (sExt==null || sExt.length()==0) sExt = "labels";
            
            ArrayList<MarkupSet> labels = cri.getLabels(gdata.classes);            
            int nMarks = labels.size();
            MSGeneral saver = new MSGeneral();
            for(int i=0; i<nMarks; i++){
               String sFile = String.format("%s%s_%03d.%s", sPath,sTitle,i+1,sExt);
               MarkupSet marks = labels.get(i);
               if (!saver.save(marks, sFile)){
                  System.err.printf("Error: failed to save label file\n (%s)\n", sFile);
                  JOptionPane.showMessageDialog(this, "Failed to save label files", "Save Failed",
                        JOptionPane.ERROR_MESSAGE);
                  return;
               }
            }
            
            String sFileDef = String.format("%s%s.def", sPath,sTitle);
            try{
               PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(sFileDef)));               
               out.printf("!(for,i,1,1,%d)\n", nMarks);
               out.println("\n<data>");
               out.println(" name = Sequence!(%d,i)");
               out.printf(" labels = %s_!(%%03d,i).%s\n", sTitle, sExt);
               out.println(" labelLoader = MLGeneral");
               out.println("</data>");
               out.println("\n!(endfor)");
               out.close();
            }
            catch(IOException ioe){
               System.err.printf("Error: failed to save def file\n (%s)\n", sFileDef);
               JOptionPane.showMessageDialog(this, "Failed to save def file", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
               return;
            }
            
            JOptionPane.showMessageDialog(this, String.format("Saved %d label files\n (with %s)", nMarks, Library.getFileName(sFileDef)), "Save Successful",
                  JOptionPane.INFORMATION_MESSAGE);
         }
      }
      else if (src == btSaveModels){
         JFileChooser fc = Library.buildFileChooser("oates,hmm", "Serialized Probabilistic Model");
         if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION){
            train(false);
            int iModel = cpModel.getSelectedIndex();
            File f = fc.getSelectedFile();
            try{
               ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(f));
               assert (SupTest.model == SupTest.getModelFromIndex(iModel));
               out.writeObject(SupTest.model);
               if (isHMM(iModel)){
                  out.writeObject(SupTest.hmmTrain);
                  out.writeInt(SupTest.nHmmStates);
                  out.writeInt(SupTest.nHmmSkip);
               }
               out.writeObject(SupTest.scoreMethod);
               int nClasses = models[iModel].length;
               out.writeInt(nClasses);
               for(int i = 0; i < nClasses; i++){
                  out.writeObject(models[iModel][i]);
                  out.writeObject(pScore[iModel][i]);
               }
               out.close();
            } catch (Exception ex){
               System.err.println("Error: failed to save serialized model in file:");
               System.err.println(" " + Library.getCanonical(f));
               System.err.println(ex);
               JOptionPane.showMessageDialog(this, "Failed to save serialized model", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
         }
      }
      else if (src == btLoadModels){
         JFileChooser fc = Library.buildFileChooser("oates,hmm", "Serialized Probabilistic Model");
         if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
            fDataSrc = fc.getSelectedFile();
            try{
               ObjectInputStream in = new ObjectInputStream(new FileInputStream(fDataSrc));
               SupTest.model = (SupTest.Model)in.readObject();
               int iModelFile = SupTest.getIndexFromModel(SupTest.model);
               if (isHMM(iModelFile)){
                  SupTest.hmmTrain = (SupTest.HmmTrain)in.readObject();
                  SupTest.nHmmStates = in.readInt();
                  SupTest.nHmmSkip = in.readInt();
               }
               SupTest.scoreMethod = (SupTest.ScoreMethod)in.readObject();
               int nClasses = in.readInt();
               models[iModelFile] = new ProbSeqModel[nClasses];
               pScore[iModelFile] = new GaussianDyn1D[nClasses];
               for(int i = 0; i < nClasses; i++){
                  models[iModelFile][i] = (ProbSeqModel)in.readObject();
                  pScore[iModelFile][i] = (GaussianDyn1D)in.readObject();
               }
               in.close();

               cpModel.setSelectedIndex(iModelFile);

               String sFile;
               try{
                  sFile = fDataSrc.getCanonicalPath();
               } catch (IOException ex){
                  ex.printStackTrace();
                  sFile = fDataSrc.getAbsolutePath();
               }
               System.err.printf("load models: %s\n", sFile);
               String sName = String.format("%s (%s)", Library.getFileName(sFile), Library.getPath(sFile));
               labLoadData.put(sName, null);
               labLoadQData.put(sName, null);
               cbModels.addItem(sName);
               cbModels.setSelectedItem(sName);

               TrainInfo ti = new TrainInfo(sName, SupTest.model, SupTest.hmmTrain, SupTest.nHmmStates,
                     SupTest.nHmmSkip, SupTest.bNorm, SupTest.scoreMethod);
               ti.setDetails(models[iModelFile], pScore[iModelFile]);
               modelCache.put(sName, ti);
            } catch (Exception ex){
               System.err.println("Error: failed to load serialized HMM in file:");
               System.err.println(" " + Library.getCanonical(fDataSrc));
               ex.printStackTrace();
               JOptionPane.showMessageDialog(this, "Failed to load serialized HMM", "Load Failed",
                     JOptionPane.ERROR_MESSAGE);
            }
         }
      }
      else if (src == btSaveHtk){
         int nProblems = 0;
         int iModel = cpModel.getSelectedIndex();
         JFileChooser fc = Library.buildFileChooser("def", "Serialized HMM (*.def)");
         if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION){
            File f = fc.getSelectedFile();
            String sf = f.getAbsolutePath();
            if (!sf.endsWith(".def")) sf += ".def";
            String format = String.format("%s%s_%%03d.%s", Library.getPath(sf), Library.getTitle(sf),
                  Library.getExt(sf));
            System.err.println("format: " + format);
            for(int i = 0; i < models[iModel].length; i++){
               String sFile = String.format(format, i + 1);
               System.err.printf("saving model %d: %s\n", i + 1, sFile);
               if (!HtkHmm.save((HmmLR)models[iModel][i], new File(sFile))){
                  System.err.println("Error: failed to save HMMs in HTK format:");
                  System.err.println(" " + Library.getCanonical(f));
                  nProblems++;
               }
            }
         }
         if (nProblems > 0){
            JOptionPane.showMessageDialog(this, String.format("Failed to save %d/%d HMMs", nProblems,
                  models[iModel].length), "Save Failed", JOptionPane.ERROR_MESSAGE);
         }
      }
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();
   }

   protected void selectLabData()
   {
      curData = labData;
      curQData = labQData;
      collectSeedMotifs();
   }

   protected void selectLoadData(String sSrc)
   {
      curData = labLoadData.get(sSrc);
      curQData = labLoadQData.get(sSrc);
      collectSeedMotifs();
   }

   /**
    * Build a seed motif list from the currently selected data (labeled or loaded).
    */
   protected void collectSeedMotifs()
   {
      // extract labeled "motif" info from tree
      seedMotifs = new ArrayList<MotifInfo>();
      Iterator<String> it = curData.keySet().iterator();
      while(it.hasNext()){
         String sClass = it.next();
         ArrayList<Sequence> examples = curData.get(sClass);
         ArrayList<WindowLocation> wlocs = new ArrayList<WindowLocation>();
         for(Sequence seq : examples)
            wlocs.add(seq.getWindowLoc());
         MotifInfo motif = new MotifInfo(null, 0, 0, 0.0, wlocs);
         seedMotifs.add(motif);
      }
   }

   public void vizChanged(ScoreGraphComplex src, int _nDisp)
   {
      nDisp = _nDisp;
      displayResults(cri, false);
   }

   public void mouseClicked(MouseEvent e)
   {
      Object src = e.getSource();

      if (src == topoViz){
         int iState = topoViz.getState(e.getPoint());
         topoViz.setHighlight(iState);
         pcviz.setHighlight(iState);
      }
   }

   public void mousePressed(MouseEvent e)
   {}

   public void mouseReleased(MouseEvent e)
   {}

   public void mouseEntered(MouseEvent e)
   {}

   public void mouseExited(MouseEvent e)
   {}
}
