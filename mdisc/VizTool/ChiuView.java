package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.regression.SimpleRegression;

import kdm.tools.SupTest;
import kdm.data.*;
import kdm.data.transform.*;
import kdm.gui.*;
import kdm.io.*;
import kdm.io.DataSaver.*;
import kdm.metrics.*;
import kdm.mlpr.*;
import kdm.models.misc.*;
import kdm.util.*;
import kdm.util.SparseMatrix.MatrixEntry;

import static kdm.models.misc.ContRecInfo.*;

/**
 * Visual implementation of Chiu's algorithm.
 */
public class ChiuView extends AbstractDiscView implements ActionListener, ChangeListener
{
   protected LabeledSlider lsWindow, lsMotif, lsPaaSegs, lsSaxSymbols, lsMinMotifs;
   protected JComboBox cbFind;
   protected JTextField tfDist, tfInput, tfProportion;
   protected JButton btDiscover, btReal, btScanR, btScanL, btSaveLabels;
   protected JCheckBox btEstR, btGraphR;
   protected LineGraph lgraphAcc, lgraphGyr;

   // protected MetricSeq metseq = new SumSeqDist(new EuclideanFV(false), true);
   // protected MetricSeq metseq = new SumSeqDist(new AbsoluteDistFV(), true);
   public MetricSeq metseq = new DTWSimple(.1);
   public ArrayList<ArrayList<WindowLocation>> occList;
   public ArrayList<Pair<Sequence, Sequence>> seeds;
   public boolean bStopAfterOneMotif = false;
   protected double estR;

   public ChiuView(GlobalData gdata)
   {
      super(gdata);
      bVerbose = true;
      gdata.chiuView = this;

      // right panel
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(220));
      rightp.setBorder(BorderFactory.createEmptyBorder(8, 4, 8, 4));
      this.add(new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER), BorderLayout.EAST);

      lsWindow = new LabeledSlider("Window Length: %d", 3, 200);
      rightp.add(lsWindow);
      lsPaaSegs = new LabeledSlider("# PAA Segments: %d", 2, 20);
      rightp.add(lsPaaSegs);
      lsSaxSymbols = new LabeledSlider("# SAX Symbols: %d", 3, 9);
      rightp.add(lsSaxSymbols);
      lsMinMotifs = new LabeledSlider("Min # Motifs: %d", 2, 50);
      rightp.add(lsMinMotifs);
      lsMotif = new LabeledSlider("Show Motif: %d", 1, 1);
      lsMotif.setEnabled(false);
      lsMotif.addChangeListener(this);
      rightp.add(lsMotif);

      lsWindow.setValue(14);
      lsPaaSegs.setValue(4);
      lsSaxSymbols.setValue(5);
      lsMinMotifs.setValue(12);

      tfDist = new JTextField(String.format("%.3f", 9.5), 8);
      btReal = new JButton("Real");
      btReal.addActionListener(this);
      if (gdata.labData == null) btReal.setEnabled(false);
      JPanel p = new JPanel();
      p.add(new JLabel("R: "));
      p.add(tfDist);
      p.add(btReal);
      rightp.add(p);

      btDiscover = new JButton("Discover!");
      btDiscover.addActionListener(this);
      rightp.add(btDiscover);
      
      rightp.add(Box.createVerticalStrut(4));
      
      btSaveLabels = new JButton("Save Labels...");
      btSaveLabels.setEnabled(false);
      btSaveLabels.addActionListener(this);
      rightp.add(btSaveLabels);      
      
      lgraphAcc = new LineGraph();
      lgraphAcc.setShowAll(true);
      lgraphAcc.setRenderMouseInfo(true, false, false);
      lgraphAcc.setDimColor(0, Color.red);
      lgraphAcc.setDimColor(1, Color.green);
      lgraphAcc.setDimColor(2, Color.blue);
      rightp.add(Box.createVerticalStrut(8));
      rightp.add(new Constrainer(lgraphAcc, 200, 60));

      lgraphGyr = new LineGraph();
      lgraphGyr.setShowAll(true);
      lgraphGyr.setRenderMouseInfo(true, false, false);
      lgraphGyr.setDimColor(3, Color.red);
      lgraphGyr.setDimColor(4, Color.green);
      lgraphGyr.setDimColor(5, Color.blue);
      rightp.add(Box.createVerticalStrut(8));
      rightp.add(new Constrainer(lgraphGyr, 200, 60));

      // left panel
      VerticalScrollPanel leftp = new VerticalScrollPanel(new VerticalLayout(220));
      leftp.setBorder(BorderFactory.createEmptyBorder(8, 4, 8, 4));
      this.add(new JScrollPane(leftp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER), BorderLayout.WEST);

      tfInput = new JTextField(12);
      leftp.add(tfInput);
      leftp.add(Box.createVerticalStrut(6));
      p = new JPanel();
      tfProportion = new JTextField(12);
      tfProportion.setText("10.0");
      p.add(new JLabel("Prop: "));
      p.add(tfProportion);
      leftp.add(p);
      leftp.add(Box.createVerticalStrut(6));
      btScanR = new JButton("Scan: R");
      btScanR.addActionListener(this);
      leftp.add(btScanR);
      leftp.add(Box.createVerticalStrut(6));
      btScanL = new JButton("Scan: Length");
      btScanL.addActionListener(this);
      leftp.add(btScanL);
      leftp.add(Box.createVerticalStrut(6));
      cbFind = new JComboBox(new String[] { "Temporal Greedy", "Look-ahead Greedy", "Dist. Greedy" });
      leftp.add(cbFind);
      btEstR = new JCheckBox("Estimate R", true);
      leftp.add(btEstR);
      btGraphR = new JCheckBox("Show Graph R", false);
      leftp.add(btGraphR);
   }

   /**
    * Discover according to current GUI settings
    */
   protected void run()
   {
      int wlen = lsWindow.getValue();
      int nPAA = lsPaaSegs.getValue();
      int nSAX = lsSaxSymbols.getValue();      
      boolean bEstR = btEstR.isSelected();
      boolean bShowGraph = btGraphR.isSelected();
      double R = bEstR ? 0 : Double.parseDouble(tfDist.getText());
      int kFindOccs = cbFind.getSelectedIndex();
      double fProp = Double.parseDouble(tfProportion.getText()) / 100.0;
      int nMinOccs = lsMinMotifs.getValue();

      run(wlen, nPAA, nSAX, bEstR, bShowGraph, R, fProp, kFindOccs, nMinOccs);
      
      scoreMotifs(occList, gdata.labData);
      lsMotif.setText(String.format("Motif: %%d / %d", occList.size()));
      lsMotif.setValues(1, 1, occList.size());
      btSaveLabels.setEnabled(true);
   }
   
   /**
    * Discover according to given settings
    * @param wlen window length
    * @param nPAA number of PAA segments
    * @param nSAX number of SAX symbols
    * @param bEstR true for neighborhood radius estimation
    * @param bShowGraph true to visually display graph during R estimation (ignored if bEstR is false)
    * @param R fixed size neighborhood radius (ignored if bEstR is true) 
    * @param fProp proportion of distances to use for R estimation (ignored if bEstR is false)
    * @param kFindOccs method for finding occurrences (0, 1, or 2)
    * @param nMinOccs minimum number of occurrences to be considered a valid motif
    */
   public void run(int wlen, int nPAA, int nSAX, boolean bEstR, boolean bShowGraph, double R, double fProp, int kFindOccs, int nMinOccs)
   {
      // double MeanLinearThresh = 0.035; // TODO include linear (only flat?) section removal

      TimerMS timer = new TimerMS();
      TimerMS timer2 = new TimerMS();
      
      int nSeqs = tseries.size();
      int nDims = gdata.getNumDims();
      int nTotal = gdata.nTotalSeriesLength;      
      
      if (bVerbose)
         System.err.printf(
               "wlen=%d  #paa=%d  #sax=%d  #seqs=%d  #dims=%d  #total=%d  min#occs=%d  est=%b  R=%.3f\n",
               wlen, nPAA, nSAX, nSeqs, nDims, nTotal, nMinOccs, bEstR, R);

      // calculate data density
      /*
       * System.err.print("Calculating data density... "); timer.reset(); Sequence all = new Sequence();
       * for(Sequence seq : tseries) all.append(seq, false); density = getDensity(all);
       * System.err.printf("done (%dms).\n", timer.time());
       */

      // build a full span list
      SpanList[] spans = new SpanList[nSeqs];
      for(int i = 0; i < nSeqs; i++)
         spans[i] = new SpanList(0, gdata.seqLens[i] - wlen, true);

      // first, we should detect and ignore linear segments
      // System.err.printf("Removing linear segments... ");
      // timer.reset();
      // WindowLocation wloc1 = new WindowLocation(-1, -1, wlen);
      // WindowLocation wloc2 = new WindowLocation(-1, -1, wlen);
      // SimpleRegression reg = new SimpleRegression();
      // timer.reset();
      // double[] y = new double[wlen];
      // for(int iSeq = 0; iSeq < nSeqs; iSeq++)q
      // {
      // wloc1.iSeries = iSeq;
      // Sequence seq1 = tseries.get(iSeq);
      // for(int ix = 0; ix + wlen <= gdata.seqLens[iSeq]; ix++)
      // {
      // wloc1.iStart = ix;
      // boolean bLinear = true;
      // seq: for(int d = 0; d < nDims; d++)
      // {
      // // extract data in this dimension
      // for(int i = 0; i < wlen; i++)
      // y[i] = tseries.get(iSeq).get(i + ix, d);
      //
      // // normalize the data
      // double ymean = StatUtils.mean(y);
      // double sdev = Math.sqrt(StatUtils.variance(y));
      // for(int i = 0; i < wlen; i++)
      // y[i] = (y[i] - ymean) / sdev;
      //
      // // compute the residual error (mse)
      // reg.clear();
      // for(int i = 0; i < wlen; i++)
      // reg.addData(i, y[i]);
      // double residual = reg.getMeanSquareError();
      // double meanRes = residual / wlen;
      // // System.err.printf("%d.%d.%d: %.2f\n", iSeq, ix, d, meanRes);
      // if (meanRes > MeanLinearThresh)
      // {
      // bLinear = false;
      // break seq;
      // }
      // }
      // if (bLinear) spans[iSeq].sub(wloc1.getFirstIndex() - wlen + 1, wloc1.getLastIndex());
      // }
      // }

      // report # spots removed, time
      int nRemoved = 0;
      // for(int i = 0; i < spans.length; i++)
      // nRemoved += spans[i].getNumRemoved();
      // System.err.printf("done (removed %d spots (%.1f%%) in %dms).\n", nRemoved, 100.0 * nRemoved / nTotal,
      // timer.time());
      int nAfterLinear = nTotal - nRemoved;

      // RLE over SAX strings
      // double[] paa = new double[nPAA];
      // char[] symbols = new char[nPAA];
      // timer.reset();
      // int nDups = 0;
      // for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      // {
      // String sPrev = null;
      // Sequence seq = tseries.get(iSeq);
      // double[][] raw = seq.toSeqArray();
      // SpanList spanRemove = new SpanList(spans[iSeq].getRangeMin(), spans[iSeq].getRangeMax(), false);
      // SpanIterator it = spans[iSeq].iterator();
      // while(it.hasMore())
      // {
      // int t = it.next();
      // StringBuffer sb = new StringBuffer();
      // for(int d = 0; d < nDims; d++)
      // {
      // assert (t + wlen <= raw[d].length) : String.format("t=%d wlen=%d len=%d last=%d", t, wlen,
      // raw[d].length, spans[iSeq].get(spans[iSeq].size() - 1));
      // SAX.genPAA(raw[d], t, wlen, paa);
      // for(int i = 0; i < nPAA; i++)
      // symbols[i] = SAX.raw2sym(paa[i], nSAX);
      // sb.append(symbols);
      // }
      // String s = sb.toString();
      // if (sPrev != null && s.equals(sPrev)) spanRemove.add(t);
      // else sPrev = s;
      // }
      // //spans[iSeq].sub(spanRemove); // TODO: for now, don't really remove them
      // nDups += spanRemove.size();
      // }
      // System.err.printf("SAX string RLE: %d dups (%dms)\n", nDups, timer.time());

      // we have to build the collision matrix
      if (bVerbose){
         int nRounds = nPAA * (nPAA - 1) / 2;
         System.err.printf("Building Collision Matrix (%d rounds)... ", nRounds);
      }
      timer.reset();
      SparseMatrix<MutableInteger> mat = buildCollisionMatrix(nPAA, nSAX, wlen, spans);
      if (bVerbose) System.err.printf("done (%d, %dms).\n", mat.size(), timer.time());

      long n = (long)nTotal * (nTotal - 1) / 2;
      int n2 = nAfterLinear * (nAfterLinear - 1) / 2;
      int hit = mat.size();
      if (bVerbose)
         System.err.printf(" hit: %d / %d = %.2f%%  (%d / %d = %.2f%%)\n", hit, n, 100.0 * hit / n, hit, n2,
               100.0 * hit / n2);

      // now that we have the collision matrix, we can start enumerating motifs
      if (bVerbose)
         System.err.printf("Enumerating motifs (R=%s)... ", bEstR ? "est" : String.format("%.2f", R));
      timer.reset();
      SpanList spanColl = new SpanList(0, nTotal - 1, true);
      occList = new ArrayList<ArrayList<WindowLocation>>();
      seeds = new ArrayList<Pair<Sequence, Sequence>>();
      lsMotif.setEnabled(true);
      int iMotif = 0;
      int[] aSeed = new int[2];
      while(true){
         SparseMatrix<MutableInteger>.MatrixEntry me = findBestCollision(mat, spanColl);
         if (me == null) break;
         aSeed[0] = me.row();
         aSeed[1] = me.column();
         if (aSeed[0] < 0 || aSeed[1] < 0) break;

         assert (aSeed[0] < aSeed[1]);
         int nHits = me.get().getValue();
         mat.remove(aSeed[0], aSeed[1]); // we never want to reprocess
         if (Math.abs(aSeed[0] - aSeed[1]) < wlen) continue; // no trivial matches

         // if (bVerbose) System.err.printf("aSeed: %d, %d (%d)\n", aSeed[0], aSeed[1], nHits);
         if (bEstR) estR = R = estimateR(aSeed, fProp, wlen, bShowGraph, iMotif);

         // make sure dist between these is < R
         int[] ofs = Library.getArrayOffset(aSeed[0], gdata.seqLens, null);
         WindowLocation wloc1 = new WindowLocation(ofs[0], ofs[1], wlen);
         Sequence seed1 = wloc1.getSeq(tseries);
         Library.getArrayOffset(aSeed[1], gdata.seqLens, ofs);
         WindowLocation wloc2 = new WindowLocation(ofs[0], ofs[1], wlen);
         Sequence seed2 = wloc2.getSeq(tseries);
         double dist = metseq.dist(seed1, seed2) / wlen;
         if (bVerbose)
            System.err.printf("Dist between seeds: %.3f  (R=%.3f)  %s  %s\n", dist, R, wloc1, wloc2);
         if (dist > R){
            // System.err .printf("Warning: seed pair too far apart (%d<->%d: %.1f)\n", aSeed[0], aSeed[1],
            // dist);
            spanColl.sub(aSeed[0]);
            spanColl.sub(aSeed[1]);
            continue;
         }

         ArrayList<WindowLocation> occs = null;
         switch(kFindOccs){
         case 0:
            occs = findOccs(wloc1, wloc2, R, spans, metseq, true);
            break;
         case 1:
            occs = findOccs2(wloc1, wloc2, R, spans);
            break;
         case 2:
            occs = findOccs3(wloc1, wloc2, R, spans);
            break;
         }
         int nOccs = occs.size();

         // remove these occs from the collision matrix
         for(int i = 0; i < nOccs; i++){
            WindowLocation wloc = occs.get(i);
            // System.err.printf("Occ %d: %s\n", i+1, wloc);
            spans[wloc.iSeries].sub(wloc.iStart - wlen + 1, wloc.getLastIndex());
            int indexA = wloc.getLinearIndex(gdata.seqLens);
            spanColl.sub(indexA - wlen + 1, indexA + wlen - 1);
         }

         if (occs != null && occs.size() >= nMinOccs){
            occList.add(occs);
            seeds.add(new Pair<Sequence, Sequence>(seed1, seed2));
            double mdl = 0;// calcDL(occs); // TODO
            if (bVerbose){
               System.err.printf("Motif %d: %d occs, R=%.3f, DL=%.2f\n", iMotif + 1, occs.size(), R, mdl);
               //for(WindowLocation wloc : occs) System.err.printf(" %s\n", wloc);
            }
            iMotif++;
            if (bStopAfterOneMotif) break;
         }
      }
      if (bVerbose) System.err.printf("done (%dms).\n", timer.time());
      if (bVerbose)
         System.err.printf("Total time for discovery (%d motifs): %dms\n", occList.size(), timer2.time());      
   }

   /** Find exact matches between multi-dim SAX strings */
   HashMap<String, MyIntList> findExactMatches(int nPAA, int nSAX, int wlen, SpanList[] spans)
   {
      HashMap<String, MyIntList> map = new HashMap<String, MyIntList>();
      int nSeqs = tseries.size();
      int nDims = gdata.getNumDims();
      double[] paa = new double[nPAA];
      char symbols[] = new char[nPAA];

      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = tseries.get(iSeq);
         double[][] raw = seq.toSeqArray();
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            int index = Library.getIndexFromArrayOffset(iSeq, t, gdata.seqLens);

            // build the SAX string
            StringBuffer sb = new StringBuffer();
            for(int d = 0; d < nDims; d++){
               SAX.genPAA(raw[d], t, wlen, paa);
               for(int i = 0; i < nPAA; i++)
                  symbols[i] = SAX.raw2sym(paa[i], nSAX);
               sb.append(symbols);
            }
            String s = sb.toString();
            MyIntList list = map.get(s);
            if (list == null){
               list = new MyIntList();
               map.put(s, list);
            }
            list.add(index);
         }
      }
      return map;
   }

   /** build the collision matrix */
   protected SparseMatrix<MutableInteger> buildCollisionMatrix(int nPAA, int nSAX, int wlen, SpanList[] spans)
   {
      int nSeqs = tseries.size();
      int nDims = gdata.getNumDims();
      double[] paa = new double[nPAA];
      char symbols[][] = new char[nDims][nPAA];

      // TODO if #round is too large (> 100?, >50?), then do a random sampling
      int nRounds = nPAA * (nPAA - 1) / 2;
      ArrayList<HashMap<String, MyIntList>> buckets = new ArrayList<HashMap<String, MyIntList>>();
      for(int iRound = 0; iRound < nRounds; iRound++)
         buckets.add(new HashMap<String, MyIntList>());

      // now we can build the matrix
      SparseMatrix<MutableInteger> mat = new SparseMatrix<MutableInteger>();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = tseries.get(iSeq);
         double[][] raw = seq.toSeqArray();
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            int index = Library.getIndexFromArrayOffset(iSeq, t, gdata.seqLens);

            // build the SAX string
            for(int d = 0; d < nDims; d++){
               SAX.genPAA(raw[d], t, wlen, paa);
               for(int i = 0; i < nPAA; i++)
                  symbols[d][i] = SAX.raw2sym(paa[i], nSAX);
            }

            int iRound = 0;
            for(int iCol = 0; iCol < nPAA; iCol++)
               for(int jCol = iCol + 1; jCol < nPAA; jCol++){
                  StringBuffer sb = new StringBuffer();
                  for(int d = 0; d < nDims; d++){
                     sb.append(symbols[d][iCol]);
                     sb.append(symbols[d][jCol]);
                  }
                  String s = sb.toString();
                  MyIntList list = buckets.get(iRound).get(s);
                  if (list == null){
                     list = new MyIntList();
                     buckets.get(iRound).put(s, list);
                  }
                  list.add(index);
                  iRound++;
               }
         }
      }

      // increment all collisions
      for(int iRound = 0; iRound < nRounds; iRound++){
         Iterator<String> it = buckets.get(iRound).keySet().iterator();
         while(it.hasNext()){
            String s = it.next();
            MyIntList list = buckets.get(iRound).get(s);
            int n = list.size();
            for(int i = 0; i < n; i++)
               for(int j = i + 1; j < n; j++){
                  int a = list.get(i);
                  int b = list.get(j);
                  assert (a < b); // true due to construction method
                  MutableInteger mi = mat.get(a, b);
                  if (mi == null){
                     mi = new MutableInteger(1);
                     mat.set(a, b, mi);
                  }
                  else mi.inc();
               }
         }
      }

      return mat;
   }

   /**
    * Find the matrix entry that is contained in the span list with the largest value
    * 
    * @param mat collision matrix
    * @param span allowed locations
    * @return index of largest allowed row or -1 if none found
    */
   protected SparseMatrix<MutableInteger>.MatrixEntry findBestCollision(SparseMatrix<MutableInteger> mat,
         SpanList span)
   {
      SparseMatrix<MutableInteger>.MatrixEntry meBest = null;
      SpanIterator it = span.iterator();
      while(it.hasMore()){
         int i = it.next();
         HashMap<Integer, MutableInteger> row = mat.getRow(i);
         if (row == null) continue;
         Iterator<Integer> jt = row.keySet().iterator();
         while(jt.hasNext()){
            int j = jt.next();
            SparseMatrix<MutableInteger>.MatrixEntry me = mat.getEntry(i, j);
            if (meBest == null || me.get().getValue() > meBest.get().getValue()) meBest = me;
         }
      }
      return meBest;
   }

   /**
    * Find the matrix row that is contained in the span list with the largest sum
    * 
    * @param coll collision matrix (sum of rows)
    * @param span allowed locations
    * @return index of largest allowed row or -1 if none found
    */
   protected int findBestCollision(int[] coll, SpanList span)
   {
      int iBest = -1;
      SpanIterator it = span.iterator();
      while(it.hasNext()){
         int i = it.next();
         if (coll[i] < 2) continue;
         if (iBest < 0 || coll[i] > coll[iBest]) iBest = i;
      }
      return iBest;
   }

   /**
    * Currently does NOT work! Temporally greedy method with a look-ahead. Should give better results than
    * findOccs() with similar run-time characteristics.
    * 
    * @param wlSeed1 seed subsequence
    * @param wlSeed2 seed subsequence
    * @param R max distance
    * @param spans allowed spots for new occurrences
    * @return list of occurrences
    */
   protected ArrayList<WindowLocation> findOccs2(WindowLocation wlSeed1, WindowLocation wlSeed2, double R,
         SpanList[] spans)
   {
      // TODO doesn't work
      assert false : "doesn't work!";

      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      Sequence seed = wlSeed1.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         int ixd = 0, cut = 0;
         double d = Double.NaN;
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            Sequence win = tseries.get(iSeq).subseq(t, t + wlen);
            double r = metseq.dist(seed, win) / wlen;

            if (Double.isNaN(d)){
               d = r;
               ixd = t;
               cut = t + wlen;
            }
            else{
               if (t >= cut){
                  if (d <= R) occs.add(new WindowLocation(iSeq, ixd, wlen));
                  d = r;
                  ixd = t;
                  cut = t + wlen;
               }
               else if (r > d){
                  if (d <= R) occs.add(new WindowLocation(iSeq, ixd, wlen));
                  d = Double.NaN;
                  it.jump(ixd + wlen);
               }
            }
         }
         if (!Double.isNaN(d)){
            assert (ixd + wlen <= gdata.seqLens[iSeq]);
            if (d <= R) occs.add(new WindowLocation(iSeq, ixd, wlen));
         }
      }

      return occs;
   }

   /**
    * Greedy method in terms of distance; uses a priority queue to get the closest value first
    * 
    * @param wlSeed1 seed subsequence
    * @param wlSeed2 seed subsequence
    * @param R max distance
    * @param spans allowed spots for new occurrences
    * @return list of occurrences
    */
   protected ArrayList<WindowLocation> findOccs3(WindowLocation wlSeed1, WindowLocation wlSeed2, double R,
         SpanList[] spans)
   {
      PriorityQueue<WindowLocation> q = new PriorityQueue<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      assert (wlen == wlSeed2.length());
      Sequence seed1 = wlSeed1.getSeq(tseries);
      Sequence seed2 = wlSeed1.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = tseries.get(iSeq);
         int T = gdata.seqLens[iSeq];
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            Sequence win = seq.subseq(t, t + wlen);
            double dist1 = metseq.dist(seed1, win) / wlen;
            double dist2 = Library.INF;
            if (dist1 < R || (dist2 = metseq.dist(seed2, win) / wlen) < R){
               q.add(new WindowLocation(iSeq, t, wlen));
               it.jump(t + wlen);
            }
         }
      }

      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      WindowLocation win;
      scored: while((win = q.poll()) != null){
         for(WindowLocation occ : occs)
            if (occ.overlaps(win)) continue scored;
         occs.add(win);
      }
      return occs;

   }

   protected double calcDL(Pair<Sequence, Sequence> seeds, ArrayList<WindowLocation> occs)
   {
      double err = 0;
      for(WindowLocation wloc : occs){
         Sequence seq = wloc.getSeq(tseries);
         double d1 = metseq.dist(seeds.first, seq);
         double d2 = metseq.dist(seeds.second, seq);
         err += Math.max(d1, d2);
      }
      err /= (occs.size() * seeds.first.length());
      return err;
   }

   /**
    * Compute the description length of the given occurrences using KDE
    * 
    * @param occs motif occurrences
    * @return measure of information in the occurrences
    */
   protected double calcDL(ArrayList<WindowLocation> occs)
   {
      if (occs == null || occs.isEmpty()) return 0;

      TimerMS timer = new TimerMS();

      // build average motif
      int nOccs = occs.size();
      Sequence motif = new Sequence(occs.get(0).getSeq(tseries));
      int wlen = motif.length();
      int nFrames = nOccs * wlen;
      int nDims = motif.getNumDims();
      for(int i = 1; i < nOccs; i++)
         motif.addFrames(occs.get(1).getSeq(tseries));
      motif._div(nOccs);
      double dlMean = 0; // TODO: need to calc this (use density array...)

      // calc mean err (dist from mean)
      double meanErr = 0;
      for(WindowLocation occ : occs){
         Sequence seq = occ.getSeq(tseries);
         for(int i = 0; i < wlen; i++){
            FeatureVec fv = seq.get(i).sub(motif.get(i));
            meanErr += fv._mul(fv).sum();
         }
      }
      //System.out.printf("DL) wlen: %d  nFrames: %d  nDims: %d  nOccs: %d  total err: %.2f  (%.2f)\n", wlen,
      //      nFrames, nDims, nOccs, meanErr, meanErr / nOccs); // TODO
      meanErr /= (nFrames * nDims);
      return meanErr;

      // calc DL of errors
      /*
       * double[] denErr = getDensity(seqErrors); double dlError = 0; for(int i=0; i<denErr.length; i++)
       * dlError -= denErr[i] * Math.log(denErr[i]); dlError /= Math.log(2.0); // we want bits
       * 
       * System.err.printf("#occs: %d #frames: %d dlMean=%.3f dlError=%.3f\n", nOccs, nFrames, dlMean,
       * dlError);
       * 
       * return (dlMean + dlError) / nFrames;
       */
   }

   /**
    * Calculate the density at each of the given points
    * 
    * @return density of each point in the sequence
    */
   protected double[] getDensity(Sequence data)
   {
      try{
         TimerMS timer = new TimerMS();
         File file = File.createTempFile("kdm_", null);
         DSRaw saver = new DSRaw();
         saver.save(data, file.getAbsolutePath());

         String cmd = String.format("kde -scaling none -noloo -data %s", file.getAbsoluteFile());
         // System.err.printf("Executing cmd: \"%s\"\n", cmd);
         Runtime rt = Runtime.getRuntime();
         Process proc = rt.exec(cmd);
         proc.waitFor();

         String sKde = Library.getPath(file.getAbsolutePath()) + Library.getTitle(file.getName())
               + ".expt1.dens";
         // System.err.printf(" output file: %s\n exist: %b\n", sKde, new File(sKde).exists());

         // now read the density back
         double[] a = new double[data.length()];
         BufferedReader in = new BufferedReader(new FileReader(sKde));
         for(int i = 0; i < a.length; i++)
            a[i] = Double.parseDouble(in.readLine());
         in.close();
         // System.err.printf("Time to get KDE: %dms\n", timer.time());
         return a;

      } catch (IOException e){
         System.err.printf("Error: failed to create temp file for kde (%s)\n", e);
         assert false : "failed to create temp file for kde";
         return null;
      } catch (InterruptedException e){
         System.err.printf("Error: kde execution was interrupted (%s)\n", e);
         assert false : "kde exec failed";
         return null;
      }
   }

   /** Dump stats on the real distance between labeled motif members */
   public void dumpRealR()
   {
      Iterator<String> it = gdata.labData.keySet().iterator();
      while(it.hasNext()){
         String sClass = it.next();
         ArrayList<Sequence> list = gdata.labData.get(sClass);
         int N = list.size();
         double[] R = new double[N * (N - 1) / 2];
         int ixR = 0;
         for(int i = 0; i < N; i++){
            for(int j = i + 1; j < N; j++)
               R[ixR++] = metseq.dist(list.get(i), list.get(j))
                     / Math.max(list.get(i).length(), list.get(j).length());
         }

         System.err.printf("%-20s: [%.2f  %.2f  %.2f]\n", sClass, StatUtils.min(R), StatUtils.mean(R),
               StatUtils.max(R));
      }
   }

   /** calculate fake tanaka MDL under non-compressed assumption */
   protected double calcTanakaMDL(ArrayList<WindowLocation> occs)
   {
      int nSeqs = tseries.size();
      SpanList[] spans = new SpanList[nSeqs];
      for(int i = 0; i < nSeqs; i++)
         spans[i] = new SpanList(0, gdata.seqLens[i] - 1, true);

      // remove the spots from the sequences
      for(WindowLocation wloc : occs)
         spans[wloc.iSeries].sub(wloc.getRange());

      // calc # segments
      int m = occs.size();
      for(int i = 0; i < nSeqs; i++)
         m += spans[i].getNumSpans();

      // DL_3 is easy
      double dl3 = m * Math.log(gdata.nTotalSeriesLength);

      // DL_1 and DL_2 are the same if we don't "compress" the string (which we can't for Chiu)
      double dl1 = 0;
      for(WindowLocation wloc : occs)
         dl1 += Math.log(wloc.length());
      for(int i = 0; i < nSeqs; i++){
         int nr = spans[i].getNumSpans();
         for(int j = 0; j < nr; j++)
            dl1 += Math.log(spans[i].getRange(j).length());
      }

      return 2 * dl1 + dl3;
   }

   /** Estimate R for the given seeds */
   protected double estimateR(int[] aSeed, double fProp, int wlen, boolean bShowGraph, int iMotif)
   {
      TimerNS timer = new TimerNS();      
      int nSeqs = tseries.size();
      Sequence seed1, seed2 = null;
      WindowLocation wlSeed1, wlSeed2 = null;

      int[] ofs = Library.getArrayOffset(aSeed[0], gdata.seqLens, null);
      wlSeed1 = new WindowLocation(ofs[0], ofs[1], wlen);
      seed1 = tseries.get(ofs[0]).subseq(ofs[1], ofs[1] + wlen);
      if (aSeed.length > 1){
         ofs = Library.getArrayOffset(aSeed[1], gdata.seqLens, null);
         wlSeed2 = new WindowLocation(ofs[0], ofs[1], wlen);
         seed2 = tseries.get(ofs[0]).subseq(ofs[1], ofs[1] + wlen);
      }
      int dj = Math.max(2, (int)Math.ceil(wlen / 10.0)); // TODO
      //System.err.printf(" wlen=%d  dj=%d\n", wlen, dj);
      MyDoubleList list = new MyDoubleList();
      WindowLocation wloc = new WindowLocation(-1, -1, wlen);
      for(int i = 0; i < nSeqs; i++){
         int ixd = 0;
         double d = Double.NaN;
         wloc.iSeries = i;
         for(int j = 0; j + wlen <= gdata.seqLens[i]; j += dj){
            wloc.iStart = j;
            Sequence win = tseries.get(i).subseq(j, j + wlen);
            if (wlSeed1.overlaps(wloc) || (wlSeed2 != null && wlSeed2.overlaps(wloc))) continue;
            double r = metseq.dist(seed1, win) / wlen;
            if (seed2 != null){
               double r2 = metseq.dist(seed2, win)/ wlen;
               if (r2 < r) r = r2;
            }
            list.add(r);
         }
      }

      if (list.isEmpty()) return Double.NaN;

      double[] rAll = list.toArray();
      int n = (int)Math.round(fProp * rAll.length);
      int iMax = Library.select(rAll, 0, rAll.length - 1, n);
      double rmax = rAll[iMax];
      double[] r = new double[n];
      int ixr = 0;
      for(int i = 0; i < rAll.length; i++)
         if (rAll[i] <= rmax) r[ixr++] = rAll[i];
      assert (ixr == r.length) : String.format("ixr=%d n=%d", ixr, n);
      
      double knee = Library.findKnee(r, true);
      double rad = r[(int)Math.round(knee)];
      if (bShowGraph){
         Sequence seqr = new Sequence("R", r);  
         seqr.setStartMS(0);
         System.err.printf("seqr: start=%d  end=%d  knee=%.2f  rad=%.2f\n", seqr.getStartMS(), seqr.getEndMS(), knee, rad);
         LineGraph lg = new LineGraph(seqr);         
         lg.setVertLine(Math.round(knee*seqr.getPeriod()*1000.0f), Color.yellow, 1.0f);
         JFrame frame = new JFrame(String.format("R - motif: %d", iMotif+1));
         frame.getContentPane().add(new GraphComplexLite(lg));
         frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
         frame.setSize(400,300);
         frame.setLocation(120, 80);
         frame.setVisible(true);
      }

      //System.err.printf(" est r: %.3f (%.1fms)\n", rad, timer.timeMS());
      return rad;
   }

   /** try different R values */
   protected void scanR()
   {
      String s = tfInput.getText();
      StringTokenizer st = new StringTokenizer(s, " :");
      if (st.countTokens() != 3){
         System.err.printf("Error: input string has incorrect format for scanR (<start>:<step>:<stop>)\n");
         return;
      }
      double start = Double.parseDouble(st.nextToken());
      double step = Double.parseDouble(st.nextToken());
      double stop = Double.parseDouble(st.nextToken());

      System.err.printf("ScanR (%.2f:%.2f:%.2f)  event.(acc,prec,recall)  frame.(acc,prec,recall)\n", start, step, stop);
      
      String sROrig = tfDist.getText();
      boolean bEstOrig = btEstR.isSelected();
      btEstR.setSelected(false);

      // bStopAfterOneMotif = true;
      boolean bVerboseOrig = bVerbose;
      bVerbose = false;

      for(double r = start; (r - stop) < 0.000001; r += step){
         // System.out.printf("Discovering first motif for R=%.3f...\n", r);
         tfDist.setText(String.format("%.3f", r));
         System.err.printf("%.2f  ", r);
         run();
         /*
          * if (cri == null || cri.getNumSpots() == 0) { System.out.println(" No motifs."); continue; }
          *  // which class is best? int iBestClass = 0; for(int i = 1; i < gdata.getNumClasses(); i++) { if
          * (cri.getClassStats(i)[CORRECT] > cri.getClassStats(iBestClass)[CORRECT]) iBestClass = i; }
          * 
          * assert (seeds.size() == 1) : String.format("# seeds: %d", seeds.size()); assert (occList.size() ==
          * 1) : String.format("# occ lists: %d\n", occList.size());
          *  // calc MDL value ArrayList<WindowLocation> occs = occList.get(0); double mdl = calcDL(occs);
          * 
          * int[] a = cri.getClassStats(iBestClass);
          *  // System.out.printf( " Discovered class: %d (%s) -- CIDS: %d, %d, %d, %d apr: %.1f, %.1f,
          * %.1f\n", // iBestClass + 1, gdata.classes[iBestClass], a[CORRECT], a[INSERTION], a[DELETION],
          * a[SUBSTITUTION], // 100*cri.getClassAcc(iBestClass), 100*cri.getClassPrec(iBestClass), 100*cri //
          * .getClassRecall(iBestClass));
          * 
          * //System.out.printf("%.6f %.2f %.2f %.2f %.4f\n", r, 100 * cri.getClassAcc(iBestClass), 100 * cri //
          * .getClassPrec(iBestClass), 100 * cri.getClassRecall(iBestClass), mdl);
          */}

      btEstR.setSelected(bEstOrig);
      tfDist.setText(sROrig);
      bStopAfterOneMotif = false;
      bVerbose = bVerboseOrig;
   }

   /** try different motif lengths */
   protected void scanL()
   {
      String s = tfInput.getText();
      StringTokenizer st = new StringTokenizer(s, " :");
      int nToks = st.countTokens();
      int start, step, stop;
      if (nToks == 2){
         start = Integer.parseInt(st.nextToken());
         step = 1;
         stop = Integer.parseInt(st.nextToken());
      }
      else if (nToks == 3){
         start = Integer.parseInt(st.nextToken());
         step = Integer.parseInt(st.nextToken());
         stop = Integer.parseInt(st.nextToken());
      }
      else{
         System.err
               .printf("Error: input string has incorrect format for scanR (<start>:<step>:<stop> or <start>:<stop>)\n");
         return;
      }

      int lenOrig = lsWindow.getValue();
      boolean bEstOrig = btEstR.isSelected();
      //btEstR.setSelected(true);

      bStopAfterOneMotif = true;
      boolean bVerboseOrig = bVerbose;
      bVerbose = false;

      System.err.printf("ScanL (%d:%d:%d)  event.(acc,prec,recall)  frame.(acc,prec,recall)\n", start, step, stop);
      
      for(int wlen = start; wlen <= stop; wlen += step){
         lsWindow.setValue(wlen);
         System.err.printf("%d:  ", wlen);
         run();

         if (cri == null || cri.getNumSpots() == 0){
            System.out.println(" No motifs.");
            continue;
         }

         // which class is best?
         int iBestClass = 0;
         for(int i = 1; i < gdata.getNumClasses(); i++){
            if (cri.getClassStats(i)[CORRECT] > cri.getClassStats(iBestClass)[CORRECT]) iBestClass = i;
         }

         assert (seeds.size() == 1) : String.format("# seeds: %d", seeds.size());
         assert (occList.size() == 1) : String.format("# occ lists: %d\n", occList.size());

         // calc MDL value
         ArrayList<WindowLocation> occs = occList.get(0);
         double mdl = calcDL(occs);
         double mdl2 = calcDL(seeds.get(0), occs);

         int[] a = cri.getClassStats(iBestClass);

          /*System.out.printf( " Discovered class: %d (%s) -- CIDS: %d, %d, %d, %d apr: %.1f, %.1f, %.1f\n",
          iBestClass + 1, gdata.classes[iBestClass], a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION],
          100*cri.getClassAcc(iBestClass), 100*cri.getClassPrec(iBestClass), 100*cri
          .getClassRecall(iBestClass));*/

         //System.out.printf("%d: bestClass=%d  %.2f %.2f %.2f [%d,%d,%d,%d] R=%.3f  %.3f\n", wlen,
               //iBestClass, 100 * cri.getClassAcc(iBestClass), 100 * cri.getClassPrec(iBestClass), 100 * cri
                 //    .getClassRecall(iBestClass), a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION],
               //estR, mdl2);
      }

      lsWindow.setValue(lenOrig);
      btEstR.setSelected(bEstOrig);
      bStopAfterOneMotif = false;
      bVerbose = bVerboseOrig;
   }
   
   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();

      if (src == btDiscover) run();
      else if (src == btReal) dumpRealR();
      else if (src == btScanR) scanR();
      else if (src == btScanL) scanL();
      else if (src == btSaveLabels) saveLabels();
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();

      if (src == lsMotif){
         int iMotif = lsMotif.getValue() - 1;
         if (iMotif >= 0 && !occList.isEmpty()){
            lgraphAcc.setData(seeds.get(iMotif).first);
            lgraphGyr.setData(seeds.get(iMotif).first);
            select(occList.get(iMotif), true);
            // System.err.printf("Motif %d: %d occs\n", iMotif + 1, occList.get(iMotif).size());
         }
      }
   }
}
