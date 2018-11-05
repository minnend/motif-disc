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
import kdm.models.Gaussian1D;
import kdm.models.misc.*;
import kdm.util.*;
import kdm.util.SparseMatrix.MatrixEntry;

import static kdm.models.misc.ContRecInfo.*;

/**
 * Visual implementation of Chiu's algorithm.
 */
public class SubDimView extends AbstractDiscView implements ActionListener, ChangeListener
{
   protected LabeledSlider lsWindow, lsShow, lsPaaSegs, lsSaxSymbols, lsMinMotifs;
   protected JComboBox cbFind;
   protected JTextField tfDist, tfInput, tfProportion;
   protected JButton btDiscover, btReal, btScanR, btScanL, btAddWhite, btAddWalk, btResetData;
   protected JCheckBox btEstR;
   protected LineGraph lgraphAcc, lgraphGyr;

   // protected MetricSeq metseq = new SumSeqDist(new EuclideanFV(false), true);
   // protected MetricSeq metseq = new SumSeqDist(new AbsoluteDistFV(), true);
   protected MetricSeq metseq = new DTWSimple(.1);
   protected ArrayList<ArrayList<WindowLocation>> occList;
   protected ArrayList<Pair<Sequence, Sequence>> seeds;
   protected boolean bStopAfterOneMotif = false;
   protected boolean bVerbose = true;
   protected double[] density;
   protected ContRecInfo cri;
   protected double estR;

   public SubDimView(GlobalData gdata)
   {
      super(gdata);
      gdata.subdimView = this;

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
      lsShow = new LabeledSlider("Show Motif: %d", 1, 1);
      lsShow.setEnabled(false);
      lsShow.addChangeListener(this);
      rightp.add(lsShow);

      lsWindow.setValue(14);
      lsPaaSegs.setValue(4);
      lsSaxSymbols.setValue(5);
      lsMinMotifs.setValue(12);

      tfDist = new JTextField(String.format("%.3f", 0.85), 8);
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

      lgraphAcc = new LineGraph();
      lgraphAcc.setShowAll(true);
      lgraphAcc.setDimColor(0, Color.red);
      lgraphAcc.setDimColor(1, Color.green);
      lgraphAcc.setDimColor(2, Color.blue);
      rightp.add(Box.createVerticalStrut(8));
      rightp.add(new Constrainer(lgraphAcc, 200, 60));

      lgraphGyr = new LineGraph();
      lgraphGyr.setShowAll(true);
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
      p = new JPanel();
      btAddWhite = new JButton("White");
      btAddWhite.addActionListener(this);
      p.add(btAddWhite);
      btAddWalk = new JButton("Walk");
      btAddWalk.addActionListener(this);
      p.add(btAddWalk);
      leftp.add(p);
      btResetData = new JButton("Reset Data");
      btResetData.addActionListener(this);
      leftp.add(btResetData);
   }

   /**
    * Discover according to current GUI settings
    */
   protected void run()
   {
      TimerMS timer = new TimerMS();
      TimerMS timer2 = new TimerMS();

      int wlen = lsWindow.getValue();
      int nPAA = lsPaaSegs.getValue();
      int nSAX = lsSaxSymbols.getValue();
      int nSeqs = tseries.size();
      int nDims = tseries.get(0).getNumDims();
      int nTotal = gdata.nTotalSeriesLength;
      int nMinOccs = lsMinMotifs.getValue();
      boolean bEstR = btEstR.isSelected();
      double R = bEstR ? 0 : Double.parseDouble(tfDist.getText());

      // double MeanLinearThresh = 0.035; // TODO

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
      if (bVerbose)
      {
         int nRounds = nPAA * (nPAA - 1) / 2;
         System.err.printf("Building Collision Matrix (%d rounds)... ", nRounds);
      }
      timer.reset();
      HashMap<String, MyIntList> map = findExactMatches(nPAA, nSAX, wlen, 0, spans); // TODO: other dims
      System.err.printf("# exact matches: %d  (%dms)\n", map.size(), timer.time());
      {
         int nHits = 0;
         int n = 0;
         Iterator<String> it = map.keySet().iterator();
         while(it.hasNext())
         {
            String s = it.next();
            int x = map.get(s).size();
            if (x > 1)
            {
               n += x;
               if (x > 1) System.err.printf(" %s: %d\n", s, x);
               nHits += x * (x - 1) / 2;
            }
         }
         System.err.printf(" total number of entries: %d => %d hits\n", n, nHits);
      }
      timer.reset();
      int nCols = 2;
      SparseMatrix<int[]> mat = null;
      while(true)
      {
         try
         {
            System.err.printf("building collision matrix with %d columns\n", nCols);
            mat = buildCollisionMatrix(nPAA, nSAX, wlen, nCols, spans);
         } catch (ColMatrixException e)
         {
            nCols++;
            continue;
         }
         break;
      }
      if (mat == null)
      {
         System.err.printf("Error: failed to create collision matrix\n");
         return;
      }
      if (bVerbose) System.err.printf("done (%d, %dms).\n", mat.size(), timer.time());

      int n = nTotal * (nTotal - 1) / 2;
      int n2 = nAfterLinear * (nAfterLinear - 1) / 2;
      int hit = mat.size();
      if (bVerbose)
         System.err.printf(" hit: %d / %d = %.3f%%  (%d / %d = %.2f%%)\n", hit, n, 100.0 * hit / n, hit, n2,
               100.0 * hit / n2);

      // now that we have the collision matrix, we can start enumerating motifs
      if (bVerbose)
         System.err.printf("Enumerating motifs (R=%s)... ", bEstR ? "est" : String.format("%.2f", R));
      timer.reset();
      SpanList spanColl = new SpanList(0, nTotal - 1, true);
      occList = new ArrayList<ArrayList<WindowLocation>>();
      seeds = new ArrayList<Pair<Sequence, Sequence>>();
      lsShow.setEnabled(true);
      int iMotif = 0;
      int[] aSeed = new int[2];
      while(true)
      {
         SparseMatrix<int[]>.MatrixEntry me = findBestCollision(mat, spanColl);
         if (me == null) break;
         aSeed[0] = me.row();
         aSeed[1] = me.column();
         if (aSeed[0] < 0 || aSeed[1] < 0) break;

         assert (aSeed[0] < aSeed[1]);
         int[] hits = me.get();
         System.err.printf("hits (%d): ", hits.length);
         for(int i = 0; i < nDims; i++)
            System.err.printf("%d ", hits[i]);
         System.err.println();
         System.exit(0); // TODO: debug

         mat.remove(aSeed[0], aSeed[1]); // we never want to reprocess
         if (Math.abs(aSeed[0] - aSeed[1]) < wlen) continue; // no trivial matches

         // if (bVerbose) System.err.printf("aSeed: %d, %d (%d)\n", aSeed[0], aSeed[1], nHits);
         if (bEstR) estR = R = estimateR(aSeed, wlen) / wlen;

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
         if (dist > R)
         {
            // System.err .printf("Warning: seed pair too far apart (%d<->%d: %.1f)\n", aSeed[0], aSeed[1],
            // dist);
            spanColl.sub(aSeed[0]);
            spanColl.sub(aSeed[1]);
            continue;
         }

         ArrayList<WindowLocation> occs = null;
         switch(cbFind.getSelectedIndex()){
         case 0:
            occs = findOccs(wloc1, wloc2, R, spans);
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
         for(int i = 0; i < nOccs; i++)
         {
            WindowLocation wloc = occs.get(i);
            // System.err.printf("Occ %d: %s\n", i+1, wloc);
            spans[wloc.iSeries].sub(wloc.iStart - wlen + 1, wloc.getLastIndex());
            int indexA = wloc.getLinearIndex(gdata.seqLens);
            spanColl.sub(indexA - wlen + 1, indexA + wlen - 1);
         }

         if (occs != null && occs.size() >= nMinOccs)
         {
            occList.add(occs);
            seeds.add(new Pair<Sequence, Sequence>(seed1, seed2));
            double mdl = 0;// calcDL(occs); // TODO
            if (bVerbose)
               System.err.printf("Motif %d: %d occs, R=%.3f, DL=%.2f\n", iMotif + 1, occs.size(), R, mdl);
            iMotif++;
            if (bStopAfterOneMotif) break;
         }
      }
      if (bVerbose) System.err.printf("done (%dms).\n", timer.time());
      if (bVerbose) System.err.printf("Total time for discovery: %dms\n", timer2.time());
      scoreMotifs(occList, gdata.labData);
      lsShow.setValues(1, 1, occList.size());
   }

   /** Find exact matches between SAX strings */
   HashMap<String, MyIntList> findExactMatches(int nPAA, int nSAX, int wlen, int d, SpanList[] spans)
   {
      HashMap<String, MyIntList> map = new HashMap<String, MyIntList>();
      int nSeqs = tseries.size();
      double[] paa = new double[nPAA];
      char symbols[] = new char[nPAA];

      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = tseries.get(iSeq);
         double[] raw = seq.extractDim(d);
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            int index = Library.getIndexFromArrayOffset(iSeq, t, gdata.seqLens);

            // build the SAX string
            SAX.genPAA(raw, t, wlen, paa);
            for(int i = 0; i < nPAA; i++)
               symbols[i] = SAX.raw2sym(paa[i], nSAX);
            String s = new String(symbols);
            MyIntList list = map.get(s);
            if (list == null)
            {
               list = new MyIntList();
               map.put(s, list);
            }
            list.add(index);
         }
      }
      return map;
   }

   /**
    * build the multidim collision rows
    * 
    * @return [D][index] sum
    */
   protected int[][] buildCollisionRows(int nPAA, int nSAX, int wlen, int nCols, SpanList[] spans)
   {
      int nSeqs = tseries.size();
      int nDims = tseries.get(0).getNumDims();
      double[] paa = new double[nPAA];
      char symbols[][] = new char[nDims][nPAA];

      // build the buckets
      int nRounds = (int)Library.choose(nPAA, 2);
      ArrayList<ArrayList<HashMap<String, MyIntList>>> buckets = new ArrayList<ArrayList<HashMap<String, MyIntList>>>();
      for(int d = 0; d < nDims; d++)
      {
         ArrayList<HashMap<String, MyIntList>> a = new ArrayList<HashMap<String, MyIntList>>();
         buckets.add(a);
         for(int iRound = 0; iRound < nRounds; iRound++)
            a.add(new HashMap<String, MyIntList>());
      }

      // now we can fill the buckets      
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = tseries.get(iSeq);
         double[][] raw = seq.toSeqArray();
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            int index = Library.getIndexFromArrayOffset(iSeq, t, gdata.seqLens);

            // build the SAX string
            for(int d = 0; d < nDims; d++)
            {
               SAX.genPAA(raw[d], t, wlen, paa);
               for(int i = 0; i < nPAA; i++)
                  symbols[d][i] = SAX.raw2sym(paa[i], nSAX);
            }

            int iRound = 0;
            for(int iCol = 0; iCol < nPAA; iCol++)
               for(int jCol = iCol + 1; jCol < nPAA; jCol++)
               {
                  for(int d = 0; d < nDims; d++)
                  {
                     StringBuffer sb = new StringBuffer();
                     sb.append(symbols[d][iCol]);
                     sb.append(symbols[d][jCol]);
                     String s = sb.toString();
                     MyIntList list = buckets.get(d).get(iRound).get(s);
                     if (list == null)
                     {
                        list = new MyIntList();
                        buckets.get(d).get(iRound).put(s, list);
                     }
                     list.add(index);
                  }
                  iRound++;
               }           
         }
      }
      
      // now we can build the row sum array
      int[][] srow = new int[nDims][gdata.nTotalSeriesLength];
      for(int d = 0; d < nDims; d++)
      {
         ArrayList<HashMap<String, MyIntList>> array = buckets.get(d);
         for(HashMap<String, MyIntList> h : array)
         {
            Iterator<String> it = h.keySet().iterator();
            while(it.hasNext())
            {
               MyIntList list = h.get(it.next());
               int n = list.size();
               for(int i=0; i<n; i++)
               {
                  int x = list.get(i);
                  srow[d][x] += (n-1);
               }
            }
         }
      }
      return srow;
   }

   /** build the multidim collision matrix */
   protected SparseMatrix<int[]> buildCollisionMatrix(int nPAA, int nSAX, int wlen, int nCols,
         SpanList[] spans) throws ColMatrixException
   {
      int nSeqs = tseries.size();
      int nDims = tseries.get(0).getNumDims();
      double[] paa = new double[nPAA];
      char symbols[][] = new char[nDims][nPAA];

      int nMaxRounds = 60;
      double nMaxHitPercent = 0.02;
      int nPairs = (int)Library.choose(gdata.nTotalSeriesLength, 2);
      int nMaxHits = (int)Math.ceil(nMaxHitPercent * nPairs);

      // store the columns to use
      int nAllRounds = (int)Library.choose(nPAA, nCols);
      int nRounds = (nAllRounds > nMaxRounds ? nMaxRounds : nAllRounds);
      System.err.printf("#rounds=%d  (all=%d, max=%d)  max # hits: %d  (total: %d)\n", nRounds, nAllRounds,
            nMaxRounds, nMaxHits, nPairs);
      ArrayList<ArrayList<HashMap<String, MyIntList>>> buckets = new ArrayList<ArrayList<HashMap<String, MyIntList>>>();
      for(int d = 0; d < nDims; d++)
      {
         ArrayList<HashMap<String, MyIntList>> a = new ArrayList<HashMap<String, MyIntList>>();
         buckets.add(a);
         for(int iRound = 0; iRound < nRounds; iRound++)
            a.add(new HashMap<String, MyIntList>());
      }
      System.err.printf(" Finished building buckets (#rounds: %d)\n", nRounds);

      // now we can fill the buckets
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = tseries.get(iSeq);
         double[][] raw = seq.toSeqArray();
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            int index = Library.getIndexFromArrayOffset(iSeq, t, gdata.seqLens);

            // build the SAX string
            for(int d = 0; d < nDims; d++)
            {
               SAX.genPAA(raw[d], t, wlen, paa);
               for(int i = 0; i < nPAA; i++)
                  symbols[d][i] = SAX.raw2sym(paa[i], nSAX);
            }

            Integer[] items = new Integer[nPAA];
            for(int i = 0; i < nPAA; i++)
               items[i] = new Integer(i);
            int[] a = new int[nCols];
            CombinationGenerator gen = new CombinationGenerator(items, nCols);
            for(int iRound = 0; iRound < nRounds; iRound++)
            {
               Object[] b = (Object[])gen.nextElement();
               for(int i = 0; i < nCols; i++)
                  a[i] = (Integer)b[i];
               for(int d = 0; d < nDims; d++)
               {
                  StringBuffer sb = new StringBuffer();
                  for(int i = 0; i < nCols; i++)
                     sb.append(symbols[d][a[i]]);
                  String s = sb.toString();
                  MyIntList list = buckets.get(d).get(iRound).get(s);
                  if (list == null)
                  {
                     list = new MyIntList();
                     buckets.get(d).get(iRound).put(s, list);
                  }
                  list.add(index);
               }
            }
         }
      }
      System.err.printf(" finished filling buckets\n");

      // count number of hits
      boolean bTooBig = false;
      for(int d = 0; d < nDims; d++)
      {
         int nHits = 0;
         for(int iRound = 0; iRound < nRounds; iRound++)
         {
            Iterator<String> it = buckets.get(d).get(iRound).keySet().iterator();
            while(it.hasNext())
            {
               String s = it.next();
               MyIntList list = buckets.get(d).get(iRound).get(s);
               int n = list.size();
               nHits += n * (n - 1) / 2;
            }
         }
         System.err.printf(" Dim: %d  #hits: %d\n", d + 1, nHits);
         bTooBig = (bTooBig || nHits > nMaxHits);
      }
      if (bTooBig) throw new ColMatrixException(ColMatrixException.tooMany);

      // increment all collisions
      int nHits = 0;
      SparseMatrix<int[]> mat = new SparseMatrix<int[]>();
      for(int iRound = 0; iRound < nRounds; iRound++)
      {
         // for(int d = 0; d < nDims; d++)
         int d = 0; // TODO: debug
         {
            Iterator<String> it = buckets.get(d).get(iRound).keySet().iterator();
            while(it.hasNext())
            {
               String s = it.next();
               MyIntList list = buckets.get(d).get(iRound).get(s);
               int n = list.size();
               for(int i = 0; i < n; i++)
                  for(int j = i + 1; j < n; j++)
                  {
                     int a = list.get(i);
                     int b = list.get(j);
                     assert (a < b); // true due to construction method
                     int[] ii = mat.get(a, b);
                     if (ii == null)
                     {
                        ii = new int[nDims];
                        ii[d] = 1;
                        mat.set(a, b, ii);
                        nHits++;
                        if (nCols < nPAA && nHits > nMaxHits)
                           throw new ColMatrixException(ColMatrixException.tooMany);
                     }
                     else ii[d]++;
                  }
            }
         }
      }

      System.err.printf("finished collision matrix (nCols=%d, nHits=%d)\n", nCols, nHits);
      return mat;
   }

   /**
    * Find the matrix entry that is contained in the span list with the largest value
    * 
    * @param mat collision matrix
    * @param span allowed locations
    * @return index of largest allowed row or -1 if none found
    */
   protected SparseMatrix<int[]>.MatrixEntry findBestCollision(SparseMatrix<int[]> mat, SpanList span)
   {
      SparseMatrix<int[]>.MatrixEntry meBest = null;
      int xBest = -1;
      SpanIterator it = span.iterator();
      while(it.hasMore())
      {
         int i = it.next();
         HashMap<Integer, int[]> row = mat.getRow(i);
         if (row == null) continue;
         Iterator<Integer> jt = row.keySet().iterator();
         while(jt.hasNext())
         {
            int j = jt.next();
            SparseMatrix<int[]>.MatrixEntry me = mat.getEntry(i, j);
            int[] ii = me.get();
            int ix = Library.maxi(ii);
            if (meBest == null || ii[ix] > xBest)
            {
               meBest = me;
               xBest = ii[ix];
            }
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
      while(it.hasNext())
      {
         int i = it.next();
         if (coll[i] < 2) continue;
         if (iBest < 0 || coll[i] > coll[iBest]) iBest = i;
      }
      return iBest;
   }

   /**
    * Find all occurrences given a seed location and a max distance
    * 
    * @param wlSeed1 first base sequence that forms center of motif
    * @param wlSeed2 second base sequence that forms center of motif
    * @param R distance from seed within which other motif occs must lie
    * @param spans span list for each sequence specifying allowable motif (start) locations
    * @return list of occurrences
    */
   protected ArrayList<WindowLocation> findOccs(WindowLocation wlSeed1, WindowLocation wlSeed2, double R,
         SpanList[] spans)
   {
      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      assert (wlen == wlSeed2.length());
      Sequence seed1 = wlSeed1.getSeq(tseries);
      Sequence seed2 = wlSeed2.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = tseries.get(iSeq);
         int T = gdata.seqLens[iSeq];
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            Sequence win = seq.subseq(t, t + wlen);
            double dist1 = metseq.dist(seed1, win) / wlen;
            if (dist1 < R || metseq.dist(seed2, win) / wlen < R)
            {
               occs.add(new WindowLocation(iSeq, t, wlen));
               it.jump(t + wlen);
            }
         }
      }
      return occs;
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
      // TODO: doesn't work
      assert false : "doesn't work!";

      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      Sequence seed = wlSeed1.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         int ixd = 0, cut = 0;
         double d = Double.NaN;
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            Sequence win = tseries.get(iSeq).subseq(t, t + wlen);
            double r = metseq.dist(seed, win) / wlen;

            if (Double.isNaN(d))
            {
               d = r;
               ixd = t;
               cut = t + wlen;
            }
            else
            {
               if (t >= cut)
               {
                  if (d <= R) occs.add(new WindowLocation(iSeq, ixd, wlen));
                  d = r;
                  ixd = t;
                  cut = t + wlen;
               }
               else if (r > d)
               {
                  if (d <= R) occs.add(new WindowLocation(iSeq, ixd, wlen));
                  d = Double.NaN;
                  it.jump(ixd + wlen);
               }
            }
         }
         if (!Double.isNaN(d))
         {
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
      PriorityQueue<ScoredWindow> q = new PriorityQueue<ScoredWindow>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      assert (wlen == wlSeed2.length());
      Sequence seed1 = wlSeed1.getSeq(tseries);
      Sequence seed2 = wlSeed1.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = tseries.get(iSeq);
         int T = gdata.seqLens[iSeq];
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore())
         {
            int t = it.next();
            Sequence win = seq.subseq(t, t + wlen);
            double dist1 = metseq.dist(seed1, win) / wlen;
            double dist2 = Library.INF;
            if (dist1 < R || (dist2 = metseq.dist(seed2, win) / wlen) < R)
            {
               q.add(new ScoredWindow(iSeq, t, wlen, Math.min(dist1, dist2)));
               it.jump(t + wlen);
            }
         }
      }

      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      ScoredWindow win;
      scored: while((win = q.poll()) != null)
      {
         for(WindowLocation occ : occs)
            if (occ.overlaps(win)) continue scored;
         occs.add(win);
      }
      return occs;

   }

   /** Dump stats on the real distance between labeled motif members */
   protected void dumpRealR()
   {
      Iterator<String> it = gdata.labData.keySet().iterator();
      while(it.hasNext())
      {
         String sClass = it.next();
         ArrayList<Sequence> list = gdata.labData.get(sClass);
         int N = list.size();
         double[] R = new double[N * (N - 1) / 2];
         int ixR = 0;
         for(int i = 0; i < N; i++)
         {
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
      for(int i = 0; i < nSeqs; i++)
      {
         int nr = spans[i].getNumSpans();
         for(int j = 0; j < nr; j++)
            dl1 += Math.log(spans[i].getRange(j).length());
      }

      return 2 * dl1 + dl3;
   }

   /** Estimate R for the given seeds */
   protected double estimateR(int[] aSeed, int wlen)
   {
      double fProp = Double.parseDouble(tfProportion.getText()) / 100.0;
      int nSeqs = tseries.size();
      Sequence seed1, seed2 = null;
      WindowLocation wlSeed1, wlSeed2 = null;

      int[] ofs = Library.getArrayOffset(aSeed[0], gdata.seqLens, null);
      wlSeed1 = new WindowLocation(ofs[0], ofs[1], wlen);
      seed1 = tseries.get(ofs[0]).subseq(ofs[1], ofs[1] + wlen);
      if (aSeed.length > 1)
      {
         ofs = Library.getArrayOffset(aSeed[1], gdata.seqLens, null);
         wlSeed2 = new WindowLocation(ofs[0], ofs[1], wlen);
         seed2 = tseries.get(ofs[0]).subseq(ofs[1], ofs[1] + wlen);
      }
      int dj = 2;// (int)Math.ceil(wlen/10.0); // TODO
      MyDoubleList list = new MyDoubleList();
      WindowLocation wloc = new WindowLocation(-1, -1, wlen);
      for(int i = 0; i < nSeqs; i++)
      {
         int ixd = 0, cut = 0;
         double d = Double.NaN;
         wloc.iSeries = i;
         for(int j = 0; j + wlen <= gdata.seqLens[i]; j += dj)
         {
            wloc.iStart = j;
            Sequence win = tseries.get(i).subseq(j, j + wlen);
            if (wlSeed1.overlaps(wloc) || (wlSeed2 != null && wlSeed2.overlaps(wloc))) continue;
            double r1 = metseq.dist(seed1, win);
            double r2 = (seed2 != null ? metseq.dist(seed2, win) : r1);
            // double r = (r1 + r2) / 2;
            double r = Math.min(r1, r2);
            list.add(r); // TODO

            /*
             * if (Double.isNaN(d)) { d = r; ixd = j; cut = j + wlen; } else { if (j == cut) { list.add(d); d =
             * r; ixd = j; cut = j + wlen; } else if (r > d) { list.add(d); d = Double.NaN; assert (ixd + wlen -
             * 1 >= j) : String.format("ixd=%d %d j=%d", ixd, ixd + wlen - 1, j); j = ixd + wlen - 1; } }
             */
         }
         /*
          * if (!Double.isNaN(d)) { assert (ixd + wlen <= gdata.seqLens[i]); list.add(d); }
          */
      }

      if (list.isEmpty()) return Double.NaN;

      double[] rAll = list.toArray();

      int n = (int)Math.round(fProp * rAll.length);
      int iMax = Library.select(rAll, 0, rAll.length - 1, n);
      // System.err.printf("all=%d n=%d iMax=%d\n", rAll.length, n, iMax);
      double rmax = rAll[iMax];
      double[] r = new double[n];
      int ixr = 0;
      for(int i = 0; i < rAll.length; i++)
         if (rAll[i] <= rmax) r[ixr++] = rAll[i];
      assert (ixr == r.length) : String.format("ixr=%d n=%d", ixr, n);

      /*
       * { // TODO DSRaw saver = new DSRaw(); saver.save(new Sequence("foo", r),
       * String.format("/home/dminn/foo%d.txt", occList.size()+1)); }
       */

      double knee = Library.findKnee(r, true);

      return r[(int)Math.round(knee)];
   }

   /** try different R values */
   protected void scanR()
   {
      String s = tfInput.getText();
      StringTokenizer st = new StringTokenizer(s, " :");
      if (st.countTokens() != 3)
      {
         System.err.printf("Error: input string has incorrect format for scanR (<start>:<step>:<stop>)\n");
         return;
      }
      double start = Double.parseDouble(st.nextToken());
      double step = Double.parseDouble(st.nextToken());
      double stop = Double.parseDouble(st.nextToken());

      String sROrig = tfDist.getText();
      boolean bEstOrig = btEstR.isSelected();
      btEstR.setSelected(false);

      // bStopAfterOneMotif = true;
      boolean bVerboseOrig = bVerbose;
      bVerbose = false;

      for(double r = start; (r - stop) < 0.000001; r += step)
      {
         // System.out.printf("Discovering first motif for R=%.3f...\n", r);
         tfDist.setText(String.format("%.3f", r));
         System.err.printf("%.2f  ", r);
         run();
         /*
          * if (cri == null || cri.getNumSpots() == 0) { System.out.println(" No motifs."); continue; } //
          * which class is best? int iBestClass = 0; for(int i = 1; i < gdata.getNumClasses(); i++) { if
          * (cri.getClassStats(i)[CORRECT] > cri.getClassStats(iBestClass)[CORRECT]) iBestClass = i; }
          * 
          * assert (seeds.size() == 1) : String.format("# seeds: %d", seeds.size()); assert (occList.size() ==
          * 1) : String.format("# occ lists: %d\n", occList.size()); // calc MDL value ArrayList<WindowLocation>
          * occs = occList.get(0); double mdl = calcDL(occs);
          * 
          * int[] a = cri.getClassStats(iBestClass); // System.out.printf( " Discovered class: %d (%s) --
          * CIDS: %d, %d, %d, %d apr: %.1f, %.1f, %.1f\n", // iBestClass + 1, gdata.classes[iBestClass],
          * a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION], // 100*cri.getClassAcc(iBestClass),
          * 100*cri.getClassPrec(iBestClass), 100*cri // .getClassRecall(iBestClass));
          * 
          * //System.out.printf("%.6f %.2f %.2f %.2f %.4f\n", r, 100 * cri.getClassAcc(iBestClass), 100 * cri //
          * .getClassPrec(iBestClass), 100 * cri.getClassRecall(iBestClass), mdl);
          */}

      btEstR.setSelected(bEstOrig);
      tfDist.setText(sROrig);
      bStopAfterOneMotif = false;
      bVerbose = bVerboseOrig;
   }

   /**
    * Calculate the expected number of hits in the collision matrix, based on Chiu's paper
    * 
    * @param nStrings number of strings being compared (k)
    * @param wlen length of window (w); note this is the length of the string, not the raw window size
    * @param nMaxErrs maximum number of allowed errors per window (d)
    * @param projLen (t) length of projected string (i.e., number of "columns" selected)
    * @param nAlpha (a) size of the alphabet
    * @return expected number of occurrences per round of lsh
    */
   public double calcExpected(int nStrings, int wlen, int nMaxErrs, int projLen, int nAlpha)
   {
      if (nMaxErrs > wlen) nMaxErrs = wlen;

      // calc k from Chiu's paper
      double x = 0;
      for(int i = 0; i <= nMaxErrs; i++)
      {
         double x1 = Math.pow(1.0 - (double)i / wlen, projLen);
         long x2 = Library.choose(wlen, i);
         double x3 = Math.pow((nAlpha - 1.0) / nAlpha, i);
         double x4 = Math.pow(1.0 / nAlpha, wlen - i);
         x += x1 * x2 * x3 * x4;
      }
      x *= Library.choose(nStrings, 2);
      return x;
   }

   /** try different motif lengths */
   protected void scanL()
   {
      String s = tfInput.getText();
      StringTokenizer st = new StringTokenizer(s, " :");
      int nToks = st.countTokens();
      int start, step, stop;
      if (nToks == 2)
      {
         start = Integer.parseInt(st.nextToken());
         step = 1;
         stop = Integer.parseInt(st.nextToken());
      }
      else if (nToks == 3)
      {
         start = Integer.parseInt(st.nextToken());
         step = Integer.parseInt(st.nextToken());
         stop = Integer.parseInt(st.nextToken());
      }
      else
      {
         System.err
               .printf("Error: input string has incorrect format for scanR (<start>:<step>:<stop> or <start>:<stop>)\n");
         return;
      }

      int lenOrig = lsWindow.getValue();
      boolean bEstOrig = btEstR.isSelected();
      btEstR.setSelected(true);

      bStopAfterOneMotif = true;
      boolean bVerboseOrig = bVerbose;
      bVerbose = false;

      for(int wlen = start; wlen <= stop; wlen += step)
      {
         lsWindow.setValue(wlen);
         run();

         if (cri == null || cri.getNumSpots() == 0)
         {
            System.out.println(" No motifs.");
            continue;
         }

         // which class is best?
         int iBestClass = 0;
         for(int i = 1; i < gdata.getNumClasses(); i++)
         {
            if (cri.getClassStats(i)[CORRECT] > cri.getClassStats(iBestClass)[CORRECT]) iBestClass = i;
         }

         assert (seeds.size() == 1) : String.format("# seeds: %d", seeds.size());
         assert (occList.size() == 1) : String.format("# occ lists: %d\n", occList.size());

         // calc MDL value
         ArrayList<WindowLocation> occs = occList.get(0);

         int[] a = cri.getClassStats(iBestClass);

         // System.out.printf( " Discovered class: %d (%s) -- CIDS: %d, %d, %d, %d apr: %.1f, %.1f, %.1f\n",
         // iBestClass + 1, gdata.classes[iBestClass], a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION],
         // 100*cri.getClassAcc(iBestClass), 100*cri.getClassPrec(iBestClass), 100*cri
         // .getClassRecall(iBestClass));

         System.out.printf("%d: bestClass=%d  %.2f %.2f %.2f [%d,%d,%d,%d] R=%.3f\n", wlen, iBestClass,
               100 * cri.getClassAcc(iBestClass), 100 * cri.getClassPrec(iBestClass), 100 * cri
                     .getClassRecall(iBestClass), a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION],
               estR);
      }

      lsWindow.setValue(lenOrig);
      btEstR.setSelected(bEstOrig);
      bStopAfterOneMotif = false;
      bVerbose = bVerboseOrig;
   }

   protected void addWhiteNoise()
   {
      System.err.printf("Adding white noise channel (current: %dD)\n", tseries.get(0).getNumDims());

      Gaussian1D gauss = new Gaussian1D(0, 10);
      int nSeqs = tseries.size();
      ArrayList<Sequence> seqs = tseries;
      tseries = new ArrayList<Sequence>();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = seqs.get(iSeq);
         int T = seq.length();
         Sequence noise = new Sequence("noise", seq.getFreq());
         for(int t = 0; t < T; t++)
            noise.add(gauss.sample());
         tseries.add(seq.addDims(noise));
      }
      System.err.printf(" New dimensionality: %dD (%d)\n", tseries.get(0).getNumDims(), gdata.getNumDims());
   }

   protected void addWalkNoise()
   {
      System.err.printf("Adding random walk channel (current: %dD)\n", tseries.get(0).getNumDims());

      Gaussian1D gauss = new Gaussian1D(0, 10);
      int nSeqs = tseries.size();
      ArrayList<Sequence> seqs = tseries;
      tseries = new ArrayList<Sequence>();
      for(int iSeq = 0; iSeq < nSeqs; iSeq++)
      {
         Sequence seq = seqs.get(iSeq);
         int T = seq.length();
         Sequence noise = new Sequence("noise", seq.getFreq());
         double v = (Library.random() - 0.5) * 10.0;
         for(int t = 0; t < T; t++)
         {
            noise.add(new FeatureVec(1, v));
            double dv = 0.1 + 2.0 * Library.random();
            if (Library.random() < 0.5) v -= dv;
            else v += dv;
         }
         tseries.add(seq.addDims(noise));
      }
      System.err.printf(" New dimensionality: %dD (%d)\n", tseries.get(0).getNumDims(), gdata.getNumDims());
   }

   protected void resetData()
   {
      tseries = new ArrayList<Sequence>();
      for(Sequence seq : gdata.tseries)
         tseries.add(seq);
      System.err.printf("Reset data: %d series, %dD\n", tseries.size(), gdata.getNumDims());
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();

      if (src == btDiscover) run();
      else if (src == btReal) dumpRealR();
      else if (src == btScanR) scanR();
      else if (src == btScanL) scanL();
      else if (src == btAddWhite) addWhiteNoise();
      else if (src == btAddWalk) addWalkNoise();
      else if (src == btResetData) resetData();
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();

      if (src == lsShow)
      {
         int iMotif = lsShow.getValue() - 1;
         if (iMotif >= 0)
         {
            lgraphAcc.setData(seeds.get(iMotif).first);
            lgraphGyr.setData(seeds.get(iMotif).first);
            select(occList.get(iMotif), true);
            // System.err.printf("Motif %d: %d occs\n", iMotif + 1, occList.get(iMotif).size());
         }
      }
   }
}
