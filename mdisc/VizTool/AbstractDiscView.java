package mdisc.VizTool;

import static kdm.models.misc.ContRecInfo.CORRECT;
import static kdm.models.misc.ContRecInfo.DELETION;
import static kdm.models.misc.ContRecInfo.INSERTION;
import static kdm.models.misc.ContRecInfo.NUMBER;
import static kdm.models.misc.ContRecInfo.SUBSTITUTION;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import javax.swing.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.io.MSGeneral;
import kdm.tools.SupTest;
import kdm.util.*;
import kdm.metrics.*;
import kdm.models.misc.ContRecInfo;
import kdm.models.misc.WordSpot;

/**
 * abstract base class for discovery views. this class setups up a border layout with the discrete sequences
 * in the lower portion of the panel.
 */
public abstract class AbstractDiscView extends JPanel
{
   protected GlobalData gdata;
   protected ArrayList<Sequence> tseries;
   protected ArrayList<DiscreteSeq> qseries;
   protected int nSymbols;
   protected Color[] colors;
   protected QuantDataView[] qv;
   protected SpanList[] spanAvail;
   protected ContRecInfo cri;
   public boolean bVerbose = true;

   public AbstractDiscView(GlobalData gdata)
   {
      super(new BorderLayout());
      this.gdata = gdata;
      tseries = gdata.tseries;
      qseries = gdata.qseries;
      nSymbols = gdata.nSymbols;

      // bottom panel
      if (qseries != null){
         JPanel p = new JPanel(new GridFillLayout(qseries.size(), 1));
         p.setOpaque(true);
         p.setBackground(Color.black);
         colors = Library.generateColors(nSymbols);
         qv = new QuantDataView[qseries.size()];
         for(int i = 0; i < qseries.size(); i++){
            qv[i] = new QuantDataView(qseries.get(i), colors, true, gdata.nAvgSeriesLength > 300 ? 0 : 1);
            qv[i].setToolTipText(String.format("Series %d: %s", i + 1, gdata.getSentence(i)));
            p.add(qv[i]);
         }
         Constrainer qcons = new Constrainer(new JScrollPane(p), 0, Integer.MAX_VALUE, 0, 300);
         this.add(qcons, BorderLayout.SOUTH);
      }

      resetAvail();
   }

   /** Highlight all occurences in the given occurence list */
   protected void select(ArrayList<? extends WindowLocation> occs, boolean bClear)
   {
      if (bClear) for(int i = 0; i < qv.length; i++)
         qv[i].clear();
      for(WindowLocation wloc : occs)
         qv[wloc.iSeries].highlight(wloc.getFirstIndex(), wloc.getLastIndex());
   }

   /** reset visualization so that all locations are available */
   protected void resetAvail()
   {
      if (qseries==null) return;
      // setup the span lists -- to start, everything is available
      spanAvail = new SpanList[qseries.size()];
      for(int i = 0; i < spanAvail.length; i++){
         spanAvail[i] = new SpanList(0, qseries.get(i).length() - 1, true);
         qv[i].clear();
      }
   }

   /**
    * Find all occurrences given a seed location and a max distance. The algorithm assumes that all positions
    * contained in the span lists are valid starting positions. It does *not* verify this assumption.
    * 
    * @param wlSeed1 first base sequence that forms center of motif
    * @param wlSeed2 second base sequence that forms center of motif (can be null)
    * @param R distance from seed within which other motif occs must lie
    * @param spans span list for each sequence specifying allowable motif (start) locations
    * @param metseq metric to use for comparing sequences
    * @param bNorm normalize distance by window length?
    * @return list of occurrences
    */
   protected ArrayList<WindowLocation> findOccs(WindowLocation wlSeed1, WindowLocation wlSeed2, double R,
         SpanList[] spans, MetricSeq metseq, boolean bNorm)
   {
      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      assert (wlSeed2 == null || wlen == wlSeed2.length());
      if (bNorm) R *= wlen;
      Sequence seed1 = wlSeed1.getSeq(tseries);
      Sequence seed2 = wlSeed2 == null ? null : wlSeed2.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = tseries.get(iSeq);
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            Sequence win = seq.subseq(t, t + wlen); // TODO don't explicitly extract subseq
            double dist1 = metseq.dist(seed1, win);
            double dist2 = (seed2 == null ? Library.INF : metseq.dist(seed2, win));
            if (dist1 < R || dist2 < R){
               occs.add(new WindowLocation(iSeq, t, wlen));
               it.jump(t + wlen);
            }
         }
      }
      return occs;
   }

   /**
    * Find all occurrences given a seed location and a max distance, while being careful not to overlap an
    * unavailable position
    * 
    * @param wlSeed1 first base sequence that forms center of motif
    * @param wlSeed2 second base sequence that forms center of motif (can be null)
    * @param R distance from seed within which other motif occs must lie
    * @param spans span list for each sequence specifying allowable motif (start) locations
    * @param metseq metric to use for comparing sequences
    * @param bNorm normalize distance by window length?
    * @return list of occurrences
    */
   protected ArrayList<WindowLocation> findOccsSafe(WindowLocation wlSeed1, WindowLocation wlSeed2,
         double R, SpanList[] spans, MetricSeq metseq, boolean bNorm)
   {
      Range r = new Range();
      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
      int nSeqs = tseries.size();
      int wlen = wlSeed1.length();
      assert (wlen > 0) : wlen;
      assert (wlSeed2 == null || wlen == wlSeed2.length());
      if (!bNorm) R /= wlen;
      Sequence seed1 = wlSeed1.getSeq(tseries);
      Sequence seed2 = wlSeed2 == null ? null : wlSeed2.getSeq(tseries);
      for(int iSeq = 0; iSeq < nSeqs; iSeq++){
         Sequence seq = tseries.get(iSeq);
         int T = seq.length();
         SpanIterator it = spans[iSeq].iterator();
         while(it.hasMore()){
            int t = it.next();
            r.set(t, t + wlen - 1);
            if (r.b >= T) break; // don't run past the end
            if (!spans[iSeq].contains(r)) continue;
            Sequence win = seq.subseq(t, t + wlen); // TODO don't explicitly extract subseq
            double dist1 = metseq.dist(seed1, win) / wlen;
            double dist2 = seed2 == null ? Library.INF : metseq.dist(seed2, win) / wlen;
            if (dist1 < R || dist2 < R){
               double d = Math.min(dist1, dist2);
               occs.add(new WindowLocation(iSeq, t, wlen));
               it.jump(t + wlen);
            }
         }
      }
      return occs;
   }

   /**
    * Compute an objective evaluation of the motifs
    */
   public void scoreMotifs(ArrayList<ArrayList<WindowLocation>> motifs,
         TreeMap<String, ArrayList<Sequence>> labData)
   {
      /*
       * { System.err.printf("#motifs: %d\n", motifs.size()); for(int i=0; i<motifs.size(); i++) {
       * System.err.printf("%d: %d\n", i+1, motifs.get(i).size()); } }
       */

      if (labData == null) return;
      if (motifs == null || motifs.isEmpty()){
         System.err.println("Warning: no motifs!");
         return;
      }

      TimerMS timer = new TimerMS();
      int nClasses = labData.size();
      int nMotifs = motifs.size();

      // setup the continuous rec info structure
      cri = new ContRecInfo(gdata.getNumClasses());
      Iterator<String> itClass = labData.keySet().iterator();
      for(int iClass = 0; iClass < labData.size(); iClass++){
         ArrayList<Sequence> subs = labData.get(itClass.next());
         cri.nLabeledWords += subs.size();
         for(Sequence sub : subs)
            cri.nLabeledFrames += sub.length();
      }

      // add the discovered motifs
      for(int iMotif = 0; iMotif < nMotifs; iMotif++){
         ArrayList<WindowLocation> occs = motifs.get(iMotif);
         for(WindowLocation occ : occs){
            WordSpot spot = new WordSpot(occ, 1.0, iMotif);
            cri.add(spot);
         }
      }

      if (bVerbose) System.err.print("Computing initial performance stats... ");
      timer.reset();
      int[][] confm = cri.scoreWordSpotWordsFast(labData);
      cri.scoreWordSpotFrames(tseries, labData);
      if (bVerbose) System.err.printf("done (%dms).\n", timer.time());

      // SupTest.dumpConfMatrix(confm, "Motif Confusion Matrix", gdata.classes, gdata.nExamples);

      // calc best class permutation and remap spots
      if (bVerbose) System.err.print("Calculating best motif <-> class mapping... ");
      timer.reset();
      short[] map = SupTest.calcBestMappingFromConfMatrix(confm, false);
      cri.remapSpotsOld(map);
      if (bVerbose) System.err.printf("done (%dms).\n", timer.time());

      // now remove spots that don't correspond to known classes
      cri.removeUnknown(gdata.getNumClasses()); // TODO

      // recalc stats with permuted classes
      if (bVerbose) System.err.print("Recomputing performance stats... ");
      timer.reset();
      confm = cri.scoreWordSpotWordsFast(labData);

      cri.scoreWordSpotFrames(tseries, labData); // TODO: should just be able to move previous numbers
      // around...
      if (bVerbose) System.err.printf("done (%dms).\n", timer.time());

      if (bVerbose) System.err.println(cri.getDumpString());
      
      if (bVerbose){
         System.err.println("Event-based Accuracy, Precision, Recall (CIDS)...");
         for(int i = 0; i < nClasses; i++){
            int[] a = cri.getClassStats(i);
            System.err.printf(" Class %d: %.1f, %.1f, %.1f  (%s) [%d,%d,%d,%d]\n", i + 1, 100.0 * cri
                  .getClassAcc(i), 100.0 * cri.getClassPrec(i), 100.0 * cri.getClassRecall(i),
                  gdata.classes[i], a[CORRECT], a[INSERTION], a[DELETION], a[SUBSTITUTION]);
         }
      }

      if (bVerbose) System.err.print("Overall Event-based Accuracy, Precision, Recall: ");
      System.err.printf(" %.1f, %.1f, %.1f  ", 100.0 * cri.getWordAcc(), 100.0 * cri.getWordPrec(),
            100.0 * cri.getWordRecall());
      if (bVerbose) System.err.println();
      int[] a = cri.getRawWordCounts();
      if (bVerbose)
         System.err.printf("  CIDS: %d, %d, %d, %d   N=%d\n", a[CORRECT], a[INSERTION], a[DELETION],
               a[SUBSTITUTION], a[NUMBER]);
      if (bVerbose) System.err.print("Overall Frame-based Accuracy, Precision, Recall: ");
      System.err.printf(" %.1f, %.1f, %.1f\n", 100.0 * cri.getFrameAcc(), 100.0 * cri.getFramePrec(),
            100.0 * cri.getFrameRecall());
      a = cri.getRawFrameCounts();
      if (bVerbose)
         System.err.printf("  CIDS: %d, %d, %d, %d   N=%d\n", a[CORRECT], a[INSERTION], a[DELETION],
               a[SUBSTITUTION], a[NUMBER]);

      if (bVerbose){
         System.err.flush();
         // Library.sleep(300); // try to keep STDERR,STDOUT separate
         SupTest.dumpConfMatrix(confm, "Motif Confusion Matrix", gdata.classes, gdata.nExamples, false);
      }
   }

   protected void saveLabels()
   {
      JFileChooser fc = Library.buildFileChooser("labels", "Label Files (*.labels)");
      if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION){
         File f = fc.getSelectedFile();
         String sPath = Library.getPath(f.getAbsolutePath());
         String sTitle = Library.getTitle(f.getAbsolutePath());
         String sExt = Library.getExt(f.getAbsolutePath());
         if (sExt == null || sExt.length() == 0) sExt = "labels";

         ArrayList<MarkupSet> labels = cri.getLabels(gdata.classes);
         int nMarks = labels.size();
         MSGeneral saver = new MSGeneral();
         for(int i = 0; i < nMarks; i++){
            String sFile = String.format("%s%s_%03d.%s", sPath, sTitle, i + 1, sExt);
            MarkupSet marks = labels.get(i);
            if (!saver.save(marks, sFile)){
               System.err.printf("Error: failed to save label file\n (%s)\n", sFile);
               JOptionPane.showMessageDialog(this, "Failed to save label files", "Save Failed",
                     JOptionPane.ERROR_MESSAGE);
               return;
            }
         }

         String sFileDef = String.format("%s%s.def", sPath, sTitle);
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
         } catch (IOException ioe){
            System.err.printf("Error: failed to save def file\n (%s)\n", sFileDef);
            JOptionPane.showMessageDialog(this, "Failed to save def file", "Save Failed",
                  JOptionPane.ERROR_MESSAGE);
            return;
         }

         JOptionPane.showMessageDialog(this, String.format("Saved %d label files\n (with %s)", nMarks,
               Library.getFileName(sFileDef)), "Save Successful", JOptionPane.INFORMATION_MESSAGE);
      }
   }
}
