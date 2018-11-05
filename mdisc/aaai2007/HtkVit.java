package mdisc.aaai2007;

import kdm.data.*;
import kdm.io.*;
import kdm.io.DataLoader.*;
import kdm.io.DataSaver.*;
import kdm.tools.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.util.*;
import kdm.mlpr.htk.*;
import kdm.gui.*;

import java.util.*;
import java.io.*;

import javax.swing.*;
import java.awt.*;

/** test interaction with HTK continuous recognition */
public class HtkVit
{
   public static void main(String[] args) throws Exception
   {
      final String BG = "bg";

      //String sFile = "/home/dminn/research/kdm/data/ah-data/exercise/sub8/accgyr.def";
      String sFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/dim13/ae.def";
      //String sFile = "/home/dminn/research/kdm/data/speech/tidigits/brett/5spkrs/ae/ae.def";
      //String sFile = "/home/dminn/research/kdm/data/asl-joseph/asl16D.def";
      //String sFile = "/home/dminn/research/kdm/data/asl-joseph/aslxyd.def";
      TreeMap<String, ArrayList<Sequence>> labData = LabeledDataLoader.load(new File(sFile));
      ArrayList<Sequence> data = LabeledDataLoader.tseries;
      ArrayList<String> sentences = LabeledDataLoader.sentences;
      System.err.printf("Found %d series in \"%s\"\n", data.size(), sFile);
      
      PrintWriter out;
      String cmd;

      ArrayList<AbstractHMM> hmms = new ArrayList<AbstractHMM>();
      int nDims = data.get(0).getNumDims();

      // create the background model
      HmmLR hmmBG = new HmmLR(1, nDims);
      hmmBG.setName(BG);
      hmmBG.setPiLeave(0, Math.log(0.2));
      GaussianDiagonal dg = (GaussianDiagonal)hmmBG.getState(0);
      GaussianDyn1D gmd[] = new GaussianDyn1D[nDims];
      for(int i = 0; i < nDims; i++)
         gmd[i] = new GaussianDyn1D();
      for(Sequence seq : data){
         int T = seq.length();
         for(int t = 0; t < T; t++)
            for(int d = 0; d < nDims; d++)
               gmd[d].add(seq.get(t, d), false);
      }
      FeatureVec fvMean = new FeatureVec(nDims);
      FeatureVec fvVar = new FeatureVec(nDims);
      for(int i = 0; i < nDims; i++){
         gmd[i].update();
         fvMean.set(i, gmd[i].getMean());
         fvVar.set(i, gmd[i].getVar());
      }
      dg.setMean(fvMean);
      dg.setVar(fvVar);
      System.err.printf(" mean: %s\n", fvMean);
      System.err.printf(" var: %s\n", fvVar);
      hmms.add(hmmBG);

      // create HMM for each class
      TimerMS timer = new TimerMS();
      int nClasses = labData.size();
      String[] classes = new String[nClasses];
      HashMap<String, Integer> c2i = new HashMap<String, Integer>();
      Iterator<String> it = labData.keySet().iterator();
      for(int ic = 0; ic < nClasses; ic++){
         classes[ic] = it.next();
         c2i.put(classes[ic], ic);
         String name = classes[ic].replace(' ', '_');
         ArrayList<Sequence> seqs = labData.get(classes[ic]);

         int vmin = seqs.get(0).length();
         int vmax = vmin;
         int vmean = vmin;
         int nSeqs = seqs.size();
         for(int i = 1; i < nSeqs; i++){
            int len = seqs.get(i).length();
            if (len < vmin) vmin = len;
            else if (len > vmax) vmax = len;
            vmean += len;
         }
         vmean = (vmean + nSeqs - 1) / nSeqs;

         // int nStates = Library.max((int)Math.round(vmean*.5), (int)Math.round(vmin*.75), 3);
         int nStates = vmin * 2 - 1; // with skip = 1         
         //int nStates = vmin; // with skip = 0
         HmmLR hmm = new HmmLR(nStates, 1, nDims);
         hmm.setName(name);  
         // TODO overlap method
         //hmm.init_segk(seqs);
         hmm.init_segk_overlap(seqs, 0.25);
         hmm.train_bw(seqs);
         double[][] tran = hmm.getFullTransMatrix();
         double meanTran = 0;
         int nTran = 0;
         for(int i = 0; i < tran.length - 1; i++){
            double p = tran[i][i];
            if (!Double.isNaN(p)){
               meanTran += Math.exp(p);
               nTran++;
            }
         }
         if (nTran > 0) meanTran = Math.max(meanTran / nTran, 0.5);
         else meanTran = 0.5;         
         hmm.setPiLeave(nStates - 1, Math.log(meanTran));
         System.err.printf(" %s (%d): %d, %d, %d (%d, %.2f)\n", classes[ic], seqs.size(), vmin, vmean, vmax,
               nStates, meanTran);
         hmms.add(hmm);
      }
      System.err.printf("Finished training HMMs (%dms)\n", timer.time());

      // perform continuous recognition
      System.err.print("Running continuous recognizer... ");
      timer.reset();
      ContRecRet crr = HTK.contRec(hmms, data);
      if (crr == null){
         System.err.println("Error: continuous recognition failed.");
         System.err.printf("\nError Stream:\n----------------------------------\n%s", crr.sErr);
         System.err.printf("\nOutput Stream:\n----------------------------------\n%s", crr.sOut);
         System.exit(1);
      }
      System.err.printf("done (%dms).\n", timer.time());
      
      // display label view
      /*ContRecViewer crv = new ContRecViewer(LabeledDataLoader.marks, crr.markupSets, new LabelEditor() {
         public String adjustLabel(String sLabel)
         {
            return sLabel.replace('_', ' ');
         }

         public boolean keepLabel(String sLabel)
         {
            return !sLabel.equals(BG);
         }
      });
      JFrame frame = new JFrame("Continuous Recognition Results");
      frame.setSize(Library.getScreenSize().width, 200);
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.getContentPane().setLayout(new BorderLayout());
      frame.getContentPane().add(crv, BorderLayout.CENTER);
      frame.setVisible(true);*/      

      // dump result stats
      ContRecInfo cri = new ContRecInfo(labData.size());
      it = labData.keySet().iterator();
      while(it.hasNext()){
         String sClass = it.next();
         ArrayList<Sequence> seqs = labData.get(sClass);
         cri.nLabeledWords += seqs.size();
         for(Sequence seq : seqs)
            cri.nLabeledFrames += seq.length();
      }

      // add the detected spots to the
      int nMarkupSets = crr.markupSets.size();
      for(int iSeq = 0; iSeq < nMarkupSets; iSeq++){
         MarkupSet marks = crr.markupSets.get(iSeq);
         int nMarks = marks.size();
         for(int i = 0; i < nMarks; i++){
            TimeMarker tm = marks.get(i);
            if (tm.getTag().equals(BG)) continue;
            String sClass = tm.getTag().replace('_', ' ');
            WindowLocation wloc = new WindowLocation(iSeq, tm.getStartIndex(), (int)tm.length());
            WordSpot spot = new WordSpot(wloc, 0.0, c2i.get(sClass));
            cri.add(spot);
         }
      }

      int[][] confm = cri.scoreWordSpotWordsFull(labData);
      cri.scoreWordSpotFrames(data, labData);

      // cri.dump();
      cri.showResults("Results", classes);
      // SupTest.dumpConfMatrix(confm, "Word Confusion Matrix", classes, null, false);
   }
}
