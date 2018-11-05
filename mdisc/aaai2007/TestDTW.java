package mdisc.aaai2007;

import java.io.*;
import java.util.*;

import kdm.data.*;
import kdm.io.*;
import kdm.metrics.*;
import kdm.util.*;

/** run some tests on DTW and lower bounding */
public class TestDTW
{
   public static void main(String[] args){
      TreeMap<String, ArrayList<Sequence>> data = LabeledDataLoader.load(new File(args[0]));
      if (data == null) System.err.printf("DiscApp) Failed to load labeled data\n (%s)\n", args[0]);
      
      // find the longest sequence
      int nLongest = 0;
      Iterator<String> it = data.keySet().iterator();
      while(it.hasNext()){
         for(Sequence seq : data.get(it.next()))
            if (seq.length() > nLongest) nLongest = seq.length();
      }
      
      MetricSeq metseq = new DTWSimple(0.1, MetricSeq.LengthPrep.none); 
      
      // build list of all resampled sequences with LBI
      ArrayList<Sequence> seqs = new ArrayList<Sequence>();
      ArrayList<LBInfo> lbi = new ArrayList<LBInfo>();
      it = data.keySet().iterator();
      while(it.hasNext()){
         for(Sequence seq : data.get(it.next())){
            seq = seq.resample(nLongest);
            seqs.add(seq);
            lbi.add(metseq.calcLBInfo(seq));
         }
      }
      
      // now run the all-pairs test
      TimerMS tt = new TimerMS();
      TimerMS timer = new TimerMS();
      long tmDist = 0, tmLB = 0;
      int nSeqs = seqs.size();
      for(int i=0; i<nSeqs; i++){
         Sequence seqA = seqs.get(i);
         LBInfo lbiA = lbi.get(i);
         //System.err.printf("i=%d/%d\n", i+1, nSeqs);
         for(int j=0; j<nSeqs; j++){
            
            Sequence seqB = seqs.get(j);
            timer.reset();
            double dist = metseq.dist(seqA, seqB);
            tmDist += timer.time();
            timer.reset();
            double lb = metseq.lowerBound(lbiA, seqB);
            tmLB += timer.time();
            assert(lb<=dist);
         }
      }
      
      System.err.printf("Time (ms): dist=%d  lb=%d  total=%d\n", tmDist, tmLB, tt.time());
   }
}
