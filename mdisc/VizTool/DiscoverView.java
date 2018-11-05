package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.border.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.tools.SupTest.HmmTrain;
import kdm.util.*;
import kdm.mlpr.*;
import kdm.mlpr.suffix_tree.*;
import kdm.models.*;
import kdm.io.*;
import mdisc.io.*;
import mdisc.data.*;

import org.apache.commons.collections.primitives.*;

public class DiscoverView extends AbstractDiscView implements ActionListener, MouseListener
{
   public static final double MIN_FOR_AFFIX = 0.60;
   
   public static final String[] sPresets = new String[]{ "Custom", "Exercise (sub8)", "Speech (ae,10ms)" };
   public static final int iCustom = 0;
   public static final int iExercise = 1;
   public static final int iSpeech = 2;
   
   protected JButton btDisc, btSave, btLoad, btTimeMerge, btTimeGrow, btSplit;
   protected JButton btAutoSplit, btGroup, btAutoGroup, btShowAvail, btAutoDisc;
   protected JTextField tfMDL, tfGrow, tfSplit;
   protected JComboBox cbPreset;
   protected HashView hashView;
   protected JPanel pCenter;
   protected VerticalScrollPanel pEntries;
   protected HashViewEntry prevHVE;
   protected ArrayList<ArrayList<WindowLocation>> motifs;
   protected AgglomSeqCluster asc;
   protected final JFileChooser fc = new JFileChooser();  
   protected Thread threadDisc;
   protected double mdlThresh = -200;//-500;
   protected double growThresh = 1.2;
   protected double splitThresh = 0.45;
   
   
   public DiscoverView(GlobalData gdata)
   {
      super(gdata);      
      gdata.discoverView = this;
      hashView = gdata.hashView;
      JPanel p;

      btDisc = new JButton("Discover!");
      btDisc.addActionListener(this);
      btSave = new JButton("Save");
      btSave.addActionListener(this);
      btLoad = new JButton("Load");
      btLoad.addActionListener(this);
      btTimeMerge= new JButton("Temporal Merge");
      btTimeMerge.addActionListener(this);
      btTimeGrow = new JButton("Temporal Grow");
      btTimeGrow.addActionListener(this);
      btSplit= new JButton("Split Info");
      btSplit.addActionListener(this);
      btAutoSplit= new JButton("Auto");
      btAutoSplit.addActionListener(this);
      btGroup= new JButton("Merge Info");
      btGroup.addActionListener(this);
      btAutoGroup= new JButton("Auto");
      btAutoGroup.addActionListener(this);
      btShowAvail = new JButton("Show Available Frames");
      btShowAvail.addActionListener(this);
      btAutoDisc = new JButton("Auto Discover");
      btAutoDisc.addActionListener(this);
      
      // right panel
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(220));
      rightp.setBorder(BorderFactory.createEmptyBorder(8, 4, 8, 4));
      this.add(new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER), BorderLayout.EAST);
      rightp.add(btDisc);
      p = new JPanel();
      p.add(btSave);
      btSave.setEnabled(false);
      p.add(btLoad);
      rightp.add(p);
      rightp.add(btTimeMerge);
      btTimeMerge.setEnabled(false);
      rightp.add(Box.createVerticalStrut(4));
      rightp.add(btTimeGrow);
      btTimeGrow.setEnabled(false);
      rightp.add(Box.createVerticalStrut(4));
      p = new JPanel();
      p.add(btSplit);
      p.add(btAutoSplit);
      rightp.add(p);
      btSplit.setEnabled(false);
      btAutoSplit.setEnabled(false);
      rightp.add(Box.createVerticalStrut(4));
      p = new JPanel();      
      p.add(btGroup);
      p.add(btAutoGroup);
      rightp.add(p);
      btGroup.setEnabled(false);
      btAutoGroup.setEnabled(false);
      rightp.add(Box.createVerticalStrut(4));
      rightp.add(btShowAvail);      
      rightp.add(Box.createVerticalStrut(4));
      rightp.add(btAutoDisc);
      
      // left panel
      VerticalScrollPanel leftp = new VerticalScrollPanel(new VerticalLayout(220));
      leftp.setBorder(BorderFactory.createEmptyBorder(8, 4, 8, 4));
      this.add(new JScrollPane(leftp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER), BorderLayout.WEST);
      
      leftp.add(Box.createVerticalStrut(4));
      tfMDL = new JTextField(""+mdlThresh, 12);
      p = new JPanel();
      p.add(new JLabel("MDL: "));
      p.add(tfMDL);
      leftp.add(p);
      leftp.add(Box.createVerticalStrut(4));
      tfGrow = new JTextField(""+growThresh, 12);
      p = new JPanel();
      p.add(new JLabel("Grow: "));
      p.add(tfGrow);
      leftp.add(p);
      leftp.add(Box.createVerticalStrut(4));
      tfSplit = new JTextField(""+splitThresh, 12);
      p = new JPanel();
      p.add(new JLabel("Split: "));
      p.add(tfSplit);
      leftp.add(p);
      leftp.add(Box.createVerticalStrut(4));
      cbPreset = new JComboBox(sPresets);
      cbPreset.addActionListener(this);
      leftp.add(cbPreset);
      
      // center panel setup
      pCenter = new JPanel(new BorderLayout());
      pEntries = new VerticalScrollPanel(10, 50, new VerticalLayout(-1, 2));
      pEntries.setBackground(Color.black);
      this.add(pCenter, BorderLayout.CENTER);     

      // setup file choose
      fc.addChoosableFileFilter(new FileFilter() {
         public boolean accept(File f)
         {
            if (f.isDirectory()) return true;
            String ext = Library.getExt(f.getName());
            if (ext == null) return false;
            return (ext.equals("seeds"));
         }

         public String getDescription()
         {
            return "Discovered Seed Locations (*.seeds)";
         }
      });

   }

   protected void removeAvail(HashViewEntry hve)
   {
      for(WindowLocation wloc : hve.getOccs())
      {
         assert(spanAvail[wloc.iSeries].contains(wloc.getRange()));
         spanAvail[wloc.iSeries].sub(wloc.getRange());
      }
   }
   
   protected void addAvail(HashViewEntry hve)
   {
      for(WindowLocation wloc : hve.getOccs())
         spanAvail[wloc.iSeries].add(wloc.getRange());
   }
   
   protected void discover()
   {
      mdlThresh = Double.parseDouble(tfMDL.getText());
      splitThresh = Double.parseDouble(tfSplit.getText());
      growThresh = Double.parseDouble(tfGrow.getText());
      
      btDisc.setEnabled(false);
      final DiscoverView discoverView = this;
      threadDisc = new Thread()
      {
         public void run()
         {
            resetAvail();
            pCenter.removeAll();
            pCenter.add(new JScrollPane(pEntries), BorderLayout.CENTER);
            pEntries.removeAll();
            prevHVE = null;
            hashView.reset();
            motifs = new ArrayList<ArrayList<WindowLocation>>();
            while(true)
            {
               hashView.runHash();         
               if (hashView.pEntries.getComponentCount() == 0) break;
               HashViewEntry hve = new HashViewEntry((HashViewEntry)hashView.pEntries.getComponent(0));
               hve.addMouseListener(discoverView);         
               System.err.printf("hve mdl=%.2f  #occs=%d\n", hve.getMDL(), hve.getCountDTW());
               if (hve.getMDL() > mdlThresh){
                  System.err.printf("MDL too high (%.2f vs. %.2f)\n", hve.getMDL(), mdlThresh);
                  break;
               }
               removeAvail(hve);
               motifs.add(hve.getOccs());
               hashView.removeMotif(hve.getOccs());               
               pEntries.add(hve);
            }            
            btDisc.setEnabled(true);
         }
      };
      threadDisc.start();
   }
   
   protected void groupInfo()
   {
      if (prevHVE==null)
      {
         System.err.println("Please select a motif.");
         return;
      }
      
      SeqDist seqd = SeqDistFactory.createNormSSD();
      int nMotifs = motifs.size();
      assert(nMotifs == pEntries.getComponentCount());
         
      ArrayList<Sequence> examples = WindowLocation.getExamples(prevHVE.getOccs(), gdata.tseries);
      
      for(int iMotif = 0; iMotif<nMotifs; iMotif++)
      {
         ArrayList<Sequence> occs = new ArrayList<Sequence>();
         for(WindowLocation wloc : motifs.get(iMotif))
         {
            Sequence seq = gdata.tseries.get(wloc.iSeries);
            Sequence sub = seq.subseq(wloc.iStart, wloc.getLastIndex() + 1, wloc.iSeries);
            occs.add(sub);
         }
                  
         double dmin = Library.INF;
         double dmax = Library.NEGINF;
         double dsum = 0;         
         
         for(int i=0; i<occs.size(); i++)
            for(int j=i+1; j<occs.size(); j++)
            {
               double d = seqd.dist(occs.get(i), occs.get(j));
               if (d<dmin) dmin = d;
               if (d>dmax) dmax = d;
               dsum += d;
            }
         int nsum = (occs.size() * (occs.size()-1))/2;
         double davg = dsum / nsum;
         
         System.err.printf("iMotif=%d  min=%.3f    max=%.4f   avg=%.3f\n", iMotif, dmin, dmax, davg);
      }
   }
   
   protected boolean canMerge(double dMerge, double dChild)
   {
      final double a = 50;
      final double b = 1000;
      double y = Math.sqrt(a*dChild+b);
      return (dMerge<y);
   }
   
   protected void autoGroup()
   {
      final boolean bDendrogram = false; // TODO - check box in UI
      final boolean bLoad = false;
      btAutoGroup.setEnabled(false);
      final DiscoverView discoverView = this;
      threadDisc = new Thread()
      {
         public void run()
         {
            Component[] comps = pEntries.getComponents();
            int N = comps.length;
            assert(N == motifs.size());
            SeqDist seqd = SeqDistFactory.createNormSSD();
            
            // create lists of examples and dist array
            ArrayList<ArrayList<Sequence>> examples = new ArrayList<ArrayList<Sequence>>(); 
            double dmap[][] = new double[N][];
            for(int i=0; i<N; i++)
            {
               dmap[i] = new double[i+1];
               ArrayList<Sequence> occs = WindowLocation.getExamples(motifs.get(i), gdata.tseries);
               examples.add(occs);
            }
            
            // fill in the dmap
            if (bLoad)
            {
               try{
                  BufferedReader in = new BufferedReader(new FileReader("dmap.txt"));
                  int n = Integer.parseInt(in.readLine());
                  assert(n == dmap.length);
                  for(int i=0; i<dmap.length; i++)
                  {
                     String line = in.readLine();
                     StringTokenizer st = new StringTokenizer(line, " \r\n");
                     for(int j=0; j<dmap[i].length; j++) dmap[i][j] = Double.parseDouble(st.nextToken());
                  }
                  in.close();
               } catch(Exception e){}
               pEntries.removeAll();
            }
            else{
               System.err.println("Building distance map...");
               TimerMS timer = new TimerMS();
               for(int iComp=0; iComp<N; iComp++)
               {   
                  System.err.printf(" Motif %d... ", iComp+1);
                  TimerMS timer2 = new TimerMS();               
                  ArrayList<Sequence> occsi = examples.get(iComp);
                  for(int jComp=0; jComp<=iComp; jComp++)
                  {
                     ArrayList<Sequence> occsj = examples.get(jComp);
                     
                     double dmin = Library.INF;
                     double dmax = Library.NEGINF;
                     double dsum = 0;
                     
                     for(int i=0; i<occsi.size(); i++)
                        for(int j=0; j<occsj.size(); j++)
                        {
                           double d = seqd.dist(occsi.get(i), occsj.get(j));
                           if (d<dmin) dmin = d;
                           if (d>dmax) dmax = d;
                           dsum += d;
                        }
                     double davg = dsum / (occsi.size() * occsj.size());
                     dmap[iComp][jComp] = davg; // min, max, or avg?
                  }
                  System.err.printf("done. (%dms)\n", timer2.time());               
                  
                  pEntries.remove(0);
                  revalidate();
               }
               System.err.printf("done. (%dms)\n", timer.time());
            }
            
            /*try{
               PrintWriter out = new PrintWriter(new FileWriter("dmap.txt"));
               out.printf("%d\n", dmap.length);
               for(int i=0; i<dmap.length; i++)
               {
                  for(int j=0; j<dmap[i].length; j++) out.printf("%.6f  ", dmap[i][j]);
                  out.println();
               }
               out.close();
            } catch(Exception e){}*/
            
            AgglomCluster ac = new AgglomCluster();
            ac.cluster(dmap, DoubleCompFactory.createAvg());
            if (bDendrogram)
            {
               pCenter.removeAll();            
               Dendrogram dend = new Dendrogram(ac.getRoot());
               dend.addMouseListener(discoverView);
               pCenter.add(dend);
            }
            else{
               ArrayList<HashViewEntry> hves = new ArrayList<HashViewEntry>();
               for(int i=0; i<comps.length; i++) hves.add((HashViewEntry)comps[i]);
               ArrayList<AgglomInfo> path = ac.getPath();
               for(int iPath=0; iPath<path.size(); iPath++)
               {
                  AgglomInfo ai = path.get(iPath);
                  assert(ai.hasKids());
                  double d = ai.getDist();
                  double d1 = ai.getChild(1).getDist();
                  double d2 = ai.getChild(2).getDist();
                  //System.err.printf("possible merge: %.3f vs %.3f, %.3f\n", d, d1, d2);
                  if (canMerge(d, d1) && canMerge(d, d2))
                  {
                     // merge the two HVEs
                     int i = ai.geti();
                     int j = ai.getj();
                     assert(i>j);
                     
                     //System.err.printf("merged: %d with %d\n", i, j);
                     
                     HashViewEntry hvei = hves.get(i);
                     HashViewEntry hvej = hves.get(j);
                     ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
                     occs.addAll(hvei.getOccs());
                     occs.addAll(hvej.getOccs());
                     HashViewEntry hveNew = new HashViewEntry(hvej.getHash(), 0, occs.size(), 0, 0, 0, colors);
                     hveNew.setOccs(occs);
                     hveNew.addMouseListener(discoverView);
                     hves.remove(i);
                     hves.set(j, hveNew);
                  }
                  else break;
               }
               
               // rebuild motif list and entries panel
               motifs.clear();
               for(int i=0; i<hves.size(); i++)
               {
                  pEntries.add(hves.get(i));
                  motifs.add(hves.get(i).getOccs());
               }
            }
            revalidate();
         }
      };
      threadDisc.start();
   }
   
   protected void save()
   {
      if (fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION)
      {
         File f = fc.getSelectedFile();
         if (!f.getName().endsWith(".seeds"))
            f = new File(f.getParentFile(), f.getName()+".seeds");
         try{
            // TODO use MotifIO
            PrintWriter out = new PrintWriter(new FileWriter(f));
            out.printf("%d\n", motifs.size());
            for(int i=0; i<motifs.size(); i++)
            {
               HashViewEntry hve = (HashViewEntry)pEntries.getComponent(i);
               ArrayList<WindowLocation> occs = motifs.get(i);
               out.printf("%s %d %d %f\n", hve.getHash(), hve.getCountFixed(), occs.size(), hve.getMDL());
               for(int j=0; j<occs.size(); j++)
                  out.println(occs.get(j).toText());
            }
            out.close();
         }
         catch(IOException e){
            System.err.println("Error: failed to save seed location information in file:");
            System.err.println(" " + Library.getCanonical(f));
            JOptionPane.showMessageDialog(this, "Failed to save seed location info", "Save Failed",
                  JOptionPane.ERROR_MESSAGE);
            return;
         }                  
         JOptionPane.showMessageDialog(this, "Saved seed location info", "Save Successful",
               JOptionPane.INFORMATION_MESSAGE);
      }
   }
   
   protected boolean load()
   {
      if (fc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION)
      {         
         File f = fc.getSelectedFile();         
         pCenter.removeAll();
         pCenter.add(new JScrollPane(pEntries), BorderLayout.CENTER);
         pEntries.removeAll();
         prevHVE = null;
         resetAvail();
         ArrayList<Motif> mis = MotifIO.loadSeeds(f, new MotifInfo());
         if (mis == null)
         {
            JOptionPane.showMessageDialog(this, "Failed to load seed motif locations", "Load Failed",
                  JOptionPane.ERROR_MESSAGE);            
         }
         else{
            int nMotifs = mis.size();
            HashViewEntry[] hves = new HashViewEntry[nMotifs];
            motifs = new ArrayList<ArrayList<WindowLocation>>();
            for(int i=0; i<nMotifs; i++)
            {
               MotifInfo mi = (MotifInfo)mis.get(i);
               motifs.add(mi.getOccs());
               hves[i] = new HashViewEntry(mi, i+1, nMotifs, colors);               
               hves[i].setOccs(mi.getOccs());
               removeAvail(hves[i]);
               hves[i].addMouseListener(this);
               pEntries.add(hves[i]);
               revalidate();
            }
         }         
         return true;
      }
      return false;
   }
   
   /**
    * Search for motifs that follow each other and affix them
    */
   protected boolean affixMotifs()
   {
      
      int nMotif = motifs.size();
      int nMaxGap = 2;

      // build the affix map
      int[][] map = new int[nMotif][nMotif];
         
      // calc map entries
      for(int i=0; i<nMotif; i++)
      {
         for(int ii=0; ii<motifs.get(i).size(); ii++)
         {
            WindowLocation iloc = motifs.get(i).get(ii);
            for(int j=0; j<nMotif; j++)
            {
               for(int jj=0; jj<motifs.get(j).size(); jj++)
               {
                  if (i==j && ii==jj) continue;
                  WindowLocation jloc = motifs.get(j).get(jj);
                  if (iloc.iSeries != jloc.iSeries) continue;
                  int td = jloc.getFirstIndex() - iloc.getLastIndex();
                  if (td>0 && td<=nMaxGap) map[i][j]++;
               }
            }
         }
      }
      
      // search for motifs to affix      
      int iBest = -1;
      int jBest = -1;
      double vBest = 0;
      for(int i=0; i<nMotif; i++)
      {
         for(int j=0; j<nMotif; j++)
         {
            if (i==j) continue;
            int nClose = map[i][j];
            int niMotifs = motifs.get(i).size();
            int njMotifs = motifs.get(j).size();
            double percent = (double)nClose / Math.max(niMotifs, njMotifs);
            if (percent > vBest)
            {
               vBest = percent;
               iBest = i;
               jBest = j;
            }
         }
      }
      
      // do we have any merge possibilities?
      if (iBest<0) return false;      
      
      //System.err.printf("best merge candidate: %d -> %d (%d / %d,%d = %.4f)\n", iBest, jBest, map[iBest][jBest], motifs.get(iBest).size(), motifs.get(jBest).size(), vBest);
      
      // is this good enough for a merge?      
      if (vBest > MIN_FOR_AFFIX)
      {
         // remove old motifs
         motifs.remove(iBest);
         motifs.remove(iBest<jBest ? jBest-1 : jBest);
         
         HashViewEntry hve1 = (HashViewEntry)pEntries.getComponent(iBest);
         HashViewEntry hve2 = (HashViewEntry)pEntries.getComponent(jBest);
         addAvail(hve1);
         addAvail(hve2);
         pEntries.remove(iBest);
         pEntries.remove(iBest<jBest ? jBest-1 : jBest);
         
         // create the new motif
         String sHash = hve1.getHash()+hve2.getHash();
         //System.err.printf("sHash=%s\n", sHash);
         HashMap<OccInfo, MyIntList> hashFixed = gdata.sufTree.findWarpedOccs(sHash, 1, 0);
         HashMap<OccInfo, MyIntList> hashDTW = gdata.sufTree.findWarpedOccs(sHash, 3, 1);
         int nCountFixed = hashFixed.size();
         int nCountDTW = hashDTW.size();
         double mdl = 0.0;    
         HashViewEntry hveNew = new HashViewEntry(sHash, nCountFixed, nCountDTW, mdl, 0, 0, colors);
         hveNew.addMouseListener(this);         
         TreeMap<OccInfo, MyIntList> sortedOccs = HashView.buildTreeMap(sHash, hashDTW);
         ArrayList<WindowLocation> occs = HashView.extractOccs(sortedOccs, tseries, spanAvail);
         hveNew.setOccs(occs);
         removeAvail(hveNew);
         pEntries.add(hveNew);
         motifs.add(occs);
         revalidate();
         return true;
      }
      
      return false;
   }
   
   protected HashViewEntry grow(HashViewEntry hve, boolean bFront, SpanList[] spanAvail)
   {
      String sMotif;
      if (bFront) sMotif = "?"+hve.getHash();
      else sMotif = hve.getHash() + "?";
      HashViewEntry ret = new HashViewEntry(sMotif, hve.getCountFixed(), hve.getCountDTW(), hve.getMDL(),0, 0, colors);
      ArrayList<WindowLocation> occs = hve.getOccs();
      for(WindowLocation occ: occs)
      {
         occ.nLength++;
         if (bFront)
         {            
            occ.iStart--;
            assert(spanAvail[occ.iSeries].contains(occ.iStart));
            spanAvail[occ.iSeries].sub(occ.iStart);
         }
         else{
            assert(spanAvail[occ.iSeries].contains(occ.getLastIndex()));
            spanAvail[occ.iSeries].sub(occ.getLastIndex());
         }
      }
      ret.setOccs(occs);      
      ret.addMouseListener(this);      
      return ret;
   }
   
   protected boolean grow()
   {      
      int nMotif = motifs.size();
      assert(pEntries.getComponentCount() == nMotif) : String.format("comps=%d  motifs=%d", pEntries.getComponentCount(), nMotif);
      
      // compute info relevant to growing
      MotifGrowInfo mgi[] = new MotifGrowInfo[nMotif];
      for(int iMotif=0; iMotif<nMotif; iMotif++)
      {         
         HashViewEntry hve = (HashViewEntry)pEntries.getComponent(iMotif);
         mgi[iMotif] = new MotifGrowInfo(iMotif, hve.getOccs(), gdata, spanAvail);         
      }
      
      // find the most promising extension = smallest since we're talking about variances
      int iBest = -1;
      boolean bBestFront = false;
      FeatureVec fvBest = null;
      for(int iMotif=0; iMotif<nMotif; iMotif++)
      {
         if (mgi[iMotif].canGrowFront() && (fvBest==null || MotifGrowInfo.compareVecs(mgi[iMotif].fvSDevFront, fvBest)<0))
         {
            iBest = iMotif;
            fvBest = mgi[iMotif].fvSDevFront;
            bBestFront = true;
         }
         if (mgi[iMotif].canGrowBack() && (fvBest==null || MotifGrowInfo.compareVecs(mgi[iMotif].fvSDevBack, fvBest)<0))
         {
            iBest = iMotif;
            fvBest = mgi[iMotif].fvSDevBack;
            bBestFront = false;
         }
      }
      
      // did we find something to grow?      
      if (iBest >= 0)
      {         
         // TODO is projection the right thing to do?  could just compare each dim individually 
         double lenBest = fvBest.projLen(FeatureVec.ones(fvBest.getNumDims()));
         //System.err.printf("best to grow: imotif=%d  front=%b  vec=%s (len=%.3f)\n", iBest, bBestFront, fvBest, lenBest);
         if (lenBest < growThresh) // is it good enough to grow?
         {
            HashViewEntry hve = (HashViewEntry)pEntries.getComponent(iBest);
            pEntries.remove(iBest);
            hve = grow(hve, bBestFront, spanAvail);
            pEntries.add(hve, iBest);            
            motifs.set(iBest, hve.getOccs());
            revalidate();
            return true;
         }
         //else System.err.printf("Not good enough to grow (%.3f)\n", lenBest);
      }
      //else System.err.printf("nothing to grow\n");
      return false;
   }
   
   protected AgglomSeqCluster cluster(HashViewEntry hve)
   {
      ArrayList<Sequence> examples = WindowLocation.getExamples(hve.getOccs(), gdata.tseries);      
      asc = new AgglomSeqCluster(SeqDistFactory.createNormSSD());
      asc.cluster(examples, DoubleCompFactory.createMax());
      return asc;
   }
   
   protected void autoSplit()
   {
      btAutoSplit.setEnabled(false);
      final DiscoverView discoverView = this;
      threadDisc = new Thread()
      {
         public void run()
         {
            Component[] comps = pEntries.getComponents();
            pEntries.removeAll();
            motifs.clear();
            revalidate();
            for(int iComp=0; iComp<comps.length; iComp++)
            {
               HashViewEntry hve = (HashViewEntry)comps[iComp];
               AgglomSeqCluster asc = cluster(hve);
               AgglomInfo aiRoot = asc.getRoot();
               double d = aiRoot.getDist();
               double d1 = aiRoot.getChild(1).getDist();
               double d2 = aiRoot.getChild(2).getDist();
               double gap = 1.0 - Math.max(d1, d2) / d;
               System.err.printf("iComp=%d  d=%.3f  gap=%.3f%%\n", iComp, d, gap*100.0);
               if (gap > splitThresh)
               {
                  ArrayList<WindowLocation> occsAll = hve.getOccs();
                  
                  AgglomInfo ai = aiRoot.getChild(1);
                  ArrayIntList list = ai.collectMembers();                  
                  ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
                  for(int i=0; i<list.size(); i++) occs.add(occsAll.get(list.get(i)));
                  HashViewEntry hve1 = new HashViewEntry(hve.getHash(), 0, occs.size(), 0.0, 0, 0, colors);
                  hve1.setOccs(occs);
                  hve1.addMouseListener(discoverView);
                  
                  ai = aiRoot.getChild(2);
                  list = ai.collectMembers();
                  occs = new ArrayList<WindowLocation>();
                  for(int i=0; i<list.size(); i++) occs.add(occsAll.get(list.get(i)));
                  HashViewEntry hve2 = new HashViewEntry(hve.getHash(), 0, occs.size(), 0.0, 0, 0, colors);
                  hve2.setOccs(occs);
                  hve2.addMouseListener(discoverView);
                  
                  pEntries.add(hve1);
                  motifs.add(hve1.getOccs());
                  pEntries.add(hve2);        
                  motifs.add(hve2.getOccs());
               }
               else{
                  pEntries.add(hve);
                  motifs.add(hve.getOccs());
               }
               revalidate();
            }            
         }
      };
      threadDisc.start();
   }
   
   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();
      
      if (src == btDisc)
      {
         discover();
         btSave.setEnabled(true);
         btTimeGrow.setEnabled(true);
         btTimeMerge.setEnabled(true);
         btSplit.setEnabled(false);
         btAutoSplit.setEnabled(true);
         btGroup.setEnabled(false);
         btAutoGroup.setEnabled(true);
      }
      else if (src == btSave) save();
      else if (src == btLoad)
      {
         if (!load()) return;
         btSave.setEnabled(true);
         btTimeGrow.setEnabled(true);
         btTimeMerge.setEnabled(true);
         btSplit.setEnabled(false);
         btAutoSplit.setEnabled(true);
         btGroup.setEnabled(false);
         btAutoGroup.setEnabled(true);
      }
      else if (src == btTimeMerge)
      {         
         while(affixMotifs()){}
      }
      else if (src == btTimeGrow)
      {         
         while(grow()){}
      }
      else if (src == btSplit)
      {
         if (prevHVE == null)
         {
            System.err.println("you must select a motif before clustering!");
         }
         else{
            AgglomSeqCluster asc = cluster(prevHVE);
            if (asc != null)
            {
               pCenter.removeAll();
               Dendrogram dend = new Dendrogram(asc.getRoot());
               dend.addMouseListener(this);
               pCenter.add(dend);            
               revalidate();
            }
         }
      }
      else if (src == btAutoSplit)
      {               
         autoSplit();
      }
      else if (src == btGroup)
      {         
         groupInfo();
      }
      else if (src == btAutoGroup)
      {
         autoGroup();
      }
      else if (src == btShowAvail)
      {
         if (prevHVE != null) prevHVE.setActive(false);
         for(int i = 0; i < qv.length; i++)
         {
            qv[i].clear();
            qv[i].highlight(spanAvail[i]);
         }
      }
      else if (src == btAutoDisc)
      {
         btAutoDisc.setEnabled(false);
         // simulate the various button press
         Thread thread = new Thread()
         {
            public void run()
            {
               ActionEvent e2 = new ActionEvent(btAutoDisc, 0, null);
               e2.setSource(btDisc);
               TimerMS timerA = new TimerMS();
               TimerMS timer = new TimerMS();
               actionPerformed(e2);               
               waitForDiscThread();
               System.err.printf("** seed discovery: %dms\n", timer.time());
               e2.setSource(btAutoSplit);
               timer.reset();
               actionPerformed(e2);
               waitForDiscThread();
               System.err.printf("** motif split: %dms\n", timer.time());
               e2.setSource(btAutoGroup);
               timer.reset();
               actionPerformed(e2);
               waitForDiscThread();
               System.err.printf("** motif merge: %dms\n", timer.time());
               e2.setSource(btTimeMerge);
               timer.reset();
               actionPerformed(e2);
               System.err.printf("** temporal merge: %dms\n", timer.time());
               e2.setSource(btTimeGrow);
               timer.reset();
               actionPerformed(e2);
               System.err.printf("** temporal grow: %dms\n", timer.time());
               System.err.printf("** total time (%d motifs): %dms\n", motifs.size(), timerA.time());               
               e2.setSource(btSave);               
               actionPerformed(e2);   
               
               btAutoDisc.setEnabled(true);
            }
         };
         thread.start();
      }
      else if (src == cbPreset){
         int ix = cbPreset.getSelectedIndex();
         if (ix == iExercise){
            gdata.hashView.lslWLen.setValue(6);
            gdata.hashView.lslMinSymbols.setValue(2);
            gdata.hashView.lslMinTrans.setValue(2);
            gdata.hashView.lslMinCountFixed.setValue(1);
            gdata.hashView.lslMinCountDTW.setValue(30);
            gdata.hashView.cbSortMotifs.setSelectedIndex(HashView.SORTBY_MDL);
            gdata.hashView.lslMaxStretch.setValue(3);
            gdata.hashView.lslMaxDontCare.setValue(1);
         }
         else if (ix == iSpeech){
            gdata.hashView.lslWLen.setValue(20);
            gdata.hashView.lslMinSymbols.setValue(2);
            gdata.hashView.lslMinTrans.setValue(2);
            gdata.hashView.lslMinCountFixed.setValue(1);
            gdata.hashView.lslMinCountDTW.setValue(4);
            gdata.hashView.cbSortMotifs.setSelectedIndex(HashView.SORTBY_MDL);
            gdata.hashView.lslMaxStretch.setValue(2);
            gdata.hashView.lslMaxDontCare.setValue(0);
         }
      }
   }
   
   protected void waitForDiscThread()
   {
      synchronized(threadDisc)
      {
         try{
            threadDisc.wait();
         } catch (InterruptedException e){
            e.printStackTrace();
         }
      }
   }

   public void mouseClicked(MouseEvent e)
   {
   }
   
   public void mousePressed(MouseEvent e)
   {
      Object src = e.getSource();
      
      if (src instanceof HashViewEntry)
      {
         HashViewEntry hve = (HashViewEntry)src;
         select(hve.getOccs(), true);
         if (prevHVE != null) prevHVE.setActive(false);
         hve.setActive(true);         
         prevHVE = hve;         
         btSplit.setEnabled(true);
         btGroup.setEnabled(true);
      }
      else if (src instanceof Dendrogram)
      {
         Dendrogram dend = (Dendrogram)src;
         AgglomInfo ai = dend.getMouseAI();
         ArrayIntList list = ai.collectMembers();
         if (prevHVE == null)
         {
            System.err.print("Members: ");
            for(int i=0; i<list.size(); i++) System.err.printf("(%d:%d) ", i, list.get(i));
            System.err.println();
            select(motifs.get(list.get(0)), true);
            for(int i=1; i<list.size(); i++) select(motifs.get(list.get(i)), false);
         }
         else{
            ArrayList<WindowLocation> occsAll = prevHVE.getOccs();
            ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();
            for(int i=0; i<list.size(); i++) occs.add(occsAll.get(list.get(i)));
            select(occs, true);
         }
      }
   }

   public void mouseReleased(MouseEvent e)
   {
   }

   public void mouseEntered(MouseEvent e)
   {
   }

   public void mouseExited(MouseEvent e)
   {
   }
}
