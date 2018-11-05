package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.util.*;
import kdm.mlpr.*;
import kdm.mlpr.suffix_tree.*;

import org.apache.commons.collections.primitives.*;


/**
 * Displays information about the data based around hashing subsequences with certain
 * properties.
 */
public class HashView extends AbstractDiscView implements ActionListener, ChangeListener, MouseListener, ItemListener,
      ListSelectionListener
{
   public static final int SORTBY_NUM_OCCS = 0;
   public static final int SORTBY_NUM_TRANS = 1;
   public static final int SORTBY_MDL = 2;
   public static final int SORTBY_MDL_RLE = 3;
   public static final int SORTBY_STRING = 4;
   
   protected SuffixTree sufTree;
   protected ArrayList<WindowLocation> occs;
   protected JPanel pCenter;
   protected LabeledSlider lslWLen, lslMinSymbols, lslMinTrans, lslMinCountFixed;
   protected LabeledSlider lslMinCountDTW, lslMaxStretch, lslMaxDontCare;
   protected JButton btHash, btGrow, btGrab, btShowAvail, btReset;
   protected JButton btGrowMotif, btExitGrow;
   protected JList listMotifs;
   protected JComboBox cbSortMotifs;
   protected HashViewEntry prevHVE = null;
   protected HashViewEntry motifHVE = null;
   protected HashAffixEntry prevHAE = null;
   protected VerticalScrollPanel pEntries;   
   protected boolean bGrow = false;   
   
   public HashView(GlobalData gdata)
   {
      super(gdata);
      gdata.hashView = this;
      sufTree = gdata.sufTree;

      // right panel
      VerticalScrollPanel rightp = new VerticalScrollPanel(new VerticalLayout(220));
      rightp.setBorder(BorderFactory.createEmptyBorder(8, 4, 8, 4));
      this.add(new JScrollPane(rightp, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER), BorderLayout.EAST);

      int wlen = 6;
      int nMinSymbols = 2;
      int nMinTrans = 2;
      int nMinCountFixed = 1;
      int nMinCountDTW = 30;
      int nMaxStretch = 3;
      int nMaxDontCare = 1;
      int iSortBy = SORTBY_MDL;

      lslWLen = new LabeledSlider("Subsequence Length: %d", 2, 32);
      lslWLen.getSlider().setPaintTicks(true);
      lslWLen.getSlider().setMinorTickSpacing(2);
      lslWLen.getSlider().setMajorTickSpacing(4);
      lslWLen.addChangeListener(this);
      
      lslMinSymbols = new LabeledSlider("Min # Symbols: %d", 1, wlen);
      lslMinSymbols.getSlider().setPaintTicks(true);
      lslMinSymbols.getSlider().setMinorTickSpacing(1);
      lslMinSymbols.getSlider().setMajorTickSpacing(2);
      
      lslMinTrans = new LabeledSlider("Min # Transitions: %d", 0, wlen-1);
      lslMinTrans.getSlider().setPaintTicks(true);
      lslMinTrans.getSlider().setMinorTickSpacing(1);
      lslMinTrans.getSlider().setMajorTickSpacing(2);      
      
      lslMinCountFixed = new LabeledSlider("Min # Occs (Fixed): %d", 1, 50);
      lslMinCountFixed.getSlider().setPaintTicks(true);
      lslMinCountFixed.getSlider().setMinorTickSpacing(5);
      lslMinCountFixed.getSlider().setMajorTickSpacing(10);      
      
      lslMinCountDTW = new LabeledSlider("Min # Occs (DTW): %d", 1, 50);
      lslMinCountDTW.getSlider().setPaintTicks(true);
      lslMinCountDTW.getSlider().setMinorTickSpacing(5);
      lslMinCountDTW.getSlider().setMajorTickSpacing(10);
      
      lslMaxStretch = new LabeledSlider("Max Stretch: %d", 1, 5);
      lslMaxStretch.getSlider().setPaintTicks(true);
      lslMaxStretch.getSlider().setMinorTickSpacing(1);
      lslMaxStretch.getSlider().setMajorTickSpacing(0);      
      
      lslMaxDontCare = new LabeledSlider("Max Don't Care: %d", 0, wlen-1);
      lslMaxDontCare.getSlider().setPaintTicks(true);
      lslMaxDontCare.getSlider().setMinorTickSpacing(1);
      lslMaxDontCare.getSlider().setMajorTickSpacing(0);
      
      cbSortMotifs = new JComboBox();
      cbSortMotifs.addItem(" # Occurrences");
      cbSortMotifs.addItem(" # Transitions");
      cbSortMotifs.addItem(" MDL");
      cbSortMotifs.addItem(" MDL (RLE)");
      cbSortMotifs.addActionListener(this);      

      btHash = new JButton("Go!");
      btHash.addActionListener(this);

      // set value of sliders
      lslWLen.setValue(wlen);
      lslMinSymbols.setValue(nMinSymbols);
      lslMinTrans.setValue(nMinTrans);
      lslMinCountFixed.setValue(nMinCountFixed);
      lslMinCountDTW.setValue(nMinCountDTW);
      lslMaxStretch.setValue(nMaxStretch);
      lslMaxDontCare.setValue(nMaxDontCare);
      cbSortMotifs.setSelectedIndex(iSortBy);

      btGrab = new JButton("Grab Motif");
      btGrab.setEnabled(false);
      btGrab.addActionListener(this);
      
      btGrow = new JButton("Grow Motif");
      btGrow.setEnabled(false);
      btGrow.addActionListener(this);

      listMotifs = new JList(new DefaultListModel());
      listMotifs.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
      listMotifs.setLayoutOrientation(JList.VERTICAL);
      listMotifs.setVisibleRowCount(6);
      listMotifs.addListSelectionListener(this);

      btShowAvail = new JButton("Show Available Frames");
      btShowAvail.addActionListener(this);
      btReset = new JButton("Reset");
      btReset.addActionListener(this);

      JPanel p;

      Box box = Box.createVerticalBox();
      rightp.add(box, BorderLayout.EAST);
      
      box.add(lslWLen);
      box.add(Box.createVerticalStrut(8));
      box.add(lslMinSymbols);
      box.add(Box.createVerticalStrut(8));
      box.add(lslMinTrans);
      box.add(Box.createVerticalStrut(8));
      box.add(lslMinCountFixed);
      box.add(Box.createVerticalStrut(8));
      box.add(lslMinCountDTW);
      box.add(Box.createVerticalStrut(8));
      p = new JPanel();
      p.add(new JLabel("Sort by: "));
      p.add(cbSortMotifs);
      box.add(p);

      box.add(Box.createVerticalStrut(8));
      box.add(lslMaxStretch);
      box.add(Box.createVerticalStrut(8));
      box.add(lslMaxDontCare);
      
      rightp.add(Box.createVerticalStrut(16));
      rightp.add(btHash);

      rightp.add(Box.createVerticalStrut(32));
      rightp.add(btGrab);
      rightp.add(Box.createVerticalStrut(4));
      rightp.add(new JScrollPane(listMotifs));
      rightp.add(Box.createVerticalStrut(4));
      rightp.add(btGrow);

      rightp.add(Box.createVerticalStrut(32));
      rightp.add(btShowAvail);
      rightp.add(Box.createVerticalStrut(8));
      rightp.add(btReset);

      // center panel setup
      pCenter = new JPanel(new BorderLayout());
      pEntries = new VerticalScrollPanel(10, 50, new VerticalLayout(-1, 2));
      pEntries.setBackground(Color.black);
      pCenter.add(new JScrollPane(pEntries), BorderLayout.CENTER);
      this.add(pCenter, BorderLayout.CENTER);      
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();
      
      if (src == lslWLen)
      {
         int wlen = lslWLen.getValue();
         lslMinSymbols.getSlider().setMaximum(Math.min(wlen, nSymbols));
         lslMinTrans.getSlider().setMaximum(Math.min(wlen-1, nSymbols));
         lslMaxDontCare.getSlider().setMaximum(wlen-1);
      }
   }

   /**
    * @return number of unique entries in the given array
    */
   protected int countUnique(int[] a)
   {
      HashSet<Integer> hash = new HashSet<Integer>();
      for(int i = 0; i < a.length; i++)
         hash.add(new Integer(a[i]));
      return hash.size();
   }
   
   /**
    * @return number of adjacent entries that differ
    */
   protected int countTrans(int[] a)
   {
      int n = 0;
      for(int i = 0; i < a.length-1; i++)
         if (a[i]!=a[i+1]) n++;
      return n;
   }

   public void runHash()
   {
      TimerMS timer = new TimerMS();
      int wlen = lslWLen.getValue();
      int nMinSymbols = lslMinSymbols.getValue();
      int nMinTrans= lslMinTrans.getValue();
      int nMinCountFixed = lslMinCountFixed.getValue();
      int nMinCountDTW = lslMinCountDTW.getValue();
      int nMaxStretch = lslMaxStretch.getValue();
      int iSortBy = cbSortMotifs.getSelectedIndex();

      // TODO: faster to rebuild suffix tree after motif is removed?       
      
      HashMap<StretchyString, OccHashEntry> hash = new HashMap<StretchyString, OccHashEntry>();
      for(int iSeq = 0; iSeq < spanAvail.length; iSeq++)
      {
         DiscreteSeq dseq = qseries.get(iSeq);
         spanAvail[iSeq].itReset();
         while(spanAvail[iSeq].itMore())
         {
            int iStart = spanAvail[iSeq].itNext();
            Range r = new Range(iStart, iStart + wlen - 1);
            if (!spanAvail[iSeq].getItRange().contains(r)) continue;
            int[] qdata = dseq.extract(iStart, wlen);
            int nUniq = countUnique(qdata);
            if (nUniq < nMinSymbols) continue;
            int nTrans = countTrans(qdata);
            if (nTrans < nMinTrans) continue;            
            String s = Library.intArr2Str(qdata);
            StretchyString ss = new StretchyString(s, nMaxStretch);
            OccHashEntry ohe = hash.get(ss);
            if (ohe == null) hash.put(ss, new OccHashEntry(s, 1, nUniq, nTrans));
            else ohe.nOccsFixed++;
         }
      }
      
      // remove entries with too few (fixed) members
      StretchyString[] keys = hash.keySet().toArray(new StretchyString[0]);
      for(int i = 0; i < keys.length; i++)
         if (hash.get(keys[i]).nOccsFixed < nMinCountFixed) hash.remove(keys[i]);
      
      // remove entries with too few (DTW) members
      keys = hash.keySet().toArray(new StretchyString[0]);
      //System.err.printf("# hash after fixed count pruning: %d\n", keys.length);
      for(int i = 0; i < keys.length; i++)
      {
         OccHashEntry oce = hash.get(keys[i]);
         findOccs(oce.sMotif, sufTree);
         oce.occs = occs;
         if (oce.getNumOccsDTW() < nMinCountDTW) hash.remove(keys[i]);
      }
      
      // calculate MDL for each remaining entry
      int nFramesAvail = 0;
      for(int i=0; i<spanAvail.length; i++) nFramesAvail += spanAvail[i].size();
      double mdlBase = Math.log(nFramesAvail) + nFramesAvail * Math.log(nSymbols);
      keys = hash.keySet().toArray(new StretchyString[0]);
      for(int i = 0; i < keys.length; i++)
      {
         OccHashEntry oce = hash.get(keys[i]);
         int nFramesUsed = oce.getNumFramesUsed();
         int nFramesLeft = nFramesAvail - nFramesUsed;
         int nOccs = oce.getNumOccsDTW();
         //  MDL = (bits per symbol) x (num frames) + nOccs x (bits per symbol in full string + 1 new one)
         if (iSortBy == SORTBY_MDL_RLE)
         {
            oce.mdl = 0;
            for(WindowLocation wloc : oce.occs)            
            {
               DiscreteSeq dseq = qseries.get(wloc.iSeries);
               // TODO: could be more efficient (no conv to string, bring log() out of inner loop)
               String sOcc = Library.intArr2Str(dseq.extract(wloc.iStart, wloc.nLength));
               int[] runs = StretchyString.getRuns(sOcc);
               for(int iRun=0; iRun<runs.length; iRun++)
                  oce.mdl -= Math.log(oce.nUniq) + Math.log(runs[iRun]);
            }
         }
         // TODO !! negative sign??
         else oce.mdl = -Math.log(oce.nUniq)*nFramesUsed + nOccs*Math.log(nSymbols+1);
         //else oce.mdl = Math.log(oce.nUniq)*nFramesUsed + nOccs*Math.log(nSymbols+1);
      }      

      // sort the remaining entries according to the current settings      
      OccHashEntry[] values = hash.values().toArray(new OccHashEntry[0]);
      Arrays.sort(values, OccHashEntry.getComparator(iSortBy));
      
      // display the results
      pEntries.removeAll();
      for(int i = 0; i < values.length; i++)
      {
         HashViewEntry hve = new HashViewEntry(values[i].sMotif, values[i].nOccsFixed, values[i].getNumOccsDTW(), 
               values[i].mdl, i + 1, values.length, colors);
         hve.setOccs(values[i].occs);
         hve.addMouseListener(this);
         pEntries.add(hve);
      }
      //System.err.printf("time to hash: %dms\n", timer.time());
      revalidate();
   }

   /**
    * Determine the time series and local offset from a global offset and window length
    * @param ix global offset of position
    * @param wlen window length
    * @return local window location
    */
   public static WindowLocation getLocation(ArrayList<Sequence> tseries, int ix, int wlen)
   {
      int iBase = 0;
      int iSeries = -1, iStart = -1;
      for(int i=0; i<tseries.size(); i++)
      {
         int len = tseries.get(i).length()+1; // add 1 to account for fake terminal character         
         if (ix < iBase+len)
         {
            iSeries = i;
            iStart = ix - iBase;
            break;
         }
         iBase += len;
      }     
      if (iSeries < 0) return null;
      else return new WindowLocation(iSeries, iStart, wlen);
   }
   
   /**
    * Find occurrences of the given subsequence in the available data and remove it from
    * the corresponding span list.
    */
   protected void removeMotif(ArrayList<WindowLocation> occs)
   {
      for(WindowLocation wloc : occs)
      {
         assert(spanAvail[wloc.iSeries].contains(wloc.getRange()));
         spanAvail[wloc.iSeries].sub(wloc.getRange());
      }
      
      pEntries.removeAll();
      revalidate();
      prevHVE = null;
   }

   /**
    * Highlight all occurrences of the given subsequence in the quantized view
    */
   protected void select(int[] seq)
   {
      String sQuery = Library.intArr2Str(seq);      
      findOccs(sQuery, sufTree);
      select(occs, true);
   }
   
   /**
    * Access method to start the process of finding matches to a given query
    * @param sQuery query string to search for
    * @param tree suffix tree in which to search
    */
   protected void findOccs(String sQuery, SuffixTree tree)
   {               
      int nMaxStretch = lslMaxStretch.getValue();
      int nMaxDontCare = lslMaxDontCare.getValue();

      HashMap<OccInfo, MyIntList> hash = tree.findWarpedOccs(sQuery, nMaxStretch, nMaxDontCare);
      
      // build a sorted tree to explore the "best" occurrences first
      TreeMap<OccInfo,MyIntList> sortedOccs = buildTreeMap(sQuery, hash);      
      occs = extractOccs(sortedOccs, tseries, spanAvail);
   }
   
   /**
    * Convert the given hashmap into a treemap
    */
   public static TreeMap<OccInfo,MyIntList> buildTreeMap(String sQuery, HashMap<OccInfo, MyIntList> hash)
   {
      // build a sorted tree to explore the "best" occurrences first
      String sAverage = OccInfo.calcAverageString(sQuery, hash.keySet());      
      TreeMap<OccInfo,MyIntList> sortedOccs = new TreeMap<OccInfo,MyIntList>();           
      Iterator<OccInfo> it = hash.keySet().iterator();
      while(it.hasNext())
      {
         OccInfo occi = it.next();         
         MyIntList coors = hash.get(occi); // must retrieve data before modifying occi         
         occi.score = StretchyString.calcEditDist(sAverage, occi.s);
         sortedOccs.put(occi, coors);
      }
      return sortedOccs;
   }
 
   /**
    * Extract the occurrences from the given tree
    * @param sortedOccs tree of occurrences which may have overlaps
    * @return list of occurrences extracted from tree (no overlaps)
    */
   public static ArrayList<WindowLocation> extractOccs(TreeMap<OccInfo,MyIntList> sortedOccs,
         ArrayList<Sequence> tseries, SpanList[] spanAvail)
   {
      SpanList[] spanCopy = new SpanList[spanAvail.length];
      for(int i=0; i<spanCopy.length; i++) spanCopy[i] = new SpanList(spanAvail[i]);         
      ArrayList<WindowLocation> occs = new ArrayList<WindowLocation>();      
      Iterator<OccInfo> it = sortedOccs.keySet().iterator();
      while(it.hasNext())
      {
         OccInfo occi = it.next();         
         MyIntList coors = sortedOccs.get(occi);
         
         for(int i = 0; i < coors.size(); i++)
         {
            int ixGlobal = coors.get(i);
            WindowLocation wloc = getLocation(tseries, ixGlobal, occi.s.length());
            if (spanCopy[wloc.iSeries].contains(wloc.getRange()))
            {
               occs.add(wloc);
               spanCopy[wloc.iSeries].sub(wloc.getRange());
            }
         }
      }
      
      return occs;
   }
   
   /**
    * Setup the display for the given hash view entry to allow "growing" the motif
    */
   protected void setupGrow(HashViewEntry hve)
   {
      JPanel p;
      int nPrefixTotal, nSuffixTotal;
      int[] nPrefix = new int[nSymbols];
      int[] nSuffix = new int[nSymbols];
      pCenter.removeAll();
      
      motifHVE = hve;
      motifHVE.setActive(true);
      
      // header is just the hve
      pCenter.add(hve, BorderLayout.NORTH);      
                  
      // setup center area
      p = new JPanel(new GridFlexLayout(2,2));
      p.setBorder(BorderFactory.createEmptyBorder(12,0,4,0));
      GridFlexLayout gfl = (GridFlexLayout)p.getLayout();
      gfl.setRow(0, GridFlexLayout.Style.pref);      
      pCenter.add(p, BorderLayout.CENTER);
      
      // setup column headers ("Prefix" and "Suffix")
      JLabel lab = new JLabel("Prefix", JLabel.CENTER);
      Font font = lab.getFont().deriveFont(18.0f);
      lab.setFont(font);      
      p.add(lab);
      lab = new JLabel("Suffix", JLabel.CENTER);
      lab.setFont(font);
      p.add(lab);
      
      // compute prefix/suffix counts
      occs = hve.getOccs();
      System.err.printf("# occs: %d\n", occs.size());
      for(WindowLocation wloc : occs)
      {         
         DiscreteSeq dseq = qseries.get(wloc.iSeries);
         if (spanAvail[wloc.iSeries].contains(wloc.getFirstIndex()-1)) nPrefix[dseq.geti(wloc.getFirstIndex()-1)]++;
         if (spanAvail[wloc.iSeries].contains(wloc.getLastIndex()+1)) nSuffix[dseq.geti(wloc.getLastIndex()+1)]++;
      }
      
      // compute prefix/suffix totals
      nPrefixTotal = nSuffixTotal = 0;
      for(int i=0; i<nSymbols; i++)
      {
         nPrefixTotal += nPrefix[i];
         nSuffixTotal += nSuffix[i];
      }
      
      // setup prefix entries
      JPanel pPrefix = new JPanel(new VerticalLayout(-1, 2));
      pPrefix.setBackground(Color.black);      
      p.add(new JScrollPane(pPrefix, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
      for(int i=0; i<nSymbols; i++)
      {
         if (nPrefix[i] == 0) continue;
         String sPrefix = ""+(char)(i+'a') + hve.getHash();
         HashAffixEntry hvePrefix = new HashAffixEntry(HashAffixEntry.PREFIX, sPrefix, nPrefix[i], nPrefixTotal, colors);
         hvePrefix.addMouseListener(this);
         pPrefix.add(hvePrefix);
      }      
      
      // setup suffix entries
      JPanel pSuffix = new JPanel(new VerticalLayout(-1, 2));
      pSuffix.setBackground(Color.black);
      p.add(new JScrollPane(pSuffix, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
      for(int i=0; i<nSymbols; i++)
      {
         if (nSuffix[i] == 0) continue;
         String sSuffix = hve.getHash() + (char)(i+'a');
         HashAffixEntry hveSuffix = new HashAffixEntry(HashAffixEntry.SUFFIX, sSuffix, nSuffix[i], nSuffixTotal, colors);
         hveSuffix.addMouseListener(this);
         pSuffix.add(hveSuffix);
      }      
      
      // setup bottom panel (exit button)
      p = new JPanel(new FlowLayout(FlowLayout.RIGHT));
      p.setBorder(BorderFactory.createMatteBorder(1,0,0,0,Color.gray));
      if (btExitGrow == null)
      {
         btExitGrow = new JButton("Close Grow View");
         btExitGrow.addActionListener(this);
         
         btGrowMotif = new JButton("Grow Motif");
         btGrowMotif.addActionListener(this);         
      }
      btGrowMotif.setEnabled(false);      
      p.add(btGrowMotif);
      p.add(btExitGrow);
      pCenter.add(p, BorderLayout.SOUTH);
            
      revalidate(); // make sure changes show up
   }
   
   public void reset()
   {
      DefaultListModel dlm = (DefaultListModel)listMotifs.getModel();
      dlm.clear();

      prevHVE = null;
      btGrab.setEnabled(false);

      resetAvail();
      
      pEntries.removeAll();         
      if (bGrow) actionPerformed(new ActionEvent(btExitGrow, 0, null));         
      revalidate();
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();

      if (src == btHash)
      {
         prevHVE = null;
         btGrab.setEnabled(false);
         runHash();         
      }
      else if (src == btGrab)
      {
         btGrab.setEnabled(false);         
         DefaultListModel dlm = (DefaultListModel)listMotifs.getModel();
         prevHVE.setActive(false);
         dlm.addElement(prevHVE);
         removeMotif(prevHVE.getOccs());
         listMotifs.setSelectedIndex(dlm.getSize() - 1);
      }
      else if (src == btGrow)
      {
         HashViewEntry hve = (HashViewEntry)listMotifs.getSelectedValue();
         listMotifs.setEnabled(false);
         btHash.setEnabled(false);
         btGrow.setEnabled(false);
         bGrow = true;
         setupGrow(hve); 
      }
      else if (src == btExitGrow)
      {
         pCenter.removeAll();
         pCenter.add(new JScrollPane(pEntries), BorderLayout.CENTER);
         btHash.setEnabled(true);
         listMotifs.setEnabled(true);
         btGrow.setEnabled(true);
         bGrow = false;
         motifHVE = null;
         revalidate();
      }
      else if (src == btGrowMotif)
      {         
         // TODO: doesn't work!
         HashViewEntry hve = new HashViewEntry(prevHAE.getHash(), prevHAE.getCount(), -1, Double.NaN, 0, 0, colors);
         //hve.setOccs(prevHAE.getO);
         setupGrow(hve); 
      }
      else if (src == btShowAvail)
      {
         listMotifs.clearSelection();
         for(int i = 0; i < qv.length; i++)
         {
            qv[i].clear();
            qv[i].highlight(spanAvail[i]);
         }
      }
      else if (src == btReset)
      {
         reset();
      }      
      else if (src == cbSortMotifs)
      {
         if (pEntries!=null && pEntries.getComponents().length > 0)
            actionPerformed(new ActionEvent(btHash, 1001, btHash.getText()));
      }
   }

   public void mouseClicked(MouseEvent e)
   {}

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
         
         if (bGrow)
         {
            if (prevHAE != null) prevHAE.setActive(false); 
            btGrowMotif.setEnabled(false);
         }
         else{            
            btGrab.setEnabled(true);
            listMotifs.clearSelection();
         }
      }
      else if (src instanceof HashAffixEntry)
      {
         HashAffixEntry hae = (HashAffixEntry)src;         
         select(hae.getQData());
         if (prevHAE != null) prevHAE.setActive(false);
         hae.setActive(true);
         prevHAE = hae;
         motifHVE.setActive(false);
         //btGrowMotif.setEnabled(true); // TODO: uncomment when you fix grow code!
      }
   }

   public void mouseReleased(MouseEvent e)
   {}

   public void mouseEntered(MouseEvent e)
   {}

   public void mouseExited(MouseEvent e)
   {}

   public void itemStateChanged(ItemEvent e)
   {
      ItemSelectable src = e.getItemSelectable();
      boolean bOn = (e.getStateChange() == ItemEvent.SELECTED);

   }

   public void valueChanged(ListSelectionEvent e)
   {
      Object src = e.getSource();
      
      if (src == listMotifs)
      {
         if (listMotifs.isSelectionEmpty())
         {
            btGrow.setEnabled(false);
         }
         else {
            if (prevHVE != null)
            {               
               prevHVE.setActive(false);
               prevHVE = null;
            }
         
            HashViewEntry hve = (HashViewEntry)listMotifs.getSelectedValue();
            occs = hve.getOccs();
            assert(occs!=null && !occs.isEmpty());
            select(hve.getOccs(), true);
            btGrow.setEnabled(true);
            btGrab.setEnabled(false);
         }
      }
      
   }

}
