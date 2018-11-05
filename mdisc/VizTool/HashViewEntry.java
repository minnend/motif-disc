package mdisc.VizTool;

import java.awt.*;
import java.util.ArrayList;

import javax.swing.*;
import kdm.data.*;
import kdm.gui.*;

/**
 * Visual representation of a hash entry computed when searching for matching subsequences
 */
public class HashViewEntry extends JPanel
{
   protected MotifInfo motif;
   protected int iRank, nEntries;
   protected Color colorDef;
   protected Color[] colorq;
   protected static Color colorActive = new Color(120, 200, 230);
   
   public HashViewEntry(HashViewEntry hve)
   {
      this(hve.getHash(), hve.getCountFixed(), hve.getCountDTW(), hve.getMDL(), hve.iRank, hve.nEntries, hve.colorq);
      setOccs(hve.getOccs());
   }
   
   public HashViewEntry(MotifInfo mi, int _iRank, int _nEntries, Color[] _colorq)
   {
      super(new GridFlexLayout(1, 6, 12, 0));
      GridFlexLayout gfl = (GridFlexLayout)getLayout();
      gfl.setRow(0, GridFlexLayout.Style.pref);
      gfl.setColumn(0, GridFlexLayout.Style.fixed, 90); // hash string
      gfl.setColumn(1, GridFlexLayout.Style.fixed, 100); // quant view
      gfl.setColumn(2, GridFlexLayout.Style.fixed, 15);  // empty
      gfl.setColumn(3, GridFlexLayout.Style.fixed, 120); // count
      gfl.setColumn(4, GridFlexLayout.Style.fixed, 150); // mdl
      gfl.setColumn(5, GridFlexLayout.Style.fixed, 150); // rank
      
      motif = mi;
      iRank = _iRank;
      nEntries = _nEntries;
      colorq = _colorq;
      
      setOpaque(true);
      colorDef = getBackground();
      
      setBorder(BorderFactory.createCompoundBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.gray),
            BorderFactory.createEmptyBorder(8,16,8,4)));
      
      DiscreteSeq dseq = new DiscreteSeq(motif.getHash(), getQData(), colorq.length);
      QuantDataView qv = new QuantDataView(dseq, colorq, false, 1);
      qv.setBorder(BorderFactory.createMatteBorder(1,1,1,1,Color.darkGray));
      
      add(new JLabel(String.format("| %s |", motif.getHash())));
      add(qv);
      add(new EmptyComponent());
      add(new JLabel(String.format("Count: %d (%d)", motif.getCountFixed(), motif.getCountDTW())));
      add(new JLabel(String.format("MDL: %.4f", motif.getMDL())));
      if (nEntries > 0) add(new JLabel(String.format("Rank: %d / %d", iRank, nEntries)));
      else add(new EmptyComponent());
   }
   
   public HashViewEntry(String sHash, int nCountFixed, int nCountDTW, double mdl, int _iRank, int _nEntries, Color[] _colorq)
   {
      this(new MotifInfo(sHash, nCountFixed, nCountDTW, mdl, null), _iRank, _nEntries, _colorq);
   }
   
   public void setOccs(ArrayList<WindowLocation> occs){ motif.setOccs(occs); }
   public ArrayList<WindowLocation> getOccs(){ return motif.getOccs(); }
   public int getCountFixed(){ return motif.getCountFixed(); }
   public int getCountDTW(){ return motif.getCountDTW(); }
   public double getMDL(){ return motif.getMDL(); }      
   public String getHash(){ return motif.getHash(); }
   public int getLength(){ return motif.getLength(); }
   
   public int[] getQData()
   {
      int n = getLength();
      String sHash = motif.getHash();
      int[] qdata = new int[n];
      for(int i=0; i<n; i++) qdata[i] = (int)sHash.charAt(i) - 'a';
      return qdata;
   }
   
   public void setActive(boolean bActive)
   {
      if (bActive) setBackground(colorActive);
      else setBackground(colorDef);
      repaint();
   }
   
   public String toString(){ return motif.toString(); }
}
