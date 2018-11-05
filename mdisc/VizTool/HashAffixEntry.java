package mdisc.VizTool;

import java.awt.*;
import java.util.*;
import javax.swing.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.util.MutableInteger;
import kdm.util.Range;

import static kdm.gui.GridFlexLayout.*;

public class HashAffixEntry extends JPanel
{
   public static final int PREFIX = 0;
   public static final int SUFFIX = 1;
   
   protected int affix;
   protected String sHash;
   protected int nCount;
   protected Color colorDef;
   protected static Color colorActive = new Color(120, 200, 230);
   protected JLabel lbTotal;
   
   public HashAffixEntry(int _affix, String _sHash, int _nCount, int nTotal, Color[] colorq)
   {
      super(new GridFlexLayout(1, 3, 12, 0));
      GridFlexLayout gfl = (GridFlexLayout)getLayout();
      gfl.setRow(0, Style.pref);
      gfl.setColumn(0, Style.fixed, 80);
      gfl.setColumn(1, Style.fixed, 30);
      gfl.setColumn(2, Style.fixed, 100);      
      
      affix = _affix;
      sHash = _sHash;
      nCount = _nCount;
      
      setOpaque(true);
      colorDef = getBackground();
      
      setBorder(BorderFactory.createCompoundBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.gray),
            BorderFactory.createEmptyBorder(8,16,8,4)));
      
      add(new JLabel(String.format("| %s |", sHash)));      
      JPanel p = new JPanel();
      p.setBackground(colorq[getQData()[affix==PREFIX ? 0 : sHash.length()-1]]);
      p.setBorder(BorderFactory.createMatteBorder(1,1,1,1,Color.darkGray));
      add(p);
      if (nCount > 0) add(new JLabel(String.format("%3d = %.1f%%", nCount, 100.0*nCount/nTotal)));
      else add(new EmptyComponent());
   }   
   
   public String getHash(){ return sHash; }
   
   public int getCount(){ return nCount; }
   
   public int[] getQData()
   {
      int n = sHash.length();
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
   
   public String toString(){ return String.format("%s (%d)", sHash, nCount); }
}
