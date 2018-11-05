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

/**
 * Displays a subsequence tree of the quantized data and shows where in the original data
 * certain leves occur.
 */
public class TreeView extends AbstractDiscView implements ItemListener, ChangeListener, SubseqTreeSelectListener
{
   protected SubseqTreeView tview;
   protected SubseqTree tree;
   protected JCheckBox btVertical;
   protected JLabel lbWLen, lbMinTrans, lbMinCount;
   protected JSlider slWLen, slMinTrans, slMinCount;

   public TreeView(GlobalData gdata)
   {
      super(gdata);
      gdata.treeView = this;

      int wlen = 8;      
      int iMinTrans = 1;
      int iMinCount = 2;
      boolean bVertical = true;

      // center panel (the tree)
      tree = new SubseqTree(qseries, wlen);
      tree.update(iMinTrans, iMinCount);
      tview = new SubseqTreeView(tree);
      tview.setVertical(bVertical);
      tview.addListener(this);
      this.add(tview, BorderLayout.CENTER);

      // right panel      
      JPanel rightp = new JPanel(new VerticalLayout(200));
      rightp.setBorder(BorderFactory.createMatteBorder(0,1,0,0,Color.gray));
      this.add(rightp, BorderLayout.EAST);      
      
      lbWLen = new JLabel("", JLabel.LEFT);
      slWLen = new JSlider(2, 16);
      slWLen.setPaintTicks(true);
      slWLen.setMinorTickSpacing(2);
      slWLen.setMajorTickSpacing(4);
      slWLen.addChangeListener(this);
      
      lbMinTrans = new JLabel("", JLabel.LEFT);
      slMinTrans = new JSlider(0, wlen-1);
      slMinTrans.setPaintTicks(true);
      slMinTrans.setMinorTickSpacing(1);
      slMinTrans.setMajorTickSpacing(2);
      slMinTrans.addChangeListener(this);
      
      lbMinCount = new JLabel("", JLabel.LEFT);
      slMinCount = new JSlider(1, 24);
      slMinCount.setPaintTicks(true);
      slMinCount.setMinorTickSpacing(4);
      slMinCount.setMajorTickSpacing(8);
      slMinCount.addChangeListener(this);
      
      slWLen.setValue(wlen);
      slMinTrans.setValue(iMinTrans);      
      slMinCount.setValue(iMinCount);
      
      btVertical = new JCheckBox("Veritcal Tree", bVertical);
      btVertical.addItemListener(this);
      
      Box box = Box.createVerticalBox();
      box.setBorder(BorderFactory.createEmptyBorder(4,4,2,2));
      rightp.add(box, BorderLayout.EAST);
      box.add(lbWLen);
      box.add(slWLen);
      box.add(Box.createVerticalStrut(8));
      box.add(lbMinTrans);
      box.add(slMinTrans);
      box.add(Box.createVerticalStrut(8));
      box.add(lbMinCount);
      box.add(slMinCount);
      box.add(Box.createVerticalStrut(8));     
      box.add(btVertical);      
   }

   public void itemStateChanged(ItemEvent e)
   {
      ItemSelectable src = e.getItemSelectable();
      boolean bOn = (e.getStateChange() == ItemEvent.SELECTED);

      if (src == btVertical) tview.setVertical(bOn);
   }   

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();      
      
      int iMinTrans = slMinTrans.getValue();
      int iMinCount = slMinCount.getValue();
      
      if (src == slMinTrans || src==slMinCount)
      {         
         lbMinTrans.setText("Min # Transitions: "+iMinTrans);
         lbMinCount.setText("Min # Occurrences: "+iMinCount);
         tree.update(iMinTrans, iMinCount);
         tview.repaint();   
      }
      else if (src == slWLen)
      {
         int wlen = slWLen.getValue();
         lbWLen.setText("Subsequence Length: "+wlen);
         slMinTrans.setMaximum(wlen-1);
         
         tree = new SubseqTree(qseries, wlen);
         tree.update(iMinTrans, iMinCount);
         tview.setTree(tree);
      }
   }
   
   public void select(int[] seq)
   {
      // TODO: display selection string somewhere
      for(int i = 0; i < qv.length; i++)
      {
         qv[i].clear();
         qv[i].highlight(seq);
      }
   }
}
