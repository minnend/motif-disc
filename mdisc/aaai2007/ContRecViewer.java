package mdisc.aaai2007;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.io.*;

/** combines a series selector and a label view */
public class ContRecViewer extends JPanel implements ActionListener
{
   protected ArrayList<MarkupSet> gt, detected;
   protected JComboBox cbSeq;
   protected LabelView lv;   
   protected LabelEditor led;
   
   public ContRecViewer(ArrayList<MarkupSet> gt, ArrayList<MarkupSet> detected, LabelEditor led)
   {
      super(new BorderLayout());
      
      this.gt = gt;
      this.detected = detected;
      this.led = led;
      
      cbSeq = new JComboBox();
      cbSeq.addActionListener(this);
      for(MarkupSet ms : gt) cbSeq.addItem(ms.getName());
      
      JPanel p = new JPanel(new FlowLayout(FlowLayout.RIGHT, 8, 2));
      p.add(new JLabel("Select Sequence:"));
      p.add(cbSeq);
      add(p, BorderLayout.NORTH);
   }

   public void actionPerformed(ActionEvent e)
   {
      assert(e.getSource() == cbSeq);      
      int ix = cbSeq.getSelectedIndex();
      ArrayList<MarkupSet> labels = new ArrayList<MarkupSet>();
      labels.add(gt.get(ix));
      labels.add(detected.get(ix));  
      if (lv == null){
         lv = new LabelView(labels, led);
         lv.setGroundTruth(0);
         add(lv, BorderLayout.CENTER);
      }
      else lv.setLabels(labels, led);
   }
   
   
}
