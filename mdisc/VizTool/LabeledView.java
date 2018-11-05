package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.border.*;
import kdm.data.*;
import kdm.gui.*;
import kdm.util.*;

/**
 * Shows the labeled data in both the continuous and quantized domains
 */
public class LabeledView extends JPanel implements MouseListener, ActionListener, ItemListener
{
   protected GlobalData gdata;
   protected JTextField tfInfo;
   protected JCheckBox btShowCont, btShowQuant;
   protected TreeMap<String, ArrayList<LabeledViewEntry>> entries;
   
   public LabeledView(GlobalData gdata)
   {
      super(new BorderLayout());
      gdata.labeledView = this;
      this.gdata = gdata;
      entries = new TreeMap<String, ArrayList<LabeledViewEntry>>();
      
      JPanel p;
      Box box;
      
      int nClasses = gdata.labData.size();
      JPanel labelp = new JPanel(new GridLayout(1, nClasses)); // panel for labels
      JPanel classp = new JPanel(new GridLayout(1, nClasses)); // panel for all classes
      
      p = new JPanel(new BorderLayout()); // p contains labelp and classp
      p.add(labelp, BorderLayout.NORTH);
      p.add(classp, BorderLayout.CENTER);
      JScrollPane sp = new JScrollPane(p);
      this.add(sp, BorderLayout.CENTER);
      
      // fill in the panels
      Set<String> labels = gdata.labData.keySet();
      Iterator<String> it = labels.iterator();
      while(it.hasNext())
      {         
         String label = it.next();
         ArrayList<Sequence> dlist = gdata.labData.get(label);
         ArrayList<DiscreteSeq> qlist = (gdata.labQData != null ? gdata.labQData.get(label) : null);
         ArrayList<LabeledViewEntry> elist = new ArrayList<LabeledViewEntry>();
         entries.put(label, elist);
         labelp.add(new JLabel(label, JLabel.CENTER)); // add the label

         int nSeq = dlist.size();
         //JPanel exp = new JPanel(new GridFillLayout(nSeq, 1, 0, 1)); // example panel (vertical)         
         JPanel exp = new JPanel(new GridFlexLayout(nSeq, 1, 0, 1));
         ((GridFlexLayout)exp.getLayout()).setRows(GridFlexLayout.Style.pref);
         //((GridFlexLayout)exp.getLayout()).setColumns(GridFlexLayout.Style.fixed, 100);
         exp.setBackground(Color.black);
         exp.setOpaque(true);
         exp.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED, Color.gray, Color.darkGray));         
         classp.add(exp);
         for(int iSeq = 0; iSeq < nSeq; iSeq++)
         {
            Sequence seq = dlist.get(iSeq);
            DiscreteSeq qseq = (gdata.labQData != null ? qlist.get(iSeq) : null);
            LabeledViewEntry lve = new LabeledViewEntry(this, seq, qseq, gdata.vizgrp);
            lve.setMaximumSize(new Dimension(-1, 60));
            lve.showLineGraphs(false);
            exp.add(lve);            
            elist.add(lve);
         }
      }

      // setup the info at bottom of screen
      JPanel bottomp = new JPanel(new BorderLayout());
      
      p = new JPanel();
      tfInfo = new JTextField(40);
      tfInfo.setEditable(false);
      p.add(tfInfo);
      bottomp.add(p, BorderLayout.CENTER);

      box = Box.createVerticalBox();
      
      btShowCont = new JCheckBox("Show Continuous", false);
      btShowCont.addItemListener(this);
      box.add(btShowCont);
      btShowQuant = new JCheckBox("Show Quantized", true);
      btShowQuant.addItemListener(this);
      btShowQuant.setEnabled(false);
      box.add(btShowQuant);
      
      if (gdata.labQData!= null)
      {         
         double[] qvals = new double[gdata.nSymbols];
         for(int i = 0; i < gdata.nSymbols; i++)
            qvals[i] = (double)i;
         QuantDataView qview = new QuantDataView(new DiscreteSeq("QV", qvals, gdata.nSymbols), Library
               .generateColors(gdata.nSymbols));
         box.add(qview);
      }

      bottomp.add(box, BorderLayout.EAST);

      this.add(bottomp, BorderLayout.SOUTH);
   }
   
   public void mouseClicked(MouseEvent e)
   {
      Object src = e.getSource();
      if (src instanceof LineGraph)
      {
         LineGraph lgraph = (LineGraph)src;
         Sequence seq = lgraph.getData();
         tfInfo.setText(String.format("Sequence: %d.%d (%d)", seq.getParentIndex(), seq
               .getParentOffset(), seq.length()));

      }
   }

   public void mousePressed(MouseEvent e)
   {}

   public void mouseReleased(MouseEvent e)
   {}

   public void mouseEntered(MouseEvent e)
   {}

   public void mouseExited(MouseEvent e)
   {}

   public void actionPerformed(ActionEvent e)
   {      
   }

   public void itemStateChanged(ItemEvent e)
   {
      ItemSelectable src = e.getItemSelectable();
      boolean bOn = (e.getStateChange() == ItemEvent.SELECTED);
      
      if (src == btShowCont) 
      {         
         Iterator<ArrayList<LabeledViewEntry>> it = entries.values().iterator();
         while(it.hasNext())
         {
            Iterator<LabeledViewEntry> jt = it.next().iterator();
            while(jt.hasNext()) jt.next().showLineGraphs(bOn);
         }
         btShowQuant.setEnabled(bOn);
         revalidate();
      }
      else if (src == btShowQuant)
      {
         Iterator<ArrayList<LabeledViewEntry>> it = entries.values().iterator();
         while(it.hasNext())
         {
            Iterator<LabeledViewEntry> jt = it.next().iterator();
            while(jt.hasNext()) jt.next().showQuant(bOn);
         }
         btShowCont.setEnabled(bOn);
         revalidate();         
      }
   }
}
