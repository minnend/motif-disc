package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import javax.swing.*;
import javax.swing.border.*;
import kdm.data.*;
import kdm.gui.*;
import kdm.mlpr.*;
import kdm.util.*;

/**
 * View of all of the raw and quantized data
 */
public class AllDataView extends JPanel implements ActionListener
{
   protected GlobalData gdata;
   protected ArrayList<Sequence> tseries;
   protected ArrayList<DiscreteSeq> qseries;
   protected QuantDataView[] qv;
   protected GraphScrollPane gsp;
   protected JComboBox cbSeries;
   protected JTextField tfSentence;

   public AllDataView(GlobalData gdata)
   {
      super(new BorderLayout());
      gdata.allDataView = this;
      this.gdata = gdata;
      tseries = gdata.tseries;
      qseries = gdata.qseries;

      setOpaque(true);

      if (tseries != null)
      {
         constructLineGraph(tseries.get(0));
         JPanel topp = new JPanel(new BorderLayout());
         JPanel p = new JPanel();
         p.add(new JLabel("View: "));
         cbSeries = new JComboBox();
         for(int i = 0; i < tseries.size(); i++)
            cbSeries.addItem(String.format(" Series %d", i + 1));
         cbSeries.addActionListener(this);
         p.add(cbSeries);
         tfSentence = new JTextField(gdata.getSentence(0), 30);
         tfSentence.setEditable(false);
         tfSentence.setCaretPosition(0);
         p.add(tfSentence);
         topp.add(p, BorderLayout.WEST);
         this.add(topp, BorderLayout.NORTH);
      }

      if (qseries != null)
      {
         Color[] colors = Library.generateColors(gdata.nSymbols);
         JPanel p = new JPanel(new GridFillLayout(qseries.size(), 1));
         p.setOpaque(true);
         p.setBackground(Color.black);
         colors = Library.generateColors(gdata.nSymbols);
         qv = new QuantDataView[qseries.size()];
         for(int i = 0; i < qseries.size(); i++)
         {
            qv[i] = new QuantDataView(qseries.get(i), colors, true, gdata.nAvgSeriesLength>300 ? 0 : 1);
            qv[i].setToolTipText(String.format("Series %d: %s", i+1, gdata.getSentence(i)));
            p.add(qv[i]);
         }
         Constrainer qcons = new Constrainer(new JScrollPane(p), 0, Integer.MAX_VALUE, 0, 260);
         this.add(qcons, BorderLayout.SOUTH);
      }
   }

   /**
    * Constructe the line graph view for the given sequence
    * 
    * @param seq sequence to display
    */
   protected void constructLineGraph(Sequence seq)
   {
      if (gsp != null)
      {
         this.remove(gsp);
         gsp = null;
      }

      GraphComplexLite[] gcl;
      if (gdata.vizgrp != null)
      {        
         int[][] vizgrp = gdata.vizgrp;
         gcl = new GraphComplexLite[vizgrp.length];
         
         for(int ig=0; ig<gcl.length; ig++)
         {
            Color[] dimColor = Library.generateColors(vizgrp[ig].length);
            
            LineGraph lgraph = new LineGraph(seq);
            gcl[ig] = new GraphComplexLite(lgraph);
            lgraph.clearDimColors();
            for(int i = 0; i < vizgrp[ig].length; i++)
               lgraph.setDimColor(vizgrp[ig][i], dimColor[i]);
         }
      }
      else{
         final int nGroup = 3; // group in sets of three by default
         int nGraphs = (int)Math.ceil((double)seq.getNumDims() / nGroup);
         gcl = new GraphComplexLite[nGraphs];
         
         Color[] dimColor = { Color.red, Color.green, Color.blue };
         
         LineGraph lgraph = null;
         for(int i = 0; i < seq.getNumDims(); i++)
         {
            if (i % nGroup == 0)
            {
               lgraph = new LineGraph(seq);            
               gcl[i / nGroup] = new GraphComplexLite(lgraph);
               lgraph.clearDimColors();
            }
            
            lgraph.setDimColor(i, dimColor[i % nGroup]);
         }
      }
      
      // TODO: show labels
      gsp = GraphScrollPane.construct(gcl, gcl[0], false, true, true, null, new ZoomReal(gcl[0]));    
      this.add(gsp, BorderLayout.CENTER);
      revalidate();
   }

   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();
      if (src == cbSeries)
      {
         int i = cbSeries.getSelectedIndex();
         constructLineGraph(tseries.get(i));
         tfSentence.setText(gdata.getSentence(i));
         tfSentence.setCaretPosition(0);
      }
   }
}
