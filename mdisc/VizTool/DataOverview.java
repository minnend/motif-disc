package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.util.*;
import kdm.models.*;

/**
 * Displays an overview of the data, including mean/variance of each dimension
 */
public class DataOverview extends JPanel implements ActionListener, ChangeListener
{
   protected GlobalData gdata;
   protected JButton btTrain;
   protected JComboBox cbClass;
   protected LabeledSlider lsZoom;
   protected int nDims;
   protected DimensionOverview dimov[];

   public DataOverview(GlobalData gdata)
   {
      super(new BorderLayout());
      setOpaque(true);
      gdata.dataOverview = this;
      this.gdata = gdata;
      
      // setup the top (labels for graphs)
      Box pTop = Box.createHorizontalBox();
      pTop.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));      
      pTop.add(new JLabel("[Distrib: all data]"));
      pTop.add(Box.createHorizontalStrut(4));
      pTop.add(new JLabel("[Gauss per class]"));
      pTop.add(Box.createHorizontalStrut(4));
      pTop.add(new JLabel("[Distrib: var per state]"));
      pTop.add(Box.createHorizontalStrut(4));      
      pTop.add(new JLabel("[Distrib: per state ...  ]"));
      pTop.add(Box.createHorizontalStrut(10));
      cbClass = new JComboBox();
      for(int i=0; i<gdata.getNumClasses(); i++) cbClass.addItem(" "+gdata.classes[i]);
      JPanel p = new JPanel();
      p.add(cbClass);
      pTop.add(p);
      pTop.add(Box.createHorizontalStrut(10));
      btTrain = new JButton("Viz per state info");
      btTrain.setToolTipText("Select model from Word Spotting Tab");
      btTrain.addActionListener(this);
      pTop.add(btTrain);
      pTop.add(Box.createHorizontalStrut(10));
      lsZoom = new LabeledSlider("Zoom: %d", 50, 300);
      lsZoom.setValue(DimensionOverview.DEF_VIEW_SIZE);
      lsZoom.getSlider().setMinorTickSpacing(50);
      lsZoom.getSlider().setMajorTickSpacing(100);
      lsZoom.getSlider().setSnapToTicks(true);
      lsZoom.getSlider().setPaintTicks(true);
      lsZoom.addChangeListener(this);
      pTop.add(lsZoom);
            
      this.add(pTop, BorderLayout.NORTH);
      
      // set up the center (the actual graphs)
      JPanel pDims = new JPanel(new VerticalLayout(-1, 2));
      pDims.setBackground(Color.black);
      this.add(new JScrollPane(pDims), BorderLayout.CENTER);
      
      nDims = gdata.tseries.get(0).getNumDims();
      dimov = new DimensionOverview[nDims];
      for(int d=0; d<nDims; d++)
      {
         dimov[d] = new DimensionOverview(d, gdata);
         pDims.add(dimov[d]);
      }
   }

   protected void updateModelViews()
   {
      gdata.wordSpottingView.train(false); // train models based on wspot params
      int iModel = gdata.wordSpottingView.getCurModelType();
      ProbSeqModel[] models = gdata.wordSpottingView.models[iModel];      
      int iClass = cbClass.getSelectedIndex();
      for(int d=0; d<nDims; d++) dimov[d].calcModelViews(iModel, iClass, models);
   }
   
   public void actionPerformed(ActionEvent e)
   {
      Object src = e.getSource();
      
      if (src == btTrain) updateModelViews();      
   }

   public void stateChanged(ChangeEvent e)
   {
      Object src = e.getSource();
      if (src == lsZoom)
      {
         if (!lsZoom.getSlider().getValueIsAdjusting())
            for(int d=0; d<nDims; d++) dimov[d].setViewSize(lsZoom.getValue());
      }
   }
}
