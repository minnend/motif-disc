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
 * A visual entry in a "labeled view" -- basically, two line graphs side by side with a
 * quantum view beneath them.
 */
public class LabeledViewEntry extends JPanel
{
   public static final int WGRAPH = 40;
   public static final int HGRAPH = 40;
   
   protected JPanel lgp;
   protected LineGraph[] lgs;
   protected QuantDataView qv;
   protected int nMaxW = 80;
   protected Dimension maxSize;
   protected int nd;
   protected boolean bShowCont, bShowQuant;

   public LabeledViewEntry(LabeledView parent, Sequence seq, DiscreteSeq qseq, int[][] vizgrp)
   {
      super(new BorderLayout());
      setOpaque(true);
      int nQuant = qseq.getNumSymbols();      
      nd = seq.getNumDims();

      bShowCont = true;
      bShowQuant = true;

      if (vizgrp != null)
      {
         // create the continuous displays
         lgp = new JPanel(new GridFillLayout(1, vizgrp.length, 1, 0));
         lgp.setOpaque(true);
         lgp.setBackground(Color.black);
         this.add(lgp, BorderLayout.CENTER);
                  
         lgs = new LineGraph[vizgrp.length];
         for(int ig=0; ig<vizgrp.length; ig++)
         {            
            lgs[ig] = new LineGraph(seq);
            lgs[ig].setRenderGridLines(false, false);
            lgs[ig].setRenderMouseInfo(false, false, false);
            lgs[ig].addMouseListener(parent);
            lgs[ig].setShowAll(true);
            Color[] dimColor = Library.generateColors(vizgrp[ig].length);
            for(int i=0; i<vizgrp[ig].length; i++)
               lgs[ig].setDimColor(vizgrp[ig][i], dimColor[i]);            
            //lgp.add(new Constrainer(lgs[ig], WGRAPH, HGRAPH)); // TODO: why doesn't this work
            lgp.add(lgs[ig]);
         }
      }
      else{
         // create the continuous displays
         lgp = new JPanel(nd > 3 ? new GridFillLayout(1, 2, 1, 0) : new BorderLayout());
         lgp.setBackground(null);
         this.add(lgp, BorderLayout.CENTER);
         
         Dimension dim;
         lgs = new LineGraph[2];
         lgs[0] = new LineGraph(seq);
         lgs[0].setDimColor(0, new Color(1.0f, 0.4f, 0.4f));
         if (nd > 1) lgs[0].setDimColor(1, new Color(0.2f, 1.0f, 0.2f));
         if (nd > 2) lgs[0].setDimColor(2, new Color(0.4f, 0.4f, 1.0f));
         lgs[0].setRenderGridLines(false, false);
         lgs[0].setRenderMouseInfo(false, false, false);
         lgs[0].addMouseListener(parent);
         dim = lgs[0].getPreferredSize();
         if (dim.width > nMaxW)
            lgs[0].setVirtualWidth(nMaxW);
         lgp.add(lgs[0]);
         
         if (nd > 3)
         {
            lgs[1] = new LineGraph(seq);
            lgs[1].setDimColor(3, new Color(1.0f, 0.4f, 0.4f));
            if (nd > 4) lgs[1].setDimColor(4, new Color(0.2f, 1.0f, 0.2f));
            if (nd > 5) lgs[1].setDimColor(5, new Color(0.4f, 0.4f, 1.0f));
            lgs[1].setRenderGridLines(false, false);
            lgs[1].setRenderMouseInfo(false, false, false);
            lgs[1].addMouseListener(parent);
            dim = lgs[1].getPreferredSize();
            if (dim.width > nMaxW)
               lgs[1].setVirtualWidth(nMaxW);
            lgp.add(lgs[1]);         
         }
      }

      // create the quantized display
      if (qseq != null)
      {
         qv = new QuantDataView(qseq, Library.generateColors(nQuant), false, qseq.length() > 20?0:1);         
         this.add(qv, BorderLayout.SOUTH);
      }
   }
   
   public void setMaximumSize(Dimension dim)
   {
      maxSize = dim;
      revalidate();
   }
   
   public Dimension getPreferredSize()
   {
      Dimension dim = super.getPreferredSize();
      if (maxSize==null) return dim;
      if (maxSize.width > 0) dim.width = Math.min(dim.width, maxSize.width);
      if (maxSize.height> 0) dim.height= Math.min(dim.height, maxSize.height);
      return dim;
   }

   public void showLineGraphs(boolean bShow)
   {
      bShowCont = bShow;
      updateLayout();
   }

   public void showQuant(boolean bShow)
   {
      bShowQuant = bShow;
      updateLayout();
   }

   protected void updateLayout()
   {
      removeAll();
      
      if (bShowCont)
      {
         add(lgp, BorderLayout.CENTER);
         for(int i=0; i<lgs.length; i++) lgp.add(lgs[i]);
         if (bShowQuant) add(qv, BorderLayout.SOUTH);
      }
      else
      { // quant only
         add(qv, BorderLayout.CENTER);
      }
   }   
}
