package mdisc.peruse;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.models.*;
import kdm.io.*;
import kdm.gui.*;
import kdm.util.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

import static java.awt.BorderLayout.*;
import static javax.swing.BoxLayout.*;

public class DOGui extends JFrame
{
   protected int nbw = 3;
   protected LineGraph lgex[];
   protected JLabel lbex[];

   public DOGui()
   {
      super("Discover GUI");

      setLocation(120, 80);
      setSize(800, 600);

      setLayout(new BorderLayout());
      Box vbox = Box.createVerticalBox();
      JPanel row = new JPanel();
      row.setLayout(new GridLayout(1, nbw));
      add(vbox, BorderLayout.CENTER);
      vbox.add(row);

      lgex = new LineGraph[nbw];
      lbex = new JLabel[nbw];
      for(int i = 0; i < nbw; i++){
         lgex[i] = new LineGraph();
         lgex[i].setShowAll(true);
         GraphComplexLite gcl = new GraphComplexLite(lgex[i]);
         lbex[i] = new JLabel();
         JPanel p = new JPanel();
         p.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
         p.setLayout(new BorderLayout());
         p.add(gcl, BorderLayout.CENTER);
         p.add(lbex[i], BorderLayout.SOUTH);
         row.add(p);
      }
   }

   public void setBestWindow(OatesWindow win, ArrayList<Sequence> vdata)
   {
      for(int i = nbw - 1; i > 0; i--)
         if (lgex[i - 1].getData() != null){
            lgex[i].setData(lgex[i - 1].getData());
            lgex[i].setShowAll(true);
            lgex[i].setDimColor(0, Color.red);
            lgex[i].setDimColor(1, Color.green);
            lgex[i].setDimColor(2, Color.blue);
            // lgex[i].setDimColor(3, Color.magenta);
            // lgex[i].setDimColor(4, Color.cyan);
            // lgex[i].setDimColor(5, Color.yellow); // TODO
            lbex[i].setText(lbex[i - 1].getText());
         }

      Sequence seq = vdata.get(win.loc.iSeries);
      int ia = win.loc.start();
      int ib = win.loc.end();
      Sequence sub = seq.subseq(ia, ib);

      lgex[0].setData(sub);
      lgex[0].setShowAll(true);
      lgex[0].setDimColor(0, Color.red);
      lgex[0].setDimColor(1, Color.green);
      lgex[0].setDimColor(2, Color.blue);
      // lgex[0].setDimColor(3, Color.magenta);
      // lgex[0].setDimColor(4, Color.cyan);
      // lgex[0].setDimColor(5, Color.yellow); // TODO
      lbex[0].setText(String.format("%d.%d (%d)  %.4f", win.loc.iSeries + 1, win.loc.iStart,
            win.loc.nLength, win.score));
   }
}
