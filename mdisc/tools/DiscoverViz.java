package mdisc.tools;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.io.*;
import kdm.io.DataLoader.*;
import kdm.io.Def.*;
import kdm.gui.*;
import kdm.models.*;
import kdm.metrics.*;
import kdm.util.*;
import kdm.mlpr.*;
import kdm.mlpr.dataTree.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

import gnu.getopt.*;

import static java.awt.BorderLayout.*;
import static javax.swing.BoxLayout.*;

public class DiscoverViz implements ActionListener, AdjustmentListener, MouseListener
{
   protected boolean bVerbose = false;
   protected JFrame frame;
   protected GraphComplexLite graphc;
   protected LineGraph lgraph, lgExample;
   protected EventBar eventBar;
   protected JComboBox comboSeries;
   protected JLabel lbSeries, lbWidth, lbRadius, lbK, lbStatus;
   protected JScrollBar sbWidth, sbRadius, sbK;
   protected JButton btMakeTree;
   protected JPanel knnPanel, kdePanel;
   protected ArrayList<Sequence> tseries;
   protected ArrayList<MarkupSet> marks;
   protected ArrayList<DefData> ddefs;
   protected int nDims = 0;
   protected int iSeries = -1;
   protected WindowLocation example;
   protected int seriesStartWindowIndex[];
   protected Color colorSeries = Color.green;
   protected Color colorMatch = Color.red;
   protected Color colorMatchEx = new Color(.5f, .2f, .1f);
   protected Color colorExample = Color.yellow;
   protected KnnSpTree knn;
   protected ArrayList<KdeResult> kde;
   protected JButton[] btKde;

   public void usage()
   {
      System.err.println("USAGE: java kdm.tools.DiscoverViz [options] <data def file>");
      System.err.println();
      System.err.println(" Options:");
      System.err.println("  -?                      display this help message");
      System.err.println("  -v                      display verbose output");
      System.err.println();
   }

   public boolean layoutTop(JPanel parent)
   {
      parent.setLayout(new BorderLayout());

      // layout the line graph
      lgraph = new LineGraph();
      lgraph.addMouseListener(this);
      graphc = new GraphComplexLite(lgraph);
      GraphScrollPane gsp = GraphScrollPane.construct(new GraphComplexLite[] { graphc }, graphc, true, true, true, null, new ZoomReal(graphc));
      eventBar = new EventBar(lgraph);
      eventBar.addGraph(graphc);
      eventBar.setScrollPane(gsp);
      gsp.addBorderComponent(eventBar, BorderLayout.SOUTH);

      // layout the header
      JPanel top = new JPanel();
      top.setLayout(new FlowLayout(FlowLayout.LEFT));
      comboSeries = new JComboBox();
      comboSeries.addActionListener(this);
      for(int i = 0; i < tseries.size(); i++)
         comboSeries.addItem("Series " + (i + 1));
      top.add(comboSeries);
      lbSeries = new JLabel();
      top.add(lbSeries);
      comboSeries.setSelectedIndex(0);

      parent.add(top, BorderLayout.NORTH);
      parent.add(gsp, BorderLayout.CENTER);
      parent.setMinimumSize(new Dimension(32, 100));
      return true;
   }

   public boolean layoutRight(JPanel parent)
   {
      Box box = Box.createVerticalBox();
      JPanel p;

      parent.setBorder(BorderFactory.createEmptyBorder(8, 8, 2, 4));
      parent.setLayout(new BorderLayout());
      parent.add(box, BorderLayout.NORTH);

      // layout width slider
      p = new JPanel(new BorderLayout(4, 4));
      p.add(new JLabel("Width:"), BorderLayout.WEST);
      sbWidth = new JScrollBar(JScrollBar.HORIZONTAL, 25, 0, 5, 100);
      sbWidth.addAdjustmentListener(this);
      p.add(sbWidth, BorderLayout.CENTER);
      lbWidth = new JLabel("" + sbWidth.getValue());
      lbWidth.setPreferredSize(new Dimension(24, 0));
      p.add(lbWidth, BorderLayout.EAST);
      box.add(p);

      box.add(Box.createVerticalStrut(4));

      // layout R slider
      p = new JPanel(new BorderLayout(4, 4));
      p.add(new JLabel("Radius:"), BorderLayout.WEST);
      sbRadius = new JScrollBar(JScrollBar.HORIZONTAL, 10, 0, 1, 100);
      sbRadius.addAdjustmentListener(this);
      p.add(sbRadius, BorderLayout.CENTER);
      lbRadius = new JLabel("" + sbRadius.getValue());
      lbRadius.setPreferredSize(new Dimension(24, 0));
      p.add(lbRadius, BorderLayout.EAST);
      box.add(p);

      box.add(Box.createVerticalStrut(8));

      // layout K slider
      p = new JPanel(new BorderLayout(4, 4));
      p.add(new JLabel("K:"), BorderLayout.WEST);
      sbK = new JScrollBar(JScrollBar.HORIZONTAL, 6, 0, 1, 20);
      sbK.addAdjustmentListener(this);
      p.add(sbK, BorderLayout.CENTER);
      lbK = new JLabel("" + sbK.getValue());
      lbK.setPreferredSize(new Dimension(24, 0));
      p.add(lbK, BorderLayout.EAST);
      box.add(p);

      box.add(Box.createVerticalStrut(8));

      // layout exemplar
      lgExample = new LineGraph();
      lgExample.setDimColor(0, colorExample);
      lgExample.setShowAll(true);
      GraphComplexLite gcl = new GraphComplexLite(lgExample);
      gcl.setBorder(BorderFactory.createLoweredBevelBorder());
      box.add(gcl);

      // layout "construct tree" button
      p = new JPanel(new BorderLayout());
      p.setBorder(BorderFactory.createEmptyBorder(0, 4, 18, 4));
      btMakeTree = new JButton("Go!");
      btMakeTree.addActionListener(this);
      btMakeTree.setEnabled(false);
      p.add(btMakeTree, BorderLayout.CENTER);
      parent.add(p, BorderLayout.SOUTH);

      return true;
   }

   public boolean create(String args[])
   {
      int c;
      LongOpt[] longopts = new LongOpt[2];
      longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');
      longopts[0] = new LongOpt("v", LongOpt.NO_ARGUMENT, null, 'v');

      Getopt g = new Getopt("SupTest", args, "?", longopts, true);
      while((c = g.getopt()) != -1)
      {
         switch(c){
         case '?':
         case 'h': // help
            usage();
            System.exit(0);
            break;
         case 'v': // verbose
            bVerbose = true;
            break;
         default:
            System.err.println("unrecognized command line option: " + c);
            System.exit(1);
            break;
         }
      }

      // make sure that a data def file was specified
      if (g.getOptind() >= args.length)
      {
         System.err.println("Error: no data definition file specified!");
         System.err.println();
         usage();
         return false;
      }

      // load the data def file
      if (bVerbose) System.err.print("Loading data def file... ");
      DataDefLoader ddl = new DataDefLoader();
      if (!ddl.loadf(args[g.getOptind()]))
      {
         System.err.println("Error: failed to load data definition file");
         System.exit(1);
      }
      if (bVerbose) System.err.println("done.");

      tseries = new ArrayList<Sequence>();
      marks = new ArrayList<MarkupSet>();
      ddefs = new ArrayList<DefData>();
      for(int iData = 0; iData < ddl.data.size(); iData++)
      {
         DefData ddef = ddl.data.get(iData);
         ddefs.add(ddef);
         Sequence seq;

         // load the data and associated labels
         try
         {
            // first load the whole data sequence
            Class cls = Library.getClass(ddef.sDataLoader, "kdm.io.DataLoader");
            DataLoader dloader = (DataLoader)cls.newInstance();
            if (!dloader.config(ddl.getBasePath(), ddef.sDataParams))
            {
               System.err.println("Error: failed to configure data loader!");
               System.err.println(" class = " + ddef.sDataLoader);
               System.err.println(" params = \"" + ddef.sDataParams + "\"");
               System.exit(1);
            }
            seq = dloader.load(ddef.sDataFile);
            if (nDims == 0)
            {
               nDims = seq.getNumDims(); // remember dimensionality
               System.err.printf("Dimensions: %d\n", nDims);
               if (nDims > 1)
               {
                  System.err.println("Error: currently DiscoverViz only supports univariate sequences");
                  return false;
               }
            }
            tseries.add(seq);

            // now load the labels
            cls = Library.getClass(ddef.sLabelLoader, "kdm.io");
            MarkupLoader mloader = (MarkupLoader)cls.newInstance();
            if (!mloader.config(ddl.getBasePath(), ddef.sLabelParams))
            {
               System.err.println("Error: failed to configure markup loader!");
               System.err.println(" class = " + ddef.sLabelLoader);
               System.err.println(" params = \"" + ddef.sLabelParams + "\"");
               System.exit(1);
            }
            MarkupSet mark = mloader.load(ddef.sLabelFile);
            marks.add(mark);

            if (bVerbose)
            {
               System.err.print("Data File: " + Library.getFileName(ddef.sDataFile));
               System.err.println("   Labels: " + Library.getFileName(ddef.sLabelFile));
               System.err.print(" Summary: " + seq.length() + " frames, " + seq.getNumDims() + " dims, ");
               System.err.println(mark.size() + " labels");
            }
         } catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }
      }

      frame = new JFrame("DiscoverViz: " + args[0]);
      frame.setLocation(120, 80);
      frame.setSize(1000, 800);
      frame.addWindowListener(new WindowAdapter() {
         public void windowClosing(WindowEvent e)
         {
            close();
         }
      });

      frame.setLayout(new BorderLayout());

      JPanel top = new JPanel();
      frame.add(top, BorderLayout.NORTH);

      JTabbedPane tab = new JTabbedPane(JTabbedPane.BOTTOM, JTabbedPane.SCROLL_TAB_LAYOUT);
      frame.add(tab, BorderLayout.CENTER);

      // status bar
      lbStatus = new JLabel("Click on the graph to select an exemplar", JLabel.CENTER);
      lbStatus.setBorder(BorderFactory.createEmptyBorder(4, 2, 2, 2));
      frame.add(lbStatus, BorderLayout.SOUTH);

      // knn panel
      knnPanel = new JPanel();
      VerticalScrollPanel vsp = new VerticalScrollPanel(4, 20);
      vsp.setLayout(new BorderLayout());
      vsp.add(knnPanel, BorderLayout.CENTER);
      JScrollPane sp = new JScrollPane(vsp, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      tab.addTab("KNN", sp);

      // kde panel
      kdePanel = new JPanel();
      vsp = new VerticalScrollPanel(4, 20);
      vsp.setLayout(new BorderLayout());
      vsp.add(kdePanel, BorderLayout.CENTER);
      sp = new JScrollPane(vsp, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      tab.addTab("KDE", sp);

      JPanel right = new JPanel();
      frame.add(right, BorderLayout.EAST);

      if (!layoutTop(top)) return false;
      if (!layoutRight(right)) return false;
      right.setPreferredSize(new Dimension(200, 100));
      top.setPreferredSize(new Dimension(1, 250));

      frame.setVisible(true);
      return true;
   }

   public void close()
   {
      frame.setVisible(false);
      frame.dispose();
      System.exit(0);
   }

   public void actionPerformed(ActionEvent e)
   {
      String cmd = e.getActionCommand();
      Object target = e.getSource();

      if (target == comboSeries)
      {
         JComboBox pd = (JComboBox)target;
         iSeries = pd.getSelectedIndex();
         Sequence seq = tseries.get(iSeries);
         graphc.setData(seq);
         graphc.setDimColor(0, colorSeries);
         if (lbSeries != null) lbSeries.setText(String.format(" Length: %d frames", seq.length()));
         frame.repaint();
      }
      else if (target == btMakeTree)
      {
         double tau = (double)sbRadius.getValue();
         double rho = 0.7;

         double[][] data = slideWindow(example.nLength);

         TimerMS timer = new TimerMS();
         knn = new KnnSpTree(data, tau, rho);
         System.err.println("KNN build time: " + timer.time());
         System.err.println("KNN Spill: " + knn.root.numSpill() + " / " + knn.root.numNodes());
         displayKNN();
         btMakeTree.setEnabled(false);

         {
            kde = new ArrayList<KdeResult>();

            timer.reset();
            int n = data.length;
            double dist[][] = new double[n][n];
            double den[] = new double[n];
            double kw = tau / 3.0;
            double a = 1.0 / (Library.SQRT_2PI * Math.sqrt(kw));
            for(int i = 0; i < n; i++)
               for(int j = i + 1; j < n; j++)
               {
                  double y = Library.dist(data[i], data[j]);
                  dist[i][j] = dist[j][i] = y;
                  double x = Math.exp(-0.5 * y * y / kw);
                  den[i] += x;
                  den[j] += x;
               }
            System.err.println("All KDE direct: " + timer.time());

            // find the kde maxima
            timer.reset();
            while(true)
            {
               // find max
               int iBest = -1;
               double vBest = Library.NEGINF;
               for(int i = 0; i < n; i++)
                  if (den[i] > vBest)
                  {
                     vBest = den[i];
                     iBest = i;
                  }

               if (iBest < 0) break; // exit if no more maxima
               double denx = den[iBest];

               // remove all NN (including this point) from the list
               MyIntList nbrs = new MyIntList();
               for(int i = 0; i < n; i++)
               {
                  if (dist[iBest][i] <= tau)
                  {
                     den[i] = Library.NEGINF;
                     nbrs.add(i);
                  }
               }

               if (nbrs.size() < 3) break; // can't have too few neighbors

               // add the new maxima
               kde.add(new KdeResult(iBest, nbrs, data, a * denx));
            }
            System.err.println("Find maxima time: " + timer.time());

            displayKDE();
         }
      }
      else if (btKde != null)
      {
         assert false;
         /*
          * for(int i=0; i<btKde.length; i++) { if (target != btKde[i]) continue;
          * ArrayList<WindowLocation> locs = new ArrayList<WindowLocation>(); for(int
          * j=0; j<kde.get(i).nbrs.length; j++)
          * locs.add(getLocation(kde.get(i).nbrs[j])); VideoFrame vf = new
          * VideoFrame(ddefs, locs); vf.setVisible(true); }
          */
      }
      else System.err.println("Unhandled Action: " + cmd);
   }

   public void adjustmentValueChanged(AdjustmentEvent e)
   {
      Object target = e.getSource();

      if (target == sbWidth)
      {
         lbWidth.setText("" + sbWidth.getValue());
         if (example != null) setExample(example.iSeries, example.iStart + example.nLength / 2);
      }
      else if (target == sbRadius)
      {
         lbRadius.setText("" + sbRadius.getValue());
         if (example != null) btMakeTree.setEnabled(true);
      }
      else if (target == sbK)
      {
         lbK.setText("" + sbK.getValue());
         displayKNN();
      }
      else System.err.println("changed: " + target);
   }

   public void mouseClicked(MouseEvent e)
   {
      Object target = e.getSource();

      if (target == lgraph)
      {
         TimeX tx = lgraph.getTimeX(e.getX());
         setExample(iSeries, tx.index);
         if (example != null) lbStatus.setText(" ");
         knn = null;
         displayKNN();
      }
   }

   protected void displayKNN()
   {
      if (knn == null)
      {
         knnPanel.removeAll();
         knnPanel.revalidate();
         knnPanel.repaint();
         return;
      }
      assert example != null;

      int k = sbK.getValue();
      Sequence seq = tseries.get(example.iSeries);
      double[] x = seq.extractDim(0, example.iStart, example.nLength);
      int[] iknn = knn.findi(x, k);
      int nCols = 3;

      knnPanel.removeAll();
      knnPanel.setLayout(new BoxLayout(knnPanel, BoxLayout.Y_AXIS));
      knnPanel.add(Box.createVerticalStrut(8));
      Box hbox = null;

      for(int i = 0; i < k; i++)
      {
         if ((i % nCols) == 0)
         {
            hbox = Box.createHorizontalBox();
            hbox.add(Box.createHorizontalStrut(8));
            knnPanel.add(hbox);
            knnPanel.add(Box.createVerticalStrut(8));
         }

         Sequence seq1 = new Sequence("NN-" + (i + 1), knn.data[iknn[i]]);
         Sequence seqn = seq1.addDims(lgExample.getData());
         LineGraph lg = new LineGraph(seqn);
         lg.setDimColor(0, colorMatch);
         lg.setDimColor(1, colorMatchEx);
         lg.setShowAll(true);
         GraphComplexLite gcl = new GraphComplexLite(lg);
         gcl.setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, Color.gray));

         JPanel p = new JPanel(new BorderLayout(4, 4));
         p.setBackground(Color.lightGray);
         p.setBorder(BorderFactory.createLoweredBevelBorder());
         p.add(gcl, BorderLayout.CENTER);
         double dist = Library.dist(x, knn.getData(iknn[i]));
         WindowLocation loc = getLocation(iknn[i]);
         JLabel label = new JLabel(String.format("%d.%d (%d)  Dist: %.3f", loc.iSeries + 1, loc.iStart,
               loc.nLength, dist), JLabel.CENTER);
         label.setBorder(BorderFactory.createEmptyBorder(0, 2, 2, 2));
         p.add(label, BorderLayout.SOUTH);

         hbox.add(p);
         hbox.add(Box.createHorizontalStrut(8));
      }

      // add filler components to keep columns consistent
      while(hbox.getComponentCount() < (1 + nCols * 2))
      {
         // TODO: why is there a shift when third comp is added?
         hbox.add(Box.createHorizontalStrut(10000));
         hbox.add(Box.createHorizontalStrut(8));
      }

      knnPanel.revalidate();
      knnPanel.repaint();
   }

   protected void displayKDE()
   {
      if (kde == null || kde.isEmpty())
      {
         kdePanel.removeAll();
         kdePanel.revalidate();
         kdePanel.repaint();
         return;
      }

      int nCols = 3;

      int k = kde.size();
      kdePanel.removeAll();
      kdePanel.setLayout(new BoxLayout(kdePanel, BoxLayout.Y_AXIS));
      kdePanel.add(Box.createVerticalStrut(8));
      Box hbox = null;
      btKde = new JButton[k];

      for(int i = 0; i < k; i++)
      {
         if ((i % nCols) == 0)
         {
            hbox = Box.createHorizontalBox();
            hbox.add(Box.createHorizontalStrut(8));
            kdePanel.add(hbox);
            kdePanel.add(Box.createVerticalStrut(8));
         }

         KdeResult kr = kde.get(i);
         Sequence seq = new Sequence("KDE: " + (i + 1), kr.data[kr.index]);
         LineGraph lg = new LineGraph(seq);
         lg.setDimColor(0, colorMatch);
         lg.setShowAll(true);
         GraphComplexLite gcl = new GraphComplexLite(lg);
         gcl.setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, Color.gray));

         JPanel p = new JPanel(new BorderLayout(4, 4));
         p.setBackground(Color.lightGray);
         p.setBorder(BorderFactory.createLoweredBevelBorder());
         p.add(gcl, BorderLayout.CENTER);
         JPanel q = new JPanel();
         WindowLocation loc = getLocation(kr.index);
         JLabel label = new JLabel(String.format("%d.%d  Den: %.3f  N: %d", loc.iSeries + 1, loc.iStart,
               kr.den, kr.nbrs.length), JLabel.CENTER);
         label.setBorder(BorderFactory.createEmptyBorder(0, 2, 2, 2));
         q.add(label);
         btKde[i] = new JButton("Video");
         btKde[i].addActionListener(this);
         q.add(btKde[i]);
         p.add(q, BorderLayout.SOUTH);

         hbox.add(p);
         hbox.add(Box.createHorizontalStrut(8));
      }

      // add filler components to keep columns consistent
      while(hbox.getComponentCount() < (1 + nCols * 2))
      {
         // TODO: why is there a shift when third comp is added?
         hbox.add(Box.createHorizontalStrut(10000));
         hbox.add(Box.createHorizontalStrut(8));
      }

      kdePanel.revalidate();
      kdePanel.repaint();
   }

   protected void setExample(int iSeries, int ix)
   {
      int w = sbWidth.getValue();
      int iStart = ix - w / 2;
      int iEnd = iStart + w;
      Sequence seq = tseries.get(iSeries);
      if (iStart >= 0 && iEnd <= seq.length())
      {
         lgExample.setData(seq.subseq(iStart, iEnd));
         lgExample.setShowAll(true);
         example = new WindowLocation(iSeries, iStart, w);
         btMakeTree.setEnabled(true);
      }
      else
      {
         example = null;
         lgExample.setData(null);
         btMakeTree.setEnabled(false);
      }
   }

   protected WindowLocation getLocation(int ix)
   {
      int iSeries = 0;
      while(iSeries + 1 < tseries.size() && ix > seriesStartWindowIndex[iSeries + 1])
         iSeries++;
      int iStart = (ix - seriesStartWindowIndex[iSeries]) * getSkip(example.nLength);
      return new WindowLocation(iSeries, iStart, example.nLength);
   }

   protected int getSkip(int w)
   {
      return Library.max(1, w / 2);
   }

   protected double[][] slideWindow(int w)
   {
      double[][] data;
      int nSkip = getSkip(w);

      // allocate space for data
      int n = 0;
      for(int iSeq = 0; iSeq < tseries.size(); iSeq++)
         n += Library.getNumSlidingWindowSites(tseries.get(iSeq).length(), w, nSkip);
      System.err.println("Num window sites: " + n);
      data = new double[n][w];

      if (seriesStartWindowIndex == null) seriesStartWindowIndex = new int[tseries.size()];

      // extract the data points
      int id = 0;
      for(int iSeq = 0; iSeq < tseries.size(); iSeq++)
      {
         seriesStartWindowIndex[iSeq] = id;
         Sequence seq = tseries.get(iSeq);
         int len = seq.length();
         for(int i = 0; i + w <= len; i += nSkip)
         {
            for(int j = 0; j < w; j++)
               data[id][j] = seq.get(i + j, 0);
            id++;
         }
      }

      return data;
   }

   public void mouseEntered(MouseEvent e)
   {}

   public void mousePressed(MouseEvent e)
   {}

   public void mouseReleased(MouseEvent e)
   {}

   public void mouseExited(MouseEvent e)
   {}

   public static void main(String args[])
   {
      DiscoverViz view = new DiscoverViz();
      view.create(args);
   }

}

// ////////////////////////////////////////////////////////////////////

class KdeResult
{
   public int index;
   public double data[][];
   public double den;
   public int nbrs[];

   public KdeResult(int _index, MyIntList _nbrs, double _data[][], double _den)
   {
      index = _index;
      data = _data;
      den = _den;
      nbrs = _nbrs.toArray();
   }
}
