package mdisc.tools;

import kdm.data.*;
import kdm.data.transform.*;
import kdm.io.*;
import mdisc.io.*;
import kdm.io.DataLoader.DataLoader;
import kdm.io.Def.*;
import kdm.util.*;
import kdm.gui.*;
import kdm.models.*;
import java.util.*;
import java.util.zip.*;
import java.io.*;
import javax.imageio.*;
import java.awt.*;
import java.awt.image.*;
import javax.swing.*;

public class VizPats
{
   /**
    * Loads the appropriate data subsequence given a pattern and data def.
    */
   public static Sequence loadData(PatternFileLoader.PatData pd, DefData dd)
   {
      StringBuffer sb = new StringBuffer();
      sb.append("<data>\n");
      if (dd.name != null) sb.append("name = " + dd.name + "\n");
      assert (dd.sDataFile != null);
      sb.append("file = " + dd.sDataFile + "\n");
      assert (dd.sDataLoader != null);
      sb.append("loader = " + dd.sDataLoader + "\n");
      if (dd.sDataParams != null) sb.append("params = " + dd.sDataParams + "\n");
      sb.append("</data>\n");

      GuiViewDefLoader gdef = new GuiViewDefLoader();
      gdef.setLoadData(false);
      if (!gdef.loadFromString(null, sb.toString()))
      {
         System.err.println("Error: failed to parse generated data block");
         assert false;
         return null;
      }

      assert (gdef.data.size() == 1);
      DefData ddef = gdef.data.get(0);
      try
      {
         Class cls = Library.getClass(ddef.sDataLoader, "kdm.io.DataLoader");
         DataLoader loader = (DataLoader)cls.newInstance();
         if (!loader.config(gdef.getBasePath(), ddef.sDataParams))
         {
            System.err.println("Error: failed to configure data loader!");
            System.err.println(" class = " + ddef.sDataLoader);
            System.err.println(" params = \"" + ddef.sDataParams + "\"");
            System.exit(1);
         }
         Sequence data = loader.load(ddef.sDataFile);
         for(DataTransform tran : ddef.trans)
            data = tran.transform(data);
         return data.subseq(pd.pos, pd.pos + pd.len);
      } catch (Exception e)
      {
         e.printStackTrace();
         return null;
      }
   }

   /**
    * Generates an image of the given model..
    */
   public static BufferedImage genModelViz(OatesModelUSamp om, double dsig, int kSig)
   {
      Sequence bounds = om.bounds(dsig);
      Sequence data = null;
      // if (kSig == 1) data = bounds.selectDims(new int[]{ 0,1,4,5,8,9 });
      // else data = bounds.selectDims(new int[]{ 2,3,6,7,10,11 });
      if (kSig == 1) data = bounds.selectDims(new int[] { 0, 1, 2, 3, 4, 5 });
      else data = bounds.selectDims(new int[] { 6, 7, 8, 9, 10, 11 });
      LineGraph lgraph = new LineGraph(data);
      lgraph.setDoubleBuffered(false);
      lgraph.setDimColor(0, new Color(1.0f, 0.4f, 0.4f));
      lgraph.setDimColor(1, new Color(1.0f, 0.4f, 0.4f));
      lgraph.setDimColor(2, new Color(0.2f, 1.0f, 0.2f));
      lgraph.setDimColor(3, new Color(0.2f, 1.0f, 0.2f));
      lgraph.setDimColor(4, new Color(0.4f, 0.4f, 1.0f));
      lgraph.setDimColor(5, new Color(0.4f, 0.4f, 1.0f));

      int w = 240; // size of image
      int h = 120;
      int w2 = 48; // width of value legend

      // setup and render the line graph
      lgraph.setVirtualWidth(w);
      BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_3BYTE_BGR);
      Graphics g = img.getGraphics();
      lgraph.setSize(w - w2, h);
      // lgraph.paint(g);
      lgraph.printAll(g);

      // add a value legend to the right side of the image
      ValueLegend valueLegend = new ValueLegend(lgraph);
      valueLegend.setDoubleBuffered(false);
      BufferedImage img2 = new BufferedImage(w2, h, BufferedImage.TYPE_3BYTE_BGR);
      valueLegend.setSize(w2, h);
      // valueLegend.paint(img2.getGraphics());
      valueLegend.printAll(img2.getGraphics());
      g.drawImage(img2, w - w2, 0, lgraph);

      return img;
   }

   /**
    * Generates an image of the signal representing the pd & dd.
    */
   public static BufferedImage genDataViz(Sequence data, PatternFileLoader.PatData pd, DefData dd, int kSig)
   {
      // System.err.println("Gen Data Viz: "+data);
      long ms = data.getLengthMS();

      LineGraph lgraph = new LineGraph(data);
      Color colors[] = { Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.yellow,
            Color.gray, Color.orange, Color.pink, Color.white };
      if (kSig == 1)
      {
         lgraph.setDimColor(0, new Color(1.0f, 0.4f, 0.4f));
         lgraph.setDimColor(1, new Color(0.2f, 1.0f, 0.2f));
         lgraph.setDimColor(2, new Color(0.4f, 0.4f, 1.0f));
      }
      else
      {
         lgraph.setDimColor(3, new Color(1.0f, 0.4f, 0.4f));
         lgraph.setDimColor(4, new Color(0.2f, 1.0f, 0.2f));
         lgraph.setDimColor(5, new Color(0.4f, 0.4f, 1.0f));
      }

      int w = 200; // size of image
      int h = 100;
      int w2 = 32; // width of value legend

      // setup and render the line graph
      BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_3BYTE_BGR);
      Graphics g = img.getGraphics();
      lgraph.setVirtualWidth(w);
      lgraph.setSize(w - w2, h);
      lgraph.paint(g);

      // add a value legend to the right side of the image
      ValueLegend valueLegend = new ValueLegend(lgraph);
      BufferedImage img2 = new BufferedImage(w2, h, BufferedImage.TYPE_3BYTE_BGR);
      valueLegend.setSize(w2, h);
      valueLegend.paint(img2.getGraphics());
      g.drawImage(img2, w - w2, 0, lgraph);

      return img;
   }

   public static void main(String args[])
   {
      if (args.length < 3) // make sure we have enough parameters
      {
         System.err.println();
         System.err
               .println("Usage: java kdm.tools.VizPats <pattern file> <data file> <output dir> <optional model file>");
         System.err.println();
         System.exit(1);
      }

      // load the data def file
      System.err.print("Loading data def file... ");
      DataDefLoader ddl = new DataDefLoader();
      if (!ddl.loadf(args[1]))
      {
         System.err.println("Error: failed to load data definition file");
         System.exit(1);
      }
      System.err.println("done.");

      // load the pattern file
      System.err.print("Loading pattern file... ");
      PatternFileLoader pfl = new PatternFileLoader();
      if (!pfl.load(args[0]))
      {
         System.err.println("Error: failed to load pattern file");
         System.exit(1);
      }
      System.err.println("done.");

      // now we can generate the html file
      try
      {
         OatesModelUSamp om = null;
         int nCols = 4;

         if (args.length == 4)
         {
            System.err.print("Loading model... ");
            FileInputStream fis = new FileInputStream(args[3]);
            GZIPInputStream gzis = new GZIPInputStream(fis);
            ObjectInputStream in = new ObjectInputStream(gzis);
            om = (OatesModelUSamp)in.readObject();
            in.close();
            System.err.println("done.");
         }

         String path = Library.ensurePathSep(args[2]);
         File dir = new File(path);
         dir.mkdirs();
         PrintWriter out = new PrintWriter(new FileWriter(path + "index.html"));
         out.println("<html><head><title>Pattern Vizualization</title><style>");
         out.println("div.pat { border: 2px solid #448; background-color: #ccf;");
         out.println("padding: 8px;  margin-top: 24px; }");
         out.println("div.occ { border: 1px solid black; background-color: #cde;");
         out.println("padding: 6px; margin: 4px; text-align: center; font-family: sans-serif; }");
         out.println("body { background-color: #000; }");
         out.println("</style></head><body>");

         for(int iPat = 0; iPat < pfl.pats.size(); iPat++)
         {
            StringBuffer sbSignal1 = new StringBuffer();
            ArrayList<PatternFileLoader.PatData> occs = pfl.get(iPat);
            System.err.println("Processing Pattern " + (iPat + 1) + " / " + pfl.pats.size() + "  ("
                  + occs.size() + ")");
            sbSignal1.append("<div class=\"pat\"><b><font size=\"+2\">Pattern " + (iPat + 1)
                  + "</font></b> (" + occs.size() + " occurrences)<br><br>\n");

            sbSignal1.append("<table width=\"100%\">\n");
            for(int i = 0; i < occs.size(); i++)
            {
               PatternFileLoader.PatData pd = occs.get(i);
               DefData dd = ddl.data.get(pd.iSeries - 1);
               Sequence data = loadData(pd, dd);
               String sScore = "";
               if (om != null) sScore = String.format("%.3f", om.eval(data));

               // dump a visualization of the signal
               BufferedImage img = genDataViz(data, pd, dd, 1);
               String sImg = "pat-viz-" + (iPat + 1) + "-" + (i + 1) + "a.jpg";
               ImageIO.write(img, "jpg", new File(path + sImg));
               if ((i % nCols) == 0) sbSignal1.append("<tr><td>\n");
               else sbSignal1.append("<td>\n");
               sbSignal1.append("<div class=\"occ\">Occurrence " + (i + 1) + ": " + sScore + "<br><i>"
                     + dd.name + " (" + pd.pos + " -> " + (pd.pos + pd.len) + ")</i><br>\n");
               sbSignal1.append("<img src=\"" + sImg + "\"></div>\n");
               if ((i % nCols) == (nCols - 1)) sbSignal1.append("</td></tr>\n");
               else sbSignal1.append("</td>\n");
            }

            // finish off the last row of the table
            int iOcc = (occs.size() % nCols);
            if (iOcc > 0) sbSignal1.append("<td columnspan=\"" + (nCols - iOcc) + "\"/></tr>\n");
            sbSignal1.append("</table></div>\n");
            out.println(sbSignal1);
         }

         // display the model if we have one
         if (args.length == 4)
         {
            double dsig = 1.0;
            out.println(String.format("Model (%.2f sigma bounds):<br><br>", dsig));

            BufferedImage img = genModelViz(om, dsig, 1);
            String sImg = "model-viz-a.jpg";
            ImageIO.write(img, "jpg", new File(path + sImg));
            out.println("<img src=\"" + sImg + "\">");

            /*
             * img = genModelViz(om, dsig, 2); sImg = "model-viz-b.jpg";
             * ImageIO.write(img, "jpg", new File(path+sImg)); out.println("<img
             * src=\""+sImg+"\">");
             */

            out.println("<br><br><hr><br><br>");
         }

         out.println("</body></html>");
         out.close();
      } catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }
}
