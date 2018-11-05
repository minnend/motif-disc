package mdisc.VizTool;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import kdm.data.*;
import kdm.gui.*;
import kdm.models.*;
import kdm.util.*;
import org.apache.commons.collections.primitives.*;

public class DimensionOverview extends JPanel
{
   public static final int DEF_VIEW_SIZE = 100;

   protected int VIEW_SIZE = DEF_VIEW_SIZE;
   protected int d, wPref;
   protected GlobalData gdata;
   protected Color cBorder;
   protected SeqStatsView ssvModelVar, ssvModelData[];

   public DimensionOverview(int _d, GlobalData _gdata)
   {
      d = _d;
      gdata = _gdata;
      int nClasses = gdata.getNumClasses();
      setLayout(new GridFlexLayout(1, 4, 8, 0));
      setBackground(Color.darkGray);
      cBorder = getBackground().brighter();
      setupGrid();
      setBorder(BorderFactory.createEmptyBorder(4, 8, 4, 8));

      // text overview / info
      JTextArea ta = new JTextArea(10, 10);
      ta.setEditable(false);
      ta.setBackground(Color.lightGray);
      ta.setBorder(BorderFactory.createLineBorder(Color.black));
      this.add(ta);

      // all data view
      double[] x = new double[gdata.nTotalSeriesLength];
      int ix = 0;
      for(Sequence seq : gdata.tseries)
         for(int i = 0; i < seq.length(); i++)
            x[ix++] = seq.get(i, d);
      SeqStatsView ssv = new SeqStatsView(x);
      ssv.setBorder(BorderFactory.createLineBorder(cBorder));
      this.add(ssv);
      ta.setText(String.format("\n Dimension: %d\n Mean:\n   %.6f\n Variance:\n   %.6f", d, ssv.getMean(),
            ssv.getVar()));

      // all class gaussians
      ssv = new SeqStatsView(x);
      ssv.setRenderHist(false);
      ssv.setBorder(BorderFactory.createLineBorder(cBorder));
      Color colors[] = Library.generateColors(nClasses);
      if (gdata.labData != null)
      {
         Iterator<String> itClasses = gdata.labData.keySet().iterator();
         for(int iClass = 0; iClass < nClasses; iClass++)
         {
            String sClass = itClasses.next();
            ArrayList<Sequence> seqs = gdata.labData.get(sClass);
            GaussianDyn1D gmd = new GaussianDyn1D();
            for(Sequence seq : seqs)
               gmd.add(seq.extractDim(d));
            gmd.setReport(ProbFVModel.Report.prob);
            ssv.addGauss(gmd, colors[iClass]);
         }
      }
      this.add(ssv);

      // distr of variances
      ssvModelVar = new SeqStatsView();
      ssvModelVar.setBorder(BorderFactory.createLineBorder(Color.darkGray.brighter()));
      this.add(ssvModelVar);
   }

   protected void setupGrid()
   {
      GridFlexLayout gfl = (GridFlexLayout)getLayout();
      gfl.setColumns(GridFlexLayout.Style.fixed, VIEW_SIZE);
      gfl.setColumn(0, GridFlexLayout.Style.pref);
      gfl.setRows(GridFlexLayout.Style.fixed, VIEW_SIZE);
   }

   public void setViewSize(int vs)
   {
      VIEW_SIZE = vs;
      setupGrid();
      revalidate();
   }

   /**
    * Remove all ssvModelData views from this container and clear the array.
    */
   protected void clearModelDataSSVs()
   {
      if (ssvModelData == null) return;
      for(int i = 0; i < ssvModelData.length; i++)
         remove(ssvModelData[i]);
      ssvModelData = null;
   }

   /**
    * Generate views that depend on trained, aligned models
    * 
    * @param iModelType type of model (oates, hmm, etc.) defined in Word Spotting View
    * @param iClass index of the class to display
    * @param models
    */
   public void calcModelViews(int iModelType, int iClass, ProbSeqModel[] models)
   {
      clearModelDataSSVs();
      int nDims = gdata.getNumDims();
      switch(iModelType){
      case 0: // oates
      {
         // calc total length of all models
         int N = 0;
         for(int i = 0; i < models.length; i++)
         {
            OatesModelUSamp oates = (OatesModelUSamp)models[i];
            N += oates.length();
         }

         // get the variances in each state in each model
         double[][] x = new double[nDims][N];
         int ix = 0;
         for(int i = 0; i < models.length; i++)
         {
            OatesModelUSamp oates = (OatesModelUSamp)models[i];
            for(int j = 0; j < oates.length(); j++)
            {
               for(int d = 0; d < nDims; d++)
                  x[d][ix] = oates.getModel(j, d).getVar();
               ix++;
            }
         }
         assert (ix == N);
         for(int d = 0; d < nDims; d++)
         {
            ssvModelVar.setData(x[d]);
            ssvModelVar.setHistColor(new Color(50, 10, 180));
         }

         // figure out which data goes in each state
         OatesModelUSamp oates = (OatesModelUSamp)models[iClass];
         int nStates = oates.length();
         ArrayDoubleList[] stateData = new ArrayDoubleList[nStates];
         for(int i = 0; i < nStates; i++)
            stateData[i] = new ArrayDoubleList();
         ArrayList<Sequence> examples = gdata.labData.get(gdata.classes[iClass]);
         for(Sequence seq : examples)
         {
            OatesMapping omap = oates.align(seq);
            for(int i = 0; i < nStates; i++)
               stateData[i].add(seq.get(omap.imap[i], d));
         }

         // display the results
         GridFlexLayout gfl = (GridFlexLayout)getLayout();
         gfl.setNumCols(4 + nStates, GridFlexLayout.Style.fixed, VIEW_SIZE);
         ssvModelData = new SeqStatsView[nStates];
         for(int i = 0; i < nStates; i++)
         {
            ssvModelData[i] = new SeqStatsView();
            ssvModelData[i].setBorder(BorderFactory.createLineBorder(Color.darkGray.brighter()));
            this.add(ssvModelData[i]);
         }
         for(int i = 0; i < nStates; i++)
         {
            ssvModelData[i].setData(stateData[i].toArray());
            ssvModelData[i].setToolTipText(String.format("kurt: %.2f", ssvModelData[i].getKurtosis()));
         }
         revalidate();
      }
         break;
      case 1: // hmm
         System.err.println("Error: hmm not currently supported");
         return;
      case 2: // dhmm
         System.err.println("Error: dhmm not currently supported");
         return;
      default:
         System.err.printf("Error: unsupported model type (%d)\n", iModelType);
         return;
      }
   }
}
