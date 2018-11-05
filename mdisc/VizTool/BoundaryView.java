package mdisc.VizTool;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.filechooser.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import kdm.data.*;
import kdm.gui.*;
import kdm.mlpr.*;
import kdm.models.*;
import kdm.models.misc.*;
import kdm.util.*;
import kdm.tools.*;
import kdm.tools.SupTest.HmmEval;
import kdm.tools.SupTest.HmmTrain;

/**
 * View unsupervised info about each series, such as potential boundary locations 
 */
public class BoundaryView extends JPanel
{
   protected GlobalData gdata;
   
   public BoundaryView(GlobalData gdata)
   {
      super(new BorderLayout());
      gdata.boundaryView = this;
      setOpaque(true);

      this.gdata = gdata;
   }
   
}