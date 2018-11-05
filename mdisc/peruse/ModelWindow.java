package mdisc.peruse;

import kdm.models.*;
import java.util.*;

public class ModelWindow
{
    public ProbSeqModel model;
    public OatesWindow window;
    
    public ModelWindow()
    {
        model = null;
        window = null;
    }

    public ModelWindow(ProbSeqModel _model, OatesWindow _window)
    {
        model = _model;
        window = _window;
    }
}
