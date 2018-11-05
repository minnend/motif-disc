package mdisc.peruse;

import kdm.data.*;
import kdm.io.*;
import kdm.util.*;
import kdm.models.misc.*;

public class OatesWindow
{
    public WindowLocation loc;
    public double score;
    public NBestOMaps maps;

    public OatesWindow()
    {
        set(-1, -1, 0, Double.NEGATIVE_INFINITY, null);
    }

    public OatesWindow(int _iSeries, int _iw, int _nw, double _score, NBestOMaps _maps)
    {
        set(_iSeries, _iw, _nw, _score, _maps);
    }

    public void set(int _iSeries, int _iw, int _nw, double _score, NBestOMaps _maps)
    {
        loc = new WindowLocation(_iSeries, _iw, _nw);
        score = _score;
        maps = _maps;
    }

    public String toString()
    {
        return String.format("[OWin: %d.%d %d %.2f | %d]",
                             loc.iSeries+1, loc.iStart, loc.nLength, score, maps.size());
        //return "[OWin: "+(loc.iSeries+1)+"."+loc.iStart+" "+loc.nLength+" "+score+" |"+maps.size()+"]";
    }
}
