package mdisc.peruse;

import java.util.*;
import java.io.*;

import kdm.data.*;
import kdm.models.*;

public class Params
{
    public static final int VER = 1001;

    public static final int MODEL_HMM = 1;
    public static final int MODEL_OATES = 2;

    public static final int HMM_VITERBI = 1;
    public static final int HMM_FORWARD = 2;
    public static final int HMM_BW = 3;

    public static final int SEARCH_GIVEN = 0;
    public static final int SEARCH_EXHAUSTIVE = 1;
    public static final int SEARCH_EXHAUST_ONE = 2;
    public static final int SEARCH_UNIFORM = 3;
    public static final int SEARCH_SUFFIX_TREE = 4;

    public int nDims = 0;
    public int nMinMatches = 3;
    public int nMaxInst = -1; // no max
    public int nMaxPats = -1; // no max
    public double aTemporal = 0.05;
    public double aSeries = 0.05;
    public int nSamples = 1000;
    public String sPath = ".";
    public int nwLen = 30;
    public int nwSkip = 10;
    public int nwGrow = 5;
    public int nMarkers = 3;
    public int kSearch = SEARCH_EXHAUSTIVE;
    public ArrayList<String> filePats;
    public int iExSeries = -1;
    public int iExStart = -1;
    public int knn = 3;
    public int splitW = 0;
    public double splitT = .01;
    public int splitG = 0;
    public int splitO = 0;
    public int model = MODEL_OATES;
    public int nGarbage = 0;
    public int nHmmStates = 5;
    public int nHmmSkip = 1;
    public int hmmEval = HMM_VITERBI;
    public int hmmTrain = HMM_BW;
    public boolean bUpdateVar = true;
    public Gaussian1D gmSeries[];
    public FeatureVec initv, minv;
    public WindowLocation given;

    public Params()
    {
        filePats = new ArrayList<String>();
    }

    public void addFilePat(String filePat)
    {
        filePats.add(filePat);
    }
}
