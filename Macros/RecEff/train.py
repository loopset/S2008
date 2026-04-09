import numpy as np
from numpy.typing import NDArray
import ROOT

ROOT.EnableImplicitMT()

rdf = ROOT.RDataFrame("GETTree", ["../../RootFiles/Cluster/Clusters_Run_0031.root"])  # type: ignore

ROOT.gInterpreter.Declare("""
    std::vector<std::array<float,3>> define_simple(ActRoot::TPCData& tpc){
        return {{-1, -1, -1}};
    }
""")
df = rdf.Define("SimpleData", "define_simple(TPCData)")

asnp = df.AsNumpy(["SimpleData"])
