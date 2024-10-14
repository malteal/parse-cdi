import pyrootutils

root = pyrootutils.setup_root(search_from=__file__, pythonpath=True)

import numpy
from cpplot import compare, stderr
from matplotlib import figure
from CDI import sffile
from glob import glob
from tools.misc import load_json
from tools.visualization import atlas_utils
import matplotlib.pyplot as plt

def binning(xs):
  return list(map(lambda x : x[0][0], xs)) + [ xs[-1][0][1] ]

def get_WP(value):
    if value==85:
        return "0_67"
    elif value==77:
        return "2_2"
    elif value==70:
        return "3_25"
    elif value==60:
        return "4_57"
    else:
        return f"cut at {value}"

# this is bad.
# it assumes all keys are there properly in all bins...
def tohists(xs):
  nmax = len(xs)
  allkeys = xs[0][2].keys()

  centralvalue = [ x[1][0] for x in xs ]
  allvars = { k : [ xs[i][2][k] for i in range(nmax) ] for k in allkeys }


  return centralvalue, allvars


def dictmap(f, d):
  return { k : f(x) for k , x in d.items() }

# MODEL_PATH = "/home/users/a/algren/scratch/trained_networks/ftag_calib/paper_runs/calibrations/ttbar/2024_06_03_22_49_18_190896_lights_z_plus_jet/plots/Flow_uniform/integral/dl1r_b/"
MODEL_PATH = "/home/users/a/algren/scratch/trained_networks/ftag_calib/paper_runs/calibrations/ttbar/2024_06_03_22_49_18_190896_lights_z_plus_jet/plots/Flow_all/integral/dl1r_b/"
GLOB_PATH = "/home/users/a/algren/work/root_stuff/parse_cdi"
if __name__ == "__main__":
    wp_lst = [60, 70, 77, 85]
    addition_name= "_uniform" if "uniform" in MODEL_PATH else ""

    for wp in wp_lst:
        vals_wp = get_WP(wp)

        wp_ot = glob(f"{MODEL_PATH}/only_sig_{vals_wp}*.json")[0]
        
        otsfs = load_json(wp_ot)
        
        transport_name = [i for i in otsfs if "widehat{T}" in i][0]
        nominal_name = [i for i in otsfs if ("widehat{T}" not in i) & ("$b$-jets" in i)][0]


        #   wp = int(argv[1])
        fname =  \
        f"{GLOB_PATH}/data/btag_ttbarPDF_mc16ade_v1.0_21-2-93_DL1r_AntiKt4EMPFlowJets_BTagging201903_{wp}.txt"
        with open(fname) as reader:
            s = reader.read()

        cdi = sffile(s, 0).value
        nominal, vars = tohists(cdi)
        bins = numpy.array(binning(cdi))

        def app(h):
            return (1 + numpy.array(h))*nominal

        nominal = numpy.array(nominal)

        dvars = dictmap(app, vars)

        plotvarnames = \
        [ "Flavor_Composition"
        , "Flavor_Response"
        , "EtaIntercalibration_Modelling"
        , "FT_EFF_ttbar_PowHW7"
        , "ttbar_mc_rad"
        , "Pileup_OffsetMu"
        , "Pileup_RhoTopology"
        , "JER_DataVsMC_MC16"
        , "JER_EffectiveNP_1"
        ]

        vars = [ dvars[k] for k in dvars ]
        statvars = [ dvars[k] for k in dvars if "stat" in k ]
        plotvars = [ dvars[k] for k in plotvarnames ]
        systvars = \
        [ dvars[k] for k in dvars
            if "stat" not in k
                and "singletop" not in k
                and "PDF" not in k
                and "Light" not in k
        ]

        alluncerts = stderr(nominal, vars)
        withstatuncerts = stderr(nominal, statvars)
        withsysuncerts = stderr(nominal, systvars)
        otcalibuncerts = stderr(nominal, plotvars)

        from maltesfs import sfs # csvsfs

        # for whatever reason, numpy removes "-" from recarray key names...
        uncerts = \
        [ "Flavor_Compositiondown"
        , "Pileup_OffsetMudown"
        , "JER_EffectiveNP_1down"
        , "Flavor_Responsedown"
        , "JER_DataVsMC_MC16down"
        , "EtaIntercalibrationdown"
        , "Pileup_RhoTopologydown"
        , "20230622_145054082514_hdamp"
        , "20230828_144908518131_lights_on_nominal"
        , "20230831_113912209609_herwig_pythia_MC"
        , "stats"
        ]


        # otsfs = sfs[transport_name]
        xs = numpy.array(otsfs["bins"])
        xerrs = numpy.diff(otsfs["bins"])/2
        ys = numpy.array(otsfs[transport_name])/numpy.array(otsfs[nominal_name])
        xscenters = (xs[1:] + xs[:-1]) / 2.0

        # yuncerts = { k : numpy.sqrt(otsfs[k].T) for k in uncerts }
        # yerrs = numpy.sqrt(numpy.sum([ otsfs[k].T for k in uncerts ], axis=0))
        # yerrs = numpy.stack([yerrs]*2)
        
        yerrs_stat = numpy.array(otsfs["stat_unc"])/numpy.array(otsfs[nominal_name])
        yerrs_total = numpy.array(otsfs["full_unc"])/numpy.array(otsfs[nominal_name])


        bins = numpy.array(bins)
        bincenters = (bins[1:] + bins[:-1]) / 2.0
        binerrs = (bins[1:] - bins[:-1]) / 2.0



        xzeros = xerrs
        # yzeros = numpy.zeros_like(yerrs_stat)
        curves = \
        [ 
        ( ys , yerrs_stat),
        ( ys , yerrs_total )
        # , ( ys + yuncerts["20230831_113912209609_herwig_pythia_MC"] , yzeros )
        # , ( ys + yuncerts["20230828_144908518131_lights_on_nominal"] , yzeros )
        # , ( ys + yuncerts["Flavor_Compositiondown"] , yzeros )
        # , ( ys + yuncerts["Pileup_RhoTopologydown"] , yzeros )
        # , ( ys + yuncerts["JER_EffectiveNP_1down"] , yzeros )
        ]

        curvelabels = \
        [ transport_name
        , "Stat. + Syst."
        # , "parton shower uncertainty"
        # , "light calib uncertainty"
        # , "JES flavor composition"
        # , "JES rho topology"
        # , "JER effective NP 1"
        ]

        nextracurves = len(curves)-2

        fig = figure.Figure((4.5, 4.5))
        plt = fig.add_subplot(111)
        compare \
        ( plt
        # , [ ( xscenters , xzeros ) ] * len(curves) + [ ( bincenters , binerrs ) ]
        , [ ( xs , xzeros ) ] * len(curves) + [ ( bincenters , binerrs ) ]
        , curves + [ alluncerts ]
        , curvelabels + [ "Likelihood fit calibration" ]
        , "$p_\\mathrm{T}$ [GeV]"
        , "efficiency scale factor"
        , alphas = [0.6,  0.3, 1 ]
        , errorfills = [ True ] * len(curves) + [ False ]
        , markers = [ None ] * len(curves) + [ "s" ]
        , linewidths = [ 2 , 1 ] + [ 2 ] * nextracurves + [ 0 ]
        , colors= \
            [ "red" ] * 2
            + [ "blue" , "green" , "gray" , "magenta" , "orange" ][:nextracurves]
            + [ "black" ]
        )

        plt.set_xscale("log")
        plt.set_ylim(0.75, 1.05)
        plt.set_xlim(20.1, 399)
        xloc=57.5
        yloc=0.92
        plt.text(xloc, yloc+0.015, "ATLAS", weight="bold", style="italic", size="medium")
        # plt.text(xloc, yloc, "Internal", style="italic", size="medium")
        plt.text(xloc, yloc, "Preliminary", style="italic", size="medium")
        plt.text(xloc, yloc-0.02, "$\\sqrt{s} = 13$ TeV", size="medium")
        
        plt.legend(loc="lower center", title=f'{atlas_utils.dl1rb} {wp} %',
                   frameon=False)

        fig.tight_layout()
        fig.savefig(f"figs/{wp}WP{addition_name}.pdf")