import pyrootutils

root = pyrootutils.setup_root(search_from=__file__, pythonpath=True)

from parsec import *
import numpy as np
from cpplot import compare, stderr
from matplotlib import figure
import matplotlib.pyplot as plt
from tools import misc
import matplotlib
from tools.visualization.atlas_utils import ATLAS_setup


font = {#'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


def tostring(p):
  return p.parsecmap("".join)

newline = string("\n")
char = none_of("")
line = many(none_of("\n")) << newline
digits = tostring(many(digit()))
integer = digits.parsecmap(int)
comma = string(",")


def floating():
  neg = string("-")
  start = digits
  p = string(".")
  rest = digits
  q = \
    joint(neg , start , p , rest) \
    ^ joint(start , p , rest) \
    ^ joint(start , p) \
    ^ start.parsecmap(lambda x: [x])

  return q.parsecmap("".join).parsecmap(float)


metaline = spaces() >> string("meta") >> line

begin = \
  ( string("Analysis")
    >> line
    >> many(metaline)
  )


ptbinning = \
  joint \
  ( spaces()
    >> string("bin(")
    >> floating()
    << string("<pt<")
  , floating() << line
  )


cv = \
  joint \
  ( spaces()
    >> string("central_value(")
    >> floating()
  , comma >> floating() << string(")")
  )


sys = \
  joint \
  ( spaces()
    >> string("sys(")
    >> tostring(many(none_of(",")))
    << comma
  , floating().parsecmap(lambda x : x * 0.01) << string("%)")
  )


sfbin = \
  joint \
  ( ptbinning << spaces() << string("{") << spaces()
  , cv
  , (spaces() >> many(sys) << spaces() << string("}")).parsecmap(dict)
  )
  

sffile = begin >> many(sfbin)


def binning(xs):
  return list(map(lambda x : x[0][0], xs)) + [ xs[-1][0][1] ]


# this is bad.
# it assumes all keys are there properly in all bins...
def tohists(xs):
  nmax = len(xs)
  allkeys = xs[0][2].keys()

  centralvalue = [ x[1][0] for x in xs ]
  allvars = { k : [ xs[i][2][k] for i in range(nmax) ] for k in allkeys }


  return centralvalue, allvars


def distmap(f, d):
  return { k : f(x) for k , x in d.items() }


if __name__ == "__main__":
  # from sys import argv
  save_fig = True

  # wp = 60
  for wp in [60, 70, 77, 85]:
    fname = \
      "data/btag_ttbarPDF_mc16ade_v1.0_21-2-93_DL1r_AntiKt4EMPFlowJets_BTagging201903_%s.txt" \
      % wp

    with open(fname) as reader:
      s = reader.read()

    cdi = sffile(s, 0).value
    nominal, vars = tohists(cdi)
    bins = np.array(binning(cdi))

    def app(h):
      return (1 + np.array(h))*nominal

    nominal = np.array(nominal)

    dvars = distmap(app, vars)
    
    bins = np.array(bins)
    bincenters = (bins[1:] + bins[:-1]) / 2.0
    binerrs = (bins[1:] - bins[:-1]) / 2.0

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

    # vars = {k: dvars[k] for k in dvars if "FT" not in k}
    vars = {k: dvars[k] for k in dvars if "stat" not in k and "singletop" not in k }
    # statvars = [ dvars[k] for k in dvars if "stat" in k ]
    # statvars_key = {k:dvars[k] for k in dvars if "stat" in k}
    # plotvars = [ dvars[k] for k in plotvarnames ]
    # systvars = \
    #   [ dvars[k] for k in dvars
    #       if "stat" not in k
    #         and "singletop" not in k
    #         and "PDF" not in k
    #         and "Light" not in k
    #   ]

    # allvars = {k: dvars[k] for k in dvars}
    # sigma_total
    alluncerts = stderr(nominal, vars.values())
    
    # sigma_vars / sigma_total
    fracuncerts = { k : stderr(nominal, [vars[k]])[1] / alluncerts[1] for k in vars }

    # uncerts_contribution = {}

    # for i,j in vars.items():
      
    #   uncerts_contribution[i] = np.sqrt(alluncerts[-1]**2-stderr(nominal, [j])[-1]**2)/alluncerts[-1]
      
    # diff_to_max = {i: np.max(np.abs(1-j)) for i,j in uncerts_contribution.items() if np.max(np.abs(1-j)) >= 0.001}

    # ranking_nr = np.argsort(list(diff_to_max.values()))[::-1]
    
    # ranking_keys = np.array(list(diff_to_max.keys()))
    
    # nr_sys= 15 # addtional systematics to plot

    # scale=2
    # fig, ax1 = plt.subplots(1,1, figsize = (8*scale, 6*scale), sharey=True)
    # st = fig.suptitle(f"WP: {wp} %", fontsize="x-large")
    # for nr in ranking_nr[:nr_sys]:
    #   i = ranking_keys[nr]
    #   ax1.plot(bincenters, uncerts_contribution[i], label=f"{i}, max ratio: {np.round(diff_to_max[i],3)}", marker="s", linestyle="--")
    # ax1.legend(frameon=False)
    # ax1.set_ylabel("$\sigma_{n-1}$/$\sigma_{nominal}$", fontsize=30)
    # # ax1.set_ylabel("SF$_{sys}$/SF$_{nominal}$")
    # ax1.set_xlabel("pT", fontsize=30)
    # plt.tight_layout()
    # if save_fig:
    #   misc.save_fig(fig, f"sys_uncertainty_contribution_{wp}WP.pdf")
      
    diff_to_max = {i: np.max(np.abs(j)) for i,j in fracuncerts.items()}# if np.max(np.abs(j)) >= 0.0031}

    ranking_nr = np.argsort(list(diff_to_max.values()))[::-1]
    
    ranking_keys = np.array(list(diff_to_max.keys()))
    
    nr_sys= 15 # addtional systematics to plot

    scale=1.1
    fig, ax1 = plt.subplots(1,1, figsize = (9*scale, 8*scale), sharey=True)
    # st = fig.suptitle(f"WP: {wp} %", fontsize="x-large")
    total_selected_unc = None
    for nr in ranking_nr[:nr_sys]:
      i = ranking_keys[nr]
      # ax1.plot(bincenters, np.abs(fracuncerts[i]), label=f"{i}, max ratio: {np.round(diff_to_max[i],3)}", marker="s", linestyle="--")
      label=i.replace("_", " ")
      ax1.errorbar(bincenters, np.abs(fracuncerts[i]),
                   label=f"{label}", # max ratio: {np.round(diff_to_max[i],3)}",
                  #  label=f"{label} max ratio: {np.round(diff_to_max[i],3)}",
                   xerr = binerrs)

      if total_selected_unc is None:
        total_selected_unc= np.abs(fracuncerts[i])**2
      else:
        total_selected_unc+= np.abs(fracuncerts[i])**2
        
    ax1.errorbar(bincenters, np.sqrt(total_selected_unc), label=r"$\sum \sigma_{sys}$", marker="s", linestyle="--", color="black", xerr = binerrs)
    ax1.errorbar(bincenters, np.ones_like(bincenters), linestyle="-", color="black",
                 xerr = binerrs)
    ax1.legend(frameon=False, ncol=2, loc="center right", bbox_to_anchor=(1, 0.675))
    ax1.set_ylabel("$\sigma_{sys}$/$\sigma_{total}$", fontsize=30)
    ax1.set_ylim([0, 1.5])
    ax1.set_xlim([20, 400])
    
    ATLAS_setup(ax=ax1, xlabel=r"$p_\mathbf{T}$ [GeV]",
                ylabel=r"$\sigma_{variation}/\sigma_{total}$",
                atlas_kw={"font_size": 20, "sub_font_size": 20, "label_font_size": 20,
                          "line_spacing": 1.5},
                text=f"WP: {wp} %, "
                )

    plt.tight_layout()
    if save_fig:
      misc.save_fig(fig, f"contribution_of_{nr_sys}_{wp}WP.pdf")
      
    import sys
  sys.exit()

  # from maltesfs import csvsfs
  
  # for whatever reason, np removes "-" from recarray key names...
  # uncerts = \
  #   [ "Flavor_Compositiondown"
  #   , "Pileup_OffsetMudown"
  #   , "JER_EffectiveNP_1down"
  #   , "Flavor_Responsedown"
  #   , "JER_DataVsMC_MC16down"
  #   , "EtaIntercalibrationdown"
  #   , "Pileup_RhoTopologydown"
  #   , "20230622_145054082514_hdamp"
  #   , "20230828_144908518131_lights_on_nominal"
  #   , "20230831_113912209609_herwig_pythia_MC"
  #   , "stats"
  #   ]


  # otsfs = csvsfs[wp]
  # xs = otsfs["bins"]
  # xerrs = np.stack([otsfs["width"]]*2)
  # ys = otsfs["center"]

  # yuncerts = { k : np.sqrt(otsfs[k].T) for k in uncerts }
  # yerrs = np.sqrt(np.sum([ otsfs[k].T for k in uncerts ], axis=0))
  # yerrs = np.stack([yerrs]*2)

  # xzeros = np.zeros_like(xerrs)
  # yzeros = np.zeros_like(yerrs)
  # curves = \
  #   [ ( ys , yerrs )
  #   , ( ys , np.stack([yuncerts["stats"]]*2) )
  #   # , ( ys + yuncerts["20230831_113912209609_herwig_pythia_MC"] , yzeros )
  #   # , ( ys + yuncerts["20230828_144908518131_lights_on_nominal"] , yzeros )
  #   # , ( ys + yuncerts["Flavor_Compositiondown"] , yzeros )
  #   # , ( ys + yuncerts["Pileup_RhoTopologydown"] , yzeros )
  #   # , ( ys + yuncerts["JER_EffectiveNP_1down"] , yzeros )
  #   ]

  # curvelabels = \
  #   [ "optimal transport"
  #   , "stat uncertainty"
  #   # , "parton shower uncertainty"
  #   # , "light calib uncertainty"
  #   # , "JES flavor composition"
  #   # , "JES rho topology"
  #   # , "JER effective NP 1"
  #   ]

  # nextracurves = len(curves)-2

  fig = figure.Figure((4.5, 4.5))
  plt = fig.add_subplot(111)

 
  compare \
    ( plt
    ,  [ ( bincenters , binerrs ) ]
    ,  [ alluncerts ]
    ,  [ "standard calib" ]
    , "jet $p_\\mathrm{T}$ [GeV]"
    , "efficiency scale factor"
    # , alphas = [0.5] + [ 1.0 ] * len(curves)
    # , errorfills = [ True ] * len(curves) + [ False ]
    # , markers = [ None ] * len(curves) + [ "s" ]
    # , linewidths = [ 2 , 0 ] + [ 2 ] * nextracurves + [ 0 ]
    , colors= \
        [ "red" ] * 2
      + [ "black" ]
    )

  plt.set_xscale("log")
  plt.set_ylim(0.75, 1.05)
  plt.text(60, 0.90, "ATLAS", weight="bold", style="italic", size="large")
  plt.text(60, 0.885, "Internal", style="italic", size="large")
  plt.text(60, 0.87, "$\\sqrt{s} = 13.6$ TeV", size="large")
  plt.legend(loc="lower center")

  fig.tight_layout()
  fig.savefig("%2dWP.pdf" % wp)
