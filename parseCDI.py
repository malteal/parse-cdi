from parsec import *


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


import matplotlib.figure as figure
import numpy
from cpplot.cpplot import comparehist, zeroerr, stderr


if __name__ == "__main__":
  from sys import argv

  wp = int(argv[1])
  fname =  \
    "btag_ttbarPDF_mc16ade_v1.0_21-2-93_DL1r_AntiKt4EMPFlowJets_BTagging201903_%s.txt" \
    % argv[1]

  with open(fname) as reader:
    s = reader.read()

  cdi = sffile(s, 0).value
  nominal, vars = tohists(cdi)
  bins = numpy.array(binning(cdi))

  def app(h):
    return (1 + numpy.array(h))*nominal

  nominal = numpy.array(nominal)

  dvars = distmap(app, vars)

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

  from maltesfs import sfs

  otsfs = numpy.array(sfs[wp][0]) , numpy.array(sfs[wp][1])

  fig = \
    comparehist \
    ( [ alluncerts, otsfs ]
    , numpy.array(bins)
    , [ "standard calib" , "MALTE's AWESOME CALIB\nDON'T NEED NO UNCERTAINTIES" ]
    , xlabel="jet $p_T$ / GeV"
    , ylabel="efficiency scale factor"
    , alphas=[0.7, 0.8]
    )

  plt = fig.axes[0]
  plt.legend()
  plt.set_xscale("log")
  fig.suptitle("%2d%% OP" % wp)
  fig.savefig("test.pdf")
