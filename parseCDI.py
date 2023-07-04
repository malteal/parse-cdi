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


import matplotlib.figure as figure
import numpy
from cpplot.cpplot import comparehist, zeroerr


def plotvars(binning, nominal, variations):
  labels, vars = zip(*list(variations.items()))
  labels = list(labels)

  nominal = numpy.array(nominal)
  vars = list(map(lambda h: (1 + numpy.array(h))*nominal, vars))

  nominal = zeroerr(nominal)
  vars = list(map(zeroerr, vars))

  binning = numpy.array(binning)
  
  return comparehist([nominal] + vars[:5], binning, ["nominal"] + labels, xlabel="", ylabel="")


if __name__ == "__main__":
  from sys import argv
  with open(argv[1]) as reader:
    s = reader.read()

  cdi = sffile(s, 0).value
  nominal, vars = tohists(cdi)
  bins = binning(cdi)

  fig = plotvars(bins, nominal, vars)
  fig.savefig("test.pdf")
