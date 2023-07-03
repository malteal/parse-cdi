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
    >> spaces()
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
  ( spaces()
    >> string("central_value(")
    >> (floating() + (comma >> floating()))
    << string(")")
  )


sys = \
  joint \
  ( spaces()
    >> string("sys(")
    >> tostring(many(none_of(",")))
    << comma
  , floating()
    << string("%)")
  )


sfbin = \
  joint \
  ( ptbinning << spaces() << string("{") << spaces()
  , cv
  , spaces() >> many(sys) << spaces() << string("}")
  )
  

sffile = begin >> many(sfbin)


if __name__ == "__main__":
  from sys import argv
  with open(argv[1]) as reader:
    s = reader.read()

  print(sffile(s, 0))
