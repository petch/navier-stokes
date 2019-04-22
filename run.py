# from problems.channelflow.simulation import *
# from problems.kanflow.simulation import *
# from problems.lucasflow.simulation import *
# from problems.poiseuilleflow.simulation import *
# from problems.shiftedvortex.simulation import *
from problems.taylorvortex.simulation import *

# from problems.dynamicvortex.simulation import *
# from problems.taylordynamicvortex.simulation import *

Domain().generate()
Nonlinear(Problem()).solve()

# Simulation().run()

