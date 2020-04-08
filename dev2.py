import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)

def ConctoDen(x):
    return x/2
def DentoConc(x):
    return x

x = range(50)
y = [n**10 for n in x]

fig, ax = plt.subplots()
ax.plot(x, y)

ax2 = ax.twinx()
ax2.plot(x, [ConctoDen(n) for n in y])
ax.figure.canvas.draw()

offset = ax2.yaxis.get_major_formatter().get_offset()
ax2.set(ylabel="Species density (" + offset + "cm$^{-3}$)")
print(offset, type(offset))

scale_x = float("1e"+offset[11:13])
print(scale_x)
ticks_x = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x/scale_x))
ax2.yaxis.set_major_formatter(ticks_x)

#ax1.yaxis.set_major_formatter(ScalarFormatter(useOffset=True))
#ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#ax1.yaxis.set_major_formatter(tick.Formatter().set_scientific((-3,4)))
#ScalarFormatter(useOffset=False)


#ax.yaxis.offsetText.set_visible(False)
#ax.yaxis.set_major_formatter()
#tick.FormatStrFormatter("%f")
#ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=1e16))

#, useOffset=-1e16)
#ax.yaxis.set_major_formatter(ScalarFormatter(use_scientific(False)))



#secax = ax.secondary_yaxis('right', functions=(ConctoDen,DentoConc))


#ax.yaxis.set_major_formatter(ticks_x)
#secax.yaxis.set_major_locator(AutoMinorLocator())


#ax.ticklabel_format(axis='y', style='plain')
#o = ax.yaxis.get_major_formatter().get_offset()


plt.show()