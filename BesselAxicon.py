from diffractio import degrees, mm, nm, np, plt, sp, um

from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.scalar_fields_XY import Scalar_field_XY

from diffractio.vector_sources_XY import Vector_source_XY
from diffractio.vector_masks_XY import Vector_mask_XY
from diffractio.vector_fields_XY import Vector_field_XY
import math
x0 = np.linspace(-125 * um, 125 * um, 256)
y0 = np.linspace(-125 * um, 125 * um, 256)

wavelength = 0.6328 * um

u0 = Scalar_source_XY(x0, y0, wavelength)
u0.bessel_beam(A = 1, r0=(0, 0), alpha=math.degrees(5), n =0, theta=0. * degrees, phi=0 * degrees, z0=0)

print(u0)


EM0 = Vector_source_XY(x0, y0, wavelength)
EM0.azimuthal_wave(256, r0=(0, 0), radius=200)
EM0.draw(kind='ellipses')
plt.title('Before mask')
plt.savefig('usage12.png')
