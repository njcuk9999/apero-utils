import numpy as np
from astropy.io import fits

def fix_shift(im):
    # reading a file
    #
    # Status:
    #  status = 0 -> all good, no shit
    #  status = 1 -> normal shift in expected direction when we have a controler issue,
    #                  we apply correction
    #  status = -1 -> opposite direction of known shift pattern, we apply correction
    #  status = -2 -> unknown shift, we need to look at the image

    # getting image size
    number_amps = 32
    width_image = image.shape[0]

    # get width of amplifier ribbon
    width_amp = width_image/number_amps

    # look for a discontinuity between consecutive columns, this traces reference vs science pixels
    gap = np.zeros(15)
    for i in range(15):
        gap[i] = np.nanmedian(im[:, i] - im[:, i + 1])

    # If all is fine, we have 4 ref pixels and the 3rd difference shows a glitch
    imax = np.argmax(np.abs(gap))

    # When images are bad, we have a value of 2 and know we need to offset odd/even amplifier
    if imax == 2:
        print('we have a shift')

        for i in range(number_amps):
            offset = 1-(i % 2)*2
            im[:,i*width_amp:i*width_amp+width_amp] = np.roll(im[:,i*width_amp:i*width_amp+width_amp],offset,axis = 1)

        status = 1

    if imax == 3:
        print('all good, no shift')

        status = 0

    if imax == 4:
        print('we have a shift, but in the direction opposite to what we saw in 20200909 dataset')

        for i in range(number_amps):
            offset = (i % 2)*2-1
            im[:,i*width_amp:i*width_amp+width_amp] = np.roll(im[:,i*width_amp:i*width_amp+width_amp],offset,axis = 1)

        status = -1

    if (imax !=2)*(imax!=3)*(image!=4):
        print('really bad, rips the fabric of the Universe!')

        status = -2

    return im, status