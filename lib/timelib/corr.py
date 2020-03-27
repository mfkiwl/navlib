# Time from reference epoch
# Correction for beginning or end of week crossovers in GNSS systems
def dt(t, t0):

    t = t - t0

    if t > 302400:
        t = t - 604800
    elif t < -302400:
        t = t + 604800

    return t
