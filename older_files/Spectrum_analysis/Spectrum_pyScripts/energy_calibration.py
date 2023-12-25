def energy_calibration(channel,slope=0.293676,increment=1.605054):
    '''
    Return corresponding energy for a given channel, with energy calibration slope and increment specified

    :param channel: --mandatory
    :param slope: --Optional
    :param increment --Optional
    :return energy in keV
    '''
    return channel*slope+increment