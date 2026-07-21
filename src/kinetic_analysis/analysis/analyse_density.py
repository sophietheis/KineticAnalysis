import numpy as np


def calculate_ribosome_density(tab_single, tab_polysome, name_single,
                               name_polysome, L_poi, L_tag):
    '''

    :param tab_single:
    :type tab_single:
    :param tab_polysome:
    :type tab_polysome:
    :param L_poi:
    :type L_poi:
    :param L_tag:
    :type L_tag:
    :return:
    :rtype:
    '''

    l_prime = L_poi+0.5*L_tag
    m_intensity_single = np.mean(tab_single[name_single])
    tab_polysome["ribosome_density"] = tab_polysome[name_polysome]/(
            m_intensity_single*l_prime)
    return m_intensity_single, tab_polysome

def estimate_density(mean_n_ribosome, prot_length, suntag_length, footprint=10):
    """
    rho_bar = mean occupancy density of the mRNA, from a measured mean
    number of active ribosomes (e.g. y_start_prot.mean() from simulation,
    or mean_ribosome_from_intensity() for real data).
    """
    L = prot_length + suntag_length
    return (mean_n_ribosome * footprint) / L