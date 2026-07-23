import numpy as np


def calculate_ribosome_density(tab_single, tab_polysome, name_single,
                               name_polysome, L_poi, L_tag):
    '''
    Calculate ribosome density based on the mean intensity of single protein and the length of the protein of interest and the tag.

    Parameters
    ----------
    tab_single : pd.DataFrame
        Dataframe containing the single protein fluorescence data.
    tab_polysome : pd.DataFrame
        Dataframe containing the polysome fluorescence data. 
    name_single : str
        Name of the column in tab_single that contains the single protein fluorescence data.
    name_polysome : str
        Name of the column in tab_polysome that contains the polysome fluorescence data.
    L_poi : int
        Length of the protein of interest in amino acids.
    L_tag : int
        Length of the tag in amino acids.
    
    Returns
    -------
    m_intensity_single : float
        Mean intensity of the single protein.
    tab_polysome : pd.DataFrame
        Updated dataframe containing the polysome fluorescence data with an additional column for ribosome density.

    '''

    l_prime = L_poi+0.5*L_tag
    m_intensity_single = np.mean(tab_single[name_single])
    tab_polysome["ribosome_density"] = tab_polysome[name_polysome]/(
        m_intensity_single*l_prime)
    return m_intensity_single, tab_polysome


def estimate_density(mean_n_ribosome, prot_length, suntag_length, footprint=10):
    """
    Generate n tracks according to one protein translation dynamics

    Parameters
    ----------
    mean_n_ribosome : int
        mean occupancy density of the mRNA, from a measured mean number of active ribosomes (e.g. y_start_prot.mean() from simulation, or mean_ribosome_from_intensity() for real data).
    prot_length : int
        length of the protein in amino acid
    suntag_length : int
        length of the suntag in amino acid
    footprint : int
        ribosome exclusion size in amino acids (default 10)
    Returns
    -------
    density : float
        estimated ribosome density in ribosome per amino acid   


    Description
    -----------
    This function generate one track of a translation site to mimic in vivo
    translation profile.
    It is based on one protein translation profile.

    """

    L = prot_length + suntag_length
    return (mean_n_ribosome * footprint) / L
