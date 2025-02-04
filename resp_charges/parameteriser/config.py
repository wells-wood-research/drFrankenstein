import dataclasses


@dataclasses
class Default_Config():
    """Default config for resp""" 
    functional = "HF"
    basis_set= "6-31G*"
    n_configs_resp  = 20 
    nconfigs_per_dihedral = 3
