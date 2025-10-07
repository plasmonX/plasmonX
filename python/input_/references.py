#!/usr/bin/env python

def define_references_to_be_printed(data, atomtypes):

    # crea l'array 'cit' come nella tua logica (per indici 0..9)
    cit = [None] * 11

    #T. Giovannini et al. arXiv.2025 
    cit[10] = """ - plasmonX 
    https://doi.org/"""

    #T. Giovannini et al. Nano Lett. 2025, 25, 10802. 
    cit[9]  = """ - wFQFMu [if picocavities]
    https://doi.org/10.1021/acs.nanolett.5c01999 
"""
    #P. Grobas Illobre et al. J. Chem. Phys. 2025, 162, 044103. 
    cit[8]  = """ - BEM
    https://doi.org/10.1063/5.0245629 
"""
    #L. Nicoli et al. Fronnt. Photonics 2023, 4, 1199598. 
    cit[7]  = """ - wFQFMu [heterostructures]
    https://doi.org/10.3389/fphot.2023.1199598 
"""
    #S. Zanotto, et al. ACS Photonics 2023, 10, 394.
    cit[6]  = """ - wFQ [if polycrystallyne graphene]
    https://doi.org/10.1021/acsphotonics.2c01157 
"""
    #T. Giovannini et al. ACS Photonics 2022, 9, 3025. 
    cit[5]  = """ - wFQFMu 
    https://doi.org/10.1021/acsphotonics.2c00761 
"""
    #P. Lafiosca et al. J. Phys. Chem. C 2021, 125, 23848. 
    cit[4]  = """ - GMRES algorithm
    https://doi.org/10.1021/acs.jpcc.1c04716 
"""
    #T. Giovannini et al. J. Phys. Chem. Lett. 2020, 11, 7595. 
    cit[3]  = """ - wFQ [graphene]
    https://doi.org/10.1021/acs.jpclett.0c02051 
"""
    #T. Giovannini, et al. Nanoscale 2019, 11, 6004. 
    cit[2]  = """ - wFQ 
    https://doi.org/10.1039/C8NR09134J 
"""
    #T. Giovannini, et al. J. Chem. Theory Comput. 2019, 15, 2233. 
    cit[1]  = """ - FQFMu forcefield
    https://doi.org/10.1021/acs.jctc.8b01149 
"""
    #T. Giovannini, et al. Chem. Soc. Rev. 2020, 49, 5664. 
    cit[0]  = """ - FQ forcefield
    https://doi.org/10.1039/C9CS00464E 
"""

    citations = []
    algorithm = (data or {}).get("algorithm", {}) or {}
    forcefield = (data or {}).get("forcefield", {}) or {}
    bem = (data or {}).get("bem", {}) or {}
    static_forcefield = (forcefield.get("static", "none") or "none").lower()
    dynamic_forcefield = (forcefield.get("dynamic", "none") or "none").lower()

    # mantengo la tua logica: se serve un singolo atomtype uso il primo

    atomtypes_lower = [str(a).strip().lower() for a in atomtypes]    

    #if static_forcefield == "fq" : 
    #    citations.append(cit[0])
    if static_forcefield == "fqfmu":
        citations.append(cit[1])

    if dynamic_forcefield == "wfq" : 
        if atomtypes_lower[0] == "c": 
            citations.append(cit[3])
            citations.append(cit[6])
        else: 
            citations.append(cit[2])
    elif dynamic_forcefield == "wfqfmu" : 
        if len(atomtypes) > 1 : #hetero
            citations.append(cit[5])
            citations.append(cit[7])
        else: 
            citations.append(cit[5])
            citations.append(cit[9])

    if bem != {}: 
        citations.append(cit[8])

    if (algorithm.get("method", "")).lower() in ["iterative", "iterative on the fly"]:
        citations.append(cit[4])

    citations.append(cit[10])

    return citations 
